import subprocess
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import os
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import argparse

import warnings
import time
from concurrent.futures import ProcessPoolExecutor
start_time = time.time()
warnings.filterwarnings("ignore")

print("\n========================================\n")
print('Setting up the environment...')

# Parse input arguments
parser = argparse.ArgumentParser(description="Run asmooth and optional for creating new cheese-mask.")
parser.add_argument("input_image", type=str, help="Path to the input image file.")
parser.add_argument("input_expmap", type=str, help="Path to the input exposure map file.")
# parser.add_argument("output_dir", type=str, help="Directory to save the output files.")
# parser.add_argument("PS_size", type=float, default=1, help="Point source extent size in arcmin.")
parser.add_argument("desired_snr", type=int, default=30, help="Desired signal-to-noise ratio for asmooth.")
parser.add_argument("cheesemask_file", type=str, help="Path to the cheese-mask file.")
parser.add_argument("--cheesemask_regions", type=str, help="Path to the cheese-mask region file. Required for creating and working with a new cheese-mask")
parser.add_argument("--detmask_file", type=str, help="Path to the detection mask file.")
parser.add_argument("--new_cheesemask", action="store_true", default=False, help="Flag to create a new cheese-mask or overwrite the existing one. Default is False.")
args = parser.parse_args()

# Input parameters
input_image = args.input_image
input_expmap = args.input_expmap
# output_dir = args.output_dir
input_cheesemask = args.cheesemask_file
input_cheesemask_regions = args.cheesemask_regions
new_cheesemask = args.new_cheesemask
detmask_file = args.detmask_file
desired_snr = args.desired_snr

log_file_path = os.path.join(os.path.dirname(input_cheesemask), "process.log")

# Read the log file and delete everything after "catprep: DONE"
if os.path.exists(log_file_path):
    with open(log_file_path, "r") as log_file:
        lines = log_file.readlines()
    
    with open(log_file_path, "w") as log_file:
        for line in lines:
            log_file.write(line)
            if "catprep: DONE" in line:
                break

print('\nEnvironment setup complete\n')
print('Starting the workflow with the following parameters:')
print(f'    Input image: {input_image}')
print(f'    Input exposure map: {input_expmap}')
print(f'    Desired signal-to-noise ratio for asmooth: {desired_snr}')
print(f'    Input Cheese-mask file: {input_cheesemask}')
if input_cheesemask_regions:
    print(f'    Remaking the cheese-mask file using the regions file: {input_cheesemask_regions}')
    print(f'    Input Detection mask file: {detmask_file}')
    if new_cheesemask:
        print(f'    new_cheesemask flag: {new_cheesemask}. Thus creating a new cheese-mask file')
    else:
        print(f'    new_cheesemask flag: {new_cheesemask}. Thus overwriting the existing cheese-mask.')
print(f'    Log file: {log_file_path}')

########## Creating a new cheese-mask ##########
if input_cheesemask_regions:
    print("\n========================================\n")
    print(f'Re-writing the cheese-mask with the new regions from {input_cheesemask_regions}...')
    if not detmask_file:
        print("ERROR: Detection mask file is required for creating a new cheese-mask.")
        exit()
    if new_cheesemask:
        cheesemask_file = input_cheesemask.replace(".fits", "_new.fits")
        print(f'New cheese-mask will be saved as {cheesemask_file}')
    else:
        cheesemask_file = input_cheesemask
        print(f'Existing cheese-mask will be overwritten.')

    hdulist = fits.open(input_image)
    ima = hdulist[0].data
    prihdr = hdulist[0].header
    pix2deg = prihdr['CDELT2']  # deg
    xsize, ysize = ima.T.shape  # transpose is required because x is RA and y is DEC

    mask_hdu = fits.open(detmask_file)
    mask = mask_hdu[0].data

    ima_wcs = wcs.WCS(prihdr, relax=False)
    ima_racen, ima_deccen = prihdr['CRVAL1'], prihdr['CRVAL2']
    ima_r = np.max((xsize, ysize)) / 2 * pix2deg  # deg
    reg = open(input_cheesemask_regions).readlines()
    ra_src, dec_src, ext_src = np.zeros(len(reg)), np.zeros(len(reg)), np.zeros(len(reg))
    
    for i in range(len(reg)):
        if 'circle(' in reg[i]:
            ra_src[i] = float(reg[i].split(',')[0].replace('fk5; circle(', ''))
            dec_src[i] = float(reg[i].split(',')[1])
            ext_src[i] = float(reg[i].split(',')[2].replace(')', ''))

    def circle(X, Y):
        x, y = np.meshgrid(X, Y)
        rho = np.sqrt(x * x + y * y)
        return rho

    x = np.arange(xsize)
    y = np.arange(ysize)

    for j in tqdm(range(len(ra_src))):
        pixim = ima_wcs.all_world2pix([[float(ra_src[j]), float(dec_src[j])]], 0)
        xp = pixim[0][0]
        yp = pixim[0][1]
        rho = circle(x - xp, y - yp) * pix2deg
        ii = np.where(rho <= ext_src[j])
        if len(ii[0]) > 0:
            mask[ii] = 0

    hdu = fits.PrimaryHDU(mask)
    hdu.header.update(ima_wcs.to_header())
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(cheesemask_file, overwrite=True)

    print(f'\nCheese-mask saved as {cheesemask_file}. Using this mask for further processing.')

########## Running asmooth ##########
print("\n========================================\n")
print('Running asmooth and filling the holes...\n')

print('Multiplying the image and exposure map with the cheese-mask...\n')
output_masked_image = input_image.replace(".fits", "_masked.fits")
output_masked_expmap = input_expmap.replace(".fits", "_masked.fits")
output_asmooth_image = input_image.replace(".fits", "_asmooth.fits")

if new_cheesemask:
    cheesemask_file = input_cheesemask.replace(".fits", "_new.fits")
else:
    cheesemask_file = input_cheesemask

sh_file_content = f"""#!/bin/bash
source /science/InitScripts/iaat-xmmsas.sh # Source the xmmsas environment here
input_image={input_image}
cheesemask={input_cheesemask}
masked_image={output_masked_image}
input_expmap={input_expmap}
masked_expmap={output_masked_expmap}
output_smooth_image={output_asmooth_image}
desiredsnr={desired_snr}

farith $input_image $cheesemask $masked_image MUL clobber=yes
farith $input_expmap $cheesemask $masked_expmap MUL clobber=yes

asmooth inset=$masked_image \
    outset=$output_smooth_image \
    weightset=$masked_expmap \
    withweightset=yes \
    withexpimageset=yes \
    expimageset=$masked_expmap \
    desiredsnr=$desiredsnr \
"""

with open('run_asmooth.sh', 'w') as file:
    file.write(sh_file_content)

print('asmooth shell script created as run_asmooth.sh\n')
print('Running asmooth...')

with open(log_file_path, "a") as log_file:
    subprocess.run(["bash", "run_asmooth.sh"], stdout=log_file, stderr=log_file)

print('\nMasked image and exposure map saved as:')
print(f'{output_masked_image} and {output_masked_expmap}')
print('\nasmooth completed successfully!')
print(f'Smoothed image saved as {output_asmooth_image}')

end_time = time.time()
time_taken = end_time - start_time

print("\n========================================\n")
if time_taken < 600:
    print(f'** All tasks completed successfully in {time_taken:.2f} seconds **')
if time_taken >= 600:
    print(f'** All tasks completed successfully in {(time_taken/60):.2f} minutes **')
if time_taken >= 3600:
    print(f'** All tasks completed successfully in {(time_taken/3600):.2f} hours **')
print("\n========================================\n")