import subprocess
import numpy as np
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy import wcs
import os
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import argparse
import warnings
warnings.filterwarnings("ignore")

print("\n========================================\n")
print('Setting up the environment...')

# Parse input arguments
parser = argparse.ArgumentParser(description="Run source detection and asmooth pipeline.")
parser.add_argument("input_image", type=str, help="Path to the input image file.")
parser.add_argument("input_expmap", type=str, help="Path to the input exposure map file.")
parser.add_argument("input_catalog", type=str, help="Path to the input catalog file.")
parser.add_argument("output_dir", type=str, help="Directory to save the output files.")
parser.add_argument("PS_size", type=float, default=1, help="Point source extent size in arcmin.")
parser.add_argument("desired_snr", type=int, default=30, help="Desired signal-to-noise ratio for asmooth.")
args = parser.parse_args()

# Input parameters
input_image = args.input_image
input_expmap = args.input_expmap
output_dir = args.output_dir
input_catalog = args.input_catalog
desired_snr = args.desired_snr

output_mask_file = os.path.join(output_dir, "detmask.fits")
output_boxlist_local = os.path.join(output_dir, "boxlist_local.fits")
output_bkgmap = os.path.join(output_dir, "BG_maps", "bkg_map.fits")
output_cheesemask = os.path.join(output_dir, "BG_maps", "cheesemask_all.fits")
output_boxlist_map = os.path.join(output_dir, "boxlist_map.fits")
output_mllist = os.path.join(output_dir, "mllist.fits")
output_sourceimage = os.path.join(output_dir, "sourceimage.fits")
output_catalog = os.path.join(output_dir, "catalog.fits")
PS_size = args.PS_size

# Create necessary directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(output_dir, "BG_maps"), exist_ok=True)

# Create log file
log_file_path = os.path.join(output_dir, "process.log")
LOG_file = open(log_file_path, "w+")

print('\nEnvironment setup complete\n')
print('Starting the workflow with the following parameters:')
print(f'Input image: {input_image}')
print(f'Input exposure map: {input_expmap}')
print(f'Input point source catalog: {input_catalog}')
print(f'Output directory: {output_dir}')
print(f'Point source extent size: {PS_size} arcmin')
print(f'Desired signal-to-noise ratio for asmooth: {desired_snr}')
print(f'Log file: {log_file_path}')

########## Defining the functions ##########

# Function to run ermask
def run_ermask(exposure_map, output_mask_file, log_file=None):
    subprocess.run(["ermask", 
                    f"expimage={exposure_map}", 
                    f"detmask={output_mask_file}",
                    ], stdout=log_file, stderr=log_file)

# Function to run erbox
def run_erbox(image_file, exposure_map, detmask_file, output_boxlist, bkg_map=None, bg_image_flag="N", ecf=1, emin=200, emax=2300, log_file=None):
    if bg_image_flag == "N":
        subprocess.run(["erbox", 
                        f"images={image_file}", 
                        f"expimages={exposure_map}",
                        f"detmasks={detmask_file}",
                        f"boxlist={output_boxlist}",
                        f"emin={emin}",
                        f"emax={emax}",
                        f"bkgima_flag={bg_image_flag}",
                        f"ecf={ecf}",
                        ], stdout=log_file, stderr=log_file)
    else:
        subprocess.run(["erbox", 
                        f"images={image_file}", 
                        f"expimages={exposure_map}",
                        f"detmasks={detmask_file}",
                        f"boxlist={output_boxlist}",
                        f"emin={emin}",
                        f"emax={emax}",
                        f"bkgimages={bkg_map}",
                        f"ecf={ecf}",
                        ], stdout=log_file, stderr=log_file)

# Function to run erbackmap
def run_erbackmap(image_file, exposure_map, detmask_file, boxlist_file, output_bkgmap, output_cheesemask, emin=200, emax=2300, log_file=None):
    subprocess.run(["erbackmap", 
                    f"image={image_file}", 
                    f"expimage={exposure_map}",
                    f"detmask={detmask_file}",
                    f"boxlist={boxlist_file}",
                    f"bkgimage={output_bkgmap}",
                    f"cheesemask={output_cheesemask}",
                    f"emax={emax}",
                    "cheesemask_flag=Y",
                    "clobber=Y",
                    ], stdout=log_file, stderr=log_file)

# Function to run ermldet
def run_ermldet(image_file, exposure_map, detmask_file, boxlist_file, bkg_map, output_mllist, output_sourceimage, emin=200, emax=2300, log_file=None):
    subprocess.run(["ermldet", 
                    f"mllist={output_mllist}", 
                    f"boxlist={boxlist_file}",
                    f"images={image_file}",
                    f"expimages={exposure_map}",
                    f"detmasks={detmask_file}",
                    f"bkgimages={bkg_map}",
                    f"srcimages={output_sourceimage}",
                    "extentmodel=gaussian",
                    f"emin={emin}",
                    f"emax={emax}",
                    "srcima_flag=Y"
                    ], stdout=log_file, stderr=log_file)

# Function to run catprep
def run_catprep(input_mllist, out_catfile, log_file=None):
    subprocess.run(["catprep", 
                    f"infile={input_mllist}", 
                    f"outfile={out_catfile}"
                    ], stdout=log_file, stderr=log_file)
    
########## Workflow ##########
print("\n========================================\n")

# Run ermask
print('1) Creating the detection mask...')

if os.path.exists(output_mask_file):
    os.remove(output_mask_file)

with open(log_file_path, "a") as log_file:
    run_ermask(input_expmap, output_mask_file, log_file=log_file)

# Run erbox (local)
print("\n========================================\n")
print('2) Running erbox in local mode...')
if os.path.exists(output_boxlist_local):
    os.remove(output_boxlist_local)

with open(log_file_path, "a") as log_file:
    run_erbox(input_image, input_expmap, output_mask_file, output_boxlist_local, log_file=log_file)

# Run erbackmap
print("\n========================================\n")
print('3) Creating background map...')
with open(log_file_path, "a") as log_file:
    run_erbackmap(input_image, input_expmap, output_mask_file, output_boxlist_local, output_bkgmap, output_cheesemask, log_file=log_file)

# Run erbox (map)
print("\n========================================\n")
print('4) Running erbox in map mode...')
if os.path.exists(output_boxlist_map):
    os.remove(output_boxlist_map)

with open(log_file_path, "a") as log_file:
    run_erbox(input_image, input_expmap, output_mask_file, output_boxlist_map, output_bkgmap, log_file=log_file)

# Run ermldet
print("\n========================================\n")
print('5) Running ermldet to identify sources...')
if os.path.exists(output_mllist):
    os.remove(output_mllist)

if os.path.exists(output_sourceimage):
    os.remove(output_sourceimage)

with open(log_file_path, "a") as log_file:
    run_ermldet(input_image, input_expmap, output_mask_file, output_boxlist_map, output_bkgmap, output_mllist, output_sourceimage, log_file=log_file)

# Run catprep
print("\n========================================\n")
print('6) Saving the final catalog...')
if os.path.exists(output_catalog):
    os.remove(output_catalog)

with open(log_file_path, "a") as log_file:
    run_catprep(output_mllist, output_catalog, log_file=log_file)

# Read the log file and count occurrences of "process done"
def check_process_completion(log_file_path):
    with open(log_file_path, "r") as log_file:
        log_content = log_file.read()
        ermask_done_count = log_content.count("ermask:DONE")
        erbox_done_count = log_content.count("erbox:DONE")
        erbackmap_done_count = log_content.count("erbackmap:DONE")
        ermldet_done_count = log_content.count("ermldet:DONE")
        catprep_done_count = log_content.count("catprep:DONE")
        
        all_tasks_successful = (ermask_done_count == 1 and 
                                erbox_done_count == 2 and 
                                erbackmap_done_count == 1 and 
                                ermldet_done_count == 1 and 
                                catprep_done_count == 1)
    return ermask_done_count, erbox_done_count, erbackmap_done_count, ermldet_done_count, catprep_done_count, all_tasks_successful

ermask_done_count, erbox_done_count, erbackmap_done_count, ermldet_done_count, catprep_done_count, all_tasks_successful = check_process_completion(log_file_path)
print(f"\n========================================\n")
if all_tasks_successful:
    print("All tasks successful")
else:
    raise RuntimeError(f'Error: Some or all tasks were not completed successfully')

# Selecting point sources and creating cheese-mask
print("\n========================================\n")
print('7) Selecting point sources and creating cheese-mask...')

detmask_file = output_mask_file
pts_cat = input_catalog
cheesemask_file = os.path.join(os.path.dirname(detmask_file), f'cheesemask_PS_{PS_size}arcmin.fits')

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
ima_coord = SkyCoord(ima_racen * u.deg, ima_deccen * u.deg, frame='icrs')

cat_src = fits.open(pts_cat)[1].data
cat_src = cat_src[(cat_src.EXT == 0) & (cat_src.DET_LIKE_0 > 12)]  # Select high S/N point sources with DET_LIKE>8
coord_src = SkyCoord(cat_src.RA * u.deg, cat_src.DEC * u.deg, frame='icrs')

# Convert RA and DEC of sources to pixel coordinates
pix_coords = ima_wcs.all_world2pix(np.column_stack((cat_src.RA, cat_src.DEC)), 0)
# Only consider pix_coords that are within the image pixel bounds
valid_pix_coords_mask = (pix_coords[:, 0] >= 0) & (pix_coords[:, 0] < (xsize - 1)) & (pix_coords[:, 1] >= 0) & (pix_coords[:, 1] < (ysize - 1))
pix_coords = pix_coords[valid_pix_coords_mask]

# Create a mask to keep only sources with mask value of 1
valid_sources_mask = mask[pix_coords[:, 1].astype(int), pix_coords[:, 0].astype(int)] == 1

# Apply the mask to cat_src
cat_src = cat_src[valid_pix_coords_mask][valid_sources_mask]

ra_src = cat_src.RA
dec_src = cat_src.DEC

# Fix the masking radius to 1arcmin. Needs to be modified for more accurate analysis.
ext_src = np.zeros(len(ra_src)) + (PS_size / 60)

def circle(X, Y):
    x, y = np.meshgrid(X, Y)
    rho = np.sqrt(x * x + y * y)
    return rho

x = np.arange(ysize)
y = np.arange(xsize)
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

with open(os.path.join(output_dir, f'Point_Sources_{PS_size}arcmin.reg'), 'w') as f:
    for i in range(len(ra_src)):
        print(f"fk5; circle({ra_src[i]},{dec_src[i]},{ext_src[i]})", file=f)

print(f'Cheese-mask saved as {cheesemask_file}')
print(f'Point sources region file saved as {os.path.join(output_dir, f"Point_Sources_{PS_size}arcmin.reg")}')

# Run asmooth
print("\n========================================\n")
print('8) Running asmooth and filling the holes...')

print('Multiplying the image and exposure map with the cheese-mask...')
output_masked_image = input_image.replace(".fits", "_masked.fits")
output_masked_expmap = input_expmap.replace(".fits", "_masked.fits")

sh_file_content = f"""#!/bin/bash
source /science/InitScripts/iaat-xmmsas.sh
input_image={input_image}
cheesemask={cheesemask_file}
masked_image={output_masked_image}
input_expmap={input_expmap}
masked_expmap={output_masked_expmap}
output_smooth_image={output_masked_image}
desiredsnr={desired_snr}

farith $input_image $cheesemask $masked_image MUL clobber=yes
farith $input_expmap $cheesemask $masked_expmap MUL clobber=yes

asmooth inset=$masked_image 
    outset=$output_smooth_image
    weightset=$masked_expmap 
    withweightset=yes
    withexpimageset=yes
    expimageset=$masked_expmap
    desiredsnr=$desiredsnr
"""

with open('run_asmooth.sh', 'w') as file:
    file.write(sh_file_content)

print('asmooth shell script created as run_asmooth.sh\n')
print('Running asmooth...')

subprocess.run(["bash", "run_asmooth.sh"])

print('asmooth completed successfully\n')
print("\n========================================")