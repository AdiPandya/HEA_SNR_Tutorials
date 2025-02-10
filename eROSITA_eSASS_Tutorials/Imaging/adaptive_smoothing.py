import subprocess
import numpy as np
from astropy.io import fits
import os
import argparse
import warnings
import time
import sys
import logging

start_time = time.time()
warnings.filterwarnings("ignore")

# Parse input arguments
parser = argparse.ArgumentParser(description="Adaptivelysmooth Images using xmm-sas or eSASS")
parser.add_argument("input_image", type=str, help="Path to the input image file.")
parser.add_argument("input_expmap", type=str, help="Path to the input exposure map file.")
parser.add_argument("desired_snr", type=int, default=30, help="Desired signal-to-noise ratio for asmooth.")
parser.add_argument("cheesemask_file", type=str, help="Path to the cheese-mask file.")
parser.add_argument("asmooth_tool", type=str, choices=["xmm-sas", "eSASS", "xmmsas", "esass"], help="Tool to use for asmooth: 'xmm-sas/xmmsas' or 'eSASS/esass'")
parser.add_argument("--boxlist_file", type=str, help="Path to the boxlist file. Required if asmooth_tool is eSASS.")
parser.add_argument("--detmask_file", type=str, help="Path to the detection mask file. Required if asmooth_tool is eSASS.")
parser.add_argument("--emin", type=float, help="Minimum energy for asmooth Image in eV. Required if asmooth_tool is eSASS.")
parser.add_argument("--emax", type=float, help="Maximum energy for asmooth Image in eV. Required if asmooth_tool is eSASS.")
parser.add_argument("--ds9", action="store_true", default=False, help="Flag to open the asmooth image in DS9. Default is False.")
args = parser.parse_args()

# Input parameters
input_image = args.input_image
input_expmap = args.input_expmap
input_cheesemask = args.cheesemask_file
desired_snr = args.desired_snr
asmooth_tool = args.asmooth_tool
boxlist_file = args.boxlist_file
detmask_file = args.detmask_file
emin = args.emin
emax = args.emax
open_ds9 = args.ds9

# Set up logging
adapative_smoothing_dir = os.path.dirname(input_image)
log_filename = os.path.join(adapative_smoothing_dir, "adaptive_smoothing.log")

if os.path.exists(log_filename):
    os.remove(log_filename)

logging.basicConfig(filename=log_filename, level=logging.INFO, format='%(message)s')
logger = logging.getLogger()

# Add a stream handler to print to console without timestamps
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console_handler)

# Log the start date and time
start_datetime = time.strftime("%d-%m-%Y %H:%M:%S", time.localtime(start_time))
logger.info(f'Start date and time: {start_datetime}')
logger.info(f'Command used: python {" ".join(sys.argv)}')
logger.info("\n========================================\n")
logger.info('Setting up the environment...\n')

logger.info('Starting the workflow with the following parameters:')
logger.info(f'    Input image: {input_image}')
logger.info(f'    Input exposure map: {input_expmap}')
logger.info(f'    Desired signal-to-noise ratio for asmooth: {desired_snr}')
logger.info(f'    Input Cheese-mask file: {input_cheesemask}')
logger.info(f'    Asmooth tool: {asmooth_tool}')
if asmooth_tool == "eSASS" or asmooth_tool == "esass":
    if boxlist_file is None:
        logger.error('Error: Boxlist file is required for eSASS asmooth.')
        exit()
    else:
        logger.info(f'    Boxlist file: {boxlist_file}')
    if detmask_file is None:
        logger.error('Error: Detection mask file is required for eSASS asmooth.')
        exit()
    else:
        logger.info(f'    Detection mask file: {detmask_file}')

    if emin is None or emax is None:
        logger.error('Error: Minimum and Maximum energy are required for eSASS asmooth.')
        exit()
    else:
        logger.info(f'    Minimum energy: {emin} eV')
        logger.info(f'    Maximum energy: {emax} eV')
logger.info(f'    Open asmooth image in DS9 flag: {open_ds9}')
logger.info(f'    Log file: {log_filename}')

########## Multiplying the point source mask with image and exposure map ##########
logger.info("\n========================================\n")
logger.info('Multiplying the image and exposure map with the cheese-mask...\n')

output_masked_image = input_image.replace(".fits", "_masked.fits")
output_masked_expmap = input_expmap.replace(".fits", "_masked.fits")

# Read the image and exposure map
image_data = fits.getdata(input_image)
expmap_data = fits.getdata(input_expmap)
cheesemask_data = fits.getdata(input_cheesemask)

# Multiply the image and exposure map with the cheese-mask
masked_image_data = image_data * cheesemask_data
masked_expmap_data = expmap_data * cheesemask_data

# Save the masked image and exposure map
masked_image_hdu = fits.PrimaryHDU(masked_image_data)
masked_image_hdu.header.update(fits.getheader(input_image))
masked_image_hdu.writeto(output_masked_image, overwrite=True)

masked_expmap_hdu = fits.PrimaryHDU(masked_expmap_data)
masked_expmap_hdu.header.update(fits.getheader(input_expmap))
masked_expmap_hdu.writeto(output_masked_expmap, overwrite=True)

logger.info('Masked image and exposure map saved as:')
logger.info(f'{output_masked_image}')
logger.info(f'{output_masked_expmap}')
logger.info("\n========================================\n")

########## Running asmooth using xmmsas ##########

if asmooth_tool == "xmm-sas" or asmooth_tool == "xmmsas":
    logger.info('Running asmooth using xmm-sas and filling the holes...\n')

    output_asmooth_image = input_image.replace(".fits", f"_asmooth_xmmsas_snr{float(desired_snr)}.fits")

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

    logger.info('asmooth shell script created as run_asmooth.sh\n')
    logger.info('Running asmooth...\n')

    with open(log_filename, "a") as log_file:
        subprocess.run(["bash", "run_asmooth.sh"], stdout=log_file, stderr=log_file)

    logger.info('\nxmm-sas Asmooth completed successfully!')
    logger.info(f'Adaptive Smoothed image saved as {output_asmooth_image}')

########## Running asmooth using eSASS ##########

elif asmooth_tool == "eSASS" or asmooth_tool == "esass":
    logger.info('Running asmooth using eSASS and filling the holes...\n')

    output_asmooth_image = input_image.replace(".fits", f"_asmooth_esass_snr{float(desired_snr)}.fits")
    
    # Run asmooth
    def erbackmap_asmooth(input_image, input_expmap, input_boxlist, input_detmask, output_image, emin, emax, snr, log_file=None):
        subprocess.run(["erbackmap",
                        f"image={input_image}",
                        f"expimage={input_expmap}",
                        f"boxlist={input_boxlist}",
                        f"detmask={input_detmask}",
                        f"bkgimage={output_image}",
                        f"emin={emin}",
                        f"emax={emax}",
                        f"snr={snr}",
                        "scut=0.001",
                        "mlmin=1.0E6",
                        "maxcut=0.5",
                        "smoothval=4.0",
                        "clobber=yes",
                        ], stdout=log_file, stderr=log_file)
    
    with open(log_filename, "a") as log_file:
        # Use input_expmap instead of output_masked_expmap because it creates holes in the output result
        erbackmap_asmooth(output_masked_image, input_expmap, boxlist_file, detmask_file, output_asmooth_image, emin, emax, desired_snr, log_file=log_file)

    with open(log_filename, "r") as log_file:
        log_content = log_file.readlines()
        erbackmap_count = sum(1 for line in log_content if 'erbackmap: DONE' in line)
        if erbackmap_count == 1:
            logger.info(f'\n Adapative smoothed image saved as {output_asmooth_image}')
        else:
            logger.info(f'Error: Adaptive smoothing failed!')
            exit()
    
    ########## Dividing the asmooth image with the exposure map ##########
    logger.info('\n Dividing the asmooth image with the exposure map to create the exposure-corrected asmooth image...\n')

    output_expcorr_image = output_asmooth_image.replace(".fits", "_expcorr.fits")

    asmooth_image_data = fits.getdata(output_asmooth_image)
    detmask_file_data = fits.getdata(detmask_file)
    expcorr_image_data = asmooth_image_data / (expmap_data*detmask_file_data)

    expcorr_image_hdu = fits.PrimaryHDU(expcorr_image_data)
    expcorr_image_hdu.header.update(fits.getheader(output_asmooth_image))
    expcorr_image_hdu.writeto(output_expcorr_image, overwrite=True)

    logger.info(f' Exposure-corrected asmooth image saved as {output_expcorr_image}')

end_time = time.time()
time_taken = end_time - start_time

logger.info("\n========================================\n")
if time_taken < 600:
    logger.info(f'** All tasks completed successfully in {time_taken:.2f} seconds **')
if time_taken >= 600:
    logger.info(f'** All tasks completed successfully in {(time_taken/60):.2f} minutes **')
if time_taken >= 3600:
    logger.info(f'** All tasks completed successfully in {(time_taken/3600):.2f} hours **')
logger.info("\n========================================\n")

if open_ds9:
    try:
        if asmooth_tool == "xmm-sas" or asmooth_tool == "xmmsas":
            view_image = output_asmooth_image
        elif asmooth_tool == "eSASS" or asmooth_tool == "esass":
            view_image = output_expcorr_image

        ds9_command = f"ds9 {view_image} &"
        subprocess.run(ds9_command, shell=True, check=True)
        logger.info(f'DS9 opened with asmooth image file {view_image}.')

    except Exception as e:
        logger.error(f'Failed to open DS9 with error: {e}')