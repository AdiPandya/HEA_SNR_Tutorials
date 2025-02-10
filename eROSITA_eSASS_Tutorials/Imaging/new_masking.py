import numpy as np
from astropy.io import fits
from astropy import wcs
from tqdm import tqdm
import argparse
import logging
import os
import time
import subprocess
import sys

start_time = time.time()

# Parse input arguments
parser = argparse.ArgumentParser(description="Create a new cheese-mask using the region file.")
parser.add_argument("cheesemask_file", type=str, help="Path to the cheese-mask file.")
parser.add_argument("cheesemask_regions", type=str, help="Path to the cheese-mask region file.")
parser.add_argument("detmask_file", type=str, help="Path to the detection mask file.")
parser.add_argument("--new_cheesemask", action="store_true", default=False, help="Flag to create a new cheese-mask or overwrite the existing one. Default is False.")
parser.add_argument("--ds9", action="store_true", default=False, help="Flag to open the cheese-mask and region file in DS9. Default is False.")

args = parser.parse_args()

# Input parameters
input_cheesemask = args.cheesemask_file
input_cheesemask_regions = args.cheesemask_regions
detmask_file = args.detmask_file
new_cheesemask = args.new_cheesemask
open_ds9 = args.ds9

# Set up logging
cheesemask_dir = os.path.dirname(input_cheesemask)
log_filename = os.path.join(cheesemask_dir, 'new_cheesemask.log')

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
logger.info(f'    Input Cheese-mask file: {input_cheesemask}')
logger.info(f'    Cheese-mask regions file: {input_cheesemask_regions}')
logger.info(f'    Detection mask file: {detmask_file}')
if new_cheesemask:
    logger.info(f'    new_cheesemask flag: {new_cheesemask}. Thus creating a new cheese-mask file')
else:
    logger.info(f'    new_cheesemask flag: {new_cheesemask}. Thus overwriting the existing cheese-mask.')
logger.info(f'    Open cheesemask file and region in DS9 flag: {open_ds9}')
logger.info(f'    Log file: {log_filename}')

########## Creating a new cheese-mask ##########
logger.info("\n========================================\n")
logger.info(f'Making the cheese-mask with the new regions from {input_cheesemask_regions}...')

if new_cheesemask:
    cheesemask_file = input_cheesemask.replace(".fits", "_new.fits")
    logger.info(f'New cheese-mask will be saved as {cheesemask_file}\n')
else:
    cheesemask_file = input_cheesemask
    logger.info(f'Existing cheese-mask will be overwritten.\n')

hdulist = fits.open(input_cheesemask)
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

logger.info(f'\nCheese-mask saved as {cheesemask_file}.')

end_time = time.time()
time_taken = end_time - start_time

logger.info("\n========================================\n")
if time_taken < 600:
    logger.info(f'** Cheese-mask creation completed successfully in {time_taken:.2f} seconds **')
if time_taken >= 600:
    logger.info(f'** Cheese-mask creation completed successfully in {(time_taken/60):.2f} minutes **')
if time_taken >= 3600:
    logger.info(f'** Cheese-mask creation completed successfully in {(time_taken/3600)::.2f} hours **')
logger.info("\n========================================\n")

# Open the cheese-mask in DS9 with the region file
if open_ds9:
    try:
        ds9_command = f"ds9 {cheesemask_file} -regions {input_cheesemask_regions} &"
        subprocess.run(ds9_command, shell=True, check=True)
        logger.info(f'DS9 opened with cheese-mask file {cheesemask_file} and regions {input_cheesemask_regions}.')
    except Exception as e:
        logger.error(f'Failed to open DS9 with error: {e}')