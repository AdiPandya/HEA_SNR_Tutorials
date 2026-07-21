import glob
import subprocess
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from concurrent.futures import ProcessPoolExecutor
import os, shutil
from tqdm import tqdm
import cv2
import warnings
import argparse
import time
import logging
import sys
start_time = time.time()
warnings.filterwarnings("ignore")

# Parse input arguments
parser = argparse.ArgumentParser(description='Extract eROSITA spectrum.')
parser.add_argument('input_events', type=str, help='Input event files from which to extract spectra')
# parser.add_argument('output_dir', type=str, help='Output directory to save the results')
parser.add_argument('src_region', type=str, help='DS9 region file for source region in pixel coordinates.')
parser.add_argument('bkg_region', type=str, help='DS9 region file for background region in pixel coordinates.')
parser.add_argument("cheesemask_file", type=str, help="Path to the cheese-mask file.")
parser.add_argument("detmask_file", type=str, help="Path to the detection mask file.")
parser.add_argument("--ds9", action="store_true", default=False, help="Flag to open the mask fits file created from the DS9 region files. Default is False.")
args = parser.parse_args()

event_file = args.input_events
src_region = args.src_region
bkg_region = args.bkg_region
cheesemask_file = args.cheesemask_file
detmask_file = args.detmask_file
open_ds9 = args.ds9

output_dir = os.path.dirname(src_region)

# Set up logging
log_filename = os.path.join(output_dir, "spectra_extraction.log")
    
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

import xspec as xs

xs.Xset.chatter = 1
xs.Xset.xsect = 'vern'
xs.Xset.abund = 'wilm'
xs.Fit.statMethod = 'cstat'
xs.Xset.cosmo = "70 0 0.7" 

logger.info('Starting the workflow with the following parameters:')
logger.info(f'    Event File: {event_file}')
logger.info(f'    Output directory: {output_dir}')
logger.info(f'    Source region: {src_region}')
logger.info(f'    Background region: {bkg_region}')
logger.info(f'    Cheese mask file: {cheesemask_file}')
logger.info(f'    Detection mask file: {detmask_file}')

if not src_region.endswith('.reg') or not bkg_region.endswith('.reg'):
    logger.error("Source and Background region file must have a '.reg' extension.")
    exit()

logger.info(f'    Open output cheese-mask and region in DS9 flag: {open_ds9}')
logger.info(f'    Log file: {log_filename}')

########## Defining the functions ##########

# Function to calculate centroid of a polygon if the input region is a polygon
def calculate_polygon_centroid(ra, dec, wcs):
    # Shoelace formula for the area of the polygon in pixel coordinates
    area_pixel = 0.5 * np.abs(np.dot(ra, np.roll(dec, 1)) - np.dot(dec, np.roll(ra, 1)))
    centroid_pixel_x = np.sum((ra + np.roll(ra, 1)) * (ra * np.roll(dec, 1) - np.roll(ra, 1) * dec)) / (6 * area_pixel)
    centroid_pixel_y = np.sum((dec + np.roll(dec, 1)) * (ra * np.roll(dec, 1) - np.roll(ra, 1) * dec)) / (6 * area_pixel)
    
    centroid_world = wcs.pixel_to_world(centroid_pixel_x, centroid_pixel_y)
    return np.array([centroid_world.ra.deg, centroid_world.dec.deg])

# Function to create a new mask from the region file and cheese mask
def create_new_mask(reg, cheese, output_file):
    reg = open(reg).readlines()
    cheese_hdu = fits.open(cheese)
    cheese_data = cheese_hdu[0].data
    wcs = WCS(cheese_hdu[0].header)
    centroid = np.zeros(2)

    if reg[0].startswith('polygon'):
        shape = 'polygon'
        reg = reg[0].replace('polygon(', '').replace(')\n', '').split(',')
        ra = [float(reg[i]) for i in range(0, len(reg), 2)]
        dec = [float(reg[i]) for i in range(1, len(reg), 2)]
        ra = np.r_[ra, ra[0]]
        dec = np.r_[dec, dec[0]]
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.fillPoly(no_mask, [np.array([ra, dec], np.int32).T], 1)

        centroid = calculate_polygon_centroid(ra, dec, wcs)
    
    elif reg[0].startswith('box'):
        shape = 'box'
        reg = reg[0].replace('box(', '').replace(')\n', '').split(',')
        ra = float(reg[0])
        dec = float(reg[1])
        width = float(reg[2])
        height = float(reg[3])
        angle = float(reg[4])
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.boxPoints(((ra, dec), (width, height), angle))
        shape_mask = cv2.fillPoly(no_mask, [shape_mask.astype(np.int32)], 1)

        centroid_coords = wcs.pixel_to_world(ra, dec)
        centroid = np.array([centroid_coords.ra.deg, centroid_coords.dec.deg])
    
    elif reg[0].startswith('circle'):
        shape = 'circle'
        reg = reg[0].replace('circle(', '').replace(')\n', '').split(',')
        ra = float(reg[0])
        dec = float(reg[1])
        radius = float(reg[2])
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.circle(no_mask, (int(ra), int(dec)), int(radius), 1, thickness=-1)
        
        icrs_coords = wcs.pixel_to_world(ra, dec)
        centroid[0] = icrs_coords.ra.deg
        centroid[1] = icrs_coords.dec.deg
        
    elif reg[0].startswith('annulus'):
        shape = 'annulus'
        reg = reg[0].replace('annulus(', '').replace(')\n', '').split(',')
        ra = float(reg[0])
        dec = float(reg[1])
        inner_radius = float(reg[2])
        outer_radius = float(reg[3])
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.circle(no_mask.copy(), (int(ra), int(dec)), int(outer_radius), 1, thickness=-1)
        shape_mask -= cv2.circle(no_mask.copy(), (int(ra), int(dec)), int(inner_radius), 1, thickness=-1)
        icrs_coords = wcs.pixel_to_world(ra, dec)
        centroid[0] = icrs_coords.ra.deg
        centroid[1] = icrs_coords.dec.deg
        
    elif reg[0].startswith('ellipse') and all('!ellipse' not in line for line in reg):
        shape = 'ellipse'
        reg = reg[0].replace('ellipse(', '').replace(')\n', '').split(',')
        ra = float(reg[0])
        dec = float(reg[1])
        a = float(reg[2])
        b = float(reg[3])
        theta = float(reg[4])
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.ellipse(no_mask, (int(ra), int(dec)), (int(a), int(b)), int(theta), 0, 360, 1, thickness=-1)
        icrs_coords = wcs.pixel_to_world(ra, dec)
        centroid[0] = icrs_coords.ra.deg
        centroid[1] = icrs_coords.dec.deg

    elif any('!ellipse' in line for line in reg):
        shape = 'ellipse annulus'
        ellipse_list = []
        for line in reg:
            ellipses = [subpart.strip() for part in line.split('&') for subpart in part.split('!') if 'ellipse(' in subpart]
            # print(len(ellipses))
            for ellipse in ellipses:
                elements = list(map(int, map(float, ellipse.replace('ellipse(', '').replace(')', '').split(','))))
                if len(elements) == 5 and elements not in ellipse_list:
                    ellipse_list.append(elements)
        ellipse_list = np.array(ellipse_list)
        outer_ellipse = ellipse_list[np.argmax(ellipse_list[:, 2:4].max(axis=1))]
        inner_ellipse = ellipse_list[np.argmin(ellipse_list[:, 2:4].max(axis=1))]
        
        ra = float(outer_ellipse[0])
        dec = float(outer_ellipse[1])
        a_outer = float(outer_ellipse[2])
        b_outer = float(outer_ellipse[3])
        theta_outer = float(outer_ellipse[4])
        
        a_inner = float(inner_ellipse[2])
        b_inner = float(inner_ellipse[3])
        theta_inner = float(inner_ellipse[4])
        
        no_mask = np.zeros_like(cheese_data)
        shape_mask = cv2.ellipse(no_mask.copy(), (int(ra), int(dec)), (int(a_outer), int(b_outer)), int(theta_outer), 0, 360, 1, thickness=-1)
        shape_mask -= cv2.ellipse(no_mask.copy(), (int(ra), int(dec)), (int(a_inner), int(b_inner)), int(theta_inner), 0, 360, 1, thickness=-1)
        icrs_coords = wcs.pixel_to_world(ra, dec)
        centroid[0] = icrs_coords.ra.deg
        centroid[1] = icrs_coords.dec.deg
        
    elif reg[0].startswith('panda'):
        shape = 'panda'
        reg = reg[0].replace('panda(', '').replace(')\n', '').split(',')
        ra = float(reg[0])
        dec = float(reg[1])
        start_angle = float(reg[2])
        finish_angle = float(reg[3])
        # reg[4] is ignored
        inner_radius = float(reg[5])
        outer_radius = float(reg[6])
        # reg[7] is ignored

        no_mask = np.zeros_like(cheese_data)

        # Generate points for the outer arc
        angles = np.linspace(start_angle, finish_angle, num=200)
        outer_x = ra + outer_radius * np.cos(np.deg2rad(angles))
        outer_y = dec + outer_radius * np.sin(np.deg2rad(angles))

        # Generate points for the inner arc (reverse order)
        inner_x = ra + inner_radius * np.cos(np.deg2rad(angles[::-1]))
        inner_y = dec + inner_radius * np.sin(np.deg2rad(angles[::-1]))

        # Combine into a polygon
        pts = np.vstack([
            np.stack([outer_x, outer_y], axis=1),
            np.stack([inner_x, inner_y], axis=1)
        ]).astype(np.int32)

        # Draw filled polygon
        shape_mask = cv2.fillPoly(no_mask.copy(), [pts], 1)

        icrs_coords = wcs.pixel_to_world(ra, dec)
        centroid[0] = icrs_coords.ra.deg
        centroid[1] = icrs_coords.dec.deg
        
    else:
        logger.info("Unsupported region format in the region file.")
        exit()

    new_mask = cheese_data * shape_mask

    # Save the new_mask to a file
    hdu = fits.PrimaryHDU(new_mask)
    hdu.header.update(cheese_hdu[0].header)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output_file, overwrite=True)
    
    logger.info(f'Region mask of shape {shape} created in: {output_file}')

    return centroid.astype(float)

# Function to run srctool
def run_srctool(event_file, centroid, out_prefix, reg_mask, TM = '8', log_file=None):
    subprocess.run(['srctool',
                f'eventfiles={event_file}',
                f'srccoord=icrs; {centroid[0]},{centroid[1]}',
                f'prefix={out_prefix}',
                'todo=SPEC ARF RMF',
                'insts=1 2 3 4 6',
                f'writeinsts= {TM}',
                f'srcreg=mask {reg_mask}', 
                f'backreg=NONE', 
                'exttype=BETA', 
                'extpars=11000 2', 
                'tstep=1', 
                'xgrid=10', 
                'psftype=NONE', 
                'pat_sel=15',
                'clobber=yes'], stdout=log_file, stderr=log_file)

# Function to perform exposure correction on source and background spectrum
def spectra_corr(src_spec, bkg_spec, src_arf=None, bkg_arf=None):
    exposure_src = fits.getval(src_spec, "exposure", ext=1)
    exposure_bkg = fits.getval(bkg_spec, "exposure", ext=1)
    backscal_src = fits.getval(src_spec, "backscal", ext=1)
    backscal_bkg = fits.getval(bkg_spec, "backscal", ext=1)
    regarea_src = fits.getval(src_spec, "regarea", ext=1)
    regarea_bkg = fits.getval(bkg_spec, "regarea", ext=1)

    r_src = regarea_src / backscal_src
    r_bkg = regarea_bkg / backscal_bkg

    if src_arf is None or bkg_arf is None:
        src_arf = src_spec.replace("SourceSpec", "ARF")
        bkg_arf = bkg_spec.replace("SourceSpec", "ARF")
        if not os.path.exists(src_spec) or not os.path.exists(bkg_spec):
            logger.info('Error: Source or background spectrum file does not exist. src_arf and bkg_arf are required.')
            exit()

    # Copy files and modify values
    src_spec_corr = src_spec.replace(".fits", "_corr.fits")
    bkg_spec_corr = bkg_spec.replace(".fits", "_corr.fits")
    src_arf_corr = src_arf.replace(".fits", "_corr.fits")
    bkg_arf_corr = bkg_arf.replace(".fits", "_corr.fits")

    short_src_arf = fits.getval(src_spec, "ANCRFILE", ext=1)
    short_bkg_arf = fits.getval(bkg_spec, "ANCRFILE", ext=1)
    short_src_arf_corr = short_src_arf.replace(".fits", "_corr.fits")
    short_bkg_arf_corr = short_bkg_arf.replace(".fits", "_corr.fits")

    shutil.copy(src_spec, src_spec_corr)
    fits.setval(src_spec_corr, "exposure", value=exposure_src / r_src, ext=1)
    fits.setval(src_spec_corr, "ANCRFILE", value=short_src_arf_corr, ext=1)
    shutil.copy(bkg_spec, bkg_spec_corr)
    fits.setval(bkg_spec_corr, "exposure", value=exposure_bkg / r_bkg * (regarea_bkg / regarea_src), ext=1)
    fits.setval(bkg_spec_corr, "ANCRFILE", value=short_bkg_arf_corr, ext=1)

    shutil.copy(src_arf, src_arf_corr)
    with fits.open(src_arf_corr) as hdu:
        hdu[1].data["SPECRESP"] *= r_src
        hdu.writeto(src_arf_corr, overwrite=True)

    shutil.copy(bkg_arf, bkg_arf_corr)
    with fits.open(bkg_arf_corr) as hdu:
        hdu[1].data["SPECRESP"] *= r_bkg
        hdu.writeto(bkg_arf_corr, overwrite=True)
    return src_spec_corr, bkg_spec_corr

# Function to group specta
def run_grouppha(infile, outfile, resp_file, group_type='snmin', group_scale=20.0):
    subprocess.run(["ftgrouppha",
                    f"infile={infile}",
                    f"outfile={outfile}",
                    f"grouptype={group_type}",
                    f"groupscale={group_scale}",
                    f"respfile={resp_file}"])

########## Workflow ##########
logger.info("\n========================================\n")
logger.info('1) Creating masks from the region files...')

src_mask = src_region.replace('.reg', '.fits')
bkg_mask = bkg_region.replace('.reg', '.fits')

src_centroid = create_new_mask(src_region, cheesemask_file, src_mask)
bkg_centroid = create_new_mask(bkg_region, cheesemask_file, bkg_mask)

logger.info("\n========================================\n")
logger.info('2) Extracting the Source spectrum...')
src_spec = False  
src_prefix = src_region.replace('.reg', '_')
with open(log_filename, "a") as log_file:
    run_srctool(event_file, src_centroid, src_prefix, src_mask, log_file=log_file)

with open(log_filename, "r") as log_file:
    log_content = log_file.readlines()
    srctool_count = sum(1 for line in log_content if 'srctool: DONE' in line)
    if srctool_count == 1:
        logger.info(f'\nSpectrum extraction for source completed successfully and files saved with prefix {src_prefix}')
        src_spec = True
    else:
        logger.info(f'\nError: Spectrum extraction for source failed!')
        exit()
      
logger.info("\n========================================\n")
logger.info('3) Extracting the Background spectrum...')

bkg_prefix = bkg_region.replace('.reg', '_')
with open(log_filename, "a") as log_file:
    run_srctool(event_file, bkg_centroid, bkg_prefix, bkg_mask, log_file=log_file)

with open(log_filename, "r") as log_file:
    log_content = log_file.readlines()
    srctool_count = sum(1 for line in log_content if 'srctool: DONE' in line)
    if src_spec == True:
        if srctool_count == 2:
            logger.info(f'\nSpectrum extraction for background completed successfully and files saved with prefix {bkg_prefix}')
        elif srctool_count == 1:
            logger.info(f'\nError: Spectrum extraction for background failed!')
            exit()
    elif src_spec == False:
        if srctool_count == 1:
            logger.info(f'\nSpectrum extraction for background completed successfully and files saved with prefix {bkg_prefix}')
        elif srctool_count == 0:
            logger.info(f'\nError: Spectrum extraction for background failed!')
            exit()

logger.info("\n========================================\n")
logger.info('4) Exposure correction between the source and background...')

src_spectrum = f'{src_prefix}820_SourceSpec_00001.fits'
bkg_spectrum  = f'{bkg_prefix}820_SourceSpec_00001.fits'

if not os.path.exists(src_spectrum) or not os.path.exists(bkg_spectrum):
    logger.error(f'Error: Source spectrum file {src_spectrum} or Background spectrum file {bkg_spectrum} does not exist.')
    exit()

src_spec_corr, bkg_spec_corr = spectra_corr(src_spectrum, bkg_spectrum)

logger.info("\n========================================\n")
logger.info('5) Grouping the spectrum...')

src_grp = src_spec_corr.replace(".fits", ".grp")
bkg_grp = bkg_spec_corr.replace(".fits", ".grp")
src_rmf = src_spectrum.replace("SourceSpec", "RMF")
bkg_rmf = bkg_spectrum.replace("SourceSpec", "RMF")

run_grouppha(src_spec_corr, src_grp, src_rmf)
run_grouppha(bkg_spec_corr, bkg_grp, bkg_rmf)
logger.info(f'Grouped source spectrum saved as: {src_grp}')
logger.info(f'Grouped background spectrum saved as: {bkg_grp}')

end_time = time.time()
time_taken = end_time - start_time

logger.info("\n========================================\n")
if time_taken < 600:
    logger.info(f'** All tasks completed successfully in {time_taken:.2f} seconds **')
if 600 <= time_taken < 3600:
    logger.info(f'** All tasks completed successfully in {(time_taken/60):.2f} minutes **')
if time_taken >= 3600:
    logger.info(f'** All tasks completed successfully in {(time_taken/3600):.2f} hours **')
logger.info("\n========================================\n")

# if open_ds9:
#     try:
#         ds9_command = f"ds9 {cheesemask_file} -regions {output_region_file} &"
#         subprocess.run(ds9_command, shell=True, check=True)
#         logger.info(f'DS9 opened with cheese-mask file {cheesemask_file} and regions {output_region_file}.')
#     except Exception as e:
#         logger.error(f'Failed to open DS9 with error: {e}')