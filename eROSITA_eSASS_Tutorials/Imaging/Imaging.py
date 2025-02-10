import subprocess
import numpy as np
from astropy.io import fits
import os
from tqdm import tqdm
import argparse
from concurrent.futures import ProcessPoolExecutor
import warnings
import time
import logging
import sys
start_time = time.time()
warnings.filterwarnings("ignore")

# Parse input arguments
parser = argparse.ArgumentParser(description='Creating Images, Exposure Maps and Exposure Corrected Images')
parser.add_argument('event_file', type=str, help='Input event file name')
parser.add_argument('output_dir', type=str, help='Output directory')
parser.add_argument('band_min', type=float, help='Minimum energy band in eV', nargs='?')
parser.add_argument('band_max', type=float, help='Maximum energy band in eV', nargs='?')
parser.add_argument('--rgb', action='store_true', help='Create RGB images')
parser.add_argument('--rgb_bands', type=float, nargs=6, metavar=('R_MIN', 'R_MAX', 'G_MIN', 'G_MAX', 'B_MIN', 'B_MAX'), help='Energy bands for RGB image creation in eV')
parser.add_argument('--ds9', action='store_true', help='Flag to open the images in DS9')
args = parser.parse_args()

# Input parameters
event_file = args.event_file
output_dir = args.output_dir
band_min = args.band_min
band_max = args.band_max
create_rgb = args.rgb
rgb_bands = args.rgb_bands
open_ds9 = args.ds9

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

if create_rgb:
    log_filename = f'{output_dir}/merged_image_RGB.log'
else:
    log_filename = f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.log'
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
logger.info('Setting up the environment...')
logger.info('Starting the workflow with the following parameters:')
logger.info(f'    Input event file: {event_file}')
logger.info(f'    Output directory: {output_dir}')
logger.info(f'    Create RGB: {create_rgb}')
if not create_rgb:
    if band_min is None or band_max is None:
        logger.error("band_min and band_max are required when --rgb is not set.")
        exit()
    else:
        logger.info(f'    Energy band: {band_min} - {band_max} eV')
if create_rgb:
    if rgb_bands is None:
        bands = [(200,700), (700, 1100), (1100, 2300)]
        logger.info("    Default RGB bands: (200-700), (700-1100), (1100-2300) eV")
    else:
        bands = [(rgb_bands[0], rgb_bands[1]), (rgb_bands[2], rgb_bands[3]), (rgb_bands[4], rgb_bands[5])]
        logger.info(f'    User RGB bands: {bands} eV')
logger.info(f'    Open images in DS9: {open_ds9}')
logger.info(f'    Log file: {log_filename}')

########## Functions ##########

# Function to run the evtool command
def run_evtool(input_name, output_name, emin, emax, gti_type='FLAREGTI', flag_type='0xe00fff30', size='auto', pattern='15', telid='1 2 3 4 5 6 7', log_file=None):
    subprocess.run(['evtool', 
                    f'eventfiles={input_name}', 
                    f'outfile={output_name}', 
                    f'gti={gti_type}', 
                    f'flag={flag_type}', 
                    f'pattern={pattern}', 
                    f'emin={emin}', 
                    f'emax={emax}',
                    f'image=yes',
                    f'events=yes',
                    f'telid={telid}',
                    f'size={size}'
                    ],
                    stdout=log_file,
                    stderr=log_file)

# Function to run the expmap command
def run_expmap(input_eventlist, input_image, output_name, emin, emax, log_file=None):
    subprocess.run(['expmap', 
                    f'inputdatasets={input_eventlist}', 
                    f'templateimage={input_image}', 
                    f'mergedmaps={output_name}', 
                    f'emin={emin}', 
                    f'emax={emax}',
                    'withvignetting=yes',
                    'withweights=yes',
                    ],
                    stdout=log_file,
                    stderr=log_file)

# Function to perform exposure correction
def exp_corr(input_image, input_expmap, output_name):
    cts = fits.open(input_image)[0].data
    exp = fits.open(input_expmap)[0].data
    hdr = fits.getheader(input_image)

    exp_corr = cts/exp

    fits.writeto(output_name, exp_corr, header=hdr, overwrite=True)

########## Image Creation ##########

if create_rgb:
    logger.info("\n========================================\n")
    logger.info("Creating RGB image with the following bands:\n")
    for band in bands:
        logger.info(f" * {band[0]} - {band[1]} eV")

    def process_band(band):
        output_image = f'{output_dir}/merged_image_{band[0]}_{band[1]}.fits'
        output_expmap = f'{output_dir}/merged_expmap_{band[0]}_{band[1]}.fits'
        output_exp_corr = f'{output_dir}/merged_exp_corr_{band[0]}_{band[1]}.fits'
    
        with open(log_filename, 'a') as log_file:
            run_evtool(event_file, output_image, band[0]/1000, band[1]/1000, log_file=log_file)
            run_expmap(event_file, output_image, output_expmap, band[0]/1000, band[1]/1000, log_file=log_file)
            exp_corr(output_image, output_expmap, output_exp_corr)

    with ProcessPoolExecutor() as executor:
        executor.map(process_band, bands)

    with open(log_filename, 'r') as log_file:
        log_content = log_file.readlines()
        evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
        expmap_count = sum(1 for line in log_content if 'expmap: DONE' in line)
        if evtool_count == 3 and expmap_count == 3:
            logger.info("\n========================================")
            logger.info("Exposure corrected RGB Image creation successful!")
    
else:
    with open(log_filename, 'a') as log_file:
        
        logger.info("\n========================================\n")
        logger.info(f"Creating image with band: {band_min} - {band_max} eV.\n")
        # Run evtool for the specified band
        run_evtool(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                    band_min/1000, band_max/1000, log_file=log_file)
        
        logger.info("\n========================================\n")
        logger.info(f"Creating exposure map for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
        # Run expmap for the specified band
        run_expmap(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                    f'{output_dir}/merged_expmap_{int(band_min)}_{int(band_max)}.fits', 
                    band_min/1000, band_max/1000, log_file=log_file)
        
        logger.info("\n========================================\n")
        logger.info(f"Creating exposure corrected image for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
        # Perform exposure correction for the specified band
        exp_corr(f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits',
                    f'{output_dir}/merged_expmap_{int(band_min)}_{int(band_max)}.fits',
                    f'{output_dir}/merged_exp_corr_{int(band_min)}_{int(band_max)}.fits')
        
    with open(log_filename, 'r') as log_file:
        log_file.seek(0)
        log_content = log_file.readlines()
        evtool_success = any('evtool: DONE' in line for line in log_content)
        expmap_success = any('expmap: DONE' in line for line in log_content)
        if evtool_success and expmap_success:
            logger.info("\n========================================")
            logger.info("Exposure corrected Image creation successful!")

end_time = time.time()
time_taken = end_time - start_time
if time_taken < 600:
    logger.info(f"** Task completed in {time_taken:.2f} seconds **")
if time_taken >= 600:
    logger.info(f"** Task completed in {time_taken/60:.2f} minutes **")
if time_taken >= 3600:
    logger.info(f"** Task completed in {time_taken/3600:.2f} hours **")
logger.info("========================================\n")

# Open the images in DS9
if open_ds9:
    try:
        if create_rgb:
            ds9_command = f"ds9 -rgb -red {output_dir}/merged_exp_corr_{bands[0][0]}_{bands[0][1]}.fits -green {output_dir}/merged_exp_corr_{bands[1][0]}_{bands[1][1]}.fits -blue {output_dir}/merged_exp_corr_{bands[2][0]}_{bands[2][1]}.fits &"
        else:
            ds9_command = f"ds9 {output_dir}/merged_exp_corr_{int(band_min)}_{int(band_max)}.fits &"
        subprocess.run(ds9_command, shell=True, check=True)
        logger.info(f'DS9 opened with images.')
    except Exception as e:
        logger.error(f'Failed to open DS9 with error: {e}')