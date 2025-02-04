import subprocess
import numpy as np
from astropy.io import fits
import os
from tqdm import tqdm
import argparse
from concurrent.futures import ProcessPoolExecutor
import warnings
import time
start_time = time.time()
warnings.filterwarnings("ignore")

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

# Main function to parse arguments and execute the workflow
if __name__ == "__main__":
    print("\n========================================\n")
    print('Setting up the environment...')
    parser = argparse.ArgumentParser(description='Creating Images, Exposure Maps and Exposure Corrected Images')
    parser.add_argument('event_file', type=str, help='Input event file name')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('band_min', type=float, help='Minimum energy band in eV', nargs='?')
    parser.add_argument('band_max', type=float, help='Maximum energy band in eV', nargs='?')
    parser.add_argument('--rgb', action='store_true', help='Create RGB images')
    parser.add_argument('--rgb_bands', type=float, nargs=6, metavar=('R_MIN', 'R_MAX', 'G_MIN', 'G_MAX', 'B_MIN', 'B_MAX'), help='Energy bands for RGB image creation in eV')

    args = parser.parse_args()

    event_file = args.event_file
    output_dir = args.output_dir
    band_min = args.band_min
    band_max = args.band_max
    create_rgb = args.rgb
    rgb_bands = args.rgb_bands

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    print('\nEnvironment setup complete\n')

    print('Starting the image creation with the following parameters:')
    print(f'    Event file: {event_file}')
    print(f'    Output directory: {output_dir}')
    print(f'    Create RGB: {create_rgb}')

    if create_rgb:
        if rgb_bands is None:
            bands = [(200,700), (700, 1100), (1100, 2300)]
        else:
            bands = [(rgb_bands[0], rgb_bands[1]), (rgb_bands[2], rgb_bands[3]), (rgb_bands[4], rgb_bands[5])]
        
        print("\n========================================\n")
        print("Creating RGB image with the following bands:")
        for band in bands:
            print(f" * {band[0]} - {band[1]} eV")
        
        log_file_path = f'{output_dir}/merged_image_RGB.log'
        with open(log_file_path, 'w+') as log_file:
            pass

        def process_band(band):
            output_image = f'{output_dir}/merged_image_{band[0]}_{band[1]}.fits'
            output_expmap = f'{output_dir}/merged_expmap_{band[0]}_{band[1]}.fits'
            output_exp_corr = f'{output_dir}/merged_exp_corr_{band[0]}_{band[1]}.fits'
        
            with open(log_file_path, 'a') as log_file:
                run_evtool(event_file, output_image, band[0]/1000, band[1]/1000, log_file=log_file)
                run_expmap(event_file, output_image, output_expmap, band[0]/1000, band[1]/1000, log_file=log_file)
                exp_corr(output_image, output_expmap, output_exp_corr)

        with ProcessPoolExecutor() as executor:
            executor.map(process_band, bands)

        with open(log_file_path, 'r') as log_file:
            log_content = log_file.readlines()
            evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
            expmap_count = sum(1 for line in log_content if 'expmap: DONE' in line)
            if evtool_count == 3 and expmap_count == 3:
                print("\n========================================")
                print("Exposure corrected RGB Image creation successful!")
        
    else:
        if band_min is None or band_max is None:
            parser.error("band_min and band_max are required when --rgb is not set.")
        
        with open(f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.log', 'w+') as log_file:
            
            print("\n========================================\n")
            print(f"Creating image with band: {band_min} - {band_max} eV.")
            # Run evtool for the specified band
            run_evtool(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                       band_min/1000, band_max/1000, log_file=log_file)
            
            print("\n========================================\n")
            print(f"Creating exposure map for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
            # Run expmap for the specified band
            run_expmap(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                       f'{output_dir}/merged_expmap_{int(band_min)}_{int(band_max)}.fits', 
                       band_min/1000, band_max/1000, log_file=log_file)
            
            print("\n========================================\n")
            print(f"Creating exposure corrected image for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
            # Perform exposure correction for the specified band
            exp_corr(f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits',
                     f'{output_dir}/merged_expmap_{int(band_min)}_{int(band_max)}.fits',
                     f'{output_dir}/merged_exp_corr_{int(band_min)}_{int(band_max)}.fits')

            log_file.seek(0)
            log_content = log_file.readlines()
            evtool_success = any('evtool: DONE' in line for line in log_content)
            expmap_success = any('expmap: DONE' in line for line in log_content)
            if evtool_success and expmap_success:
                print("\n========================================")
                print("Exposure corrected Image creation successful!")
    
    end_time = time.time()
    time_taken = end_time - start_time
    if time_taken < 600:
        print(f"** Task completed in {time_taken:.2f} seconds **")
    if time_taken >= 600:
        print(f"** Task completed in {time_taken/60:.2f} minutes **")
    if time_taken >= 3600:
        print(f"** Task completed in {time_taken/3600:.2f} hours **")
    print("========================================\n")
