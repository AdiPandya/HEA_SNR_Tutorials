import subprocess
import numpy as np
from astropy.io import fits
import os
from tqdm import tqdm
import argparse
import warnings

# Ignore warnings
warnings.filterwarnings("ignore")

# Function to run the evtool command
def run_evtool(input_name, output_name, emin, emax, gti_type='FLAREGTI', flag_type='0xe00fff30', pattern='15', telid='1 2 3 4 5 6 7', log_file=None):
    subprocess.run(['evtool', 
                    f'eventfiles={input_name}', 
                    f'outfile={output_name}', 
                    f'gti={gti_type}', 
                    f'flag={flag_type}', 
                    f'pattern={pattern}', 
                    f'emin={emin}', 
                    f'emax={emax}',
                    f'image=yes',
                    f'events=no',
                    f'telid={telid}'
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
    parser = argparse.ArgumentParser(description='Creating Images, Exposure Maps and Exposure Corrected Images')
    parser.add_argument('event_file', type=str, help='Input event file name')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('band_min', type=float, help='Minimum energy band', nargs='?')
    parser.add_argument('band_max', type=float, help='Maximum energy band', nargs='?')
    parser.add_argument('--rgb', action='store_true', help='Create RGB images')
    parser.add_argument('--rgb_bands', type=float, nargs=6, metavar=('R_MIN', 'R_MAX', 'G_MIN', 'G_MAX', 'B_MIN', 'B_MAX'), help='Energy bands for RGB image creation')

    args = parser.parse_args()

    event_file = args.event_file
    output_dir = args.output_dir
    band_min = args.band_min
    band_max = args.band_max
    create_rgb = args.rgb
    rgb_bands = args.rgb_bands

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    if create_rgb:
        if rgb_bands is None:
            bands = [(200,700), (700, 1100), (1100, 2300)]
        else:
            bands = [(rgb_bands[0], rgb_bands[1]), (rgb_bands[2], rgb_bands[3]), (rgb_bands[4], rgb_bands[5])]
        
        print("\n========================================")
        print("Creating RGB image with the following bands:")
        for band in bands:
            print(f" * {band[0]} - {band[1]} eV")
        print("========================================\n")

        with open(f'{output_dir}/merged_image_RGB.log', 'w+') as log_file:
            for band in tqdm(bands):
                output_image = f'{output_dir}/merged_image_{band[0]}_{band[1]}.fits'
                output_expmap = f'{output_dir}/merged_expmap_{band[0]}_{band[1]}.fits'
                output_exp_corr = f'{output_dir}/merged_exp_corr_{band[0]}_{band[1]}.fits'
                
                # Run evtool for each band
                run_evtool(event_file, output_image, band[0]/1000, band[1]/1000, log_file=log_file)
                # Run expmap for each band
                run_expmap(event_file, output_image, output_expmap, band[0]/1000, band[1]/1000, log_file=log_file)
                # Perform exposure correction for each band
                exp_corr(output_image, output_expmap, output_exp_corr)
            
            log_file.seek(0)
            log_content = log_file.readlines()
            evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
            expmap_count = sum(1 for line in log_content if 'expmap: DONE' in line)
            if evtool_count == 3 and expmap_count == 3:
                print("\n========================================")
                print("All evtool and expmap tasks completed successfully for RGB Image")
                print("========================================\n")
        
    else:
        if band_min is None or band_max is None:
            parser.error("band_min and band_max are required when --rgb is not set.")
        
        with open(f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.log', 'w+') as log_file:
            
            print("\n========================================")
            print(f"Creating image with band: {band_min} - {band_max} eV for the file")
            print(f"{event_file}")
            print("========================================\n")
            # Run evtool for the specified band
            run_evtool(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                       band_min/1000, band_max/1000, log_file=log_file)
            
            print("\n========================================")
            print(f"Creating exposure map for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
            print("========================================\n")
            # Run expmap for the specified band
            run_expmap(event_file, f'{output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits', 
                       f'{output_dir}/merged_expmap_{int(band_min)}_{int(band_max)}.fits', 
                       band_min/1000, band_max/1000, log_file=log_file)
            
            print("\n========================================")
            print(f"Creating exposure corrected image for the file {output_dir}/merged_image_{int(band_min)}_{int(band_max)}.fits")
            print("========================================\n")
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
                print("Image creation successful!")
                print("========================================\n")
