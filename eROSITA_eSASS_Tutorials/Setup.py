############ Setup the environment ############
import glob
import subprocess
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import warnings
import argparse

warnings.filterwarnings("ignore")

plt.style.use("default")
plt.rc('xtick', direction='in', top=True)
plt.rc('ytick', direction='in', right=True)
plt.rc('axes', linewidth=1.15)
plt.rc("mathtext", fontset="dejavuserif")

print("\n========================================\n")
print('Setting up the environment...')

parser = argparse.ArgumentParser(description='Process eROSITA data.')
parser.add_argument('input_dir', type=str, help='Input directory containing raw data')
parser.add_argument('output_dir', type=str, help='Output directory for filtered data')
parser.add_argument('timebin', type=str, default='20', help='Time bin size for lightcurve')
parser.add_argument('center_ra', type=str, help='Center RA')
parser.add_argument('center_dec', type=str, help='Center DEC')
parser.add_argument('--logs', action='store_true', default=True, help='Flag to create log files')
parser.add_argument('--ff_plots', action='store_true', default=True, help='Flag to create flare filtering plots')
parser.add_argument('--ff_proof', action='store_true', default=False, help='Flag to proof check flare filtering')
parser.add_argument('--separate_tm', action='store_true', default=False, help='Flag to separate merged event list by TMs')

args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
timebin = args.timebin
center_ra = args.center_ra
center_dec = args.center_dec
create_logs = args.logs
ff_plots = args.ff_plots
proof_check = args.ff_proof
separate_tm = args.separate_tm

elist = glob.glob(f'{input_dir}/???/???/EXP_010/e?01_??????_020_EventList_c010.fits.gz')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(f'{output_dir}/Lightcurves'):
    os.makedirs(f'{output_dir}/Lightcurves')

if not os.path.exists(f'{output_dir}/Merged'):
    os.makedirs(f'{output_dir}/Merged')

print('Created directories for output files:')
print(f'  {output_dir}')
print(f'  {output_dir}/Lightcurves')
print(f'  {output_dir}/Merged')

clean_list = np.empty(len(elist), dtype=object)
lightcurve0_list = np.empty(len(elist), dtype=object)
lightcurve_list = np.empty(len(elist), dtype=object)
filtered_list = np.empty(len(elist), dtype=object)

for i in range(len(elist)):
    clean_list[i] = f'{output_dir}/' + elist[i].split('/')[-1].replace('EventList_c010.fits.gz', 'c010_s01_CleanedEvents.fits')
    lightcurve0_list[i] = f'{output_dir}/Lightcurves/' + elist[i].split('/')[-1].replace('EventList_c010.fits.gz', f'c010_s02_LC0_tb{timebin}.fits')
    lightcurve_list[i] = f'{output_dir}/Lightcurves/' + elist[i].split('/')[-1].replace('EventList_c010.fits.gz', f'c010_s03_LC_tb{timebin}.fits')
    filtered_list[i] = f'{output_dir}/' + elist[i].split('/')[-1].replace('EventList_c010.fits.gz', 'c010_s04_FlareFilteredEvents.fits')

with open(f'{output_dir}/filtered.list', 'w') as f:
    for e in filtered_list:
        f.write(f'{e}\n')

print('\nEnvironment setup complete\n')
print('Starting the workflow with the following parameters:')
print(f'  Input directory: {input_dir}')
print(f'  Output directory: {output_dir}')
print(f'  Time bin size: {timebin}')
print(f'  Center RA: {center_ra}')
print(f'  Center DEC: {center_dec}')
print(f'  Create logs: {create_logs}')
print(f'  Create flare filtering plots: {ff_plots}')
print(f'  Proof check flare filtering: {proof_check}')
print(f'  Separate merged event list by TMs: {separate_tm}')

############ Define functions to run evtool, radec2xy and flaregti ############

def run_evtool(input_name, output_name, gti_type='GTI', flag_type='0xe00fff30', pattern='15', emin='0.2', emax='10.0', image='no', events='yes', telid='1 2 3 4 5 6 7', log_file=None):
    subprocess.run(['evtool', 
                    f'eventfiles={input_name}', 
                    f'outfile={output_name}', 
                    f'gti={gti_type}', 
                    f'flag={flag_type}', 
                    f'pattern={pattern}', 
                    f'emin={emin}', 
                    f'emax={emax}',
                    f'image={image}',
                    f'events={events}',
                    f'telid={telid}'
                    ],
                    stdout=log_file,
                    stderr=log_file)

def run_radec2xy(input_name, ra, dec, log_file=None):
    subprocess.run(['radec2xy', 
                    f'{input_name}', 
                    f'ra0={ra}', 
                    f'dec0={dec}'],
                    stdout=log_file,
                    stderr=log_file)
    
def run_flaregti(input_name, output_lightcurve, pimin='5000', source_size='150', gridsize='26', timebin=timebin, threshold='-1', log_file=None):
    subprocess.run(['flaregti', 
                    f'{input_name}', 
                    f'pimin={pimin}', 
                    f'source_size={source_size}',
                    f'gridsize={gridsize}',
                    f'lightcurve={output_lightcurve}',
                    'write_mask=no',
                    f'timebin={timebin}',
                    f'threshold={threshold}'
                    ],
                    stdout=log_file,
                    stderr=log_file)

############ Making clean Event list for all tiles ############
print("\n========================================\n")
print('1) Creating clean event list for all tiles:')
if create_logs:
    with open(f'{output_dir}/evtool_s01.log', 'w+') as log_file:    
        for tile in tqdm(range(len(elist))):   
            run_evtool(elist[tile], clean_list[tile], log_file=log_file)

        log_file.seek(0)
        log_content = log_file.readlines()
        evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
        print(f'evtool finished successfully for {evtool_count} out of {len(elist)} files ({evtool_count/len(elist)*100}%)')
    print(f'Log file saved as {output_dir}/evtool_s01.log')
else:
    for tile in tqdm(range(len(elist))):   
        run_evtool(elist[tile], clean_list[tile])

############ Extract Lightcurves ############
print("\n========================================\n")
print('2) Extracting lightcurves for all tiles:')
if create_logs:
    with open(f'{output_dir}/Lightcurves/flaregti_s02.log', 'w+') as log_file:
        for tile in tqdm(range(len(elist))):   
            run_flaregti(clean_list[tile], lightcurve0_list[tile], log_file=log_file)
        
        log_file.seek(0)
        log_content = log_file.readlines()
        flaregti_count = sum(1 for line in log_content if 'flaregti: DONE' in line)
        print(f'flaregti finished successfully for {flaregti_count} out of {len(elist)} files ({flaregti_count/len(elist)*100}%)')
    print(f'Log file saved as {output_dir}/Lightcurves/flaregti_s02.log')
else:
    for tile in tqdm(range(len(elist))):   
        run_flaregti(clean_list[tile], lightcurve0_list[tile])

############ Flare Filtering functions ############

def gaussian(x, amplitude, mean, stdev):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stdev**2))

def fit_gaussian(data, bins='auto'):
    bin_heights, bin_borders = np.histogram(data, bins=bins)
    bin_centers = (bin_borders[:-1] + bin_borders[1:]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, bin_heights, p0=[1., np.mean(data), np.std(data)])
    return popt, pcov, bin_borders

def sigma_clipping(data, popt):
    data_std = popt[2]
    data_mean = popt[1]
    sigma_threshold = data_std * 3
    lower_limit  = data_mean - sigma_threshold 
    upper_limit = data_mean + sigma_threshold
    clipped_data = data[(data >= lower_limit) & (data <= upper_limit)]
    return clipped_data, lower_limit, upper_limit

def threshold_lightcurve(input_data, ff_plots=ff_plots, output_dir=f'{output_dir}/Lightcurves/'):
    lightcurve = fits.open(input_data)
    time = lightcurve[1].data['TIME']
    rate = lightcurve[1].data['RATE']
    positive_rate = rate[rate > 0]

    popt_rate, pcov_rate, rate_borders = fit_gaussian(rate)
    popt_pos, _, _ = fit_gaussian(positive_rate)
    clipped_data, lower_limit, upper_limit = sigma_clipping(positive_rate, popt_pos)
    popt_clip, pcov_clip, _ = fit_gaussian(clipped_data, bins=rate_borders)

    if ff_plots:
        plt.rc('font', family='DejaVu Serif', size=11)

        fig, ax = plt.subplots(3, 1, figsize=(8, 7))
        fig.subplots_adjust(hspace=0.4)  

        main_color = 'tab:red'
        clipping_region_color = 'k'

        # Plot 1
        ax[0].plot((time - time[0]) / 1e3, rate, lw=1.5, color=main_color)
        ax[0].set_ylabel('Rate \n $[\\mathrm{cts\\ s^{-1}\\ deg^{-2}}]$')
        ax[0].set_xlabel('Time [ks]')
        ax[0].axhline(popt_rate[1], color=clipping_region_color, linestyle='--', label='Mean')
        ax[0].axhspan(lower_limit, upper_limit, color=clipping_region_color, alpha=0.3, label='Clipping Region')
        ax[0].legend()

        # Plot 2
        ax[1].plot(np.arange(0, len(rate)) * 10 / 1e3, rate, lw=1.5, color=main_color)
        ax[1].set_ylabel('Rate \n $[\\mathrm{cts\\ s^{-1}\\ deg^{-2}}]$')
        ax[1].set_xlabel('Time [ks]')
        ax[1].axhline(popt_rate[1], color=clipping_region_color, linestyle='--', label='Mean')
        ax[1].axhspan(lower_limit, upper_limit, color=clipping_region_color, alpha=0.3, label='Clipping Region')
        ax[1].legend()

        # Plot 3
        ax[2].hist(rate, bins=rate_borders, alpha=0.5, label='Data', color=main_color)
        x_fit_interval = np.linspace(rate_borders[0], rate_borders[-1], 100)
        ax[2].axvspan(lower_limit, upper_limit, color='steelblue', alpha=0.25)
        ax[2].hist(clipped_data, bins=rate_borders, alpha=0.75, label='Clipped Data', color='tab:blue')
        ax[2].plot(x_fit_interval, gaussian(x_fit_interval, *popt_clip), label='Fitted Gaussian (Clipped)', color='tab:red')

        ax[2].set_xlabel('Rate $[\\mathrm{cts\\ s^{-1}\\ deg^{-2}}]$')
        ax[2].set_ylabel('Counts')
        ax[2].legend()

        fig.savefig(output_dir + input_data.split('/')[-1].replace('.fits', '.png'), dpi=300)
        plt.close(fig)  # Close the figure to avoid displaying it
    
    return upper_limit

############ Flare Filtering ############

print('\n2.1) Calculating Thresholds for flare filtering:')
tile_thresholds = np.zeros(len(lightcurve0_list))
for tile in tqdm(range(len(lightcurve0_list))):   
    tile_thresholds[tile] = threshold_lightcurve(lightcurve0_list[tile])
print(f'Flare filtering plot saved in {output_dir}/Lightcurves/')

print("\n========================================\n")
print('3) Running flaregti for all tiles with the calculated thresholds:')
if create_logs:
    with open(f'{output_dir}/Lightcurves/flaregti_s03.log', 'w+') as log_file:
        for tile in tqdm(range(len(elist))):   
            run_flaregti(clean_list[tile], lightcurve_list[tile], threshold=tile_thresholds[tile], log_file=log_file)
        
        log_file.seek(0)
        log_content = log_file.readlines()
        flaregti_count = sum(1 for line in log_content if 'flaregti: DONE' in line)
        print(f'flaregti finished successfully for {flaregti_count} out of {len(elist)} files ({flaregti_count/len(elist)*100}%)')
    print(f'Log file saved as {output_dir}/Lightcurves/flaregti_s03.log')
else:
    for tile in tqdm(range(len(elist))):   
        run_flaregti(clean_list[tile], lightcurve_list[tile], threshold=tile_thresholds[tile])

print("\n========================================\n")
print('4) Running evtool for all tiles with the flare filtered lightcurves:')
if create_logs:
    with open(f'{output_dir}/evtool_s04.log', 'w+') as log_file:    
        for tile in tqdm(range(len(elist))):   
            run_evtool(clean_list[tile], filtered_list[tile], gti_type="FLAREGTI", log_file=log_file)

        log_file.seek(0)
        log_content = log_file.readlines()
        evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
        print(f'evtool finished successfully for {evtool_count} out of {len(clean_list)} files ({evtool_count/len(clean_list)*100}%)')
    print(f'Log file saved as {output_dir}/evtool_s04.log')
else:    
    for tile in tqdm(range(len(elist))):   
        run_evtool(clean_list[tile], filtered_list[tile], gti_type="FLAREGTI")

if proof_check:
    print('\n4.1) Proof checking flare filtering:')
    if not os.path.exists(f'{output_dir}/Lightcurves/Proof_check'):
        os.makedirs(f'{output_dir}/Lightcurves/Proof_check')

    pc_lightcurve_list = np.empty(len(elist), dtype=object)

    for i in range(len(filtered_list)):
        pc_lightcurve_list[i] = f'{output_dir}/Lightcurves/Proof_check/' + filtered_list[i].split('/')[-1].replace('s04_FlareFilteredEvents.fits', f's041_pcLC_tb{timebin}.fits')

    if create_logs:
        with open(f'{output_dir}/Lightcurves/Proof_check/flaregti_s041.log', 'w+') as log_file:
            for tile in tqdm(range(len(elist))):   
                run_flaregti(filtered_list[tile], pc_lightcurve_list[tile], log_file=log_file)
            
            log_file.seek(0)
            log_content = log_file.readlines()
            flaregti_count = sum(1 for line in log_content if 'flaregti: DONE' in line)
            print(f'flaregti finished successfully for {flaregti_count} out of {len(elist)} files ({flaregti_count/len(elist)*100}%)')
        print(f'Log file saved as {output_dir}/Lightcurves/Proof_check/flaregti_s041.log')
    else:
        for tile in tqdm(range(len(elist))):   
            run_flaregti(filtered_list[tile], pc_lightcurve_list[tile])

    for tile in tqdm(range(len(pc_lightcurve_list))):   
        threshold_lightcurve(pc_lightcurve_list[tile], image=True, output_dir=f'{output_dir}/Lightcurves/Proof_check/')
    print(f'Plots for proof checking saved in {output_dir}/Lightcurves/Proof_check/')

print("\n========================================\n")
print('5) Merging tiles event list:')
if create_logs:
    with open(f'{output_dir}/Merged/merged_evtool_s05.log', 'w+') as log_file:    
        run_evtool(f'@{output_dir}/filtered.list', f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', gti_type="FLAREGTI", log_file=log_file)
        run_radec2xy(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', center_ra, center_dec, log_file=log_file)

        log_file.seek(0)
        log_content = log_file.readlines()
        evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
        radec2xy_count = sum(1 for line in log_content if 'radec2xy: DONE' in line)
        if evtool_count == 1 and radec2xy_count == 1:
            print('Merged tiles eventlist successfully')
    print(f'Log file saved as {output_dir}/Merged/merged_evtool_s05.log')
else:
    run_evtool(f'@{output_dir}/filtered.list', f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', gti_type="FLAREGTI")
    run_radec2xy(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', center_ra, center_dec)

if separate_tm:
    print('\n6) Separating merged event list into TM 1 2 3 4 5 6 7 8 9:')
    TM_list = np.array([1, 2, 3, 4, 5, 6, 7])

    if create_logs:
        with open(f'{output_dir}/Merged/separate_TM_evtool_s05.log', 'w+') as log_file:
            for i in tqdm(range(len(TM_list))):
                run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_{TM_list[i]}20_s05_TM{TM_list[i]}_Events.fits', telid=f'{TM_list[i]}', log_file=log_file)

            run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_820_s05_TM8_Events.fits', telid='1 2 3 4 6', log_file=log_file)

            run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_920_s05_TM9_Events.fits', telid='5 7', log_file=log_file)

            log_file.seek(0)
            log_content = log_file.readlines()
            evtool_count = sum(1 for line in log_content if 'evtool: DONE' in line)
            print(f'evtool successfully separated file into {evtool_count} files for {len(TM_list) + 2} TMs ({evtool_count / (len(TM_list) + 2) * 100}%)')

        print(f'Log file saved as {output_dir}/Merged/separate_TM_evtool_s05.log')
    else:
        for i in tqdm(range(len(TM_list))):
            run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_{TM_list[i]}20_s05_TM{TM_list[i]}_Events.fits', telid=f'{TM_list[i]}')

        run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_820_s05_TM8_Events.fits', telid='1 2 3 4 6')

        run_evtool(f'{output_dir}/Merged/Merged_020_s05_TM0_Events.fits', f'{output_dir}/Merged/Merged_920_s05_TM9_Events.fits', telid='5 7')

print("\n========================================\n")
print('** All tasks completed successfully **')
print("\n========================================\n")