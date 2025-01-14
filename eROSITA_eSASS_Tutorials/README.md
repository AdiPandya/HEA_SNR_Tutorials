# Tutorials for eROSITA data reduction using the eSASS software.

The tutorials are structured as follows:

## 1. [Introduction and Setup](01_Introduction_and_Setup.ipynb)
This contains the introduction to the eROSITA data reduction and common steps required for Image creation and spectral analysis.

The notebook covers the following steps:

- **Make Clean Event List**:
    - Define functions to run `evtool` and `radec2xy` commands.
    - Process event files to create cleaned event lists.

- **Extract the Lightcurve**:
    - Define a function to run `flaregti` for lightcurve extraction.
    - Process cleaned event lists to extract initial lightcurves.

- **Flare Filtering**:
    - Define functions for Gaussian fitting and sigma clipping.
    - Apply thresholds to lightcurves to identify and filter out flares.

- **Re-run `flaregti` with Calculated Thresholds**:
    - Re-extract lightcurves using the calculated thresholds to ensure flares are filtered out.

- **Re-run `evtool` using the New GTIs**:
    - Re-process event files using the new Good Time Intervals (GTIs) obtained from flare filtering.

- **Combine the Tiles**:
    - Merge the filtered event lists from different tiles into a single event list.
    - Separate the merged event list into individual Telescope Modules (TMs) for further analysis.

These steps provide a comprehensive guide to preparing eROSITA data for image creation and spectral analysis.

### Using the [Setup Script](Setup.py)

This script processes eROSITA data by performing several steps mentioned above. It uses parallel processing to speed up the data reduction process. Each file is processed in parallel to filter out flares and then merged into a single event list. 

#### Usage

To run the script, use the following command:
```bash
python Setup.py <input_dir> <output_dir> <timebin> <center_ra> <center_dec> [--ff_plots] [--ff_proof] [--separate_tm]
```

- `<input_dir>`: Input directory containing raw data.
- `<output_dir>`: Output directory for filtered data.
- `<timebin>`: Time bin size for lightcurve extraction.
- `<center_ra>`: Center RA in deg.
- `<center_dec>`: Center Dec in deg.

##### Optional Flags
- `--ff_plots`: Create flare filtering plots. Default is `True`.
- `--ff_proof`: Proof check flare filtering by creating lightcurves of the filtered event lists. Default is `False`.
- `--separate_tm`: Separate merged event list by Telescope Modules (TMs) 1 2 3 4 5 6 7 8 9. Default is `False`.

#### Example Usage

```bash
python Setup.py Data/Raw_data Data/Filtered_data 20 83.63 22.39 --ff_plots --separate_tm
```

This command will process the data in the specified input directory, save the results in the output directory, use a time bin size of 20, and set the center coordinates to RA 83.63 and DEC -5.39. It will also create log files and flare filtering plots.


## 2. [Image Creation](Imaging)
### Using the Imaging Script and Notebook

The `Imaging/Imaging.py` script and the `Imaging/02_Imaging_tutorial.ipynb` notebook process eROSITA data to generate images, exposure maps, and exposure-corrected images. Below are the steps to use the script and notebook:

#### Using the Imaging Script

- **Navigate to the Imaging Directory**:
    ```bash
    cd Imaging
    ```

- **Run the Imaging Script**:
    ```bash
    python Imaging.py <event_file> <output_dir> <band_min> <band_max> [--rgb] [--rgb_bands R_MIN R_MAX G_MIN G_MAX B_MIN B_MAX]
    ```

    - `<event_file>`: Input event file name.
    - `<output_dir>`: Output directory where the results will be saved.
    - `<band_min>`: Minimum energy band in eV (required if `--rgb` is not set).
    - `<band_max>`: Maximum energy band in eV (required if `--rgb` is not set).

- **Script Options**:
    - `--rgb`: Create RGB images using predefined or specified energy bands.
    - `--rgb_bands R_MIN R_MAX G_MIN G_MAX B_MIN B_MAX`: Energy bands for RGB image creation. If not provided, default bands (200-700 eV, 700-1100 eV, 1100-2300 eV) will be used.

- **Example Usage**:
    ```bash
    python Imaging.py event.fits output_dir 200 2000
    python Imaging.py event.fits output_dir --rgb --rgb_bands 200 700 700 1100 1100 2300
    ```

- **Script Description**:
    The script performs the following steps:
    - Creates an image using the specified energy band or RGB bands.
    - Generates an exposure map for the created image.
    - Creates an exposure-corrected image.

- **Functions**:
    - `run_evtool(input_name, output_name, emin, emax, gti_type='FLAREGTI', flag_type='0xe00fff30', pattern='15', telid='1 2 3 4 5 6 7', log_file=None)`: Runs the 'evtool' command to create an image from the event file.
    - `run_expmap(input_eventlist, input_image, output_name, emin, emax, log_file=None)`: Runs the 'expmap' command to generate an exposure map for the image.
    - `exp_corr(input_image, input_expmap, output_name)`: Creates an exposure-corrected image by dividing the counts image by the exposure map.

#### Using the Imaging Notebook

The [02_Imaging_tutorial.ipynb](Imaging/02_Imaging_tutorial.ipynb) notebook provides a tutorial of the python script [Imaging.py](Imaging/Imaging.py) to process the eROSITA data. It includes the following steps:

- **Creating Image and Exposure Map**: Use functions to create images and exposure maps.
- **Exposure Correction**: Generate exposure-corrected images.
- **Viewing the Image**: Visualize the images using interactive widgets.
- **RGB Image**: Create and visualize RGB images using predefined or specified energy bands.

For more detailed information, refer to the comments and documentation within the [Imaging.py](Imaging/Imaging.py) script and the [02_Imaging_tutorial.ipynb](Imaging/02_Imaging_tutorial.ipynb) notebook.



## 3. [Spectral Analysis](Spectra)