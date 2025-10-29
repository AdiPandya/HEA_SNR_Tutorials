# hidePS (Python): Hide Point Sources in FITS Images

This directory contains:
- [hidePS.py](hidePS.py) — Python script to replace point-source regions in a FITS image with random values sampled from surrounding pixels.
- [hidePS_visualization.ipynb](hidePS_visualization.ipynb) — A tutorial notebook to visualize the hidePS workflow on a synthetic test image.

## What hidePS does

- Reads a FITS image (2D or the first slice of 3D/4D) and a DS9 region file (circle/ellipse).
- Expands each source region by a multiplicative factor to form an annulus and samples pixels from that annulus.
- Replaces source pixels with random draws from the annulus sample and writes the modified FITS to a new file.

## Setup

Install dependencies using the following command:

```bash
pip install -r requirements.txt
```

## Usage

Basic:
```bash
python hidePS.py input.fits regions.reg output.fits
```
Options:
```bash
python hidePS.py input.fits regions.reg output.fits --expand-factor FLOAT
```

Quick notes:
- Output file is overwritten if it exists.
- If annulus has no pixels, increase --expand-factor (e.g., 1.5–2.0).

### Regions (DS9, pixel coords)
- Supported: circle(x, y, r) and ellipse(x, y, a, b, angle)
- The script ignores lines starting with #, global, image, or containing fk5/icrs

## Tutorial notebook
[hidePS_visualization.ipynb](hidePS_visualization.ipynb) — creates a test FITS file with random point-sources and their region file, visualizes masks/annuli, runs the fill step-by-step, and shows before/after comparisons (plots + saved processed FITS).

## Author's Note
I hope you find the script and tutorial useful for your analysis needs. If you have any suggestions for improvements or new features, please feel free to reach out at my email id: aditya.pandya@astro.uni-tuebingen.de.