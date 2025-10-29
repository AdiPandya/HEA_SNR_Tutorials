#!/usr/bin/env python3
"""
Requirements:
- astropy
- numpy

Usage: python hidePS.py input.fits regions.reg output.fits
"""

import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import logging
import re
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# EXPAND_FACTOR = 1.2  # Expand factor for regions (times the original size)

class Region:
    """Region class for basic geometric shapes."""
    
    def __init__(self, shape_type, center_x, center_y, *params):
        self.shape_type = shape_type.lower()
        self.center_x = float(center_x)
        self.center_y = float(center_y)
        self.params = [float(p) for p in params]
    
    def contains_point(self, x, y):
        """Check if a point is inside the region."""
        dx = x - self.center_x
        dy = y - self.center_y
        
        if self.shape_type == 'circle':
            radius = self.params[0]
            return (dx*dx + dy*dy) <= radius*radius
        
        elif self.shape_type == 'ellipse':
            a = self.params[0]  # semi-major axis
            b = self.params[1]  # semi-minor axis
            angle = self.params[2] if len(self.params) > 2 else 0.0
            
            # Rotate coordinates
            cos_a = np.cos(np.radians(angle))
            sin_a = np.sin(np.radians(angle))
            dx_rot = dx * cos_a + dy * sin_a
            dy_rot = -dx * sin_a + dy * cos_a
            
            return (dx_rot*dx_rot)/(a*a) + (dy_rot*dy_rot)/(b*b) <= 1.0
        
        return False
    
    def create_expanded(self, expand_factor):
        """Create an expanded version of this region."""
        if self.shape_type == 'circle':
            new_params = [self.params[0] * expand_factor]
        elif self.shape_type == 'ellipse':
            new_params = [self.params[0] * expand_factor, 
                         self.params[1] * expand_factor]
            if len(self.params) > 2:
                new_params.append(self.params[2])  # Keep angle
        else:
            new_params = self.params
        
        return Region(self.shape_type, self.center_x, self.center_y, *new_params)


def parse_ds9_region_line(line):
    """Parse a single line from a DS9 region file."""
    # Remove comments and whitespace
    line = line.split('#')[0].strip()
    if not line:
        return None
    
    # Basic parsing for circle and ellipse
    circle_match = re.match(r'circle\s*\(\s*([^,]+),\s*([^,]+),\s*([^)]+)\)', line)
    if circle_match:
        x, y, r = circle_match.groups()
        return Region('circle', x, y, r)
    
    ellipse_match = re.match(r'ellipse\s*\(\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+)(?:,\s*([^)]+))?\)', line)
    if ellipse_match:
        x, y, a, b, angle = ellipse_match.groups()
        if angle is None:
            angle = 0.0
        return Region('ellipse', x, y, a, b, angle)
    
    logger.warning(f"Could not parse region line: {line}")
    return None


class RegionProcessor:
    """Process the region files to mask out Point sources."""
    
    def __init__(self, fits_file):
        """Initialize with a FITS file."""
        self.hdu_list = fits.open(fits_file)
        self.image_data = self.hdu_list[0].data.astype(np.float32)
        self.header = self.hdu_list[0].header
        
        # Handle different image dimensions
        if self.image_data.ndim == 2:
            self.image_2d = self.image_data
        elif self.image_data.ndim == 3:
            self.image_2d = self.image_data[0]
        elif self.image_data.ndim == 4:
            self.image_2d = self.image_data[0, 0]
        else:
            raise ValueError(f"Unsupported image dimensions: {self.image_data.ndim}")
        
        # Create output copy
        self.output_data = self.image_2d.copy()
        
        # Get image dimensions
        self.height, self.width = self.image_2d.shape
        
        # Set up WCS for coordinate transformations
        try:
            self.wcs = WCS(self.header, naxis=2)
            self.has_wcs = True
        except Exception as e:
            logger.warning(f"Could not parse WCS: {e}. Using pixel coordinates.")
            self.wcs = None
            self.has_wcs = False
    
    def pixel_to_world(self, x, y):
        """Convert pixel coordinates to world coordinates."""
        if self.wcs is None:
            return x, y
        try:
            # Convert to world coordinates
            world_coords = self.wcs.pixel_to_world(x, y)
            return world_coords.ra.deg, world_coords.dec.deg
        except:
            return x, y
    
    def world_to_pixel(self, ra, dec):
        """Convert world coordinates to pixel coordinates."""
        if self.wcs is None:
            return ra, dec
        try:
            pixel_coords = self.wcs.world_to_pixel_values(ra, dec)
            return pixel_coords[0], pixel_coords[1]
        except:
            return ra, dec
    
    def create_region_mask(self, region):
        """Create a boolean mask for pixels inside a region."""
        mask = np.zeros((self.height, self.width), dtype=bool)
        
        # Create coordinate grids
        y_coords, x_coords = np.mgrid[0:self.height, 0:self.width]
        
        # Check each pixel (this could be optimized for large images)
        for y in range(self.height):
            for x in range(self.width):
                if region.contains_point(x, y):
                    mask[y, x] = True
        
        return mask
    
    def fill_region_with_random_values(self, region):
        """Fill a region with random values from surrounding pixels."""
        # logger.info(f"Processing {region.shape_type} region at ({region.center_x:.1f}, {region.center_y:.1f})")
        
        # Create expanded region
        expanded_region = region.create_expanded(EXPAND_FACTOR)
        
        # Get masks
        region_mask = self.create_region_mask(region)
        expanded_mask = self.create_region_mask(expanded_region)
        
        # Get pixels in the expanded region but not in the original region (annulus)
        annulus_mask = expanded_mask & ~region_mask
        
        # Get values from the annulus
        annulus_values = self.image_2d[annulus_mask]
        
        if len(annulus_values) == 0:
            logger.warning(f"No pixels found in annulus for region. Skipping.")
            return
        
        # Remove NaN and infinite values
        finite_values = annulus_values[np.isfinite(annulus_values)]
        
        if len(finite_values) == 0:
            logger.warning(f"No finite pixels found in annulus for region. Skipping.")
            return
        
        # Get number of pixels to fill
        n_pixels_to_fill = np.sum(region_mask)
        
        if n_pixels_to_fill == 0:
            logger.warning(f"No pixels to fill in region. Skipping.")
            return
        
        # Generate random samples
        random_indices = np.random.randint(0, len(finite_values), n_pixels_to_fill)
        random_values = finite_values[random_indices]
        
        # Fill the region
        self.output_data[region_mask] = random_values
        
        # logger.info(f"Filled {n_pixels_to_fill} pixels with random values from {len(finite_values)} source pixels")
    
    def process_regions_file(self, region_file):
        """Process all regions in a DS9 region file."""
        try:
            with open(region_file, 'r') as f:
                lines = f.readlines()
            
            regions = []
            for line in lines:
                # Skip comments and coordinate system lines
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('global') or \
                   'fk5' in line.lower() or 'icrs' in line.lower() or 'image' in line.lower():
                    continue
                
                region = parse_ds9_region_line(line)
                if region:
                    regions.append(region)
            
            logger.info(f"Found {len(regions)} valid regions in {region_file}")
            
            # Process each region
            for i, region in tqdm(enumerate(regions), total=len(regions), desc="Processing regions"):
                # logger.info(f"Processing region {i+1}/{len(regions)}")
                self.fill_region_with_random_values(region)
                
        except Exception as e:
            logger.error(f"Error processing regions file {region_file}: {e}")
            raise
    
    def save_output(self, output_file):
        """Save the modified image to a new FITS file."""
        # Create new HDU with modified data
        if self.image_data.ndim == 2:
            new_data = self.output_data
        elif self.image_data.ndim == 3:
            new_data = self.image_data.copy()
            new_data[0] = self.output_data
        elif self.image_data.ndim == 4:
            new_data = self.image_data.copy()
            new_data[0, 0] = self.output_data
        
        # Create new HDU
        primary_hdu = fits.PrimaryHDU(data=new_data, header=self.header)
        
        # Add history
        primary_hdu.header['HISTORY'] = 'Processed with hideregions2_simple.py'
        primary_hdu.header['HISTORY'] = 'Point sources filled with random surrounding values'
        
        # Save to file
        hdu_list = fits.HDUList([primary_hdu])
        hdu_list.writeto(output_file, overwrite=True)
        logger.info(f"Saved output to {output_file}")
    
    def close(self):
        """Close the FITS file."""
        self.hdu_list.close()


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Fill regions in FITS images with random values from surrounding pixels.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Region file format:
  - DS9 region files (.reg) with circle() and ellipse() regions in pixel coordinates.
  - Example: circle(100.5, 200.3, 15.0)
  - Example: ellipse(150.2, 180.7, 20.0, 10.0, 45.0)
        """
    )
    
    parser.add_argument('input_fits', help='Input FITS image file')
    parser.add_argument('region_file', help='DS9 region file (.reg format)')
    parser.add_argument('output_fits', help='Output FITS image file')
    parser.add_argument('--expand-factor', type=float, default=1.2,
                       help='Factor to expand regions (default: 1.2)')

    args = parser.parse_args()
 
    # Set expand size
    global EXPAND_FACTOR
    EXPAND_FACTOR = args.expand_factor
    
    try:
        # Process the image
        processor = RegionProcessor(args.input_fits)
        processor.process_regions_file(args.region_file)
        processor.save_output(args.output_fits)
        processor.close()
        
        logger.info("Processing completed successfully!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
