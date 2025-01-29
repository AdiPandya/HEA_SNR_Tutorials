from matplotlib import pyplot as plt
import matplotlib
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
import numpy as np
import numpy.ma as ma
import math
import pyregion
import argparse
import os
from astropy.coordinates import Angle

def WriteRegListToFiles(sources, outroot, psradius, with_names):

    outfile = outroot + '.reg'    
    outexcl = outroot + '_exclude.reg'    
    outfits = outroot + '.fits'    


    # add num to sources
    regions_nums = [i+1 for i in range(len(sources["id"]))]

    ds9_regs = ""

    for source in sources:
        
        if psradius == "HEW" or psradius == "hew":
            #ds9_regs += 'fk5;circle(%f,%f,%f") # tag={%s} tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"], source["id"], source["EXT_LIKE"])
            ds9_regs += 'fk5;circle(%f,%f,%f) # tag={%s} tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"]/3600.0, source["id"], source["EXT_LIKE"])
        elif psradius == "cat" or psradius == "CAT":
            #ds9_regs += 'fk5;circle(%f,%f,%f") # tag={%s} tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"], source["id"], source["EXT_LIKE"])
            ds9_regs += 'fk5;circle(%f,%f,%f) # tag={%s} tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"]/3600.0, source["id"], source["EXT_LIKE"])
        else:
            # try number
            try: float(psradius)
            except:
                print("Wrong psradius format, possible are numbers (arcsec unit), hew or cat")
                break
            
            ds9_regs += 'fk5;circle(%f,%f,%f) # tag={%s} tag={extended=%i} \n' % (source["ra"], source["dec"], float(psradius)/3600.0, source["id"], source["EXT_LIKE"])


    with open(outfile, "w") as of:
        of.write(ds9_regs)

    print("Written output regfile", outfile)

    ds9_regs = ""

    for i, source in enumerate(sources):
        
        if not with_names:
            if psradius == "HEW" or psradius == "hew":
                #ds9_regs += 'fk5;-circle(%f,%f,%f") # tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"], source["EXT_LIKE"])
                ds9_regs += 'fk5;-circle(%f,%f,%f) # tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"]/3600.0, source["EXT_LIKE"])
            elif psradius == "cat" or psradius == "CAT":
                #ds9_regs += 'fk5;-circle(%f,%f,%f") # tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"], source["EXT_LIKE"])
                ds9_regs += 'fk5;-circle(%f,%f,%f) # tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"]/3600.0, source["EXT_LIKE"])
            else:
                try: float(psradius)
                except:
                    print("Wrong psradius format, possible are numbers (arcsec unit), hew or cat")
                    break
            
                ds9_regs += 'fk5;-circle(%f,%f,%f) # tag={extended=%i} \n' % (source["ra"], source["dec"], float(psradius)/3600.0, source["EXT_LIKE"])

        else:
            if psradius == "HEW" or psradius == "hew":
                #ds9_regs += 'fk5;-circle(%f,%f,%f") # text={%i} tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"], i+1, source["EXT_LIKE"])
                ds9_regs += 'fk5;-circle(%f,%f,%f) # text={%i} tag={extended=%i} \n' % (source["ra"], source["dec"], source["rhew"]/3600.0, i+1, source["EXT_LIKE"])
            elif psradius == "cat" or psradius == "CAT":
                #ds9_regs += 'fk5;-circle(%f,%f,%f") # text={%i} tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"], i+1, source["EXT_LIKE"])
                ds9_regs += 'fk5;-circle(%f,%f,%f) # text={%i} tag={extended=%i} \n' % (source["ra"], source["dec"], source["r"]/3600.0, i+1, source["EXT_LIKE"])
            else:
                try: float(psradius)
                except:
                    print("Wrong psradius format, possible are numbers (arcsec unit), hew or cat")
                    break

                ds9_regs += 'fk5;-circle(%f,%f,%f) # text={%i} tag={extended=%i} \n' % (source["ra"], source["dec"], float(psradius)/3600.0, i+1, source["EXT_LIKE"])
                      


    with open(outexcl, "w") as of:
        of.write(ds9_regs)

    print("Written output regfile", outexcl)

    columns = []

    columns.append(fits.Column(name='NUM', array=np.asarray(regions_nums), format='J'))
    columns.append(fits.Column(name='RA', array=sources['ra'], format='D'))
    columns.append(fits.Column(name='DEC', array=sources['dec'], format='D'))
    columns.append(fits.Column(name='R_ERR3', array=sources['r'], format='D', unit='arcsec'))
    columns.append(fits.Column(name='R_HEW2', array=sources['rhew'], format='D', unit='arcsec'))
    columns.append(fits.Column(name='ID', array=sources['id'], format='A40'))

    try:
        columns.append(fits.Column(name='ML_CTS_0', array=sources['cts_0'], format='D'))
        columns.append(fits.Column(name='ML_CTS_ERR_0', array=sources['cerr_0'], format='D'))
        columns.append(fits.Column(name='ML_RATE_0', array=sources['rate_0'], format='D'))
        columns.append(fits.Column(name='ML_CTS_1', array=sources['cts_1'], format='D'))
        columns.append(fits.Column(name='ML_CTS_ERR_1', array=sources['cerr_1'], format='D'))
        columns.append(fits.Column(name='ML_RATE_1', array=sources['rate_1'], format='D'))
        columns.append(fits.Column(name='ML_CTS_2', array=sources['cts_2'], format='D'))
        columns.append(fits.Column(name='ML_CTS_ERR_2', array=sources['cerr_2'], format='D'))
        columns.append(fits.Column(name='ML_RATE_2', array=sources['rate_2'], format='D'))
        columns.append(fits.Column(name='ML_CTS_3', array=sources['cts_3'], format='D'))
        columns.append(fits.Column(name='ML_CTS_ERR_3', array=sources['cerr_3'], format='D'))
        columns.append(fits.Column(name='ML_RATE_3', array=sources['rate_3'], format='D'))
    except:
        pass

    columns.append(fits.Column(name='DET_LIKE_0', array=sources["DET_LIKE_0"] , format='D'))
    try: columns.append(fits.Column(name='DET_LIKE_1', array=sources["DET_LIKE_1"] , format='D'))
    except: columns.append(fits.Column(name='DET_LIKE_P1', array=sources["DET_LIKE_P1"] , format='D'))

    try:
        columns.append(fits.Column(name='DET_LIKE_P2', array=sources["DET_LIKE_P2"] , format='D'))
        columns.append(fits.Column(name='DET_LIKE_P3', array=sources["DET_LIKE_P3"] , format='D'))
    except:
        pass
    
    columns.append(fits.Column(name='EXTENDED', array=sources['EXT_LIKE'], format='D'))

    t = fits.BinTableHDU.from_columns(columns)
    t.writeto(outfits, overwrite=True)

    print("Written output fitsfile", outfits)

def ParseRegFile(filename):

    regions = []

    if os.path.isfile(filename):
        f = open(filename, "r")

        for line in f:
            if line.startswith("#") or line.startswith("global") or line.strip() == "fk5" or line.startswith("image") or line.strip() == "":
                continue
            line = line.strip()

            reg = {}

            if line.startswith("fk5;"):
                line = line.replace("fk5;", "")

            if line.startswith("-"):
                reg["exclude"] = True
                line = line[1:]
            else:
                reg["exclude"] = False
            # get shape
            reg["shape"] = line[:line.index("(")]
            reg["comment"] = ""
            
            # get coordinates and radius with formated line
            coords = line[line.index("(")+1:].replace(")", "")
            
            if ":" in coords: sexa = True
            else: sexa = False

            coords = coords.split(",")

            # comment support
            if "#" in line:
                reg["comment"] = line[line.index("#")-1:]
                line = line[:line.index("#")]
                coords = line[line.index("(")+1:].replace(")", "")

                coords = coords.split(",")

            if reg["shape"] == "circle":
                if sexa: reg["ra"] = Angle(coords[0] + " hours").degree
                else: reg["ra"] = float(coords[0])

                if sexa: reg["dec"] = Angle(coords[1] + " degrees").degree
                else: reg["dec"] = float(coords[1])

                reg["r"] = float(coords[2].replace('"', ""))/3600.00
            elif reg["shape"] == "ellipse":

                if sexa: reg["ra"] = Angle(coords[0] + " hours").degree
                else: reg["ra"] = float(coords[0])

                if sexa: reg["dec"] = Angle(coords[1] + " degrees").degree
                else: reg["dec"] = float(coords[1])

                reg["r1"] = float(coords[2].replace('"', ""))/3600.00
                reg["r2"] = float(coords[3].replace('"', ""))/3600.00
                reg["angle"] = float(coords[4])
            elif reg["shape"] == "polygon":
                reg["ra_coords"] = []
                reg["dec_coords"] = []

                #print(reg.coord_list)
                for i in range(0, len(coords), 2):
                    if sexa: reg["ra_coords"].append(Angle(coords[i] + " hours").degree  )                      
                    else: reg["ra_coords"].append(float(coords[i]))

                    if sexa: reg["dec_coords"].append(Angle(coords[i+1] + " degrees").degree)
                    else: reg["dec_coords"].append(float(coords[i+1]))

                reg["ra"] = np.mean(reg["ra_coords"])
                reg["dec"] = np.mean(reg["dec_coords"])
                reg["r_enclosing"] = max( max(reg["ra_coords"])-min(reg["ra_coords"]), max(reg["dec_coords"])-min(reg["dec_coords"])  )/2*1.5 # radius = min max extend/2

            else:
                print("Unsupported region shape. Supported are: circles, ellipses, polygons")

            regions.append(reg)

        f.close()

    return regions

def ParseFitsTab(filename):

    print("Parsing fits source list %s" % filename)

    if os.path.isfile(filename):
        hdu_list = fits.open(filename,memmap=True)

        #print(hdu_list[1].columns)
        increase_factor = 1

        regions = Table(hdu_list[1].data)

        regions["ra"] = regions['RA']
        regions["dec"] = regions['DEC']
        regions['r'] = 3.*regions['RADEC_ERR']
        # create a column with 2*HEW/2
        regions['rhew'] = 2*28./2.*increase_factor
        try: regions['id'] = regions['DETUID']
        except: regions['id'] = regions['ID_SRC']
        try: regions['hr1'] = regions['HR_1']
        except: pass
        try: regions['hr1e'] = regions['HR_1_ERR']
        except: pass
        try: regions['hr2'] = regions['HR_2']
        except: pass
        try: regions['hr2e'] = regions['HR_2_ERR']
        except: pass
        try: regions['hr3'] = regions['HR_3']
        except: pass
        try: regions['hr3e'] = regions['HR_3_ERR']
        except: pass

        try: regions['cts_0'] = regions['ML_CTS_0']
        except: pass
        try: regions['cerr_0'] = regions['ML_CTS_ERR_0']
        except: pass
        try: regions['rate_0'] = regions['ML_RATE_0']
        except: pass

        try: regions['cts_1'] = regions['ML_CTS_1']
        except: pass
        try: regions['cerr_1'] = regions['ML_CTS_ERR_1']
        except: pass
        try: regions['rate_1'] = regions['ML_RATE_1']
        except: pass

        try: regions['cts_2'] = regions['ML_CTS_2']
        except: pass
        try: regions['cerr_2'] = regions['ML_CTS_ERR_2']
        except: pass
        try: regions['rate_2'] = regions['ML_RATE_2']
        except: pass

        try: regions['cts_3'] = regions['ML_CTS_3']
        except: pass
        try: regions['cerr_3'] = regions['ML_CTS_ERR_3']
        except: pass
        try: regions['rate_3'] = regions['ML_RATE_3']
        except: pass

        regions["DET_LIKE_0"] = regions["DET_LIKE_0"]

        try: regions["DET_LIKE_1"] = regions["DET_LIKE_1"]
        except: regions["DET_LIKE_P1"] = regions["DET_LIKE_P1"]

        try: regions["DET_LIKE_P2"] = regions["DET_LIKE_P2"]
        except: pass
        try: regions["DET_LIKE_P3"] = regions["DET_LIKE_P3"]
        except: pass

        regions["EXT_LIKE"] = regions["EXT_LIKE"]
        
        ni = np.in1d(regions['RADEC_ERR'],[0])
        #ni = np.isin(regions['RADEC_ERR'],[0])
        regions['r'][ni] = 1.0

        mask = regions['r'] <= 1.0
        regions['r'][mask] = 1.0

        ma0_arr = ma.masked_less(regions['EXT_LIKE'], 10.0)
        regions['EXT_LIKE'][ma0_arr.mask]=0

        ma1_arr = ma.masked_greater(regions['EXT_LIKE'], 9.9)
        regions['EXT_LIKE'][ma1_arr.mask]=1

        hdu_list.close()

    return regions


def select_regions(selregion, fitstable, outroot, psradius, det_like, with_name):

    if len(det_like) > 4:
        det_like = det_like[0:4]

    regs1 = ParseRegFile(selregion)
    regs2 = ParseFitsTab(fitstable)

    ra_col = []
    dec_col = []
    indices = []

    det_vals = [[] for i in range(len(det_like))]

    for i, reg in enumerate(regs2):
        ra_col.append(reg["ra"])
        dec_col.append(reg["dec"])
        indices.append(i)

        # add det_like cols
        if len(det_like) > 0:
            for d in range(len(det_like)):
                det_vals[d].append(reg["DET_LIKE_%i" % d])

    for d in range(len(det_like)):
        det_vals[d] = np.asarray(det_vals[d])

    ra_col = np.asarray(ra_col)
    dec_col = np.asarray(dec_col)

    valid_regs = None

    print("Filtering sources in source regions.")

    for reg in regs1:
        
        if reg["exclude"]:
            continue

        x_cent = reg["ra"]
        y_cent = reg["dec"]
        if reg["shape"] == "circle":
            radius = reg["r"]
        elif reg["shape"] == "ellipse":
            radius = max(reg["r1"], reg["r2"])
        elif reg["shape"] == "polygon":
            radius = reg["r_enclosing"]
            print("WARNING: Excluding all regions in a circle around polygon.")
        else:
            print("Region shape not supported. Supported are circle, ellipse, polygon")

        x_cent = np.full(ra_col.shape, x_cent)
        y_cent = np.full(dec_col.shape, y_cent)


        mask = np.rad2deg(np.arccos(np.sin(np.deg2rad(dec_col))*np.sin(np.deg2rad(y_cent))+np.cos(np.deg2rad(dec_col))
                *np.cos(np.deg2rad(y_cent))*np.cos(np.deg2rad(abs(ra_col-x_cent))))) < radius 

        if not np.any(valid_regs):
            valid_regs = np.asarray(indices)[mask]      
        else:
            valid_regs = np.append(valid_regs, np.asarray(indices)[mask])   

    valid_regs = np.sort(np.unique(valid_regs))  

    ra_good = ra_col[valid_regs]
    dec_good = dec_col[valid_regs]
    indices_good = np.asarray(indices)[valid_regs]
    

    masks_exclude = []
    # second go with exclude regions
    for reg in regs1:

        if not reg["exclude"]:
            continue

        x_cent = reg["ra"]
        y_cent = reg["dec"]
        if reg["shape"] == "circle":
            radius = reg["r"]
        elif reg["shape"] == "ellipse":
            radius = max(reg["r1"], reg["r2"])
            print("WARNING: Only circle exclude regions are exact.")
        elif reg["shape"] == "polygon":
            radius = reg["r_enclosing"]
            print("WARNING: Excluding all regions in a circle around polygon.")
        else:
            print("Region shape not supported. Supported are circle, ellipse, polygon")

        x_cent = np.full(ra_good.shape, x_cent)
        y_cent = np.full(dec_good.shape, y_cent)

        mask = np.rad2deg(np.arccos(np.sin(np.deg2rad(dec_good))*np.sin(np.deg2rad(y_cent))+np.cos(np.deg2rad(dec_good))
                *np.cos(np.deg2rad(y_cent))*np.cos(np.deg2rad(abs(ra_good-x_cent))))) > radius 
  
        masks_exclude.append(mask)

    # make master exclude region mask with logical and
    if len(masks_exclude) > 0:

        master_mask = masks_exclude[0]

        for i in range(1, len(masks_exclude)):
            master_mask = master_mask & masks_exclude[i]

        indices_good = indices_good[master_mask]

    valid_regs_good = np.sort(np.unique(indices_good))  
    
    for d in range(len(det_like)):
        det_vals[d] = det_vals[d][valid_regs_good]


    masks_det = []
    # lastly apply detection likelihood thresholds
    if len(det_like) > 0:

        for i in range(len(det_like)):

            if i > 3:
                print("Only up to 4 numbers supported.")
                break
            
            # check if number or skipping
            try:
                thresh = float(det_like[i])
                if thresh == 0:
                    raise ValueError
            except ValueError:
                print("No threshold for DET_LIKE_%i" % i)
                continue

            print("Applying threshold of %f for DET_LIKE_%i" % (thresh, i))

            thresh = np.full(det_vals[i].shape, thresh)


            mask = np.asarray(det_vals[i]) >= thresh

            masks_det.append(mask)

    # make master det_like mask with logical and
    if len(det_like) > 0:

        master_mask = masks_det[0]

        for i in range(1, len(masks_det)):
            master_mask = master_mask & masks_det[i]

        indices_good = indices_good[master_mask]


    valid_regs_good = np.sort(np.unique(indices_good))  

    good_regs = np.asarray(regs2)[valid_regs_good]

    WriteRegListToFiles(good_regs, outroot, psradius, with_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--selregion", help="Region file for the region inside which the sources should be located (circles, RA and Dec in degrees, radius in arcsec)", type=str, required=True)
    parser.add_argument("--fitstable", help="eRASS Source catalogue (FITS)", type=str, required=True)
    parser.add_argument("--outroot", required=True, type=str, help="Output file name (for the region files and for fits).")
    parser.add_argument("--psradius", required=False, type=str, help="What radius to use for the regions for the point sources. Input 'HEW' if 2 * HEW/2. Input cat to use reported positional error of catalog * 3 (Default). Input number for a fixed source radius (units arcseconds)", default="cat")
    parser.add_argument("--det_like", required=False, type=str, nargs="+", help="Detection likelihood threshold. Default is no threshold. Supports up to 4 numbers, e.g. '3 3 0 3' where 0 represents all sources. In this example sources with DET_LIKE_0, DET_LIKE_P1 and DET_LIKE_P3 greater 3  - with no threshold on DET_LIKE_P2 - would be included in the output. Separate numbers by space.", default=[])
    parser.add_argument("--with_name", required=False, action="store_true", help="Give this argument if output regions should be named with numbers.")

    args = parser.parse_args()
    select_regions(args.selregion, args.fitstable, args.outroot, args.psradius, args.det_like, args.with_name)
    
