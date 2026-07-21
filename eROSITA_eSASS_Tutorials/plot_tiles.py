"""
plot_tiles.py  -  Overlay eROSITA survey-tile footprints on an HiPS background.

All tiles from the SMAPS catalogue are drawn as thin outlines; the tiles
actually present in the input event-file list are highlighted.

Usage:
    python plot_tiles.py <skymap.fits> <event_files> <output.png>
                         [--hips HIPS_ID] [--nx NX] [--ny NY]

Called programmatically:
    from plot_tiles import plot
    plot(skymap_path, event_files, output_path)
    plot(skymap_path, event_files, output_path,
         hips_id='erosita/dr1/rate/024', nx=1800, ny=900)

Arguments:
    skymap       : path to SKYMAPS_052022_MPE.fits (tile catalogue)
    event_files  : list/iterable of raw event-file paths  -or-  path to a
                   text file whose lines are event-file paths
    output_path  : path for the output PNG image
    --hips       : HiPS identifier (default: erosita/dr1/rate/024)
    --nx / --ny  : output image dimensions in pixels (default: 1800 x 900)
"""

import sys
import os
import glob
import argparse
import logging
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
import astropy.units as u


# --------------------------------------------------------------------------- #
#  Helpers                                                                     #
# --------------------------------------------------------------------------- #


DEFAULT_HIPS = 'erosita/dr1/rate/024'

# In ZEA, the projection-plane distance at 90° from the centre is:
#   R = (180/π) * 2 * sin(45°) = (180/π) * √2  degrees
_ZEA_R90 = (180.0 / np.pi) * 2.0 * np.sin(np.pi / 4.0)  # ≈ 81.03 deg


def _make_halfsky_wcs(nx, glon_cen, glat_cen):
    """Return a Galactic ZEA WCS centred on (glon_cen, glat_cen).

    The hemisphere within 90° of the centre fills a circle of diameter *nx*
    pixels.  Image is square (nx x nx).
    """
    cdelt = _ZEA_R90 / (nx / 2.0)
    w = WCS(naxis=2)
    w.wcs.ctype = ['GLON-ZEA', 'GLAT-ZEA']
    w.wcs.crval = [float(glon_cen % 360), float(glat_cen)]
    w.wcs.crpix = [nx / 2.0 + 0.5, nx / 2.0 + 0.5]
    w.wcs.cdelt = [-cdelt, cdelt]
    w.pixel_shape = (nx, nx)
    return w


def _catalog_centroid(smaps):
    """Return (glon, glat) of the vector mean of all tile centres."""
    l = np.radians(smaps['GLON_CEN'].astype(float))
    b = np.radians(smaps['GLAT_CEN'].astype(float))
    x = np.mean(np.cos(b) * np.cos(l))
    y = np.mean(np.cos(b) * np.sin(l))
    z = np.mean(np.sin(b))
    glon = float(np.degrees(np.arctan2(y, x)) % 360)
    glat = float(np.degrees(np.arcsin(z / np.sqrt(x**2 + y**2 + z**2))))
    return glon, glat


def _fetch_hips(hips_id, wcs, logger):
    """Fetch an HiPS image via hips2fits and reproject to the target WCS.

    hips2fits silently returns all-zeros for Galactic coordinate types and
    for ZEA with fov >= ~180°.  We therefore fetch as an ICRS AIT all-sky
    image (which works reliably) and reproject into the caller's WCS.
    """
    from astroquery.hips2fits import hips2fits as _h2f
    from reproject import reproject_interp

    nx = wcs.pixel_shape[0]
    ny = wcs.pixel_shape[1]

    # ICRS AIT all-sky WCS (standard 2:1 aspect, no coordinate issues)
    fetch_nx = nx
    fetch_ny = nx // 2
    fetch_wcs = WCS(naxis=2)
    fetch_wcs.wcs.ctype = ['RA---AIT', 'DEC--AIT']
    fetch_wcs.wcs.crval = [0.0, 0.0]
    fetch_wcs.wcs.crpix = [fetch_nx / 2.0 + 0.5, fetch_ny / 2.0 + 0.5]
    fetch_wcs.wcs.cdelt = [-360.0 / fetch_nx, 180.0 / fetch_ny]
    fetch_wcs.pixel_shape = (fetch_nx, fetch_ny)

    logger.info(f'Fetching HiPS background: {hips_id} …')
    _h2f.timeout = 120
    result = _h2f.query_with_wcs(hips=hips_id, wcs=fetch_wcs, format='fits')
    raw = result.data.astype(float) if hasattr(result, 'data') else result[0].data.astype(float)

    # Reproject ICRS AIT → Galactic ZEA so pixels align with tile outlines
    data, _ = reproject_interp((raw, fetch_wcs), wcs, shape_out=(ny, nx))

    n_fin = int(np.isfinite(data).sum())
    n_pos = int((data > 0).sum())
    return data


def _skyfield_from_file(path):
    """Return the integer SKYFIELD value from an event-file primary header."""
    try:
        with fits.open(path, memmap=True) as hdul:
            return int(hdul[0].header['SKYFIELD'])
    except Exception:
        return None


def _tile_polygon_pixels(ra_min, ra_max, de_min, de_max, wcs, ny, nx, npts=8):
    """
    Convert a tile (RA/Dec rectangle) to pixel-space polygon vertices.


    Returns None for tiles that cross the Galactic anti-meridian (l ≈ 180°)
    or whose vertices project to NaN (outside the AIT ellipse).
    """
    def _edge(r0, r1, d0, d1):
        return np.linspace(r0, r1, npts), np.linspace(d0, d1, npts)

    edges = [
        _edge(ra_min, ra_max, de_min, de_min),
        _edge(ra_max, ra_max, de_min, de_max),
        _edge(ra_max, ra_min, de_max, de_max),
        _edge(ra_min, ra_min, de_max, de_min),
    ]
    ra_all  = np.concatenate([e[0] for e in edges])
    dec_all = np.concatenate([e[1] for e in edges])

    sky = SkyCoord(ra=ra_all, dec=dec_all, unit='deg', frame='icrs')
    gal = sky.galactic
    l, b = gal.l.deg, gal.b.deg

    # Normalise longitudes around the tile centre so they are contiguous
    l_c = l.copy()
    med = np.median(l_c)
    l_c[l_c - med >  180] -= 360
    l_c[l_c - med < -180] += 360

    try:
        px, py = wcs.all_world2pix(l_c, b, 0)
    except Exception:
        return None

    if np.any(np.isnan(px)) or np.any(np.isnan(py)):
        return None

    # Tiles straddling the anti-meridian (l ≈ 180°) have vertices that jump
    # from ~0 to ~nx in pixel space — skip them to avoid spurious lines.
    if np.ptp(px) > nx / 2:
        return None

    return np.column_stack([px, py])


# --------------------------------------------------------------------------- #
#  Main function                                                               #
# --------------------------------------------------------------------------- #

def plot(skymap_path, event_files, output_path,
         hips_id=DEFAULT_HIPS, nx=1800, logger=None):
    """
    Produce a PNG showing all survey tiles overlaid on an HiPS background,
    with the input tiles highlighted.

    Uses a ZEA (zenithal equal-area) projection centred on the survey centroid
    so the hemisphere of eROSITA_DE sky fills a circle.

    Parameters
    ----------
    skymap_path   : str  - path to SMAPS tile catalogue FITS
    event_files   : list of str  -or-  path to a text file containing paths
    output_path   : str
    hips_id       : str  - HiPS identifier for the background image
    nx            : int  - image size in pixels (square, default 1800)
    logger        : optional logging.Logger
    """
    log = logger or logging.getLogger('plot_tiles')

    # ------------------------------------------------------------------ #
    #  1. Load SMAPS catalogue (needed for centroid before WCS)            #
    # ------------------------------------------------------------------ #

    with fits.open(skymap_path) as hdul:
        smaps = hdul['SMAPS'].data

    srvmap_all  = smaps['SRVMAP'].astype(int)
    ra_min_all  = smaps['RA_MIN']
    ra_max_all  = smaps['RA_MAX']
    dec_min_all = smaps['DE_MIN']
    dec_max_all = smaps['DE_MAX']
    ra_cen_all  = smaps['RA_CEN']
    dec_cen_all = smaps['DE_CEN']
    n_tiles     = len(srvmap_all)

    # ------------------------------------------------------------------ #
    #  2. Build ZEA WCS centred on survey, fetch HiPS background          #
    # ------------------------------------------------------------------ #
    glon_cen, glat_cen = _catalog_centroid(smaps)

    wcs   = _make_halfsky_wcs(nx, glon_cen, glat_cen)
    ny    = nx   # square image
    image = _fetch_hips(hips_id, wcs, log)

    # ------------------------------------------------------------------ #
    #  3. Resolve event-file list  →  used SKYFIELD IDs                   #
    # ------------------------------------------------------------------ #
    if isinstance(event_files, (str, os.PathLike)):
        ev_path = str(event_files)
        if os.path.isfile(ev_path):
            with open(ev_path) as f:
                event_files = [l.strip() for l in f if l.strip()]
        else:
            event_files = glob.glob(ev_path)

    event_files = list(event_files)
    # log.info(f'Reading SKYFIELD from {len(event_files)} raw event file(s) …')

    used_ids = set()
    for path in event_files:
        sf = _skyfield_from_file(path)
        if sf is not None:
            used_ids.add(sf)

    log.info(f'Unique tiles in input: {sorted(used_ids)}')

    # ------------------------------------------------------------------ #
    #  4. Build polygon arrays                                             #
    # ------------------------------------------------------------------ #

    bg_polys   = []   # all catalogue tiles
    used_polys = []   # tiles present in the input data

    for i in range(n_tiles):
        verts = _tile_polygon_pixels(
            ra_min_all[i], ra_max_all[i],
            dec_min_all[i], dec_max_all[i],
            wcs, nx, nx   # square image
        )
        if verts is None:
            continue

        if srvmap_all[i] in used_ids:
            used_polys.append(Polygon(verts, closed=True))
        else:
            bg_polys.append(Polygon(verts, closed=True))

    # ------------------------------------------------------------------ #
    #  5. Build zoom WCS centred on input tiles                          #
    # ------------------------------------------------------------------ #
    from reproject import reproject_interp as _reproj

    zoom_glons, zoom_glats = [], []
    for sf in used_ids:
        idx_sf = np.where(srvmap_all == sf)[0]
        if len(idx_sf) == 0:
            continue
        i = idx_sf[0]
        for ra, dec in [(ra_min_all[i], dec_min_all[i]),
                        (ra_min_all[i], dec_max_all[i]),
                        (ra_max_all[i], dec_min_all[i]),
                        (ra_max_all[i], dec_max_all[i])]:
            s = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs').galactic
            zoom_glons.append(s.l.deg)
            zoom_glats.append(s.b.deg)

    if zoom_glons:
        zoom_glons = np.array(zoom_glons)
        zoom_glats = np.array(zoom_glats)
        med_l = np.median(zoom_glons)
        zoom_glons[zoom_glons - med_l >  180] -= 360
        zoom_glons[zoom_glons - med_l < -180] += 360
        z_l_cen = float(np.mean(zoom_glons) % 360)
        z_b_cen = float(np.mean(zoom_glats))
        fov = max(5.0,
                  max(float(np.ptp(zoom_glons)),
                      float(np.ptp(zoom_glats))) * 1.5)
    else:
        z_l_cen, z_b_cen, fov = 0.0, 0.0, 20.0

    zoom_nx = nx
    zoom_cdelt = fov / zoom_nx
    zoom_wcs = WCS(naxis=2)
    zoom_wcs.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']
    zoom_wcs.wcs.crval = [z_l_cen, z_b_cen]
    zoom_wcs.wcs.crpix = [zoom_nx / 2.0 + 0.5, zoom_nx / 2.0 + 0.5]
    zoom_wcs.wcs.cdelt = [-zoom_cdelt, zoom_cdelt]
    zoom_wcs.pixel_shape = (zoom_nx, zoom_nx)

    zoom_image, _ = _reproj((image, wcs), zoom_wcs,
                            shape_out=(zoom_nx, zoom_nx))

    zoom_bg_polys   = []
    zoom_used_polys = []
    for i in range(n_tiles):
        verts = _tile_polygon_pixels(
            ra_min_all[i], ra_max_all[i],
            dec_min_all[i], dec_max_all[i],
            zoom_wcs, zoom_nx, zoom_nx
        )
        if verts is None:
            continue
        if srvmap_all[i] in used_ids:
            zoom_used_polys.append(Polygon(verts, closed=True))
        else:
            zoom_bg_polys.append(Polygon(verts, closed=True))

    # Tight pixel bounds of the input tiles in the zoom frame
    if zoom_used_polys:
        _all_verts = np.vstack([p.get_xy() for p in zoom_used_polys])
        _tile_span = max(np.ptp(_all_verts[:, 0]), np.ptp(_all_verts[:, 1]))
        _tile_margin = _tile_span * 0.15
        z_xlim = (float(_all_verts[:, 0].min() - _tile_margin),
                  float(_all_verts[:, 0].max() + _tile_margin))
        z_ylim = (float(_all_verts[:, 1].min() - _tile_margin),
                  float(_all_verts[:, 1].max() + _tile_margin))
    else:
        z_xlim = (0, zoom_nx)
        z_ylim = (0, zoom_nx)

    # Find all tile IDs whose centres fall inside the cropped view
    zoom_visible_ids = set()
    for i in range(n_tiles):
        _sky = SkyCoord(ra=ra_cen_all[i], dec=dec_cen_all[i],
                        unit='deg', frame='icrs')
        _gal = _sky.galactic
        try:
            _px, _py = zoom_wcs.all_world2pix(_gal.l.deg, _gal.b.deg, 0)
        except Exception:
            continue
        if (z_xlim[0] <= float(_px) <= z_xlim[1]
                and z_ylim[0] <= float(_py) <= z_ylim[1]):
            zoom_visible_ids.add(srvmap_all[i])
    # remove the input tiles — they are already labelled at full brightness
    zoom_visible_ids -= used_ids

    # ------------------------------------------------------------------ #
    #  6. Figure — side-by-side layout                                   #
    # ------------------------------------------------------------------ #

    fig = plt.figure(figsize=(20, 10), facecolor='black')
    ax   = fig.add_axes([0.01, 0.02, 0.57, 0.96], projection=wcs)
    ax_z = fig.add_axes([0.60, 0.02, 0.39, 0.96], projection=zoom_wcs)

    # Compute a shared norm from the full-sky image so both panels use the
    # same colour scale (vmin = 0.5th percentile, vmax = 99.5th percentile).
    _main_fin = image[np.isfinite(image)]
    if _main_fin.size > 1:
        _shared_norm = simple_norm(_main_fin, stretch='asinh',
                                   min_percent=0.5, max_percent=99.5)
    else:
        from matplotlib.colors import Normalize as _Norm
        _shared_norm = _Norm(vmin=0, vmax=1)

    # -- shared render helper --
    def _render_panel(panel_ax, panel_image, panel_wcs, panel_nx,
                      panel_bg_polys, panel_used_polys,
                      label_fontsize=11, lw_used=1.4, norm=None,
                      extra_label_ids=None, extra_label_alpha=0.45,
                      show_used_labels=True):
        img_disp = np.where(np.isfinite(panel_image), panel_image, np.nan)
        if norm is None:
            _fin = img_disp[np.isfinite(img_disp)]
            if _fin.size > 1:
                norm = simple_norm(_fin, stretch='asinh',
                                   min_percent=0.5, max_percent=99.5)
            else:
                from matplotlib.colors import Normalize as _Norm
                norm = _Norm(vmin=0, vmax=1)
        panel_ax.imshow(img_disp, origin='lower', cmap='inferno', norm=norm,
                        interpolation='nearest')
        if panel_bg_polys:
            panel_ax.add_collection(PatchCollection(
                panel_bg_polys, match_original=False,
                facecolor='none', edgecolor=(0.55, 0.55, 0.55),
                linewidth=0.3, alpha=0.85))
        if panel_used_polys:
            panel_ax.add_collection(PatchCollection(
                panel_used_polys, match_original=False,
                facecolor=(0.2, 0.6, 1.0), edgecolor='none', alpha=0.2))
            panel_ax.add_collection(PatchCollection(
                panel_used_polys, match_original=False,
                facecolor='none', edgecolor=(0.2, 0.8, 1.0),
                linewidth=lw_used))
        for sf in sorted(used_ids) if show_used_labels else []:
            _idx = np.where(srvmap_all == sf)[0]
            if len(_idx) == 0:
                continue
            _i = _idx[0]
            _sky = SkyCoord(ra=ra_cen_all[_i], dec=dec_cen_all[_i],
                            unit='deg', frame='icrs')
            _gal = _sky.galactic
            try:
                _px, _py = panel_wcs.all_world2pix(_gal.l.deg, _gal.b.deg, 0)
            except Exception:
                continue
            panel_ax.text(float(_px), float(_py), f'{sf:06d}',
                          color='white', fontsize=label_fontsize,
                          ha='center', va='center', fontweight='bold',
                          path_effects=[pe.withStroke(linewidth=1.5,
                                                      foreground='black')])
        # Background tile labels (dimmer)
        if extra_label_ids:
            for sf in sorted(extra_label_ids):
                _idx = np.where(srvmap_all == sf)[0]
                if len(_idx) == 0:
                    continue
                _i = _idx[0]
                _sky = SkyCoord(ra=ra_cen_all[_i], dec=dec_cen_all[_i],
                                unit='deg', frame='icrs')
                _gal = _sky.galactic
                try:
                    _px, _py = panel_wcs.all_world2pix(_gal.l.deg, _gal.b.deg, 0)
                except Exception:
                    continue
                panel_ax.text(float(_px), float(_py), f'{sf:06d}',
                              color='white', fontsize=label_fontsize,
                              ha='center', va='center',
                              alpha=extra_label_alpha,
                              path_effects=[pe.withStroke(linewidth=1.0,
                                                          foreground='black')])
        panel_ax.set_facecolor('black')
        _ctype0 = panel_wcs.wcs.ctype[0].split('-')[0].lower()
        _ctype1 = panel_wcs.wcs.ctype[1].split('-')[0].lower()
        for _c in (panel_ax.coords[_ctype0], panel_ax.coords[_ctype1]):
            _c.set_ticks_visible(False)
            _c.set_ticklabel_visible(False)
            _c.set_axislabel('')
            _c.grid(False)
        panel_ax.coords.frame.set_color('none')

    _render_panel(ax,   image,      wcs,      nx,      bg_polys,      used_polys,
                  label_fontsize=11, lw_used=1.4, norm=_shared_norm,
                  show_used_labels=False)
    _render_panel(ax_z, zoom_image, zoom_wcs, zoom_nx, zoom_bg_polys, zoom_used_polys,
                  label_fontsize=18,  lw_used=2.0, norm=_shared_norm,
                  extra_label_ids=zoom_visible_ids, extra_label_alpha=0.45)

    # Crop the zoom panel to just the input tile area
    ax_z.set_xlim(*z_xlim)
    ax_z.set_ylim(*z_ylim)

    # ------------------------------------------------------------------ #
    #  7. Galactic reference lines — main panel                          #
    # ------------------------------------------------------------------ #
    _ref_l = np.linspace(0.0, 360.0, 1800)
    _ref_b = np.linspace(-90.0, 90.0, 900)

    def _draw_lat_line(b_deg, lw, ls, alpha):
        _px, _py = wcs.all_world2pix(_ref_l, np.full_like(_ref_l, b_deg), 0)
        _v = (np.isfinite(_px) & np.isfinite(_py)
              & (_px >= 0) & (_px < nx) & (_py >= 0) & (_py < nx))
        if _v.sum() < 2:
            return
        for _seg in np.split(np.column_stack([_px, _py, _v]),
                              np.where(np.diff(_v.astype(int)) != 0)[0] + 1):
            _m = _seg[:, 2].astype(bool)
            if _m.sum() > 1:
                ax.plot(_seg[_m, 0], _seg[_m, 1],
                        color='white', lw=lw, alpha=alpha, ls=ls,
                        transform=ax.transData)
        try:
            _lx = float(wcs.all_world2pix([1.0], [float(b_deg)], 0)[0][0])
            _ly = float(wcs.all_world2pix([1.0], [float(b_deg)], 0)[1][0])
            if 0 <= _lx < nx and 0 <= _ly < nx:
                ax.text(_lx - 10, _ly + 6, fr'${b_deg:+g}$°',
                        color='white', fontsize=16, alpha=0.85,
                        ha='center', va='bottom')
        except Exception:
            pass

    def _draw_lon_line(l_deg, lw, ls, alpha):
        _px, _py = wcs.all_world2pix(np.full_like(_ref_b, l_deg), _ref_b, 0)
        _v = (np.isfinite(_px) & np.isfinite(_py)
              & (_px >= 0) & (_px < nx) & (_py >= 0) & (_py < nx))
        if _v.sum() < 2:
            return
        ax.plot(_px[_v], _py[_v], color='white', lw=lw, alpha=alpha, ls=ls,
                transform=ax.transData)
        _idx = np.where(_v)[0][_v.sum() // 2]
        ax.text(float(_px[_idx]) + 5, float(_py[_idx]) - 10,
                fr'${l_deg}$°', color='white', fontsize=16, alpha=1,
                ha='left', va='center')

    for _b in (30, 60, 0, -30, -60):
        _draw_lat_line(_b, lw=1, ls='--', alpha=0.80)
    for _l in (0, 180, 210, 240, 270, 300, 330):
        _draw_lon_line(_l, lw=1, ls=':', alpha=0.8)

    # ------------------------------------------------------------------ #
    #  8. Reference grid on zoom panel + yellow box on main panel        #
    # ------------------------------------------------------------------ #
    _z_ref_l = np.linspace(z_l_cen - fov * 2, z_l_cen + fov * 2, 600)
    _z_ref_b = np.linspace(z_b_cen - fov * 2, z_b_cen + fov * 2, 600)
    _grid_step = max(1, int(fov / 4))  # sensible grid spacing in degrees

    for _b in range(-90, 91, _grid_step):
        if abs(_b - z_b_cen) > fov:
            continue
        _px, _py = zoom_wcs.all_world2pix(
            _z_ref_l, np.full_like(_z_ref_l, float(_b)), 0)
        _v = (np.isfinite(_px) & np.isfinite(_py)
              & (_px >= z_xlim[0]) & (_px <= z_xlim[1])
              & (_py >= z_ylim[0]) & (_py <= z_ylim[1]))
        if _v.sum() > 1:
            ax_z.plot(_px[_v], _py[_v], color='white', lw=0.6,
                      alpha=0.45, ls='--', transform=ax_z.transData)
            # label at left edge of the cropped view
            try:
                _lx = float(z_xlim[0]) + 5
                _ly = float(np.interp(_lx, _px[_v], _py[_v]))
                if z_ylim[0] <= _ly <= z_ylim[1]:
                    ax_z.text(_lx, _ly, f'b={_b:+d}°',
                              color='white', fontsize=14, alpha=0.7,
                              ha='left', va='center')
            except Exception:
                pass

    for _l in range(0, 361, _grid_step):
        _l_n = float(_l)
        if abs(((_l_n - z_l_cen + 180) % 360) - 180) > fov:
            continue
        _px, _py = zoom_wcs.all_world2pix(
            np.full_like(_z_ref_b, _l_n), _z_ref_b, 0)
        _v = (np.isfinite(_px) & np.isfinite(_py)
              & (_px >= z_xlim[0]) & (_px <= z_xlim[1])
              & (_py >= z_ylim[0]) & (_py <= z_ylim[1]))
        if _v.sum() > 1:
            ax_z.plot(_px[_v], _py[_v], color='white', lw=0.6,
                      alpha=0.45, ls=':', transform=ax_z.transData)
            # label at bottom edge of the cropped view
            try:
                _ly = float(z_ylim[0]) + 5
                _lx = float(np.interp(_ly, _py[_v], _px[_v]))
                if z_xlim[0] <= _lx <= z_xlim[1]:
                    ax_z.text(_lx, _ly, f'l={int(_l_n)}°',
                              color='white', fontsize=14, alpha=0.7,
                              ha='center', va='bottom')
            except Exception:
                pass

    # Yellow box on main panel marking the cropped tile region
    _cx = [z_xlim[0], z_xlim[1], z_xlim[1], z_xlim[0], z_xlim[0]]
    _cy = [z_ylim[0], z_ylim[0], z_ylim[1], z_ylim[1], z_ylim[0]]
    try:
        _box_world = zoom_wcs.all_pix2world(np.column_stack([_cx, _cy]), 0)
        _bx, _by = wcs.all_world2pix(_box_world[:, 0], _box_world[:, 1], 0)
        ax.plot(_bx, _by, color='yellow', lw=1.2, alpha=0.85,
                transform=ax.transData)
    except Exception:
        pass

    # Yellow border around zoom panel
    for spine in ax_z.spines.values():
        spine.set_edgecolor('yellow')
        spine.set_linewidth(1.2)

    # ------------------------------------------------------------------ #
    #  9. Legend + titles                                                  #
    # ------------------------------------------------------------------ #
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=(0.55, 0.55, 0.55), lw=1,
               label=f'All survey tiles ({n_tiles})'),
        Line2D([0], [0], color=(0.2, 0.8, 1.0), lw=2,
               label=f'Input tiles ({len(used_ids)})')]
    ax.legend(handles=legend_elements, loc='lower right',
              facecolor='#111111', edgecolor='white',
              labelcolor='white', fontsize=16, framealpha=0.8)

    ax.set_title('eROSITA-DE DR1 Band 024 (0.2-2.3 keV) Count Rate Map',
                 color='white', fontsize=22, pad=8)
    ax_z.set_title('Input tile region (zoomed)',
                   color='white', fontsize=22, pad=8)

    # -- save --
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='black')
    plt.close(fig)
    log.info(f'\nTile map saved as {output_path}')


# --------------------------------------------------------------------------- #
#  CLI entry point                                                             #
# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Visualise eROSITA survey tiles on an HiPS background.')
    parser.add_argument('skymap',
                        help='Path to SMAPS tile catalogue (SKYMAPS_052022_MPE.fits)')
    parser.add_argument('event_files',
                        help='Glob pattern or text-file list of raw event files')
    parser.add_argument('output',
                        help='Output image path (.png)')
    parser.add_argument('--hips', default=DEFAULT_HIPS, metavar='HIPS_ID',
                        help=f'HiPS identifier for the background image '
                             f'(default: {DEFAULT_HIPS})')
    parser.add_argument('--nx', type=int, default=1800,
                        help='Image size in pixels — square (default: 1800)')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(message)s', stream=sys.stdout)

    ev = args.event_files if os.path.isfile(args.event_files) else glob.glob(args.event_files)
    plot(args.skymap, ev, args.output,
         hips_id=args.hips, nx=args.nx)
