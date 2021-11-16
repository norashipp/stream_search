#!/usr/bin/env python
"""
Generic python script.
"""
import os
from os.path import abspath, dirname, realpath
import subprocess
from collections import OrderedDict as odict
import copy

import matplotlib as mpl
import __main__ as main
# if not hasattr(main, '__file__'):
#     mpl.use('Agg')

import fitsio
import healpy as hp
import numpy as np
import pylab as plt
from matplotlib.colors import LogNorm
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u
import yaml
from matplotlib.patches import Ellipse


# import skymap.survey
# from skymap.survey import DESSkymapMcBryde
from skymap.survey import DESSkymap
# from skymap.utils import gal2cel, cel2gal
# from skymap.utils import setdefaults
# from skymap.healpix import ang2disc
from ugali.utils.shell import mkdir
from ugali.utils.projector import mod2dist, angsep, dist2mod


from numpy.polynomial import polynomial
from polyfit2d import polyfit2d
from results import create_results
import results
from streams import GLOBULARS, DWARFS, GALAXIES
import elysian

# import load_stream

import galstreams
import rotation_matrix

# DATADIR = os.path.join(dirname(realpath(__file__)), '../data')
DATADIR = os.path.join(dirname(realpath(__file__)), '/home/s1/nshipp/projects/stream_search/data')
STREAMFILE = os.path.join(DATADIR, 'streams_v6.0.yaml')


def load_streams(filename=STREAMFILE):
    return yaml.load(open(STREAMFILE, 'r'))

# STREAMS = load_streams()

FRACMIN = 0.5
SIGMA = 0.2
NSIDE = 256

# try:
#     POINTINGS = np.genfromtxt('../data/pointings/stream_pointing_maps_v7.csv', dtype=None, delimiter=',')
#     # POINTINGS = np.genfromtxt('../data/stream_pointing_maps_v2.csv', dtype=None, delimiter=',')
# except:
#     POINTINGS = np.genfromtxt('/Users/nora/projects/proper_motions/data/pointings/stream_pointing_maps_v7.csv', dtype=None, delimiter=',')
#     # POINTINGS = np.genfromtxt('/Users/nora/projects/proper_motions/data/stream_pointing_maps_v2.csv', dtype=None, delimiter=',')


class DESSkymapNGC288(DESSkymap):
    """Class for plotting a zoom on DES. This is relatively inflexible."""
    # [[RA1,RA2],[DEC1,DEC2]] frame limits
    FRAME = [[9, 18], [-22, -33]]
    FIGSIZE = [3.4, 5]

    def draw_inset_colorbar(self, *args, **kwargs):
        defaults = dict(loc=1, height="3%", width="30%",
                        bbox_to_anchor=(0, -0.05, 1, 1))
        setdefaults(kwargs, defaults)
        super(DESSkymap, self).draw_inset_colorbar(*args, **kwargs)

    def get_map_range(self, hpxmap, pixel=None, nside=None):
        """ Calculate the longitude and latitude range for an implicit map. """
        return FRAME


class DESSkymapNGC1261(DESSkymap):
    """Class for plotting a zoom on DES. This is relatively inflexible."""
    # RA, DEC frame limits
    FRAME = [[47, 48], [-60, -50]]
    FIGSIZE = [3.4, 5]

    def draw_inset_colorbar(self, *args, **kwargs):
        defaults = dict(loc=3, height="4%", width="20%",
                        bbox_to_anchor=(0, 0.05, 1, 1))
        setdefaults(kwargs, defaults)
        super(DESSkymap, self).draw_inset_colorbar(*args, **kwargs)


class DESSkymapNGC1851(DESSkymap):
    """Class for plotting a zoom on DES. This is relatively inflexible."""
    # RA, DEC frame limits
    #FRAME = [[80,80],[-50,-30]]
    FRAME = [[80, 80], [-47, -33]]
    FIGSIZE = [3.4, 5]

    def draw_inset_colorbar(self, *args, **kwargs):
        defaults = dict(loc=4, height="5%", width="25%",
                        bbox_to_anchor=(0, 0.05, 1, 1))
        setdefaults(kwargs, defaults)
        super(DESSkymap, self).draw_inset_colorbar(*args, **kwargs)


class DESSkymapNGC1904(DESSkymap):
    """Class for plotting a zoom on DES. This is relatively inflexible."""
    # RA, DEC frame limits
    FRAME = [[80.5, 82.5], [-30, -21]]
    FIGSIZE = [3.4, 5]

    def draw_inset_colorbar(self, *args, **kwargs):
        defaults = dict(loc=4, height="5%", width="25%",
                        bbox_to_anchor=(0, 0.05, 1, 1))
        setdefaults(kwargs, defaults)
        super(DESSkymap, self).draw_inset_colorbar(*args, **kwargs)


def skymap_factory(proj):
    proj = proj.lower()
    if proj == 'mbtpq':
        smap = skymap.survey.DESSkymap
    elif proj == 'q1':
        smap = skymap.survey.DESSkymapQ1
    elif proj == 'q2':
        smap = skymap.survey.DESSkymapQ2
    elif proj == 'q3':
        smap = skymap.survey.DESSkymapQ3
    elif proj == 'q4':
        smap = skymap.survey.DESSkymapQ4
    elif proj == 'ngc288':
        smap = DESSkymapNGC288
    elif proj == 'ngc1261':
        smap = DESSkymapNGC1261
    elif proj == 'ngc1851':
        smap = DESSkymapNGC1851
    elif proj == 'ngc1904':
        smap = DESSkymapNGC1904
    elif proj == 'splaea':
        smap = skymap.survey.DESPolarLambert
    elif proj == 'laea':
        smap = skymap.survey.DESLambert
    elif proj == 'cyl':
        smap = skymap.survey.DESSkymapCart
    else:
        smap = skymap.survey.SurveySkymap

    return smap

DATA = odict([
    ('data', dict()),
    ('bkg', dict()),
    ('resid', dict()),
    ('frac', dict()),
    #('sfd',dict()),
])


CBAR_KWARGS = odict([
    ('mbtpq', dict(label='Counts/pixel')),
    ('splaea', dict(loc=7, height='3%', label='Counts/pixel')),
    ('laea', dict(loc=7, bbox_to_anchor=(-0.01, 0.07, 1, 1), label='Counts/healpix')),
    ('q1', dict(loc=3, bbox_to_anchor=(0, 0.03, 1, 1), label='Counts/healpix')),
    ('q2', dict(loc=2, bbox_to_anchor=(0, -0.1, 1, 1), label='Counts/healpix')),
    ('q3', dict(loc=3, bbox_to_anchor=(0, 0.12, 1, 1), label='Counts/healpix')),
    ('q4', dict(loc=4, bbox_to_anchor=(0, 0.05, 1, 1), label='Counts/healpix')),
    ('default', dict(label='Counts/healpix')),
])


def mask_lmc(nside=256):
    from skymap.constants import RA_LMC, DEC_LMC, RADIUS_LMC
    RADIUS_LMC = 5.7
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    mask[ang2disc(nside, RA_LMC, DEC_LMC, 3 * RADIUS_LMC)] = True
    return mask


def mask_milky_way_north_3(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    ra = ra - 360 * (ra > 180)
    vertices = [[-50, 90], [-130, 30], [-50, 30]]
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_milky_way_north_2(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    ra = ra - 360 * (ra > 180)
    vertices = [[-12, 36], [-22, -20], [-60, -15], [-40, 27]]
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_milky_way_north(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    ra = ra - 360 * (ra > 180)
    vertices = [[-125, -3], [-110, 35], [-88, 35], [-103, -1], [-125, -3]]
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_milky_way(nside=256):
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    #mask[(ra < 304) & (ra > 295) & (dec > -60) & (dec < -40)] = True
    mask[(ra < 315) & (ra > 295) & (dec > -70) & (dec < -35)] = True
    return mask


def mask_acs(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    # vertices = [[121, 60], [130, 60], [-70 + 360, 63], [-77 + 360, 55]]
    # vertices = [[130, 60], [135, 65], [155, 75], [-140+360, 80], [-100+360, 72], [-90+360, 65], [-80+360, 60], [-75+360, 68], [-85+360, 78], [-150+360, 86], [140, 80], [125, 75], [125, 70], [120, 60], [130, 60]]
    # vertices = [[-81.20842981166365+360, 60.02847479426181], [-103.18991410673861+360, 71.44587167686264], [-139.29088384800647+360, 77.70003482326835], [162.456378757321, 74.84730497914623], [133.21261597970803, 61.14200024319528], [119.66962482115947, 60.05509696491776], [120.47322378748359, 70.6530402996451], [159.60586655370196, 83.52806532949364], [-83.9730393866897+360, 77.83186299691168], [-73.65549049748984+360, 67.61016199435278]]
    vertices = np.array([[129.0948921713502, 62.52810768283558], [133.63357495373106, 66.79929543092777], [139.99127285541002, 71.18224005454947], [150.56539516173876, 75.29338855683497], [168.56101032015746, 79.16443785287761], [-163.2546873976258, 80.94741486188961], [-135.54629935256855, 79.87729083771019], [-111.91899764770008, 76.26405863933468], [-97.34702120529064, 71.62669968091069], [-93.82267929110368, 67.2709685408379], [-86.99439599698219, 63.61224297643058], [-78.6869084100243, 68.68932790216034], [-80.80161083732156, 72.72211970974755],
                         [-85.2655169949618, 75.79296094818497], [-92.50988663530084, 80.38692902617865], [-102.86187563750477, 82.74627019631558], [-139.84688477685475, 85.9680487649359], [177.27997933346577, 86.06424242318538], [158.10739160942848, 84.2624478173823], [139.38461451738698, 81.55584566801063], [137.38898697973735, 79.79972916219887], [132.22455296516068, 77.23410025797286], [130.36894206992676, 73.82593566527407], [128.6008401703078, 69.74585710631004], [124.44422835367168, 65.38883293725237], [124.6106442460963, 63.08164363064252]])
    vertices[:, 0][vertices[:, 0] < 0] += 360
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_sgr_north(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    # ra = ra - 360 * (ra > 180)
    # vertices = [[120, 15], [120, 23], [-125+360, -3], [-130+360, -10]]
    vertices = [[120, 15], [120, 26], [-125 + 360, 0], [-130 + 360, -10]]
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_sgr(nside=256):
    from matplotlib.path import Path
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    ra = ra - 360 * (ra > 180)
    vertices = [(0, -30), (47, -3), (47, +7), (36, +7), (5, -11)]
    mask[Path(vertices).contains_points(np.array([ra, dec]).T)] = True
    return mask


def mask_globulars(nside=256):
    """ Mask sources """
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    for key, val in GLOBULARS.items():
        #radius = 10*val['rh']/60.0
        if 'jacobi' in val:
            radius = np.degrees(
                np.arctan(1e-3 * val['jacobi'] / mod2dist(val['mod'])))
        else:
            radius = 0.5
        #print(key, radius)
        pix = ang2disc(nside, val['ra'], val['dec'], radius, inclusive=True)
        mask[pix] = True
    return mask


def mask_dwarfs(nside=256):
    """ Mask sources """
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    for key, val in DWARFS.items():
        radius = 5 * val['rh'] / 60.0
        #print(key, radius)
        pix = ang2disc(nside, val['ra'], val['dec'], radius, inclusive=True)
        mask[pix] = True
    return mask


def mask_more_globs_dwarfs(nside=256):
    # File from Adrian - McConnachie+2012, Simon+2018, DES papers for the dwarfs, and Vasiliev+2018 for the globular clusters
    """ Mask sources """
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)

    globs_dwarfs = fitsio.read('../data/MW_dwarfs_globs.fits')

    for gd in globs_dwarfs:
        radius = 5 * gd['r_h'] / 60.0
        if radius >= 2.0:
            radius = 2.0
        pix = ang2disc(nside, gd['ra'], gd['dec'], radius, inclusive=True)
        mask[pix] = True
    return mask


def mask_galaxies(nside=256):
    """ Mask sources """
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    for key, val in GALAXIES.items():
        radius = 2 * val['rh'] / 60.0
        #print(key, radius)
        pix = ang2disc(nside, val['ra'], val['dec'], radius, inclusive=True)
        mask[pix] = True
    return mask


def mask_plane(nside=256, bmax=25):
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    b = c.galactic.b.deg
    l = c.galactic.l.deg
    # | ((np.abs(b) < bmax + 20) & (np.abs(l) < 60))
    mask[(np.abs(b) < bmax)] = True
    return mask


def mask_center(nside=256, cmax=25):
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    sep = angsep(0, 0, c.galactic.l.deg, c.galactic.b.deg)
    mask[(sep < cmax)] = True
    return mask


# lmc=False, milky_way=False, sgr=False, globulars=False, dwarfs=False, galaxies=False, plane=False, bmax=25, center=False, globs_dwarfs=False, acs=False

def make_mask(nside=256, masking=True, lmc=True, milky_way=True, sgr=False, globulars=True, dwarfs=True, galaxies=True, plane=True, bmax=25, center=True, globs_dwarfs=True, acs=False, cmax=25):
    mask = np.zeros(hp.nside2npix(nside), dtype=bool)
    if masking:
        if lmc:
            mask |= mask_lmc(nside)
        if milky_way:
            mask |= mask_milky_way(nside)
            mask |= mask_milky_way_north(nside)
            mask |= mask_milky_way_north_2(nside)
            mask |= mask_milky_way_north_3(nside)
        if sgr:
            mask |= mask_sgr(nside)
            mask |= mask_sgr_north(nside)
        if globulars:
            mask |= mask_globulars(nside)
        if dwarfs:
            mask |= mask_dwarfs(nside)
        if galaxies:
            mask |= mask_galaxies(nside)
        if galaxies:
            mask |= mask_galaxies(nside)
        if plane:
            mask |= mask_plane(nside, bmax)
        if center:
            mask |= mask_center(nside, cmax)
        if globs_dwarfs:
            mask |= mask_more_globs_dwarfs(nside)
        if acs:
            mask |= mask_acs(nside)
    return mask


def degrade(nside, hpxmap):
    print("Degrading healpix map...")
    orig_nside = hp.get_nside(hpxmap)
    if nside > orig_nside:
        msg = "nside too large"
        raise Exception(msg)
    spix = superpixel(np.arange(len(hpxmap)), orig_nside, nside)
    data = np.zeros(hp.nside2npix(nside))
    np.add.at(data, spix, hpxmap)
    return data


def prepare_data(hpxmap, fracdet, fracmin=FRACMIN, clip=None, degrade=None, mask_kw=dict()):
    nside = hp.get_nside(hpxmap)
    mask = (fracdet < fracmin) | (hpxmap > np.percentile(hpxmap, clip)) | make_mask(nside, **mask_kw)
    data = np.ma.array(hpxmap, mask=mask, fill_value=np.nan)
    data /= fracdet
    return data


def get_rotmat(stream=None, ends=None, center=None):
    if stream is not None:
        print(stream)
        mw_streams = galstreams.MWStreams(verbose=False)
        ends = [(mw_streams[stream].end_f.ra.value, mw_streams[stream].end_f.dec.value),
                (mw_streams[stream].end_o.ra.value, mw_streams[stream].end_o.dec.value)]
    elif ends is not None:
        pass
    else:
        print('Need stream name or ends!')

    phi, theta, psi = results.euler_angles(
        ends[0][0], ends[0][1], ends[1][0], ends[1][1], center)

    print(phi, theta, psi)

    R = results.create_matrix(phi, theta, psi)
    return R


def get_streampix(data, stream=None, ends=None, R=None):
    if R is None:
        R = get_rotmat(stream=stream, ends=ends)
    nside = hp.get_nside(data)
    lon, lat = hp.pix2ang(nside, np.arange(len(data)), lonlat=True)
    streampix = hp.ang2pix(
        nside, *rotation_matrix.phi12_rotmat(lon, lat, np.linalg.inv(R)), lonlat=True)
    return streampix


def fit_bkg_poly(data, center=(0, 0), coords='cel', coord_stream=None, proj='cyl', sigma=0.1, percent=[2, 95], deg=5):
    """ Fit foreground/background with a polynomial """
    nside = hp.get_nside(data.mask)
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)

    vmin, vmax = np.percentile(data.compressed(), q=percent)
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)
    data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False),
                       mask=data.mask)

    smap = skymap.Skymap(parallels=False, meridians=False,
                         lon_0=center[0], lat_0=center[1], projection=proj)

    if coords == 'gal':
        galpix = hp.ang2pix(nside, *gal2cel(lon, lat), lonlat=True)
        data = data[galpix]
    elif coords == 'stream':
        if not coord_stream:
            print('Need to input coord_stream!')
        else:
            print('Converting to %s coords.' % coord_stream)
        streampix = get_streampix(data=data, stream=coord_stream)
        data = data[streampix]

    sel = ~data.mask
    x, y = smap(lon[sel], lat[sel])
    v = data[sel]
    c = polyfit2d(x, y, v, [deg, deg])

    # Evaluate the polynomial
    x, y = smap(lon, lat)
    bkg = polynomial.polyval2d(x, y, c)
    bkg = np.ma.array(bkg, mask=data.mask, fill_value=np.nan)
    return bkg


def fit_bkg_gauss(data, sigma=2.0, percent=[2, 95]):
    """ Fit foreground/background counts with using a Gaussian smoothing """
    nside = hp.get_nside(data.mask)
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)

    vmin, vmax = np.percentile(data.compressed(), q=percent)
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)
    bkg = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False),
                      mask=data.mask)

    return bkg


def fit_bkg(data, proj='mbtpq', sigma=SIGMA):
    # plt.figure()
    smap = skymap_factory(proj)()
    # plt.close()

    nside = hp.get_nside(data.mask)

    vmin, vmax = np.percentile(data[~data.mask], q=[2, 95])
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)
    data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False),
                       mask=data.mask)

    # Fit the polynomial
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    sel = ~data.mask
    x, y = smap(lon[sel], lat[sel])
    v = data[sel]
    c = polyfit2d(x, y, v, [5, 5])

    # Evaluate the polynomial
    bkg = np.zeros_like(data)
    x, y = smap(lon, lat)
    bkg = polynomial.polyval2d(x, y, c)
    bkg = np.ma.array(bkg, mask=data.mask)

    return bkg


def fit_data_bkg(hpxmap, fracdet, fracmin=FRACMIN, sigma=SIGMA, percent=[2, 95]):

    nside = hp.get_nside(hpxmap)
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    mask = (fracdet < fracmin) | make_mask(nside)
    data = np.ma.array(hpxmap, mask=mask)
    data /= fracdet
    vmin, vmax = np.percentile(data.compressed(), q=percent)
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)
    data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False),
                       mask=data.mask)

    # Background (don't draw)
    smap = skymap.McBrydeSkymap(parallels=False, meridians=False)

    sel = ~data.mask
    x, y = smap(lon[sel], lat[sel])
    v = data[sel]
    c = polyfit2d(x, y, v, [5, 5])

    # Evaluate the polynomial
    x, y = smap(lon, lat)
    bkg = polynomial.polyval2d(x, y, c)
    bkg = np.ma.array(bkg, mask=data.mask, fill_value=np.nan)

    return data, bkg


def draw_arc(smap, points, coord='cel', offset=0, **kwargs):
    import elysian

    if len(points) < 2:
        raise Exception("Must specify at least two points.")

    stream = elysian.factory(dict(ends=points))

    lons, lats = np.array(points).T
    if coord == 'gal':
        lons, lats = cel2gal(lons, lats)

    lon, lat = smap.great_circle(lons[0], lats[0], lons[1], lats[1], 'short')
    lon, lat = apply_offset(lon, lat, offset)
    smap.plot(*smap(lon, lat), **kwargs)

    lonmax, latmax = lon[np.argmax(lat)], lat[np.argmax(lat)]

    # for i in np.arange(len(points)-1):
    #    lon = lons[i:i+2]
    #    lat = lats[i:i+2]
    #
    #    perp = np.array([-(lat[1]-lat[0]),lon[1]-lon[0]])
    #    perp /= np.sqrt( (perp**2).sum() )
    #
    #    #lon += (offset/np.cos(np.radians(lat)))*perp[0]
    #    lon += offset*perp[0]
    #    lat += offset*perp[1]
    #
    #    if lat.max() > latmax:
    #        lonmax,latmax = lon[np.argmax(lat)],lat[lat.argmax()]
    #
    #    smap.draw_great_circle(lon[0],lat[0],lon[1],lat[1],'short',**kwargs)
    #    #smap.plot(*smap(lon,lat),**kwargs)

    return lonmax, latmax

# def apply_offset(lons,lats,offset):
#    out_lon = []
#    out_lat = []
#    for i in np.arange(len(lons)-1):
#        lon0,lon1 = lons[i:i+2]
#        lat0,lat1 = lats[i:i+2]
#
#        perp = np.array([(lat1-lat0),lon1-lon0])
#        perp /= np.sqrt((perp**2).sum() )
#        #perp[0] /= np.cos(np.mean([lat0,lat1])
#
#        print [lon0,lat0],[lon1,lat1]
#        print perp
#
#        #lon += (offset/np.cos(np.radians(lat)))*perp[0]
#        if i == 0:
#            out_lon.append( offset*perp[0]+lon0)
#            out_lat.append( offset*perp[1]+lat0)
#
#        out_lon.append( offset*perp[0]+lon1)
#        out_lat.append( offset*perp[1]+lat1)
#        print out_lon[-1],out_lat[-1]
#        print
#
#    return np.array(out_lon),np.array(out_lat)


def fit_polynomial(points, deg=2, ivar='lon', **kwargs):
    lons, lats = np.array(points).T

    if ivar == 'lon':
        dvar = 'lat'
        coeff = np.polyfit(lons, lats, deg, **kwargs)
        poly = np.poly1d(coeff)
    else:
        dvar = 'lon'
        coeff = np.polyfit(lats, lons, deg, **kwargs)
        poly = np.poly1d(coeff)

    polystr = '%s = ' % ivar
    for i, c in enumerate(coeff[::-1]):
        polystr += ' %+6.3f*%s^%i ' % (c, dvar, i)
    print(polystr)
    return poly


def apply_offset(lons, lats, offset):
    out_lon = []
    out_lat = []
    for i in np.arange(len(lons) - 1):
        lon0, lon1 = lons[i:i + 2]
        lat0, lat1 = lats[i:i + 2]
        phi, theta, psi = results.euler_angles(lon0, lat0, lon1, lat1)
        matrix = results.create_matrix(phi, theta, psi)
        L, B = elysian.sky_to_stream([lon0, lon1], [lat0, lat1], matrix)
        lon, lat = elysian.stream_to_sky(L, B + offset, matrix)
        if i == 0:
            out_lon.append(lon[0])
            out_lat.append(lat[0])

        out_lon.append(lon[1])
        out_lat.append(lat[1])
    return np.array(out_lon), np.array(out_lat)


def draw_stream_poly(smap, stream, coord='cel', offset=0, deg=2, ivar='lat', **kwargs):
    """ Draw a polynomial fit to Palca """
    defaults = dict(lw=2)
    setdefaults(kwargs, defaults)

    points = np.array(stream['ridge'])
    poly = fit_polynomial(points, deg, ivar=ivar)

    lat = np.linspace(points[:, 1].min(), points[:, 1].max(), 50)
    lon = poly(lat)
    lon, lat = apply_offset(lon, lat, offset)
    if coord == 'gal':
        lon, lat = cel2gal(lon, lat)

    smap.plot(*smap(lon, lat), **kwargs)
    return lon[np.argmax(lat)], lat[np.argmax(lat)]

    # smap.scatter(*smap(lon,lat),**kwargs)

# def draw_arc2(smap,stream,coord='cel',offset=0,**kwargs):
#    s = elysian.factory(stream)()
#    lon,lat = s.get_arc(s.ends[0][0],s.ends[0][1],s.ends[1][0],s.ends[1][1],
#                        offset=-offset)
#    if coord == 'gal':
#        lon,lat = cel2gal(lon,lat)
#    smap.plot(*smap(lon,lat),**kwargs)
#
#    return lon[np.argmax(lat)],lat[np.argmax(lat)]


def draw_stream_arc(smap, stream, coord='cel', offset=0, **kwargs):
    ends = np.array(stream['ends'])
    #lons,lats = apply_offset(ends[:,0],ends[:,1],offset)
    #lon,lat = smap.great_circle(lons[0],lats[0],lons[1],lats[1],'short')
    lon, lat = smap.great_circle(ends[0][0], ends[0][1], ends[1][0], ends[1][1],
                                 'short')
    lon, lat = apply_offset(lon, lat, offset)
    if coord == 'gal':
        lon, lat = cel2gal(lon, lat)
    smap.plot(*smap(lon, lat), **kwargs)

    return lon[np.argmax(lat)], lat[np.argmax(lat)]


def get_stream_label_coords(smap, stream, coord='cel', offset=0, **kwargs):
    ends = np.array(stream['ends'])
    lon, lat = smap.great_circle(ends[0][0], ends[0][1], ends[
                                 1][0], ends[1][1], 'short')
    lon, lat = apply_offset(lon, lat, offset)
    if coord == 'gal':
        lon, lat = cel2gal(lon, lat)
    return lon[np.argmax(lat)], lat[np.argmax(lat)]


def label_pointing(smap, text, lon, lat, **kwargs):
    defaults = dict(xycoords='data', fontsize=12, color='r', weight='black',
                    ha='right', va='bottom')

    if mpl.rcParams.get('text.usetex'):
        text = r"\textbf{%s}" % text.replace(' ', '\ ')

    setdefaults(kwargs, defaults)

    plt.gca().annotate(text, smap(lon - 0.8, lat - 0.8), **kwargs)


def draw_stream_ellipse(smap, stream, coord='cel', radius=1.0, offset=0., label=True, **kwargs):
    ax = plt.gca()

    pointings = np.copy(POINTINGS)
    names = pointings['f0']
    pointings = pointings[names == stream['name'].replace(' ', '_')]

    # if stream['name'] == 'Tucana III':
    #     print pointings

    lonmax, latmax = get_stream_label_coords(smap, stream, coord, offset=1.0)

    if len(pointings) == 0:
        return lonmax, latmax

    ra, dec = pointings['f2'], pointings['f3']
    x, y = smap(ra, dec)
    ra_label, dec_label = apply_offset(ra, dec, offset=offset)
    # ra_label, dec_label = ra, dec
    widths = np.abs(smap(ra + radius, dec)
                    [0] - smap(ra - radius, dec)[0]) / np.cos(np.deg2rad(dec))
    heights = np.abs((smap(ra, dec + radius)[1] - smap(ra, dec - radius)[1]))
    # ells = [Ellipse((x[i], y[i]), width=widths[i], height=heights[i], angle=0, fc=None, ec='red', lw=1, fill=False, ls='--') for i in range(len(ra))]
    try:
        observed = pointings['f6']
        print(observed)
    except:
        observed = np.ones_like(ra)

    # ells = [Ellipse((x[i], y[i]), width=widths[i], height=heights[i], angle=0, fc='cornflowerblue', ec='navy', lw=0.5, fill=observed[i], ls='-') for i in range(len(ra))]
    ells = []
    for i in range(len(ra)):
        if observed[i] == 2:
            continue
        elif observed[i] == 3:
            observed[i] = 1
            color = 'tomato'
            ecolor = 'maroon'
        else:
            color = 'cornflowerblue'
            ecolor = 'navy'
        ell = Ellipse((x[i], y[i]), width=widths[i], height=heights[
                      i], angle=0, fc=color, ec=ecolor, lw=0.5, fill=observed[i], ls='-')
        ells.append(ell)

    for i, e in enumerate(ells):
        ax.add_artist(e)
        if label:
            label_pointing(smap, str(pointings['f1'][i]), ra_label[
                           i], dec_label[i], fontsize=10, color='orange')

    return lonmax, latmax


def draw_streams(smap, mod=None, coord='cel', quad=None, streams=None, offset=0,
                 label=True, label_kw={}, pointings=False, label_pointings=True, pointing_radius=1.0, **kwargs):
    """ Draw streams and labels. """
    orig_kwargs = copy.deepcopy(kwargs)
    for key, values in STREAMS.items():
        kwargs = copy.deepcopy(orig_kwargs)

        setdefaults(kwargs, values['kwargs'])
        setdefaults(kwargs, dict(lw=2, ls='-'))

        name = values.get('name', 'key')
        if (mod is not None) and (values['mods'][0] > mod) or (values['mods'][1] < mod):
            continue
        quads = values['quadrants']
        squads = ['q%i' for i in quads]
        if (quad is not None) and (quad not in quads) and (quad not in squads):
            continue
        if (streams is not None) and (name not in streams):
            continue

        #points = np.array(values.get('ridge',values['ends']))
        print(name)

        if isinstance(offset, bool) and offset:
            off = values.get('offset', 1.0)
        else:
            off = offset
        off = np.atleast_1d(off)

        #lonmax,latmax = draw_arc(smap,points,coord=coord,offset=off,**kwargs)
        if pointings:
            lonmax, latmax = draw_stream_ellipse(
                smap, values, coord=coord, radius=pointing_radius, label=label_pointings)
        else:
            for o in off:
                if name in ['Palca', 'ATLAS']:
                    lonmax, latmax = draw_stream_poly(
                        smap, values, coord=coord, offset=o, **kwargs)
                else:
                    lonmax, latmax = draw_stream_arc(
                        smap, values, coord=coord, offset=o, **kwargs)

        kw = copy.deepcopy(label_kw)
        defaults = dict(lon=lonmax, lat=latmax, color=kwargs['color'])
        if coord == 'gal':
            setdefaults(kw, values.get('gal_text', {}))
        else:
            setdefaults(kw, values.get('cel_text', {}))
        setdefaults(kw, defaults)
        if label:
            label_stream(smap, name, **kw)

        # smap.draw_great_circle(lon[0],lat[0],lon[1],lat[1],**kwargs)


def label_stream(smap, text, lon, lat, **kwargs):
    defaults = dict(xycoords='data', fontsize=12, weight='black',
                    ha='right', va='bottom')

    if mpl.rcParams.get('text.usetex'):
        text = r"\textbf{%s}" % text.replace(' ', '\ ')
    #    text = r"{\bf %s}"%(text.replace(' ','\ '))
    #    #text = r"%s"%text.replace(' ','\ ')
    #    #defaults['fontsize'] = 12
    # else:
    #    pass
    setdefaults(kwargs, defaults)
    plt.gca().annotate(text, smap(lon, lat), **kwargs)


def draw_globulars(smap, mod, coord='cel'):
    """ Draw streams and labels. """
    for name, values in GLOBULARS.items():
        kwargs = values['kwargs']
        defaults = dict(marker='o', mec='none')
        setdefaults(kwargs, defaults)
        if (values['mods'][0] > mod) or (values['mods'][1] < mod):
            continue

        lon, lat = values['ra'], values['dec']
        if coord == 'gal':
            lon, lat = cel2gal(lon, lat)

        smap.plot(*smap(lon, lat), **kwargs)
        plt.gca().annotate(name, smap(lon, lat), xycoords='data',
                           xytext=(-3, 3), textcoords='offset points',
                           color=kwargs['color'], weight='bold', fontsize=10,
                           ha='right', va='bottom'
                           )


def draw_ngc1851_pm(smap, **kwargs):
    defaults = dict(ls='-', color='chartreuse')
    setdefaults(kwargs, defaults)
    ra, dec = 78.52817, -40.046
    pmra = 1.23 / np.cos(np.radians(dec))
    pmdec = 2.51
    pmra = -5.26
    pmdec = 0.54
    plt.plot(*smap([ra, ra + pmra], [dec, dec + pmdec]), **kwargs)


def draw_grillmair17(smap, coord='cel', offset=1.0, **kwargs):
    # Interpolate in ra from dec
    molonglo = [[345.017, -0.5843, +0.0182],
                'blue', np.linspace(-24, -12, 50)]
    murrumbidgee = [[+367.893, -0.4647, -0.00862, +0.000118, +1.2347e-6, -1.13758e-7],
                    'limegreen', np.linspace(-30, 38, 100)]
    # Interpolate in dec from ra
    orinoco = [[-25.5146, +0.1672, -0.003827, -0.0002835, -5.3133e-6],
               'red', np.linspace(-40, 20, 50)]
    kwando = [[-7.817, -2.354, +0.1202, -0.00215],
              'cyan', np.linspace(19, 33, 50)]

    for stream, c, dec in [molonglo, murrumbidgee]:
        poly = np.poly1d(stream[::-1])
        ra = poly(dec)
        if coord == 'gal':
            lon, lat = cel2gal(ra, dec)
        else:
            lon, lat = ra, dec
        smap.plot(*smap(*apply_offset(lon, lat, offset)), color=c, **kwargs)
        smap.plot(*smap(*apply_offset(lon, lat, -offset)), color=c, **kwargs)

    for stream, c, ra in [orinoco, kwando]:
        poly = np.poly1d(stream[::-1])
        dec = poly(ra)
        if coord == 'gal':
            lon, lat = cel2gal(ra, dec)
        else:
            lon, lat = ra, dec
        smap.plot(*smap(*apply_offset(lon, lat, offset)), color=c, **kwargs)
        smap.plot(*smap(*apply_offset(lon, lat, -offset)), color=c, **kwargs)


def hollywood(infiles, outfile=None, delay=40, queue='local'):
    print("Lights, Camera, Action...")
    infiles = np.atleast_1d(infiles)
    if not len(infiles):
        msg = "No input files found"
        raise ValueError(msg)

    infiles = ' '.join(infiles)
    if not outfile:
        outfile = infiles[0].replace('.png', '.gif')
    cmd = 'convert -delay %i -quality 100 %s %s' % (delay, infiles, outfile)
    if queue != 'local':
        cmd = 'csub -q %s ' % (queue) + cmd
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def quadrant_boundary(quadrant=1, steps=100):
    """ Return a quadrant boundaries in RA, DEC """
    CORNERS = [
        [(10, -68), (-71.3, -68), (-45, -38), (6, -38)],
        [(73.67, -45), (0, -45), (0, 8), (60, 8)],
        [(94.3, -68), (5, -68), (3.1, -38), (60, -38)],
        [(125.5, -55), (70, -55), (50, -15), (90, -15)]
    ]
    plt.figure()
    smap = skymap.survey.DESSkymap()
    plt.close()
    corners = CORNERS[quadrant - 1]
    ll = smap(*corners[0])
    lr = smap(*corners[1])
    ur = smap(*corners[2])
    ul = smap(*corners[3])

    left = smap(ll[0] * np.ones(steps),
                np.linspace(ul[1], ll[1], steps), inverse=True)
    bottom = smap(np.linspace(ll[0], lr[0], steps), ll[
                  1] * np.ones(steps), inverse=True)
    right = smap(lr[0] * np.ones(steps),
                 np.linspace(lr[1], ur[1], steps), inverse=True)
    top = smap(np.linspace(ur[0], ul[0], steps), ul[
               1] * np.ones(steps), inverse=True)

    x = np.hstack([left[0], bottom[0], right[0], top[0]])
    y = np.hstack([left[1], bottom[1], right[1], top[1]])
    return x, y


def draw_quadrant(smap, quadrant=1, coord='cel', **kwargs):
    """ Draw a quadrant boundary """
    lon, lat = quadrant_boundary(quadrant)
    if coord == 'gal':
        lon, lat = cel2gal(lon, lat)
    return smap.plot(*smap.proj(lon, lat), **kwargs)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
