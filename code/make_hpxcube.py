import os
import glob
import logging
from collections import OrderedDict as odict
import gc

import fitsio
import numpy as np
import healpy as hp
# import pylab as plt
import scipy.ndimage as nd
from matplotlib.path import Path

from astropy import units as u
from astropy.coordinates import SkyCoord

from utils import load_infiles
from ugali.utils import healpix
from ugali.analysis.isochrone import factory as isochrone_factory

from surveys import surveys
import filter_data


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', '--filename', default='iso_hpxcube.fits.gz')
    parser.add_argument('-n', '--nside', default=512)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-s', '--survey', default='DES_Y3')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format='%(message)s', level=level)
    nside = args.nside

    survey = args.survey
    print('Filtering %s...' %survey)
    mag = surveys[survey]['mag']
    mag_g = mag % 'G'
    mag_r = mag % 'R'
    mag_i = mag % 'I'
    ext = surveys[survey]['ext']
    if ext is not None:
        ext = surveys[survey]['ext']
        ext_g = ext % 'G'
        ext_r = ext % 'R'
        ext_i = ext % 'I'
        # fix this, not relevant for now
    minmag = surveys[survey]['minmag']
    maxmag = surveys[survey]['maxmag']
    # columns = ['RA', 'DEC', mag_g, mag_r, mag_i] # include i eventually, try mp search
    columns = ['RA', 'DEC', mag_g, mag_r]
    if stargal is not None:
        columns.append(surveys[survey]['stargal'])    

    ###################
    dmu = 0.1
    # modulii = np.arange(15, 20 + dmu, dmu)
    modulii = np.arange(surveys[survey]['modulii'][0], surveys[survey]['modulii'][1])
    age = 12.0  # from DES search, compared to 12.5, 0.0001, doesn't make much difference along main sequence
    z = 0.0002

    metal_poor = False
    ###################

    if surveys[survey]['fracdet'] is not None:
        print "Reading coverage fraction..."
        frac = fitsio.read(surveys[survey]['fracdet'])
        scale = (4096 / nside)**2
        pix = hp.nest2ring(nside, frac['PIXEL'] // scale)

        fracdet = np.zeros(hp.nside2npix(nside))
        np.add.at(fracdet, pix, frac['SIGNAL'])
        fracdet /= scale

    print "Reading catatogs..."
    dirname = surveys[survey]['data_dir']
    filenames = sorted(glob.glob(dirname + '/*.fits'))[:]

    if survey == 'PS1':
        pix = []
        for f in filenames:
            pix.append(int(f[-10:-5]))
        ang = hp.pix2ang(32, pix, nest=False, lonlat=True)
        c = SkyCoord(ang[0], ang[1], frame='icrs', unit='deg')
        b = c.galactic.b.deg
        BMIN = 20
        idx = np.abs(b) > BMIN
        filenames = np.asarray(filenames)[idx]

    data = load_infiles(filenames, columns=columns, multiproc=16)
    gc.collect()

    # Select magnitude range
    print "Selecting: %.1f < %s < %.1f" % (minmag, mag_g, maxmag)
    # data = data[(data[mag_g] < maxmag) & (data[mag_g] > minmag)]
    a1 = data[mag_g] < maxmag
    a2 = data[mag_g] > minmag
    a1 &= a2
    data = data[a1]
    gc.collect()

    mincolor = 0.
    maxcolor = 1.
    print "Selecting: %.1f < %s < %.1f" % (mincolor, 'g - r', maxcolor)
    a1 = data[mag_g] - data[mag_r] < maxcolor
    a2 = data[mag_g] - data[mag_r] > mincolor
    a1 &= a2
    data = data[a1]
    gc.collect()

    if ext is not None:
        data[mag_g] -= ext_g
        data[mag_r] -= ext_r
        data[mag_i] -= ext_i

    if stargal is not None:
        a1 = data[surveys[survey]['stargal']] <= surveys[survey]['stargal_cut']
        data = data[a1]
        gc.collect()

    hpxcube = np.zeros((hp.nside2npix(nside), len(modulii)))
    for i, mod in enumerate(modulii):
        print(" bin=%i: m-M = %.1f..." % (i, mod))

        C = surveys[survey]['C']
        E = surveys[survey]['E']
        err = surveys[survey]['err']
        gmin = 19.5 - (16.8 - mod)
        sel = filter_data.select_isochrone(data[mag_g], data[mag_r], err=err, iso_params=[mod, age, z], C=C, E=E, gmin=gmin)
        d = data[sel]

        if metal_poor:
            sel = filter_data.select_metal_poor(data[mag_g], data[mag_r], data[mag_i])
            d = data[sel]

        pixel = healpix.ang2pix(nside, d['RA'], d['DEC'])
        pix, cts = np.unique(pixel, return_counts=True)
        hpxcube[pix, i] = cts

    print("Writing %s..." % args.filename)
    header = healpix.header_odict(nside, coord='C', partial=False).values()
    f = fitsio.FITS(args.filename, 'rw', clobber=True)
    print("  Writing hpxcube...")
    f.write(hpxcube, header=header, extname='hpxcube')
    if surveys[survey]['fracdet'] is not None:
        print("  Writing fracdet...")
        f.write(fracdet, extname='fracdet', header=header)
    print("  Writing bins...")
    f.write(modulii, extname='modulus', header={'dmu': dmu})
    f.close()
