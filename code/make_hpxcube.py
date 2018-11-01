import os
import glob
import logging
from collections import OrderedDict as odict

import fitsio
import numpy as np
import healpy as hp
import pylab as plt
import scipy.ndimage as nd
from matplotlib.path import Path

from utils import load_infiles
from ugali.utils import healpix
from ugali.analysis.isochrone import factory as isochrone_factory

import surveys
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
    mag = surveys.surveys[survey]['mag']
    mag_g = mag % G
    mag_r = mag % R
    mag_i = mag % I
    ext = surveys.surveys[survey]['ext']
    if ext is not None:
        ext = surveys.surveys[survey]['ext']
        ext_g = ext % G
        ext_r = ext % R
        ext_i = ext % I
    minmag = surveys.surveys[survey]['minmag']
    maxmag = surveys.surveys[survey]['maxmag']
    columns = ['RA', 'DEC', mag_g, mag_r, mag_i]

    ###################
    dmu = 0.1
    modulii = np.arange(15, 20 + dmu, dmu)
    age = 12.0
    z = 0.0002

    metal_poor = True
    ###################

    if surveys.surveys[survey]['fracdet'] is not None:
        print "Reading coverage fraction..."
        frac = fitsio.read(surveys.surveys[survey]['fracdet'])
        scale = (4096 / nside)**2
        pix = hp.nest2ring(nside, frac['PIXEL'] // scale)

        fracdet = np.zeros(hp.nside2npix(nside))
        np.add.at(fracdet, pix, frac['SIGNAL'])
        fracdet /= scale

    print "Reading catatogs..."
    dirname = surveys.surveys[survey]['data_dir']
    filenames = sorted(glob.glob(dirname + '/*.fits'))[:]
    data = load_infiles(filenames, columns=columns, multiproc=8)

    # Select magnitude range
    print "Selecting: %.1f < %s < %.1f" % (minmag, mag_g, maxmag)
    data = data[(data[mag_g] < maxmag) & (data[mag_g] > minmag)]

    if ext is not None:
        data[mag_g] -= ext_g
        data[mag_r] -= ext_r
        data[mag_i] -= ext_i

    hpxcube = np.zeros((hp.nside2npix(nside), len(modulii)))
    for i, mod in enumerate(modulii):
        print(" bin=%i: m-M = %.1f..." % (i, mod))

        sel = filter_data.select_isochrone(data[mag_g], data[mag_r], mod, age=age, z=z, dmu=dmu, C=[0.05, 0.1], E=2)
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
    if surveys.surveys[survey]['fracdet'] is not None:
        print("  Writing fracdet...")
        f.write(fracdet, extname='fracdet', header=header)
    print("  Writing bins...")
    f.write(modulus, extname='modulus', header={'dmu': dmu})
    f.close()
