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

from multiprocessing import Pool
from multiprocessing import Process, Value, Array
from multiprocessing import sharedctypes
# import cPickle as pickle

from astropy import units as u
from astropy.coordinates import SkyCoord

from utils_alex import load_infiles
from ugali.utils import healpix
# from ugali.analysis.isochrone import factory as isochrone_factory

from surveys import surveys
import filter_data


GRZ = False
GI = False


def run(arguments):
    mod = arguments
    print("m-M = %.1f..." % (mod))

    gmin = 3.4 + mod  # abs mag cutoff

    if GI:
        print('Using g-i sel')
        sel1 = filter_data.select_isochrone_gi(data[mag_g], data[mag_i], err=err, iso_params=[mod, age, z], C=C, E=E, gmin=gmin, survey=survey)
    else:
        print('Using g-r sel')
        sel1 = filter_data.select_isochrone(data[mag_g], data[mag_r], err=err, iso_params=[mod, age, z], C=C, E=E, gmin=gmin, survey=survey)

    if GRZ:
        print('Using g-r-z sel')
        sel2 = filter_data.select_isochrone_grz(data[mag_g], data[mag_r], data[mag_z], err=err, iso_params=[mod, age, z], C=C, E=E, gmin=gmin, survey=survey)
        sel = sel1 & sel2
    else:
        sel = sel1

    d = data[sel]

    pixel = healpix.ang2pix(nside, d['RA'], d['DEC'])
    pix, cts = np.unique(pixel, return_counts=True)
    return mod, pix, cts

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-f', '--filename', default='iso_hpxcube.fits.gz')
    parser.add_argument('-n', '--nside', default=512)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-s', '--survey', default='DES_Y3')
    parser.add_argument('-mp', '--multiproc', default=0)
    parser.add_argument('-a', '--age', default=12.)
    parser.add_argument('-z', '--metallicity', default=0.0004)
    parser.add_argument('-mm', '--maxmag', default=None)
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format='%(message)s', level=level)
    nside = args.nside

    survey = args.survey
    print('Filtering %s...' % survey)
    mag = surveys[survey]['mag']
    mag_g = mag % 'G'
    mag_r = mag % 'R'
    mag_i = mag % 'I'
    mag_z = mag % 'Z'
    ext = surveys[survey]['ext']
    stargal = surveys[survey]['stargal']
    stargal_cut = surveys[survey]['stargal_cut']
    # if ext is not None:
    #     # fix this, not relevant for now
    #     ext = surveys[survey]['ext']
    #     ext_g = ext % 'G'
    #     ext_r = ext % 'R'
    #     ext_i = ext % 'I'
    # columns = ['RA', 'DEC', mag_g, mag_r, mag_i] # include i eventually, try
    # mp search
    if GRZ:
        columns = ['RA', 'DEC', mag_g, mag_r, mag_z]
    elif GI:
        columns = ['RA', 'DEC', mag_g, mag_i]
    else:
        columns = ['RA', 'DEC', mag_g, mag_r]
    if stargal is not None:
        columns.append(stargal)
    if survey == 'DELVE':
        # FOR DELVE
        print('DELVE stargal')
        columns.append('WAVG_SPREAD_MODEL_G')

    ###################
    dmu = 0.1
    # moduli = np.arange(15, 20 + dmu, dmu)
    moduli = np.arange(surveys[survey]['moduli'][0], surveys[
                       survey]['moduli'][1] + dmu, dmu)

    print('Moduli: ', moduli)
    # 12.0  # from DES search, compared to 12.5, 0.0001, doesn't make much
    # difference along main sequence
    age = float(args.age)
    z = float(args.metallicity)  # 0.0002
    maxmag = float(args.maxmag)
    minmag = surveys[survey]['minmag']
    if maxmag is None:
        maxmag = surveys[survey]['maxmag']

    metal_poor = False
    ###################

    if surveys[survey]['fracdet'] is not None:
        print("Reading coverage fraction...")
        frac = fitsio.read(surveys[survey]['fracdet'])
        scale = (4096 / nside)**2
        pix = hp.nest2ring(nside, frac['PIXEL'] // scale)

        fracdet = np.zeros(hp.nside2npix(nside))
        np.add.at(fracdet, pix, frac['SIGNAL'])
        fracdet /= scale

    dirname = surveys[survey]['data_dir']
    print("Reading catalogs from %s..." % dirname)
    filenames = sorted(glob.glob(dirname + '/*.fits'))[:]

    if survey == 'PS1':
        pix = []
        for f in filenames:
            pix.append(int(f[-10:-5]))
        ang = hp.pix2ang(32, pix, nest=False, lonlat=True)
        c = SkyCoord(ang[0], ang[1], frame='icrs', unit='deg')
        b = c.galactic.b.deg

        BMIN = 0
        BMAX = 1000
        # idx = np.abs(b) > BMIN
        # idx = np.abs(b) <= BMAX
        idx = (np.abs(b) >= BMIN) & (np.abs(b) < BMAX)

        filenames = np.asarray(filenames)[idx]

    data = load_infiles(filenames, columns=columns, multiproc=16)
    # outfile = '/data/des40.b/data/nshipp/stream_search/%s_pickle.dat' % survey
    # pickle.dump(data, outfile)
    gc.collect()

    # Select magnitude range
    print("Selecting: %.1f < %s < %.1f" % (minmag, mag_g, maxmag))
    # data = data[(data[mag_g] < maxmag) & (data[mag_g] > minmag)]
    a1 = data[mag_g] < maxmag
    a2 = data[mag_g] > minmag
    a1 &= a2
    data = data[a1]
    gc.collect()

    if GI:
        mincolor = -1
        maxcolor = 2
        print("Selecting: %.1f < %s < %.1f" % (mincolor, 'g - r', maxcolor))
        a1 = data[mag_g] - data[mag_i] < maxcolor
        a2 = data[mag_g] - data[mag_i] > mincolor
    else:
        mincolor = 0.
        maxcolor = 1.
        print("Selecting: %.1f < %s < %.1f" % (mincolor, 'g - r', maxcolor))
        a1 = data[mag_g] - data[mag_r] < maxcolor
        a2 = data[mag_g] - data[mag_r] > mincolor
    a1 &= a2
    data = data[a1]
    gc.collect()

    if ext is not None:
        data[mag_g] -= ext_g
        data[mag_r] -= ext_r
        data[mag_i] -= ext_i

    if survey == 'DELVE':
        # FOR DELVE
        print('Using DELVE stargal: abs(wavg_spread_model_g) < 0.003')
        a1 = np.abs(data['WAVG_SPREAD_MODEL_G']) < 0.003
        data = data[a1]
        gc.collect()
    elif stargal is not None:
        print('Selecting: %s <= %i' % (stargal, stargal_cut))
        a1 = data[stargal] <= stargal_cut
        data = data[a1]
        gc.collect()

    if metal_poor:
        sel = filter_data.select_metal_poor(
            data[mag_g], data[mag_r], data[mag_i])
        data = data[sel]
        gc.collect()

    C = surveys[survey]['C']
    E = surveys[survey]['E']
    err = surveys[survey]['err']

    hpxcube = np.zeros((hp.nside2npix(nside), len(moduli)))

    # data = Array('d', data)
    # hpxcube = Array('d', hpxcube)
    # data = np.ctypeslib.as_ctypes(data)
    # hpxcube = np.ctypeslib.as_ctypes(hpxcube)
    # data = sharedctypes.RawArray(data._type_, data)
    # hpxcube = sharedctypes.RawArray(hpxcube._type_, hpxcube)

    arguments = moduli

    multiproc = int(args.multiproc)
    if multiproc:
        p = Pool(multiproc, maxtasksperchild=1)
        results = p.map(run, arguments)
    else:
        results = [r for r in map(run, arguments)]

    for res in results:
        mod, pix, cts = res
        i = np.argmin(np.abs(moduli - mod))
        hpxcube[pix, i] = cts

    # for i, mod in enumerate(moduli):
    #     print(" bin=%i: m-M = %.1f..." % (i, mod))

        # C = surveys[survey]['C']
        # E = surveys[survey]['E']
        # err = surveys[survey]['err']
        # gmin = 19.5 - (16.8 - mod)

        # if args.multiproc > 0:
        #     from multiprocessing import Pool
        #     p = Pool(args.multiproc, maxtasksperchild=1)
        #     sel = p.map(run, args)
        # else:
        #     sel = [run(arg) for arg in args]

        # sel = filter_data.select_isochrone(data[mag_g], data[mag_r], err=err, iso_params=[
        # mod, age, z], C=C, E=E, gmin=gmin, survey=survey)

        # d = data[sel]

        # if metal_poor:
        #     sel = filter_data.select_metal_poor(
        #         data[mag_g], data[mag_r], data[mag_i])
        #     d = data[sel]

        # pixel = healpix.ang2pix(nside, d['RA'], d['DEC'])
        # pix, cts = np.unique(pixel, return_counts=True)
        # hpxcube[pix, i] = cts

    print("Writing %s..." % args.filename)
    header = healpix.header_odict(nside, coord='C', partial=False).values()
    f = fitsio.FITS(args.filename, 'rw', clobber=True)
    print("  Writing hpxcube...")
    # f.write(hpxcube, header=header, extname='hpxcube')
    f.write(hpxcube, extname='hpxcube') # header not working
    if surveys[survey]['fracdet'] is not None:
        print("  Writing fracdet...")
        # f.write(fracdet, extname='fracdet', header=header)
        f.write(fracdet, extname='fracdet')
    print("  Writing bins...")
    f.write(moduli, extname='modulus', header={'dmu': dmu})
    f.close()
