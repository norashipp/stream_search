from __future__ import division

import os
import glob
import subprocess

import numpy as np
from collections import OrderedDict as odict
from scipy.ndimage.filters import gaussian_filter
import healpy as hp
import fitsio

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import skymap
from skymap.utils import cel2gal, gal2cel
from skymap.utils import setdefaults

import ugali.utils.healpix as uhp
# from ugali.utils import fileio

import galstreams

import streamlib
reload(streamlib)

import surveys


def plot_pretty(dpi=175, fontsize=15, labelsize=15, figsize=(10, 8)):
    # import pyplot and set some parameters to make plots prettier

    plt.rc('savefig', dpi=dpi)
    plt.rc('text', usetex=True)
    plt.rc('font', size=fontsize)
    plt.rc('xtick.major', pad=5)
    plt.rc('xtick.minor', pad=5)
    plt.rc('ytick.major', pad=5)
    plt.rc('ytick.minor', pad=5)
    plt.rc('figure', figsize=figsize)

    mpl.rcParams['xtick.labelsize'] = labelsize
    mpl.rcParams['ytick.labelsize'] = labelsize
    mpl.rcParams.update({'figure.autolayout': True})


def load_hpxcube(filename='../data/iso_hpxcube.fits.gz'):
    print("Reading %s..." % filename)
    f = fitsio.FITS(filename)
    hpxcube = f['HPXCUBE'].read()
    try:
        fracdet = f['FRACDET'].read()
    except:
        fracdet = np.zeros_like(hpxcube[:, 0])
        fracdet[np.where(hpxcube[:, 0] > 0)] = 1
    modulus = f['MODULUS'].read()

    return hpxcube, fracdet, modulus


def prepare_hpxmap(mu, hpxcube, fracdet, modulus, fracmin=0.5, clip=100, sigma=0.2, **mask_kw):
    i = np.argmin(np.abs(mu - modulus))
    hpxmap = np.copy(hpxcube[:, i])

    data = streamlib.prepare_data(hpxmap, fracdet, fracmin=fracmin, clip=clip, mask_kw=mask_kw)
    # bkg = streamlib.fit_bkg_poly(data, sigma=sigma)
    bkg = None
    return data, bkg


def plot_streams(smap, mu, dmu=0.5, coords='cel', filename=None):
    mw_streams = galstreams.MWStreams(verbose=False)
    for stream in mw_streams.keys():
        if stream in ['Her-Aq', 'EriPhe', 'Sgr-L10']:
            continue
        mu_stream = dist2mod(mw_streams[stream].Rhel[0])
        if np.abs(mu - mu_stream) < dmu:
            # print stream, mw_streams[stream].ra.max() - mw_streams[stream].ra.min(), mw_streams[stream].dec.max(), mw_streams[stream].dec.min()
            if coords == 'gal':
                x, y = smap(mw_streams[stream].l, mw_streams[stream].b)
            else:
                x, y = smap(mw_streams[stream].ra, mw_streams[stream].dec)
            smap.plot(x, y, '.', alpha=0.5)
            plt.gca().annotate(mw_streams[stream].name.replace('_', ' '), (x.min(), y.min()), color='c', fontsize=15)

    if filename:
        plt.savefig(filename)


def plot_stream_list(smap, streams, coords='cel', filename=None):
    mw_streams = galstreams.MWStreams(verbose=False)
    for stream in streams:
        if coords == 'gal':
            x, y = smap(mw_streams[stream].l, mw_streams[stream].b)
        else:
            x, y = smap(mw_streams[stream].ra, mw_streams[stream].dec)
        smap.plot(x, y, '.', alpha=0.5)
        plt.gca().annotate(mw_streams[stream].name.replace('_', ' '), (x.min(), y.min()), color='b', fontsize=15)

    if filename:
        plt.savefig(filename)


def plot_density(data, bkg, coords='cel', center=(0, 0), proj='mbtfpq', filename=None, **kwargs):
    defaults = dict(cmap='gray_r', xsize=400, smooth=0.2)
    setdefaults(kwargs, defaults)

    nside = hp.get_nside(data)

    lon, lat = hp.pix2ang(nside, np.arange(len(data)), lonlat=True)
    galpix = hp.ang2pix(nside, *gal2cel(lon, lat), lonlat=True)

    plt.figure()
    smap = skymap.Skymap(projection=proj, lon_0=center[0], lat_0=center[1], celestial=False)

    if coords == 'gal':
        # smap.draw_hpxmap((data - bkg)[galpix], **kwargs)
        smap.draw_hpxmap(data[galpix], **kwargs)
    else:
        # smap.draw_hpxmap((data - bkg), **kwargs)
        smap.draw_hpxmap(data, **kwargs)

    if filename:
        plt.savefig(filename)

    return smap


def mod2dist(distance_modulus):
    return 10**(distance_modulus / 5. + 1.) / 1e3


def dist2mod(distance):
    return 5. * (np.log10(np.array(distance) * 1.e3) - 1.)


def make_movie(infiles, outfile=None, delay=40, queue='local'):
    print("Making movie...")
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


if __name__ == "__main__":
    plot_pretty(figsize=(18, 14))
    filename = '../data/iso_hpxcube_ps1.fits.gz'
    movdir = '/data/des40.b/data/nshipp/stream_search/plots/ps1_cap/'
    movdir_labeled = '/data/des40.b/data/nshipp/stream_search/plots/ps1_cap/labeled/'
    hpxcube, fracdet, modulus = load_hpxcube(filename)
    for mu in modulus[:-16]:
        # for lon in [0, 180, -180]:
        #     print 'Plotting m-M = %.1f...' % mu
        #     data, bkg = prepare_hpxmap(mu, hpxcube, fracdet, modulus, plane=True, center=True, sgr=False, bmax=25, cmax=40)
        #     smap = plot_density(data, bkg, vmax=15, center=(lon, 30), filename=movdir + 'density_ps1_%.2f_%i.png' % (mu, lon))
        #     plot_streams(smap, mu, filename=movdir_labeled + 'density_ps1_%.2f_%i_labeled.png' % (mu, lon))

        if os.path.exists(movdir_labeled + 'density_ps1_cap_%.2f_labeled.png' % (mu)):
            print 'Skipping m-M = %.1f' % mu
            continue
        print 'Plotting m-M = %.1f...' % mu
        data, bkg = prepare_hpxmap(mu, hpxcube, fracdet, modulus, plane=True, center=True, sgr=False, bmax=25, cmax=40)
        smap = plot_density(data, bkg, vmax=10, center=(0, -90), proj='ortho', coords='gal', filename=movdir + 'density_ps1_cap_%.2f.png' % (mu))
        plot_streams(smap, mu, coords='gal', filename=movdir_labeled + 'density_ps1_cap_%.2f_labeled.png' % (mu))
