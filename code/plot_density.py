from __future__ import division

import os
import glob

import numpy as np
from collections import OrderedDict as odict
from scipy.ndimage.filters import gaussian_filter
import healpy as hp

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import skymap
from skymap.utils import cel2gal, gal2cel
import ugali.utils.healpix as uhp
from ugali.utils import fileio

import surveys


def plot_pretty(dpi=175, fontsize=15, labelsize=15):
    # import pyplot and set some parameters to make plots prettier

    plt.rc('savefig', dpi=dpi)
    plt.rc('text', usetex=True)
    plt.rc('font', size=fontsize)
    plt.rc('xtick.major', pad=5)
    plt.rc('xtick.minor', pad=5)
    plt.rc('ytick.major', pad=5)
    plt.rc('ytick.minor', pad=5)

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
        fractdet = np.zeros_like(hpxcube)
        fracdet[np.where(hpxcube > 0)] = 1
    modulus = f['MODULUS'].read()

    return hpxcube, fracdet, modulus


def plot_density(mu, hpxmap, fracdet, modulus, sigma=0.2, fracmin=0.5, filename=None):
    i = np.argmin(np.abs(mu - modulus))
    hpxmap = np.copy(hpxcube[:, i])

    mask_kw = dict(lmc=False, milky_way=False, sgr=False, globulars=False, dwarfs=False, galaxies=False)
    data = streamlib.prepare_data(hpxmap, fracdet, fracmin=fracmin, mask_kw=mask_kw)
    bkg = streamlib.fit_bkg_poly(data, sigma=sigma)

    kwargs = dict(cmap='gray_r', xsize=400, smooth=sigma)

    plt.figure()
    smap = skymap.Skymap(projection='mbtfpq', lon_0=0, lat_0=0)
    smap.draw_hpxmap((data - bkg), **kwargs)

    if filename:
        plt.savefig(filename)


def make_movie():
    pass


if __name__ == "__main__":
    filename = 'iso_hpxcube_ps1.fits.gz'
    hpxcube, fracdet, modulus = load_hpxcube(filename)
    mu = 16.8
    plot_density(mu, hpxcube, fracdet, modulus, filename='density_ps1_%.2f.png' % mu)
