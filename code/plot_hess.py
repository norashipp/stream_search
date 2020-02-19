#!/usr/bin/env python -W ignore::DeprecationWarning

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.ndimage.filters import gaussian_filter
import healpy as hp


def plot_pretty(dpi=175, fontsize=15):
    # import pyplot and set some parameters to make plots prettier
    plt.rc("savefig", dpi=dpi)
    plt.rc('text', usetex=True)
    plt.rc('font', size=fontsize)
    plt.rc('xtick.major', pad=5)
    plt.rc('xtick.minor', pad=5)
    plt.rc('ytick.major', pad=5)
    plt.rc('ytick.minor', pad=5)


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def area_correction(data_on, data_off, nside=256):
    npix_on = float(
        len(np.unique(hp.ang2pix(nside, data_on['RA'], data_on['DEC'], lonlat=True))))
    npix_off = float(
        len(np.unique(hp.ang2pix(nside, data_off['RA'], data_off['DEC'], lonlat=True))))
    return npix_on / npix_off


def plot_hess(stream, data_on, data_off=None, gmax=23.5, gmin=16, grmax=1, grmin=0, vmin=None, vmax=None, dx=0.2 / 5., dy=1. / 6., ax=None, gband='PSF_MAG_SFD_G', rband='PSF_MAG_SFD_R', no_bkg=False, smoothing=0.75):
    print('dx, dy = ', dx, dy)
    xbins = np.arange(grmin, grmax + dx, dx)
    ybins = np.arange(gmin, gmax + dy, dy)

    g1 = data_on[gband]
    r1 = data_on[rband]
    gr1 = g1 - r1

    h_on, bx, by = np.histogram2d(gr1, g1, bins=[xbins, ybins])
    if data_off is None:
        h_off = 0
    else:
        g2 = data_off[gband]
        r2 = data_off[rband]
        gr2 = g2 - r2

        h_off, bx, by = np.histogram2d(gr2, g2, bins=[xbins, ybins])
        h_off *= area_correction(data_on, data_off)
    hdiff = h_on - h_off
    if no_bkg:
        hdiff = h_on

    smooth = gaussian_filter(hdiff, smoothing)

    if ax == None:
        plt.figure(figsize=(5, 6))
    else:
        plt.sca(ax)

    plt.title(r'$\mathrm{%s}$' % stream)
    im1 = plt.imshow(smooth.T, origin='upper', aspect='auto', extent=[
                     grmin, grmax, gmax, gmin], interpolation='none', cmap='binary', rasterized=True, vmin=vmin, vmax=vmax)
    colorbar(im1)

    plt.xlabel(r'$g-r$')
    plt.ylabel(r'$g$')

    return im1
