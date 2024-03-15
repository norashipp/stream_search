import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import healpy as hp
import os
import sys
_path = os.path.abspath('/Users/nora/code/slegs')
sys.path.append(_path)

import astropy
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import fits as fitsio
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

import gala.coordinates as gc
# import galstreams

from numpy.polynomial import polynomial
from polyfit2d import polyfit2d
sys.path.append('/Users/norashipp/code/slegs')
from stream_helpers import get_rot
import streamlib

from utils import mod2dist, dist2mod


def get_vscale(img, q=[5, 95]):
    from astropy.stats import sigma_clip

    _tmp = img.ravel()
    stdfunc = lambda x, axis: 1.5 * np.median(np.abs(x - np.median(x, axis=axis)), axis=axis)
    _tmp_clipped = sigma_clip(_tmp[_tmp != 0], stdfunc=stdfunc)

    return np.percentile(_tmp_clipped[_tmp_clipped > 0], q)


def plot_stream(hpxcube, modulus, stream=None, ends=None, mu_min=14.0, mu_max=20.0, reso=3.0, smoothing=0.15, savedir='../plots/', save=True, plot_streams=False):
    import galstreams
    print(mu_min, mu_max)

    if ends is None:
        if isinstance(stream, str):
            mw_streams = galstreams.MWStreams(verbose=False)
            stream = mw_streams[stream]
        # ends = [[stream.end_o.ra.deg, stream.end_o.dec.deg], [stream.end_f.ra.deg, stream.end_f.dec.deg]]
        ends = [[stream.end_f.ra.deg, stream.end_f.dec.deg], [stream.end_o.ra.deg, stream.end_o.dec.deg]]
        name = stream.name.replace(' ', '_')
    else:
        name = 'zoom'
    print(ends)
    fr = gc.GreatCircleICRSFrame.from_endpoints(coord.SkyCoord(ends[0][0], ends[0][1], unit='deg'),
                                                coord.SkyCoord(ends[1][0], ends[1][1], unit='deg'))
    phi, theta, psi = get_rot(fr)
    # phi, theta, psi = zyx_euler_from_endpoints(ends[0][0], ends[0][1], ends[1][0], ends[1][1])
    print(phi, theta, psi)

    if savedir is None:
        savedir = '../plots/%s/' % name
    if not os.path.isdir(savedir):
        os.mkdir(savedir)

    dmu = 0.1
    for mu in np.arange(mu_min, mu_max + dmu / 2., dmu):
        # for mu in [18.5]:
        print('m-M = %.1f' % mu)
        hpxmap = hpxcube[:, np.argmin(np.abs(mu - modulus))]
        nside = hp.npix2nside(hpxmap.shape[0])
        func = lambda x, y, z: hp.vec2pix(nside, x, y, z)

        if smoothing > 0:
            hpxmap_smooth = hp.smoothing(hpxmap, sigma=np.radians(smoothing), verbose=False)
        else:
            hpxmap_smooth = hpxmap
        proj = hp.projector.GnomonicProj(xsize=1024, ysize=800, rot=[phi, theta, psi], reso=reso)
        img = proj.projmap(hpxmap_smooth, func)

        vmin, vmax = get_vscale(img, q=[1, 90])
        print('vmin, vmax = ', vmin, vmax)

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.imshow(img, origin='bottom', vmin=vmin, vmax=vmax, cmap='Greys', extent=proj.get_extent())
        plt.title('m-M = %.1f' % mu, fontsize=22)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        if plot_streams:
            plot_stream_footprints(ax, proj, mu, dmu=100)

        if save:
            if plot_streams:
                savename = savedir + 'labeled_%s_%.1f.png' % (name, mu)
            else:
                savename = savedir + '%s_%.1f.png' % (name, mu)
            print('Saving %s...' % savename)
            plt.savefig(savename)
            plt.close('all')

    if save:
        os.system('convert -delay 15 -quality 100 %s/%s_*.png %s/%s.gif' % (savedir, name, savedir, name))
    else:
        return ax, proj


def plot_dwarfs_gcs(ax, proj, mu, dmu=0.5):
    xlim = ax.get_xlim()
    ylim = ax.set_ylim()

    dwarfs_globs = fitsio.open('../data/MW_dwarfs_globs.fits')[1].data

    x_array, y_array, name_array = [], [], []
    for dg in dwarfs_globs:
        x, y = proj.ang2xy(dg['ra'], dg['dec'], lonlat=True)
        if np.abs(mu - dg['DM']) < dmu and x > xlim[0] and x < xlim[1] and y > ylim[0] and y < ylim[1] and not np.isnan(x) and not np.isnan(y):
            x_array.append(x)
            y_array.append(y)
            name_array.append(dg['name'])

    x_array = np.asarray(x_array)
    y_array = np.asarray(y_array)
    name_array = np.asarray(name_array, dtype='str')

    ax.scatter(x_array, y_array, s=50, marker='o', facecolors='none', edgecolors='b')
    for i in range(len(x_array)):
        # name = (r'$\mathrm{%s}$' %name_array[i]).replace(' ', '\ ')
        name = name_array[i]
        plt.annotate(name, (x_array[i] + 0.01, y_array[i] + 0.01), color='b')

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])


def plot_stream_footprints(ax, proj, mu, dmu=0.5):
    import galstreams
    xlim = ax.get_xlim()
    ylim = ax.set_ylim()

    mw_streams = galstreams.MWStreams(verbose=False)
    for stream_name in mw_streams.keys():
    # for stream_name in ['ATLAS', 'Aliqa_Uma', 'Elqui', 'Chenab', 'Indus', 'Jhelum', 'Phoenix', 'Ravi', 'Turbio', 'Turranburra', 'Wambelong', 'TucanaIII', 'Willka_Yaku']:
        if stream_name in ['Her-Aq', 'EriPhe', 'Sgr-L10']:
            continue
        mu_stream = dist2mod(mw_streams[stream_name].Rhel[0])
        if np.abs(mu - mu_stream) < dmu:
            x, y = proj.ang2xy(mw_streams[stream_name].ra, mw_streams[stream_name].dec, lonlat=True)
            idx = ~np.isnan(x) & ~np.isnan(y)
            x, y = x[idx], y[idx]
            if len(x) < 1 or x.max() < xlim[0] or x.min() > xlim[1] or y.max() < ylim[0] or y.min() > ylim[1]:
                continue
            ax.plot(x, y, '.', alpha=0.1, label=mw_streams[stream_name].name.replace('_', ' '))

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    leg = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize=15, frameon=False, ncol=2)
    for lh in leg.legendHandles:
        lh._legmarker.set_alpha(1)


def get_gnomonic_proj(stream=None, ends=None, rot=None, xsize=1024, ysize=800, reso=3.0):
    if rot is not None:
        phi1, theta, psi = rot
    elif ends is not None:
        print('Using ends = ', ends)
        phi, theta, psi = get_euler(ends=ends)
    else:
        try:
            print('Using stream = ', stream.name)
        except:
            print('Using stream = ', stream)
        phi, theta, psi = get_euler(stream=stream)

    proj = hp.projector.GnomonicProj(xsize=xsize, ysize=ysize, rot=[phi, theta, psi], reso=reso)

    return proj


def get_mollweide_proj(stream=None, ends=None, xsize=1024):
    if ends is not None:
        phi, theta, psi = get_euler(ends=ends)
    else:
        phi, theta, psi = get_euler(stream=stream)

    proj = hp.projector.MollweideProj(xsize=xsize, rot=[phi, theta, psi])

    return proj


def get_ortho_proj(stream=None, ends=None, half_sky=True, **kwargs):
    if ends is not None:
        phi, theta, psi = get_euler(ends=ends)
    else:
        phi, theta, psi = get_euler(stream=stream)

    proj = hp.projector.OrthographicProj(rot=[phi, theta, psi], half_sky=half_sky, **kwargs)

    return proj


def plot_proj(proj, data, ax=None, q=[5, 90], vmin=None, vmax=None):
    nside = hp.npix2nside(len(data))
    func = lambda x, y, z: hp.vec2pix(nside, x, y, z)

    img = proj.projmap(data, func)

    if vmin is None or vmax is None:
        vmin, vmax = get_vscale(img, q=q)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    im = ax.imshow(img, origin='lower', vmin=vmin, vmax=vmax, cmap='Greys', extent=proj.get_extent())

    return ax, im
