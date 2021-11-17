from __future__ import division

import os
import glob
import subprocess
# from importlib import reload

import numpy as np
from collections import OrderedDict as odict
from scipy.ndimage.filters import gaussian_filter
import healpy as hp
import fitsio

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import astropy.coordinates as coords  # NEED TO RENAME
import astropy.units as u

import skymap
from skymap.utils import cel2gal, gal2cel
from skymap.utils import setdefaults

import ugali.utils.healpix as uhp
from ugali.utils.projector import angsep
# from ugali.utils import fileio

import galstreams

from importlib import reload
import streamlib
reload(streamlib)
import rotation_matrix
import surveys
# reload(surveys)
import results
import elysian


# HOMEDIR = '/Users/nora/'
HOMEDIR = '/home/s1/nshipp/'


def plot_pretty(dpi=175, fontsize=15, labelsize=15, figsize=(10, 8), tex=True):
    # import pyplot and set some parameters to make plots prettier

    plt.rc('savefig', dpi=dpi)
    plt.rc('text', usetex=tex)
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
        print('fracdet test', np.sum(fracdet > 0.5))
    except:
        print('Skipping fracdet...')
        fracdet = np.zeros_like(hpxcube[:, 0])
        fracdet[np.where(hpxcube[:, 0] > 0)] = 1
    try:
        modulus = f['MODULUS'].read()
    except:
        print('Error reading modulus...')
        modulus = np.array([16.])
    return hpxcube, fracdet, modulus


def prepare_hpxmap(mu, hpxcube, fracdet, modulus, fracmin=0.5, clip=100, sigma=0.2, **mask_kw):
    # fracdet = np.ones(hpxmap.size)
    # fracdet[hpxmap == hp.UNSEEN] = 0

    i = np.argmin(np.abs(mu - modulus))
    hpxmap = np.copy(hpxcube[:, i])

    data = streamlib.prepare_data(
        hpxmap, fracdet, fracmin=fracmin, clip=clip, mask_kw=mask_kw)
    # bkg = streamlib.fit_bkg_poly(data, sigma=sigma)
    # bkg = None
    return data


def fit_background(data, center=(0, 0), coords='cel', coord_stream=None, proj='cyl', sigma=0.2, percent=[2, 95], deg=5):
    bkg = streamlib.fit_bkg_poly(data, center=center, coords=coords,
                                 coord_stream=coord_stream, proj=proj, sigma=sigma, percent=percent, deg=deg)
    return bkg


def plot_dwarfs_globs(smap, data, mu, dmu=0.5, coords='cel', coord_stream=None, filename=None):
    dwarfs_globs = fitsio.read(
        '../data/MW_dwarfs_globs.fits')

    nside = hp.get_nside(data)

    ra_array, dec_array, name_array = [], [], []
    for dg in dwarfs_globs:
        pix = hp.ang2pix(nside, dg['ra'], dg['dec'], lonlat=True)
        if np.abs(mu - dg['DM']) < dmu and ~data.mask[pix]:
            ra_array.append(dg['ra'])
            dec_array.append(dg['dec'])
            name_array.append(dg['name'])
    ra_array = np.asarray(ra_array)
    dec_array = np.asarray(dec_array)
    name_array = np.asarray(name_array, dtype='str')

    if coords == 'gal':
        l, b = cel2gal(ra_array, dec_array)
        x, y = smap(l, b)
    elif coords == 'stream':
        if not coord_stream:
            print('Need to input coord_stream!')
        R = streamlib.get_rotmat(stream=coord_stream)
        phi1, phi2 = rotation_matrix.phi12_rotmat(ra_array, dec_array, R)
        x, y = smap(phi1, phi2)
    elif coords == 'cel':
        x, y = smap(ra_array, dec_array)
    else:
        x, y = smap(ra_array, dec_array)
        print('Unknown coord, using cel.')

    smap.scatter(x, y, s=50, marker='o', facecolors='none', edgecolors='b')
    for i in range(len(x)):
        plt.annotate(name_array[i], (x[i] + 0.01, y[i] + 0.01), color='b')

    if filename:
        plt.savefig(filename)


def plot_streams(smap, mu, dmu=0.5, coords='cel', coord_stream=None, filename=None):
    mw_streams = galstreams.MWStreams(verbose=False)
    for stream in mw_streams.keys():
        if stream in ['Her-Aq', 'EriPhe', 'Sgr-L10']:
            continue
        mu_stream = dist2mod(mw_streams[stream].Rhel[0])
        if np.abs(mu - mu_stream) < dmu:
            # print stream, mw_streams[stream].ra.max() -
            # mw_streams[stream].ra.min(), mw_streams[stream].dec.max(),
            # mw_streams[stream].dec.min()
            if coords == 'gal':
                x, y = smap(mw_streams[stream].l, mw_streams[stream].b)
            elif coords == 'stream':
                if not coord_stream:
                    print('Need to input coord_stream!')
                ra, dec = mw_streams[stream].ra, mw_streams[stream].dec
                R = streamlib.get_rotmat(stream=coord_stream)
                phi1, phi2 = rotation_matrix.phi12_rotmat(ra, dec, R)
                x, y = smap(phi1, phi2)
            else:
                x, y = smap(mw_streams[stream].ra, mw_streams[stream].dec)
            smap.plot(x, y, '.', alpha=0.5)
            plt.gca().annotate(mw_streams[stream].name.replace(
                '_', ' '), (x[0], y[0]), color='navy', fontsize=20)

    if filename:
        plt.savefig(filename)


def plot_stream_list(smap, streams, coords='cel', coord_stream=None, filename=None):
    mw_streams = galstreams.MWStreams(verbose=False)
    for stream in streams:
        if coords == 'gal':
            x, y = smap(mw_streams[stream].l, mw_streams[stream].b)
        elif coords == 'stream':
            if not coord_stream:
                print('Need to input coord_stream!')
            ra, dec = mw_streams[stream].ra, mw_streams[stream].dec
            R = streamlib.get_rotmat(stream=coord_stream)
            phi1, phi2 = rotation_matrix.phi12_rotmat(ra, dec, R)
            x, y = smap(phi1, phi2)
        else:
            x, y = smap(mw_streams[stream].ra, mw_streams[stream].dec)
        smap.plot(x, y, '.', alpha=0.5)
        plt.gca().annotate(mw_streams[stream].name.replace(
            '_', ' '), (x.min(), y.min()), color='b', fontsize=15)

    if filename:
        plt.savefig(filename)


def plot_density(data, bkg, coords='cel', coord_stream=None, center=(0, 0), proj='mbtfpq', filename=None, smap=None, **kwargs):
    defaults = dict(cmap='gray_r', xsize=400, smooth=0.2)
    setdefaults(kwargs, defaults)

    nside = hp.get_nside(data)

    if smap is None:
        plt.figure()
        smap = skymap.Skymap(projection=proj, lon_0=center[0], lat_0=center[1], celestial=False)

    if coords == 'gal':
        lon, lat = hp.pix2ang(nside, np.arange(len(data)), lonlat=True)
        galpix = hp.ang2pix(nside, *gal2cel(lon, lat), lonlat=True)
        # smap.draw_hpxmap((data - bkg)[galpix], **kwargs)
        # IS THIS RIGHT?
        smap.draw_hpxmap((data[galpix] - bkg), **kwargs)
        # smap.draw_hpxmap(data[galpix], **kwargs)
    elif coords == 'stream':
        if not coord_stream:
            print('Need to input coord_stream!')
        streampix = streamlib.get_streampix(data=data, stream=coord_stream)
        smap.draw_hpxmap((data[streampix] - bkg), **kwargs)

    else:
        smap.draw_hpxmap((data - bkg), **kwargs)
        # smap.draw_hpxmap(data, **kwargs)

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


def plot_stream_zoom(hpxcube, fracdet, modulus, stream=None, ends=None, mu=None, width=0.3, sigma=0.2, delta=0.2, vmin=-1, vmax=5, filename=None):
    if ends is not None:
        pass        
    elif stream is not None:
        mw_streams = galstreams.MWStreams(verbose=False)
        ends = [(mw_streams[stream].end_f.ra.value, mw_streams[stream].end_f.dec.value),
                (mw_streams[stream].end_o.ra.value, mw_streams[stream].end_o.dec.value)]
        # center = (mw_streams[stream].ra.mean(), mw_streams[stream].dec.mean())
        mu = dist2mod(mw_streams[stream].Rhel[0])
    elif ends is not None and mu is not None:
        pass
    else:
        print('Need stream name or ends and modulus!')

    length = angsep(ends[0][0], ends[0][1], ends[1][0], ends[1][1])

    data = prepare_hpxmap(mu, hpxcube, fracdet, modulus, clip=100, sigma=sigma,
                          plane=False, center=False, sgr=False, bmax=25, cmax=40)
    bkg = fit_background(data, center=(
        0, 0), coords='cel', coord_stream=None, sigma=sigma, deg=5)

    if sigma:
        data.fill_value = np.ma.median(data)
        data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False),
                           mask=data.mask, fill_value=np.nan)

    nside = hp.get_nside(data)

    height = max(1.2 * 6 * width, 5)
    length = max(1.2 * length, 1.5 * height)

    if ends is not None:
        streampix = streamlib.get_streampix(data=data, stream=None, ends=ends)
    else:
        streampix = streamlib.get_streampix(data=data, stream=stream, ends=ends)

    # delta = 0.1
    aspect = float(height) / length
    b = np.arange(-height, height, delta)
    l = np.arange(-length, length, delta)
    ll, bb = np.meshgrid((l[1:] + l[:-1]) / 2., (b[1:] + b[:-1]) / 2.)

    phi, theta, psi = results.euler_angles(
        ends[0][0], ends[0][1], ends[1][0], ends[1][1])

    if stream is None:
        stream = 'New'
    stream_cls = elysian.euler_factory(stream, phi, theta, psi, ends=ends)

    rot = stream_cls(Lambda=ll.flatten() * u.deg, Beta=bb.flatten() * u.deg)
    icrs = rot.transform_to(coords.ICRS)
    ends = coords.ICRS(*(stream_cls.ends.T * u.deg)).transform_to(stream_cls)
    ends.Lambda.wrap_angle = 180 * u.deg

    pix = hp.ang2pix(nside, icrs.ra.deg, icrs.dec.deg, lonlat=True)
    value = (data - bkg)[pix].reshape(ll.shape)
    mask = (data.mask)[pix].reshape(ll.shape)

    fig, ax = plt.subplots(1, 1, figsize=(5.0 / aspect, 5.0))
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes('right', size=1.2, pad=0.1, sharey=ax)
    cax = divider.append_axes('top', size='7%', pad=0.05)

    plt.setp(ax2.get_yticklabels(), visible=False)
    # Do the plotting...
    kw = dict(vmin=vmin, vmax=vmax, rasterized=True, cmap='gray_r')
    # if stream == 'ATLAS':
    #     kw.update(vmin=-1.5, vmax=2.0)
    im = ax.pcolormesh(l, b, value, **kw)

    plt.colorbar(im, cax, orientation='horizontal')
    cax.tick_params(axis='x', top=True, bottom=False,
                    labeltop=True, labelbottom=False)

    # ax.plot(ends.L.deg, ends.B.deg, ls='--', lw=1.5, color='deepskyblue')
    cts = value.sum(axis=1)
    area = (~mask).sum(axis=1)
    ax2.plot(cts / area, (b[1:] + b[:-1]) / 2., color='k')

    label = r'${\rm %s}$' % stream.replace(' ', '\ ')
    ax2.annotate(label, (0.95, 0.95), xycoords='axes fraction',
                 fontsize=11, ha='right', va='top')
    # And the labeling
    ax.set_xlabel(r'$\Lambda\ {\rm (deg)}$')
    ax.set_ylabel(r'${\rm B\ (deg)}$')
    ax.set_xlim(l.min(), l.max())
    ax.set_ylim(b.min(), b.max())
    # plt.suptitle(name)
    plt.subplots_adjust(bottom=0.05, top=0.98, right=0.98)

    if filename:
        plt.savefig(filename)

    return ax


def plot_stream(stream, hpxcube, fracdet, modulus):
    mw_streams = galstreams.MWStreams(verbose=False)
    center = (mw_streams[stream].ra.mean(), mw_streams[stream].dec.mean())
    mu = dist2mod(mw_streams[stream].Rhel[0])

    data = prepare_hpxmap(mu, hpxcube, fracdet, modulus, clip=100, sigma=sigma,
                          plane=False, center=False, sgr=False, bmax=25, cmax=40)
    bkg = fit_background(data, center=(
        0, 0), coords='cel', coord_stream=None, sigma=sigma, deg=5)

    smap = plot_density(data, bkg, center=center, vmax=8,
                        coords='cel', proj='ortho', xsize=600)

    smap = plot_density(data, bkg, center=center, vmax=8,
                        coords='cel', proj='ortho', xsize=600)
    plot_stream_list(smap, [stream])

    smap = plot_density(data, bkg, center=center, vmax=8,
                        coords='cel', proj='ortho', xsize=600)
    plot_streams(smap, mu, 50)


if __name__ == "__main__":
    plot_pretty(figsize=(18, 14))

    survey = 'BASS' # 'DECaLS'
    proj = 'ortho'

    # if survey == 'DECaLS':
    #     coords = 'stream'
    #     coord_stream = 'Lethe'
    #     vmin, vmax = 0, 15
    #     xsize = 1000
    #     sigma = 0.3
    #     version = 4
    #     filename = '../data/decals_dr8_iso_hpxcube_v%i.fits.gz' % version
    #     movdir = '../plots/decals/v%i/' % version

    #     # stream = 'Lethe'
    #     mw_streams = galstreams.MWStreams(verbose=False)
    #     if coords == 'cel':
    #         center = (mw_streams[stream].ra.mean(),
    #                   mw_streams[stream].dec.mean())
    #     elif coords == 'gal':
    #         center = (mw_streams[stream].l.mean(), mw_streams[stream].b.mean())
    #     elif coords == 'stream':
    #         center = 0, 0
    #     background_center = center

    if survey == 'DECaLS':
        coord_stream = None
        coords = 'cel'
        lon, lat = -5, 10 # 190, 0
        center = (lon, lat)
        # vmin, vmax = 0, 12
        vmin, vmax = 0, 15
        xsize = 1000
        sigma = 0.2
        version = 6
        filename = '../data/decals_dr8_iso_hpxcube_v%i.fits.gz' % version
        movdir = '../plots/decals/v%i/' % version

        # stream = 'Lethe'
        mw_streams = galstreams.MWStreams(verbose=False)
        if center is None:
            if coords == 'cel':
                center = (mw_streams[stream].ra.mean(),
                          mw_streams[stream].dec.mean())
            elif coords == 'gal':
                center = (mw_streams[stream].l.mean(), mw_streams[stream].b.mean())
            elif coords == 'stream':
                center = 0, 0
        background_center = center

    if survey == 'BASS':
        coord_stream = None
        coords = 'cel'
        lon, lat = 190, 80
        # coords = 'gal'
        # lon, lat = 120, 50
        vmin, vmax = -5, 3
        xsize = 1000
        sigma = 0.2
        version = 6
        filename = '../data/bass_dr8_iso_hpxcube_v%i.fits.gz' % version
        movdir = '../plots/bass/v%i/' % version
        center = (lon, lat)
        background_center = (190, 80)

    if center[0] > 180:
        center = (center[0] - 360, center[1])

    hpxcube, fracdet, modulus = load_hpxcube(filename)
    modulus = np.arange(14, 18.05, 0.1)

    for mu in modulus:
        if os.path.exists(movdir + 'density_%s_%s_%.2f.png' % (survey, coord_stream, mu)):
            print('Skipping m-M = %.1f' % mu)
            continue
        print('Plotting m-M = %.1f...' % mu)
        data = prepare_hpxmap(mu, hpxcube, fracdet, modulus, clip=100, plane=False, center=False, sgr=False, bmax=25, cmax=40, sigma=sigma)
        # bkg = 0
        bkg = fit_background(data, center=background_center, coords=coords, coord_stream=coord_stream, sigma=sigma, deg=5)
        smap = plot_density(data, bkg, coord_stream=coord_stream, center=center, vmin=vmin, vmax=vmax, coords=coords, proj=proj,
                            xsize=xsize, smooth=sigma, filename=movdir + 'density_%s_%s_%.2f.png' % (survey, coord_stream, mu))

        # plot_streams(smap, mu, dmu=2, coords=coords, coord_stream=coord_stream,
        # filename=movdir + 'density_%s_%s_%.2f_streams.png' % (survey,
        # coord_stream, mu))

        # smap = plot_density(data, bkg, coord_stream=coord_stream, center=center,
        #                     vmin=vmin, vmax=vmax, coords=coords, proj='ortho', xsize=3000)
        # plot_dwarfs_globs(smap, data, mu, coords=coords, coord_stream=coord_stream,
        # filename=movdir + 'density_%s_%s_%.2f_globs.png' % (survey,
        # coord_stream, mu))

        plt.close('all')
