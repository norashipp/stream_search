import numpy as np
import matplotlib.pyplot as plt
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
import galstreams

from numpy.polynomial import polynomial
from polyfit2d import polyfit2d
from stream_helpers import get_rot
import streamlib


def mod2dist(distance_modulus):
    return 10**(distance_modulus / 5. + 1.) / 1e3


def dist2mod(distance):
    return 5. * (np.log10(np.array(distance) * 1.e3) - 1.)


def zyx_euler_from_endpoints(lon1, lat1, lon2, lat2):
    c1 = coord.SkyCoord(lon1 * u.deg, lat1 * u.deg)
    c2 = coord.SkyCoord(lon2 * u.deg, lat2 * u.deg)
    fr = gc.GreatCircleICRSFrame.from_endpoints(c1, c2)
    origin = fr.realize_frame(coord.UnitSphericalRepresentation(0 * u.deg, 0 * u.deg))

    gc_icrs = origin.transform_to(coord.ICRS)
    R = gc.greatcircle.reference_to_greatcircle(coord.ICRS, fr)
    psi = -np.degrees(np.arctan2(R[2, 1], R[2, 2]))

    return [gc_icrs.ra.degree, gc_icrs.dec.degree, psi]


def get_vscale(img, q=[5, 95]):
    from astropy.stats import sigma_clip

    _tmp = img.ravel()
    stdfunc = lambda x, axis: 1.5 * np.median(np.abs(x - np.median(x, axis=axis)), axis=axis)
    _tmp_clipped = sigma_clip(_tmp[_tmp != 0], stdfunc=stdfunc)

    return np.percentile(_tmp_clipped[_tmp_clipped > 0], q)


def plot_stream(hpxcube, modulus, stream=None, ends=None, mu_min=14.0, mu_max=20.0, reso=3.0, smoothing=0.15, savedir='../plots/', save=True, plot_streams=False):
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

        hpxmap_smooth = hp.smoothing(hpxmap, sigma=np.radians(smoothing), verbose=False)
        proj = hp.projector.GnomonicProj(xsize=1024, ysize=800, rot=[phi, theta, psi], reso=reso)
        img = proj.projmap(hpxmap_smooth, func)

        vmin, vmax = get_vscale(img, q=[1, 90])

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.imshow(img, origin='bottom', vmin=vmin, vmax=vmax, cmap='Greys', extent=proj.get_extent())
        plt.title('m-M = %.1f' % mu, fontsize=22)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        if plot_streams:
            plot_stream_footprints(ax, proj, mu, dmu=100)

        if plot_streams:
            leg = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1, fontsize=15, frameon=False)
            for lh in leg.legendHandles:
                lh._legmarker.set_alpha(1)

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
        plt.annotate(name_array[i], (x_array[i] + 0.01, y_array[i] + 0.01), color='b')

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])


def plot_stream_footprints(ax, proj, mu, dmu=0.5):
    xlim = ax.get_xlim()
    ylim = ax.set_ylim()

    mw_streams = galstreams.MWStreams(verbose=False)
    for stream_name in mw_streams.keys():
        if stream_name in ['Her-Aq', 'EriPhe', 'Sgr-L10']:
            continue
        mu_stream = dist2mod(mw_streams[stream_name].Rhel[0])
        if np.abs(mu - mu_stream) < 100:
            x, y = proj.ang2xy(mw_streams[stream_name].ra, mw_streams[stream_name].dec, lonlat=True)
            idx = ~np.isnan(x) & ~np.isnan(y)
            x, y = x[idx], y[idx]
            if len(x) < 1 or x.max() < xlim[0] or x.min() > xlim[1] or y.max() < ylim[0] or y.min() > ylim[1]:
                continue
            ax.plot(x, y, '.', alpha=0.1, label=mw_streams[stream_name].name.replace('_', ' '))

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])


def fit_bkg(data, proj, sigma=0.1, percent=[2, 95], deg=5):
    nside = hp.get_nside(data.mask)
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)

    vmin, vmax = np.percentile(data.compressed(), q=percent)
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)

    smoothed = hp.smoothing(data, sigma=np.radians(sigma), verbose=False)

    # not sure how to use astropy convolve with healpix
    # kernel = Gaussian2DKernel(x_stddev=sigma)
    # smoothed = convolve(data, kernel)

    data = np.ma.array(smoothed, mask=data.mask)

    sel = ~data.mask
    x, y = proj.ang2xy(lon[sel], lat[sel], lonlat=True)
    
    xmin, xmax, ymin, ymax = proj.get_extent()
    sel2 = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
    # sel2 = ~np.isnan(x) & ~np.isnan(y)
    
    v = data[sel][sel2]
    x = x[sel2]
    y = y[sel2]

    c = polyfit2d(x, y, v, [deg, deg])

    # Evaluate the polynomial
    x, y = proj.ang2xy(lon, lat, lonlat=True)
    bkg = polynomial.polyval2d(x, y, c)
    bkg = np.ma.array(bkg, mask=data.mask, fill_value=np.nan)

    return bkg


def prepare_data(mu, hpxcube, modulus, fracdet, fracmin=0.5, clip=100, sigma=0.1, percent=[2, 95], **mask_kw):
    i = np.argmin(np.abs(mu - modulus))
    hpxmap = np.copy(hpxcube[:, i])

    fracdet = np.ones(hpxmap.size)
    fracdet[hpxmap == hp.UNSEEN] = 0

    data = streamlib.prepare_data(hpxmap, fracdet, fracmin=fracmin, clip=clip, mask_kw=mask_kw)
    # data = np.ma.array(hpxmap, mask=fracdet < 0.5)
    # data = np.ma.array(hpxmap, mask=hpxmap == hp.UNSEEN)
    data.fill_value = 0
    # data = hpxmap

    # vmin, vmax = np.percentile(data.compressed(), q=percent)
    # data = np.clip(data, vmin, vmax)

    # data.fill_value = np.ma.median(data)
    
    # print(data.fill_value)
    # print(np.unique(data[~data.mask][~np.isnan(data[~data.mask])]))

    # data = np.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False))
    data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False), mask=data.mask)
    # data = np.ma.array(hp.smoothing(data, sigma=np.radians(sigma), verbose=False))

    data.fill_value = np.ma.median(data)

    return data


def get_euler(stream=None, ends=None):
    if ends is None:
        if isinstance(stream, str):
            mw_streams = galstreams.MWStreams(verbose=False)
            stream = mw_streams[stream]
        # ends = [[stream.end_o.ra.deg, stream.end_o.dec.deg], [stream.end_f.ra.deg, stream.end_f.dec.deg]]
        ends = [[stream.end_f.ra.deg, stream.end_f.dec.deg], [stream.end_o.ra.deg, stream.end_o.dec.deg]]

    fr = gc.GreatCircleICRSFrame.from_endpoints(coord.SkyCoord(ends[0][0], ends[0][1], unit='deg'), coord.SkyCoord(ends[1][0], ends[1][1], unit='deg'))
    phi, theta, psi = get_rot(fr)
    # phi, theta, psi = zyx_euler_from_endpoints(ends[0][0], ends[0][1], ends[1][0], ends[1][1])
    print(phi, theta, psi)
    return phi, theta, psi


def get_gnomonic_proj(stream=None, ends=None, xsize=1024, ysize=800, reso=3.0):
    if ends is not None:
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


def plot_proj(proj, data, q=[5, 90], vmin=None, vmax=None):
    nside = hp.npix2nside(len(data))
    func = lambda x, y, z: hp.vec2pix(nside, x, y, z)

    img = proj.projmap(data, func)

    if vmin is None or vmax is None:
        vmin, vmax = get_vscale(img, q=q)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.imshow(img, origin='bottom', vmin=vmin, vmax=vmax, cmap='Greys', extent=proj.get_extent())

    return ax
