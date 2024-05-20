import numpy as np
import healpy as hp
import fitsio
import os
import healpy as hp

import streamlib

from utils import mod2dist, dist2mod, apwnorm


def prepare_data(mu, hpxcube, modulus, fracdet, fracmin=0.5, clip=100, sigma=0.1, percent=[2, 95], **mask_kw):
    i = np.argmin(np.abs(mu - modulus))
    hpxmap = np.copy(hpxcube[:, i])

    if fracdet is None:
        # fracdet = np.ones(hpxmap.size)
        # fracdet[hpxmap == hp.UNSEEN] = 0
        fracdet = np.zeros_like(hpxcube[:, 0])
        fracdet[np.where(np.sum(hpxcube, axis=1) > 0)] = 1

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


def fit_bkg(data, proj, sigma=0.1, percent=[2, 95], deg=5):
    nside = hp.get_nside(data.mask)
    lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)

    vmin, vmax = np.percentile(data.compressed(), q=percent)
    data = np.clip(data, vmin, vmax)
    data.fill_value = np.ma.median(data)

    smoothed = hp.smoothing(data, sigma=np.radians(sigma), verbose=False)
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



def get_survey_data(hpxcube, fracdet, modulus, survey, proj=None, slices=None, use_bkg=True, printing=False, version=0, gmax=23., rerun=False):
    if slices is None:
        slices = [(26, 50),
                (15, 26),
                (0, 15)]

    slice_str = '_'.join(['{}-{}'.format(x, y) for x, y in slices])
    dist = mod2dist(modulus)
    for x, y in slices:
        print(dist[x], dist[y])

    ###

    data_path = '../data/data_stack_%s_v%i_g%.1f_bkg%i.npy' % (survey.lower(), version, gmax, use_bkg)
    if os.path.exists(data_path) and not rerun:
        data_stack = np.load(data_path)
    else:
        print('Rerunning data stack...')
        data_stack = []

        for i in range(np.min(slices), np.max(slices)):
            mu = modulus[i]
            if printing:
                print('%.1f' %mu)

            data = prepare_data(mu, hpxcube, modulus, None, sigma=0.2, fracmin=0.5, clip=100, all_off=True, lmc=False, milky_way=False, sgr=False,
                                                    globulars=False, dwarfs=False, galaxies=False, plane=False, bmax=25, center=False, globs_dwarfs=False, acs=False)
            if use_bkg:
                nside = hp.npix2nside(len(hpxcube[:, 0]))
                func = lambda x, y, z: hp.vec2pix(nside, x, y, z)

                data_masked = prepare_data(mu, hpxcube, modulus, None, sigma=0.2, fracmin=0.5, clip=100, sgr=True, acs=True, globs_dwarfs=True, lmc=True)
                bkg = plot_density_healpy.fit_bkg(data_masked, proj, sigma=0.2)
                bkg.mask = data.mask
            else:
                bkg = 0

            data_stack.append(data - bkg)

        data_stack = np.asarray(data_stack)

        print('Saving data stack %s...' %data_path)
        np.save(data_path, data_stack.T)

    ###
    
    if data_stack.shape[1] > data_stack.shape[0]:
        data_stack = data_stack.T

    mask = streamlib.make_mask(nside=512, all_off=True, lmc=False, sgr=False, acs=False, globs_dwarfs=False, dwarfs=False, globulars=False, milky_way=False)
    mask |= fracdet < 1.0
    mask_stack = np.vstack([mask]*data_stack.shape[1]).T

    data_stack_masked = np.ma.array(data_stack, mask=mask_stack)
    data_stack_masked.fill_value = np.nan

    ###

    rgb = [np.sum(data_stack_masked[:, x:y], axis=1) / (y - x) for x, y in slices]
    for s in rgb:
        s.fill_value = np.nan

    return rgb, data_stack_masked
    

def renorm_rgb(rgb, pmin=1, pmax=99):
    for i in range(3):
        X = rgb[i]

        amin = np.nanpercentile(np.ma.filled(X, np.nan), pmin)
        amax = np.nanpercentile(np.ma.filled(X, np.nan), pmax)

        rgb[i] = apwnorm(X, min=amin, max=amax)
        rgb[i][X == 0] = 0.

    return rgb