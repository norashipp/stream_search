from __future__ import division

import numpy as np
from matplotlib.path import Path

from ugali.analysis.isochrone import factory as isochrone_factory

import surveys


###########
# CMD CUT #
###########

def mkpol(mu, age=12., z=0.0004, dmu=0.5, C=[0.05, 0.05], E=4., err=None):
    if err == None:
        print 'Using PS1 err!'
        err = surveys.surveys['PS1']['err']
    """ Builds ordered polygon for masking """

    # if survey == 'DES':
    #     err = lambda x: 0.0010908679647672335 + np.exp((x - 27.091072029215375) / 1.0904624484538419)  # median
    # elif survey == 'PS1':
    #     err = lambda x: 0.00363355415 + np.exp((x - 23.9127145) / 1.09685211)
    # err changes between surveys!

    # iso = ic.isochrone_factory('Dotter', age=age, distance_modulus=mu, z=z, dirname='/home/s1/nshipp/.ugali/isochrones/des/dotter2008')
    try:
        iso = isochrone_factory('Dotter', age=age, distance_modulus=mu, z=z, dirname='/home/s1/kadrlica/.ugali/isochrones/ps1/dotter2008')
    except:
        iso = isochrone_factory('Dotter2008', age=age, distance_modulus=mu, z=z, dirname='/Users/nora/projects/proper_motions/data/ps1')
    c = iso.color
    m = iso.mag
    mnear = m + mu - dmu / 2.
    mfar = m + mu + dmu / 2.
    C = np.r_[c + E * err(mfar) + C[1], c[::-1] - E * err(mnear[::-1]) - C[0]]
    M = np.r_[m, m[::-1]]
    return np.c_[C, M]


def select_isochrone(mag_g, mag_r, err, iso_params=[17.0, 12.5, 0.0001], dmu=0.5, C=[0.01, 0.01], E=2, gmin=None):
    mu, age, z = iso_params

    mk = mkpol(mu=mu, age=age, z=z, dmu=dmu, C=C, E=E, err=err)
    pth = Path(mk)
    cm = np.vstack([mag_g - mag_r, mag_g - mu]).T
    idx = pth.contains_points(cm)
    if gmin:
        idx &= (mag_g > gmin)
    return idx


def calculate_err():
    # binned statistic median
    # curve fit func
    pass


##################
# METAL POOR CUT #
##################

def calculate_locus(gr, ri):
    out = binned_statistic(gr, ri, 'median', bins=100, range=[0.2, 0.8])
    xx = (out[1][:-1] + out[1][1:]) / 2.
    yy = out[0]

    p = np.polyfit(xx, yy, deg=1)
    return p


def locus(mag_g, mag_r, mag_i, p=None):
    if np.any(p == None):
        p = calculate_locus(mag_g - mag_r, mag_r - mag_i)
    return p[0] * (mag_g - mag_r) + p[1]


def select_metal_poor(mag_g, mag_r, mag_i, survey='DES_Y3A2', dmin=0.02, dmax=0.06, p=None):
    stellar_locus = locus(data, p=p)
    idx = (mag_r - mag_i > stellar_locus + dmin) & (mag_r - mag_i < stellar_locus + dmax)
    return data[idx], idx
