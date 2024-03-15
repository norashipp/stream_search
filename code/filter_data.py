from __future__ import division

import numpy as np
from matplotlib.path import Path

from ugali.analysis.isochrone import factory as isochrone_factory

import surveys


###########
# CMD CUT #
###########

HOMEDIR = '/home/s1/nshipp/'
# HOMEDIR = '/Users/nora/'


def mkpol(mu, age=12., z=0.0004, dmu=0.5, C=[0.05, 0.05], E=4., err=None, survey='DECaLS', clip=None):
    if err == None:
        try:
            err = surveys.surveys[survey]['err']
        except:
            print('Using DES err!')
            err = surveys.surveys['DES_DR1']['err']
    """ Builds ordered polygon for masking """

    # if survey == 'DES':
    #     err = lambda x: 0.0010908679647672335 + np.exp((x - 27.091072029215375) / 1.0904624484538419)  # median
    # elif survey == 'PS1':
    #     err = lambda x: 0.00363355415 + np.exp((x - 23.9127145) / 1.09685211)
    # err changes between surveys!

    # iso = ic.isochrone_factory('Dotter', age=age, distance_modulus=mu, z=z, dirname='/home/s1/nshipp/.ugali/isochrones/des/dotter2008')
    
    if survey in ['PS1']:
        iso = isochrone_factory(
            'Dotter', survey='ps1', age=age, distance_modulus=mu, z=z)
    elif survey in ['DES_DR1', 'DES_Y3A2', 'DECaLS', 'DES_Y6', 'DES_Y6_GOLD', 'DELVE', 'DECaLS_DR9', 'DECaLS_DR10', 'DECaLS_DR10_OFF', 'DELVE_R1', 'DELVE_R2', 'DELVE_DR3']:
        iso = isochrone_factory('Dotter', survey='des',
                                age=age, distance_modulus=mu, z=z)
    elif survey in ['BASS', 'BASS_DR9']:
        try:
            print(
                HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
            iso = np.loadtxt(
                HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
        except:
            print('Error loading BASS isochrone...')
            print(mu, age, z)

        g_ps1 = iso[:, 9]
        r_ps1 = iso[:, 10]
        i_ps1 = iso[:, 11]

        g_bass = g_ps1 + 0.00464 + 0.08672 * \
            (g_ps1 - i_ps1) - 0.00668 * (g_ps1 - i_ps1) ** 2 - \
            0.00255 * (g_ps1 - i_ps1) ** 3
        r_bass = r_ps1 + 0.00110 - 0.06875 * \
            (g_ps1 - i_ps1) + 0.02480 * (g_ps1 - i_ps1)**2 - \
            0.00855 * (g_ps1 - i_ps1)**3

    else:
        print('Survey error - update isochrones.')

    if survey in ['BASS', 'BASS_DR9']:
        c = g_bass - r_bass
        m = g_bass
    else:
        c = iso.color
        m = iso.mag

    if clip is not None:
        # Clip for plotting, use gmin otherwise
        # clip abs mag
        cut = (m > clip) & ((m + mu) < 23.0) & (c > 0) & (c < 1)
        c = c[cut]
        m = m[cut]

    mnear = m + mu - dmu / 2.
    mfar = m + mu + dmu / 2.
    C = np.r_[c + E * err(mfar) + C[1], c[::-1] - E * err(mnear[::-1]) - C[0]]
    M = np.r_[m, m[::-1]]
    return np.c_[C, M]


def select_isochrone(mag_g, mag_r, err, iso_params=[17.0, 12.5, 0.0001], dmu=0.5, C=[0.01, 0.01], E=2, gmin=None, survey='DECaLS'):
    mu, age, z = iso_params

    mk = mkpol(mu=mu, age=age, z=z, dmu=dmu, C=C, E=E, err=err, survey=survey)
    pth = Path(mk)
    cm = np.vstack([mag_g - mag_r, mag_g - mu]).T
    idx = pth.contains_points(cm)
    if gmin:
        idx &= (mag_g > gmin)
    return idx


def mkpol_grz(mu, age=12., z=0.0004, dmu=0.5, C=[0.05, 0.05], E=4., err=None, survey='DECaLS', clip=None):
    if err == None:
        print('Using DES err!')
        err = surveys.surveys['DES_DR1']['err']

    if survey in ['PS1']:
        iso = np.loadtxt(
            HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))

        g_ps1 = iso[:, 9]
        r_ps1 = iso[:, 10]
        # i_ps1 = iso[:, 11]
        z_ps1 = iso[:, 12]

        c = r_ps1 - z_ps1
        m = g_ps1

    elif survey in ['DES_DR1', 'DES_Y3A2', 'DECaLS', 'DES_Y6', 'DECaLS_DR10', 'DECaLS_DR10_OFF', 'DELVE_DR3']:
        iso = np.loadtxt(
            HOMEDIR + '.ugali/isochrones/des/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))

        # u_decam = iso[:, 9]
        g_decam = iso[:, 10]
        r_decam = iso[:, 11]
        # i_decam = iso[:, 12]
        z_decam = iso[:, 13]

        c = r_decam - z_decam
        m = g_decam

    elif survey in ['BASS', 'BASS_DR9']:
        try:
            print(
                HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
            iso = np.loadtxt(
                HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
        except:
            print('Error loading BASS isochrone...')
            print(mu, age, z)

        g_ps1 = iso[:, 9]
        r_ps1 = iso[:, 10]
        i_ps1 = iso[:, 11]
        z_ps1 = iso[:, 12]

        g_bass = g_ps1 + 0.00464 + 0.08672 * \
            (g_ps1 - i_ps1) - 0.00668 * (g_ps1 - i_ps1) ** 2 - \
            0.00255 * (g_ps1 - i_ps1) ** 3
        r_bass = r_ps1 + 0.00110 - 0.06875 * \
            (g_ps1 - i_ps1) + 0.02480 * (g_ps1 - i_ps1)**2 - \
            0.00855 * (g_ps1 - i_ps1)**3
        z_bass = z_ps1 + 0.03664 - 0.11084 * \
            (g_ps1 - i_ps1) + 0.04477 * (g_ps1 - i_ps1)**2 - \
            0.01223 * (g_ps1 - i_ps1)**3

        c = r_bass - z_bass
        m = g_bass

    else:
        print('Survey error - update isochrones.')

    if clip is not None:
        # Clip for plotting, use gmin otherwise
        cut = (m > clip) & ((m + mu) < 23.0) & (c > 0) & (c < 1)
        c = c[cut]
        m = m[cut]

    mnear = m + mu - dmu / 2.
    mfar = m + mu + dmu / 2.
    C = np.r_[c + E * err(mfar) + C[1], c[::-1] - E * err(mnear[::-1]) - C[0]]
    M = np.r_[m, m[::-1]]
    return np.c_[C, M]


def select_isochrone_grz(mag_g, mag_r, mag_z, err, iso_params=[17.0, 12.5, 0.0001], dmu=0.5, C=[0.01, 0.01], E=2, gmin=None, survey='DECaLS'):
    if survey in ['BASS', 'BASS_DR9']:
        C = [0.1, 0.05]

    mu, age, z = iso_params

    mk = mkpol_grz(mu=mu, age=age, z=z, dmu=dmu,
                   C=C, E=E, err=err, survey=survey)
    pth = Path(mk)
    cm = np.vstack([mag_r - mag_z, mag_g - mu]).T
    idx = pth.contains_points(cm)
    if gmin:
        idx &= (mag_g > gmin)
    return idx

def mkpol_gi(mu, age=12., z=0.0004, dmu=0.5, C=[0.05, 0.05], E=4., err=None, survey='DECaLS', clip=None):
    if err == None:
        print('Using DES err!')
        err = surveys.surveys['DES_DR1']['err']

    if survey in ['PS1']:
        iso = np.loadtxt(
            HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))

        g_ps1 = iso[:, 9]
        # r_ps1 = iso[:, 10]
        i_ps1 = iso[:, 11]
        # z_ps1 = iso[:, 12]

        c = g_ps1 - i_ps1
        m = g_ps1

    elif survey in ['DES_DR1', 'DES_Y3A2', 'DECaLS', 'DES_Y6', 'DES_Y6_GOLD', 'DELVE_R1', 'DELVE_R2', 'DELVE', 'DELVE_DR3']:
        iso = np.loadtxt(
            HOMEDIR + '.ugali/isochrones/des/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))

        # u_decam = iso[:, 9]
        g_decam = iso[:, 10]
        # r_decam = iso[:, 11]
        i_decam = iso[:, 12]
        # z_decam = iso[:, 13]

        c = g_decam - i_decam
        m = g_decam

    # elif survey in ['BASS', 'BASS_DR9']:
    #     print('Survey error - use g, r for BASS')
    #     try:
    #         print(
    #             HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
    #         iso = np.loadtxt(
    #             HOMEDIR + '.ugali/isochrones/ps1/dotter2016/iso_a%.1f_z%.5f.dat' % (age, z))
    #     except:
    #         print('Error loading BASS isochrone...')
    #         print(mu, age, z)

    #     g_ps1 = iso[:, 9]
    #     r_ps1 = iso[:, 10]
    #     i_ps1 = iso[:, 11]
    #     z_ps1 = iso[:, 12]

    #     g_bass = g_ps1 + 0.00464 + 0.08672 * \
    #         (g_ps1 - i_ps1) - 0.00668 * (g_ps1 - i_ps1) ** 2 - \
    #         0.00255 * (g_ps1 - i_ps1) ** 3
    #     r_bass = r_ps1 + 0.00110 - 0.06875 * \
    #         (g_ps1 - i_ps1) + 0.02480 * (g_ps1 - i_ps1)**2 - \
    #         0.00855 * (g_ps1 - i_ps1)**3
    #     z_bass = z_ps1 + 0.03664 - 0.11084 * \
    #         (g_ps1 - i_ps1) + 0.04477 * (g_ps1 - i_ps1)**2 - \
    #         0.01223 * (g_ps1 - i_ps1)**3

    #     c = r_bass - z_bass
    #     m = g_bass

    else:
        print('Survey error - update isochrones.')

    if clip is not None:
        # Clip for plotting, use gmin otherwise
        cut = (m > clip) & ((m + mu) < 23.0) & (c > 0) & (c < 1)
        c = c[cut]
        m = m[cut]

    mnear = m + mu - dmu / 2.
    mfar = m + mu + dmu / 2.
    C = np.r_[c + E * err(mfar) + C[1], c[::-1] - E * err(mnear[::-1]) - C[0]]
    M = np.r_[m, m[::-1]]
    return np.c_[C, M]


def select_isochrone_gi(mag_g, mag_i, err, iso_params=[17.0, 12.5, 0.0001], dmu=0.5, C=[0.01, 0.01], E=2, gmin=None, survey='DECaLS'):
    # if survey in ['BASS', 'BASS_DR9']:
        # C = [0.1, 0.05]

    mu, age, z = iso_params

    mk = mkpol_gi(mu=mu, age=age, z=z, dmu=dmu,
                   C=C, E=E, err=err, survey=survey)
    pth = Path(mk)
    cm = np.vstack([mag_g - mag_i, mag_g - mu]).T
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
    idx = (mag_r - mag_i > stellar_locus +
           dmin) & (mag_r - mag_i < stellar_locus + dmax)
    return data[idx], idx
