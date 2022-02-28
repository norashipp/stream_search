import glob
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits as fitsio

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit

DATA_DIR = '/data/des81.b/data/tavangar/skim_y6_gold_1_1/'

MAG = 'SOF_PSF_MAG_CORRECTED_%s'
MAG_ERR = 'SOF_PSF_MAG_ERR_%s'


def load_data(nfiles=1, saving=False, loading=False):
    if loading:
        g_mag = np.load('g_mag.npy', g_mag)
        g_mag_err = np.load('g_mag_err.npy', g_mag_err)
        r_mag = np.load('r_mag.npy', r_mag)
        r_mag_err = np.load('r_mag_err.npy', r_mag_err)
        i_mag = np.load('i_mag.npy', i_mag)
        i_mag_err = np.load('i_mag_err.npy', i_mag_err)
        return (g_mag, g_mag_err), (r_mag, r_mag_err), (i_mag, i_mag_err)

    fnames = glob.glob(DATA_DIR + '*.fits')
    fnames.sort()
    idx = np.random.choice(np.arange(len(fnames)), size=nfiles, replace=False)
    fnames = [fnames[i] for i in idx]

    g_mag = []
    g_mag_err = []
    r_mag = []
    r_mag_err = []
    i_mag = []
    i_mag_err = []
    # z_mag = []
    # z_mag_err = []

    for f in fnames:
        dat = fitsio.open(f)[1].data

        g_mag.append(dat[MAG % 'G'])
        g_mag_err.append(dat[MAG % 'G'])
        r_mag.append(dat[MAG % 'R'])
        r_mag_err.append(dat[MAG % 'R'])
        i_mag.append(dat[MAG % 'I'])
        i_mag_err.append(dat[MAG % 'I'])
        # z_mag.append(dat[MAG % 'Z'])
        # z_mag_err.append(dat[MAG % 'Z'])

    if saving:
        np.save('g_mag.npy', g_mag)
        np.save('g_mag_err.npy', g_mag_err)
        np.save('r_mag.npy', r_mag)
        np.save('r_mag_err.npy', r_mag_err)
        np.save('i_mag.npy', i_mag)
        np.save('i_mag_err.npy', i_mag_err)


    return (np.hstack(g_mag), np.hstack(g_mag_err)), (np.hstack(r_mag), np.hstack(r_mag_err)), (np.hstack(i_mag), np.hstack(i_mag_err)) # , (np.hstack(z_mag), np.hstack(z_mag_err))


def fit_error(nfiles, mags=None, plotting=True):
    try:
        g_mag, r_mag, i_mag = mags
    except:
        g_mag, r_mag, i_mag = load_data(nfiles=nfiles)


    bins = np.linspace(16, 25, 50)
    bin_mids = (bins[1:] + bins[:-1])/2.
    func = lambda x, a, b, c: a + np.exp((x - b)/c) 

    # G-BAND
    g_meds, _, _ = binned_statistic(g_mag[0], g_mag[1], statistic='median', bins=bins)
    (ag, bg, cg), _ = curve_fit(func, bin_mids, g_meds, (0.001, 27.09, 1.09))
    g_func = lambda g: func(g, ag, bg, cg)

    # R-BAND
    r_meds, _, _ = binned_statistic(r_mag[0], r_mag[1], statistic='median', bins=bins)
    (ar, br, cr), _ = curve_fit(func, bin_mids, r_meds, (0.001, 27.09, 1.09))
    r_func = lambda r: func(r, ar, br, cr)

    # I-BAND
    i_meds, _, _ = binned_statistic(i_mag[0], i_mag[1], statistic='median', bins=bins)
    (ai, bi, ci), _ = curve_fit(func, bin_mids, i_meds, (0.001, 27.09, 1.09))
    i_func = lambda g: func(i, ai, bi, ci)

    # Z-BAND
    # z_meds, _, _ = binned_statistic(z_mag[0], z_mag[1], statistic='median', bins=bins)
    # (a0, b0, c0), _ = curve_fit(func, bin_mids, z_meds, (0.001, 27.09, 1.09))
    # z_func = lambda g: func(z, az, bz, cz)

    if plotting:
        plot_fit(bin_mids, [g_meds, r_meds, i_meds], [g_mag, r_mag, i_mag], [g_func, r_func, i_func])

    return (ag, bg, cg), (ar, br, cr), (ai, bi, ci)


def plot_fit(bins, meds, mags, funcs):
    print('Plotting')
    for i in range(len(funcs)):
        plt.figure()
        plt.scatter(mags[i][0], mags[i][1], s=1, alpha=0.2)
        plt.plot(bins, meds, 'k--', lw=2)
        plt.plot(bins, funcs[i](bins), 'b-', lw=3)



if __name__ == '__main__':
    nfiles = 1
    print('Loading data')
    g_mag, r_mag, i_mag = load_data(nfiles, saveing=True)
    print('Fitting data')
    out = fit_error(nfiles, [g_mag, r_mag, i_mag], plotting=True)






