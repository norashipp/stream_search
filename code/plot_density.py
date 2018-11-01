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


def load_hpxcube():
    pass


def plot_density():
    pass


def make_movie():
    pass
