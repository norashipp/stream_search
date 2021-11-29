import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import astropy.io.fits as fitsio
from astropy.coordinates import SkyCoord
import astropy.units as u
import galstreams
import healpy as hp

import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

import warnings
warnings.filterwarnings('ignore')

from importlib import reload

import sys
# sys.path.append('/Users/nora/projects/stream_search/code')
sys.path.append('../code')

# import plot_density
# reload(plot_density)
import plot_density_healpy
reload(plot_density_healpy)
import streamlib
reload(streamlib)
import load_data
reload(load_data)


plot_density_healpy.plot_pretty(figsize=(12, 12))

