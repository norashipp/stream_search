from collections import OrderedDict as odict
import numpy as np

surveys = odict([
    ('DES_Y3A2',
     dict(
         mag='PSF_MAG_SFD_%s',
         ext=None,
         data_dir='/home/s1/kadrlica/projects/field_of_streams/v8/skim/',
         fracdet='/data/des40.b/data/y3a2/maps/y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz',
         minmag=16.,
         maxmag=23.5,
         err=lambda x: 0.0010908679647672335 +
             np.exp((x - 27.091072029215375) / 1.0904624484538419),
         C=[0.05, 0.1],
         E=2.,
         moduli=[15, 20])),
    ('DES_DR1',
     dict(
         mag='WAVG_MAG_PSF_%s_DERED',
         ext=None,
         data_dir='/data/des40.b/data/nshipp/skim518/',
         fracdet='/data/des40.b/data/y3a2/maps/y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz',
         minmag=16.,
         maxmag=23.5,
         err=lambda x: 0.0010908679647672335 +
             np.exp((x - 27.091072029215375) / 1.0904624484538419),
         C=[0.05, 0.1],
         E=2.,
         moduli=[15, 20])),
    ('DES_Y6',
     dict(
         mag='SOF_PSF_MAG_CORRECTED_%s',
         ext=None,
         data_dir='/data/des81.b/data/mmcnanna/y6a1/skim_y6_gold_1_1/',
         fracdet=None,
         minmag=16.,
         maxmag=24.0,
         stargal='EXT_SOF',
         stargal_cut=1,
         err=lambda x: 0.0010908679647672335 + np.exp((x - 27.091072029215375) / 1.0904624484538419),  # Y3
         C=[0.05, 0.1],
         E=2.,
         moduli=[15, 20])),
    ('PS1',
     dict(
         mag='%sFPSFMAG_SFD',  # '%sPSFMAG',
         ext=None,  # 'EXTSFD_%s',
         # '/data/des40.b/data/pan-starrs/dr1/healpix/',
         data_dir='/home/s1/kadrlica/projects/ps1/dr1/v0/skim_ext_0_0/',
         fracdet=None,  # '/home/s1/smau/projects/panstarrs/simple_v3/panstarrs_pseudo_fracdet.fits.gz',
         minmag=14.,  # 12.,
         maxmag=21.5,
         stargal=None,
         stargal_cut=None,
         err=lambda x: 0.00363355415 + np.exp((x - 23.9127145) / 1.09685211),
         C=[0.05, 0.05],
         E=1.,
         moduli=[13, 18])),
    ('DECaLS',
     dict(
         mag='MAG_SFD_%s',
         ext=None,
         # '/data/des40.b/data/decals/dr7/skim/',
         data_dir='/data/des40.b/data/decals/dr8/south_skim/',
         # data_dir='/data/des51.b/data/tavangar/south_skim_zcut/',
         # data_dir='/data/des51.b/data/tavangar/south_skim_gal/',
         fracdet=None,  # '/data/des40.b/data/nshipp/projects/stream_search/data/decals_dr7_pseudo_fracdet.fits.gz',
         minmag=16.,
         maxmag=23.0,
         stargal='EXTENDED_CLASS',
         # stargal=None,
         stargal_cut=0,
         err=lambda x: 0.0010908679647672335 + \
             np.exp((x - 27.091072029215375) / 1.0904624484538419),
         C=[-0.04, 0.1],  # [0.01, 0.075], # [0.05, 0.1],
         E=2.,  # 3., # 4., # 2.,
         moduli=[15., 20.])),  # [14, 20])),
    ('BASS',
     dict(
         mag='MAG_SFD_%s',
         ext=None,
         data_dir='/data/des40.b/data/decals/dr8/north_skim/',
         fracdet=None,
         minmag=16.,
         maxmag=23.0,
         stargal='EXTENDED_CLASS',
         stargal_cut=0,
         err=lambda x: 6.31366550e-04 + \
             np.exp((x - 2.57279516e+01) / 1.15917318e+00),
         C=[0.05, 0.5],
         E=2.,
         moduli=[14, 20])),
    ('SDSS_DR13',  # star galaxy separation?
     dict(
         mag='MAG_PSF_SFD_%s',
         ext=None,
         data_dir='/data/des40.b/data/sdss/dr13/healpix/',
         fracdet=None,
         minmag=16.,
         maxmag=22.0))
])
