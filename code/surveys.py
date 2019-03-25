from collections import OrderedDict as odict

surveys = odict([
    ('DES_Y3A2',
     dict(
         mag='PSF_MAG_SFD_%s',
         ext=None,
         data_dir='/home/s1/kadrlica/projects/field_of_streams/v8/skim/',
         fracdet='/data/des40.b/data/y3a2/maps/y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz',
         minmag=16.,
         maxmag=23.5)),
    ('DES_DR1',
     dict(
         mag='WAVG_MAG_PSF_%s_DERED',
         ext=None,
         data_dir='/data/des40.b/data/nshipp/skim518/',
         fracdet='/data/des40.b/data/y3a2/maps/y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz',
         minmag=16.,
         maxmag=23.5)),
    ('PS1',
     dict(
         mag='%sFPSFMAG',  # '%sPSFMAG',
         ext=None,  # 'EXTSFD_%s',
         data_dir='/home/s1/kadrlica/projects/ps1/dr1/v0/skim_ext_0_0/',  # '/data/des40.b/data/pan-starrs/dr1/healpix/',
         fracdet=None,
         minmag=12.,
         maxmag=23.3)),
    ('SDSS_DR13',  # star galaxy separation?
     dict(
         mag='MAG_PSF_SFD_%s',
         ext=None,
         data_dir='/data/des40.b/data/sdss/dr13/healpix/',
         fracdet=None,
         minmag=16.,
         maxmag=22.0))
])
