import numpy as np
import healpy as hp
import fitsio


SURVEY = ['BASS', 'DECaLS', 'DES_Y6', 'BASS_DR2', 'DECaLS_DR9', 'DES_Y6_GOLD', 'DELVE']
Z = np.array([0.0001, 0.0004, 0.0010])
AGE = np.array([11.0, 12.0, 13.0, 13.5])
GMAX = np.array([22.0, 23.0])


def load_data(survey='DECaLS', z=0.0001, age=13.5, gmax=23.0, gmin=16., filename=None):
    if filename is None:
        if survey.lower() not in [s.lower() for s in SURVEY]:
            print('Survey options are BASS or DECaLS. Using DECaLS.')
            survey = 'DECaLS'

        z = Z[np.argmin(np.abs(z - Z))]
        age = AGE[np.argmin(np.abs(age - AGE))]
        gmax = GMAX[np.argmin(np.abs(gmax - GMAX))]

        filename = '/Users/nora/projects/stream_search/data/%s_iso_hpxcube_z%.4f_a%.1f_gmax%i_gmin%i.fits.gz' % (survey, z, age, gmax, gmin)

        if survey == 'DES_Y6':
            # filename = '/Users/nora/projects/stream_search/data/iso_hpxcube_desy6.fits.gz'
            # filename = '/Users/nora/projects/stream_search/data/iso_hpxcube_desy6_phoenix.fits.gz'
            # filename = '../data/DES_Y6_iso_hpxcube_z0.0001_a11.0_gmax24_gmin3p5_E2_C0p050p1.fits.gz'
            filename = '../data/DES_Y6_iso_hpxcube_NGC1851.fits.gz'

    # print('Loading %s...' %filename)
    hpxcube, fracdet, modulus = load_hpxcube(filename)

    if fracdet is None:
        fracdet = load_fracdet(survey)
    if fracdet is None:
        fracdet = np.zeros_like(hpxcube[:,0])
        fracdet[np.where(np.sum(hpxcube, axis=1) > 0)] = 1

    return hpxcube, fracdet, modulus


def load_fracdet(survey):
    if survey == 'DECaLS':
        print('Reading %s' % '../data/DECaLS_nside128.fits.gz...')
        hpxmap128 = fitsio.read('../data/DECaLS_nside128.fits.gz')
        fracdet128 = hpxmap128 > 0
        fracdet512 = convert_fracdet(fracdet128, 512)
        return fracdet512
    elif survey == 'BASS':
        print('Reading %s' % '../data/BASS_nside256.fits.gz...')
        hpxmap128 = fitsio.read('../data/BASS_nside256.fits.gz')
        fracdet128 = hpxmap128 > 0
        fracdet512 = convert_fracdet(fracdet128, 512)
        return fracdet512
    # elif survey in ['DES_Y6', 'DES_Y6_GOLD']:
    #     print('Reading %s' % 'y6a2_g_o.4096_t.32768_frac_EQU.fits.fz')
    #     fracdet4096 = np.zeros(hp.nside2npix(4096))
    #     f = fitsio.read('../data/2021/y6a2_g_o.4096_t.32768_frac_EQU.fits.fz')
    #     fracdet4096[f['PIXEL']] = f['SIGNAL']
    #     fracdet512 = convert_fracdet(fracdet4096, 512)
    #     return fracdet512
    else:
        return None


def load_hpxcube(filename='../data/iso_hpxcube.fits.gz'):
    print("Reading %s..." % filename)
    f = fitsio.FITS(filename)
    hpxcube = f['HPXCUBE'].read()
    try:
        fracdet = f['FRACDET'].read()
        print('fracdet test', np.sum(fracdet > 0.5))
    except:
        # print('Skipping fracdet...')
        # fracdet = np.zeros_like(hpxcube[:, 0])
        # fracdet[np.where(hpxcube[:, 0] > 0)] = 1
        fracdet = None
    try:
        modulus = f['MODULUS'].read()
    except:
        print('Error reading modulus...')
        modulus = np.array([16.])
    return hpxcube, fracdet, modulus


def convert_fracdet(fracdet, nside):
    nside1 = nside
    fracdet2 = fracdet
    nside2 = hp.npix2nside(len(fracdet2))

    pix1 = np.arange(hp.nside2npix(nside1))
    ra1, dec1 = hp.pix2ang(nside1, pix1, lonlat=True)
    pix12 = hp.ang2pix(nside2, ra1, dec1, lonlat=True)

    fracdet1 = np.zeros(hp.nside2npix(nside1))
    fracdet1 = fracdet2[pix12]

    return fracdet1
