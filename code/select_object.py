import numpy as np
import os
from astropy import table

from ugali.utils import healpix
from ugali.utils.projector import angsep
from utils_alex import load_infiles

# import streamlib
from streamlib_new import get_rotmat
from rotation_matrix import phi12_rotmat

# import galstreams


def get_object(ra, dec, survey='DECaLS', outfile='cutout.fits', radius=0.1):
    if survey == 'DECaLS':
        filename = '/data/des40.b/data/decals/dr8/south_skim/decals-dr8-sweep_%0.5d.fits'
    elif survey == 'BASS':
        filename = '/data/des40.b/data/decals/dr8/north_skim/decals-dr8-sweep_%0.5d.fits'

    filenames = [filename % i for i in healpix.ang2disc(32, ra, dec, np.maximum(5, radius * 2)) if os.path.exists(filename % i)]
    columns = ['RA', 'DEC', 'MAG_G', 'MAG_R', 'MAG_Z', 'MAG_SFD_G', 'MAG_SFD_R', 'MAG_SFD_Z', 'EXTENDED_CLASS']

    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['EXTENDED_CLASS'] == 0) & (data['MAG_SFD_G'] < 23.5) & (data['MAG_SFD_G'] > 16) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] > 0) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] < 1)
    sep = angsep(ra, dec, data['RA'], data['DEC'])
    sel &= (sep < radius)

    tab = table.Table(data[sel])
    tab.write(outfile)


def get_object_desy6(ra, dec, outfile='cutout.fits', radius=0.1):
    filename = '/data/des81.b/data/mmcnanna/y6a1/skim_y6_gold_1_1/y6_gold_1_0_%0.5d.fits'

    filenames = [filename % i for i in healpix.ang2disc(32, ra, dec, np.maximum(5, radius * 2)) if os.path.exists(filename % i)]
    columns = ['RA', 'DEC', 'SOF_PSF_MAG_CORRECTED_R', 'SOF_PSF_MAG_CORRECTED_R', 'EXT_SOF']

    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['EXT_SOF'] <= 1) & (data['SOF_PSF_MAG_CORRECTED_G'] < 24.5) & (data['SOF_PSF_MAG_CORRECTED_G'] > 16) & (
        data['SOF_PSF_MAG_CORRECTED_G'] - data['SOF_PSF_MAG_CORRECTED_R'] > 0) & (data['SOF_PSF_MAG_CORRECTED_G'] - data['SOF_PSF_MAG_CORRECTED_R'] < 1)
    sep = angsep(ra, dec, data['RA'], data['DEC'])
    sel &= (sep < radius)

    tab = table.Table(data[sel])
    tab.write(outfile)


def get_object_ps1(ra, dec, outfile='cutout.fits', radius=0.1):
    filename = '/home/s1/kadrlica/projects/ps1/dr1/v0/skim_ext_0_0/ps1_dr1_%0.5d.fits'

    filenames = [filename % i for i in healpix.ang2disc(32, ra, dec, np.maximum(5, radius * 2)) if os.path.exists(filename % i)]
    columns = ['RA', 'DEC', 'GFPSFMAG_SFD', 'RFPSFMAG_SFD']

    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['GFPSFMAG_SFD'] < 21.5) & (data['GFPSFMAG_SFD'] > 14) & (data['GFPSFMAG_SFD'] - data['RFPSFMAG_SFD'] > 0) & (data['GFPSFMAG_SFD'] - data['RFPSFMAG_SFD'] < 1)
    sep = angsep(ra, dec, data['RA'], data['DEC'])
    sel &= (sep < radius)

    tab = table.Table(data[sel])
    tab.write(outfile)


def get_stream(ends, survey='DECaLS', outfile='cutout.fits'):
    if survey == 'DECaLS':
        filename = '/data/des40.b/data/decals/dr8/south_skim/decals-dr8-sweep_%0.5d.fits'
    elif survey == 'BASS':
        filename = '/data/des40.b/data/decals/dr8/north_skim/decals-dr8-sweep_%0.5d.fits'

    length = angsep(ends[0][0], ends[0][1], ends[1][0], ends[1][1])

    filenames = [filename % i for i in healpix.ang2disc(32, np.mean([ends[0][0], ends[1][0]]), np.mean([ends[0][1], ends[1][1]]), length * 1.5) if os.path.exists(filename % i)]
    columns = ['RA', 'DEC', 'MAG_G', 'MAG_R', 'MAG_Z', 'MAG_SFD_G', 'MAG_SFD_R', 'MAG_SFD_Z', 'MAG_G', 'MAG_R', 'MAG_Z', 'EXTENDED_CLASS']

    print('Loading data...')
    data = load_infiles(filenames, columns=columns, multiproc=32)

    print('Selecting data...')
    sel = (data['EXTENDED_CLASS'] == 0) & (data['MAG_SFD_G'] < 23.5) & (data['MAG_SFD_G'] > 16) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] > 0) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] < 1)
    data = data[sel]

    print('Converting coordinates...')
    R = get_rotmat(ends=ends)
    phi1, phi2 = phi12_rotmat(data['RA'], data['DEC'], R)

    sel = (np.abs(phi1) < (length / 2. + 5)) & (np.abs(phi2) < 5)

    print('Writing data...')
    tab = table.Table(data[sel])
    tab.write(outfile)


def get_stream_desy6(ends, outfile='cutout.fits'):
    filename = '/data/des81.b/data/mmcnanna/y6a1/skim_y6_gold_1_1/y6_gold_1_0_%0.5d.fits'

    length = angsep(ends[0][0], ends[0][1], ends[1][0], ends[1][1])

    filenames = [filename % i for i in healpix.ang2disc(32, np.mean([ends[0][0], ends[1][0]]), np.mean([ends[0][1], ends[1][1]]), length * 1.5) if os.path.exists(filename % i)]
    columns = ['RA', 'DEC', 'SOF_PSF_MAG_CORRECTED_G', 'SOF_PSF_MAG_CORRECTED_R', 'EXT_SOF']

    print('Loading data...')
    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['EXT_SOF'] <= 1) & (data['SOF_PSF_MAG_CORRECTED_G'] < 24.5) & (data['SOF_PSF_MAG_CORRECTED_G'] > 16) & (
        data['SOF_PSF_MAG_CORRECTED_G'] - data['SOF_PSF_MAG_CORRECTED_R'] > 0) & (data['SOF_PSF_MAG_CORRECTED_G'] - data['SOF_PSF_MAG_CORRECTED_R'] < 1)

    data = data[sel]

    print('Converting coordinates...')
    R = get_rotmat(ends=ends)
    phi1, phi2 = phi12_rotmat(data['RA'], data['DEC'], R)
    sel = (np.abs(phi1) < (length / 2. + 5)) & (np.abs(phi2) < 5)

    print('Writing file %s...' %outfile)
    tab = table.Table(data[sel])
    tab.write(outfile)



if __name__ == '__main__':
    # stream = sys.argv[1]
    # mw_streams = galstreams.MWStreams(verbose=False)
    # stream = mw_streams[stream]
    # ends = [[19.465113557599395, -26.584615187212712], [31.04372386479431, -32.98118501241838]]
    # ends = [[stream.end_f.ra.deg, stream.end_f.dec.deg], [stream.end_o.ra.deg, stream.end_o.dec.deg]]
    # ends = [[31.04372386, -32.98118501], [20.082460505880235, -56.996486198871246]]
    # get_stream(ends, survey='DECaLS', outfile='/data/des40.b/data/nshipp/stream_search/data/cutouts/Phoenix_cutout.fits')
    # ends = [[-15.565158278830129, 9.145015179988334], [-9.804650620987337, 17.42797474620764]]
    # get_stream(ends, survey='DECaLS', outfile='/data/des40.b/data/nshipp/stream_search/data/cutouts/Pal13_2.fits.gz')
    # ra_ngc1851, dec_ngc1851 =  78.528, -40.047
    # get_object_desy6(ra_ngc1851, dec_ngc1851, radius=30, outfile='../data/NGC1851_DES_Y6.fits.gz')

    # ra_sgr, dec_sgr = 283.7629167, -30.4783333
    # get_object_ps1(ra_sgr, dec_sgr, radius=10, outfile='../data/Sgr_dsph_PS1.fits.gz')

    # ra_leo2, dec_leo2 = 168.3716667, 22.1547222
    # get_object_ps1(ra_leo2, dec_leo2, radius=10, outfile='../data/LeoII_PS1.fits.gz')

    ends = [[-44.05, -56.54], [-56.98, -50.85]]
    get_stream_desy6(ends, outfile='y6_rrl_stream.fits.gz')
