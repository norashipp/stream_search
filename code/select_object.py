import numpy as np

from astropy import table

from ugali.utils import healpix
from ugali.utils import load_infiles
from utils_alex import load_infiles

import streamlib


def get_object(ra, dec, survey='DECaLS'):
    if survey == 'DECaLS':
        filename = '/data/des40.b/data/decals/dr8/south_skim/decals-dr8-sweep_%0.5d.fits'
    elif survey == 'BASS':
        filename = '/data/des40.b/data/decals/dr8/north_skim/decals-dr8-sweep_%0.5d.fits'

    filenames = [filename % i for i in healpix.ang2disc(32, ra, dec, 5)]
    columns = ['RA', 'DEC', 'MAG_SFD_G', 'MAG_SFD_R', 'MAG_SFD_Z', 'EXTENDED_CLASS']

    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['EXTENDED_CLASS'] == 0) & (data['MAG_SFD_G'] < 23.5) & (data['MAG_SFD_G'] > 16) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] > 0) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] < 1)
    sep = angsep(ra, dec, data['RA'], data['DEC'])
    sel &= (sep < 0.1)

    tab = table.Table(data[sel])
    tab.write(filename)


def get_stream(ends, survey='DECaLS', filename='cutout.fits'):
    if survey == 'DECaLS':
        filename = '/data/des40.b/data/decals/dr8/south_skim/decals-dr8-sweep_%0.5d.fits'
    elif survey == 'BASS':
        filename = '/data/des40.b/data/decals/dr8/north_skim/decals-dr8-sweep_%0.5d.fits'

    length = angsep(ends[0][0], ends[0][1], ends[1][0], ends[1][1])

    filenames = [filename % i for i in healpix.ang2disc(32, np.mean([ends[0][0], ends[1][0]]), np.mean([ends[0][1], ends[1][1]]), length * 1.5)]
    columns = ['RA', 'DEC', 'MAG_SFD_G', 'MAG_SFD_R', 'MAG_SFD_Z', 'EXTENDED_CLASS']

    data = load_infiles(filenames, columns=columns, multiproc=32)

    sel = (data['EXTENDED_CLASS'] == 0) & (data['MAG_SFD_G'] < 23.5) & (data['MAG_SFD_G'] > 16) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] > 0) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] < 1)

    R = streamlib.get_rotmat(ends=ends)
    phi1, phi2 = phi12_rotmat(data['RA'], data['DEC'], R)

    sel &= (np.abs(phi1) < (length / 2. + 5)) & (np.abs(phi2) < 5)

    tab = table.Table(data[sel])
    tab.write(filename)


if __name__ == '__main__':
    stream = 'ATLAS'
    ends = [[19.465113557599395, -26.584615187212712], [31.04372386479431, -32.98118501241838]]
    get_stream(ends, survey='DECaLS', filename='/data/des40.b/data/nshipp/stream_search/data/cutouts/ATLAS_cutout.fits')
