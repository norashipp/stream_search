import numpy as np
import glob

import astropy.io.fits as fitsio
from astropy import table

from utils_alex import load_infiles
from surveys import surveys
from stream_data import stream_matrices
import rotation_matrix


def load_data(survey='DECaLS_DR9', ra_range=None):
    gmin, gmax = surveys[survey]['minmag'], surveys[survey]['maxmag']
    data_dir = surveys[survey]['data_dir']
    mag = surveys[survey]['mag']
    columns = ['RA', 'DEC', mag % 'G', mag % 'R']
    filenames = glob.glob(data_dir + '/*.fits')
    if ra_range is not None:
        filenames = [f for f in filenames if (f[75:78] < ra_range[1]) and (f[83:86] > ra_range[0])]
    filenames = filenames[:5]
    data = load_infiles(filenames, columns=columns)
    data['RA'] = data['RA'] - 360 * (data['RA'] > 180)
    idx = (data[mag % 'G'] < gmax) & (data[mag % 'G'] > gmin) & (data[mag % 'R'] < gmax) & (data[mag % 'R'] > gmin)
    return data[idx]


def select_cutout(data, extent):
    data = (data['RA'] > extent[0]) & (data['RA'] < extent[1]) & (data['DEC'] > extent[2]) & (data['DEC'] < extent[3])


def select_stream_cutout(data, stream, extent):
    R = stream_matrices[stream]
    phi1, phi2 = rotation_matrix.phi12_rotmat(data['RA'], data['DEC'], R)
    data['PHI1'] = phi1
    data['PHI2'] = phi2
    data = (data['PHI1'] > extent[0]) & (data['PHI1'] < extent[1]) & (data['PHI2'] > extent[2]) & (data['PHI2'] < extent[3])
    return data


def save_data(data, filename='cutout.fits', data_dir='/data/des40.b/data/nshipp/stream_cutouts/v3.0/'):
    data = table.Table(data)
    data.write(data_dir + filename)


if __name__ == '__main__':
    data = load_data(survey='DECaLS_DR9', ra_range=[0,60])
    data = select_stream_cutout(data=data, extent=[-30, 30, -10, 10], stream='ATLAS')
    save_data(data, filename='AAU_LS9_cutout.fits')
