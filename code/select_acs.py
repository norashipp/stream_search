import numpy as np
from utils import load_infiles
from astropy import table
import os


pix = np.loadtxt('../data/acs_pix.dat')
columns = ['RA', 'DEC', 'MAG_SFD_G', 'MAG_SFD_R', 'EXTENDED_CLASS']

data = load_infiles([filename % i for i in pix if os.path.exists(
    filename % i)], columns=columns, multiproc=16)
sel = (data['EXTENDED_CLASS'] == 0) & (data['MAG_SFD_G'] < 23) & (data['MAG_SFD_G'] > 16) & (
    data['MAG_SFD_G'] - data['MAG_SFD_R'] > 0) & (data['MAG_SFD_G'] - data['MAG_SFD_R'] < 1)
data = data[sel]

phi1, phi2 = rotation_matrix.phi12_rotmat(data['RA'], data['DEC'], R)
sel = (phi2 < 6) & (phi2 > -10) & (phi1 > -70)
data = data[sel]

tab = table.Table(data)
tab.write('../data/acs_data.fits_3.gz')
