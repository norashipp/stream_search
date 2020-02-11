import os

for survey in ['BASS', 'DECaLS']:
	for z in [0.0001, 0.0004, 0.001]:
		print('Running python makepxcube.py -f ../data/%s_iso_hpxcube_z%.4f.fits.gz -mp 62 -a 12 -z %.4f -s %s...' %(survey, z, z, survey))
		os.system('python makepxcube.py -f ../data/%s_iso_hpxcube_z%.4f.fits.gz -mp 62 -a 12 -z %.4f -s %s' %(survey, z, z, survey))