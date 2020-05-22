import os

for survey in ['DECaLS']:
    for z in [0.001]:
        for age in [11.0]:
            print('Running python make_hpxcube.py -f ../data/%s_iso_hpxcube_z%.4f_a%.1f_gmax23_gmin20.2.fits.gz -mp 31 -a %.1f -z %.4f -s %s...' % (survey, z, age, age, z, survey))
            os.system('python make_hpxcube.py -f ../data/%s_iso_hpxcube_z%.4f_a%.1f_gmax23_gmin20.2.fits.gz -mp 31 -a %.1f -z %.4f -s %s' % (survey, z, age, age, z, survey))
