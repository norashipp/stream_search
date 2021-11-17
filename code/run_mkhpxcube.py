import os

counter = 0
# for survey in ['BASS_DR9']:
for survey in ['DECaLS_DR9', 'BASS_DR9', 'DELVE', 'DES_Y6_GOLD']:
    for z in [0.0001]:
        for age in [12.0]:
            print('Running python make_hpxcube.py -f ../data/%s_iso_hpxcube_z%.4f_a%.1f_gmax23_gmin16_b0.fits.gz -mp 10 -a %.1f -z %.4f -s %s...' % (survey, z, age, age, z, survey))
            os.system('csub -o log_%i python make_hpxcube.py -f ../data/2021/%s_iso_hpxcube_z%.4f_a%.1f_gmax23_gmin16_b0.fits.gz -mp 10 -a %.1f -z %.4f -s %s' % (counter, survey, z, age, age, z, survey))
            counter += 1
