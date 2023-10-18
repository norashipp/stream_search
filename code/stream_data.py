stream_dists = {'Tucana III': 25.1,
                'ATLAS': 22.9,
                'Phoenix': 19.1,
                'Indus': 16.6,
                'Jhelum': 13.2,
                'Chenab': 39.8,
                'Elqui': 50.1,
                'Aliqa Uma': 28.8,
                'Turranburra': 27.5,
                '300S':18,
                'Ravi':22.9,
                'Herc':132}

stream_rvs = {'ATLAS': -99.8,
              'Phoenix': 47.9,
              'Indus': -47.2,
              'Jhelum': (0.20, -34.0),
              'Chenab': -145.5,
              'Elqui': -33.9,
              'Aliqa Uma': -34.0,
              'Ravi':0.,
              'Herc':46.81}


stream_masses = {'Tucana III': 2e-7,
                 'ATLAS': 2e-6,
                 'Phoenix': 2e-6,
                 'Indus': 0.001,
                 'Jhelum': 0.002,
                 'Chenab': 0.001,
                 'Elqui': 0.0001,
                 'Aliqa Uma': 2e-6,
                 'Turranburra': 0.0001,
                 '300S':2e-6,
                 'Ravi':0.0001,
                 'Herc':2e-6}

stream_sigmas = {'Tucana III': 0.05,
                 'ATLAS': 0.01,
                 'Phoenix': 0.01,
                 'Indus': 0.1,
                 'Jhelum': 0.1,
                 'Chenab': 0.5,
                 'Elqui': 0.1,
                 'Aliqa Uma': 0.01,
                 'Turranburra': 0.1,
                 '300S':0.01,
                 'Ravi':0.01,
                 'Herc':0.01}


stream_matrices = {'Aliqa Uma': [[0.66315359, 0.48119409, -0.57330582], [0.74585903, -0.36075668, 0.5599544], [-0.06262284, 0.79894109, 0.59814004]],
                   'ATLAS': [[0.83697865, 0.29481904, -0.4610298], [0.51616778, -0.70514011, 0.4861566], [0.18176238, 0.64487142, 0.74236331]],
                   'Chenab': [[0.51883185, -0.34132444, -0.78378003], [-0.81981696, 0.06121342, -0.56934442], [-0.24230902, -0.93795018, 0.2480641]],
                   'Elqui': [[0.74099526, 0.20483425, -0.63950681], [0.57756858, -0.68021616, 0.45135409], [0.34255009, 0.70381028, 0.62234278]],
                   'Indus': [[0.47348784, -0.22057954, -0.85273321], [0.25151201, -0.89396596, 0.37089969], [0.84412734, 0.39008914, 0.3678036]],
                   'Jhelum': [[0.60334991, -0.20211605, -0.7714389], [-0.13408072, -0.97928924, 0.15170675], [0.78612419, -0.01190283, 0.61795395]],
                   'Phoenix': [[0.5964467, 0.27151332, -0.75533559], [-0.48595429, -0.62682316, -0.60904938], [0.63882686, -0.73032406, 0.24192354]],
                   'Tucana III': [[0.505715, -0.007435, -0.862668], [-0.078639, -0.996197, -0.037514], [0.859109, -0.086811, 0.504377]],
                   'Turranburra': [[0.36111266, 0.85114984, -0.38097455], [0.87227667, -0.16384562, 0.46074725], [-0.32974393, 0.49869687, 0.80160487]],
                   'Jet': [[-0.69798645,  0.61127501, -0.37303856], [-0.62615889, -0.26819784,  0.73211677], [0.34747655,  0.744589,  0.56995374]],
                   '300S': [[-0.88197819,  0.38428506,  0.27283596], [ 0.43304104,  0.88924457,  0.14737555], [-0.18598367,  0.24813119, -0.95070552]],
                   'Ravi': [[0.57336113, -0.22475898, -0.78787081], [0.57203155, -0.57862539, 0.58135407], [0.58654661, 0.78401279, 0.20319208]],
                   'Herc': [[-0.36919252, -0.90261961,  0.22130235], [-0.75419217,  0.43013129,  0.49616655], [-0.54303872,  0.01627647, -0.8395499 ]]}

stream_lengths = {'Aliqa Uma': 10.0,
                  'ATLAS': 22.6,
                  'Chenab': 18.5,
                  'Elqui': 9.4,
                  'Indus': 20.3,
                  'Jhelum': 29.2,
                  'Phoenix': 13.6,
                  'Tucana III': 4.8,
                  'Turranburra': 16.9,
                  'Jet':20.0,
                  '300S':17.0,
                  'Ravi':60.,
                  'Herc':1.}

stream_widths = {'Aliqa Uma': 0.26,
                 'ATLAS': 0.24,
                 'Chenab': 0.71,
                 'Elqui': 0.54,
                 'Indus': 0.83,
                 'Jhelum': 1.16,
                 'Phoenix': 0.16,
                 'Tucana III': 0.18,
                 'Turranburra': 0.60,
                 'Jet':0.2,
                 '300S':0.47,
                 'Ravi':0.72,
                 'Herc':6.3/60.}

stream_vrs = {'Aliqa Uma': 0,
              'ATLAS': 0,
              'Chenab': -150,
              'Elqui': 0,
              'Indus': 0,
              'Jhelum': 0,
              'Phoenix': 0,
              'Tucana III': 0,
              'Turranburra': 0,
              'Jet':272.5,
              '300S':300,
              'Ravi': 0,
              'Herc':46.81}

stream_vr_widths = {'Aliqa Uma': 0,
                    'ATLAS': 0,
                    'Chenab': 4.32,
                    'Elqui': 8.40,
                    'Indus': 5.76,
                    'Jhelum': 13.30,
                    'Phoenix': 0,
                    'Tucana III': 0,
                    'Turranburra': 0,
                    'Jet':0,
                    '300S':0,
                    'Ravi':0,
                    'Herc':0}

stream_phi2s = {'Aliqa Uma': 0,
                'ATLAS': 0.66,
                'Chenab': 0,
                'Elqui': 0,
                'Indus': 0,
                'Jhelum': 0,
                'Phoenix': 0,
                'Tucana III': 0,
                'Turranburra': 0,
                'Jet':0,
                '300S':0,
                'Ravi':0,
                'Herc':0}

stream_mids = {'Aliqa Uma': (35.96519575304428, -34.98107408877647),
               'ATLAS': (19.40440473586367, -27.453578354852883),
               'Chenab': (-33.339759611694205, -51.60798647784715),
               'Elqui': (15.452462557149355, -39.75505320589948),
               'Indus': (-24.978978038471325, -58.510205017947726),
               'Jhelum': (-18.520322938617255, -50.483277585344595),
               'Phoenix': (24.475886399510763, -49.054701430221975),
               'Tucana III': (-1.4880868099542681, -59.641037041735935),
               'Turranburra': (67.01021150038846, -22.394061648029155),}

stream_phi12_pms = {'Aliqa Uma': {'pm1': 0.9803, 'e_pm1': 0.0350, 'pm2': -0.3416, 'e_pm2': 0.0277, 'grad_pm1': -0.0224, 'e_grad_pm1': 0.0209, 'grad_pm2': -0.0362, 'e_grad_pm2': 0.0212},
                    'ATLAS': {'pm1': 1.6602, 'e_pm1': 0.0428, 'pm2': -0.1537, 'e_pm2': 0.0351, 'grad_pm1': 0.0155, 'e_grad_pm1': 0.0049, 'grad_pm2': -0.0179, 'e_grad_pm2': 0.0044},
                    'Chenab': {'pm1': 1.0336, 'e_pm1': 0.0454, 'pm2': -0.5975, 'e_pm2': 0.0287, 'grad_pm1': 0.0440, 'e_grad_pm1': 0.0130, 'grad_pm2': -0.0213, 'e_grad_pm2': 0.0084},
                    'Elqui': {'pm1': 0.5584, 'e_pm1': 0.0606, 'pm2': -0.0280, 'e_pm2': 0.0491, 'grad_pm1': -0.0270, 'e_grad_pm1': 0.0199, 'grad_pm2': -0.0433, 'e_grad_pm2': 0.0141},
                    'Indus': {'pm1': -3.0886, 'e_pm1': 0.0319, 'pm2': 0.2053, 'e_pm2': 0.0285, 'grad_pm1': 0.0542, 'e_grad_pm1': 0.0043, 'grad_pm2': 0.0436, 'e_grad_pm2': 0.0041},
                    # 'Jhelum': {'pm1': -5.9330, 'e_pm1': 0.9163, 'pm2': -0.7612, 'e_pm2': 0.7987, 'grad_pm1': 0.0258, 'e_grad_pm1': 0.1969, 'grad_pm2': 0.0347, 'e_grad_pm2': 0.1505},
                    'Jhelum': {'pm1': -5.9330, 'e_pm1': 0.03, 'pm2': -0.7612, 'e_pm2': 0.05, 'grad_pm1': 0.0, 'e_grad_pm1': 0.0, 'grad_pm2': 0.0, 'e_grad_pm2': 0.0},
                    'Phoenix': {'pm1': -1.9439, 'e_pm1': 0.0216, 'pm2': -0.3649, 'e_pm2': 0.0227, 'grad_pm1': -0.0091, 'e_grad_pm1': 0.0062, 'grad_pm2': 0.0088, 'e_grad_pm2': 0.0068},
                    'Tucana III': {'pm1': 1.0835, 'e_pm1': 0.0311, 'pm2': -0.0260, 'e_pm2': 0.0343, 'grad_pm1': 0.1200, 'e_grad_pm1': 0.0309, 'grad_pm2': -0.0618, 'e_grad_pm2': 0.0319},
                    'Turranburra': {'pm1': 0.6922, 'e_pm1': 0.0455, 'pm2': -0.2223, 'e_pm2': 0.0436, 'grad_pm1': 0.0016, 'e_grad_pm1': 0.0159, 'grad_pm2': -0.0287, 'e_grad_pm2': 0.0138}}

stream_radec0_pms = {'Aliqa Uma': {'pmra': 0.2465, 'e_pmra': 0.0330, 'pmdec': -0.7073, 'e_pmdec': 0.0517},
                     'ATLAS': {'pmra': 0.0926, 'e_pmra': 0.0326, 'pmdec': -0.8783, 'e_pmdec': 0.0328},
                     'Chenab': {'pmra': 0.3223, 'e_pmra': 0.0365, 'pmdec': -2.4659, 'e_pmdec': 0.0434},
                     'Elqui': {'pmra': 0.1311, 'e_pmra': 0.0387, 'pmdec': -0.3278, 'e_pmdec': 0.0923},
                     'Phoenix': {'pmra': 2.7572, 'e_pmra': 0.0217, 'pmdec': -0.0521, 'e_pmdec': 0.0222},
                     'Tucana III': {'pmra': -0.0995, 'e_pmra': 0.0390, 'pmdec': -1.6377, 'e_pmdec': 0.0373},
                     'Turranburra': {'pmra': 0.4348, 'e_pmra': 0.0386, 'pmdec': -0.8875, 'e_pmdec': 0.0426}}


stream_phi120_pms = {'Aliqa Uma': {'pm1': -0.6634, 'pm2': -0.3479},
                     'ATLAS': {'pm1': -0.5586, 'pm2': -0.6841},
                     'Chenab': {'pm1': 2.1318, 'pm2': -1.2805},
                     'Elqui': {'pm1': -0.2986, 'pm2': -0.1883},
                     'Phoenix': {'pm1': -0.9694, 'pm2': -2.5817},
                     'Indus': {'pm1': -6.33, 'pm2': -1.34},
                     'Jhelum': {'pm1': -8.04, 'pm2': -3.98},
                     'Tucana III': {'pm1': 0.2048, 'pm2': -1.6279},
                     'Turranburra': {'pm1': -0.8193, 'pm2': -0.5528}}

stream_peri_apo = {'Aliqa Uma': (15.82, 46.04),
                   'ATLAS': (14.62, 47.64),
                   'Chenab': (33.09, 81.53),
                   'Elqui': (7.25, 66.37),
                   'Phoenix': (11.40, 17.33),
                   'Indus': (11.63, 19.25),
                   'Jhelum': (8.83, 34.17),
                   'Jet': (8.94, 38.98)}
