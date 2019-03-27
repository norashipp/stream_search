#!/usr/bin/env python
"""
Generic python script.
"""
import os
from os.path import join, abspath,dirname,realpath
from collections import OrderedDict as odict
import numpy as np
import yaml

from results import create_results

### #print __file__
### #print realpath(__file__)
### DATADIR = join(dirname(realpath(__file__)),'../data')
### #STREAMFILE = join(DATADIR,'streams_v1.0.yaml')
### #STREAMFILE = join(DATADIR,'streams_v2.1.yaml')
### STREAMFILE = join(DATADIR,'streams_v3.0.yaml')
### STREAMS = yaml.load(open(STREAMFILE,'r'))
### RESULTS = create_results(STREAMS)
###  
### def load_streams(filename=STREAMFILE):
###     return yaml.load(open(STREAMFILE,'r'))

# Mostly from harris
# rh : half-light radius in arcmin
GLOBULARS = odict([
    ('NGC288',
     dict(ra=13.1885, dec=-26.5826, mods=[14.6,15.1], 
          mod = 14.85, rh=2.23, # arcmin
          jacobi = 76.4, #pc
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1261',
     dict(ra=48.0675, dec=-55.2162, mods=[15.8,16.3], #
          mod = 16.10, rh=0.68,
          jacobi = 146.4, #pc
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1851',
     dict(ra=78.52817, dec=-40.046, mods=[15.2,15.8], #
          mod = 15.5, rh=0.51,
          jacobi = 166.5, #pc
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1904',
     dict(ra=81.0462, dec=-24.52472, mods=[15.2,15.9], #
          mod = 15.56, rh=0.65,
          jacobi = 153.8, #pc
          kwargs=dict(color='chartreuse')
      )),
    ('NGC 7089',
     dict(ra=-36.637, dec=-0.823, mods=[15.3,15.7], #
          mod = 15.48, rh=1.06,
          kwargs=dict(color='chartreuse')
      )),
    ('Reticulum',
     dict(ra=69.0375, dec=-58.85833, mods=[18.1,18.7], #
          mod = 18.43, rh = np.nan,
          kwargs=dict(color='chartreuse')
      )),
    ('Whiting 1',
     dict(ra=30.7375, dec=-3.25277, mods=[17.1,17.7], #
          mod = 17.45, rh=0.22,
          kwargs=dict(color='chartreuse')
      )),
    ('AM-1',
     dict(ra=58.759583, dec=-49.615227, mods=[20.1,20.7], #
          mod = 20.41, rh=0.41, 
          kwargs=dict(color='chartreuse')
      )),
    ('Eridanus',
     dict(ra=66.18541667, dec=-21.18694, mods=[19.5,20.1], #
          mod = 19.78, rh=0.46,
          kwargs=dict(color='chartreuse')
      )),

    ])
GLOBULARS['M2'] = GLOBULARS['NGC 7089']

# From NED & McConnachie
DWARFS = odict([
    ('Fornax',
     dict(ra=39.99720, dec=-34.44919, mods=[-np.inf,np.inf], 
          mod = 20.84, rh = 16.60, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('Sculptor',
     dict(ra=15.03898, dec=-33.70903, mods=[-np.inf,np.inf], 
          mod = 19.67 , rh = 11.30, #arcmin
          kwargs=dict(color='magenta')
      )),
])

# From NED
GALAXIES = odict([
    ('Tucana',
     dict(ra=340.45667, dec=-64.41944, mods=[-np.inf,np.inf], 
          mod = 24.73, rh = 2.9, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('Phoenix',
     dict(ra=27.77642, dec=-44.44469, mods=[-np.inf,np.inf],
          # NED gave rh = 4.9, which wasn't large enough
          mod = 23.10, rh = 10.0, #arcmink
          kwargs=dict(color='magenta')
      )),
    ('IC 1613',
     dict(ra=16.19913, dec=2.11778, mods=[-np.inf,np.inf],
          mod = 24.29, rh = 16.2, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 55', # Very elliptical galaxy
     dict(ra=3.723333, dec=-39.196639, mods=[-np.inf,np.inf],
          mod = 26.32, rh = 19.0, # average arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 247',
     dict(ra=11.785625, dec=-20.760389, mods=[-np.inf,np.inf], 
          mod = 27.53, rh = 21.4, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 253',
     dict(ra=11.888002, dec=-25.288220, mods=[-np.inf,np.inf], 
          mod = 27.52, 
          #rh = 27.5, #arcmin from NED
          rh = 10.0, #just making this up
          kwargs=dict(color='magenta')
      )),
    ('NGC 300',
     dict(ra=13.722833, dec=-37.684389, mods=[-np.inf,np.inf], 
          mod = 26.43, rh = 21.9, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 1395',
     dict(ra=54.623955, dec=-23.027524, mods=[-np.inf,np.inf], 
          mod = 31.74, rh = 5.9, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 1399',
     dict(ra=54.620941, dec=-35.450657, mods=[-np.inf,np.inf], 
          mod = 31.23, rh = 6.9, #arcmin
          kwargs=dict(color='magenta')
      )),
    ('NGC 1407',
     dict(ra=55.049417, dec=-18.580111, mods=[-np.inf,np.inf], 
          mod = 31.81, rh = 4.6, #arcmin
          kwargs=dict(color='magenta')
      )),
])

""" 
ADW: This should all be deprecated
OTHERS = odict([
    ('Murrumbidgee',
     dict(name = 'Murrumbidgee',
          altname = None,
          ref = 'Grillmair:2017',
          ends=[(123.3,-65.), (358.7,16.3)], # ra, dec
          modulus=16.6,
          width_deg = 0.37,
          richness=np.nan,
          stellar_mass=np.nan,
          kwargs=dict(color='orange'),
          mods=[16, 18.4],
          quadrants = [2],
          cel_text = dict(lon=11,lat=-21,ha='center',va='bottom'),
          offset = 1.25,
          status = 0,
      )),
    ('Orinoco',
     dict(name = 'Orinoco',
          altname = None,
          ref = 'Grillmair:2017',
          ends=[(11.7, -21.5), (31.1, -33.2)], # ra, dec
          modulus=16.6,
          width_deg = 0.67,
          richness=np.nan,
          stellar_mass=np.nan,
          kwargs=dict(color='orange'),
          mods=[16, 18.4],
          quadrants = [2],
          cel_text = dict(lon=11,lat=-21,ha='center',va='bottom'),
          offset = 1.25,
          status = 0,
      )),
    ('Murrumbidgee',
     dict(name = 'Murrumbidgee',
          altname = None,
          ref = 'Grillmair:2017',
          ends=[(11.7, -21.5), (31.1, -33.2)], # ra, dec
          modulus=16.6,
          width_deg = 0.3,
          richness=np.nan,
          stellar_mass=np.nan,
          kwargs=dict(color='orange'),
          mods=[16, 18.4],
          quadrants = [2],
          cel_text = dict(lon=11,lat=-21,ha='center',va='bottom'),
          offset = 1.25,
          status = 0,
      )),
])

STREAMS2 = odict([
    ('ATLAS',
     dict(name = 'ATLAS',
          altname = None,
          ref = 'Koposov et al. 2014',
          #ends=[(11.7, -21.5), (31.3, -33.9)], # ra, dec
          ends=[(11.7, -21.5), (31.1, -33.2)], # ra, dec
          pole=(162.4, 3.7), # l,b
          modulus=16.7,
          distance=21.9, # kpc
          width=(133, 0.35), # pc, deg
          length=(7.5, 21.3), # kpc, deg
          richness=2.0e4,
          stellar_mass=4.7e3,
          kwargs=dict(color='orange'),
          mods=[16, 18.4],
          quadrants = [2],
          cel_text = dict(lon=11,lat=-21,ha='center',va='bottom'),
          offset = 1.25,
          status = 0,
      )),
    ('Molonglo',
     dict(name = 'Molonglo',
          altname = None,
          ref = 'Grillmair 2017',
          ends=[(6.4, -24.4), (20.9, -32.8)],
          pole=(336.7, 1.2),
          modulus=16.3,
          distance=18.2,
          width=(89, 0.28),
          length=(5, 15.2),
          stellar_mass=1.7e3,
          kwargs=dict(color='orange'),
          mods=[16.3, 16.9],
          quadrants = [2],
          cel_text = dict(lon=6,lat=-24,ha='left',va='bottom'),
          offset = 1.0,
          status = 0,
      )),
    ('Phoenix',
     dict(name = 'Phoenix',
          altname = None,
          ref = 'Balbinot et al. 2015',
          ends=[(20.0, -55.7), (28.3, -42.3)],
          #ends=[(20.8, -53.9), (28.0, -42.9)],
          pole=None,
          modulus=16.4,
          distance=19.0,
          width=(133, 0.40),
          length=(4.0, 12.0),
          stellar_mass=1.6e3,
          kwargs=dict(color='orange'),
          mods=[15.8, 17.7],
          quadrants = [3],
          cel_text = dict(lon=18,lat=-55,ha='left',va='bottom'),
          offset = 1.5,
          status = 0,
      )),
    ('Tucana III',
     dict(name = 'Tucana III',
          altname = 'Tucana',
          ref = 'Drlica-Wagner et al. 2015',
          ends=[(-6.3, -59.9), (3.5, -59.5)],
          pole=(285.2, 30.2),
          modulus=16.7,
          distance=21.9,
          width=(76, 0.2),
          length=(1.6, 4.1),
          stellar_mass=5.6e3,
          kwargs=dict(color='orange'),
          mods=[16.5, 18],
          quadrants = [1],
          cel_text = dict(ha='left',va='bottom'),
          offset = 1.0,
          status = 0,
      )),
    ('Caelum',
     dict(name = 'Caelum',
          altname = None,
          ref = None,
          ends=[(90.5, -45.6), (79.3, -34.3)],
          pole=(201.7, 50.7),
          modulus=15.7,
          distance=13.1,
          width=(107, 0.47),
          length=(3.3, 14.1),
          stellar_mass=1.5e3,
          kwargs=dict(color='cyan'),
          mods=[15.7, 16.7],
          quadrants = [1],
          offset = 1.0,
          status = 1,
      )),
    ('Eridanus',
     dict(name = 'Eridanus',
          altname = None,
          ref = None,
          ends=[(65.6, -21.2), (75.2, -26.4)],
          pole=(167.2, 35.3),
          modulus=16.9,
          distance=20.4,
          width=(155, 0.38),
          length=(4.3, 10.2),
          stellar_mass=2.6e3,
          kwargs=dict(color='cyan'),
          mods=[17.1, 18.0],
          quadrants = [1],
          offset = 1.0,
          status = 1,
      )),
    ('Fornax',
     dict(name = 'Fornax',
          altname = None,
          ref = None,
          ends=[(31.7, -31.5), (40.6, -38.3)],
          pole=(176.5, 9.7),
          modulus=17.3,
          distance=31.6,
          width=(254, 0.46),
          length=(5.6, 10.0),
          stellar_mass=3.1e3,
          kwargs=dict(color='cyan'),
          mods=[17.0, 18.0],
          quadrants = [2],
          offset = 1.0,
          status = 1,
      )),
    ('Indus',
     dict(name = 'Indus',
          altname = 'Ind 1',
          ref = None,
          ends=[(-37.1,-52.5),(-8.0, -64.8)],
          pole=(304.8, 39.1),
          modulus=15.8,
          distance=14.5,
          width=(291, 1.15),
          length=(4.8, 18.3),
          stellar_mass=4.7e3,
          kwargs=dict(color='cyan'),
          mods=[15.0, 16.8],
          quadrants = [1],
          cel_text = dict(ha='left',va='bottom'),
          offset = 1.5,
          status = 1,
      )),
    ('Jhelum (Indus2)',
     dict(name = 'Jhelum',
          altname = 'Ind 2',
          ref = None,
          ends=[(-38.8, -45.1), (4.73, -51.7)],
          pole=(285.8, 21.4),
          modulus=16.0,
          distance=15.8,
          width=(226, 0.82),
          length=(4.8, 18.3),
          stellar_mass=1.6e4,
          kwargs=dict(color='cyan'),
          mods=[15.0, 16.0],
          quadrants = [1],
          cel_text = dict(ha='left',va='bottom'),
          offset = 1.5,
          status = 1,
      )),
    ('Chenab (Indus3)',
     dict(name = 'Chenab',
          altname = 'Ind 3',
          ref = None,
          ends=[(-40.7, -59.9), (-28.3, -43.0)],
          pole=(34.0, 30.6),
          modulus=17.9,
          distance=38.0,
          width=(385, 0.58),
          length=(12.4, 18.5),
          stellar_mass=8.5e3,
          kwargs=dict(color='cyan'),
          mods=[17.3, 18.5],
          quadrants = [1],
          cel_text = dict(ha='left',va='bottom'),
          offset = 1.5,
          status = 1,
      )),
    ('Ravi (Indus4)',
     dict(name = 'Ravi',
          altname = 'Ind 4',
          ref = None,
          ends=[(-25.2, -44.1), (-16.0, -59.7)],
          pole=(353.5, 34.9),
          modulus=16.4,
          distance=19.1,
          width=(180, 0.54),
          length=(5.7, 16.6),
          stellar_mass=2.8e3,
          kwargs=dict(color='cyan'),
          mods=[16.5, 17.6],
          quadrants = [1],
          status = 1,
      )),
    ('Elqui (Smudge)',
     dict(name = 'Elqui',
          altname = 'Smudge',
          ref = None,
          ends=[(10.7, -36.9), (20.6, -42.4)],
          pole=(341.7, 8.7),
          modulus=18.5,
          distance=50.1,
          width=(324, 0.37),
          length=(10.1, 11.5),
          stellar_mass=5.9e3,
          kwargs=dict(color='cyan'),
          mods=[17.5, 19.6],
          quadrants = [2],
          offset = 1.5,
          status = 1,
      )),
    ('Reticulum',
     dict(name = 'Reticulum',
          altname = None,
          ref = None,
          ends=[(36.1, -64.6), (38.4, -58.3)],
          pole=(233.9, 26.6),
          modulus=17.3,
          distance=28.8,
          width=(166, 0.33),
          length=(3.2, 6.4),
          stellar_mass=2.4e3,
          kwargs=dict(color='cyan'),
          mods=[17.3, 18.2],
          quadrants = [3],
          offset = 1.0,
          status = 1,
      )),
    ('Loa',
     dict(name = 'Loa',
          altname = 'BigOne',
          ref = None,
          ends=[(28.,-60.), (27.9,-45.), (27.8,-30.),(25.7,-25.3), (16.7, 1.6)],
          pole=None,
          modulus=17.8,
          distance=None,
          width=(None, 2.0),
          length=(None, 60.0),
          stellar_mass=np.nan,
          kwargs=dict(color='cyan',dashes=(5,2)),
          mods=[17.2, 18.4],
          quadrants = [2,3],
          offset = 0,
          status = 1,
      )),
    # ADW: Derived characteristics are fabricated
    ('NGC 1851',
     dict(name = 'NGC 1851',
          altname = None,
          ref = None,
          ends=[(76.0,-35.3),(80.8, -45.0)],
          pole=None,
          modulus=15.5,
          distance=None,
          width=(None, 0.2),
          length=(None, 10),
          stellar_mass=np.nan,
          kwargs=dict(color='cyan'),
          mods=[15.1, 15.5],
          quadrants = [1],
          cel_text = dict(lon=81,lat=-45,ha='left',va='top'),
          offset = 1.5,
          status = 1,
      )),
    # ADW: Derived characteristics are fabricated
    ('Comet',
     dict(name = 'Comet',
          altname = None,
          ref = None,
          ends=[(56.1, -50.4),(77.4, -40.7)],
          pole=None,
          modulus=17.1,
          distance=None,
          width=(None, 0.2),
          length=(None, 40),
          stellar_mass=np.nan,
          kwargs=dict(color='cyan'),
          mods=[17.3, 17.9],
          quadrants = [1],
          cel_text = dict(ha='left',va='bottom'),
          offset = 1.5,
          status = 1,
      )),
    # ADW: Derived characteristics are fabricated
    ('Turbio',
     dict(name = 'Turbio',
          altname = None,
          ref = None,
          ends=[(28., -61.),(27.9, -46.0)],
          pole=None,
          modulus=16.8,
          distance=None,
          width=(None, 0.3),
          length=(None, 15),
          stellar_mass=np.nan,
          kwargs=dict(color='cyan'),
          mods=[16.1, 17.4],
          quadrants = [1],
          offset = 1.0,
          status = 1,
      )),

])

OLD_STREAMS = odict([
    #('Sagittarius',
    # dict(ends=[(25., -7.), (40., 2.)],
    #      kwargs=dict(color='orange'),
    #      mods=[15.6, 21],
    #      )),
    #('Magellanic',
    # dict(ends = [(80,-72),(355,-2.5)],
    #      kwargs = dict(color = 'orange'),
    #      mods = [15,21],
    #  )),
    ('ATLAS',
     dict(ends=[(11.7, -21.5), (31.3, -33.9)],
         #ends=[(12.0, -22.1), (29.9, -32.7)],
          kwargs=dict(color='orange'),
          mods=[16, 18.3],
          )),
    ('Phoenix',
     dict(ends=[(20.0, -55.7), (28.3, -42.3)],
          #ends=[(20.8, -53.9), (28.0, -42.9)],
          #ends=[(26.4, -45.8), (20.7, -54.2)],
          kwargs=dict(color='orange'),
          mods=[16.0, 17.2],
          )),
    ('Tucana III',
     dict(ends=[(-4.0, -59.8), (4.0, -59.5)],
          kwargs=dict(color='orange'),
          mods=[16.5, 18],
          )),
    ('Alpheus', # Grillmair
     dict(ends=[(21.6, -69.0), (27.7, -45.0)],
          #ends=[(19.55, -33.19), (6.19, -25.30)],
          #ends=[(20.859, -32.791), (6.357, -24.412)],
          kwargs=dict(color='orange'),
          mods=[14.0, 20.0],
          )),
    ('Molonglo',
     dict(ends=[(6.357, -24.412), (20.859, -32.791)],
          #ends=[(19.55, -33.19), (6.19, -25.30)],
          #ends=[(20.859, -32.791), (6.357, -24.412)],
          kwargs=dict(color='orange'),
          mods=[16.3, 16.9],
          )),
    ('Indus',
     dict(ends=[(-39.3,-53.3),(-10.253, -64.107)],
          #ends=[(-10.253, -64.107), (-37.947, -51.887)],
          #ends=[(-28.7, -57.1),(-14.8, -62.5)],
          kwargs=dict(color='cyan'),
          mods=[15.0, 16.8],
          )),
    ('Ravi (Indus2)',
     dict(ends=[(-8.7, -51.3), (-37.5, -47.3)],
          kwargs=dict(color='cyan'),
          mods=[15.0, 16.0],
          )),
    ('Chenab (Indus3)',
     dict(ends=[(-40.684, -59.939), (-28.3, -43.0)],
          #ends=[(-38.4, -60.3),(-27.3, -42.0)]
          kwargs=dict(color='cyan'),
          mods=[17.3, 18.5],
          )),
    ('Indus4',
     dict(ends=[(-25.2, -44.1), (-16.029, -59.665)],
          kwargs=dict(color='cyan'),
          mods=[16.5, 17.6],
          )),
    ('Elqui (Smudge)',
     dict(ends=[(10.7, -36.9), (22.9, -43.6)],
          kwargs=dict(color='cyan'),
          mods=[17.5, 19.6],
          )),
    #('Smudge',
    # dict(ends=[(10.6, -37.3), (24.1, -43.8)],
    #      kwargs=dict(color='cyan'),
    #      mods=[18.6, 19.2],
    #      )),
    #('Smudge2',
    # dict(ends=[(14.44, -38.80), (22.82, -44.39)],
    #      kwargs=dict(color='cyan'),
    #      mods=[17.1, 18.3],
    #      )),
    ('Eridanus',
     dict(ends=[(65.60, -21.19), (75.19, -26.44)],
          kwargs=dict(color='cyan'),
          mods=[16.9, 18.0],
          )),
    ('Fornax',
     dict(ends=[(31.7, -31.5), (40.6, -38.3)],
          #ends=[(27.907, -28.488), (43.593, -39.096)],
          #ends=[(27.907, -28.488), (43.593, -39.096)],
          #ends = [(39.5, -36.6), (27.9, -28.9)],
          kwargs=dict(color='cyan'),
          mods=[17.0, 17.7],
          )),
    ('Reticulum',
     dict(ends=[(36.06, -64.62), (38.37, -58.30)],
          #ends=[(35.9, -64.7), (38.7, -58.3)],
          kwargs=dict(color='cyan'),
          mods=[17.3, 18.2],
          )),
    # ADW: Don't really see this one
    # ('Reticulum2',
    #  dict(ends=[(46.056, -64.471), (47.006, -56.298)],
    #       kwargs=dict(color='cyan'),
    #       mods=[15.6, 16.2],
    #       )),
    # ADW: Don't really see this one
    #('Horologium',
    # dict(ends=[(63.01, -56.30), (44.44, -44.33)],
    #      #ends=[(44.103, -43.310), (65.348, -56.563)],
    #      kwargs=dict(color='cyan'),
    #      mods=[16.2, 17.1],
    #      )),
    # ADW: Don't really see this one
    #('Cetus',
    # dict(ends=[(25.85, -25.86), (21.39, -20.88)],
    #      #ends=[(26.403, -26.527), (20.939, -16.933)],
    #      kwargs=dict(color='cyan'),
    #      mods=[17.3, 18.0],
    #      )),
    ('Caelum',
     dict(ends=[(90.51, -45.55), (79.25, -34.34)],
          #ends=[(87.94, -43.73), (77.11, -29.22)],
          kwargs=dict(color='cyan'),
          mods=[15.3, 16.5],
          )),
    #('Pictor',
    # dict(ends=[(78.9, -41.2), (80.9, -52.2)],
    #      #ends=[(86.19, -54.98), (80.49, -45.55)],
    #      #ends=[(79.87, -43.13), (82.25, -51.13)],
    #      kwargs=dict(color='cyan'),
    #      mods=[14.6, 15.5],
    #      )),

    ('Loa (BigOne)',
     dict(ends=[(28., -65.), (27,-45), (26,-30),(25,-20),(15, 3.0)],
          #ends=[(86.19, -54.98), (80.49, -45.55)],
          #ends=[(79.87, -43.13), (82.25, -51.13)],
          kwargs=dict(color='cyan'),
          mods=[17.0, 19.0],
      )),

    #('NGC1851',
    # dict(ends=[(76.1, -35.3), (80.8, -44.8)],
    #      kwargs=dict(color='cyan'),
    #      mods=[14.9, 15.2],
    #      )),
    ('NGC 1851',
     dict(ends=[(76.01,-35.29),(80.77, -44.97)],
          kwargs=dict(color='cyan'),
          mods=[15.1, 15.5],
      )),

    ('Comet',
     dict(ends=[(56.1, -50.4),(77.4, -40.7)],
          kwargs=dict(color='cyan'),
          mods=[17.3, 17.9],
      )),

    #
    #('NGC288',
    # dict(ends=[(12.94, -31.67), (12.99, -31.5)],
    #      kwargs=dict(color='cyan'),
    #      mods=[14.9, 14.9],
    #      )),

    #('RA=30',
    # dict(ends=[(30., -70.), (30., 7.)],
    #      #ends=[(86.19, -54.98), (80.49, -45.55)],
    #      #ends=[(79.87, -43.13), (82.25, -51.13)],
    #      kwargs=dict(color='red',dashes=(2,2)),
    #      mods=[18.0, 20.0],
    #      )),
    #('RA=60',
    # dict(ends=[(60.,-70.), (60., 7.)],
    #      #ends=[(86.19, -54.98), (80.49, -45.55)],
    #      #ends=[(79.87, -43.13), (82.25, -51.13)],
    #      kwargs=dict(color='red',dashes=(2,2)),
    #      mods=[18.0, 20.0],
    #      )),

])


GLOBULARS = odict([
    ('NGC288',
     dict(ra=13.1885, dec=-26.5826, mods=[14.6,15.1], #
          mod = 14.85, rh=2.23,
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1261',
     dict(ra=48.0675, dec=-55.2162, mods=[15.8,16.3], #
          mod = 16.10, rh=0.68,
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1851',
     dict(ra=78.52817, dec=-40.046, mods=[15.2,15.8], #
          mod = 15.5, rh=0.51,
          kwargs=dict(color='chartreuse')
      )),
    ('NGC1904',
     dict(ra=81.0462, dec=-24.52472, mods=[15.2,15.9], #
          mod = 15.56, rh=0.65,
          kwargs=dict(color='chartreuse')
      )),
    ('NGC 7089',
     dict(ra=-36.637, dec=-0.823, mods=[15.3,15.7], #
          mod = 15.48, rh=1.06
          kwargs=dict(color='chartreuse')
      )),
    ('Reticulum',
     dict(ra=69.0375, dec=-58.85833, mods=[18.1,18.7], #
          mod = 18.43, 
          kwargs=dict(color='chartreuse')
      )),
    ('Whiting 1',
     dict(ra=30.7375, dec=-3.25277, mods=[17.1,17.7], #
          mod = 17.45, rh=0.22,
          kwargs=dict(color='chartreuse')
      )),
    ('AM-1',
     dict(ra=58.759583, dec=-49.615227, mods=[20.1,20.7], #
          mod = 20.41, rh=0.41, 
          kwargs=dict(color='chartreuse')
      )),
    ('Eridanus',
     dict(ra=66.18541667, dec=-21.18694, mods=[19.5,20.1], #
          mod = 19.78, rh=0.46,
          kwargs=dict(color='chartreuse')
      )),

    ])
GLOBULARS['M2'] = GLOBULARS['NGC 7089']

DWARFS = odict([
    ('Fornax',
     dict(ra=39.99720, dec=-34.44919, mods=[-np.inf,np.inf], 
          mod = 20.84, rh = 710,
          kwargs=dict(color='magenta')
      )),
    ('Sculptor',
     dict(ra=15.03898, dec=-33.70903, mods=[-np.inf,np.inf], 
          mod = 19.67 , rh = 283,
          kwargs=dict(color='magenta')
      )),
])

SELECTION_PARAMS = odict([
    ('ATLAS', ('y', 0.5 / 2., 1.0)),  # 'y', 1.5, 3
    ('Phoenix', ('x', 0.5 / 2., 1.0)),  # ('x', 1.5, 10.)
    ('TucIII', ('y', 0.75 / 2., 1.5)),  # ('y', 1., 4.)
    ('Vertical', ('x', 1.5, 5.)),
    ('New', ('y', 0.8, 4.)),
    ('Sagittarius', ('y', 5., 10.)),
    ('Smudge', ('y', 1.0 / 2., 2.0)),
    ('Eridanus', ('x', 2.0 / 2., 2.0)),
    ('Phoenix3', ('y', 0.75 / 2., 1.0)),
    ('Indus', ('x', 3. / 2., 3.5)),  # 1: left
    ('Indus2', ('x', 3. / 2., 3.5)),
    ('Indus3', ('x', 2.5 / 2., 3.0)),  # 2: right
    ('Indus4', ('x', 1.0 / 2., 2.0)),
    ('Fornax', ('x', 1.0 / 2, 10.)), # 1: left
    ('Reticulum', ('y', 0.5 / 2, 5.0)),
    ('Reticulum2', ('x', 0.5 / 2, 1.0)),
    ('Horologium', ('y', 1.0 / 2, 1.0)),
    ('Cetus', ('x', 1.0 / 2, 2.0)), # 0.5 / 2, 1.0
    ('Caelum', ('y', 1.5 / 2., 2.0)),
    ('Pictor', ('x', 1.0 / 2, 5.0))
])


# r, mu, age, Z (age, Z not fit)
ISO_PARAMS = odict([
    ('ATLAS', (2.0e4, 16.7, 13., 0.0004)), # ??
    ('Phoenix', (6.8e3, 16.4, 11.5, 0.0004)),
    ('Tucana', (5.5e4, 16.5, 13., 0.0001)), # !! 16.8 looks way better
    ('Sagittarius', (2.1e6, 17.3, 13., 0.0004)),
    ('Smudge', (2.5e4, 18.5, 13., 0.0004)),
    ('Eridanus', (1.1e4, 16.9, 13., 0.0004)),
    ('Phoenix3', (7.1e3, 16.3, 13., 0.0004)),
    ('Indus', (2e4, 15.8, 13., 0.0004)),
    ('Indus2', (1.8e4, 15.7, 13., 0.0004)),
    ('Indus3', (3.6e4, 17.9, 13., 0.0004)),
    ('Indus4', (1.0e4, 15.7, 13., 0.0004)), # !!
    ('Fornax', (1.3e4, 17.3, 13., 0.0004)),
    ('Reticulum', (1.0e4, 17.3, 13., 0.0004)), # ??
    ('Horologium', (7.1e3, 15.8, 13., 0.0004)), # !!
    ('Caelum', (6.4e3, 15.7, 13., 0.0004)), # ??
    ('Pictor', (3.8e3, 15.5, 13., 0.0004)), # ??
    # ('Reticulum2', ()),
    # ('Cetus', ())
])



def select_stream(p1, p2, ra, dec, w=0.8, axis='x'):
    m = (p1[1] - p2[1]) / (p1[0] - p2[0])
    b = p1[1] - m * p1[0]

    #print 'line parameters = %.2f, %.2f' % (m, b)

    if axis == 'x':
        x = (dec - b) / m

        dec_min = np.minimum(p1[1], p2[1])
        dec_max = np.maximum(p1[1], p2[1])

        idx = (np.abs(ra - x) < w) & (dec > dec_min) & (dec < dec_max)

    elif axis == 'y':
        y = m * ra + b

        ra_min = np.minimum(p1[0], p2[0])
        ra_max = np.maximum(p1[0], p2[0])

        idx = (np.abs(dec - y) < w) & (ra > ra_min) & (ra < ra_max)

    else:
        raise IOError('Valid axes are x, y.')

    return idx


def select_off_stream(p1, p2, ra, dec, d=1.5, w=0.8, dec_min=-57, dec_max=-45, axis='x'):
    m = (p1[1] - p2[1]) / (p1[0] - p2[0])
    b = p1[1] - m * p1[0]

    if axis == 'x':
        x = (dec - b) / m

        dec_min = np.minimum(p1[1], p2[1])
        dec_max = np.maximum(p1[1], p2[1])

        idx1 = (np.abs(ra - (x + d)) < w) & (dec > dec_min) & (dec < dec_max)
        idx2 = (np.abs(ra - (x - d)) < w) & (dec > dec_min) & (dec < dec_max)

    elif axis == 'y':
        y = m * ra + b

        ra_min = np.minimum(p1[0], p2[0])
        ra_max = np.maximum(p1[0], p2[0])

        idx1 = (np.abs(dec - (y + d)) < w) & (ra > ra_min) & (ra < ra_max)
        idx2 = (np.abs(dec - (y - d)) < w) & (ra > ra_min) & (ra < ra_max)

    else:
        raise IOError('Valid axes are x, y.')

    return idx1, idx2

"""
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()



