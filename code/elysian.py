#!/usr/bin/env python
"""
Various Local Group coordinate frames. 

Adapted from:
http://docs.astropy.org/en/stable/generated/examples/coordinates/plot_sgr-coordinate-frame.html
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict as odict

from astropy.coordinates import frame_transform_graph
#from astropy.coordinates.angles import rotation_matrix
from astropy.coordinates.matrix_utilities import rotation_matrix
import astropy.coordinates as coord
import astropy.units as u
import healpy as hp

#from ugali.utils.projector import gal2cel,cel2gal
from results import create_result
# These should be in this module, but want `results` to be standalone
from results import find_pole, euler_angles, create_matrix

__all__ = ['Stream']

#def find_pole(lon1,lat1,lon2,lat2):
#    vec = np.cross(hp.ang2vec(lon1,lat1,lonlat=True),
#                   hp.ang2vec(lon2,lat2,lonlat=True))
#    lon,lat = hp.vec2ang(vec,lonlat=True)
#    return [np.asscalar(lon),np.asscalar(lat)]

class Stream(coord.BaseCoordinateFrame):
    """
    A Heliocentric spherical coordinate system defined by the orbit
    of the Sagittarius dwarf galaxy, as described in
        http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
    and further explained in
        http://www.astro.virginia.edu/~srm4n/Sgr/.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    Lambda : `Angle`, optional, must be keyword
        The longitude-like angle corresponding to Sagittarius' orbit.
    Beta : `Angle`, optional, must be keyword
        The latitude-like angle corresponding to Sagittarius' orbit.
    distance : `Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    """
    default_representation = coord.SphericalRepresentation

    frame_specific_representation_info = {
        'spherical': [coord.RepresentationMapping('lon', 'Lambda'),
                      coord.RepresentationMapping('lat', 'Beta'),
                      coord.RepresentationMapping('distance', 'distance')],
        'unitspherical': [coord.RepresentationMapping('lon', 'Lambda'),
                          coord.RepresentationMapping('lat', 'Beta')]
    }

    @property 
    def L(self): return self.Lambda

    @property 
    def B(self): return self.Beta
        
    def get_arc(self,lon0,lat0,lon1,lat1,offset=0):
        if not hasattr(self,'ends'):
            return AttributeError()
        # The ends are stored in ICRS
        ends = coord.SkyCoord([lon0,lon1],[lat0,lat1],
                              unit=(u.deg,u.deg),
                              frame='icrs').transform_to(self)
        ends.Lambda.wrap_angle = 180*u.deg
        #print(ends.Lambda)
        L = np.linspace(ends.Lambda.min(), ends.Lambda.max(),100)
        B = np.zeros_like(L) + offset*u.deg
        # Return in ICRS
        icrs = coord.SkyCoord(L,B,frame=self).transform_to(coord.ICRS)
        return icrs.ra.value,icrs.dec.value

def sky_to_stream(lon, lat, matrix):
    """ Compute the transformation from spherical sky coordinates
        (ICRS, Galactic, etc. ) to heliocentric stream coordinates.
    Parameters:
    -----------
    lon : sky longitude (deg)
    lat : sky latitude (deg)
    matrix : rotation matrix

    Returns:
    --------
    L,B : stream coordinates (deg)
    """
    lon = np.atleast_1d(np.radians(lon))
    lat = np.atleast_1d(np.radians(lat))

    X = np.cos(lat)*np.cos(lon)
    Y = np.cos(lat)*np.sin(lon)
    Z = np.sin(lat)

    # Calculate X,Y,Z,distance in the stream system
    Xs, Ys, Zs = matrix.dot(np.array([X, Y, Z]))
    Zs = -Zs

    # Calculate the angular coordinates lambda,beta
    L = np.arctan2(Ys,Xs)
    B = np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs))
    return np.degrees(L),np.degrees(B)

def stream_to_sky(L, B, matrix):
    """ 
    Compute the transformation from heliocentric stream coordinates to
    spherical sky coordinates.

    Parameters:
    -----------
    L : stream longitude (deg)
    B : stream latitude (deg)
    matrix : rotation matrix

    Returns:
    --------
    lon,lat : sky coordinates (deg)
    """
    L = np.atleast_1d(np.radians(L))
    B = np.atleast_1d(np.radians(B))

    Xs = np.cos(B)*np.cos(L)
    Ys = np.cos(B)*np.sin(L)
    Zs = np.sin(B)
    Zs = -Zs

    X, Y, Z = matrix.T.dot(np.array([Xs, Ys, Zs]))

    lon = np.arctan2(Y,X)
    lat = np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z))
    #lon[lon<=0] += 2*np.pi*u.radian
    return np.degrees(lon),np.degrees(lat)

def euler_factory(name, phi, theta, psi, ends=None, pole=None, flip=False):
    cls = type(name, (Stream,),{'__init__':Stream.__init__})
    cls.ends = ends if ends is None else np.array(ends)
    cls.pole = pole if pole is None else np.array(pole)
    cls.euler_angles = phi,theta,psi
    cls.flip = flip

    MATRIX = create_matrix(phi,theta,psi)

    @frame_transform_graph.transform(coord.FunctionTransform, coord.ICRS, cls)
    def icrs_to_stream(sky_coord, stream_frame, cls=cls, matrix=MATRIX):
        """ Compute the transformation from ICRS spherical to
            heliocentric stream coordinates.
        """
        Lambda,Beta=sky_to_stream(sky_coord.ra.deg,sky_coord.dec.deg,matrix)
        Lambda[Lambda < 0] +=  360.
        Lambda *= (-1)**cls.flip
        return cls(Lambda=Lambda*u.deg, Beta=Beta*u.deg,
                   distance=sky_coord.distance)

    @frame_transform_graph.transform(coord.FunctionTransform, cls, coord.ICRS)
    def stream_to_icrs(stream_coord, sky_frame, cls=cls, matrix=MATRIX):
        """ Compute the transformation from heliocentric stream coordinates to
            spherical ICRS.
        """
        lon,lat=stream_to_sky((-1)**cls.flip * stream_coord.Lambda.deg,
                              stream_coord.Beta.deg,matrix)
        lon[lon<=0] += 360.
        return coord.ICRS(ra=lon*u.deg, dec=lat*u.deg,
                          distance=stream_coord.distance)
    return cls

def stream_factory(stream):
    result = create_result(stream)
    name = result['name'].replace(' ','')

    ends = np.array(result['ends'])
    pole = np.array(result['pole'])
    phi,theta,psi = np.array(result['euler_angles'])
    return euler_factory(name, phi, theta, psi, ends=ends, pole=pole)

factory = stream_factory

# There is some craziness here because Nidever et al. 2008 performs a mirror flip
# "As in the Magellanic coordinate system of Wakker, L_MS decreases along the MS (toward its tail)."
Magellanic  = euler_factory("Magellanic", 84.295-90, 90+17.7848, 180-268.9597,
                            flip=True)

Sagittarius = euler_factory("Sagittarius", 180+3.75, 90-13.46, 180+14.111534)

#def factory(stream):
#    """
#    Create a stream coordinate system.
#
#    Parameters:
#    -----------
#    stream : the stream dictionary including name and ends
#
#    Returns:
#    --------
#    cls    : Sream coordinate class
#    """
#
#    result = create_result(stream)
#    name = result['name'].replace(' ','')
#
#    # Create the new class dynamically from the Stream base class
#    cls = type(name, (Stream,),{'__init__':Stream.__init__})
#    cls.ends = np.array(result['ends'])
#    cls.pole = np.array(result['pole'])
#    cls.euler_angles = phi,theta,psi = np.array(result['euler_angles'])
#
#    MATRIX = create_matrix(phi,theta,psi)
#
#    # Generate the rotation matrix using the x-convention (see Goldstein)
#    #D = rotation_matrix(np.radians(phi),   "z", unit=u.radian)
#    #C = rotation_matrix(np.radians(theta), "x", unit=u.radian)
#    #B = rotation_matrix(np.radians(psi),   "z", unit=u.radian)
#    #MATRIX = np.array(B.dot(C).dot(D))
#
#    @frame_transform_graph.transform(coord.FunctionTransform, coord.ICRS, cls)
#    def icrs_to_stream(sky_coord, stream_frame, cls=cls, matrix=MATRIX):
#        """ Compute the transformation from ICRS spherical to
#            heliocentric stream coordinates.
#        """
#        Lambda, Beta = sky_to_stream(sky_coord.ra.deg,sky_coord.dec.deg,matrix)
#        Lambda[Lambda < 0] +=  360.
#        return cls(Lambda=Lambda*u.deg, Beta=Beta*u.deg,
#                   distance=sky_coord.distance)
#
#        #lon = np.atleast_1d(sph_coord.ra.radian)
#        #lat = np.atleast_1d(sph_coord.dec.radian)
#        #
#        #X = np.cos(lat)*np.cos(lon)
#        #Y = np.cos(lat)*np.sin(lon)
#        #Z = np.sin(lat)
#        #
#        ## Calculate X,Y,Z,distance in the stream system
#        #Xs, Ys, Zs = matrix.dot(np.array([X, Y, Z]))
#        #Zs = -Zs
#        #
#        ## Calculate the angular coordinates lambda,beta
#        #Lambda = np.arctan2(Ys,Xs)*u.radian
#        #Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2*np.pi*u.radian
#        #Beta = np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs))*u.radian
#        #return cls(Lambda=Lambda, Beta=Beta,distance=sph_coord.distance,
#
#    @frame_transform_graph.transform(coord.FunctionTransform, cls, coord.ICRS)
#    def stream_to_icrs(stream_coord, sky_frame, cls=cls, matrix=MATRIX):
#        """ Compute the transformation from heliocentric stream coordinates to
#            spherical ICRS.
#        """
#        lon,lat = stream_to_sky(stream_coord.Lambda.deg,stream_coord.Beta.deg,matrix)
#        lon[lon<=0] += 360.
#        return coord.ICRS(ra=lon*u.deg, dec=lat*u.deg,
#                          distance=stream_coord.distance)
#
#        #L = np.atleast_1d(stream_coord.Lambda.radian)
#        #B = np.atleast_1d(stream_coord.Beta.radian)
#        #
#        #Xs = np.cos(B)*np.cos(L)
#        #Ys = np.cos(B)*np.sin(L)
#        #Zs = np.sin(B)
#        #Zs = -Zs
#        #
#        #X, Y, Z = matrix.T.dot(np.array([Xs, Ys, Zs]))
#        #
#        #lon = np.arctan2(Y,X)*u.radian
#        #lat = np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z))*u.radian
#        #
#        #lon[lon<=0] += 2*np.pi*u.radian
#        #return coord.ICRS(ra=lon, dec=lat, distance=stream_coord.distance)
#
#    return cls

##class Sagittarius(Stream): pass
##class Magellanic(Stream): pass
# 
#COORDS = odict([
#        #('Sagittarius',
#        # # Euler angles (from Law & Majewski 2010)
#        # dict(pole = (5.6193, -14.0666),
#        #      euler = (180 + 3.75, 90 - 13.46, 180 + 14.111534),
#        #      ends  = [(45,5),(22.5,-7)],
#        #      )
#        # ),
#        ('Sagittarius',
#         # Euler angles and orbital pole from Majewski et al. 2003
#         # https://arxiv.org/abs/astro-ph/0304198
#         # Orbital pole:
#         #  l,b    = 273.8, -13.5
#         # Our convention has the opposite sign
#         # (right-handed vs left-handed coordinate systems?):
#         #  l,b    = 93.8, 13.5
#         dict(pole  = (273.8, -13.5),
#              euler = (273.8 - 90, -13.5 + 90, 180 + 14.1),
#              ends  = [(45,5),(22.5,-7)],
#              )
#         ),
#        ('Magellanic',
#         # Orbital pole from Nidever et al. 2008
#         # https://arxiv.org/abs/0706.1578
#         # l,b = 188.5, -7.5
#         # Euler angles and ends derived by hand...
#         dict(pole  = (188.5, -7.5),
#              euler = (188.5 - 90, -7.5 + 90, 180 + 32.72421421),
#              ends  = [(80,-72),(355,-2.5)],
#              )
#         ),
#        ('ATLAS',
#         # End points from galstreams (https://github.com/cmateu/galstreams)
#         # ra,dec = (19.3549, -26.6374),(31.1024, -32.8941)
#         # Orbital pole from Koposov et al. 2015
#         # ra,dec = 77.16, 46.92
#         # l,b = 161.03, 3.98
#         # Euler angles derived by hand...
#         dict(pole  = (161.03, 3.98),
#              euler = None,
#              #phi   = -90 + 161.03,
#              #theta = 90 + 3.98,
#              #psi   = 0,
#              ends  = [(12.0, -22.1), (29.9, -32.7)]
#              )
#         ),
#        ('GD-1',
#         # Endpoints eyeballed from Grillmair & Diontos 2006 [0604332]:
#         # ra,dec = (141,26.5),(214, 58.5),
#         # Orbital pole derived from:
#         # Koposov, Rix, Hogg (2010) [0907.1085]
#         # ra,dec = (34.5987, 29.7331)
#         # l,b = (144.6,-29.45)
#         # Then flipped to be at positive b
#         # l,b = (324.6,29.45)
#         dict(pole = (144.6, -29.45),
#              euler = None,
#              #phi   = -90 + 240.57,
#              #theta = 90 + 17.72,
#              #psi   = 0,
#              ends  = [(26.4, -45.8), (20.7, -54.2)]
#              )
#         ),
#        ('Pal5',
#         # Endpoints eyeballed from :
#         # Odenkirchen et al. 2003
#         # ra,dec = (227.5, -2.6), ( 233.45, 3.05)
#         # ra,dec = (140.23, 46.44)
#         # l,b = (173.25,44.64)
#         dict(pole = (173.25,44.64),
#              euler = None,
#              #phi   = -90 + 240.57,
#              #theta = 90 + 17.72,
#              #psi   = 0,
#              ends  = [(26.4, -45.8), (20.7, -54.2)]
#              )
#         ),
#        ('Styx',
#         # Endpoints from:
#         # Grillmair 2009 [0811.3965]
#         # ra,dec = (194,20), (259,21)
#         # ra,dec = (48.9,66,1)
#         # l,b = (136.8,7.2)
#         dict(pole = (136.8,7.2),
#              euler = None,
#              #phi   = -90 + 240.57,
#              #theta = 90 + 17.72,
#              #psi   = 0,
#              ends  = [(26.4, -45.8), (20.7, -54.2)]
#              )
#         ),
#        ('TriPsc',
#         # Endpoints from:
#         # Bonaca et al. 2012
#         # ra,dec = (21.346, 34.977),(23.975, 23.201)
#         # ra,dec = (118.1,9.6))
#         # l,b = (211.15,17.86)
#         dict(pole = (211.15,17.86),
#              euler = None,
#              #phi   = -90 + 240.57,
#              #theta = 90 + 17.72,
#              #psi   = 0,
#              ends  = [(26.4, -45.8), (20.7, -54.2)]
#              )
#         ),
#        ('Phoenix',
#         # Orbital pole derived from endpoints:
#         # RA,DEC = [(20.7, -54.2), (26.4, -45.8)]
#         # ra,dec = (132.1,-14.74)
#         # l,b = (240.57, 17.72)
#         # From Balbinot et al. (2015) [1509.04283]
#         # 72.247, -29.698
#         # End position
#         # Euler angles derived by hand...
#         dict(pole = (240.57, 17.72),
#              euler = None,
#              #phi   = -90 + 240.57,
#              #theta = 90 + 17.72,
#              #psi   = 0,
#              ends  = [(26.4, -45.8), (20.7, -54.2)]
#              )
#         ),
#        ('Tucana',
#         # Orbital pole derived from endpoints:
#         # RA,DEC =  [(4.0, -59.5), (-4.0, -59.8)]
#         # ra,dec = (175.09,-30.20)
#         # l,b = (285.16, 30.21)
#         # End position
#         # Euler angles derived by hand...
#         dict(ends  = [(-4.0, -59.8), (4.0, -59.5)],
#              )
#         ),
#        ('Indus',
#         # Orbital pole derived from endpoints.
#         dict(ends  = [(-28.7, -57.1),(-14.8, -62.5)],
#              )
#         ),
#        ('Eridanus',
#         # Orbital pole derived from endpoints.
#         dict(ends  = [(65.60, -21.19),(75.19, -26.44)],
#              )
#         ),
#        ('Phoenix2',
#         # Orbital pole derived from endpoints.
#         dict(ends  = [(14.44, -38.80),(22.82, -44.39)],
#              )
#         ),
# 
#        ])
# 
## Need to bind `cls` early:
## https://www.python-course.eu/global_vs_local_variables.php
#for name,coords in COORDS.items():
#    # Create the new class dynamically from the Stream base class
#    cls = type(name, (Stream,),{'__init__':Stream.__init__})
#    ends = np.array(coords['ends'])
#    cls.ends = ends
#    globals()[name] = cls
# 
#    # Pole in Galactic Coordinates
#    if not coords.get('pole'):
#        ra,dec = find_pole(ends[0][0],ends[0][1],ends[1][0],ends[1][1])
#        glon,glat = cel2gal(ra,dec)
#    else:
#        glon,glat = coords.get('pole')
# 
#    # Euler Angles
#    if not coords.get('euler'):
#        phi   = coords.get('phi'  , glon - 90)
#        theta = coords.get('theta', glat + 90)
#        psi   = coords.get('psi'  , 0)
#    else:
#        phi,theta,psi = coords['euler']
# 
# 
#    #cls = coords['cls']
# 
#    # Generate the rotation matrix using the x-convention (see Goldstein)
#    D = rotation_matrix(np.radians(phi),   "z", unit=u.radian)
#    C = rotation_matrix(np.radians(theta), "x", unit=u.radian)
#    B = rotation_matrix(np.radians(psi),   "z", unit=u.radian)
#    MATRIX = np.array(B.dot(C).dot(D))
# 
#    @frame_transform_graph.transform(coord.FunctionTransform, coord.Galactic, cls)
#    def galactic_to_stream(gal_coord, sgr_frame, cls=cls, matrix=MATRIX):
#        """ Compute the transformation from Galactic spherical to
#            heliocentric Sgr coordinates.
#        """
# 
#        l = np.atleast_1d(gal_coord.l.radian)
#        b = np.atleast_1d(gal_coord.b.radian)
# 
#        X = np.cos(b)*np.cos(l)
#        Y = np.cos(b)*np.sin(l)
#        Z = np.sin(b)
# 
#        # Calculate X,Y,Z,distance in the stream system
#        Xs, Ys, Zs = matrix.dot(np.array([X, Y, Z]))
#        Zs = -Zs
# 
#        # Calculate the angular coordinates lambda,beta
#        Lambda = np.arctan2(Ys,Xs)*u.radian
#        Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2.*np.pi*u.radian
#        Beta = np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs))*u.radian
#        #import pdb; pdb.set_trace()
#        return cls(Lambda=Lambda, Beta=Beta,
#                   distance=gal_coord.distance)
# 
#    @frame_transform_graph.transform(coord.FunctionTransform, cls, coord.Galactic)
#    def stream_to_galactic(sgr_coord, gal_frame, cls=cls, matrix=MATRIX):
#        """ Compute the transformation from heliocentric Sgr coordinates to
#            spherical Galactic.
#        """
#        L = np.atleast_1d(sgr_coord.Lambda.radian)
#        B = np.atleast_1d(sgr_coord.Beta.radian)
# 
#        Xs = np.cos(B)*np.cos(L)
#        Ys = np.cos(B)*np.sin(L)
#        Zs = np.sin(B)
#        Zs = -Zs
# 
#        X, Y, Z = matrix.T.dot(np.array([Xs, Ys, Zs]))
# 
#        l = np.arctan2(Y,X)*u.radian
#        b = np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z))*u.radian
# 
#        l[l<=0] += 2*np.pi*u.radian
#        return coord.Galactic(l=l, b=b, distance=sgr_coord.distance)




if __name__ == '__main__':
    # Sagittarius and LMC positions
    SGR_GLON,SGR_GLAT = 5.6193,-14.0666
    LMC_GLON,LMC_GLAT = 280.47, -32.75

    # Pole of Magellanic stream
    MS_GLON,MS_GLAT = 188.5, -7.5

    # Random spot on the sky?
    icrs = coord.ICRS(280.161732*u.degree, 11.91934*u.degree)
    print(icrs.transform_to(Sagittarius))
    print(icrs.transform_to(Magellanic))

    # Sgr Dwarf
    icrs = coord.Galactic(SGR_GLON*u.degree, SGR_GLAT*u.degree)
    print(icrs.transform_to(Sagittarius))

    # Magellanic Stream Pole
    icrs = coord.Galactic(MS_GLON*u.degree, MS_GLAT*u.degree)
    print(icrs.transform_to(Magellanic))
     
    # LMC
    icrs = coord.Galactic(LMC_GLON*u.degree, LMC_GLAT*u.degree)
    print(icrs.transform_to(Magellanic))

