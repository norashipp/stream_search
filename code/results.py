#!/usr/bin/env python
"""
Calculate augmented results from a stream yaml file
"""
import os
import copy
from collections import OrderedDict as odict

import yaml
import numpy as np
import healpy as hp

from astropy.coordinates import SkyCoord, Galactocentric
from astropy.coordinates.matrix_utilities import rotation_matrix
import astropy.units as u


def create_result(stream):
    """ Create an augmented result dictionary from a single stream.

    Parameters:
    -----------
    stream : the dictionary for a single stream

    Returns:
    result : the augmented result dictionary
    """
    result = copy.deepcopy(stream)

    points = np.array(stream['ends'])[[0, -1]]
    result['endpoints'] = points.tolist()
    ra, dec = points.T

    # Pole in celestial coordinates
    pole_cel = find_pole(ra[0], dec[0], ra[-1], dec[-1])
    result['pole'] = pole_cel

    # Euler angles to stream coordinates
    euler_cel = euler_angles(ra[0], dec[0], ra[-1], dec[-1])
    result['euler_angles'] = euler_cel

    # Calculate the physical ditance
    distance_kpc = float(mod2dist(stream['modulus']))
    result['distance'] = distance_kpc
    # Pole in celestial coordinates
    coord = SkyCoord(ra * u.deg, dec * u.deg,
                     distance=distance_kpc * np.ones(len(ra)) * u.kpc,
                     frame='icrs')

    galacto = Galactocentric(galcen_distance=8.3 * u.kpc, z_sun=0. * u.pc)
    gal_xyz = coord.transform_to(galacto).cartesian.xyz.value.T

    result['gal_xyz'] = gal_xyz.tolist()
    v1 = gal_xyz[0]
    v2 = gal_xyz[-1]
    rgc1 = float(np.sqrt(np.sum(v1**2)))
    rgc2 = float(np.sqrt(np.sum(v2**2)))
    result['gal_dist'] = [rgc1, rgc2]
    result['gal_mean_dist'] = float(np.mean([rgc1, rgc2]))

    # phi = azimuthal
    # psi = theta = polar
    phi, psi = hp.vec2ang(np.cross(v1, v2), lonlat=True)

    phi[psi < 0] = np.mod(phi[psi < 0] + 180, 360)
    psi[psi < 0] *= -1
    psi = np.asscalar(90 - psi)
    phi = np.asscalar(phi)  # this is not the same as Erkal 2016, which had a sign error.

    result['gal_pole'] = [phi, psi]

    # Convert to physical parameters
    length_deg = float(angsep(ra[0], dec[0], ra[-1], dec[-1]))
    length_kpc = float(distance_kpc * np.radians(length_deg))
    width_deg = float(stream['width_deg'])
    width_kpc = float(distance_kpc * np.radians(width_deg))
    result['length_deg'] = length_deg
    result['length_kpc'] = length_kpc
    result['width_deg'] = width_deg
    result['width_pc'] = width_kpc * 1e3

    iso_kw = dict(name='Dotter2008', distance_modulus=stream['modulus'],
                  age=stream['age'],
                  metallicity=stream['metallicity'])
    richness = nstars2richness(stream['nstars'], **iso_kw)

    # Isochrone-derived parameters
    iso_params = isochrone_params(richness=richness, **iso_kw)
    iso_params['richness'] = richness
    abs_mag = iso_params['absolute_magnitude']
    mu = float(surface_brightness(abs_mag, distance_kpc, length_deg, width_deg))
    iso_params['surface_brightness'] = mu

    for k, v in iso_params.items():
        result[k] = v

    # Progenitor mass calculation
    result['progenitor_mass'] = float(progenitor_mass(stream))

    return result


def create_results(streams):
    """
    Create a results dictionary from the raw streams dictionary.

    Parameters:
    -----------
    streams : input streams dictionary or yaml filename

    Returns:
    --------
    results : dictionary of augmented stream results
    """

    if isinstance(streams, str):
        streams = yaml.load(open(streams, 'r'))
    results = copy.deepcopy(streams)

    for key, stream in results.items():
        results[key] = create_result(stream)

    return results


def mod2dist(distance_modulus):
    """
    Return distance (kpc) for a given distance modulus.
    """
    distance_modulus = np.array(distance_modulus)
    return 10**((0.2 * distance_modulus) - 2.)


def dist2mod(distance):
    """
    Return distance modulus for a given distance (kpc).
    """
    distance = np.array(distance)
    return 5. * (np.log10(distance * 1.e3) - 1.)


def angsep(lon1, lat1, lon2, lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance

    Parameters:
    -----------
    lon1 : longitude of the first point (deg)
    lat1 : latitude of the first point (deg)
    lon2 : longitude of the second point (deg)
    lat2 : latitude of the second point (deg)

    Returns:
    --------
    sep  : angular separation (deg)
    """
    lon1, lat1 = np.radians([lon1, lat1])
    lon2, lat2 = np.radians([lon2, lat2])

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1, num2), denominator))


def find_pole(lon1, lat1, lon2, lat2):
    """ Find the pole of a great circle orbit between two points.

    Parameters:
    -----------
    lon1 : longitude of the first point (deg)
    lat1 : latitude of the first point (deg)
    lon2 : longitude of the second point (deg)
    lat2 : latitude of the second point (deg)

    Returns:
    --------
    lon,lat : longitude and latitude of the pole
    """
    vec = np.cross(hp.ang2vec(lon1, lat1, lonlat=True),
                   hp.ang2vec(lon2, lat2, lonlat=True))
    lon, lat = hp.vec2ang(vec, lonlat=True)
    return [np.asscalar(lon), np.asscalar(lat)]


def create_matrix(phi, theta, psi):
    """ Create the transformation matrix.
    """
    # Generate the rotation matrix using the x-convention (see Goldstein)
    D = rotation_matrix(np.radians(phi),   "z", unit=u.radian)
    C = rotation_matrix(np.radians(theta), "x", unit=u.radian)
    B = rotation_matrix(np.radians(psi),   "z", unit=u.radian)
    return np.array(B.dot(C).dot(D))


def euler_angles(lon1, lat1, lon2, lat2, center=None):
    """ Calculate the Euler angles for spherical rotation using the x-convention
    (see Goldstein).

    Parameters:
    -----------
    lon1 : longitude of the first point (deg)
    lat1 : latitude of the first point (deg)
    lon2 : longitude of the second point (deg)
    lat2 : latitude of the second point (deg)

    Returns:
    --------
    phi,theta,psi : rotation angles around Z,X,Z
    """
    pole = find_pole(lon1, lat1, lon2, lat2)

    # Initial rotation
    phi = pole[0] - 90.
    theta = pole[1] + 90.
    psi = 0.

    matrix = create_matrix(phi, theta, psi)
    # Generate the rotation matrix using the x-convention (see Goldstein)
    #D = rotation_matrix(np.radians(phi),   "z", unit=u.radian)
    #C = rotation_matrix(np.radians(theta), "x", unit=u.radian)
    #B = rotation_matrix(np.radians(psi),   "z", unit=u.radian)
    #MATRIX = np.array(B.dot(C).dot(D))

    if center is not None:
        lon = np.radians([center[0]])
        lat = np.radians([center[1]])

        X = np.cos(lat) * np.cos(lon)
        Y = np.cos(lat) * np.sin(lon)
        Z = np.sin(lat)

        # Calculate X,Y,Z,distance in the stream system
        Xs, Ys, Zs = matrix.dot(np.array([X, Y, Z]))
        Zs = -Zs
        # print('no z flip')

        # Calculate the transformed longitude
        Lambda = np.arctan2(Ys, Xs)
        Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2. * np.pi
        psi = float(np.mean(np.degrees(Lambda)))

    else:
        lon = np.radians([lon1, lon2])
        lat = np.radians([lat1, lat2])

        X = np.cos(lat) * np.cos(lon)
        Y = np.cos(lat) * np.sin(lon)
        Z = np.sin(lat)

        # Calculate X,Y,Z,distance in the stream system
        Xs, Ys, Zs = matrix.dot(np.array([X, Y, Z]))
        Zs = -Zs

        # Calculate the transformed longitude
        Lambda = np.arctan2(Ys, Xs)
        Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2. * np.pi

        psi = float(np.mean(np.degrees(Lambda)))

    return [phi, theta, psi]


def nstars2richness(nstars, gmax=23.5, **kwargs):
    """ Calculate the richness from the number of stars in the MS selection

    Parameters:
    -----------
    nstars: number of (area corrected, not isochrone corrected) stars in the MS

    Returns:
    --------
    richness
    """

    richness = np.nan

    try:
        from ugali.isochrone import factory
        iso = factory(**kwargs)
        # Ok, hack back the representer added by model.py
        yaml.Dumper.yaml_representers.pop(odict, None)
    except ImportError:
        return params

    # Sample the isochrone/IMF
    init, pdf, act, mag_1, mag_2 = iso.sample(mass_steps=1e4)

    # NStars is calculated from an isochrone selection with 67% efficiency
    nstars *= (1 / 0.67)

    # Select just MS stars
    MSTO = 3.5
    gmin = MSTO + iso.distance_modulus
    g = mag_1 + iso.distance_modulus

    # Calculate the richness the fraction of the total number of stars
    sel = (g > gmin) & (g < gmax)
    richness = nstars / (np.sum(pdf[sel]) / pdf.sum())

    return float(richness)


def isochrone_params(richness, **kwargs):
    """ Calculate isochrone parameters.

    NOTE: This has a dependence on ugali:
    https://github.com/DarkEnergySurvey/ugali

    Parameters:
    -----------
    richness : isochrone richness (normalization) parameter
    kwargs   : keyword arguments (age, metallicity, distance_modulus, etc.)

    Returns:
    --------
    params   : dictionary of isochrone-derived parameters
    """

    params = dict(stellar_mass=np.nan,
                  nstars_new=np.nan,
                  absolute_magnitude=np.nan)

    try:
        from ugali.isochrone import factory
        iso = factory(**kwargs)
        # Ok, hack back the representer added by model.py
        yaml.Dumper.yaml_representers.pop(odict, None)
    except ImportError:
        return params

    params['stellar_mass'] = float(richness * iso.stellar_mass())
    mass_init, mass_pdf, mass_act, mag_1, mag_2 = iso.sample()

    # This calculates just main sequence stars
    MSTO = 3.5
    gmin = MSTO + iso.distance_modulus
    gmax = 23.5
    g = mag_1 + iso.distance_modulus
    sel = (g > gmin) & (g < gmax)
    nstars = richness * np.sum(mass_pdf[sel]) / mass_pdf.sum()
    params['nstars_new'] = float(nstars)
    params['absolute_magnitude'] = float(iso.absolute_magnitude(richness))

    return params


def progenitor_mass(stream):
    """ Returns progenitor mass in units of Msun. """
    try:
        from progenitor_masses import progenitor_masses_worker
    except ImportError:
        return np.nan

    return progenitor_masses_worker(stream)


def surface_brightness(abs_mag, distance, length_deg, width_deg):
    """ Calculate the surface brightness of a stream

    Parameters:
    -----------
    abs_mag     : V-band absolute magnitude (Mv)
    distance    : distance to the stream (kpc)
    length_kpc  : full length (deg)
    width_deg   : standard deviation of best-fit Gaussian (deg)

    Returns:
    --------
    surface_brightness : surface brightness (mag arcsec^-2)
    """
    length_asec = length_deg * 3600.
    width_asec = width_deg * 3600.

    # We are assuming that width is the std of a gaussian fit to the stream
    # and that 34.1% of the flux is contained within this angle (half-width)
    c_v = -2.5 * np.log10(0.341)

    return abs_mag + dist2mod(distance) + c_v + 2.5 * np.log10(length_asec * width_asec)


def onstream(length, width):
    verts = [
        [length / 2., -2 * width],
        [length / 2.,  2 * width],
        [-length / 2.,  2 * width],
        [-length / 2., -2 * width],
        [length / 2., -2 * width],
    ]
    return [np.array(verts)]


def offstream(length, width):
    on = onstream(length, width)[0]

    off1 = np.copy(on)
    off1[:, 1] += 4 * width

    off2 = np.copy(on)
    off2[:, 1] -= 4 * width

    return [off1, off2]


def angsep(lon1, lat1, lon2, lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1, lat1 = np.radians([lon1, lat1])
    lon2, lat2 = np.radians([lon2, lat2])

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1, num2), denominator))


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filename')
    parser.add_argument('outfile', nargs='?')
    parser.add_argument('-s', '--stream', default=None)
    args = parser.parse_args()

    results = create_results(args.filename)
    if args.stream:
        results = results.get(args.stream)

    if args.outfile:
        with open(args.outfile, 'w') as out:
            out.write(yaml.dump(results))
    else:
        print(yaml.dump(results))
