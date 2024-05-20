import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

def plot_pretty(dpi=175, fontsize=15, labelsize=15, figsize=(10, 8), tex=True):
    # import pyplot and set some parameters to make plots prettier

    plt.rc('savefig', dpi=dpi)
    plt.rc('text', usetex=tex)
    plt.rc('font', size=fontsize)
    plt.rc('xtick.major', pad=5)
    plt.rc('xtick.minor', pad=5)
    plt.rc('ytick.major', pad=5)
    plt.rc('ytick.minor', pad=5)
    plt.rc('figure', figsize=figsize)

    mpl.rcParams['xtick.labelsize'] = labelsize
    mpl.rcParams['ytick.labelsize'] = labelsize
    mpl.rcParams.update({'figure.autolayout': True})

def mod2dist(distance_modulus):
    return 10**(distance_modulus / 5. + 1.) / 1e3


def dist2mod(distance):
    return 5. * (np.log10(np.array(distance) * 1.e3) - 1.)

def get_euler(stream=None, ends=None):
    import galstreams
    if ends is None:
        if isinstance(stream, str):
            mw_streams = galstreams.MWStreams(verbose=False)
            stream = mw_streams[stream]
        # ends = [[stream.end_o.ra.deg, stream.end_o.dec.deg], [stream.end_f.ra.deg, stream.end_f.dec.deg]]
        ends = [[stream.end_f.ra.deg, stream.end_f.dec.deg], [stream.end_o.ra.deg, stream.end_o.dec.deg]]

    fr = gc.GreatCircleICRSFrame.from_endpoints(coord.SkyCoord(ends[0][0], ends[0][1], unit='deg'), coord.SkyCoord(ends[1][0], ends[1][1], unit='deg'))
    phi, theta, psi = get_rot(fr)
    # phi, theta, psi = zyx_euler_from_endpoints(ends[0][0], ends[0][1], ends[1][0], ends[1][1])
    print(phi, theta, psi)
    return phi, theta, psi

def zyx_euler_from_endpoints(lon1, lat1, lon2, lat2):
    c1 = coord.SkyCoord(lon1 * u.deg, lat1 * u.deg)
    c2 = coord.SkyCoord(lon2 * u.deg, lat2 * u.deg)
    fr = gc.GreatCircleICRSFrame.from_endpoints(c1, c2)
    origin = fr.realize_frame(coord.UnitSphericalRepresentation(0 * u.deg, 0 * u.deg))

    gc_icrs = origin.transform_to(coord.ICRS)
    R = gc.greatcircle.reference_to_greatcircle(coord.ICRS, fr)
    psi = -np.degrees(np.arctan2(R[2, 1], R[2, 2]))

    return [gc_icrs.ra.degree, gc_icrs.dec.degree, psi]


def apwnorm(x, min=None, max=None):
    if min is None:
        min = np.nanmin(x[~np.isinf(x)])
    if max is None:
        max = np.nanmax(x[~np.isinf(x)])
    return np.clip((x - min) / (max - min), 0, 1)


def func_vec2pix(nside):
    return lambda x, y, z: hp.vec2pix(nside, x, y, z)