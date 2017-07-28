import numpy as np

def coord2cartesian(r, lat, lon):
    """

    :param r: [km]
    :param lat: [deg]
    :param lon: [deg]
    :return:
    """
    DR = np.pi/180
    x = r * np.cos(lat * DR) * np.cos(lon * DR)
    y = r * np.cos(lat * DR) * np.sin(lon * DR)
    z = r * np.sin(lat * DR)

    return x, y, z

def cartesian2coord(x, y, z):
    """

    :param x: [km]
    :param y: [km]
    :param z: [km]
    :return:
    """
    r = np.sqrt(x**2 + y**2 + z**2)

    xy = np.sqrt(x**2 + y**2)

    lat = np.arctan(z / xy) * 180 / np.pi
    lon = np.arctan2(y, x) * 180 / np.pi

    return r, lat, lon



def interpolate(z1, z2, z3, z4, x):
    """
    Third Order Lagrange Interpolation function
    Reference: Section 2.5.7.1 of GSA's "Ionospheric Correction
    Algorithm for Galileo Single Frequency Users"

    :param z1:
    :param z2:
    :param z3:
    :param z4:
    :param x:
    :return:
    """
    # if abs(2 * x) < 10 ** -10:
    #     return z2

    delta = 2 * x - 1
    g1 = z3 + z2
    g2 = z3 - z2
    g3 = z4 + z1
    g4 = (z4 - z1) / 3.0

    a0 = 9 * g1 - g3
    a1 = 9 * g2 - g4
    a2 = g3 - g1
    a3 = g4 - g2

    return 1 / 16.0 * (a0 + a1 * delta + a2 * delta ** 2 + a3 * delta ** 3)

def interpolate2d(Z, x, y):
    assert (np.shape(Z) == (4,4))

    deltax = 2 * x - 1
    deltay = 2 * y - 1
    # Interpolate horizontally first

    G1 = Z[2,:] + Z[1,:]
    G2 = Z[2,:] - Z[1,:]
    G3 = Z[3,:] + Z[0,:]
    G4 = (Z[3,:] - Z[0,:]) / 3.0

    A0 = 9 * G1 - G3
    A1 = 9 * G2 - G4
    A2 = G3 - G1
    A3 = G4 - G2

    z = 1 / 16.0 * (A0 + A1 * deltay + A2 * deltay ** 2 + A3 * deltay ** 3)

    g1 = z[2] + z[1]
    g2 = z[2] - z[1]
    g3 = z[3] + z[0]
    g4 = (z[3] - z[0]) / 3.0

    a0 = 9 * g1 - g3
    a1 = 9 * g2 - g4
    a2 = g3 - g1
    a3 = g4 - g2

    return 1 / 16.0 * (a0 + a1 * deltax + a2 * deltax ** 2 + a3 * deltax ** 3)

def epstein(peak_amp, peak_height, thickness, H):
    return peak_amp * NeqClipExp((H - peak_height) / thickness) / np.power((1 + NeqClipExp((H - peak_height) / thickness)), 2)


def NeqJoin(dF1, dF2, dAlpha, dX):
    """
    Allows smooth joining of functions f1 and f2 (i.e. continuous first derivatives) at origin.
    Alpha determines width of transition region. Calculates value of joined functions at x.
    :param dF1:
    :param dF2:
    :param dAlpha:
    :param dX:
    :return:
    """
    ee = NeqClipExp(dAlpha * dX)
    return (dF1 * ee + dF2) / (ee + 1)


def NeqClipExp(dPower):
    """

    :param dPower: Power for exponential function [double]
    :return:
    """

    assert(not np.any(np.isnan(dPower)))

    mask1 = np.logical_and(dPower < 80, dPower > -80)
    mask2 = dPower > 80
    mask3 = dPower < -80
    out = np.exp(dPower, where=mask1)
    if type(out) == np.ndarray:
        out[mask2] = 5.5406 * 10 ** 34
        out[mask3]= 1.8049 * 10 ** -35
    else:
        if mask2:
            out = 5.5406 * 10 ** 34
        elif mask3:
            out = 1.8049 * 10 ** -35


    assert( not np.any(out < 0))
    return out


def NeqCriticalFreqToNe(f0):
    """

    :param f0: peak plasma frequency of layer [MHz]
    :return:
    """
    return 0.124 * f0 ** 2


if __name__ == "__main__":
    # unit testing
    assert NeqClipExp(-100) == 1.8049 * 10 ** -35
    assert NeqClipExp(100) == 5.5406 * 10 ** 34

    assert np.all(NeqClipExp(np.array([-100, 0, 100])) == np.array([1.8049 * 10 ** -35, 1., 5.5406 * 10 ** 34]))
