import csv
import pprint

def compute_MODIP(latitude, longitude, stModip):
    """
    Reference: Section 2.5.4.3 of GSA's "Ionospheric Correction
    Algorithm for Galileo Single Frequency Users"

    :param latitude: [deg]
    :param longitude: [deg]
    :param stModip: array of MODIP vales
    :return: MODIP [deg]
    """

    if (latitude > 90) or (latitude < -90):
        mu = 90
        return mu

    l = int((longitude + 180)/10) - 2
    if l<0:
        l = l + 36
    if l>33:
        l = l - 36

    a = (latitude + 90) / 5.0 + 1
    x = a = int(a)
    i = int(a) - 2

    z = []
    for j in range(1,5):
        z_row = []
        for k in range(1,5):
            z_jk = stModip[i+j,l+k]
            z_row.append(z_jk)
        z.append(z_row)

    z_k = []
    for k in range(1,5):
        args = [z[0][k], z[1][k], z[2][k], z[3][k]]
        z_k.append(interpolate(args[0], args[1], args[2], args[3], x))

    b = (longitude + 180) / 10.0
    y = b = int(b)

    mu = interpolate(z_k[0], z_k[1], z_k[2], z_k[3], y)

    return mu

def interpolate(z1,z2,z3,z4,x):
    """
    Third Order Interpolation function
    Reference: Section 2.5.7.1 of GSA's "Ionospheric Correction
    Algorithm for Galileo Single Frequency Users"

    :param z1:
    :param z2:
    :param z3:
    :param z4:
    :param x:
    :return:
    """
    if abs(2*x) < 10**-10:
        return z2

    delta = 2*x - 1
    g1 = z3 + z2
    g2 = z3 - z2
    g3 = z4 + z1
    g4 = (z4 - z1) / 3.0

    a0 = 9 * g1 - g3
    a1 = 9 * g2 - g4
    a2 = g3 - g1
    a3 = g4 - g2

    return 1/16.0 * (a0 + a1*delta + a2*delta**2 + a3*delta**3)

def read_stMoDIP(path):
    with open(path) as f:
        reader = csv.reader(f, delimiter=' ')
        data = list(rec for rec in csv.reader(f, delimiter=','))
        data = [map(float, row) for row in data]
        return data


