"""
This file contains functions that computes the 13 parameters used by NeQuick-G model.
NeQuick-G is an ionospheric model proposed for Galileo GNSS to correct for ionospheric time delay
"""

import numpy as np
import csv
f = open('dump','w')
dumpwriter = csv.writer(f, delimiter = ',')

def main():
    ### Local system
    localpath = '/home/tpl/Documents/Airbus/Project/Papers/Nequick/CCIR_MoDIP/'

    # STAGE 0 ##################################################################
    # Receiver knowledge
    lat = 30
    longitude = 0
    UT = 12
    mth = 10 # mth = {01, 02, 03, ... 11, 12}
    # Broadcast correction parameters
    ai0 = 80
    ai1 = 0
    ai2 = 0

    stage0_data = [lat, longitude, UT, mth, ai0, ai1, ai2]

    # TODO: find example broadcast correction parameters

    # STAGE 1 ##################################################################
    stMODIP = read_stMoDIP(localpath)
    MODIP = compute_MODIP(lat, longitude, stMODIP)
    Az = effective_ionization(ai0, ai1, ai2, MODIP)
    Azr = effective_sunspot_number(Az)
    solarcosine, solarsine = solar_declination(mth, UT)
    LT = localtime(UT, longitude)
    chi = solar_zenith(lat, LT, solarcosine, solarsine)
    chi_eff = effective_solar_zenith(chi)

    stage1_data = [MODIP, Az, Azr, chi_eff]
    print "stage1_data = [MODIP, Az, Azr, chi_eff]"
    print stage1_data
    # STAGE 2 ##################################################################

    F2, Fm3 = readccirXXfiles(mth, localpath)
    AF2, Am3 = interpolate_AZR(F2, Fm3, Azr)
    cf2, cm3 = F2fouriertimeseries(UT, AF2, Am3)
    foF2, M3000F2, NmF2 = legendre_calculation(MODIP, lat, longitude, cf2, cm3)
    foE, NmE = ELayer(lat, Az, chi_eff, mth)
    foF1, NmF1 = F1Layer(foE, NmF2)

    stage2a_data = [foE, foF1, foF2]
    stage2b_data = [NmE, NmF1, NmF2]
    stage2c_data = [M3000F2]

    stage2_data = stage2a_data + stage2b_data + stage2c_data

    print "stage 2a_data"
    print stage2a_data

    # STAGE 3 ##################################################################

    hmE = get_hmE()
    hmF2 = get_hmF2(foE, foF2, M3000F2)
    hmF1 = get_hmF1(hmF2, hmE)

    stage3_data = [hmE, hmF1, hmF2]

    # STAGE 4 ##################################################################
    B2bot = get_B2bot(NmF2, foF2, M3000F2)
    B1top = get_B1top(hmF1, hmF2)
    B1bot = get_B1bot(hmF1, hmE)
    BEtop = get_BEtop(B1bot)
    BEbot = get_BEbot()

    stage4_data = [BEtop, BEtop, B1top, B1bot, B2bot]

    # STAGE 5 ##################################################################

    A1 = get_A1(NmF2)
    A2, A3 = get_A2A3(NmE, NmF1, A1, hmF2, hmF1, hmE, BEtop, B1bot, B2bot, foF1)
    k = shape_parameter(mth, NmF2, hmF2, B2bot, Azr)
    H0 = get_H0(B2bot, k)

    stage5_data = [k, A1, A2, A3, H0]

    # OUTPUT ##################################################################

    bottomside_para = [hmE, hmF1, hmF2, BEtop, BEbot, B1top, B1bot, B2bot, A1, A2, A3]
    topside_para = [NmF2, hmF2, H0]

    print "Nm:"
    print [NmE, NmF1, NmF2]

    print "hm:"
    print [hmE, hmF1, hmF2]

    print "B:"
    print [BEtop, BEbot, B1top, B1bot, B2bot]

    print "A:"
    print [A1, A2, A3]

    return bottomside_para, topside_para, Azr


################################################################
# STAGE 1:
def read_stMoDIP(path):
    path = path + 'modipNeQG_wrapped.txt'
    with open(path) as f:
        data = list(rec for rec in csv.reader(f, delimiter=','))
        data = [map(float, row) for row in data]
        dumpwriter.writerows(data)
        return data


def interpolate(z1, z2, z3, z4, x):
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
    if abs(2 * x) < 10 ** -10:
        return z2

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

    l = int((longitude + 180) / 10) - 2
    if l < 0:
        l = l + 36
    if l > 33:
        l = l - 36
    assert (type(l) == int)
    a = (latitude + 90) / 5.0 + 1
    x = a - int(a)
    i = int(a) - 2

    z = []
    for j in range(4):
        z_row = []
        for k in range(4):
            z_jk = stModip[i + j][l + k]
            z_row.append(z_jk)
        z.append(z_row)

    z_k = []
    for k in range(4):
        args = [z[0][k], z[1][k], z[2][k], z[3][k]]
        z_k.append(interpolate(args[0], args[1], args[2], args[3], x))

    b = (longitude + 180) / 10.0
    y = b - int(b)
    mu = interpolate(z_k[0], z_k[1], z_k[2], z_k[3], y)

    return mu


def effective_ionization(ai0, ai1, ai2, MoDIP):
    """
    :param ai0:
    :param ai1:
    :param ai2:
    :param MoDIP: [deg]
    :return: Effective ionoiszation level
    Remark: Az is equivalent to F10.7 in climatological NeQuick

    """
    if (ai0 == 0) and (ai1 == 0) and (ai2 == 0):
        Az = 63.7
    else:
        Az = ai0 + ai1 * MoDIP + ai2 * MoDIP ** 2

    # Reference Section 3.3
    assert Az <= 400
    assert Az >= 0

    return Az


def effective_sunspot_number(Az):
    """
    This parameter is equivalent to R 12 in climatological NeQuick
    Reference: Section 2.5.4.5
    :param Az:
    :return:

    """
    return np.sqrt(167273 + (Az - 63.7) * 1123.6) - 408.99


def solar_declination(month, universal_time):
    """
    Compute sin(delta_Sun ), cos(delta_Sun ), the sine and cosine of the solar declination.
    :param month: [mth]
    :param universal_time: [hours and decimals]
    :return:(Cosine, Sine)
    """
    # Compute day of year at the middle of the month
    dy = 30.5 * month - 15
    # Compute time [days]:
    t = dy + (18 - universal_time) / 24
    # Compute the argument
    a_m = (0.9856 * t - 3.289) * (np.pi / 180)  # radians
    a_l = a_m + (1.916 * np.sin(a_m) + 0.020 * np.cos(2 * a_m) + 282.634) * (np.pi / 180)  # radians

    # Compute sine and cosine of solar declination

    Sine = 0.39782 * np.sin(a_l)
    Cosine = np.sqrt(1 - Sine ** 2)
    return (Cosine, Sine)


def localtime(universal_time, longitude):
    """

    :param universal_time:
    :return:
    """
    return universal_time + longitude / 15.0


def solar_zenith(latitude, localtime, solarcosine, solarsine):
    """
    Reference : Section 2.5.4.7
    :param latitude: [deg]
    :param localtime: [hours and decimals]
    :param solarcosine: of solar declination
    :param solarsine: of solar declination
    :return: solar zenith angle [deg]
    """
    coschi = np.sin(latitude * np.pi / 180) * solarsine + np.cos(latitude * np.pi / 180) * solarcosine * np.cos(
        np.pi / 12 * (12 - localtime))
    chi = np.arctan2(np.sqrt(1 - coschi ** 2), coschi)

    return chi


def effective_solar_zenith(chi, chi0=86.23292796211615):
    """
    Reference: Section 2.5.4.8
    :param chi: solar zenith angle [deg]
    :param chi0: default solar zenith angle [deg]
    :return: effective solar zenith angle[deg]
    """
    # numerator = chi + (90 - 0.24 * np.exp(20 - 0.2 * chi)) * np.exp(12 * (chi - chi0))
    # denominator = 1.0 + np.exp(12 * (chi - chi0))
    # return numerator / denominator

    return NeqJoin(90 - 0.24 * NeqClipExp(20 - 0.2 * chi), chi, 12, chi - chi0)


#########################################################################

# STAGE 2:
def ELayer(latitude, Az, chi_eff, month):
    """
    Reference: 2.5.5.1
    :param latitude: [deg]
    :param Az: Effective Ionisation Level
    :param chi_eff: effective_solar_zenith [deg]
    :param month: mth
    :return: E layer critical frequency foE [MHz], NmE [10^11 m^-3]
    """
    # Define the seas parameter as a function of the month of the year as follows
    if (month in [1, 2, 11, 12]):
        seas = -1
    elif (month in [3, 4, 9, 10]):
        seas = 0
    elif (month in [5, 6, 7, 8]):
        seas = 1
    else:
        raise ValueError('Month must be an integer between 1 and 12')

    # Introduce the latitudinal dependence
    ee = NeqClipExp(0.3 * latitude)
    seasp = seas * (ee - 1) / (ee + 1)

    foE = np.sqrt((1.112 - 0.019 * seasp) ** 2 * np.sqrt(Az) * np.cos(chi_eff * np.pi / 180) ** 0.6 + 0.49)

    NmE = NeqCriticalFreqToNe(foE)

    return foE, NmE


def readccirXXfiles(month, path):
    """
    Reference: Section 2.5.3.2
    :param month:
    :param path:
    :return:f2_ijk array and fm3_ijk
    """

    data = []
    with open(path + '/ccir' + str(month + 10) + '.txt') as f:
        for row in csv.reader(f, delimiter=' '):
            row = [num for num in row if num != '']  # filter
            data = data + [float(num) for num in row]
    assert (len(data) == 2858)

    F2data = data[:1976]
    F2 = np.reshape(np.array(F2data), (2, 76, 13))

    Fm3data = data[1976:]
    Fm3 = np.reshape(np.array(Fm3data), (2, 49, 9))


    return F2, Fm3


def interpolate_AZR(F2, Fm3, Azr):
    """

    :param F2:
    :param Fm3:
    :param Azr:
    :return:
    """

    AF2 = F2[0, :, :] * (1 - Azr / 100.0) + F2[1, :, :] * Azr / 100.0
    Am3 = Fm3[0, :, :] * (1 - Azr / 100.0) + Fm3[1, :, :] * Azr / 100.0

    return AF2, Am3


def F2fouriertimeseries(UT, AF2, Am3):
    """

    :param UT: Universal Time UT [hours],
    :param AF2:
    :param Am3:
    :return:CF2, Cm3 vectors of coefficients for Legendre calculation for foF2 and M(3000)F2
    """

    # Compute the time argument
    T = (15 * UT - 180) * np.pi / 180.0

    # calculate the Fourier time series for foF2
    cos = np.cos(np.arange(1, 7) * T)
    sin = np.sin(np.arange(1, 7) * T)
    x = np.empty(12)
    x[0::2] = sin
    x[1::2] = cos
    y = np.ones(13)
    y[1:] = x

    cf2 = np.sum(AF2 * y, axis=1)

    # calculate the Fourier time series for M(3000)F2
    cos = np.cos(np.arange(1, 5) * T)
    sin = np.sin(np.arange(1, 5) * T)
    x = np.empty(8)
    x[0::2] = sin
    x[1::2] = cos
    y = np.ones(9)
    y[1:] = x

    cm3 = np.sum(Am3 * y, axis=1)
    dumpwriter.writerow(cf2)
    dumpwriter.writerow(cm3)

    return cf2, cm3


def legendre_calculation(modip, latitude, longitude, CF2, Cm3):
    """

    :param modip:
    :param latitude:
    :param longitude:
    :param CF2:
    :param Cm3:
    :return: foF2 [MHz], M(3000)F2
    """

    M = np.ones(12)
    P = np.ones(9)
    S = np.ones(9)
    C = np.ones(9)

    # Compute MODIP coefficients
    factor = np.sin(modip * np.pi / 180)
    for i in range(1, 12):
        M[i] = factor * M[i - 1]

    # Compute latitude and longitude coefficients
    factor = np.cos(latitude * np.pi / 180)
    for i in range(1, 9):
        P[i] = factor * P[i - 1]
        S[i] = np.sin(i * longitude * np.pi / 180)
        C[i] = np.cos(i * longitude * np.pi / 180)
    # Compute foF2
    # Order 0 term
    foF2_1 = np.sum(CF2[:12] * M)
    # Legendre grades
    Q = np.array([12, 12, 9, 5, 2, 1, 1, 1, 1])
    K = np.empty(9 , dtype = np.int)
    K[0] = -Q[0]
    for i   in range(1, 9):
        K[i] = K[i - 1] + 2 * Q[i - 1]  # [-12.,  12.,  36.,  54.,  64.,  68.,  70.,  72.,  74.]

    foF2_n = np.zeros(9)
    foF2_n[0] = foF2_1
    for i in range(1, 9):  # check! there is a bug in the indices
        for j in range(Q[i]):
            foF2_n[i] += CF2[K[i] + 2 * j] * C[i] * M[j] * P[i]
            # print CF2[K[i] + 2 * j - 1]
            foF2_n[i] += CF2[K[i] + 2 * j+1] * S[i] * M[j] * P[i]
    # print foF2_n
    foF2 = np.sum(foF2_n)

    # compute 0 order term
    M3000F2_1 = np.sum(Cm3[:7] * M[:7])

    R = np.array([7, 8, 6, 3, 2, 1, 1])
    H = np.empty(7, dtype = np.int)
    H[0] = -R[0]
    for i in range(1, 7):
        H[i] = H[i - 1] + 2 * R[i - 1]  # [ -7.,   7.,  23.,  35.,  41.,  45.,  47.])

    M3000F2_n = np.zeros(7)
    M3000F2_n[0] = M3000F2_1

    for i in range(1, 7):
        for j in range(R[i]):
            M3000F2_n[i] += (Cm3[H[i] + 2 * j - 1] * C[i] + Cm3[H[i] + 2 * j] * S[i]) * (M[j] * P[j])

    M3000F2 = np.sum(M3000F2_n)

    NmF2 = NeqCriticalFreqToNe(foF2)

    return foF2, M3000F2, NmF2


def F1Layer(foE, foF2):
    """

    :param foE: E layer critical frequency [MHz]
    :param foF2: F2 layer critical frequency [MHz]
    :return: (foF1, NmF1 [10^11 m^-3] )
    """
    if foE >= 2.0:
        foF1 = 1.4 * foE  # Titheridge's Formula
    else:
        foF1 = 0

    foF1 = NeqJoin(1.4 * foE, 0, 1000.0, foE - 2)
    foF1 = NeqJoin(0, foF1, 1000.0, foE - foF1)
    foF1 = NeqJoin(foF1, 0.85 * foF1, 60.0, 0.85 * foF2 - foF1)

    if foF1 < 10 ** -6:
        foF1 = 0

    # F1 layer maximum density
    if (foF1 <= 0) and (foE > 2):
        NmF1 = NeqCriticalFreqToNe(foE + 0.5)
    else:
        NmF1 = NeqCriticalFreqToNe(foF1)

    return foF1, NmF1


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
    if dPower > 80:
        return 5.5406 * 10 ** 34
    elif dPower < -80:
        return 1.8049 * 10 ** -35
    else:
        return np.exp(dPower)

def NeqCriticalFreqToNe(f0):
    """

    :param f0: peak plasma frequency of layer [MHz]
    :return:
    """
    return 0.124 * f0**2


#########################################################################

def get_hmE():
    return 120  # [km]


def get_hmF1(hmF2, hmE):
    """

    :param hmF2: [km]
    :param hmE: [km]
    :return:
    """
    hmF1 = (hmF2 + hmE) / 2.0

    return hmF1


def get_hmF2(foE, foF2, M3000F2):
    """

    :param foE: [Mhz]
    :param foF2: [Mhz]
    :param M3000F2:
    :return: [km]
    Based on Dudeney (1983)
    """

    numerator = 1490 * M3000F2 * np.sqrt((0.0196 * M3000F2 ** 2 + 1) / (1.2967 * M3000F2 ** 2 - 1))
    if foE < 10 ** -30:
        deltaM = - 0.012
    else:
        r = float(foF2) / foE
        top = r * np.exp(20 * (r - 1.75)) + 1.75
        bottom = np.exp(20 * (r - 1.75)) + 1.0
        rho = top / bottom
        deltaM = 0.253 / (rho - 1.215) - 0.012

    denominator = M3000F2 + deltaM

    hmF2 = numerator / denominator - 176

    return hmF2


#########################################################################

def get_B2bot(NmF2, foF2, M3000F2):
    """

    :param NmF2:[10^11 m^-3 ]
    :param foF2: [MHz]
    :param M3000F2:
    :return: [km]
    """

    top = 0.385 * NmF2
    bottom = 0.01 * np.exp(-3.467 + 0.857 * np.log(foF2 ** 2) + 2.02 * np.log(M3000F2))

    return top / bottom


def get_B1top(hmF1, hmF2):
    """

    :param hmF1: [km]
    :param hmF2: [km]
    :return: [km]
    """
    return 0.3 * (hmF2 - hmF1)


def get_B1bot(hmF1, hmE):
    """

    :param hmF1: [km]
    :param hmE: [km]
    :return: [km]
    """
    return 0.5 * (hmF1 - hmE)


def get_BEtop(B1bot):
    """

    :param B1bot: [km]
    :return: [km]
    """
    return max(B1bot, 7)


def get_BEbot():
    """

    :return: [km]
    """
    return 5.0


#########################################################################

def get_A1(NmF2):
    """

    :param NmF2:  [10^11 m^-3 ]
    :return: F2 layer amplitude
    """
    return 4 * NmF2


def get_A2A3(NmE, NmF1, A1, hmF2, hmF1, hmE, BEtop, B1bot, B2bot, foF1):
    """

    :param NmE:
    :param NmF1:
    :param A1:
    :param hmF2:
    :param hmF1:
    :param hmE:
    :param BEtop:
    :param B1bot:
    :param B2bot:
    :param foF1:
    :return:
    """
    if foF1 < 0.5:
        A2 = 0.0
        A3 = 4.0 * (NmE - epstein(A1, hmF2, B2bot, hmE))
    else:
        A3a = 4.0 * NmE
        for i in range(5):
            A2a = 4.0 * (NmF1 - epstein(A1, hmF2, B2bot, hmF1) - epstein(A3a, hmE, BEtop, hmF1))
            # A2a = (A2a * np.exp(A2a - 0.8 * NmF1) + 0.8 * NmF1) / (1 + np.exp(A2a - 0.8 * NmF1))
            A2a = NeqJoin(A2a, 0.8 * NmF1, 1, A2a - 0.8 * NmF1)
            A3a = 4.0 * (NmE - epstein(A2a, hmF1, B1bot, hmE) - epstein(A1, hmF2, B2bot, hmE))
        A2 = A2a
        # A3 = (A3a * np.exp(60 * (A3a - 0.005)) + 0.05) / (1 + np.exp(60 * (A3a - 0.005)))
        A3 = NeqJoin(A3a, 0.05, 60.0, A3a - 0.005)
    return A2, A3


def epstein(peak_amp, peak_height, thickness, H):
    return peak_amp * np.exp((H - peak_height) / thickness) / (1 + np.exp((H - peak_height) / thickness)) ** 2


def shape_parameter(mth, NmF2, hmF2, B2bot, Azr):
    """

    :param mth:
    :param NmF2:
    :param hmF2:
    :param B2bot:
    :param Azr:
    :return:
    """
    if mth in [4, 5, 6, 7, 8, 9]:
        ka = 6.705 - 0.014 * Azr - 0.008 * hmF2
    elif mth in [1, 2, 3, 10, 11, 12]:
        ka = -7.77 + 0.097 * (hmF2 / B2bot) ** 2 + 0.153 * NmF2
    else:
        raise ValueError("Invalid Month")
    kb = (ka * np.exp(ka - 2) + 2) / (1 + np.exp(ka - 2))
    k = (8 * np.exp(kb - 8) + kb) / (1 + np.exp(kb - 8))
    return k


def get_H0(B2bot, k):
    """

    :param B2bot: [km]
    :param k:
    :return: topside thickness parameter [km]
    """

    Ha = k * B2bot
    x = (Ha - 150.0) / 100.0
    v = (0.041163 * x - 0.183981) * x + 1.424472
    Ho = Ha / v
    return Ho
