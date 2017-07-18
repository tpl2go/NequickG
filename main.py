import numpy as np
import csv
def group_delay(f, STEC):
    """

    :param f: frequency [Hz]
    :param STEC: Slanted Total Electron Content [electrons / m^2]
    :return: [m]
    """
    return 40.3 / f ** 2 * STEC


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
    :param universal_time: [hour]
    :return:(Cosine, Sine)
    """
    # Compute day of year at the middle of the month
    dy = 30.5*month -15
    # Compute time [days]:
    t = dy + (18-universal_time)/24
    # Compute the argument
    a_m = (0.9856 * t - 3.289) * (np.pi/180) # radians
    a_l = a_m + (1.916 * np.sin(a_m) + 0.020 * np.cos(2*a_m) + 282.634) * (np.pi/180) # radians

    # Compute sine and cosine of solar declination

    Sine = 0.39782 * np.sin(a_l)
    Cosine = np.sqrt(1-Sine**2)
    return (Cosine, Sine)

def solar_zenith(latitude, localtime, solarsine, solarcosine):
    """
    Reference : Section 2.5.4.7
    :param latitude: [deg]
    :param localtime: [hours]
    :param solarsine: of solar declination
    :param solarcosine: of solar declination
    :return: solar zenith angle [deg]
    """
    coschi = np.sin(latitude * np.pi/180) * solarsine + np.cos(latitude * np.pi/180) * solarcosine * np.cos(np.pi/12 * (12 - localtime))
    chi = np.arctan2(np.sqrt(1 - coschi**2), coschi)
    return chi

def effective_solar_zenith(chi, chi0=86.23292796211615):
    """
    Reference: Section 2.5.4.8
    :param chi: solar zenith angle [deg]
    :param chi0: default solar zenith angle [deg]
    :return: effective solar zenith angle[deg]
    """
    numerator = chi + (90 -0.24 *np.exp(20 - 0.2 * chi )) * np.exp(12*(chi-chi0))
    denominator = 1.0 + np.exp(12 * (chi - chi0))
    return numerator / denominator

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
    if (month in [1,2,11,12]):
        seas = -1
    elif (month in [3,4,9,10]):
        seas = 0
    elif (month in [5,6,7,8]):
        seas = 1
    else:
        raise ValueError('Month must be an integer between 1 and 12')

    #Introduce the latitudinal dependence
    ee = np.exp(0.3 * latitude)
    seasp = seas * (ee - 1) / (ee + 1)

    foE = np.sqrt((1.112 - 0.019*seasp)**2 * np.sqrt(Az) * np.cos(chi_eff * np.pi/180)**0.6 + 0.49)

    NmE = 0.124 * foE**2

    return foE, NmE

def F1Layer(foE, foF2):
    """

    :param foE: E layer critical frequency [MHz]
    :param foF2: F2 layer critical frequency [MHz]
    :return: (foF1, NmF1 [10^11 m^-3] )
    """
    if foE >= 2.0:
        foF1 = 1.4 * foE
    else:
        foF1 = 0
    # TODO: tweak foF1 with NeQJOIN (Section F.2.12.1)
    if foF1 < 10**-6:
        foF1 = 0

    #F1 layer maximum density
    if (foF1<=0) and (foE>2):
        NmF1 = 0.124 * (foE + 0.5)**2
    else:
        NmF1 = 0.124 * foF1**2

    return foF1, NmF1

def F2Layer(month):
    """
    Reference: 2.5.5.3
    :param month:
    :return:
    """
    # Read ccirXX.asc values
    #TODO
    pass
    F2, Fm3 = readccirXXfiles(month, '/home/tpl/Documents/Airbus/Project/Papers/Nequick/CCIR_MoDIP')

def readccirXXfiles(month, path):
    """
    Reference: Section 2.5.3.2
    :param month:
    :param path:
    :return:f2_ijk array and fm3_ijk
    """

    F2 = np.zeros((2,76, 13))
    Fm3 = np.zeros()

    data = []
    with open(path + '/ccir' + str(month + 10) + '.txt') as f:
        for row in csv.reader(f, delimiter = ' '):
            data = data + [float(num) for num in row]

    F2data = data[:1976]
    F2 = np.reshape(np.array(F2data),(2,76,13))

    Fm3data = data[1976:]
    Fm3 = np.reshape(np.array(Fm3data),(2,49,9))

    return (F2, Fm3)

def interpolate_AZR(F2, Fm3, Azr):
    """
    TODO: I DONT UNDERSTAND WHAT THIS DOES AT ALL
    :return:
    """
    AF2 = F2[0,:,:] * (1 - Azr/100.0) + F2[1,:,:] * Azr/100.0
    Am3 = Fm3[0,:,:] * (1 - Azr/100.0) + Fm3[1,:,:] * Azr/100.0

    return (AF2,Am3)

def F2fouriertimeseries(time, AF2, Am3):
    """

    :param time: Universal Time UT [hours],
    :param AF2:
    :param Am3:
    :return:CF2, Cm3 vectors of coefficients for Legendre calculation for foF2 and M(3000)F2
    """

    # Compute the time argument
    T = (15 * time - 180 ) * np.pi/180

    #calculate the Fourier time series for foF2
    cos = np.cos(np.arange(1,7) * T)
    sin = np.sin(np.arange(1,7) * T)
    x = np.empty(12)
    x[0::2] = sin
    x[1::2] = cos
    y = np.ones(13)
    y[1:] = x

    cf2 = AF2 * y

    #calculate the Fourier time series for M(3000)F2
    cos = np.cos(np.arange(1,5) * T)
    sin = np.sin(np.arange(1,5) * T)
    x = np.empty(8)
    x[0::2] = sin
    x[1::2] = cos
    y = np.ones(9)
    y[1:] = x

    cm3 = Am3 * y

    return (cf2, cm3)

def legendre_calculation(modip, latitude, longitude, CF2, Cm3):
    """

    :param modip:
    :param latitude:
    :param longitude:
    :param CF2:
    :param Cm3:
    :return: foF2 [MHz], M(3000)F2
    """

    assert np.size(CF2) == 12

    M = np.ones(12)
    P = np.ones(9)
    S = np.ones(9)
    C = np.ones(9)

    # Compute MODIP coefficients
    factor = np.sin(modip * np.pi/180)
    for i in range(1,12):
        M[i] = factor * M[i-1]

    # Compute latitude and longitude coefficients
    factor = np.cos(latitude * np.pi / 180)
    for i in range(1,9):
        P[i] = factor * P[i-1]
        S[i] = np.sin(i * longitude * np.pi / 180)
        C[i] = np.cos(i * longitude * np.pi / 180)

    # Compute foF2
    # Order 0 term
    foF2_1 = np.sum(CF2 * M)
    # Legendre grades
    Q = np.array([12,12,9,5,2,1,1,1,1])
    K = np.empty(9)
    K[0] = -Q[0]
    for i in range(1,9):
        K[i] = K[i-1] + 2 * Q[i-1] # [-12.,  12.,  36.,  54.,  64.,  68.,  70.,  72.,  74.]

    foF2_n = np.zeros(9)
    foF2_n[0] = foF2_1
    for i in range(1,9): # check! there is a bug in the indices
        for j in range(Q[i]):
            foF2_n[i] += (CF2[K[i] + 2 * j - 1] * C[i] ) * M[j] * P[i]
            foF2_n[i] += CF2[K[i] + 2 * j ] * S[i] * M[j] * P[i]

    foF2 = np.sum(foF2_n)

    # compute 0 order term
    M3000F2_1 = np.sum(Cm3[:7] * M[:7])

    R = np.array([7,8,6,3,2,1,1])
    H = np.empty(7)
    H[0] = -R[0]
    for i in range(1,7):
        H[i] = H[i-1] + 2 * R[i-1] # [ -7.,   7.,  23.,  35.,  41.,  45.,  47.])

    M3000F2_n = np.zeros(7)
    M3000F2_n[0] = M3000F2_1

    for i in range(1,7):
        for j in range(R[i]):
            M3000F2_n[i] += (Cm3[H[i] + 2*j - 1] * C[i] + Cm3[H[i] + 2*j ] * S[i] ) * (M[j]*P[j])

    M3000F2 = np.sum(M3000F2_n)

    NmF2 = 0.124 * foF2**2

    return foF2, M3000F2, NmF2

def get_hmE():
    return 120 # [km]

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
    """

    numerator = 1490 * M3000F2 * np.sqrt((0.0196 * M3000F2**2 + 1)/(1.2967 * M3000F2**2 -1))
    if foE<10**-30:
        deltaM = - 0.012
    else:
        r = float(foF2) / foE
        top = r * np.exp(20 * (r - 1.75)) + 1.75
        bottom = np.exp(20 * (r - 1.75))  + 1.0
        rho = top/bottom
        deltaM = 0.253/ (rho - 1.215) - 0.012

    denominator = M3000F2 + deltaM

    hmF2 = numerator/denominator - 176

    return hmF2

# THICKNESS PARAMETERS
def get_B2bot(NmF2, foF2, M3000F2):
    """

    :param NmF2:[10^11 m^-3 ]
    :param foF2: [MHz]
    :param M3000F2:
    :return: [km]
    """

    top = 0.385 * NmF2
    bottom = 0.01*np.exp(-3.467 + 0.857 * np.log(foF2**2) + 2.02 * np.log(M3000F2))

    return top/bottom

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
    return max(B1bot,7)

def get_BEbot():
    """

    :return: [km]
    """
    return 5.0

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
    if foF1  < 0.5:
        A2 = 0.0
        A3 = 4.0 * (NmE - Epst(A1, hmF2, B2bot, hmE))
    else:
        A3a = 4.0 * NmE
        for i in range(5):
            A2a = 4.0 * (NmF1 - Epst(A1, hmF2, B2bot, hmF1) - Epst(A3a, hmE, BEtop, hmF1))
            A2a = (A2a * np.exp(A2a - 0.8 * NmF1) + 0.8*NmF1) / (1 + np.exp(A2a - 0.8 * NmF1))
            A3a = 4.0 * (NmE - Epst(A2a, hmF1, B1bot, hmE) - Epst(A1, hmF2, B2bot, hmE))
        A2 = A2a
        A3 = (A3a * np.exp(60 * (A3a - 0.005)) + 0.05) / (1 + np.exp(60 * (A3a - 0.005)))

    return A2,A3

def Epst(peak_amp,peak_height,thickness,H):
    return peak_amp * np.exp((H - peak_height)/thickness) / (1 + np.exp((H - peak_height)/thickness)) ** 2

def shape_parameter(mth, NmF2, hmF2, B2bot, Azr):
    """

    :param mth:
    :param NmF2:
    :param hmF2:
    :param B2bot:
    :param Azr:
    :return:
    """
    if mth in [4,5,6,7,8,9]:
        ka = 6.705 - 0.014 * Azr  - 0.008 * hmF2
    elif mth in [1,2,3,10,11,12]:
        ka = -7.77 + 0.097 * (hmF2 / B2bot)**2 + 0.153 * NmF2
        kb  = (ka * np.exp(ka-2) + 2) / (1 + np.exp(ka - 2))
        k = (8 * np.exp(kb - 8) + kb) / (1 + np.exp(kb - 8))
    return k

def get_H0(B2bot, k):
    """

    :param B2bot: [km]
    :param k:
    :return: topside thickness parameter [km]
    """

    Ha = k * B2bot
    x = (Ha - 150.0)/ 100.0
    v = (0.041163 * x - 0.183981) * x + 1.424472
    Ho = Ha / v
    return Ho

def electron_density(h, lat, long, a0, a1, a2, mth, UT):
    """

    :param h:
    :param lat:
    :param long:
    :param a0:
    :param a1:
    :param a2:
    :param mth:
    :param UT:
    :return:
    """
    # TODO: implement chain of processing to get parameters
    hmF2 = get_hmF2(foE, foF2, M3000F2)

    if h < hmF2: #bottomside
        return bottomside_electrondensity(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot,B1top, B1bot, BEtop, BEbot)
    else:
        return topside_electrondensity(h, NmF2, hmF2, H0)


def bottomside_electrondensity(h, A1, A2, A3, hmF2, hmF1, hmE, B2bot,B1top, B1bot, BEtop, BEbot):
    """

    :param h:
    :param A1:
    :param A2:
    :param A3:
    :param hmF2:
    :param hmF1:
    :param hmE:
    :param B2bot:
    :param B1top:
    :param B1bot:
    :param BEtop:
    :param Bebot:
    :return:
    """
    if h > hmE:
        BE = BEtop
    else:
        BE = BEbot

    if h > hmF1:
        BF1 = B1top
    else:
        BF1 = B1bot

    #Compute the exponential arguments for each layer
    if h < 100:
        h = 100
    alpha1 = (h - hmF2) / B2bot
    alpha2 = (h - hmF1) / BF1 * np.exp(10 / (1 + abs(h - hmF2)))
    alpha3 = (h - hmE) / BE * np.exp(10 / (1 + abs(h - hmF2)))
    alpha = np.array([alpha1, alpha2, alpha3])

    A = np.array([A1,A2,A3])

    S = np.empty(3)
    mask1 = np.abs(alpha) > 25
    mask2 = np.logical_not(mask1)
    S[mask1] = 0
    S[mask2] = np.dot( A[mask2], np.exp(alpha[mask2])/np.power(1 + np.exp(alpha[mask2]), 2) )

    if h >= 100:
        N = (np.sum(S)) * 10**11
    else:
        B = np.array([B2bot, BF1, BE])
        ds = np.empty(3)
        ds[mask1] = 0
        ds[mask2] = (1 - np.exp(alpha[mask2])) / (B[mask2]* (1 + np.exp(alpha[mask2])) )

    # Chapman parameters
    BC = 1 - 10 * np.sum(S * ds) / np.sum(S)
    z = (h - 100) / 10.0

    # electron density
    N = np.sum(S) * np.exp(1 - BC * z - np.exp(-z)) * 10**11

    return N

def topside_electrondensity(h, NmF2, hmF2, H0):
    """

    :param h:
    :param NmF2:
    :param hmF2:
    :param H0:
    :return:
    """
    g = 0.125
    r = 100

    deltah = h - hmF2
    z = deltah / (H0 * (1 + r*g*deltah/ (r*H0 + g * deltah)) )

    ea = np.exp(z)
    if ea > 10**11:
        N = 4 * NmF2/ ea * 10**11
    else:
        N = 4 * NmF2 * ea * 10**11 / (1 + ea) ** 2

    return N

def vertical_integration(h1, h2, epsilon, A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot, NmF2, H0):
    """

    :param h1:
    :param h2:
    :param epsilon:
    :param A1:
    :param A2:
    :param A3:
    :param hmF2:
    :param hmF1:
    :param hmE:
    :param B2bot:
    :param B1top:
    :param B1bot:
    :param BEtop:
    :param BEbot:
    :param NmF2:
    :param H0:
    :return:
    """
    n = 8
    deltan = (h2 - h1) / n
    g = 0.5773502691896 *   deltan # why this strange number
    g1 =  0.5773502691896 * (h2 - h1) # my guess
    y = g1 + (deltan - g)/2

    GN2 = 0
    for i in range(n):
        GN2 += deltan/2.0 * electron_density(h, lat, )
