from aux import *
import CCIR_MoDIP.ccir_fm3
import CCIR_MoDIP.ccir_f2
import CCIR_MoDIP.modip
import numpy as np

# np.seterr(invalid = 'raise')

class NEQTime:
    def __init__(self, mth, universal_time):
        self.mth = int(mth)  # {01, 02, 03 ... 11, 12}
        self.universal_time = universal_time  # hours and decimals

    def __hash__(self):
        return hash((self.mth, self.universal_time))

    def __eq__(self, other):
        return (self.mth, self.universal_time) == (other.mth, other.universal_time)


class Position:
    def __init__(self, latitude, longitude):
        self.latitude = latitude  # degrees
        self.longitude = longitude  # degrees

    def __hash__(self):
        return hash((int(self.latitude), int(self.longitude)/2))

    def __eq__(self, other):
        # arbitrary design decision: decimal places in position have no significance
        return (int(self.latitude), int(self.longitude)/2) == (int(other.latitude), int(other.longitude)/2)


class GalileoBroadcast:
    def __init__(self, ai0, ai1, ai2):
        self.ai0 = ai0
        self.ai1 = ai1
        self.ai2 = ai2




class NequickG:
    def __init__(self, parameters):
        # self.Para = parameters
        self.hmF2 = parameters.hmF2
        topside_para = parameters.topside_para()
        bottomside_para = parameters.bottomside_para()
        self.topside = NequickG_topside(*topside_para)
        self.bottomside = NequickG_bottomside(*bottomside_para)

    def electrondensity(self, h):
        """

        :param h: [km]
        :return: electron density [m^-3]
        """
        h = np.array(h)

        mask1 = h < self.hmF2
        mask2 = np.logical_not(mask1)

        h_bot = h[mask1]
        h_top = h[mask2]

        N = np.empty(np.shape(h))

        N[mask1] = self.bottomside.electrondensity(h_bot)
        N[mask2] = self.topside.electrondensity(h_top)

        assert (not np.any(N < 0))

        return N

    def vTEC(self, h1, h2, tolerance=None):
        """
        Vertical TEC numerical Integration
        :param h1: integration lower endpoint
        :param h2: integration higher endpoint
        :param tolerance:
        :return:
        """

        assert (h2 > h1)

        if tolerance == None:
            if h1 < 1000:
                tolerance = 0.001
            else:
                tolerance = 0.01

        n = 8

        GN1 = self.__single_quad(h1, h2, n)
        n *= 2
        GN2 = self.__single_quad(h1, h2, n)  # TODO: there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 20:
            GN1 = GN2
            n *= 2
            GN2 = self.__single_quad(h1, h2, n)
            count += 1

        if count == 20:
            print "vTEC integration did not converge"

        return (GN2 + (GN2 - GN1) / 15.0)


    def vTEC_ratio(self):
        bot = self.vTEC(0, self.hmF2)
        top = self.vTEC(self.hmF2, 20000)

        return top / bot

    def __single_quad(self, h1, h2, n):

        delta = float(h2 - h1) / n

        g = .5773502691896 * delta  # delta / sqrt(3)
        y = h1 + (delta - g) / 2.0

        h = np.empty(2 * n)
        I = np.arange(n)
        h[0::2] = y + I * delta
        h[1::2] = y + I * delta + g
        N = self.electrondensity(h)
        GN = delta / 2.0 * np.sum(N)

        return GN




class NequickG_parameters:
    def __init__(self, pos, broadcast, time):
        self.Position = pos  # Nequick position object
        self.Broadcast = broadcast  # Nequick broadcast object
        self.Time = time  # Nequick time object
        self.compute_parameters()

    def compute_parameters(self):
        # Stage 1
        self.__compute_MODIP__()
        self.__effective_ionization__()
        self.__effective_sunspot_number__()
        self.__solar_declination__()
        self.__solar_zenith__()
        self.__effective_solar_zenith__()

        # Stage 2
        self.__readccirXXfiles__()
        self.__interpolate_AZR__()
        self.__F2fouriertimeseries__()
        self.F2Layer()
        self.ELayer()
        self.F1Layer()

        # Stage 3
        self.get_hmE()
        self.get_hmF2()
        self.get_hmF1()

        # Stage 4
        self.get_B2bot()
        self.get_B1top()
        self.get_B1bot()
        self.get_BEtop()
        self.get_BEbot()

        # Stage 5
        self.get_A1()
        self.get_A2A3()
        self.shape_parameter()
        self.get_H0()

    def topside_para(self):
        return [self.NmF2, self.hmF2, self.H0]

    def bottomside_para(self):
        # ensure that the order agrees with the argument order of bottomside class
        return [self.hmE, self.hmF1, self.hmF2, self.BEtop, self.BEbot, self.B1top, self.B1bot, self.B2bot, self.A1,
                self.A2, self.A3]

    ############################ STAGE 1####################################
    def __compute_MODIP__(self):
        """
        Reference: Section 2.5.4.3 of GSA's "Ionospheric Correction
        Algorithm for Galileo Single Frequency Users"

        :param latitude: [deg]
        :param longitude: [deg]
        :param stModip: array of MODIP vales
        :return: MODIP [deg]
        """

        latitude = self.Position.latitude
        longitude = self.Position.longitude

        sq , lon_excess, lat_excess= self.__square(latitude, longitude)
        mu2 = self.__interpolate2d(sq, lon_excess, lat_excess)
        self.modip = mu2

        return self.modip

    def __square(self, latitude, longitude):
        """
        extract a 4x4 array of "nearby" modip values from the standard Modip table
        :param latitude:
        :param longitude:
        :return:
        --x--x--x--x--
        --------------
        --------------
        --x--x--x--x--
        ------*-------
        --------------
        --x--x--x--x--
        --------------
        --------------
        --x--x--x--x--
        """
        # what if longitude is [0,360]
        # what if longitude is [-180, 180]
        num_division = 36
        lon_division = 10.0
        lat_division = 5.0

        # Longitude
        lon = (longitude + 180)/lon_division
        lon_start = int(lon) - 2 # range : -2 to 34 or 16 to 52 or -20 to 16
        lon_excess = lon - int(lon)

        # max permissible lon_start is 34
        # min permissible lon_start is 0
        # so this hack is needed
        if (lon_start < 0):
            lon_start = lon_start + num_division
        if (lon_start > (num_division - 3)):
            lon_start = lon_start - num_division

        lat = (latitude + 90.0) / lat_division + 1
        lat_start = int(lat - 1e-6) - 2 # why?
        lat_excess = lat - lat_start - 2

        stModip = np.array(CCIR_MoDIP.modip.stModip)
        square =stModip[lat_start+1:lat_start+5,lon_start+2:lon_start+6]
        return square, lon_excess, lat_excess

    def __interpolate2d(self, Z, x, y):
        """

        :param Z: 4x4 2D array of modip values
        :param x: normalised longitude displacement [0<x<1]
        :param y: normalised latitude displacement [0<y<1]
        :return: lagrange 3-order interpolated modip
        """
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


    def __effective_ionization__(self):
        """

        :param ai0:
        :param ai1:
        :param ai2:
        :param MoDIP: [deg]
        :return: Effective ionoiszation level
        Remark: Az is equivalent to F10.7 in climatological NeQuick
        It is a mystery to me why a solar index is a designed to vary with geographical location
        """
        ai0 = self.Broadcast.ai0
        ai1 = self.Broadcast.ai1
        ai2 = self.Broadcast.ai2

        MoDIP = self.modip

        if (ai0 == 0) and (ai1 == 0) and (ai2 == 0):
            Az = 63.7
        else:
            Az = ai0 + ai1 * MoDIP + ai2 * MoDIP ** 2

        if Az < 0:
            Az = 0
        elif Az > 400:
            Az = 400
        # Reference Section 3.3
        assert Az <= 400
        assert Az >= 0
        self.Az = Az
        return self.Az

    def __effective_sunspot_number__(self):
        """
        This parameter is equivalent to R 12 in climatological NeQuick
        Reference: Section 2.5.4.5
        :param Az:
        :return:

        """
        self.Azr = np.sqrt(167273 + (self.Az - 63.7) * 1123.6) - 408.99
        # min = -99.636351241819057

        # CCIR recommends that Azr saturates at 150
        # does a negative Azr even make sense??
        # assert(not np.isnan(self.Azr))
        # self.Azr_unclip = self.Azr
        # if self.Azr < 0:
        #     self.Azr = 0
        # assert(self.Azr > 0)

        # assert(self.Azr < 150)
        # if self.Azr > 150:
        #     self.Azr = 150

        return self.Azr

    def __solar_declination__(self, ):
        """
        Compute sin(delta_Sun ), cos(delta_Sun ), the sine and cosine of the solar declination.
        :param month: [mth]
        :param universal_time: [hours and decimals]
        :return:(Cosine, Sine)
        """

        month = self.Time.mth
        universal_time = self.Time.universal_time

        # Compute day of year at the middle of the month
        dy = 30.5 * month - 15
        # Compute time [days]:
        t = dy + (18 - universal_time) / 24.0
        # Compute the argument
        a_m = (0.9856 * t - 3.289) * (np.pi / 180)  # radians
        a_l = a_m + (1.916 * np.sin(a_m) + 0.020 * np.sin(2 * a_m) + 282.634) * (np.pi / 180)  # radians

        # Compute sine and cosine of solar declination
        Sine = 0.39782 * np.sin(a_l)
        Cosine = np.sqrt(1 - Sine ** 2)

        self.solarsine = Sine
        self.solarcosine = Cosine

        return (Cosine, Sine)

    def __localtime(self, universal_time, longitude):
        return universal_time + longitude / 15.0 # 15 degrees per hour

    def __solar_zenith__(self):
        """
        Reference : Section 2.5.4.7
        :param latitude: [deg]
        :param localtime: [hours and decimals]
        :param solarcosine: of solar declination
        :param solarsine: of solar declination
        :return: solar zenith angle [deg]
        """

        latitude = self.Position.latitude
        LT = self.__localtime(self.Time.universal_time, self.Position.longitude)
        solarsine = self.solarsine
        solarcosine = self.solarcosine

        coschi = np.sin(latitude * np.pi / 180) * solarsine + np.cos(
            latitude * np.pi / 180) * self.solarcosine * np.cos(
            np.pi / 12 * (12 - LT))

        self.chi = np.arctan2(np.sqrt(1 - coschi ** 2), coschi) * 180 / np.pi

        return self.chi

    def __effective_solar_zenith__(self):
        """
        Reference: Section 2.5.4.8
        :param chi: solar zenith angle [deg]
        :param chi0: default solar zenith angle [deg]
        :return: effective solar zenith angle[deg]
        """

        chi0 = 86.23292796211615
        self.chi_eff = NeqJoin(90 - 0.24 * NeqClipExp(20 - 0.2 * self.chi), self.chi, 12, self.chi - chi0)
        return self.chi_eff

    def get_stage1para(self):
        name = ['modip', 'effective_ionization', 'effective_sunspot_number', 'solarcosine', 'solarsine', 'solar_zenith', 'effective_solar_zenith']
        return name, [self.modip, self.Az, self.Azr, self.solarcosine, self.solarsine, self.chi, self.chi_eff]
    ############################ STAGE 2####################################
    def ELayer(self):
        """
        Reference: 2.5.5.1
        :param latitude: [deg]
        :param Az: Effective Ionisation Level
        :param chi_eff: effective_solar_zenith [deg]
        :param month: mth
        :return: E layer critical frequency foE [MHz], NmE [10^11 m^-3]
        """
        latitude = self.Position.latitude
        Az = self.Az
        chi_eff = self.chi_eff
        month = self.Time.mth

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
        self.seasp = seas * (ee - 1) / (ee + 1)
        seasp = self.seasp
        self.foE = np.sqrt((1.112 - 0.019 * seasp) ** 2 * np.sqrt(Az) * np.cos(chi_eff * np.pi / 180) ** 0.6 + 0.49)
        self.NmE = NeqCriticalFreqToNe(self.foE)


        # #checking
        # chin = NeqJoin(90 - 0.24*NeqClipExp(20-0.2*self.chi), self.chi, 12, self.chi - 86.23292796211615)
        # sfac = (1.112 - 0.019 * seasp) * np.sqrt(np.sqrt(self.Az))
        # fa = sfac * NeqClipExp(np.log(np.cos(chin*np.pi/180)) * 0.3)
        # foE = np.sqrt(fa * fa + 0.49)
        #
        # print foE - self.foE
        if (self.foE < 0):
            print self.foE
        assert(self.foE >= 0)
        return self.foE, self.NmE

    def __readccirXXfiles__(self):
        """
        Reference: Section 2.5.3.2
        :param month:
        :param path:
        :return:f2_ijk array and fm3_ijk
        """

        self.F2 = np.array(CCIR_MoDIP.ccir_f2.F2[self.Time.mth])
        self.Fm3 = np.array(CCIR_MoDIP.ccir_fm3.Fm3[self.Time.mth])

        return self.F2, self.Fm3

    def __interpolate_AZR__(self):
        """
        AZR is a type of solar index
        CCIR provides "spherical harmonic" coefficients of foF2 and M3000F2 at high and low AZR
        Linearly interpolate the coefficients based on current AZR
        :param F2:
        :param Fm3:
        :param Azr:
        :return:
        """
        F2 = self.F2
        Fm3 = self.Fm3
        Azr = self.Azr

        self.AF2 = F2[0, :, :] * (1 - Azr / 100.0) + F2[1, :, :] * Azr / 100.0
        self.Am3 = Fm3[0, :, :] * (1 - Azr / 100.0) + Fm3[1, :, :] * Azr / 100.0

        return self.AF2, self.Am3

    def __F2fouriertimeseries__(self):
        """

        :param UT: Universal Time UT [hours],
        :param AF2:
        :param Am3:
        :return:CF2, Cm3 vectors of coefficients for Legendre calculation for foF2 and M(3000)F2
        """
        UT = self.Time.universal_time
        AF2 = self.AF2
        Am3 = self.Am3

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

        self.CF2 = np.sum(AF2 * y, axis=1)

        # calculate the Fourier time series for M(3000)F2
        cos = np.cos(np.arange(1, 5) * T)
        sin = np.sin(np.arange(1, 5) * T)
        x = np.empty(8)
        x[0::2] = sin
        x[1::2] = cos
        y = np.ones(9)
        y[1:] = x

        self.Cm3 = np.sum(Am3 * y, axis=1)

        return self.CF2, self.Cm3

    def __geographical_variation__(self, Q):
        """

        :param Q: array of length of modip expansion
        :return:
        Refer to page 18 of "Advances in Ionospheric Mapping by Numerical Methods (1966)"
        """
        modip = self.modip
        latitude = self.Position.latitude
        longitude = self.Position.longitude

        G = []
        m = np.sin(modip * np.pi / 180)
        p = np.cos(latitude * np.pi / 180)
        for i in range(Q[0]):
            G.append(m**i)

        for i in range(1,len(Q)):
            for j in range(Q[i]):
                G.append(m**j * p**i * np.cos(i * longitude * np.pi / 180) )
                G.append(m**j * p**i * np.sin(i * longitude * np.pi / 180) )

        return np.array(G)


    def F2Layer(self):
        """
        legendre_calculation
        :param modip:
        :param latitude:
        :param longitude:
        :param CF2:
        :param Cm3:
        :return: foF2 [MHz], M(3000)F2
        """
        CF2 = self.CF2
        Cm3 = self.Cm3

        G = self.__geographical_variation__([12, 12, 9, 5, 2, 1, 1, 1, 1])
        assert (len(G) == 76)
        foF2 = np.sum(CF2 * G)
        self.foF2 = foF2


        G = self.__geographical_variation__([7, 8, 6, 3, 2, 1, 1])
        assert (len(G) == 49)
        M3000F2 = np.sum(Cm3 * G)
        self.M3000F2 = M3000F2
        assert ( self.M3000F2 > 0 )

        self.NmF2 = NeqCriticalFreqToNe(self.foF2)

        # if foF2 < 0:
        #     print foF2, self.Azr, self.Az, self.NmF2, self.Position.latitude, self.Position.longitude
        # assert (foF2 > 0) # this will fail

        return self.foF2, self.M3000F2, self.NmF2

    def F1Layer(self):
        """

        :param foE: E layer critical frequency [MHz]
        :param foF2: F2 layer critical frequency [MHz]
        :return: (foF1, NmF1 [10^11 m^-3] )
        """
        foE = self.foE
        foF2 = self.foF2

        # In the day, foF1 = 1.4foE. In the night time, foF1 = 0
        # use NeqJoin for a smooth day-night transition
        # gradient factor of 1000 is arbitrary and large so that neqjoin can approx a step function
        # why is foE = 2 a threshold for day -night boundary?

        # print '1.4 * foE', 1.4 * foE
        foF1 = NeqJoin(1.4 * foE, 0, 1000.0, foE - 2)
        # print 'foF1: ', foF1
        foF1 = NeqJoin(0, foF1, 1000.0, foE - foF1)
        # print 'foF1: ', foF1
        foF1 = NeqJoin(foF1, 0.85 * foF1, 60.0, 0.85 * foF2 - foF1)
        # print 'foF1: ', foF1

        if foF1 < 10 ** -6:
            foF1 = 0

        # F1 layer maximum density
        if (foF1 <= 0) and (foE > 2): # how can foF1 be negative??
            NmF1 = NeqCriticalFreqToNe(foE + 0.5)
        else:
            NmF1 = NeqCriticalFreqToNe(foF1)

        assert (foF1 >= 0)
        assert (NmF1 >= 0)
        self.foF1 = foF1
        self.NmF1 = NmF1

        return foF1, NmF1

    def get_stage2para(self):
        name = ['foE', 'foF1', 'foF2', 'M3000F2']
        return name , [self.foE, self.foF1, self.foF2, self.M3000F2]
    ############################ STAGE 3####################################
    def get_hmE(self):
        self.hmE = 120
        return 120  # [km]

    def get_hmF1(self):
        """

        :param hmF2: [km]
        :param hmE: [km]
        :return:
        """
        self.hmF1 = (self.hmF2 + self.hmE) / 2.0

        return self.hmF1

    def get_hmF2(self):
        """

        :param foE: [Mhz]
        :param foF2: [Mhz]
        :param M3000F2:
        :return: [km]
        Based on Dudeney (1983)
        """
        foE = self.foE
        foF2 = self.foF2
        M3000F2 = self.M3000F2
        numerator = 1490 * M3000F2 * np.sqrt((0.0196 * M3000F2 ** 2 + 1) / (1.2967 * M3000F2 ** 2 - 1))
        assert (not np.any(np.isnan(numerator)))
        if foE < 10 ** -30: # avoid divide by zero
            deltaM = - 0.012
        else:
            r = float(foF2) / foE
            top = r * np.exp(20 * (r - 1.75)) + 1.75
            bottom = np.exp(20 * (r - 1.75)) + 1.0
            rho = top / bottom
            deltaM = 0.253 / (rho - 1.215) - 0.012
            assert rho > 1.73
        denominator = M3000F2 + deltaM

        self.hmF2 = numerator / denominator - 176

        return self.hmF2

    ############################ STAGE 4####################################

    def get_B2bot(self):
        """

        :param NmF2:[10^11 m^-3 ]
        :param foF2: [MHz]
        :param M3000F2:
        :return: [km]
        """

        top = 0.385 * self.NmF2
        bottom = 0.01 * np.exp(-3.467 + 0.857 * np.log(self.foF2 ** 2) + 2.02 * np.log(self.M3000F2))

        self.B2bot = top / bottom
        return self.B2bot

    def get_B1top(self):
        """

        :param hmF1: [km]
        :param hmF2: [km]
        :return: [km]
        """
        self.B1top = 0.3 * (self.hmF2 - self.hmF1)
        return self.B1top

    def get_B1bot(self):
        """

        :param hmF1: [km]
        :param hmE: [km]
        :return: [km]
        """
        self.B1bot = 0.5 * (self.hmF1 - self.hmE)
        return self.B1bot

    def get_BEtop(self):
        """

        :param B1bot: [km]
        :return: [km]
        """
        self.BEtop = max(self.B1bot, 7)
        return self.BEtop

    def get_BEbot(self):
        """

        :return: [km]
        """
        self.BEbot = 5.0
        return 5.0

    ############################ STAGE 5####################################
    def get_A1(self):
        """

        :param NmF2:  [10^11 m^-3 ]
        :return: F2 layer amplitude
        """
        self.A1 = 4 * self.NmF2
        return self.A1

    def get_A2A3(self):
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
        if self.foF1 < 0.5:
            self.A2 = 0.0
            self.A3 = 4.0 * (self.NmE - epstein(self.A1, self.hmF2, self.B2bot, self.hmE))
            return self.A2, self.A3
        else:
            A3a = 4.0 * self.NmE
            # print 'A3a: ', A3a
            for i in range(5):
                A2a = 4.0 * (self.NmF1 -
                             epstein(self.A1, self.hmF2, self.B2bot, self.hmF1) -
                             epstein(A3a, self.hmE, self.BEtop, self.hmF1))
                # print 'A2a: ', A2a
                A2a = NeqJoin(A2a, 0.8 * self.NmF1, 1, A2a - 0.8 * self.NmF1)
                # print 'A2a', A2a
                A3a = 4.0 * (self.NmE -
                             epstein(A2a, self.hmF1, self.B1bot, self.hmE) -
                             epstein(self.A1, self.hmF2, self.B2bot, self.hmE))

            self.A2 = A2a
            self.A3 = NeqJoin(A3a, 0.05, 60.0, A3a - 0.005)
            # print 'A3: ', self.A3
        return self.A2, self.A3

    def shape_parameter(self):
        """

        :param mth:
        :param NmF2:
        :param hmF2:
        :param B2bot:
        :param Azr:
        :return:
        """
        mth = self.Time.mth
        if mth in [4, 5, 6, 7, 8, 9]:
            ka = 6.705 - 0.014 * self.Azr - 0.008 * self.hmF2
        elif mth in [1, 2, 3, 10, 11, 12]:
            ka = -7.77 + 0.097 * (self.hmF2 / self.B2bot) ** 2 + 0.153 * self.NmF2
        else:
            raise ValueError("Invalid Month")
        # kb = (ka * np.exp(ka - 2) + 2) / (1 + np.exp(ka - 2))
        # self.k = (8 * np.exp(kb - 8) + kb) / (1 + np.exp(kb - 8))
        kb = NeqJoin(ka,2,1,ka - 2.0)
        kb = NeqJoin(8, kb, 1, kb - 8)

        self.k = kb
        return self.k

    def get_H0(self):
        """

        :param B2bot: [km]
        :param k:
        :return: topside thickness parameter [km]
        """

        Ha = self.k * self.B2bot
        x = (Ha - 150.0) / 100.0
        v = (0.041163 * x - 0.183981) * x + 1.424472
        self.H0 = Ha / v
        return self.H0


class NequickG_bottomside:
    def __init__(self, hmE, hmF1, hmF2, BEtop, BEbot, B1top, B1bot, B2bot, A1, A2, A3):
        self.hmF2 = hmF2
        self.hmF1 = hmF1
        self.hmE = hmE


        self.BEtop = BEtop
        self.BEbot = BEbot
        self.B1top = B1top
        self.B1bot = B1bot
        self.B2bot = B2bot

        self.AmpF2 = A1
        self.AmpF1 = A2
        self.AmpE = A3

        assert (hmF2 > 0)
        assert (hmF1 > 0)
        assert (hmE > 0)
        assert (BEtop > 0)
        assert (BEbot > 0)
        assert (B1top > 0)
        assert (B1bot > 0)
        assert (B2bot > 0)
        assert (A1 >= 0)
        # assert (A2 >= 0)
        # assert (A3 >= 0) # this will fail. Why? Becuse E layer does not exist at night

    def electrondensity(self, h):
        assert not np.any(h > self.hmF2)

        if type(h) != np.ndarray:
            h = np.array(h)

        BE = np.empty(np.shape(h))
        mask = h > self.hmE
        BE[mask] = self.BEtop
        BE[np.logical_not(mask)] = self.BEbot

        BF1 = np.empty(np.shape(h))
        mask = h > self.hmF1
        BF1[mask] = self.B1top
        BF1[np.logical_not(mask)] = self.B1bot

        # Compute the exponential arguments for each layer

        h[h < 100] = 100

        # thickness parameter with fade out exponent for E and F1 layer
        thickF2 = self.B2bot
        thickF1 = BF1 / np.exp(10 / (1 + np.abs(h - self.hmF2)))
        thickE = BE / np.exp(10 / (1 + np.abs(h - self.hmF2)))

        EpstF2 = epstein(self.AmpF2, self.hmF2, thickF2, h)
        EpstF1 = epstein(self.AmpF1, self.hmF1, thickF1, h)
        EpstE = epstein(self.AmpE, self.hmE, thickE, h)

        # suppress small values in epstein layer
        diffF2 = (h - self.hmF2)
        diffF1 = (h - self.hmF1)
        diffE = (h - self.hmE)

        alphaF2 = diffF2 / thickF2
        alphaF1 = diffF1 / thickF1
        alphaE = diffE / thickE

        EpstF2[np.abs(alphaF2) > 25] = 0
        EpstF1[np.abs(alphaF1) > 25] = 0
        EpstE[np.abs(alphaE) > 25] = 0

        # sum the 3 semi eptstein layers
        N = np.empty(np.shape(h))

        mask1 = h >= 100  # no corrections needed
        S = EpstF2 + EpstF1 + EpstE
        N[mask1] = S[mask1] * 10 ** 11

        mask2 = np.logical_not(mask1)  # chapman corrections needed for heights less than 100km

        dsF2 = (1 - np.exp(alphaF2)) / (thickF2 * (1 + np.exp(alphaF2)))
        dsF2[np.abs(alphaF2) > 25] = 0
        assert( not np.any(np.isnan(dsF2)))
        dsF1 = (1 - np.exp(alphaF1)) / (thickF1 * (1 + np.exp(alphaF1)))
        dsF1[np.abs(alphaF1) > 25] = 0
        assert( not np.any(np.isnan(dsF1)))
        dsE = (1 - np.exp(alphaE)) / (thickE * (1 + np.exp(alphaE)))
        dsE[np.abs(alphaE) > 25] = 0
        assert( not np.any(np.isnan(dsE)))
        BC = 1 - 10 * (EpstF2 * dsF2 + EpstF1 * dsF1 + EpstE * dsE) / S
        z = (h - 100) / 10.0

        N[mask2] = S[mask2] * np.exp(1 - BC[mask2] * z[mask2] - np.exp(-z[mask2])) * 10 ** 11
        assert( np.all(N > 0) )
        return N


class NequickG_topside:
    def __init__(self, NmF2, hmF2, H0):
        self.hmF2 = hmF2
        self.NmF2 = NmF2
        self.H0 = H0

    def electrondensity(self, h):
        assert not np.any(h <= self.hmF2)
        g = 0.125
        r = 100

        deltah = h - self.hmF2
        z = deltah / (self.H0 * (1 + r * g * deltah / (r * self.H0 + g * deltah)))
        ea = np.exp(z)

        if type(h) == np.ndarray:
            mask1 = ea > 10 ** 11
            mask2 = np.logical_not(mask1)
            N = np.empty(np.size(h))
            N[mask1] = 4 * self.NmF2 / ea[mask1] * 10 ** 11
            N[mask2] = 4 * self.NmF2 * ea[mask2] * 10 ** 11 / (1 + ea[mask2]) ** 2

        else:
            if ea > 10 ** 11:
                N = 4 * self.NmF2 / ea * 10 ** 11
            else:
                N = 4 * self.NmF2 * ea * 10 ** 11 / (1 + ea) ** 2

        assert( np.all(N > 0) )
        return N

########################################################################################
