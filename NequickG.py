import csv
import numpy as np
import matplotlib.pyplot as plt
from aux import *


class NEQTime:
    def __init__(self, mth, universal_time):
        self.mth = mth  # {01, 02, 03 ... 11, 12}
        self.universal_time = universal_time  # hours and decimals


class Position:
    def __init__(self, latitude, longitude):
        self.latitude = latitude  # degrees
        self.longitude = longitude  # degrees


class GalileoBroadcast:
    def __init__(self, ai0, ai1, ai2):
        self.ai0 = ai0
        self.ai1 = ai1
        self.ai2 = ai2


class NequickG_global:
    def __init__(self, time, broadcast):
        """

        :param time: Nequick time object
        :param broadcast: Nequick broadcast object
        """
        self.time = time
        self.broadcast = broadcast

    def get_Nequick_local(self, position):
        """
        :param position: list of nequick position objects
        :return: list of nequick models at positions
        """
        if type(position) == list:
            models = []
            for pos in position:
                Para = NequickG_parameters(pos, self.broadcast, self.time)
                models.append(NequickG(Para))

            return models
        else:
            Para = NequickG_parameters(position, self.broadcast, self.time)
            return NequickG(Para)

    def sTEC2(self, h1, lat1, lon1, h2, lat2, lon2, tolerance=None):
        "this method avoids the 0.1km threshold in ray perigee radius"
        if tolerance == None:
            if h1 < 1000:
                tolerance = 0.001
            else:
                tolerance = 0.01
        roughpaper = SlantRayAnalysis2(h1, lat1, lon1, h2, lat2, lon2)

        n = 8

        xx, yy, zz, delta = roughpaper.sampleX(n)
        rr, latlat, lonlon = roughpaper.cartesian2coord(xx, yy, zz)
        hh = roughpaper.radius2height(rr)

        GN1 = self.__integrate2(hh, latlat, lonlon, delta)

        n *= 2

        xx, yy, zz, delta = roughpaper.sampleX(n)
        rr, latlat, lonlon = roughpaper.cartesian2coord(xx, yy, zz)
        hh = roughpaper.radius2height(rr)

        GN2 = self.__integrate2(hh, latlat, lonlon, delta)  # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 5:

            GN1 = GN2

            n *= 2

            xx, yy, zz, delta = roughpaper.sampleX(n)
            rr, latlat, lonlon = roughpaper.cartesian2coord(xx, yy, zz)
            hh = roughpaper.radius2height(rr)

            GN2 = self.__integrate2(hh, latlat, lonlon, delta)

            count += 1

        return (GN2 + (GN2 - GN1) / 15.0)

    def __integrate2(self, hh, lats, lons, delta):
        NEQs = []
        electrondensity = []

        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            NEQ = self.get_Nequick_local(pos)
            NEQs.append(NEQ)
            electrondensity.append(NEQ.electrondensity(hh[i]))
        GN = delta / 2.0 * np.sum(electrondensity)

        return GN

    def sTEC(self, h1, lat1, lon1, h2, lat2, lon2, tolerance=None):
        """

        :param h1:
        :param lat1:
        :param lon1:
        :param h2:
        :param lat2:
        :param lon2:
        :param tolerance:
        :return:
        """

        if tolerance == None:
            if h1 < 1000:
                tolerance = 0.001
            else:
                tolerance = 0.01

        roughpaper = SlantRayAnalysis(h1, lat1, lon1, h2, lat2, lon2)

        if roughpaper.rp < 0.1:
            neq = self.get_Nequick_local(Position(lat1, lon1))
            return neq.vTEC(h1, h2)

        s1, s2 = roughpaper.ray_endpoints()

        n = 8

        GN1 = self.__integrate(s1, s2, n, roughpaper)
        n *= 2
        GN2 = self.__integrate(s1, s2, n, roughpaper)  # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 5:
            GN1 = GN2
            n *= 2
            GN2 = self.__integrate(s1, s2, n, roughpaper)
            count += 1

        return (GN2 + (GN2 - GN1) / 15.0)


    def __sampleX(self, x1, x2, n):
        """
        Decides which 2n points to sample for integration between x1 and x2
        :param x1:
        :param x2:
        :param n:
        :return:
        """
        delta = float(x2 - x1) / n
        g = .5773502691896 * delta  # delta / sqrt(3)
        y = x1 + (delta - g) / 2.0

        ss = np.empty(2 * n)
        I = np.arange(n)
        ss[0::2] = y + I * delta
        ss[1::2] = y + I * delta + g

        return ss, delta

    def __integrate(self, s1, s2, n, roughpaper):
        ss, delta = self.__sampleX(s1, s2, n)

        positions = []
        NEQs = []
        electrondensity = []
        heights, lats, lons = roughpaper.ray_coords(ss)
        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            positions.append(pos)
            NEQ = self.get_Nequick_local(pos)
            NEQs.append(NEQ)
            electrondensity.append(NEQ.electrondensity(heights[i]))
        GN = delta / 2.0 * np.sum(electrondensity)
        return GN


class NequickG:
    def __init__(self, parameters):
        self.Para = parameters
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

        mask1 = h < self.Para.hmF2
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
        GN2 = self.__single_quad(h1, h2, n)  # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 5:
            GN1 = GN2
            n *= 2
            GN2 = self.__single_quad(h1, h2, n)
            count += 1

        return (GN2 + (GN2 - GN1) / 15.0)

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


class SlantRayAnalysis:
    """
    This class automatically computes all geometric quanties related to slant ray.
    Think of it as an intern who does all the tedious calculation

    Intermediate computation uses Ray Perigee as construct
    Reference: 2.5.8.2
    """

    def __init__(self, h1, lat1, lon1, h2, lat2, lon2):

        self.h1 = h1
        self.h2 = h2
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2

        self.zenith = self.__zenith_angle(h1, lat1, lon1, h2, lat2, lon2)
        self.rp = self.__ray_perigee_radius(h1, h2, self.zenith)
        self.latp, self.lonp = self.__ray_perigee_coord(lat1, lon1, lat2, lon2, self.zenith)
        self.great_angle = self.__great_circle_angle(lat2, lon2, self.latp, self.lonp)
        self.sigma_p = self.__ray_perigee_azimuth(lat2, lon2, self.latp, self.lonp, self.great_angle)

    def rayperigee_position(self):
        return self.rp, self.latp, self.lonp

    def ray_endpoints(self):
        r1 = 6371.2 + self.h1
        r2 = 6371.2 + self.h2

        s1 = np.sqrt(r1 ** 2 - self.rp ** 2)
        s2 = np.sqrt(r2 ** 2 - self.rp ** 2)

        return s1, s2

    def ray_coords(self, s):

        hs = np.sqrt(np.power(s, 2) + self.rp ** 2) - 6371.2

        rp, latp, lonp, sigma_p = self.rp, self.latp, self.lonp, self.sigma_p

        DR = np.pi / 180
        sine_sigma_p, cosine_sigma_p = (np.sin(sigma_p * DR), np.cos(sigma_p * DR))
        # great circle parameters
        # perigee triangle angle
        tan_delta = np.divide(s, rp)
        cosine_delta = 1.0 / np.sqrt(1 + tan_delta ** 2)

        sine_delta = tan_delta * cosine_delta
        # lat
        sin_lats = np.sin(latp * DR) * cosine_delta + np.cos(latp * DR) * sine_delta * cosine_sigma_p
        cos_lats = np.sqrt(1 - sin_lats ** 2)
        lats = np.arctan2(sin_lats, cos_lats) * 180 / np.pi

        # lon
        sin_lons = sine_delta * sine_sigma_p * np.cos(latp * DR)
        cos_lons = cosine_delta - np.sin(latp * DR) * sin_lats
        lons = np.arctan2(sin_lons, cos_lons) * 180 / np.pi + lonp

        return hs, lats, lons

    def __greatcircle_delta(self, lat1, lon1, lat2, lon2):
        DR = np.pi / 180.0

        cosine_delta = np.sin(lat1 * DR) * np.sin(lat2 * DR) + np.cos(lat1 * DR) * np.cos(lat2 * DR) * np.cos(
            (lon2 - lon1) * DR)
        sine_delta = np.sqrt(1 - cosine_delta ** 2)

        return np.arctan2(sine_delta, cosine_delta) * 180 / np.pi

    def __zenith_angle(self, h1, lat1, lon1, h2, lat2, lon2):
        """

        :param lat1: [deg]
        :param lon1: [deg]
        :param lat2: [deg]
        :param lon2: [deg]
        :return: [deg]
        """
        delta = self.__greatcircle_delta(lat1, lon1, lat2, lon2)

        # cosine = np.sin(lat1) * np.sin(lat2) +np.cos(lat1) * np.cos(lat2) *np.cos(lon2 - lon1)
        # sine = np.sqrt(1 - cosine**2)
        DR = np.pi / 180
        cosine = np.cos(delta * DR)
        sine = np.sin(delta * DR)

        r1 = 6371.2 + h1
        r2 = 6371.2 + h2

        zenith = np.arctan2(sine, cosine - r1 / r2) * 180.0 / np.pi

        return zenith

    # Ray-perigee computation
    def __ray_perigee_radius(self, h1, h2, zenith):

        r1 = 6371.2 + h1
        r2 = 6371.2 + h2

        DR = np.pi / 180.0

        rp = r1 * np.sin(zenith * DR)

        return rp

    def __ray_perigee_coord(self, lat1, lon1, lat2, lon2, zenith):

        assert not np.any(np.array(lat1) > 90)
        assert not np.any(np.array(lat2) > 90)
        assert not np.any(np.array(lon1) > 90)
        assert not np.any(np.array(lon2) > 90)
        assert not np.any(np.array(zenith) > 180)

        assert not np.any(np.array(lat1) < -90)
        assert not np.any(np.array(lat2) < -90)
        assert not np.any(np.array(lon1) < -90)
        assert not np.any(np.array(lon2) < -90)
        assert not np.any(np.array(zenith) < -180)

        DR = np.pi / 180.0
        delta = self.__greatcircle_delta(lat1, lon1, lat2, lon2)
        cosine_delta = np.cos(delta * DR)
        sine_delta = np.sin(delta * DR)

        sine_sigma = np.sin((lon2 - lon1) * DR) * np.cos(lat2 * DR) / sine_delta
        cosine_sigma = (np.sin(lat2 * DR) - cosine_delta * np.sin(lat1 * DR)) / (sine_delta * np.cos(lat1 * DR))

        assert not np.any(np.array(sine_sigma) > 1)
        assert not np.any(np.array(cosine_sigma) > 1)

        delta_p = 90 - zenith

        sine_latp = np.sin(lat1 * DR) * np.cos(delta_p * DR) - np.cos(lat1 * DR) * np.sin(delta_p * DR) * cosine_sigma
        cosine_latp = np.sqrt(1 - sine_latp ** 2)

        assert not np.any(np.array(sine_latp) > 1)
        assert not np.any(np.array(cosine_latp) > 1)

        if abs(abs(lat1) - 90) < 10 ** -10 * 180 / np.pi:
            if lat1 < 0:
                latp = -zenith
            else:
                latp = zenith

        else:
            latp = np.arctan2(sine_latp, cosine_latp) * 180 / np.pi

        if abs(abs(lat1) - 90) < 10 ** -10:
            if zenith < 0:
                lonp = lon2
            else:
                lonp = lon2 + 180
        else:
            sine_lonp = - sine_sigma * np.sin(delta_p * DR) / np.cos(latp * DR)
            cosine_lonp = (np.cos(delta_p * DR) - np.sin(lat1 * DR) * np.sin(latp * DR)) / (
            np.cos(lat1 * DR) * np.cos(latp * DR))
            lonp = np.arctan2(sine_lonp, cosine_lonp) * 180 / np.pi + lon1

        return latp, lonp

    def __great_circle_angle(self, lat2, lon2, latp, lonp):
        DR = np.pi / 180
        if abs(abs(latp) - 90) < 10 ** 10 * 180 / np.pi:
            return abs(lat2 - latp)
        else:
            cosine = np.sin(latp * DR) * np.sin(lat2 * DR) + np.cos(latp * DR) * np.cos(lat2 * DR) * np.cos(
                (lon2 - lonp) * DR)
            sine = np.sqrt(1 - cosine ** 2)

            return np.arctan2(sine, cosine) * 180 / np.pi

    def __ray_perigee_azimuth(self, lat2, lon2, latp, lonp, great_angle):

        if abs(abs(latp) - 90) < 10 ** -10 * 180 / np.pi:
            if latp < 0:
                return 0
            else:
                return 180
        else:
            DR = np.pi / 180.0

            sine_sigmap = np.sin((lon2 - lonp) * DR) * np.cos(lat2 * DR) / np.sin(great_angle * DR)
            cosine_sigmap = (np.sin(lat2 * DR) - np.cos(great_angle * DR) * np.sin(latp)) / (
            np.sin(great_angle * DR) * np.cos(latp * DR))
            return np.arctan2(sine_sigmap, cosine_sigmap) * 180 / np.pi

class SlantRayAnalysis2:
    """
    Similar to SlantRayAnalysis. Except that intermediate computation is done in cartesian which is much simpler
    """
    def __init__(self, h1, lat1, lon1, h2, lat2, lon2):

        self.h1 = h1
        self.h2 = h2
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2

        self.x1, self.y1, self.z1 = self.coord2cartesian(6371.2 + h1, lat1, lon1)
        self.x2, self.y2, self.z2 = self.coord2cartesian(6371.2 + h2, lat2, lon2)

    def coord2cartesian(self, r, lat, lon):
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

    def cartesian2coord(self, x, y, z):
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

    def radius2height(self, r):
        return r - 6371.2

    def height2radius(self, h):
        return h + 6371.2

    def intermediatepoints(self, n):
        xx = np.linspace(self.x1, self.x2, n)
        yy = np.linspace(self.y1, self.y2, n)
        zz = np.linspace(self.z1, self.z2, n)

        rr, latlat, lonlon = self.cartesian2coord(xx, yy, zz)

        return rr, latlat, lonlon

    def sampleX(self, n):
        """
        Decides which 2n points to sample for integration between x1 and x2
        :param n:
        :return:
        """
        deltax = float(self.x2 - self.x1) / n
        x_g = .5773502691896 * deltax  # delta / sqrt(3)
        x_y = self.x1 + (deltax - x_g) / 2.0

        xx = np.empty(2 * n)
        I = np.arange(n)
        xx[0::2] = x_y + I * deltax
        xx[1::2] = x_y + I * deltax + x_g

        deltay = float(self.y2 - self.y1) / n
        y_g = .5773502691896 * deltay  # delta / sqrt(3)
        y_y = self.y1 + (deltax - y_g) / 2.0

        yy = np.empty(2 * n)
        I = np.arange(n)
        yy[0::2] = y_y + I * deltay
        yy[1::2] = y_y + I * deltay + y_g

        deltaz = float(self.z2 - self.z1) / n
        z_g = .5773502691896 * deltaz  # delta / sqrt(3)
        z_y = self.z1 + (deltaz - y_g) / 2.0

        zz = np.empty(2 * n)
        I = np.arange(n)
        zz[0::2] = z_y + I * deltaz
        zz[1::2] = z_y + I * deltaz + z_g

        delta = np.sqrt(deltax**2 + deltay**2 + deltaz**2)


        return xx, yy, zz, delta


class NequickG_parameters:
    def __init__(self, pos, broadcast, time):
        self.Position = pos  # Nequick position object
        self.Broadcast = broadcast  # Nequick broadcast object
        self.Time = time  # Nequick time object
        self.stmodip_path = '/home/tpl/Documents/Airbus/Project/Papers/Nequick/CCIR_MoDIP/modipNeQG_wrapped.txt'
        self.CCIR_path = '/home/tpl/Documents/Airbus/Project/Papers/Nequick/CCIR_MoDIP/ccir'
        self.compute_parameters()

    def compute_parameters(self):
        # Stage 1
        self.__read_stMoDIP__()
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
    def __read_stMoDIP__(self):
        with open(self.stmodip_path) as f:
            data = list(rec for rec in csv.reader(f, delimiter=','))
            data = [map(float, row) for row in data]
            return data

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
        stModip = self.__read_stMoDIP__()

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
        self.modip = mu
        return mu

    def __effective_ionization__(self):
        """
        :param ai0:
        :param ai1:
        :param ai2:
        :param MoDIP: [deg]
        :return: Effective ionoiszation level
        Remark: Az is equivalent to F10.7 in climatological NeQuick

        """
        ai0 = self.Broadcast.ai0
        ai1 = self.Broadcast.ai1
        ai2 = self.Broadcast.ai2

        MoDIP = self.modip

        if (ai0 == 0) and (ai1 == 0) and (ai2 == 0):
            Az = 63.7
        else:
            Az = ai0 + ai1 * MoDIP + ai2 * MoDIP ** 2

        # Reference Section 3.3
        assert Az <= 400
        assert Az >= 0

        self.Az = Az
        return Az

    def __effective_sunspot_number__(self):
        """
        This parameter is equivalent to R 12 in climatological NeQuick
        Reference: Section 2.5.4.5
        :param Az:
        :return:

        """
        self.Azr = np.sqrt(167273 + (self.Az - 63.7) * 1123.6) - 408.99
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
        t = dy + (18 - universal_time) / 24
        # Compute the argument
        a_m = (0.9856 * t - 3.289) * (np.pi / 180)  # radians
        a_l = a_m + (1.916 * np.sin(a_m) + 0.020 * np.cos(2 * a_m) + 282.634) * (np.pi / 180)  # radians

        # Compute sine and cosine of solar declination
        Sine = 0.39782 * np.sin(a_l)
        Cosine = np.sqrt(1 - Sine ** 2)

        self.solarsine = Sine
        self.solarcosine = Cosine

        return (Cosine, Sine)

    def __localtime(self, universal_time, longitude):
        return universal_time + longitude / 15.0

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
        self.chi = np.arctan2(np.sqrt(1 - coschi ** 2), coschi)

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
        seasp = seas * (ee - 1) / (ee + 1)

        self.foE = np.sqrt((1.112 - 0.019 * seasp) ** 2 * np.sqrt(Az) * np.cos(chi_eff * np.pi / 180) ** 0.6 + 0.49)
        self.NmE = NeqCriticalFreqToNe(self.foE)

        return self.foE, self.NmE

    def __readccirXXfiles__(self):
        """
        Reference: Section 2.5.3.2
        :param month:
        :param path:
        :return:f2_ijk array and fm3_ijk
        """
        path = self.CCIR_path
        month = self.Time.mth
        data = []
        with open(path + str(month + 10) + '.txt') as f:
            for row in csv.reader(f, delimiter=' '):
                row = [num for num in row if num != '']  # filter
                data = data + [float(num) for num in row]
        assert (len(data) == 2858)

        F2data = data[:1976]
        self.F2 = np.reshape(np.array(F2data), (2, 76, 13))

        Fm3data = data[1976:]
        self.Fm3 = np.reshape(np.array(Fm3data), (2, 49, 9))

        return self.F2, self.Fm3

    def __interpolate_AZR__(self):
        """

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
        modip = self.modip
        latitude = self.Position.latitude
        longitude = self.Position.longitude
        CF2 = self.CF2
        Cm3 = self.Cm3

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
        K = np.empty(9, dtype=np.int)
        K[0] = -Q[0]
        for i in range(1, 9):
            K[i] = K[i - 1] + 2 * Q[i - 1]  # [-12,  12,  36,  54,  64,  68,  70,  72,  74]

        foF2_n = np.zeros(9)
        foF2_n[0] = foF2_1
        for i in range(1, 9):  # check! there might be a bug in the indices
            for j in range(Q[i]):
                foF2_n[i] += CF2[K[i] + 2 * j] * C[i] * M[j] * P[i]
                foF2_n[i] += CF2[K[i] + 2 * j + 1] * S[i] * M[j] * P[i]
        self.foF2 = np.sum(foF2_n)

        # compute 0 order term
        M3000F2_1 = np.sum(Cm3[:7] * M[:7])

        R = np.array([7, 8, 6, 3, 2, 1, 1])
        H = np.empty(7, dtype=np.int)
        H[0] = -R[0]
        for i in range(1, 7):
            H[i] = H[i - 1] + 2 * R[i - 1]  # [ -7,   7,  23,  35,  41,  45,  47])

        M3000F2_n = np.zeros(7)
        M3000F2_n[0] = M3000F2_1

        for i in range(1, 7):
            for j in range(R[i]):
                M3000F2_n[i] += (Cm3[H[i] + 2 * j - 1] * C[i] + Cm3[H[i] + 2 * j] * S[i]) * (M[j] * P[j])

        self.M3000F2 = np.sum(M3000F2_n)

        self.NmF2 = NeqCriticalFreqToNe(self.foF2)

        return self.foF2, self.M3000F2, self.NmF2

    def F1Layer(self):
        """

        :param foE: E layer critical frequency [MHz]
        :param foF2: F2 layer critical frequency [MHz]
        :return: (foF1, NmF1 [10^11 m^-3] )
        """
        foE = self.foE
        foF2 = self.foF2

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

        self.foF1 = foF1
        self.NmF1 = NmF1

        return foF1, NmF1

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
        if foE < 10 ** -30:
            deltaM = - 0.012
        else:
            r = float(foF2) / foE
            top = r * np.exp(20 * (r - 1.75)) + 1.75
            bottom = np.exp(20 * (r - 1.75)) + 1.0
            rho = top / bottom
            deltaM = 0.253 / (rho - 1.215) - 0.012

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
            A2 = 0.0
            A3 = 4.0 * (self.NmE - epstein(self.A1, self.hmF2, self.B2bot, self.hmE))
        else:
            A3a = 4.0 * self.NmE
            for i in range(5):
                A2a = 4.0 * (
                    self.NmF1 - epstein(self.A1, self.hmF2, self.B2bot, self.hmF1) - epstein(A3a, self.hmE, self.BEtop,
                                                                                             self.hmF1))
                A2a = NeqJoin(A2a, 0.8 * self.NmF1, 1, A2a - 0.8 * self.NmF1)
                A3a = 4.0 * (
                    self.NmE - epstein(A2a, self.hmF1, self.B1bot, self.hmE) - epstein(self.A1, self.hmF2, self.B2bot,
                                                                                       self.hmE))
            self.A2 = A2a
            self.A3 = NeqJoin(A3a, 0.05, 60.0, A3a - 0.005)
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
        kb = (ka * np.exp(ka - 2) + 2) / (1 + np.exp(ka - 2))
        self.k = (8 * np.exp(kb - 8) + kb) / (1 + np.exp(kb - 8))
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

        self.AmpF1 = A2
        self.AmpF2 = A1
        self.AmpE = A3

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

        mask2 = np.logical_not(mask1)  # chapman corrections needed

        dsF2 = (1 - np.exp(alphaF2)) / (thickF2 * (1 + np.exp(alphaF2)))
        dsF2[np.abs(alphaF2) > 25] = 0
        dsF1 = (1 - np.exp(alphaF1)) / (thickF1 * (1 + np.exp(alphaF1)))
        dsF1[np.abs(alphaF1) > 25] = 0
        dsE = (1 - np.exp(alphaE)) / (thickE * (1 + np.exp(alphaE)))
        dsE[np.abs(alphaE) > 25] = 0

        BC = 1 - 10 * (EpstF2 * dsF2 + EpstF1 * dsF1 + EpstE * dsE) / S
        z = (h - 100) / 10.0

        N[mask2] = S[mask2] * np.exp(1 - BC[mask2] * z[mask2] - np.exp(-z[mask2])) * 10 ** 11

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

        return N

########################################################################################
