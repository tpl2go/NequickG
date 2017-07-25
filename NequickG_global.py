import csv
import numpy as np
import matplotlib.pyplot as plt
from aux import *
import time
from NequickG import Position, GalileoBroadcast, NEQTime, NequickG, NequickG_parameters
import spheretrig


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
            return NequickG(Para), Para

    def ray(self, h1, lat1, lon1, h2, lat2, lon2):
        with open('Slant_Ray_Logger.dat', 'w') as f:
            writer = csv.writer(f, delimiter=',')

            roughpaper = SlantRayAnalysis2(h1, lat1, lon1, h2, lat2, lon2)

            n = 200

            xx, yy, zz, delta = roughpaper.sampleX(n)
            rr, latlat, lonlon = roughpaper.cartesian2coord(xx, yy, zz)
            hh = roughpaper.radius2height(rr)
            writer.writerow(['lon', 'lat', 'height', 'el density'])
            for i in range(len(latlat)):
                pos = Position(latlat[i], lonlon[i])
                NEQ, para = self.get_Nequick_local(pos)
                writer.writerow([lonlon[i], latlat[i], hh[i], NEQ.electrondensity(hh[i])])

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
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 20:
            GN1 = GN2

            n *= 2

            xx, yy, zz, delta = roughpaper.sampleX(n)
            rr, latlat, lonlon = roughpaper.cartesian2coord(xx, yy, zz)
            hh = roughpaper.radius2height(rr)

            GN2 = self.__integrate2(hh, latlat, lonlon, delta)

            count += 1

        if count == 20:
            print "Integration2 did not converge"

        return (GN2 + (GN2 - GN1) / 15.0) * 1000

    def __integrate2(self, hh, lats, lons, delta):
        NEQs = []
        electrondensity = []

        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            NEQ, para = self.get_Nequick_local(pos)
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

        roughpaper = Ray(h1, lat1, lon1, h2, lat2, lon2)

        if roughpaper.p_radius < 0.1:
            neq, para = self.get_Nequick_local(Position(lat1, lon1))
            return neq.vTEC(h1, h2)

        s1, s2 = roughpaper.ray_endpoints()

        n = 8

        GN1 = self.__integrate(s1, s2, n, roughpaper)
        n *= 2
        GN2 = self.__integrate(s1, s2, n, roughpaper)  # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 20:
            GN1 = GN2
            n *= 2
            GN2 = self.__integrate(s1, s2, n, roughpaper)
            count += 1

        if count == 20:
            print "Integration did not converge"

        return (GN2 + (GN2 - GN1) / 15.0) * 1000

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
            NEQ, para = self.get_Nequick_local(pos)
            NEQs.append(NEQ)
            electrondensity.append(NEQ.electrondensity(heights[i]))
        GN = delta / 2.0 * np.sum(electrondensity)
        return GN

    def map(self, lat1, lon1, lat2, lon2, resolution=40):
        lats = np.linspace(lat1, lat2, resolution)
        lons = np.linspace(lon1, lon2, resolution)

        lonlon, latlat = np.meshgrid(lons, lats)

        vtec = np.empty([resolution, resolution])

        for i in range(resolution):
            for j in range(resolution):
                lat = lats[i]
                lon = lons[j]

                pos = Position(lat, lon)
                neq, para = self.get_Nequick_local(pos)
                vtec[i, j] = neq.vTEC(100, 1000)
                # is this [i j] or [j i]
        return latlat, lonlon, vtec



class Ray:
    def __init__(self, h1, lat1, lon1, h2, lat2, lon2):
        "By convention, point 1 is at a lower height"
        self.ob_h = h1
        self.ob_lat = lat1
        self.ob_lon = lon1

        self.sat_h = h2
        self.sat_lat = lat2
        self.sat_lon = lon2

        self.ob_x, self.ob_y, self.ob_z = self.coord2cartesian(6371.2 + h1, lat1, lon1)
        self.sat_x, self.sat_y, self.sat_z = self.coord2cartesian(6371.2 + h2, lat2, lon2)

        self.length = np.sqrt(
            (self.sat_x - self.ob_x) ** 2 + (self.sat_y - self.ob_y) ** 2 + (self.sat_z - self.ob_z) ** 2)

        self.ob_zenith = spheretrig.zenith(h1, lat1, lon1, h2, lat2, lon2)
        self.ob_azimuth = spheretrig.azimuth(lat1, lon1, lat2, lon2)
        self.greatcircle = spheretrig.greatcircle(lat1, lon1, lat2, lon2)

        self.p_radius = self.perigee_radius()
        self.p_lat, self.p_lon = self.perigee_coords()
        self.p_azimuth = self.perigee_azimuth()


    def isvalid(self):
        # TODO: implement a method to test if ray passes through the earth
        pass

    def coord2cartesian(self, r, lat, lon):
        """

        :param r: [km]
        :param lat: [deg]
        :param lon: [deg]
        :return:
        """
        DR = np.pi / 180
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
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        xy = np.sqrt(x ** 2 + y ** 2)

        lat = np.arctan(z / xy) * 180 / np.pi
        lon = np.arctan2(y, x) * 180 / np.pi

        return r, lat, lon

    def radius2height(self, r):
        return r - 6371.2

    def height2radius(self, h):
        return h + 6371.2

    def linspace(self, n):
        xs = np.linspace(self.ob_x, self.sat_x, n)
        ys = np.linspace(self.ob_y, self.sat_y, n)
        zs = np.linspace(self.ob_z, self.sat_z, n)

        rs, lats, lons = self.cartesian2coord(xs, ys, zs)

        return rs, lats, lons

    # Ray-perigee computation
    def perigee_radius(self):

        r1 = 6371.2 + self.ob_h
        DR = np.pi / 180.0
        rp = r1 * np.sin(self.ob_zenith * DR)

        return rp

    def perigee2ob_greatcircle(self):
        return 90 - self.ob_zenith

    def perigee_coords(self):

        sigma = self.ob_azimuth
        delta_p = self.perigee2ob_greatcircle()
        zeta = self.ob_zenith

        DR = np.pi / 180

        # Perigee Latitude
        # spherical cosine law on spherical triangle between pole, observer and perigee
        sine_latp = np.sin(self.ob_lat * DR) * np.cos(delta_p * DR) - np.cos(self.ob_lat * DR) * np.sin(
            delta_p * DR) * np.cos(sigma * DR)
        cosine_latp = np.sqrt(1 - sine_latp ** 2)

        # if ray pierces pole
        if abs(abs(self.ob_lat) - 90) < 10 ** -10 * 180 / np.pi:
            if self.ob_lat < 0:
                latp = -zeta
            else:
                latp = zeta

        else:
            latp = np.arctan2(sine_latp, cosine_latp) * 180 / np.pi

        # Perigee Longitude
        sine_lonp = - np.sin(sigma * DR) * np.sin(delta_p * DR) / np.cos(latp * DR)
        cosine_lonp = (np.cos(delta_p * DR) - np.sin(self.ob_lat * DR) * np.sin(latp * DR)) / (
            np.cos(self.ob_lat * DR) * np.cos(latp * DR))
        # if ray pierces pole
        if abs(abs(self.ob_lat) - 90) < 10 ** -10 * 180 / np.pi:
            if zeta < 0:
                lonp = self.sat_lon
            else:
                lonp = self.sat_lon + 180
        else:
            lonp = np.arctan2(sine_lonp, cosine_lonp) * 180 / np.pi + self.ob_lon

        return latp, lonp

    def perigee2sat_greatcircle(self):
        latp = self.p_lat
        lonp = self.p_lon

        DR = np.pi / 180
        if abs(abs(latp) - 90) < 10 ** 10 * 180 / np.pi:
            return abs(self.sat_lat - latp)
        else:
            cosine = np.sin(latp * DR) * np.sin(self.sat_lat * DR) + np.cos(latp * DR) * np.cos(
                self.sat_lat * DR) * np.cos((self.sat_lon - lonp) * DR)
            sine = np.sqrt(1 - cosine ** 2)

            return np.arctan2(sine, cosine) * 180 / np.pi

    def perigee_azimuth(self):

        if abs(abs(self.p_lat) - 90) < 10 ** -10 * 180 / np.pi:
            if self.p_lat < 0:
                return 0
            else:
                return 180
        else:
            return spheretrig.azimuth(self.p_lat, self.p_lon, self.sat_lat, self.sat_lon)

    def perigee_distance(self):
        r1 = 6371.2 + self.ob_h  # radius of earth
        r2 = 6371.2 + self.sat_h

        s1 = np.sqrt(r1 ** 2 - self.p_radius ** 2)
        s2 = np.sqrt(r2 ** 2 - self.p_radius ** 2)

        return s1, s2

    def perigeedistance2coords(self, s):

        hs = np.sqrt(np.power(s, 2) + self.p_radius ** 2) - 6371.2

        rp, latp, lonp, sigma_p = self.p_radius, self.p_lat, self.p_lon, self.p_azimuth

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

        self.length = np.sqrt((self.x2 - self.x1) ** 2 + (self.y2 - self.y1) ** 2 + (self.z2 - self.z1) ** 2)

    def coord2cartesian(self, r, lat, lon):
        """

        :param r: [km]
        :param lat: [deg]
        :param lon: [deg]
        :return:
        """
        DR = np.pi / 180
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
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        xy = np.sqrt(x ** 2 + y ** 2)

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
        Decides which 2n points to sample for Gauss-Legendre integration between x1 and x2
        :param n:
        :return:
        """
        xx, deltax = gaussquadrature2_sampler(n, self.x1, self.x2)

        yy, deltay = gaussquadrature2_sampler(n, self.y1, self.y2)

        zz, deltaz = gaussquadrature2_sampler(n, self.z1, self.z2)

        delta = np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)

        return xx, yy, zz, delta



def gaussquadrature2_sampler(n, x1, x2):
    """returns array of x points to sample for second order Gauss Lagrange quadrature"""
    deltax = float(x2 - x1) / n
    x_g = .5773502691896 * deltax  # delta / sqrt(3)
    x_y = x1 + (deltax - x_g) / 2.0

    xx = np.empty(2 * n)
    I = np.arange(n)
    xx[0::2] = x_y + I * deltax
    xx[1::2] = x_y + I * deltax + x_g

    return xx, deltax
