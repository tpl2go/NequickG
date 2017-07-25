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
        # """
        # :param position: list of nequick position objects
        # :return: list of nequick models at positions
        # """
        # if type(position) == list:
        #     models = []
        #     for pos in position:
        #         Para = NequickG_parameters(pos, self.broadcast, self.time)
        #         models.append(NequickG(Para))
        #
        #     return models
        # else:
        #     Para = NequickG_parameters(position, self.broadcast, self.time)
        #     return NequickG(Para), Para
        Para = NequickG_parameters(position, self.broadcast, self.time)
        return NequickG(Para), Para

    def table(self, h1, lat1, lon1, h2, lat2, lon2):
        with open('Slant_Ray_Logger.dat', 'w') as f:
            writer = csv.writer(f, delimiter=',')

            ray= Ray(h1, lat1, lon1, h2, lat2, lon2)
            n = 200

            hs, lats, lons, delta = self._segment2(n, ray)
            writer.writerow(['lon', 'lat', 'height', 'foE', 'foF1', 'foF2', 'eldensity'])
            for i in range(len(lats)):
                pos = Position(lats[i], lons[i])
                NEQ, para = self.get_Nequick_local(pos)
                out = [lons[i], lats[i], hs[i]]
                out = out + [para.foE, para.foF1, para.foF2]
                out.append( NEQ.electrondensity(hs[i]))

                writer.writerow(out)

    def _segment2(self, n, ray):
        xx, deltax = gaussquadrature2_segment(n, ray.ob_x, ray.sat_x)
        yy, deltay = gaussquadrature2_segment(n, ray.ob_y, ray.sat_y)
        zz, deltaz = gaussquadrature2_segment(n, ray.ob_z, ray.sat_z)
        delta = np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)
        rr, latlat, lonlon = cartesian2coord(xx, yy, zz)
        hh = radius2height(rr)

        return hh, latlat, lonlon, delta

    def _segment(self, n, ray):
        s1, s2 = ray.perigee_distancelimits()
        ss, delta = gaussquadrature2_segment(n, s1, s2)
        hs, lats, lons = ray.perigeedistance2coords(ss)

        return hs, lats, lons, delta

    def sTEC(self, h1, lat1, lon1, h2, lat2, lon2, tolerance=None):
        seg = self._segment
        if tolerance == None:
            if h1 < 1000:
                tolerance = 0.001
            else:
                tolerance = 0.01

        ray = Ray(h1, lat1, lon1, h2, lat2, lon2)
        if ray.p_radius < 0.1:
            print "calulcating vTEC instead"
            pos = Position(lat1, lon1)
            neq, para = self.get_Nequick_local(pos)
            return neq.vTEC(h1, h2)

        n = 8
        hs, lats, lons, delta = seg(n, ray)
        GN1 = self._integrate(hs, lats, lons, delta)

        n *= 2

        hs, lats, lons, delta = seg(n, ray)
        GN2 = self._integrate(hs, lats, lons, delta)  # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 20:
            GN1 = GN2

            n *= 2

            hs, lats, lons, delta = seg(n, ray)
            GN2 = self._integrate(hs, lats, lons, delta)

            count += 1

        if count == 20:
            print "Integration2 did not converge"


        return (GN2 + (GN2 - GN1) / 15.0) * 1000

    def _integrate(self, hh, lats, lons, delta):
        electrondensity = np.empty(len(hh))

        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            NEQ, para = self.get_Nequick_local(pos)
            electrondensity[i] = NEQ.electrondensity(hh[i])
        GN = delta / 2.0 * np.sum(electrondensity)

        return GN

    def map_vTEC(self, lat1, lon1, lat2, lon2, resolution=40):
        lats = np.linspace(lat1, lat2, resolution)
        lons = np.linspace(lon1, lon2, resolution)

        lonlon, latlat = np.meshgrid(lons, lats)

        vtec = np.empty([resolution, resolution])

        for i in range(resolution):
            for j in range(resolution):
                pos = Position(latlat[i, j], lonlon[i, j])
                neq, para = self.get_Nequick_local(pos)
                vtec[i, j] = neq.vTEC(100, 20000)

        return latlat, lonlon, vtec

    def map_parameters(self,attrs, lat1, lon1, lat2, lon2, resolution=40):

        lats = np.linspace(lat1, lat2, resolution)
        lons = np.linspace(lon1, lon2, resolution)

        lonlon, latlat = np.meshgrid(lons, lats)

        outs = []
        for i in range(len(attrs)):
            outs.append(np.empty([resolution, resolution]))
        for i in range(resolution):
            for j in range(resolution):
                pos = Position(latlat[i, j], lonlon[i, j])
                neq, para = self.get_Nequick_local(pos)
                for k in range(len(attrs)):
                    out = outs[k]
                    out[i, j] = getattr(para, attrs[k])

        return latlat, lonlon, outs


class Ray:
    def __init__(self, h1, lat1, lon1, h2, lat2, lon2):
        "By convention, point 1 is at a lower height"
        self.ob_h = h1
        self.ob_lat = lat1
        self.ob_lon = lon1

        self.sat_h = h2
        self.sat_lat = lat2
        self.sat_lon = lon2

        self.ob_x, self.ob_y, self.ob_z = coord2cartesian(6371.2 + h1, lat1, lon1)
        self.sat_x, self.sat_y, self.sat_z = coord2cartesian(6371.2 + h2, lat2, lon2)

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

    def linspace(self, n):
        xs = np.linspace(self.ob_x, self.sat_x, n)
        ys = np.linspace(self.ob_y, self.sat_y, n)
        zs = np.linspace(self.ob_z, self.sat_z, n)

        rs, lats, lons = cartesian2coord(xs, ys, zs)

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

    def perigee_distancelimits(self):
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

def gaussquadrature2_segment(n, x1, x2):
    """returns array of x points to sample for second order Gauss Lagrange quadrature"""
    deltax = float(x2 - x1) / n
    x_g = .5773502691896 * deltax  # delta / sqrt(3)
    x_y = x1 + (deltax - x_g) / 2.0

    xx = np.empty(2 * n)
    I = np.arange(n)
    xx[0::2] = x_y + I * deltax
    xx[1::2] = x_y + I * deltax + x_g

    return xx, deltax

def coord2cartesian(r, lat, lon):
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

def cartesian2coord(x, y, z):
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

def radius2height(r):
    return r - 6371.2

def height2radius(h):
        return h + 6371.2

def segment2( n, ray):
    xx, deltax = gaussquadrature2_segment(n, ray.ob_x, ray.sat_x)
    yy, deltay = gaussquadrature2_segment(n, ray.ob_y, ray.sat_y)
    zz, deltaz = gaussquadrature2_segment(n, ray.ob_z, ray.sat_z)
    delta = np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)
    rr, latlat, lonlon = cartesian2coord(xx, yy, zz)
    hh = radius2height(rr)

    return hh, latlat, lonlon, delta

def segment(n, ray):
    s1, s2 = ray.perigee_distancelimits()
    ss, delta = gaussquadrature2_segment(n, s1, s2)
    hs, lats, lons = ray.perigeedistance2coords(ss)

    return hs, lats, lons, delta

if __name__ == "__main__":
    ray = Ray(0.07811, 82.49, 297.66, 20281.54618, 54.29, 8.23)

    hs, lats, lons, delta = segment(8, ray)
    print hs
    print lats
    print lons
    print delta

    hs, lats, lons, delta = segment2(8, ray)
    print hs
    print lats
    print lons
    print delta
