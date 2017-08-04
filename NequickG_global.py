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
        Nequick model is fully defined by time and solar activity.
        But numerical values can only be evaluated by constructing local NequickG model
        :param time: Nequick time object
        :param broadcast: Nequick broadcast object
        """
        self.time = time
        self.broadcast = broadcast

    def get_Nequick_local(self, position):
        Para = NequickG_parameters(position, self.broadcast, self.time)
        return NequickG(Para), Para

    def slant_table(self, n, ray, path, ex_attr=None):
        """writes a table of electron density along slant ray. optional paramenters allowed.
        Function inspired by International Reference Ionosphere model's interface"""

        # This is unpythonic and hence commented out
        # but at least it advertises the attributes available

        # allowed_attrs = set(['foF1', 'foF2', 'foE', 'M3000F2', 'NmF2', 'NmF1', 'NmE', 'hmE', 'hmF1', 'hmF2', 'modip',
        #      'Az','Azr', 'solarsine', 'solarcosine', 'chi', 'chi_eff', 'H0', 'B1bot', 'B1top', 'B2bot', 'BEtop', 'BEbot'
        #      ,'A1', 'A2', 'A3', 'k'])
        # if not set(ex_attr).issubset(allowed_attrs):
        #     raise ValueError('Invalid attribute present')

        with file(path, 'w') as f:
            writer = csv.writer(f, delimiter=',')
            rs, lats, lons, delta = ray.linspace(n)
            hs = ray.radius2height(rs)
            header = ['height', 'lat', 'lon', 'el_density'] + ex_attr
            writer.writerow(header)
            for i in range(len(lats)):
                # create local nequick object
                pos = Position(lats[i], lons[i])
                NEQ, para = self.get_Nequick_local(pos)

                # write values
                out = [hs[i], lats[i], lons[i]]
                out.append(NEQ.electrondensity(hs[i]))
                for attr in ex_attr:
                    out.append(getattr(para, attr))
                writer.writerow(out)

    def _gaussspace2(self, n, ray):
        """Segment a ray for Gauss quadrature by cartesian distance"""
        # same result as _segment(...)
        xx, deltax = gaussquadrature2_segment(n, ray.ob_x, ray.sat_x)
        yy, deltay = gaussquadrature2_segment(n, ray.ob_y, ray.sat_y)
        zz, deltaz = gaussquadrature2_segment(n, ray.ob_z, ray.sat_z)
        delta = np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)
        rr, latlat, lonlon = cartesian2coord(xx, yy, zz)
        hh = radius2height(rr)

        return hh, latlat, lonlon, delta

    def _gassspace(self, n, ray):
        """Segment a ray for Gauss quadrature by perigee distance"""
        # same result as _segment2(...)
        s1, s2 = ray.perigee_distancelimits()
        ss, delta = gaussquadrature2_segment(n, s1, s2)
        hs, lats, lons = ray.perigeedistance2coords(ss)

        return hs, lats, lons, delta

    def _integrate(self, hh, lats, lons, delta):
        electrondensity = np.empty(len(hh))

        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            NEQ, para = self.get_Nequick_local(pos)
            electrondensity[i] = NEQ.electrondensity(hh[i])
        GN = delta / 2.0 * np.sum(electrondensity)

        return GN

    def sTEC(self, ray, tolerance=None):
        """
        slant TEC integration.
        Doubles the number of intervals until change falls within tolerance
        :param ray:
        :param tolerance:
        :return:
        """
        # implementation can be optimised further

        # assert (ray.isvalid())
        seg = self._gassspace
        if tolerance == None:
            if ray.ob_h < 1000:
                tolerance = 0.001
            else:
                tolerance = 0.01

        if ray.p_radius < 0.1:
            print "calulcating vTEC instead"
            pos = Position(ray.ob_lat, ray.ob_lon)
            neq, para = self.get_Nequick_local(pos)
            return neq.vTEC(ray.ob_h, ray.sat_h)

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
            print "Warning: Integration2 did not converge"

        return (GN2 + (GN2 - GN1) / 15.0) * 1000

    def sTEC2(self, ray):
        """
        slant TEC integration with increased tolerance in topside
        implemented as directed by Galileo reference document to save computation
        not needed in this implementation because computational power is not scare
        """
        ht, latt, lont = ray.height2coords(1000)
        ray1 = Ray(ray.ob_h, ray.ob_lat, ray.ob_lon, ht, latt, lont)
        ray2 = Ray(ht, latt, lont, ray.sat_h, ray.sat_lat, ray.sat_lon)
        stec1 = self.sTEC(ray1, tolerance=0.001)
        stec2 = self.sTEC(ray2, tolerance=0.01)

        return stec1 + stec2

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

    def map_parameters(self, attrs, lat1, lon1, lat2, lon2, resolution=40):
        """
        Two points (lat1,lon1) (lat2,lon2) defines a rectangle in geographical grid
        attrs
        :param attrs: list of Nequick parameter names [strings] to map
        :param lat1:
        :param lon1:
        :param lat2:
        :param lon2:
        :param resolution:
        :return: at least three 2D grids. one of longitude, one for latitude.
        the rest are maps of paramenters
        """

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
                    try:
                        out[i, j] = getattr(para, attrs[k])
                    except AttributeError:
                        if attrs[k] == 'vTEC':
                            out[i, j] = neq.vTEC(100, 20000)
                        else:
                            raise AttributeError

        return latlat, lonlon, outs


class Ray:
    # TODO: make construction lazy
    def __init__(self, h1, lat1, lon1, h2, lat2, lon2):
        "By convention, point 1 is at a lower height"
        self.ob_h = h1
        self.ob_radius = 6371.2 + h1
        self.ob_lat = lat1
        self.ob_lon = lon1

        self.sat_h = h2
        self.sat_radius = 6371.2 + h2
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
        if self.p_radius > self.ob_radius and self.p_radius < self.sat_radius:
            r = self.p_radius
        else:
            r = self.ob_radius

        if r < 6371.2:
            print r
            print self.p_radius
            print self.ob_radius
            return False
        else:
            return True

    def arange(self, step):
        xs = np.arange(self.ob_x, self.sat_x, step)
        ys = np.arange(self.ob_y, self.sat_y, step)
        zs = np.arange(self.ob_z, self.sat_z, step)

        rs, lats, lons = cartesian2coord(xs, ys, zs)

        return rs, lats, lons

    def linspace(self, n):
        xs = np.linspace(self.ob_x, self.sat_x, n)
        ys = np.linspace(self.ob_y, self.sat_y, n)
        zs = np.linspace(self.ob_z, self.sat_z, n)

        rs, lats, lons = cartesian2coord(xs, ys, zs)

        delta = np.sqrt((xs[1] - xs[0]) ** 2 + (ys[1] - ys[0]) ** 2 + (ys[1] - ys[0]) ** 2)

        return rs, lats, lons, delta

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

    def height2coords(self, h):
        h = np.array(h)
        s = np.sqrt((6371.2 + h) ** 2 + self.p_radius ** 2)

        return self.perigeedistance2coords(s)

    @staticmethod
    def radius2height(r):
        return r - 6371.2


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


def segment2(n, ray):
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
