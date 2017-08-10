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
        r1 = 6371.2 + self.h1  # radius of earth
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
        assert not np.any(np.array(lon1) > 360)
        assert not np.any(np.array(lon2) > 360)
        assert not np.any(np.array(zenith) > 180)

        assert not np.any(np.array(lat1) < -90)
        assert not np.any(np.array(lat2) < -90)
        assert not np.any(np.array(lon1) < -360)
        assert not np.any(np.array(lon2) < -360)
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
        xx, deltax = gaussquadrature2_segment(n, self.x1, self.x2)

        yy, deltay = gaussquadrature2_segment(n, self.y1, self.y2)

        zz, deltaz = gaussquadrature2_segment(n, self.z1, self.z2)

        delta = np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)

        return xx, yy, zz, delta

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

        ray = Ray(h1, lat1, lon1, h2, lat2, lon2)

        if ray.p_radius < 0.1:
            neq, para = self.get_Nequick_local(Position(lat1, lon1))
            return neq.map_vTEC(h1, h2)

        s1, s2 = ray.perigee_distancelimits()

        n = 8
        ss, delta = gaussquadrature2_segment(n, s1, s2)
        hs, lats, lons = ray.perigeedistance2coords(ss)
        GN1 = self.__integrate2(hs, lats, lons, delta)
        n *= 2
        ss, delta = gaussquadrature2_segment(n, s1, s2)
        hs, lats, lons = ray.perigeedistance2coords(ss)
        GN2 = self.__integrate2(hs, lats, lons, delta) # there is repeated work here. can be optimized

        count = 1
        while (abs(GN2 - GN1) > tolerance * abs(GN1)) and count < 20:
            GN1 = GN2
            n *= 2
            ss, delta = gaussquadrature2_segment(n, s1, s2)
            hs, lats, lons = ray.perigeedistance2coords(ss)
            GN2 = self.__integrate2(hs, lats, lons, delta)
            count += 1

        if count == 20:
            print "Integration did not converge"

        return (GN2 + (GN2 - GN1) / 15.0) * 1000

    def __integrate(self, s1, s2, n, ray):
        ss, delta = gaussquadrature2_segment(n, s1, s2)

        positions = []
        NEQs = []
        electrondensity = []
        heights, lats, lons = ray.perigeedistance2coords(ss)
        for i in range(len(lats)):
            pos = Position(lats[i], lons[i])
            positions.append(pos)
            NEQ, para = self.get_Nequick_local(pos)
            NEQs.append(NEQ)
            electrondensity.append(NEQ.electrondensity(heights[i]))
        GN = delta / 2.0 * np.sum(electrondensity)
        return GN
