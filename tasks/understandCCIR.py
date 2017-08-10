"""
Context:
Each CCIR data file contains the spherical harmonic (I think) coefficients for the
atlas of foF2 and M(3000)F2 under two conditions:
1) low solar activity (sunspot number R12 (or Azr) = 0
2) high solar activity (sunspot number R12 (or Azr) = 100

Purpose:
This script visualises on a map what the CCIR data file encodes.
Relevant code is copied from main script to avoid unnecessary computation

"""
import numpy as np
import CCIR_MoDIP.ccir_fm3
import CCIR_MoDIP.ccir_f2
import CCIR_MoDIP.modip
from NequickG import Position, NEQTime
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


class CCIR_Unpacker:
    def __init__(self, pos, time):
        self.Position = pos  # Nequick position object
        self.Time = time  # Nequick time object

    def compute(self, solar):
        # Stage 1
        self.__compute_MODIP__()
        # Stage 2
        self.__readccirXXfiles__()
        self.__interpolate_AZR__(solar)
        self.__F2fouriertimeseries__()

        return self.F2Layer()

    def __compute_MODIP__(self):

        latitude = self.Position.latitude
        longitude = self.Position.longitude

        sq, lon_excess, lat_excess = self.Square(latitude, longitude)
        mu2 = self.interpolate2d(sq, lon_excess, lat_excess)
        self.modip = mu2

        return self.modip

    def Square(self, latitude, longitude):
        # what if longitude is [0,360]
        # what if longitude is [-180, 180]
        num_division = 36
        lon_division = 10.0
        lat_division = 5.0

        # Longitude
        lon = (longitude + 180) / lon_division
        lon_start = int(lon) - 2  # range : -2 to 34 or 16 to 52 or -20 to 16
        lon_excess = lon - int(lon)

        # max permissible lon_start is 34
        # min permissible lon_start is 0
        # so this hack is needed
        if (lon_start < 0):
            lon_start = lon_start + num_division
        if (lon_start > (num_division - 3)):
            lon_start = lon_start - num_division

        lat = (latitude + 90.0) / lat_division + 1
        lat_start = int(lat - 1e-6) - 2  # why?
        lat_excess = lat - lat_start - 2

        stModip = np.array(CCIR_MoDIP.modip.stModip)
        square = stModip[lat_start + 1:lat_start + 5, lon_start + 2:lon_start + 6]
        return square, lon_excess, lat_excess

    def interpolate2d(self, Z, x, y):
        assert (np.shape(Z) == (4, 4))

        deltax = 2 * x - 1
        deltay = 2 * y - 1
        # Interpolate horizontally first

        G1 = Z[2, :] + Z[1, :]
        G2 = Z[2, :] - Z[1, :]
        G3 = Z[3, :] + Z[0, :]
        G4 = (Z[3, :] - Z[0, :]) / 3.0

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

    def __readccirXXfiles__(self):

        self.F2 = np.array(CCIR_MoDIP.ccir_f2.F2[self.Time.mth])
        self.Fm3 = np.array(CCIR_MoDIP.ccir_fm3.Fm3[self.Time.mth])

        return self.F2, self.Fm3

    def __interpolate_AZR__(self, solar):
        F2 = self.F2
        Fm3 = self.Fm3

        if solar == 'high':
            self.AF2 = F2[1, :, :]
            self.Am3 = Fm3[1, :, :]
        if solar == 'low':
            self.AF2 = F2[0, :, :]
            self.Am3 = Fm3[0, :, :]

        return self.AF2, self.Am3

    def __F2fouriertimeseries__(self):

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
        modip = self.modip
        latitude = self.Position.latitude
        longitude = self.Position.longitude

        G = []
        m = np.sin(modip * np.pi / 180)
        p = np.cos(latitude * np.pi / 180)
        for i in range(Q[0]):
            G.append(m ** i)

        for i in range(1, len(Q)):
            for j in range(Q[i]):
                G.append(m ** j * p ** i * np.cos(i * longitude * np.pi / 180))
                G.append(m ** j * p ** i * np.sin(i * longitude * np.pi / 180))

        return np.array(G)

    def F2Layer(self):
        CF2 = self.CF2
        Cm3 = self.Cm3

        G = self.__geographical_variation__([12, 12, 9, 5, 2, 1, 1, 1, 1])
        foF2 = np.sum(CF2 * G)
        self.foF2 = foF2

        G = self.__geographical_variation__([7, 8, 6, 3, 2, 1, 1])
        M3000F2 = np.sum(Cm3 * G)
        self.M3000F2 = M3000F2
        return self.foF2, self.M3000F2


def plotmap(lonlon, latlat, mapmap, title, path):
    # Plotting
    plt.figure()
    mapp = Basemap(projection='cyl', llcrnrlat=-90., urcrnrlat=90., \
                   resolution='c', llcrnrlon=-180., urcrnrlon=180.)
    # draw coastlines, country boundaries, fill continents.
    mapp.drawcoastlines()
    mapp.drawstates()
    mapp.drawcountries()

    # -- create and draw meridians and parallels grid lines
    mapp.drawparallels(np.arange(-90., 90., 30.), labels=[1, 0, 0, 0], fontsize=10)
    mapp.drawmeridians(np.arange(-180., 180., 30.), labels=[0, 0, 0, 1], fontsize=10)

    xx, yy = mapp(lonlon, latlat)
    cs = mapp.contourf(xx, yy, mapmap)
    plt.title(title)
    plt.colorbar()
    plt.savefig(path)
    plt.close()
    plt.figure()


res = 100
lats = np.linspace(-70, 70, res)
lons = np.linspace(-180, 180, res)

lonlon, latlat = np.meshgrid(lons, lats)

for mth in range(1, 13):

    if not os.path.exists('foF2/' + str(mth) + '/'):
        os.makedirs('foF2/' + str(mth) + '/')
    if not os.path.exists('M3000F2/' + str(mth) + '/'):
        os.makedirs('M3000F2/' + str(mth) + '/')

    for UT in range(0,23):

        foF2map_high = np.ones([res, res], dtype=np.float32)
        foF2map_low = np.ones([res, res], dtype=np.float32)
        M3000map_low = np.ones([res, res], dtype=np.float32)
        M3000map_high = np.ones([res, res], dtype=np.float32)

        for i in range(res):
            for j in range(res):
                lat = latlat[j, i]
                lon = lonlon[j, i]
                pos = Position(lat, lon)
                time = NEQTime(mth, UT)
                ccir = CCIR_Unpacker(pos, time)
                foF2map_low[j, i], M3000map_low[j, i] = ccir.compute('low')
                foF2map_high[j, i], M3000map_high[j, i] = ccir.compute('high')

        path = os.path.join('foF2', str(mth), "high" + '{:02d}'.format(UT) + '.png')
        plotmap(lonlon, latlat, foF2map_high, 'foF2map: mth= ' + str(mth) + " UT=" + str(UT),path)
        path = os.path.join('foF2', str(mth), "low" + '{:02d}'.format(UT) + '.png')
        plotmap(lonlon, latlat, foF2map_low, 'foF2map: mth= ' + str(mth) + " UT=" + str(UT), path)
        path = os.path.join('M3000F2', str(mth), 'high' + '{:02d}'.format(UT) + '.png')
        plotmap(lonlon, latlat, M3000map_high, 'M3000F2map: mth= ' + str(mth) + " UT=" + str(UT),path)
        path = os.path.join('M3000F2', str(mth), 'low' + '{:02d}'.format(UT) + '.png')
        plotmap(lonlon, latlat, M3000map_low, 'M3000F2map: mth=' + str(mth) + " UT=" + str(UT), path)
