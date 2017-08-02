import numpy as np
import CCIR_MoDIP.ccir_fm3
import CCIR_MoDIP.ccir_f2
import CCIR_MoDIP.modip
from aux import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
def __compute_MODIP__(latitude, longitude):
    """
    Reference: Section 2.5.4.3 of GSA's "Ionospheric Correction
    Algorithm for Galileo Single Frequency Users"

    :param latitude: [deg]
    :param longitude: [deg]
    :param stModip: array of MODIP vales
    :return: MODIP [deg]
    """


    stModip = CCIR_MoDIP.modip.stModip

    lngp = 36
    dlatp = 5
    dlngp = 10

    if (latitude > 90):
        modip = 90
        return 90
    elif (latitude < -90):
        modip = -90
        return -90

    lng1 = (longitude + 180.0) / dlngp
    sj = int(lng1) - 2
    dj= lng1 - int(lng1)

    if (sj < 0):
        sj = sj + lngp
    if (sj > (lngp - 3)):
        sj = sj - lngp

    lat1 = (latitude + 90.0) / dlatp + 1
    si = int(lat1 - 1e-6) - 2
    di = lat1 - si - 2

    z = []
    for k in range(1,5):
        z_row = []
        for j in range(1,5):
            z_jk = stModip[si + j][sj+k + 1]
            z_row.append(z_jk)
        z_row.append(di)
        z.append(interpolate(*z_row))
    z.append(dj)
    mu =  interpolate(*z)
    return mu



for mth in range(1,13):
    for UT in range(12,13):
        Azr = 100
        # UT = 12
        res = 100
        F2 = np.array(CCIR_MoDIP.ccir_f2.F2[mth])
        Fm3 = np.array(CCIR_MoDIP.ccir_fm3.Fm3[mth])

        AF2 = F2[0, :, :] * (1 - Azr / 100.0) + F2[1, :, :] * Azr / 100.0
        Am3 = Fm3[0, :, :] * (1 - Azr / 100.0) + Fm3[1, :, :] * Azr / 100.0

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

        CF2 = np.sum(AF2 * y, axis=1)

        # calculate the Fourier time series for M(3000)F2
        cos = np.cos(np.arange(1, 5) * T)
        sin = np.sin(np.arange(1, 5) * T)
        x = np.empty(8)
        x[0::2] = sin
        x[1::2] = cos
        y = np.ones(9)
        y[1:] = x

        Cm3 = np.sum(Am3 * y, axis=1)


        lats = np.linspace(-70,70, res)
        lons = np.linspace(-180,180, res)

        lonlon, latlat = np.meshgrid(lons, lats)

        foF2map = np.ones([res,res], dtype=np.float32)
        M3000map = np.ones([res,res], dtype=np.float32)


        for lat in range(res):
            for lon in range(res):
                latitude = latlat[lat,lon]
                longitude = lonlon[lat,lon]
                modip = __compute_MODIP__(latitude, longitude)

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
                foF2 = np.sum(foF2_n)
##############################################################

                foF2_orders = []
                foF2_orders.append(np.sum(CF2[:12] * M)) # 0th order
                total = 12
                MM = np.empty(24)
                MM[0::2] = M
                MM[1::2] = M
                for i in range(1,9):
                    q = [12, 12, 9, 5, 2, 1, 1, 1, 1][i]

                    CS = np.empty(24)
                    CS[0::2] = C[i]
                    CS[1::2] = S[i]

                    foF2_orders.append(np.sum(CF2[total:total + 2*q] * MM[:2*q] * CS[:2*q]))
                    total +=2*q
                foF2_orders = np.array(foF2_orders)
                foF2_o = np.sum(foF2_orders * P)

#########################################################################
                foF2map[lat,lon] = foF2_o


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
                        M3000F2_n[i] += Cm3[H[i] + 2 * j] * C[i] * (M[j] * P[i])
                        M3000F2_n[i] += Cm3[H[i] + 2 * j + 1] * S[i] * (M[j] * P[i])

                M3000F2 = np.sum(M3000F2_n)

                M3000map[lat,lon] = M3000F2



        plt.figure()
        mapp = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
                          resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
        # draw coastlines, country boundaries, fill continents.
        mapp.drawcoastlines()
        mapp.drawstates()
        mapp.drawcountries()

        #-- create and draw meridians and parallels grid lines
        mapp.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
        mapp.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)

        xx, yy = mapp(lonlon, latlat)
        cs = mapp.contour(xx, yy, foF2map, linewidths=1.5)
        plt.title('foF2map: '+ str(mth) + " " + str(UT))
        plt.savefig('foF2/' + str(mth) + "_" + str(UT) + '.png')
        plt.colorbar()
        plt.close()
        plt.figure()

        mapp = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
                          resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
        # draw coastlines, country boundaries, fill continents.
        mapp.drawcoastlines()
        mapp.drawstates()
        mapp.drawcountries()

        #-- create and draw meridians and parallels grid lines
        mapp.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
        mapp.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)

        xx, yy = mapp(lonlon, latlat)
        cs = mapp.contour(xx, yy, M3000map, linewidths=1.5)
        plt.title('M3000F2: ' + str(mth) + " " + str(UT))
        plt.colorbar()
        plt.savefig('M3000F2/' + str(mth) + "_" + str(UT)+ '.png')
        plt.close()


