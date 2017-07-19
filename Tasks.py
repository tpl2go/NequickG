import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from NequickG import NEQTime, Position, GalileoBroadcast, NequickG, NequickG_global


def Task1():
    # TASK 1: Plot NequickG at 8 different locations

    #mth = 10, universal_time = 12, latitude = 30, longitude = 0
    time = [10, 12]
    positions = [
        [30, 0],
        [45, 0],
        [60, 0],
        [75, 0]]
    broadcasts = [[80,0,0],[193,0,0]]
    TX = NEQTime(*time)
    for bx in broadcasts:
        BX = GalileoBroadcast(*bx)
        for pos in positions:

            RX = Position(*pos)
            NEQ_global = NequickG_global(TX,BX)
            NEQ = NEQ_global.get_Nequick_local(RX)

            hmin = 100
            hmax = 1000
            h = np.arange(hmin, hmax)

            N = NEQ.electrondensity(h)
            #plt.figure()
            plt.plot(h, N)
            plt.xlabel('Height(km)')
            plt.ylabel("Electron Density (m^-3)")
            Azr = NEQ.Para.Azr
            label = "Azr" + str(int(Azr)) + " Oct 12UT " + str(pos[0]) + "N " + str(pos[1]) + "E"
            plt.title("Nequick-G:\n" + label)
            plt.grid()
            #plt.savefig(label)
    plt.savefig("all")


def Task2():
    #test vTEC
    TX = NEQTime(10,12)
    RX = Position(40,0)
    BX = GalileoBroadcast(80,0,0)
    NEQ_global = NequickG_global(TX, BX)
    NEQ = NEQ_global.get_Nequick_local(RX)
    print NEQ.vTEC(100,1000)



def Task3():
    # vTEC map
    BX = GalileoBroadcast(80,0,0)
    vTEC = []
    for lat in np.arange(40, 60, 0.5):
        vTEC_lon = []
        for lon in range(-20, 20, 1):
            pos = [lat, lon]
            TX = NEQTime(10, 12)
            RX = Position(*pos)
            NEQ_global = NequickG_global(TX, BX)
            NEQ = NEQ_global.get_Nequick_local(RX)
            vTEC_lon.append(NEQ.vTEC(100,1000))
        vTEC.append(vTEC_lon)
    plt.figure()
    plt.contour(vTEC)
    plt.show()


def Task4():
    # set up orthographic map projection with
    # perspective of satellite looking down at 50N, 100W.
    # use low resolution coastlines.
    plt.figure()
    mapp = Basemap(projection='ortho',lat_0=43,lon_0=1.4,resolution='l')
    # draw coastlines, country boundaries, fill continents.
    mapp.drawcoastlines(linewidth=0.25)
    mapp.drawcountries(linewidth=0.25)
    mapp.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    mapp.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    mapp.drawmeridians(np.arange(0,360,30))
    mapp.drawparallels(np.arange(-90,90,30))

    BX = GalileoBroadcast(80,0,0)
    vTEC = []
    lats = np.arange(40, 60, 0.5)
    lons = np.arange(-20, 20, 1)
    for lat in lats:
        vTEC_lon = []
        for lon in lons:

            pos = [lat, lon]
            TX = NEQTime(10, 12)
            RX = Position(*pos)
            NEQ_global = NequickG_global(TX, BX)
            NEQ = NEQ_global.get_Nequick_local(RX)

            vTEC_lon.append(NEQ.vTEC(100,1000))
        vTEC.append(vTEC_lon)
    x, y = np.meshgrid(lons, lats)
    xx, yy = mapp(x, y)
    cs = mapp.contour(xx, yy, np.array(vTEC), linewidths=1.5)
    plt.title('vertical Total Electron Count over Toulouse')
    plt.show()




def Task5():
    #test sTEC
    TX = NEQTime(10,12)
    BX = GalileoBroadcast(80,0,0)
    NEQ_global = NequickG_global(TX, BX)
    stec = NEQ_global.sTEC(100,40,0,1000,40,0)
    print stec

Task2()
Task5()
