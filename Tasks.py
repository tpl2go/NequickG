import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from NequickG import Position, GalileoBroadcast, NequickG


def Task1():
    # TASK 1: Plot NequickG at 8 different locations

    #mth = 10, universal_time = 12, latitude = 30, longitude = 0
    positions = [
        [10, 12, 30, 0],
        [10, 12, 45, 0],
        [10, 12, 60, 0],
        [10, 12, 75, 0]]
    broadcasts = [[80,0,0],[193,0,0]]

    for pos in positions:
        for bx in broadcasts:
            RX = Position(*pos)
            BX = GalileoBroadcast(*bx)
            NEQ = NequickG(RX,BX)

            hmin = 100
            hmax = 1000
            h = np.arange(hmin, hmax)

            N = NEQ.vertical_electrondensity(h)
            #plt.figure()
            plt.plot(h, N)
            plt.xlabel('Height(km)')
            plt.ylabel("Electron Density (m^-3)")
            Azr = NEQ.Para.Azr
            label = "Azr" + str(int(Azr)) + " Oct 12UT " + str(pos[2]) + "N " + str(pos[3]) + "E"
            plt.title("Nequick-G:\n" + label)
            plt.grid()
            #plt.savefig(label)
    plt.savefig("all")
# Task1()

def Task2():
    #test vTEC
    RX = Position(10,12,30,0)
    BX = GalileoBroadcast(80,0,0)
    NEQ = NequickG(RX, BX)
    print NEQ.vTEC(100,1000)

# Task2()

def Task3():
    # vTEC map
    BX = GalileoBroadcast(80,0,0)
    vTEC = []
    for lat in np.arange(40, 60, 0.5):
        vTEC_lon = []
        for lon in range(-20, 20, 1):
            pos = [10, 12] + [lat, lon]
            RX = Position(*pos)
            NEQ = NequickG(RX,BX)
            vTEC_lon.append(NEQ.vTEC(100,1000))
        vTEC.append(vTEC_lon)

    plt.contour(vTEC)
    plt.show()

# Task3()

def Task4():
    # set up orthographic map projection with
    # perspective of satellite looking down at 50N, 100W.
    # use low resolution coastlines.
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
            pos = [10, 12] + [lat, lon]
            RX = Position(*pos)
            NEQ = NequickG(RX,BX)
            vTEC_lon.append(NEQ.vTEC(100,1000))
        vTEC.append(vTEC_lon)
    x, y = np.meshgrid(lons, lats)
    xx, yy = mapp(x, y)
    cs = mapp.contour(xx, yy, np.array(vTEC), linewidths=1.5)
    plt.title('contour lines over filled continent background')
    plt.show()

Task4()



