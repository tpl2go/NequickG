""" Plot vertical electron density profile at a given location and time """
import matplotlib.pyplot as plt
import numpy as np
from NequickG import NEQTime, Position, GalileoBroadcast
from NequickG_global import NequickG_global


""" Input Parameters here """
mth = 4
UT = 12
lat = 40
lon = 100
Az = 64

""" Processing Below """
# Create input objects
TX = NEQTime(mth,UT)
BX = GalileoBroadcast(Az,0,0)
RX = Position(lat,lon)

# Create Nequick models
NEQ_global = NequickG_global(TX,BX)
NEQ, Para = NEQ_global.get_Nequick_local(RX)

# Extract information from model
hmin = 100
hmax = 1000
hs = np.arange(hmin, hmax)
N = NEQ.electrondensity(hs)
Azr = Para.Azr

# Plotting
label = ' lat:' + str(lat) + ' lon:' + str(lon)
plt.plot(hs, N, label = label)
plt.plot(Para.hmF2,Para.NmF2 * 10 ** 11, marker = 'o', markersize=5, linestyle='None', label = 'F2 anchor pt')
plt.plot(Para.hmF1,Para.NmF1* 10 ** 11, marker = 'o', markersize=5, linestyle='None', label = 'F1 anchor pt')
plt.plot(Para.hmE,Para.NmE* 10 ** 11, marker = 'o', markersize=5, linestyle='None', label = 'E anchor pt')
plt.xlabel('Height(km)')
plt.ylabel("Electron Density (m^-3)")
plt.title("Nequick-G:\n" + 'mth:' + str(mth) + ' UT:' + str(UT) + ' Azr:' + str(int(Azr)))
plt.grid()
plt.legend()
plt.show()

