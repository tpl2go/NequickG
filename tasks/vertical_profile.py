""" Plot vertical electron density profile at a given location and time """
import matplotlib.pyplot as plt
import numpy as np
from NequickG import NEQTime, Position, GalileoBroadcast
from NequickG_global import NequickG_global


# Create input objects
TX = NEQTime(4,12)
BX = GalileoBroadcast(236.831,0,0)
RX = Position(40,0)

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
label = "Azr" + str(int(Azr)) + " Apr 12UT 40N 0E"
plt.plot(hs, N, label = label)
plt.xlabel('Height(km)')
plt.ylabel("Electron Density (m^-3)")
plt.title("Nequick-G:\n" + " Apr 12UT ")
plt.grid()
plt.legend()
plt.show()
