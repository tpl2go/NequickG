"""
Time and Solar activity fully defines a global Nequick model.
From the global model, Nequick processing parameters are computed at each geographical location.
This task visualises the spatial variation in each Nequick processing parameter.
Useful for understanding the physical meaning behind each processing parameter and as a debugging tool
"""
import matplotlib.pyplot as plt
import numpy as np
from NequickG import NEQTime, Position, GalileoBroadcast
from NequickG_global import NequickG_global
from mpl_toolkits.basemap import Basemap
import os

def plotparameters(TX, BX, path):
    # processing parameters to visualise
    attrs = ['foF1', 'foF2', 'foE', 'M3000F2', 'NmF2', 'NmF1', 'NmE', 'hmE', 'hmF1', 'hmF2', 'modip',
             'Az','Azr', 'solarsine', 'solarcosine', 'chi', 'chi_eff', 'H0', 'B1bot', 'B1top', 'B2bot', 'BEtop', 'BEbot'
             ,'A1', 'A2', 'A3', 'k', 'vTEC', 'seasp']


    NEQ_global = NequickG_global(TX, BX)

    latlat, lonlon, outs = NEQ_global.map_parameters(attrs, -70, -180, 70, 180, resolution=150)


    for i in range(len(outs)):
        plt.figure()
        mapp = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
                  resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
        #-- draw coastlines, state and country boundaries, edge of map
        mapp.drawcoastlines()
        mapp.drawstates()
        mapp.drawcountries()

        #-- create and draw meridians and parallels grid lines
        mapp.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
        mapp.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)
        xx, yy = mapp(lonlon, latlat)

        # filled contour
        cs = mapp.contourf(xx, yy, outs[i])

        plt.title(attrs[i])
        mapp.colorbar(cs)

        if not os.path.exists(path):
            os.makedirs(path)

        plt.savefig(os.path.join(path,  attrs[i]+'.png'))
        plt.close()


# plot parameters for Validation test scenarios
solar = 'Medium' # or 'Medium' or 'Low'

if solar=='Low':
    BX = GalileoBroadcast(2.580271,0.127628236,0.0252748384) # Low solar activity
elif solar=='Medium':
    BX = GalileoBroadcast(121.129893,0.351254133,0.0134635348)# Medium solar activity
elif solar == 'High':
    BX = GalileoBroadcast(236.831641, -0.39362878, 0.00402826613)# High solar activity
else:
    raise ValueError

for time in [0, 4, 8, 12, 16, 20]:
    TX = NEQTime(4, time)
    plotparameters(TX, BX, os.path.join('maps', solar, str(time)))
