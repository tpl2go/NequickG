from NequickG_global import Ray
import csv
import os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

def run(table_type):
    """
    Plot all trajectory of specified in validation table
    """
    print "Plotting trajectories in", table_type, "Solar Activity data table"
    with open(os.path.join('Validation', table_type + '_reference.dat')) as infile:
        reader = csv.reader(infile, delimiter = ' ')
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
        for row in reader:
            row = map(float, row)

            lon1 = row[2]
            lat1 = row[3]
            h1 = row[4] / 1000 # change m to km

            lon2 = row[5]
            lat2 = row[6]
            h2 = row[7] / 1000 # change m to km

            ray = Ray(h1, lat1, lon1, h2, lat2, lon2)

            rs, lats, lons, delta = ray.linspace(100)

            x, y = mapp(lons, lats)
            mapp.plot(x,y, marker='.', linestyle='None')

        plt.savefig(table_type+'.png')

run('Medium')
run('Low')
run('High')
