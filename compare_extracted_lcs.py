# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 15:51:23 2018

@author: pnolan86
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from numpy import genfromtxt
ndata = np.genfromtxt('nowindage.csv', delimiter=',')
wdata = np.genfromtxt('windage0,019.csv', delimiter=',')


lon_min = -70.916465753
lon_max = -70.2909317079
lat_min = 41.0250244184
lat_max = 41.4548797564
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )


m.plot(ndata[:,0],ndata[:,1], latlon=True,color='r')
m.plot(wdata[:,0],wdata[:,1], latlon=True,color='b')
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

ax = plt.gca()
def format_coord(x, y):
    return 'x=%.4f, y=%.4f'%(m(x, y, inverse = True))
ax.format_coord = format_coord
plt.show()