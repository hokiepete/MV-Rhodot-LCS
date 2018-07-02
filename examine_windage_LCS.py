#import h5py as hp
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import calendar
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ydim=(240-1)*2+1
xdim=(263-1)*2+1
timestep=0
dx = 200
dy = 200
thresh = 6
#with hp.File('withwindage.nc','r') as londfile:
#root = Dataset('MV_FTLE_-6hrs_NoWindage.nc','r')
root = Dataset('MV_FTLE_-6hrs_Windage=0,019.nc','r')
loadfile = root.variables
print loadfile.keys()
#print londfile.keys()
lat = loadfile['initial_lat'][:].reshape([ydim,xdim])
lon = loadfile['initial_lon'][:].reshape([ydim,xdim])
ftle = loadfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
cg = loadfile['CauchyGreen'][:,timestep,:,:].reshape([ydim,xdim,2,2])
#grad = londfile['VelocityGradient'][:]
t = loadfile['time'][:]
print time.gmtime(t[0]*24*60*60+tstart)
#nftle[nftle<0]=0
ftle = ma.masked_where(ftle==999,ftle)


dfdy,dfdx = np.gradient(ftle,dy/2.0,dx/2.0,edge_order=2)
dfdydy,dfdydx = np.gradient(dfdy,dy/2.0,dx/2.0,edge_order=2)
dfdxdy,dfdxdx = np.gradient(dfdx,dy/2.0,dx/2.0,edge_order=2)

dirdiv = np.ma.empty([ydim,xdim])
concav = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in  range(xdim):
        if (dfdx[i,j] and dfdy[i,j] and dfdxdy[i,j] and dfdydy[i,j] and dfdxdx[i,j] and dfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(cg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            dirdiv[i,j] = np.dot([dfdx[i,j],dfdy[i,j]],eigenVectors[:,0])
            concav[i,j] = np.dot(np.dot([[dfdxdx[i,j],dfdxdy[i,j]],[dfdydx[i,j],dfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        #print nconcav[i,j]
        else:
            dirdiv[i,j] = np.ma.masked
            concav[i,j] = np.ma.masked




dirdiv = np.ma.masked_where(concav>0,dirdiv)
dirdiv = np.ma.masked_where(ftle<=thresh,dirdiv)


plt.close('all')
lon_min = lon.min()
lon_max = lon.max()
lat_min = lat.min()
lat_max = lat.max()
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
'''
plt.figure(1)
ridge = m.contour(lon,lat,dirdiv,levels =[0],colors='blue',latlon=True,alpha=0.6)
m.drawcoastlines()

parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax = plt.gca()
def format_coord(x, y):
    return 'x=%.4f, y=%.4f'%(m(x, y, inverse = True))
ax.format_coord = format_coord
plt.show()
'''

matplotlib.use('Agg')
ridgelines = m.contour(lon,lat,dirdiv,levels =[0], latlon=True)
pp = ridgelines.collections[0].get_paths()
plt.close('all')


###View contours
plt.close('all')
for p in range(len(pp)):
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    if x.size > 5:
        m.plot(x,y)#, latlon=True)
        m.drawcoastlines()
        plt.title('{:03d}'.format(p))
        plt.savefig('{:03d}.png'.format(p))
        plt.close('all')
        
ridgelines = m.contour(lon,lat,dirdiv,levels =[0], latlon=True)
plt.show()        
'''
"""
import numpy as np
##OUTPUT desired contours
#from operator import itemgetter 
#index = [104,75,31,32,26] #no windage
index = [78,54,16,17,10,19] # windage = 0.019
v = np.empty([0,2])
for i in range(len(index)):
    if i<3:
        v = np.concatenate((v,list(reversed(pp[index[i]].vertices))))
    else:
        v = np.concatenate((v,pp[index[i]].vertices))
    #v = nvpp[index[i]].vertices

x,y = m(v[:,0],v[:,1],inverse=True)
#x = v[:,0]
#y = v[:,1]
m.plot(x,y, latlon=True)
m.drawcoastlines()


del v
import csv

v = np.stack((x,y),axis=1)
with open("windage0,019.csv", "w") as f:
    writer = csv.writer(f)
    for row in v:
        writer.writerow(row)
    f.close()










