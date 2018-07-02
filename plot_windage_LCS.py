#import h5py as hp
from netCDF4 import Dataset
import matplotlib.pyplot as plt
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
nthresh = 3
wthresh = 4
nthresh = 5
wthresh = 5
root = Dataset('MV_FTLE_-6hrs_NoWindage.nc','r')
londfile = root.variables
print londfile.keys()
nlat = londfile['initial_lat'][:].reshape([ydim,xdim])
nlon = londfile['initial_lon'][:].reshape([ydim,xdim])
nftle = londfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
ncg = londfile['CauchyGreen'][:,timestep,:,:].reshape([ydim,xdim,2,2])
t = londfile['time'][:]
print time.gmtime(t[0]*24*60*60+tstart)
nftle = ma.masked_where(nftle==999,nftle)


root = Dataset('MV_FTLE_-6hrs_Windage=0,019.nc','r')
londfile = root.variables
wlat = londfile['initial_lat'][:].reshape([ydim,xdim])
wlon = londfile['initial_lon'][:].reshape([ydim,xdim])
wftle = londfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
wcg = londfile['CauchyGreen'][:,timestep,:,:].reshape([ydim,xdim,2,2])

wftle = ma.masked_where(wftle==999,wftle)


ndfdy,ndfdx = np.gradient(nftle,dy/2.0,dx/2.0,edge_order=2)
ndfdydy,ndfdydx = np.gradient(ndfdy,dy/2.0,dx/2.0,edge_order=2)
ndfdxdy,ndfdxdx = np.gradient(ndfdx,dy/2.0,dx/2.0,edge_order=2)
wdfdy,wdfdx = np.gradient(wftle,dy/2.0,dx/2.0,edge_order=2)
wdfdydy,wdfdydx = np.gradient(wdfdy,dy/2.0,dx/2.0,edge_order=2)
wdfdxdy,wdfdxdx = np.gradient(wdfdx,dy/2.0,dx/2.0,edge_order=2)

ndirdiv = np.ma.empty([ydim,xdim])
wdirdiv = np.ma.empty([ydim,xdim])
nconcav = np.ma.empty([ydim,xdim])
wconcav = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        if (ndfdx[i,j] and ndfdy[i,j] and ndfdxdy[i,j] and ndfdydy[i,j] and ndfdxdx[i,j] and ndfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(ncg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            ndirdiv[i,j] = np.dot([ndfdx[i,j],ndfdy[i,j]],eigenVectors[:,0])
            nconcav[i,j] = np.dot(np.dot([[ndfdxdx[i,j],ndfdxdy[i,j]],[ndfdydx[i,j],ndfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        else:
            ndirdiv[i,j] = np.ma.masked
            nconcav[i,j] = np.ma.masked

        if (wdfdx[i,j] and wdfdy[i,j] and wdfdxdy[i,j] and wdfdydy[i,j] and wdfdxdx[i,j] and wdfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(wcg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            wdirdiv[i,j] = np.dot([wdfdx[i,j],wdfdy[i,j]],eigenVectors[:,0])
            wconcav[i,j] = np.dot(np.dot([[wdfdxdx[i,j],wdfdxdy[i,j]],[wdfdydx[i,j],wdfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        else:
           wdirdiv[i,j] = np.ma.masked
           wconcav[i,j] = np.ma.masked



ndirdiv = np.ma.masked_where(nconcav>0,ndirdiv)
ndirdiv = np.ma.masked_where(nftle<=nthresh,ndirdiv)
wdirrdiv = np.ma.masked_where(wconcav>0,wdirdiv)
wdirdiv = np.ma.masked_where(wftle<=wthresh,wdirdiv)




nftle = ma.masked_where(nftle<nthresh,nftle)
wftle = ma.masked_where(wftle<wthresh,wftle)

lon_min = wlon.min()
lon_max = wlon.max()
lat_min = wlat.min()
lat_max = wlat.max()
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

nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],latlon=True)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],latlon=True)

'''
npp = nridge.collections[0].get_paths()
wpp = wridge.collections[0].get_paths()
nindex = [104,75,31,32,26] #no windage
windex = [78,54,16,17,10,19] # windage = 0.019
nv = np.empty([0,2])
wv = np.empty([0,2])
plt.close('all')
for i in range(len(nindex)):
    if i<3:
        nv = np.concatenate((nv,list(reversed(npp[nindex[i]].vertices))))
    else:
        nv = np.concatenate((nv,npp[nindex[i]].vertices))

for i in range(len(windex)):
    if i<3:
        wv = np.concatenate((wv,list(reversed(wpp[windex[i]].vertices))))
    else:
        wv = np.concatenate((wv,wpp[windex[i]].vertices))
        
nx,ny = m(nv[:,0],nv[:,1],inverse=True)
wx,wy = m(wv[:,0],wv[:,1],inverse=True)
m.plot(nx,ny, latlon=True,color='r')
m.plot(wx,wy, latlon=True,color='b')
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
'''
'''
plt.figure(2)
plt.subplot(121)
m.pcolormesh(wlon,wlat,wftle,latlon=True,vmin=0, vmax=wftle.max(),cmap='Blues')#,alpha=0.4)
m.pcolormesh(nlon,nlat,nftle,latlon=True,vmin=0, vmax=nftle.max(),cmap='Reds')#,alpha=0.4)
#lon, lat = np.meshgrid(lon,lat,indexing='ij')
#m.pcolormesh(lon,lat,ftle,latlon=True,shading='gourand')
#m.contourf(lon,lat,ftle,latlon=True,levels=np.linspace(0,ftle.max(),3001))
#plt.colorbar()
m.drawcoastlines()
#m.drawrivers()
#m.drawstates()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title("Red: No windage, Blue: windage = 0.019; both fields masked under 3 days^{-1}")

plt.subplot(122)
m.pcolormesh(nlon,nlat,nftle,latlon=True,vmin=0, vmax=nftle.max(),cmap='Reds')#,alpha=0.4)
m.pcolormesh(wlon,wlat,wftle,latlon=True,vmin=0, vmax=wftle.max(),cmap='Blues')#,alpha=0.4)
#lon, lat = np.meshgrid(lon,lat,indexing='ij')
#m.pcolormesh(lon,lat,ftle,latlon=True,shading='gourand')
#m.contourf(lon,lat,ftle,latlon=True,levels=np.linspace(0,ftle.max(),3001))
#plt.colorbar()
m.drawcoastlines()
#m.drawrivers()
#m.drawstates()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
'''

'''
plt.figure(1)
plt.subplot(111)

m.drawcoastlines()
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],colors='blue',latlon=True,alpha=0.6)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
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
















