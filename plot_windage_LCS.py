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
wthresh = 3
#with hp.File('withwindage.nc','r') as loadfile:
root = Dataset('MV_FTLE_-6hrs_NoWindage.nc','r')
loadfile = root.variables
print loadfile.keys()
#print loadfile.keys()
nlat = loadfile['initial_lat'][:].reshape([ydim,xdim])
nlon = loadfile['initial_lon'][:].reshape([ydim,xdim])
nftle = loadfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
ncg = loadfile['CauchyGreen'][:,timestep,:,:].reshape([ydim,xdim,2,2])
#grad = loadfile['VelocityGradient'][:]
t = loadfile['time'][:]
#print loadfile['time'].units
#epoch 2017-06-01 00:00:00 UTC    
#tstart = tt.mktime(tt.strptime("01.06.2017 00:00:00", "%d.%m.%Y %H:%M:%S")) UTC
#tstart = tt.mktime(tt.strptime("31.05.2017 22:00:00", "%d.%m.%Y %H:%M:%S"))
#print tt.gmtime(tstart)
print time.gmtime(t[0]*24*60*60+tstart)
nftle[nftle<0]=0
nftle = ma.masked_where(nftle==999,nftle)


#with hp.File('MV_FTLE_-6hrs_Windage=0,019.nc','r') as loadfile:
root = Dataset('MV_FTLE_-6hrs_Windage=0,019.nc','r')
loadfile = root.variables
#print loadfile.keys()
wlat = loadfile['initial_lat'][:].reshape([ydim,xdim])
wlon = loadfile['initial_lon'][:].reshape([ydim,xdim])
wftle = loadfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
wcg = loadfile['CauchyGreen'][:,timestep,:,:].reshape([ydim,xdim,2,2])
#grad = loadfile['VelocityGradient'][:]
#time = loadfile['time'][:]

wftle[wftle<0]=0
wftle = ma.masked_where(wftle==999,wftle)


ndfdy,ndfdx = np.gradient(nftle,dy/2.0,dx/2.0,edge_order=2)
wdfdy,wdfdx = np.gradient(wftle,dy/2.0,dx/2.0,edge_order=2)


ndfdy,ndfdx = np.gradient(nftle,dy/2.0,dx/2.0,edge_order=2)
wdfdy,wdfdx = np.gradient(wftle,dy/2.0,dx/2.0,edge_order=2)


ndirdiv = np.ma.empty([ydim,xdim])
wdirdiv = np.ma.empty([ydim,xdim])
for i in range(ydim):
    for j in range(xdim):
        if (ndfdx[i,j] and ndfdy[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(ncg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            ndirdiv[i,j] = np.dot([ndfdx[i,j],ndfdy[i,j]],eigenVectors[:,0])
            
        else:
            ndirdiv[i,j] = np.ma.masked
        if (wdfdx[i,j] and wdfdy[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(wcg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            wdirdiv[i,j] = np.dot([wdfdx[i,j],wdfdy[i,j]],eigenVectors[:,0])
            
        else:
            wdirdiv[i,j] = np.ma.masked

ndirdiv = np.ma.masked_where(nftle<=nthresh,ndirdiv)
wdirdiv = np.ma.masked_where(wftle<=wthresh,wdirdiv)


plt.close('all')
plt.figure(1)
ax=plt.subplot(121)
plt.title("no windage")
plt.contourf(nlon,nlat,nftle,levels=np.linspace(0,nftle.max(),301),cmap='viridis')
plt.colorbar()
ax.set_aspect('equal', adjustable='box', anchor='C')
ax=plt.subplot(122)
plt.title("windage = 0.019")
plt.contourf(wlon,wlat,wftle,levels=np.linspace(0,wftle.max(),301),cmap='viridis')
ax.set_aspect('equal', adjustable='box', anchor='C')
plt.colorbar()
'''
plt.subplot(223)
plt.title("Red: No windage, Blue: windage = 0.019; both fields masked under 3 days^{-1}")
'''
nftle = ma.masked_where(nftle<3,nftle)
#plt.pcolormesh(nlon,nlat,nftle,vmin=0, vmax=nftle.max(),cmap='Reds')
wftle = ma.masked_where(wftle<3,wftle)
#plt.pcolormesh(wlon,wlat,wftle,vmin=0, vmax=wftle.max(),cmap='Blues')
#plt.colorbar()
plt.figure(2)
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
            resolution = 'c',
            area_thresh=0.,
            )
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
'''
parallels = np.arange(41.1,lat_max+0.1,0.1)
meridians = np.arange(-70.2,lon_min-0.1,-0.1)
'''
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
'''
parallels = np.arange(41.1,lat_max+0.1,0.1)
meridians = np.arange(-70.2,lon_min-0.1,-0.1)
'''
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

plt.show()


'''
plt.subplot(224)
plt.title("windage = 0.19")
plt.pcolor(glon,glat,gftle,vmin=0,cmap='viridis')
plt.colorbar()
'''