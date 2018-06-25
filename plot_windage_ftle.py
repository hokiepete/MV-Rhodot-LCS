#import h5py as hp
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import calendar
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
ydim=240
xdim=263
timestep=0
#with hp.File('withwindage.nc','r') as loadfile:
root = Dataset('MV_FTLE_-6hrs_NoWindage.nc','r')
loadfile = root.variables
#print loadfile.keys()
nlat = loadfile['initial_lat'][:].reshape([ydim,xdim])
nlon = loadfile['initial_lon'][:].reshape([ydim,xdim])
nftle = loadfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
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
llat = loadfile['initial_lat'][:].reshape([ydim,xdim])
llon = loadfile['initial_lon'][:].reshape([ydim,xdim])
lftle = loadfile['FTLE'][:,timestep,0].reshape([ydim,xdim])
#grad = loadfile['VelocityGradient'][:]
#time = loadfile['time'][:]

lftle[lftle<0]=0
lftle = ma.masked_where(lftle==999,lftle)

plt.close('all')
plt.figure(1)
plt.subplot(221)
plt.title("no windage")
plt.contourf(nlon,nlat,nftle,levels=np.linspace(0,nftle.max(),301),cmap='viridis')
plt.colorbar()
plt.subplot(222)
plt.title("windage = 0.019")
plt.contourf(llon,llat,lftle,levels=np.linspace(0,lftle.max(),301),cmap='viridis')
plt.colorbar()
plt.subplot(223)
plt.title("Red: No windage, Blue: windage = 0.019; both fields masked under 3 days^{-1}")
nftle = ma.masked_where(nftle<3,nftle)
plt.pcolormesh(nlon,nlat,nftle,vmin=0, vmax=nftle.max(),cmap='Reds')
lftle = ma.masked_where(lftle<3,lftle)
plt.pcolormesh(llon,llat,lftle,vmin=0, vmax=lftle.max(),cmap='Blues')
#plt.colorbar()

plt.figure(2)
lon_min = llon.min()
lon_max = llon.max()
lat_min = llat.min()
lat_max = llat.max()
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

m.pcolormesh(llon,llat,lftle,latlon=True,vmin=0, vmax=lftle.max(),cmap='Blues')
m.pcolormesh(nlon,nlat,nftle,latlon=True,vmin=0, vmax=nftle.max(),cmap='Reds')
#lon, lat = np.meshgrid(lon,lat,indexing='ij')
#m.pcolormesh(lon,lat,ftle,latlon=True,shading='gourand')
#m.contourf(lon,lat,ftle,latlon=True,levels=np.linspace(0,ftle.max(),3001))
#plt.colorbar()
m.drawcoastlines()
#m.drawrivers()
#m.drawstates()

#m.drawcountries()
parallels = np.linspace(lat_min,lat_max,6)
#parallels = np.arange(lat_min,lat_max,6)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.linspace(lon_min,lon_max,6)
#meridians = np.arange(lon_min,lon_max,6)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title("Red: No windage, Blue: windage = 0.019; both fields masked under 3 days^{-1}")
plt.show()


'''
plt.subplot(224)
plt.title("windage = 0.19")
plt.pcolor(glon,glat,gftle,vmin=0,cmap='viridis')
plt.colorbar()
'''