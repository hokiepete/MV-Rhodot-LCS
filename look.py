from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def cot(th):
    return 1.0/np.tan(th)

def sec(th):
    return 1.0/np.cos(th)

def deg2rad(deg):
    return deg*np.pi/180.0

def rad2deg(rad):
    return rad*180.0/np.pi

def lonlat2km(reflon,reflat,lon,lat ):
    #LONLAT2KM Summary of this function goes here
    #   Uses Lambert Conformal Projection
    #stdlat1  =  deg2rad(30)
    #stdlat2  =  deg2rad(60)
    stdlat1  =  deg2rad(40)
    stdlat2  =  deg2rad(42)
    R=6371
    reflon = deg2rad(reflon)
    reflat = deg2rad(reflat)
    lon = deg2rad(lon)
    lat = deg2rad(lat)
    n = np.log(np.cos(stdlat1)*sec(stdlat2)) / np.log(np.tan(0.25*np.pi+0.5*stdlat2)*cot(0.25*np.pi+0.5*stdlat1))
    F=(np.cos(stdlat1)*(np.tan(0.25*np.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*np.pi+0.5*reflat)**n)
    p = R*F*(cot(0.25*np.pi+0.5*lat)**n)
    th = n*(lon-reflon)
    x=p*np.sin(th)
    y=p0-p*np.cos(th)
    return x,y

ncfile = "MIT_nsf_alpha200m_surf_vel_2017081200_2017081412_2017081612_01h_r01.nc"
#ncfile = "nam.t00z.conusnest.hiresf00.tm00.nc"
root = Dataset(ncfile,'r')
vars = root.variables
lat = vars["lat"][:]
lon = vars["lon"][:]
latorg = 0.5*(lat[0]+lat[-1])
lonorg = 0.5*(lon[0]+lon[-1])
timeo = vars['time'][:]*86400
#lon, lat = np.meshgrid(lon,lat)
yy = np.linspace(-23.9,23.9,24)
xx = np.linspace(-26.2,26.2,26)
#yy = np.linspace(-99.9,99.9,240)
#xx = np.linspace(-100.2,100.2,260)
xx,yy = np.meshgrid(xx,yy)

import time
import calendar
tt = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
tp = timeo[0]+tt
print tp
print time.gmtime(tp)
'''
time
shape = (49,)
dimension[time] = 49
number of attributes = 6
attribute[standard_name] = time
attribute[long_name] = time
attribute[units] = days since 2017-06-01 00:00:00 UTC
attribute[calendar] = gregorian
attribute[_CoordinateAxisType] = Time
attribute[_FillValue] = 1e+35
'''
#print time.gmtime(timeo[0]*86400)#1346114717972/1000.)


#plt.scatter(lon,lat,color='r')

#ncfile = "MIT_nsf_alpha200m_surf_vel_2017081200_2017081412_2017081612_01h_r01.nc"
ncfile = "nam.t00z.conusnest.hiresf00.tm00.nc"
root = Dataset(ncfile,'r')
vars = root.variables
#keys=vars.keys()
#gridspacing = 3
lat = vars["latitude"][:]
lon = vars["longitude"][:]-360
y = vars["y"][:]
x = vars["x"][:]
u = vars["UGRD_10maboveground"][:]
timew = vars['time'][:]
print time.gmtime(timew)
print u.shape
dim = lat.shape
#plt.scatter(lon,lat)

tol = 0.013
#Find the index of the origin
print "want lat lon", latorg, lonorg
for i in range(dim[0]):
    for j in range(dim[1]):
        #print lat[i,j], lon[i,j]
        if np.fabs(lat[i,j]-latorg)<tol and np.fabs(lon[i,j]-lonorg)<tol:
            print 'closest we can get is', lat[i,j], lon[i,j]
            print i,j
            iorg = i
            jorg = j

xorg,yorg= lonlat2km(lonorg,latorg,lon[iorg,jorg],lat[iorg,jorg] )

x = (x - x[jorg])/1000 + xorg
y = (y - y[iorg])/1000 + yorg
print x[jorg]
print y[iorg]

x,y = np.meshgrid(x,y)
print x[iorg,jorg]
print y[iorg,jorg]

#Domain size
gridpointsI = 101
gridpointsJ = 103
#Calculate start and end points based on 
#domain size and origin
imin = int(np.floor(gridpointsI/2.0))
jmin = int(np.floor(gridpointsJ/2.0))
imax = int(np.ceil(gridpointsI/2.0))
jmax = int(np.ceil(gridpointsJ/2.0))

u = u[0,(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]#.ravel()
#v = vvar[:,(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]#.ravel()
x = x[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]
y = y[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]

#f = interpolate.interp2d(x.ravel(), y.ravel(), u.ravel(), kind='cubic')
#f = interpolate.interp2d(x, y, u, kind='cubic
#f = interpolate.RectBivariateSpline(x[0,:], y[:,0], u.T)
f = interpolate.RectBivariateSpline(y[:,0], x[0,:], u)

#uout = f(xx.ravel(), yy.ravel())
#uout = f(xx[0,:], yy[:,0])
uout = f(yy[:,0], xx[0,:])
dim2 = uout.shape

for i in range(dim2[0]):
    for j in range(dim2[1]):
        i
#uout = uout.T
#import matplotlib.pyplot as plt
u = np.ma.masked_where(x>xx.max(),u)
u = np.ma.masked_where(x<xx.min(),u)
u = np.ma.masked_where(y>yy.max(),u)
u = np.ma.masked_where(y<yy.min(),u)
#x = np.ma.masked_where(x>26.2,x)
#x = np.ma.masked_where(x<-26.2,x)
#x = x[x<=26.2]
#x = x[x>=-26.2]
#y = y[y<=23.9]
#y = y[y>=-23.9]
plt.close('all')
fig = plt.figure(1)
ax = plt.pcolormesh(x[0,:], y[:,0],u,vmin=u.min(),vmax=u.max())
plt.colorbar()

fig = plt.figure(2)
ax = plt.pcolormesh(xx[0,:], yy[:,0],uout,vmin=u.min(),vmax=u.max())
plt.colorbar()

