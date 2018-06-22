from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import time
import calendar
siavashtime = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

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
ncdir = "ncep_nam_20170814_00z/"
ncfile = "MIT_nsf_alpha200m_surf_vel_2017081200_2017081412_2017081612_01h_r01.nc"
#ncfile = "nam.t00z.conusnest.hiresf00.tm00.nc"

root = Dataset(ncdir+ncfile,'r')
vars = root.variables
lato = vars["lat"][:]
lono = vars["lon"][:]
seau = vars["East_vel"][:]
seav = vars["North_vel"][:] 
latorg = 0.5*(lato[0]+lato[-1])
lonorg = 0.5*(lono[0]+lono[-1])
timeo = vars['time'][:]*86400 + siavashtime
timeout = vars['time'][:]
root.close()

yy = np.linspace(-23.9,23.9,240)
xx = np.linspace(-26.2,26.2,263)
xx,yy = np.meshgrid(xx,yy)

windu = np.empty(seau.shape)
windv = np.empty(seav.shape)
#print time.gmtime(timew)
for tt in range(timeo.shape[0]):
    print tt
    ncfile = "nam.t00z.conusnest.hiresf{:d}.tm00.nc".format(tt+12)
    print ncfile
    root = Dataset(ncdir+ncfile,'r')
    vars = root.variables
    lat = vars["latitude"][:]
    lon = vars["longitude"][:]-360
    y = vars["y"][:]
    x = vars["x"][:]
    u = vars["UGRD_10maboveground"][:].squeeze()
    v = vars["VGRD_10maboveground"][:].squeeze()
    timew = vars['time'][:]
    dim = lat.shape
    root.close()
    if timew == timeo[tt]:
        print True
        '''
        ##FIND where to recenter grid
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
        '''
        iorg = 740
        jorg = 1639
        xorg,yorg= lonlat2km(lonorg,latorg,lon[iorg,jorg],lat[iorg,jorg] )
        x = (x - x[jorg])/1000 + xorg
        y = (y - y[iorg])/1000 + yorg
        x,y = np.meshgrid(x,y)
        #Domain size
        gridpointsI = 101
        gridpointsJ = 103
        #Calculate start and end points based on 
        #domain size and origin
        imin = int(np.floor(gridpointsI/2.0))
        jmin = int(np.floor(gridpointsJ/2.0))
        imax = int(np.ceil(gridpointsI/2.0))
        jmax = int(np.ceil(gridpointsJ/2.0))
        
        u = u[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]
        v = v[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]
        x = x[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]
        y = y[(iorg-imin):(iorg+imax),(jorg-jmin):(jorg+jmax)]
        f = interpolate.RectBivariateSpline(y[:,0], x[0,:], u)
        windu[tt,:,:] = f(yy[:,0], xx[0,:])
        del f
        f = interpolate.RectBivariateSpline(y[:,0], x[0,:], v)
        windv[tt,:,:] = f(yy[:,0], xx[0,:])
        del f
    else:
        print False
        break
        


##Write DATA
dim2 = windu.shape
dataset = Dataset('windagedata.nc', mode='w', format='NETCDF4_CLASSIC') 
lat = dataset.createDimension('lat',dim2[1])
lon = dataset.createDimension('lon',dim2[2])
time = dataset.createDimension('time', None)

times = dataset.createVariable('time',np.float64,('time',),fill_value=999)
lats = dataset.createVariable('lat',np.float64,('lat',),fill_value=999)
lons = dataset.createVariable('lon',np.float64,('lon',),fill_value=999)
uo = dataset.createVariable('eastward_vel',np.float64,('time','lat','lon',),fill_value=999)
vo = dataset.createVariable('northward_vel',np.float64,('time','lat','lon',),fill_value=999)
uw = dataset.createVariable('eastward_wind',np.float64,('time','lat','lon',),fill_value=999)
vw = dataset.createVariable('northward_wind',np.float64,('time','lat','lon',),fill_value=999) 

lons.standard_name = 'longitude'
lons.units = 'degree_east'
lons.positive = 'east'
lons._CoordinateAxisType = 'Lon'
lons.axis = 'X'
lons.coordsys = 'geographic'
lons[:] = lono

lats.standard_name = 'latitude'
lats.units = 'degree_north'
lats.positive = 'up'
lats._CoordinateAxisType = 'Lat'
lats.axis = 'Y'
lats.coordsys = 'geographic'
lats[:] = lato

times.standard_name = 'time'
times.long_name = 'time'
times.units = 'days since 2017-06-01 00:00:00 UTC'
times.calendar = 'gregorian'
times._CoordinateAxisType = 'Time'
times[:] = timeo

uo.standard_name = 'surface_eastward_sea_water_velocity'
uo.long_name = 'surface_eastward_sea_water_velocity'
uo.units = 'meter second-1'
uo.coordsys = 'geographic'
uo.positive = 'toward east'
uo.coordinates = 'Longitude Latitude datetime'
uo[:] = seau

vo.standard_name = 'surface_northward_sea_water_velocity'
vo.long_name = 'surface_northward_sea_water_velocity'
vo.units = 'meter second-1'
vo.coordsys = 'geographic'
vo.positive = 'toward north'
vo.coordinates = 'Longitude Latitude datetime'
vo[:] = seav

uw.standard_name = 'eastward_wind'
uw.long_name = 'eastward_wind'
uw.units = 'meter second-1'
uw.coordsys = 'geographic'
uw.positive = 'toward east'
uw.coordinates = 'Longitude Latitude datetime'
uw[:] = windu

vw.standard_name = 'northward_wind'
vw.long_name = 'northward_wind'
vw.units = 'meter second-1'
vw.coordsys = 'geographic'
vw.positive = 'toward north'
vw.coordinates = 'Longitude Latitude datetime'
vw[:] = windv

dataset.close()

