# driver_5Perc_MC_v3.py
#
# Plots centered-composites for VRILEs
#
# Madeline Frank
# December 2024
import numpy as np
import netCDF4
import pickle
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import datetime as dt
#from mpl_toolkits.basemap import Basemap
#from time import strptime
import time
import os.path
import struct
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from global_land_mask import globe
import cftime
import scipy.ndimage as ndimage

import nsidc_extent
import seaIceMasks
import helpers
import regrid_map
import bootstrap
import polar_project

#geoms = fiona.open(shpreader.natural_earth(resolution='110m', category='physical', name='land'))
#land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry']) for geom in geoms])
land_shp_fname = shpreader.natural_earth(resolution='110m', category='physical', name='land')
land_geom = unary_union(list(shpreader.Reader(land_shp_fname).geometries()))

land = prep(land_geom)

timeStringFormat = "%Y-%m-%d-%H"

def is_land(x, y):
  return land.contains(sgeom.Point(x, y))

# 6/22 MC: getting my vRILE locations
def vrile_data(dataDir='/raid3/maddie/Data/', iceFile='summer_vRILE_withboth_annual_bwfilter_location_closest.txt'):
    # vRILE[i] is a 3 variable list [datetimeobj, lat, lon]
    lat,lon = np.loadtxt(dataDir + iceFile,usecols=(1,2),dtype='str',delimiter = ', ',unpack = True)
    year,month,day = np.loadtxt(dataDir + iceFile,usecols=(0,1,2),dtype='str',delimiter='-',unpack=True)        
    vRILE = []
    for itr in range(len(lat)):
        tmpDate = dt.datetime(int(year[itr]),int(month[itr]),int(day[itr]))
        tmpLat = float(lat[itr]); tmpLon = float(lon[itr])
        vRILE.append([tmpDate, tmpLat, tmpLon])
    return vRILE


# 6/25 MC: for mslp, key should be 'msl' and infoKey 'sfc'.  For theta on trop, key is 'pt' and infoKey is 'tpvTracker'
def driver_values_localRef_cases(infoCases, dtLag=0, refFrame='sfcWind',isAtmo=True, key='msl', infoKey='sfc',
             rCutoff_ice_sfcCyclone=1000., rCutoff_tpv_sfcCyclone=1000., rEarth=6370.,calcMonthlyAnomaly=True):
  vals_allCases = []
  nCases = len(infoCases)
  for iCase in range(nCases):
    #unpack information for case
    #6/22 MC: changing this to match vRILE = [datetime,lat,lon] and changing tEvent to not care about the wind
    #2/11/19 MC: if *_windLoc.pkl, then infoCases = [event time, event time - 5 days, lat, lon] and you have to change the lat and lon indicies
    t0 = infoCases[iCase][0]; lat_windMax = infoCases[iCase][1]; lon_windMax=infoCases[iCase][2];
    #rotationAngle = angle_windMax
    tEvent = t0 - dt.timedelta(days=dtLag)
    if tEvent.year > 2022:
      continue
    if (refFrame=='min_DT-MSLP'):
      min_sfc_val, min_sfc_lat, min_sfc_lon, min_DT_val, min_DT_lat, min_DT_lon, angle_DT2MSLP = find_mslp_DT_min(tEvent, lat_windMax, lon_windMax, rCutoff_ice_sfcCyclone, rCutoff_tpv_sfcCyclone, rEarth, calcMonthlyAnomaly=calcMonthlyAnomaly)
      lat_windMax = min_sfc_lat; lon_windMax = min_sfc_lon; angle_windMax = angle_DT2MSLP
      #lat_windMax = min_DT_lat; lon_windMax = min_DT_lon; angle_windMax = angle_DT2MSLP
      rotationAngle = -angle_windMax
 
    if (isAtmo):
      f = helpers.get_atmo_filename_year(tEvent, info=infoKey)
      data = netCDF4.Dataset(f,'r')
      fileTimes = data.variables['time'][:]
      iTime = helpers.eraI_time2TimeInd(fileTimes,tEvent)
      #if key == 'msl':
      #  vals = data.variables[key][iTime,0,:,:]
      #else:
      vals = data.variables[key][iTime,:,:]
      lat1d = data.variables['latitude'][:]
      lon1d = data.variables['longitude'][:]
      data.close()      
      if (calcMonthlyAnomaly): #subtract off 1979-2014 mean for that month
        valsRef,valsStd,mthly_lat,mthly_lon = helpers.get_monthlyMean_andStd_ERA5(key=key,scaleFac=1.0,iMonth=tEvent.month-1)
        [minLat] = np.where(mthly_lat == np.min(lat1d))    
        step = int(len(mthly_lon)/len(lon1d))
        vals = (vals-valsRef[0:minLat[0]+1:step,::step])/valsStd[0:minLat[0]+1:step,::step]
      lon2d,lat2d = np.meshgrid(lon1d,lat1d)
    else:
      f = '/raid4/maddie/era5_flux_into_waves_daily_50N.nc'
      data = netCDF4.Dataset(f,'r')
      fileTimes = data.variables['time'][:]
      iTime = helpers.eraI_time2TimeInd(fileTimes,tEvent)
      vals = data.variables['phiaw'][iTime,:,:]
      lat1d = data.variables['latitude'][:]
      lon1d = data.variables['longitude'][:]
      data.close()
      vals = np.ma.masked_values(vals, -32767.)
      lon2d, lat2d = np.meshgrid(lon1d, lat1d)      
    gridX,gridY, mRef,vals_loc = regrid_map.driver(lat_windMax, lon_windMax, lat2d, lon2d, vals)
    
    vals_allCases.append(vals_loc)

  #return (gridX-x0, gridY-y0, vals_allCases)  
  return (gridX, gridY, vals_allCases)


def find_mslp_DT_min(t0, latRef, lonRef, rSfc, rDT, rEarth, calcMonthlyAnomaly=True):
  #return the locations and values of MSLP and DT PT minima w/in rX of reference point
  
  #get mslp and DT PT fields
  fAtmo = helpers.get_atmo_filename_year(t0, info='mslp')
  data = netCDF4.Dataset(fAtmo,'r')
  fileTimes = data.variables['time'][:]
  iTime = helpers.eraI_time2TimeInd(fileTimes,t0)
  vals_sfc = data.variables['msl'][iTime,:,:]
  lat_sfc = data.variables['latitude'][:]
  lon_sfc = data.variables['longitude'][:]
  data.close()
  
  fAtmo = helpers.get_atmo_filename_year(t0, info='2pvu_theta_wind')
  data = netCDF4.Dataset(fAtmo,'r')
  fileTimes = data.variables['time'][:]
  iTime = helpers.eraI_time2TimeInd(fileTimes,t0)
  vals_DT = data.variables['pt'][iTime,:,:]
  lat_DT = data.variables['latitude'][:]
  lon_DT = data.variables['longitude'][:]
  data.close()
  
  if (calcMonthlyAnomaly): #subtract off 1979-2014 mean for that month
      valsRef, std, latclimo, lonclimo  = helpers.get_monthlyMean_andStd_ERA5(key='msl',scaleFac=1.0e-2,iMonth=t0.month-1)
      [minLat] = np.where(latclimo == np.min(lat_sfc))
      step = int(len(lonclimo)/len(lon_sfc))
      vals_sfc -= valsRef[minLat[0]::step,::step]
      
      valsRef, std, latclimo, lonclimo = helpers.get_monthlyMean_andStd_ERA5(key='pt',scaleFac=1.0,iMonth=t0.month-1)
      [minLat] = np.where(latclimo == np.min(lat_DT))
      step = int(len(lonclimo)/len(lon_DT))
      vals_DT -= valsRef[minLat[0]::step, ::step]
  
  #find MSLP min near ice point
  lon2d,lat2d = np.meshgrid(lon_sfc,lat_sfc)
  d2r = np.pi/180.
  distMesh = helpers.calc_distSphere_multiple(rEarth, latRef*d2r, lonRef*d2r, lat2d*d2r, lon2d*d2r)
  iMin = np.argmin(vals_sfc[distMesh<rSfc])
  min_sfc_val = vals_sfc[distMesh<rSfc][iMin]
  min_sfc_lat = lat2d[distMesh<rSfc][iMin]
  min_sfc_lon = lon2d[distMesh<rSfc][iMin]
  
  #find DT PT min near MSLP point
  lon2d,lat2d = np.meshgrid(lon_DT,lat_DT)
  distMesh = helpers.calc_distSphere_multiple(rEarth, min_sfc_lat*d2r, min_sfc_lon*d2r, lat2d*d2r, lon2d*d2r)
  iMin = np.argmin(vals_DT[distMesh<rDT])
  min_DT_val = vals_DT[distMesh<rDT][iMin]
  min_DT_lat = lat2d[distMesh<rDT][iMin]
  min_DT_lon = lon2d[distMesh<rDT][iMin]
  
  dLat = min_sfc_lat-min_DT_lat; dLat *= d2r
  dLon = min_sfc_lon-min_DT_lon; dLon *= d2r
  latRef = .5*(min_DT_lat+min_sfc_lat); u = np.cos(latRef*d2r)*dLon; v = dLat
  #I don't think we need to trap u~0 since atan2 handles +/- inf appropriately.
  angle = np.arctan2(v,u)

  if (True):
    print('min,lat,lon for MSLP, DT: ', t0, min_sfc_val, min_sfc_lat, min_sfc_lon, min_DT_val, min_DT_lat, min_DT_lon, angle/d2r)
  
  return (min_sfc_val, min_sfc_lat, min_sfc_lon, min_DT_val, min_DT_lat, min_DT_lon, angle)


def calc_eventLocation_ctrMass(objFile,VRILE):
    data = np.load(objFile)
    lat = data['lat']; lon = data['lon']; obj = data['vals']; mask = data['mask']
    
    infoCases = []
    for itr in range(len(VRILE)):
        tmp_obj = obj[itr,:,:] * mask[itr,:,:]
        s = np.shape(tmp_obj)
        X_1d = np.reshape(lon,s[0]*s[1]); Y_1d = np.reshape(lat,s[0]*s[1]); obj_1d = np.reshape(np.abs(tmp_obj),s[0]*s[1])
        tmp_land = globe.is_land(Y_1d,X_1d)
        [ind] = np.where(tmp_land == True)
        obj_1d[ind] = 0
        M = np.sum(obj_1d)
        x_loc = (np.sum(X_1d * obj_1d))/M
        y_loc = (np.sum(Y_1d * obj_1d))/M
        # if x_loc > 180:
        #     x_loc = x_loc - 360
        infoCases.append([VRILE[itr][0], y_loc, x_loc])
    fEventInfo = '/raid3/maddie/Data/djf_bwfilter_3d_5percentile_dtUnique01_massLoc_1979-2023_noshift.pkl'
    output = open(fEventInfo, 'wb')
    pickle.dump(infoCases,output)
    output.close()


def calc_loc_maxWind_iceLossObject(u,v,mask):
  #Input u,v[times,lats,lons]. mask[lats,lons].
  #return iTime,iLat,iLon, and wind angle (math angle wrt x-axis) for max wind over valid mask
  
  speed = np.sqrt(u*u+v*v)
  speed = speed*mask[None,:,:]
  ind = np.argmax(speed)
  iTime,iLat,iLon = np.unravel_index(ind, speed.shape)
  
  angle = np.arctan2(v[iTime,iLat,iLon],u[iTime,iLat,iLon])
  
  return (iTime,iLat,iLon,angle)


def get_mean(data,ax):
  data = np.asarray(data)
  data = np.nanmean(data,ax)
  return data

def get_std(data,ax):
  data = np.asarray(data)
  data = np.std(data,ax)
  return data

def make_event_locations(N,minLat=55,minLon=-180,step=0.5,yrMin=1979,yrMax=2016,season='winter'):
  # This function makes fake vRILE locations by selecting a random place (over ocean or ice) at a random time (in specified
  # date range)
  Cases = []
  #Cases = np.empty((N,3))
  # this is for the land mask
  #m = Basemap(projection='npstere',boundinglat=minLat,lon_0=0,resolution='l')
 
  year_list = range(yrMin,yrMax + 1)
  if season == 'summer':
    day_list = list(range(1,31)) + list(range(1,32)) + list(range(1,32))
    mth_list = [6] * 30 + [7] * 31 + [8] * 31
  elif season == 'winter':
    day_list = list(range(1,32)) + list(range(1,32)) + list(range(1,29))
    mth_list = [12] * 31 + [1] * 31 + [2] * 28

  else:
    print('Which season?')

  lats = np.arange(minLat,90 + step,step)
  if minLon < 0:
    lons = np.arange(minLon,180 + step,step)
  else:
    lons = np.arange(0,360 + step,step)
  
  for itr in range(N):
    land = True
    # This bit just makes sure that the random point we pick is over ocean or ice and not land.
    while land == True:
      latId = np.random.randint(0,len(lats))
      lonId = np.random.randint(0,len(lons))
      tmpLat = lats[latId]; tmpLon = lons[lonId]
      #x,y = m(tmpLon,tmpLat)
      land = is_land(tmpLon,tmpLat)

    yearId = np.random.randint(0,len(year_list))
    dayId = np.random.randint(0,len(day_list))

    tmpDate = dt.datetime(year_list[yearId],mth_list[dayId],day_list[dayId])
    tmpLat = lats[latId]; tmpLon = lons[lonId]

    Cases.append([tmpDate, tmpLat, tmpLon])
    #Cases[itr,0] = tmpDate; #Cases[itr,1] = tmpLat; Cases[itr,2] = tmpLon

  # And now we save the random dates and locations we made
  fEventInfo = '/raid3/maddie/Data/' + season + '_random_location.pkl'
  output = open(fEventInfo, 'wb')
  pickle.dump(Cases,output)
  output.close()
  

def plot_mean(gridX,gridY,vals_mslp,vals_pt):
  mean_mslp = get_mean(vals_mslp,ax=0)#/100.
  #vals_pt = np.ma.masked_values(vals_pt, -32767.)
  #mean_pt = np.ma.mean(vals_pt,axis=0)
  #mean_pt = get_mean(vals_pt,ax=0)
  #print(np.min(mean_pt), np.max(mean_pt), np.min(mean_mslp), np.max(mean_mslp))
  #print(np.min(vals_mslp[i,:,:]), np.max(vals_mslp[i,:,:]), np.min(
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #levs = np.arange(0.15,1,0.1)
  levs = np.arange(-3,3.5,0.5)
  clevs = np.arange(-3,3.1,.1)
  #levs = np.arange(-16,16+1,2)
  #clevs = np.arange(-16,16.1,.1)
  
#  plt.contourf(gridX,gridY,mean_pt,levels=clevs,cmap='YlGnBu',extend='both')
  plt.contourf(gridX,gridY,mean_mslp,cmap='bwr',levels=clevs, extend='both')
  cbar = plt.colorbar()
#  cbar.set_label('DT PT Anomaly ($\sigma$)')

  Z = ndimage.gaussian_filter(mean_mslp,sigma=2,order=0) 
  cs = plt.contour(gridX,gridY,Z,colors='k',levels=levs)
  plt.clabel(cs,levels=[-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3],inline=True, fontsize=10)

  plt.xlabel('Distance (3000 km)')
  plt.ylabel('Distance (3000 km)')
  #plt.title('5 Day Lag ( Winter vRILEs 1979 - 2016)')
  ax.set_xticks([]); ax.set_yticks([])
  plt.savefig('/raid3/maddie/test.png')

def plot_sig(gridX,gridY,vals_mslp,vals_pt,BS,cutoff):
  # this shows significance of the mslp signal
  mean_mslp = get_mean(vals_mslp, ax=0)
  mean_pt = get_mean(vals_pt, ax=0)  
  r,c = bootstrap.is_sig(BS,mean_mslp,cutoff)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  clevs = np.arange(-3,3+0.1,0.1)
  levs = np.arange(-3,3.5,0.5)

  plt.contourf(gridX,gridY,mean_pt,levels=clevs,cmap='bwr',extend='both')
  cbar = plt.colorbar()
  cbar.set_label('DT PT Anomaly ($\sigma$)')
  
  cs = plt.contour(gridX,gridY,mean_mslp,colors='k',levels=levs)
  plt.clabel(cs,levels=[-3,-2,-1,0,1,2,3],inline=True, fontsize=10)

  plt.plot(gridX[r,c],gridY[r,c],'k.',markersize = 1,alpha=0.5,label=str((1-cutoff)*100) + '% Significance')
  ax.legend()

  plt.xlabel('Distance (3000 km)')
  plt.ylabel('Distance (3000 km)')
  ax.set_xticks([]); ax.set_yticks([])
  savename = '/raid3/maddie/test.png'
  plt.savefig(savename) 


def steven_dat_to_datetime_VRILE(dataDir,fname):
  tmp = np.loadtxt(dataDir + fname)
  # this should be ## x 3, where the 1st col is an int or float YYYYMMDDHH
  dateStr = [str(int(tmp[x,0])) for x in range(np.shape(tmp)[0])]
  V = []
  for itr in range(tmp.shape[0]):
    d = dateStr[itr]
    t0 = dt.datetime.strptime(d,'%Y%m%d%H')
    V.append([t0,tmp[itr,2]])
    
    if (False):
      # check and see if we have the ice concentrations files
      fStart,fEnd = helpers.get_ice_filenames(t0)
      if os.path.isfile(fStart) == False:
        print('oh no! no file')
        helpers.download_iceConcentration_forChange(t0)
  return V


def get_seaIce_objects(VRILE):
    obj = []; mask = []
    
    for itr in range(len(VRILE)):
        V = VRILE[itr]
        print(itr, V[0])
        fIceStart,fIceEnd = helpers.get_ice_filenames(V[0])
        if (os.path.isfile(fIceStart) == False):
          print('Missing: ', fIceStart)
        if (os.path.isfile(fIceEnd) == False):
          print('Missing: ', fIceEnd) 
        #if (os.path.isfile(fIceStart) == False):    
        #    one_day = dt.timedelta(days=1)
        #    fIceStart,fIceEnd = helpers.get_ice_filenames(V[0] - one_day)
        #    count = 0
        #    while ((os.path.isfile(fIceStart) == False) | (os.path.isfile(fIceEnd) == False)) & (count < 2):
         #       fIceStart,fIceEnd = helpers.get_ice_filenames(V[0] + one_day * count)
          #      count += 1
          #  if os.path.isfile(fIceStart) == False:
          #      print('oh no! no file')
                #helpers.download_iceConcentration_forChange(V[0])    
        print(fIceStart, fIceEnd)
        lat,lon,iceLoss = seaIceMasks.calc_iceLossMask_refPoint(fIceStart, fIceEnd,returnMask=True,dxiceThresh=-.1)
        if np.min(iceLoss) <= -999.:
          print('Bad Data: ', V[0])
          continue
        mask.append(iceLoss)
    
        latCon,lonCon,xice,icekey = helpers.get_iceConc(fIceEnd,iceKey='nsidc_nt_seaice_conc',canSwitch=False)
        obj.append(xice)
    
    fname = '/raid3/maddie/Data/jja_bwfilter_3d_5percentile_dtUnique01_objLoc_1979-2023_noshift.npz'
    np.savez(fname, lat=lat, lon=lon, vals=obj, mask=mask)

 
f = open('/raid3/maddie/Data/djf_bwfilter_3d_5percentile_dtUnique01_massLoc_1979-2023_noshift.pkl','rb')
infoCases = pickle.load(f)
f.close()


#get_seaIce_objects(V)

#objFile = '/raid3/maddie/Data/djf_bwfilter_3d_5percentile_dtUnique01_objLoc_1979-2023_noshift.npz'
#calc_eventLocation_ctrMass(objFile,V)


#infoCases = vrile_data(dataDir='/raid3/datasets/seaIce/VRILE/', iceFile='rapid_seaice_loss_events_int_withboth_jja_bwfilter_3d_5percentile.dat')

#refFrame is either "sfcWind" or "min_DT-MSLP"
gridX, gridY, vals_mslp = driver_values_localRef_cases(infoCases, dtLag=3, refFrame='min_DT-MSLP',isAtmo=True, key='msl', infoKey='mslp', rCutoff_ice_sfcCyclone=1000., rCutoff_tpv_sfcCyclone=1000., rEarth=6370.,calcMonthlyAnomaly=True)

gridX, gridY, vals_pt = driver_values_localRef_cases(infoCases, dtLag=3, refFrame='min_DT-MSLP',isAtmo=True, key='pt', infoKey='2pvu_theta_wind',rCutoff_ice_sfcCyclone=1000., rCutoff_tpv_sfcCyclone=1000., rEarth=6370.,calcMonthlyAnomaly=True)


#plot_mean(gridX,gridY,vals_mslp,vals_pt)
if (True):
  fCaseVals_slp = '/raid3/maddie/Data/djf_mslpCtr_massLoc_mslp_stdAnom_1979-2023_Lag3.npz'
  fCaseVals_dt = '/raid3/maddie/Data/djf_mslpCtr_massLoc_theta_stdAnom_1979-2023_Lag3.npz'
  np.savez(fCaseVals_slp, gridX=gridX, gridY=gridY, vals=vals_mslp)
  np.savez(fCaseVals_dt, gridX=gridX, gridY=gridY, vals=vals_pt)

#MCaseVals = '/raid4/maddie/Data_CESM2/djf_mslpCtr_massLoc_mslp_cesm_1960-2000.npz'
#DCaseVals = '/raid3/maddie/Data/jja_bwfilter_mslpCtr_massLoc_theta_stdAnom_1979-2022.npz'
#MCaseVals = '/raid3/maddie/Data/jja_mslpCtr_massLoc_mslp_stdAnom_1979-2023.npz'
#DCaseVals = '/raid3/maddie/Data/jja_mslpCtr_massLoc_theta_stdAnom_1979-2023.npz'
#DCaseVals = '/raid3/maddie/Data/jja_iceCtr_massLoc_ke_1988-2019.npz'

#data_BS = np.load('/raid3/maddie/summer_iceCtr_stdAnom_1979-2022_mslp.npz'); BS = data_BS['BS'][:,:,:]
#data_msl = np.load(MCaseVals); vals_mslp = data_msl['vals'][:,:,:]; gridX = data_msl['gridX'][:]; gridY = data_msl['gridY'][:]
#data_pt = np.load(DCaseVals); vals_pt = data_pt['vals'][:,:,:]

#plot_mean(gridX,gridY,vals_mslp,vals_pt)

#make_event_locations(N=5000,minLat=55,minLon=-180,step=0.5,yrMin=1979,yrMax=2021,season='summer')

  
  
