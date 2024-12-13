# make_era5_composites_fromevents_forpaper
#
# Gets ERA5 data corresponding to VRILE event files
#
# Steven Cavallo
# December 2024
###########################################
# imports
import os

import numpy as np
import netCDF4
import datetime as dt
import matplotlib.pyplot as plt
from scipy import ndimage
import os
from netCDF4 import Dataset

from mstats import *
import utilities_modules as um
import helpers_v2


#########################################
fEventInfo_dir = '/raid3/datasets/VRILEs/1979_2023/'
fEventInfo_name = 'rapid_seaice_loss_events_int_withboth_jja_bwfilter_3d_10percentile_dtUnique01_noshift.dat'

fEvent_out_dir = '/raid3/scavallo/data/VRILEs/1979_2023/'
fEvent_netcdf_out = 'ERA5_AtmosComposites_jja_bwfilter_3d_10percentile_dtUnique01_noshift_1989_2023.nc'

fpath_atmos = '/raid1/ERA5/0p5/mslp/'
fpath_trop = '/raid1/ERA5/0p5/2pvu_theta_wind/'
fpath_atmos_climo = '/raid1/ERA5/climo/0p5/mslp.nc'
fpath_trop_climo = '/raid1/ERA5/climo/0p5/2pvu_theta.nc'

atmos_info = 'sfc'
atmos_dt_info = 'DT'
nDaysOffset = 0
#year_range_reanal = [1989,1998]
#year_range_reanal = [2013,2022]
#year_range_reanal = [1989,2022]

#year_range_reanal = [1989,1998]
#year_range_reanal = [2014,2023]
year_range_reanal = [1989,2023]
#########################################

data_trop_climo = netCDF4.Dataset(fpath_trop_climo,'r')
pt_climo_times = data_trop_climo.variables['time'][:]

data_atmos_climo = netCDF4.Dataset(fpath_atmos_climo,'r')
slp_climo_times = data_atmos_climo.variables['time'][:]

fEventInfo = fEventInfo_dir + fEventInfo_name
bb = np.loadtxt(fEventInfo, skiprows=1)       
infoCases = bb[:,0].astype(int)
Ntimeranges = len(infoCases)
print(infoCases)

yyyy_climo_ref = '2020'

nsamps_sfc = 0
nsamps_DT = 0
iternum = 0
for t0 in range(0,Ntimeranges):        			
    datenow = infoCases[t0]    
    date_record_str = str(datenow)	

    yyyy_event = date_record_str[0:4]
    mm_event = date_record_str[4:6]
    dd_event = date_record_str[6:8]
    hh_event = date_record_str[8:10]
    yyyymmdd_event = date_record_str[0:8]
    yyyymmddhh_event = date_record_str[0:10]

    if ( (int(yyyy_event)< year_range_reanal[0])):
        continue

    if ( (int(yyyy_event) > year_range_reanal[1]) ):
        continue
    
    
    
    fnow = yyyy_event + '.nc'
    fPath = fpath_atmos + fnow
    fPath_DT = fpath_trop + fnow
	
    if (True):
        print('Event date %s and Surface file %s' %(date_record_str, fPath))    
        print('Tropopause file %s' %(fPath_DT))
        
    tnow = dt.datetime(int(yyyy_event),int(mm_event),int(dd_event),int(hh_event))    
    tnow_climo = dt.datetime(int(yyyy_climo_ref),int(mm_event),int(dd_event),int(hh_event))  
    timeindex_climo = um.date_to_dayofyear(yyyy_event,mm_event,dd_event)-1
    print(timeindex_climo)

    tRef = dt.datetime(int(yyyy_climo_ref),1,1,0)    
    timeVal = (tnow_climo-tRef).total_seconds()/3600.
    timeVal = int(timeVal)  
    iTime0 =  np.argmin(np.absolute(timeVal-pt_climo_times))
    stimeindex = iTime0-2
    etimeindex = iTime0+2 
    
    pt_climo = np.nanmean(data_trop_climo.variables['pt'][stimeindex:etimeindex,:,:],0).squeeze()       
    mslp_climo = (np.nanmean(data_atmos_climo.variables['msl'][stimeindex:etimeindex,:,:],0).squeeze())/100. 
    
    tRef = dt.datetime(1900,1,1,0)    
    timeVal = (tnow-tRef).total_seconds()/3600.
    timeVal = int(timeVal)      

    try:
        data_trop = netCDF4.Dataset(fPath_DT,'r') 
        pt_times = data_trop.variables['time'][:]        
        lats_trop = data_trop.variables['latitude'][:]
        lons_trop = data_trop.variables['longitude'][:]
        iTime0 =  np.argmin(np.absolute(timeVal-pt_times)) 
        stimeindex_trop = iTime0-2
        etimeindex_trop = iTime0+2 
        trth_now = np.nanmean(data_trop.variables['pt'][stimeindex_trop:etimeindex_trop,:,:],0).squeeze() 
        nsamps_DT += 1 	
        surfaceflag = True
    except:
        print('Could not find tropopause file %s' %(fPath_DT))
        surfaceflag = False     
   
    try:
        data_atmos = netCDF4.Dataset(fPath,'r') 
        atmos_times = data_atmos.variables['time'][:]        
        lats_atmos = data_atmos.variables['latitude'][:]
        lons_atmos = data_atmos.variables['longitude'][:]
        iTime0 =  np.argmin(np.absolute(timeVal-atmos_times)) 
        stimeindex_atmos = iTime0-2
        etimeindex_atmos = iTime0+2 
        mslp_now = (np.nanmean(data_atmos.variables['msl'][stimeindex_atmos:etimeindex_atmos,:,:],0).squeeze())/100.	
        nsamps_sfc += 1
        tropflag = True
    except:
        print('Could not find surface file %s' %(fPath))
        tropflag = False      


    mslp_anom = mslp_now - mslp_climo
    trth_anom = trth_now - pt_climo

    mslp_anom_std = (mslp_now - mslp_climo) / np.nanstd(mslp_climo)
    trth_anom_std = (trth_now - pt_climo) / np.nanstd(pt_climo)
      
    if iternum == 0:
        if surfaceflag ==  True: 
            mslp_sum = mslp_now
            mslp_anom_sum = mslp_anom
            mslp_anom_std_sum = mslp_anom_std
	
        if tropflag == True:
            trth_sum = trth_now
            trth_anom_sum = trth_anom
            trth_anom_std_sum = trth_anom_std
    else:
        if surfaceflag ==  True: 
            mslp_sum = mslp_sum + mslp_now	
            mslp_anom_sum = mslp_anom_sum + mslp_anom
            mslp_anom_std_sum = mslp_anom_std_sum + mslp_anom_std
	
        if tropflag == True:
            trth_sum = trth_sum + trth_now
            trth_anom_sum = trth_anom_sum + trth_anom
            trth_anom_std_sum = trth_anom_std_sum + trth_anom_std

    
    iternum += 1

data_atmos.close()
data_trop.close()    
data_trop_climo.close()
data_atmos_climo.close()

print(nsamps_sfc, nsamps_DT)
mslp_comp = mslp_sum / nsamps_sfc    
#u10_comp = u10_sum / nsamps_sfc    
#v10_comp = v10_sum / nsamps_sfc    

print("mslp_comp:")
mstats(mslp_comp)

trth_comp = trth_sum / nsamps_DT    
#tru_comp = tru_sum / nsamps_DT    
#trv_comp = trv_sum / nsamps_DT    

mslp_anom_comp = mslp_anom_sum / nsamps_sfc    
trth_anom_comp = trth_anom_sum / nsamps_DT    
#tru_anom_comp = tru_anom_sum / nsamps_DT    
#trv_anom_comp = trv_anom_sum / nsamps_DT    

mslp_anom_std_comp = mslp_anom_std_sum / nsamps_sfc    
trth_anom_std_comp = trth_anom_std_sum / nsamps_DT    
#tru_anom_std_comp = tru_anom_std_sum / nsamps_DT    
#trv_anom_std_comp = trv_anom_std_sum / nsamps_DT  


nlats_sfc = len(lats_atmos)
nlons_sfc = len(lons_atmos)
nlats_DT = len(lats_trop)
nlons_DT = len(lons_trop)

nlats = nlats_sfc
nlons = nlons_sfc

nsample_arr_sfc = np.zeros_like(mslp_comp)
nsample_arr_trop = np.zeros_like(trth_comp)
nsample_arr_sfc[:,:] = nsamps_sfc
nsample_arr_trop[:,:] = nsamps_DT


fout = fEvent_out_dir + fEvent_netcdf_out
print(fout)
ncfile = Dataset(fout,'w',format='NETCDF4') 
ncfile.createDimension('latitude',nlats)
ncfile.createDimension('longitude',nlons)
ncfile.createDimension('latitude_dt',nlats_DT)
ncfile.createDimension('longitude_dt',nlons_DT)
lats_out = ncfile.createVariable('latitude','float32',('latitude'))
lons_out = ncfile.createVariable('longitude','float32',('longitude'))
lats_out_dt = ncfile.createVariable('latitude_dt','float32',('latitude_dt'))
lons_out_dt = ncfile.createVariable('longitude_dt','float32',('longitude_dt'))
mslp_out = ncfile.createVariable('mslp_comp','float32',('latitude','longitude'))
#u10_out = ncfile.createVariable('u10_comp','float32',('latitude','longitude'))
#v10_out = ncfile.createVariable('v10_comp','float32',('latitude','longitude'))
trth_out = ncfile.createVariable('trth_comp','float32',('latitude_dt','longitude_dt'))
#tru_out = ncfile.createVariable('tru_comp','float32',('latitude','longitude'))
#trv_out = ncfile.createVariable('trv_comp','float32',('latitude','longitude'))
mslp_out_anom = ncfile.createVariable('mslp_anom_comp','float32',('latitude','longitude'))
trth_out_anom = ncfile.createVariable('trth_anom_comp','float32',('latitude_dt','longitude_dt'))
#tru_out_anom = ncfile.createVariable('tru_anom_comp','float32',('latitude','longitude'))
#trv_out_anom = ncfile.createVariable('trv_anom_comp','float32',('latitude','longitude'))
mslp_out_anom_std = ncfile.createVariable('mslp_anom_std_comp','float32',('latitude','longitude'))
trth_out_anom_std = ncfile.createVariable('trth_anom_std_comp','float32',('latitude_dt','longitude_dt'))
#tru_out_anom_std = ncfile.createVariable('tru_anom_std_comp','float32',('latitude','longitude'))
#trv_out_anom_std = ncfile.createVariable('trv_anom_std_comp','float32',('latitude','longitude'))
nsamples_sfc_out = ncfile.createVariable('nsamples_sfc','float32',('latitude','longitude'))
nsamples_trop_out = ncfile.createVariable('nsamples_trop','float32',('latitude_dt','longitude_dt'))
lats_out.units = 'degrees_north'
lons_out.units = 'degrees_east'
lats_out_dt.units = 'degrees_north'
lons_out_dt.units = 'degrees_east'
mslp_out.units = 'hPa'
#u10_out.units = 'm s-1'
#v10_out.units = 'm s-1'
trth_out.units = 'Kelvin'
#tru_out.units = 'm s-1'
#trv_out.units = 'm s-1'
mslp_out_anom.units = 'hPa'
trth_out_anom.units = 'Kelvin'
#tru_out_anom.units = 'm s-1'
#trv_out_anom.units = 'm s-1'
mslp_out_anom_std.units = 'Standard deviations'
trth_out_anom_std.units = 'Standard deviations'
#tru_out_anom_std.units = 'Standard deviations'
#trv_out_anom_std.units = 'Standard deviations'
lats_out[:] = lats_atmos
lons_out[:] = lons_atmos
lats_out_dt[:] = lats_trop
lons_out_dt[:] = lons_trop
mslp_out[:] = mslp_comp
#u10_out[:] = u10_comp
#v10_out[:] = v10_comp
trth_out[:] = trth_comp
#tru_out[:] = tru_comp
#trv_out[:] = trv_comp
mslp_out_anom[:] = mslp_anom_comp
trth_out_anom[:] = trth_anom_comp
#tru_out_anom[:] = tru_anom_comp
#trv_out_anom[:] = trv_anom_comp
mslp_out_anom_std[:] = mslp_anom_std_comp
trth_out_anom_std[:] = trth_anom_std_comp
#tru_out_anom_std[:] = tru_anom_std_comp
#trv_out_anom_std[:] = trv_anom_std_comp
nsamples_sfc_out[:] = nsample_arr_sfc
nsamples_trop_out[:] = nsample_arr_trop
#ncfile = Dataset(outfilename,'w') 
ncfile.close()

