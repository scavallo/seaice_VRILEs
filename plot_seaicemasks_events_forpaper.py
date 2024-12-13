# plot_seaicemasks_events_forpaper
#
# makes sea ice location masks for VRILE data
#
# Steven Cavallo
# December 2024
###########################################
import numpy as np
import netCDF4
import datetime as dt
import matplotlib.pyplot as plt
#import cPickle as pickle
import pickle
from scipy import ndimage
from mpl_toolkits.basemap import Basemap, addcyclic
import os

from mstats import *
import utilities_modules as um
#import helpers_v2
import seaIceMasks

ftp_seaiceconc = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02202_V3/north/daily/'
#fDir_seaiceconc = '/arctic3/datasets/seaIce/concentration/'
fDir_seaiceconc = '/data1/scavallo/data/seaice/concentration/daily/'

ftp_seaiceage = 'ftp://sidads.colorado.edu/DATASETS/nsidc0611_seaice_age_v4/data/'
fDir_seaiceage = '/data1/scavallo/data/seaice/ice_age/'


fEventInfo_dir = '/data1/scavallo/data/seaice/VRILEs/1979_2023/'
fEventInfo_name = 'rapid_seaice_loss_events_int_withboth_nonjja_bwfilter_3d_10percentile_dtUnique01_noshift.dat'
fEventInfo_outfile = 'VRILE_masks_nonjja_bwfilter_3d_5percentile_dtUnique01_noshift_2014_2023.dat'
file_numpy =  '/data1/scavallo/data/seaice/VRILEs/1979_2023/VRILE_masks_nonjja_bwfilter_3d_10percentile_dtUnique01_noshift_2014_2023.npz'
year_range = ['2014010100','2023123118']

imagedir = '/home/scavallo/scripts/python_scripts/images/'


map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
proj_latlon = [55. , 270.]
zoom = 'false'
cen_lon = 270.
label_fontsize = 22

#################################################
#year2plot = 2017
outfile = open(fEventInfo_dir+fEventInfo_outfile,'w')
outfile.write('date lat lon iceloss')
outfile.write('\n')

#fpath_iceage = fDir_seaiceage + 'iceage_nh_12.5km_2017.nc'
fpath_iceage = fDir_seaiceage + 'seaiceage_august_1984_2018.nc'
f = netCDF4.Dataset(fpath_iceage,'r')
lons_iceage = f.variables['lon'][:]
lats_iceage = f.variables['lat'][:]
years_iceage = f.variables['years'][:]
#[yind] = np.where(years_iceage == year2plot)
print(years_iceage)
iceage = np.nanmedian(f.variables['age_of_sea_ice'][:,:,:],0).squeeze()
iceage_min = np.nanmedian(f.variables['age_of_sea_ice'][-10:,:,:],0).squeeze()
iceage_max = np.nanmedian(f.variables['age_of_sea_ice'][0:10,:,:],0).squeeze()
f.close


iceage = um.filter_numeric_nans(iceage,10,float('NaN'),'high')
iceage_min = um.filter_numeric_nans(iceage_min,10,float('NaN'),'high')
iceage_max = um.filter_numeric_nans(iceage_max,10,float('NaN'),'high')
mstats(iceage)


if 1 == 1:
    fEventInfo = fEventInfo_dir + fEventInfo_name
    bb = np.loadtxt(fEventInfo, skiprows=0)       
    infoCases = bb[:,0].astype(int)
    Ntimeranges = len(infoCases)
    print(infoCases)
    
    os.chdir(fDir_seaiceconc)
    
    iceVals = 0
    niters = 2

    lat_app = []
    lon_app = []
    ice_app = []
    
    lat_pt = []
    lon_pt = []
    mask_pt = []   
    date_event_pt = [] 
    for t0 in range(0,Ntimeranges):        			
        datenow = infoCases[t0]
        #print(t0,datenow)    
        date_record_str = str(datenow)	
        for tt in range(0,niters):	
            #print('Event data:' , date_record_str)

            yyyy_event = date_record_str[0:4]
            mm_event = date_record_str[4:6]
            dd_event = date_record_str[6:8]
            hh_event = date_record_str[8:10]
            yyyymmdd_event = date_record_str[0:8]

            if int(yyyy_event) < 1987:
               ndays = 6
            else:
               ndays = 5

            if int(yyyy_event) < 1987:
                fname_in = 'seaice_conc_daily_nh_n07_' + yyyymmdd_event + '_v03r01.nc'
            elif ( (int(yyyy_event) >=1987) & (int(yyyy_event) <= 1991)):
                fname_in = 'seaice_conc_daily_nh_f08_' + yyyymmdd_event + '_v03r01.nc'  		    
            elif ( (int(yyyy_event) >=1992) & (int(yyyy_event) <= 1995)):
                fname_in = 'seaice_conc_daily_nh_f11_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >=1996) & (int(yyyy_event) <= 2007)): 
                fname_in = 'seaice_conc_daily_nh_f13_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >=2008) & (int(yyyy_event) <= 2017)):
                fname_in = 'seaice_conc_daily_nh_f17_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >= 2018) ):
                fname_in = 'seaice_conc_daily_icdr_nh_f18_'+ yyyymmdd_event + '_v01r00.nc'	

            fpath_in = fDir_seaiceconc +  fname_in
            #print(fpath_in)
	    
            fcheck = os.path.isfile(fpath_in)
            if fcheck == True :    
                #print "File exists"
                tt = niters+1
            else:	    		
                #print "Checking one day back"	   
                date_record_str = um.advance_time(date_record_str,-24) 	    
	
        date_event = date_record_str
        fpath_event = fpath_in
		
        date_record_str = um.advance_time(date_record_str,-ndays*24) 	
        for tt in range(0,niters):	
	    #print('Reference data:', date_record_str)

            yyyy_event = date_record_str[0:4]
            mm_event = date_record_str[4:6]
            dd_event = date_record_str[6:8]
            hh_event = date_record_str[8:10]
            yyyymmdd_event = date_record_str[0:8]

            if int(yyyy_event) < 1987:
                ndays = 6
            else:
                ndays = 5

            if int(yyyy_event) < 1987:
                fname_in = 'seaice_conc_daily_nh_n07_' + yyyymmdd_event + '_v03r01.nc'
            elif ( (int(yyyy_event) >=1987) & (int(yyyy_event) <= 1991)):
                fname_in = 'seaice_conc_daily_nh_f08_' + yyyymmdd_event + '_v03r01.nc'  		    
            elif ( (int(yyyy_event) >=1992) & (int(yyyy_event) <= 1995)):
                fname_in = 'seaice_conc_daily_nh_f11_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >=1996) & (int(yyyy_event) <= 2007)): 
                fname_in = 'seaice_conc_daily_nh_f13_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >=2008) & (int(yyyy_event) <= 2017) ):
                fname_in = 'seaice_conc_daily_nh_f17_' + yyyymmdd_event + '_v03r01.nc'		    
            elif ( (int(yyyy_event) >= 2018) ):
                fname_in = 'seaice_conc_daily_icdr_nh_f18_'+ yyyymmdd_event + '_v01r00.nc'	
	

            fpath_in = fDir_seaiceconc +  fname_in
	    #print(fpath_in)	    

            fcheck = os.path.isfile(fpath_in) 
            if fcheck == True :    
		#print "File exists"
                tt = niters+1
            else:	    
                #print "Checking one day back"	   
                date_record_str = um.advance_time(date_record_str,-24) 	    
		        

            date_ref = date_record_str
            fpath_ref = fpath_in
        
        print(date_event, date_ref)	    
        fIceStart = fpath_ref 
        fIceEnd = fpath_event
	
        try:
            lat_ptnow,lon_ptnow, mask_ptnow = seaIceMasks.calc_iceLossMask_refPoint(fIceStart, fIceEnd, returnMask=False, dxiceThresh=-.1)
            lat_pt.append(lat_ptnow)
            lon_pt.append(lon_ptnow)
            mask_pt.append(mask_ptnow)
            date_event_pt.append(date_event)	
        except:
            print('missed', fIceStart,fIceEnd)
	    
        print(int(date_event),int(year_range[0]),int(year_range[1]))			
        if ( (int(date_event) >= int(year_range[0]) ) and  (int(date_event) <= int(year_range[1])) ):
	#if ( int(date_event) > 2009123118 ):
	#if ( int(date_event) >= 2008123118 ):
	#if ( int(date_event) <= 2018123118 ):
            print(int(date_event),int(year_range[0]),int(year_range[1]))	
            try:
                print(fIceStart,fIceEnd)	        
                lat,lon, mask = seaIceMasks.calc_iceLossMask_refPoint(fIceStart, fIceEnd,returnMask=True, dxiceThresh=-.1)		
	        
                iceVals += mask	    
		
		#print(lat,lon)
                mstats(iceVals)
                lat_app.append(lat)
                lon_app.append(lon)
                ice_app.append(iceVals)	
		
                latv = lat
                lonv = lon	    		
							
            except:
                print('Missed')
    # lon,lat,iceVals    
np.savez(file_numpy,iceVals=iceVals,icelats=lat,icelons=lon)         		

lat_pt = np.array(lat_pt)
lon_pt = np.array(lon_pt)
mask_pt = np.array(mask_pt)
date_event_pt = np.array(date_event_pt)
for ii in range(0,len(lat_pt)):
    if ( (int(date_event_pt[ii]) >= int(year_range[0]) ) and  (int(date_event_pt[ii]) <= int(year_range[1])) ):
        outfile.write('%-10s %7.2f %7.2f %7.5f\n' %(date_event_pt[ii],lat_pt[ii],lon_pt[ii],mask_pt[ii]))    

