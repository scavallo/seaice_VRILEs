# plot_era5_composites_forpaper
#
# Plots VRILE locations/masks, SLP, and ice age composites from ERA5 data
#
# Steven Cavallo
# December 2024
###########################################
# imports
import netCDF4

import os, datetime, pylab,sys
import numpy as np
import matplotlib as mpl
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy import ndimage
from matplotlib.colors import LogNorm
import numpy.ma as ma

# Add a couple of user defined functions
import weather_modules as wm
import utilities_modules as um
from mstats import *


import warnings
warnings.filterwarnings("ignore")
###################################
# Set user options
###################################
print("Python version")
print(sys.version)


fdir = '/data1/scavallo/data/seaice/VRILEs/1979_2023/'
fname = 'ERA5_AtmosComposites_nonjja_bwfilter_3d_10percentile_dtUnique01_noshift_2014_2023.nc'

data_source = 'ERA5'

imagedir = '/home/scavallo/scripts/python_scripts/images/'
fimg_save = 'era5_comp_VRILE_bwfilter_5percentile_icemask_nonjja_2014_2023.png'

# Used if plot_option == 20 or 22 below
file_numpy = '/data1/scavallo/data/seaice/VRILEs/1979_2023/VRILE_masks_nonjja_bwfilter_3d_10percentile_dtUnique01_noshift_2014_2023.npz'

# Below is only used if plot_option == 21 below
file_vrile_locations = '/data1/scavallo/data/seaice/VRILEs/1979_2022/VRILE_locations_nonjja_bwfilter_3d_10percentile_dtUnique01_2007010100_2022123118.dat'

# Below is only used if plot_option == 22 below
file_seaiceage = '/data1/scavallo/data/seaice/ice_age/'

cbarlabel = 'Number of VRILEs (2014-2023)' 
showFig = True
num_smooth_iterations = 0
standardize_anomalies = False
plot_option = 22 # 1 for mslp, 
                 # 2 for mslp anomalies
		 # 3 for mslp and trth
		 # 4 for mslp and 10-m wind
		 # 5 for mslp gradient and mslp
		 # 6 for trth and mslp 
		 # 7 for trth anomalies and mslp 
		 # 8 for trwind and mslp
		 # 9 for tru anomalies and mslp
		 # 10 for trv anomalies and mslp
		 # 20 for VRILE density, slp, trth 
		 # 21 for VRILE locations, slp, trth
		 # 22 for mslp, slp, sea ice age (jja)
		 # 23 for mslp gradient, slp, sea ice age (non-jja)

map_projection = 'npstere' # 'ortho' for orthographic projection, 'lcc' for Lambert Conformal projection 
proj_latlon = [55. , 270.]
zoom = 'false'
cen_lon = 270.
label_fontsize = 22
minlat = 25

###################################
# END user options
###################################
#static_path = '/data2/scavallo/era_interim/long_term_means/era_static_fields.nc'
static_path = '/data2/scavallo/ERA5/lsm_monthly/ERA5_lsm_fields_jja_1979_2022.nc'

#fpath = fdir + fname
f = netCDF4.Dataset(static_path,'r')
#lsm = f.variables['lsm'][:,:].squeeze()
lsm_full = f.variables['lsm'][:,0,0::2,0::2].squeeze()
lsm = np.nanmean(lsm_full,0).squeeze()
f.close()
    

fpath = fdir + fname
f = netCDF4.Dataset(fpath,'r')
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
lons_dt = f.variables['longitude_dt'][:]
lats_dt = f.variables['latitude_dt'][:]
mslp = f.variables['mslp_comp'][:,:]

mslp_anom = f.variables['mslp_anom_comp'][:,:]
trth = f.variables['trth_comp'][:,:]
trth_anom = f.variables['trth_anom_comp'][:,:]
trth_anom_std = f.variables['trth_anom_std_comp'][:,:]

if data_source == 'ERA5':
    u10 = np.zeros_like(trth).astype('f')
    v10 = np.zeros_like(trth).astype('f')
    tru = np.zeros_like(trth).astype('f')
    trv = np.zeros_like(trth).astype('f')
    tru_anom = np.zeros_like(trth).astype('f')
    trv_anom = np.zeros_like(trth).astype('f')
    tru_anom_std = np.zeros_like(trth).astype('f')
    trv_anom_std = np.zeros_like(trth).astype('f')    
else:
    u10 = f.variables['u10_comp'][:,:]
    v10 = f.variables['v10_comp'][:,:]
    tru = f.variables['tru_comp'][:,:]
    tru_anom = f.variables['tru_anom_comp'][:,:]
    tru_anom_std = f.variables['tru_anom_std_comp'][:,:]
    trv = f.variables['trv_comp'][:,:]
    trv_anom = f.variables['trv_anom_comp'][:,:]
    trv_anom_std = f.variables['trv_anom_std_comp'][:,:]

f.close()

if plot_option == 1:
    plotvar = mslp
    contourvar = mslp_anom
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r     
    cbarlabel = 'hPa'

    vMinMax = [1000,1020]
    vcint = 1
    if standardize_anomalies == False:
        vMinMax2 = [-50,50]
        vcint2 = 1
        vMinMax3 = [-50,50]
        vcint3 = 1    
    elif standardize_anomalies == True:  
        vMinMax2 = [-4,4]
        vcint2 = 0.5
        vMinMax3 = [-2,2]
        vcint3 = 0.2        
elif plot_option == 2:
    plotvar = mslp_anom
    contourvar = trth_anom
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r    
    if standardize_anomalies == False:    
        vMinMax = [-5,5]
        vcint = 0.5
        vMinMax2 = [-50,50]
        vcint2 = 1
        vMinMax3 = [-50,50]
        vcint3 = 1   
        cbarlabel = 'hPa'     
    elif standardize_anomalies == True:    
        vMinMax = [-4,4]
        vcint = 0.5
        vMinMax2 = [-4,4]
        vcint2 = 0.5
        vMinMax3 = [-2,2]
        vcint3 = 0.2 
        cbarlabel = 'Standard deviations'       
elif plot_option == 3:
    plotvar = mslp
    contourvar = trth
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r  
    cbarlabel = 'hPa'    
    
    vMinMax = [1000,1020]
    vcint = 1
    vMinMax2 = [290,350]
    vcint2 = 2
    vMinMax3 = [-50,50]
    vcint3 = 1    

elif plot_option == 4:
    windmag = np.sqrt(u10**2.0 + v10**2.0)
    
    plotvar = windmag
    contourvar = mslp
    contourvar2 = trth_anom
    
    plotvar = ma.masked_where(lsm == 1, plotvar)
    
    cmap_opt = plt.cm.hot_r    
    cbarlabel = 'm s-1'    
    
    vMinMax = [0,5]
    vcint = 0.1
    vMinMax2 = [990,1028]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1    
elif plot_option == 5:    
    dpdy,dpdx = wm.gradient_sphere(mslp, lats, lons)
    
    plotvar = (np.sqrt(dpdx**2.0 + dpdy**2.0))*1000.*1000 # per 1000 km
    
    plotvar = ma.masked_where(lsm == 1, plotvar)
    
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.hot_r    
    cbarlabel = 'hPa / 1000 km'    
    
    vMinMax = [0,4]
    vcint = 0.1
    vMinMax2 = [990,1028]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1  
    
    #infile = np.load(file_numpy)   
    #icearr = infile['iceVals']
    #lon_numparr = infile['icelons']
    #lat_numparr = infile['icelats']

elif plot_option == 6:
    plotvar = trth
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.jet
    cbarlabel = 'hPa'    
    
    vMinMax = [300,375]
    vcint = 1
    vMinMax2 = [960,1040]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1  
elif plot_option == 7:
    if standardize_anomalies == False:
        plotvar = trth_anom
        vMinMax = [-5,5]
        vcint = 0.5	
        cbarlabel = 'Kelvin'
    else:
        plotvar = trth_anom_std
        vMinMax = [-2,2]
        vcint = 0.1		
        cbarlabel = 'Standard deviations'   
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r      
    
    vMinMax2 = [960,1040]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1          
elif plot_option == 8:
    plotvar = np.sqrt(tru**2.0 + trv**2.0)
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.hot_r
    cbarlabel = 'm s-1'    
    
    vMinMax = [5,30]
    vcint = 1
    vMinMax2 = [960,1040]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1    
elif plot_option == 9:
    if standardize_anomalies == False:
        plotvar = tru_anom
        vMinMax = [-10,10]
        vcint = 0.5	
        cbarlabel = 'm s-1'
    else:
        plotvar = tru_anom_std
        vMinMax = [-2,2]
        vcint = 0.1		
        cbarlabel = 'Standard deviations'   
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r      
    
    vMinMax2 = [960,1040]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1    
elif plot_option == 10:
    if standardize_anomalies == False:
        plotvar = trv_anom
        vMinMax = [-10,10]
        vcint = 0.5	
        cbarlabel = 'm s-1'   	
    else:
        plotvar = trv_anom_std
        vMinMax = [-2,2]
        vcint = 0.1		
        cbarlabel = 'Standard deviations'   
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.RdBu_r
         
    vMinMax = [-10,10]
    vcint = 0.5
    vMinMax2 = [960,1040]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1    
elif ( plot_option == 20 ):
    infile = np.load(file_numpy)   
    icearr = infile['iceVals']
    lon_numparr = infile['icelons']
    lat_numparr = infile['icelats']


    plotvar = icearr
    contourvar = mslp
    contourvar2 = trth
    
    cmap_opt = plt.cm.hot_r
    cbarlabel = 'Number of VRILEs'    
    
    vMinMax = [1,16]
    vcint = 1
    vMinMax2 = [960,1040]
    vcint2 = 1
       
    [linds] = np.where(lats>=30)
    minvalnow = np.ceil(np.min(contourvar2[linds,:]))
    print(minvalnow)
    
    vMinMax3 = [minvalnow,minvalnow+6]
    vcint3 = 2  
elif ( plot_option == 21 ):
    infile = np.loadtxt(file_vrile_locations,skiprows = 0)   
    datelist_vrile = infile[:,1]
    vrile_lons = infile[:,3]
    vrile_lats = infile[:,2]
    
    periods2plot = [1989,1999,2009]
    vrile_lats_period1 = []
    vrile_lats_period2 = []
    vrile_lats_period3 = []
    vrile_lons_period1 = []
    vrile_lons_period2 = []
    vrile_lons_period3 = []    
    
    ncases_vrile = np.size(datelist_vrile)
    tt = 0
    while tt < ncases_vrile:
        datestr_vrile = str(datelist_vrile[tt])
        yyyynow = int(datestr_vrile[0:4])
        if ( (yyyynow >= periods2plot[0]) and (yyyynow <periods2plot[1]) ):
            vrile_lats_period1.append(vrile_lats[tt])
            vrile_lons_period1.append(vrile_lons[tt])
        elif ( (yyyynow >= periods2plot[1]) and (yyyynow <periods2plot[2]) ): 
            vrile_lats_period2.append(vrile_lats[tt])
            vrile_lons_period2.append(vrile_lons[tt])	
        elif (yyyynow >= periods2plot[2]):
            vrile_lats_period3.append(vrile_lats[tt])
            vrile_lons_period3.append(vrile_lons[tt])    
    
        tt+=1
    dpdy,dpdx = wm.gradient_sphere(mslp, lats, lons)
    
    plotvar = (np.sqrt(dpdx**2.0 + dpdy**2.0))*1000.*1000 # per 1000 km
    
    plotvar = ma.masked_where(lsm == 1, plotvar)
    
    contourvar = mslp
    contourvar2 = trth_anom
    
    cmap_opt = plt.cm.hot_r    
    cbarlabel = 'hPa / 1000 km'    
    
    vMinMax = [0,4]
    vcint = 0.1
    vMinMax2 = [990,1028]
    vcint2 = 1
    vMinMax3 = [-50,50]
    vcint3 = 1  
elif ( (plot_option == 22) or (plot_option == 23) ):
    iceage_arr = []
    iceage_arr_min = []
    iceage_arr_max = []
    if plot_option == 22:
        for mm in range(0,3):        
            if mm == 0:	    
                fpath_iceage = file_seaiceage + 'seaiceage_june_1984_2021.nc'
            elif mm == 1:	    
                fpath_iceage = file_seaiceage + 'seaiceage_july_1984_2021.nc'	    
            elif mm == 2:	    
                fpath_iceage = file_seaiceage + 'seaiceage_august_1984_2021.nc'	    
	    
            f = netCDF4.Dataset(fpath_iceage,'r')
            lons_iceage = f.variables['lon'][:]
            lats_iceage = f.variables['lat'][:]
            years_iceage = f.variables['years'][:]
            #[yind] = np.where(years_iceage == year2plot)
            print(years_iceage)
	    
            iceage = np.nanmedian(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
            iceage_min = np.nanmin(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
            iceage_max = np.nanmax(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
	    
            #iceage = np.nanmedian(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()
            #iceage_max = np.nanmax(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()
            #iceage_min = np.nanmin(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()	    	    
	    
            f.close
            
            iceage = um.filter_numeric_nans(iceage,10,float('NaN'),'high')
            iceage_min = um.filter_numeric_nans(iceage_min,10,float('NaN'),'high')
            iceage_max = um.filter_numeric_nans(iceage_max,10,float('NaN'),'high')

            iceage_arr.append(iceage)
            iceage_arr_min.append(iceage_min)
            iceage_arr_max.append(iceage_max)
    else:
        for mm in range(0,9):        
            print(mm)
            if mm == 0:	    
                fpath_iceage = file_seaiceage + 'seaiceage_january_1984_2021.nc'
            elif mm == 1:	    
                fpath_iceage = file_seaiceage + 'seaiceage_february_1984_2021.nc'	    
            elif mm == 2:	    
                fpath_iceage = file_seaiceage + 'seaiceage_march_1984_2021.nc'	    	
            if mm == 3:	    
                fpath_iceage = file_seaiceage + 'seaiceage_april_1984_2021.nc'
            elif mm == 4:	    
                fpath_iceage = file_seaiceage + 'seaiceage_may_1984_2021.nc'	    
            elif mm == 5:	    
                fpath_iceage = file_seaiceage + 'seaiceage_september_1984_2021.nc'	    	
            if mm == 6:	    
                fpath_iceage = file_seaiceage + 'seaiceage_october_1984_2021.nc'
            elif mm == 7:	    
                fpath_iceage = file_seaiceage + 'seaiceage_november_1984_2021.nc'	    
            elif mm == 8:	    
                fpath_iceage = file_seaiceage + 'seaiceage_december_1984_2021.nc'	    	
	
            f = netCDF4.Dataset(fpath_iceage,'r')
            lons_iceage = f.variables['lon'][:]
            lats_iceage = f.variables['lat'][:]
            years_iceage = f.variables['years'][:]
            #[yind] = np.where(years_iceage == year2plot)
            print(years_iceage)
	    
            #iceage = np.nanmedian(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
            #iceage_max = np.nanmax(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
            #iceage_min = np.nanmin(f.variables['age_of_sea_ice'][-15:,:,:],0).squeeze()
	    
            iceage = np.nanmedian(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()
            iceage_max = np.nanmax(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()
            iceage_min = np.nanmin(f.variables['age_of_sea_ice'][0:25,:,:],0).squeeze()	    
	    
            f.close
    
            iceage = um.filter_numeric_nans(iceage,10,float('NaN'),'high')
            iceage_min = um.filter_numeric_nans(iceage_min,10,float('NaN'),'high')
            iceage_max = um.filter_numeric_nans(iceage_max,10,float('NaN'),'high')

            iceage_arr.append(iceage)
            iceage_arr_min.append(iceage_min)
            iceage_arr_max.append(iceage_max)	    
    
    iceage_comp = np.nanmedian(iceage_arr,0).squeeze()
    iceage_comp_min = np.nanmin(iceage_arr_min,0).squeeze()
    iceage_comp_max =  np.nanmax(iceage_arr_max,0).squeeze()
    
    #iceage_comp_min = np.nanmedian(iceage_arr_min,0).squeeze()
    #iceage_comp_max = np.nanmedian(iceage_arr_max,0).squeeze()
    iceage_comp_absmax = np.nanmax(iceage_arr_max,0).squeeze()
    iceage_comp_absmin = np.nanmin(iceage_arr_min,0).squeeze()
    
    contourvar2 = trth_anom

    vMinMax3 = [-50,50]
    vcint3 = 1  

    infile = np.load(file_numpy)   
    icearr = infile['iceVals']
    lon_numparr = infile['icelons']
    lat_numparr = infile['icelats']


    plotvar = icearr
    contourvar = mslp
    contourvar2 = trth
    
    cmap_opt = plt.cm.hot_r
    #cbarlabel = 'Number of VRILEs (2010-2021)'    
    
    vMinMax = [0.25,4]
    vcint = 0.5
    vMinMax2 = [960,1040]
    vcint2 = 1
    
mstats(plotvar)

if num_smooth_iterations > 0: 
    smooth_amp = 0.75
    for ii in range(0,num_smooth_iterations):
        contourvar  = ndimage.gaussian_filter(contourvar,smooth_amp)
        contourvar2  = ndimage.gaussian_filter(contourvar2,smooth_amp)
    
    
if ( (plot_option == 20) or (plot_option == 22) or (plot_option == 23) ):
    lonin = lons      
    contourvar, lons = um.addcyclic(contourvar, lonin)
    lonin = lons_dt
    contourvar2, lons_dt = um.addcyclic(contourvar2, lonin)  
else:
    lonin = lons
    plotvar, lons = um.addcyclic(plotvar, lonin)
    contourvar, dummy = um.addcyclic(contourvar, lonin)
    lonin = lons_dt
    contourvar2, lons_dt = um.addcyclic(contourvar2, lonin)

cflevs =  np.arange(vMinMax[0],vMinMax[1]+(vcint/2), vcint)


cflevs_contour1 =  np.arange(vMinMax[0],vMinMax[1]+(vcint/2), vcint)
cflevs_contour2 =  np.arange(vMinMax2[0],vMinMax2[1]+(vcint2/2), vcint2)
cflevs_contour3 =  np.arange(vMinMax3[0],vMinMax3[1]+(vcint3/2), vcint3)

golden = (np.sqrt(5)+1.)/2.
####################################################
# Figure 1
####################################################
fig = plt.figure(figsize=(8., 16./golden), dpi=128)  
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
if map_projection == 'ortho':
   if zoom == 'false':   
      m = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
               resolution = 'l', area_thresh = 1000.,ax=ax1)
   else:
      m1 = Basemap(projection='ortho', lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],
               resolution = 'l', area_thresh = 1000.,ax=ax1)      		           

      width = m1.urcrnrx - m1.llcrnrx
      height = m1.urcrnry - m1.llcrnry

      #coef = 0.5
      coef = 0.7
      width = width*coef
      height = height*coef
      m = Basemap(projection='ortho',lat_0 = proj_latlon[0], lon_0 = proj_latlon[1],resolution='l',\
          llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)		  

elif map_projection == 'lcc':
    m = Basemap(llcrnrlon=-120.0,llcrnrlat=20.,urcrnrlon=-60.0,urcrnrlat=50.0,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax1)	
elif map_projection == 'npstere':
    #m = Basemap(projection='npstere',boundinglat=minlat,lon_0=cen_lon,resolution='l')         
    m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
F = plt.gcf()  # Gets the current figure

#m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
#m.drawcountries(linewidth=1, color='#444444', zorder=5)
#m.drawstates(linewidth=0.66, color='#444444', zorder=4)
m.drawstates(color='#444444', linewidth=1.25)
m.drawcoastlines(color='#444444')
m.drawcountries(color='#444444', linewidth=1.25)
m.drawmapboundary
if ( (plot_option == 20) or (plot_option == 22) or (plot_option == 23)):
    m.fillcontinents(color='Wheat',zorder=1) 
# draw lat/lon grid lines every 30 degrees.
#m.drawmeridians(np.arange(0, 360, 30))
m.drawparallels(np.arange(-90, 90, 10))
m.drawmeridians(np.arange(0, 360, 30),labels=[True,True,True,True])

X, Y = np.meshgrid(lons, lats)   
x, y = m(X, Y)

X_dt, Y_dt = np.meshgrid(lons_dt, lats_dt)   
x_dt, y_dt = m(X_dt, Y_dt)

mstats(plotvar)


col = '0.35'   
if ( (plot_option == 20) or (plot_option == 22) or (plot_option == 23)):
    #Xice, Yice = np.meshgrid(lon_numparr, lat_numparr)   
    xice, yice = m(lon_numparr,lat_numparr)
    mstats(xice)
    CS1 = m.contourf(xice,yice, plotvar, cmap=cmap_opt,levels=cflevs_contour1, extend='both',zorder=1)
    cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
else:
    CS1 = m.contourf(x,y, plotvar, cmap=cmap_opt,levels=cflevs_contour1, extend='both',zorder=1)
    cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',extend='both',pad=0.05)
if plot_option == 6:
    CS2 = m.contour(x, y, plotvar, cflevs_contour1, colors='k', linewidths=0.25)      
    col = 'k'

mstats(contourvar)
CS2 = m.contour(x, y, contourvar, cflevs_contour2, colors=col, linewidths=2.0)      
plt.clabel(CS2, cflevs_contour2, fmt = '%4d', inline=True, fontsize=10)    
if plot_option == 20:
    CS3 = m.contour(x_dt, y_dt, contourvar2, cflevs_contour3, colors='b', linewidths=1.5)      
    plt.clabel(CS3, cflevs_contour3, fmt = '%4.1f', inline=True, fontsize=10)    
if plot_option == 21:
    #xpt,ypt = m(vrile_lons,vrile_lats)
    #CSx = m.scatter(xpt,ypt,c='m',s=80,zorder=20,marker='*')   
    xpt1,ypt1 = m(vrile_lons_period1,vrile_lats_period1)
    xpt2,ypt2 = m(vrile_lons_period2,vrile_lats_period2)
    xpt3,ypt3 = m(vrile_lons_period3,vrile_lats_period3)
    CSx1 = m.scatter(xpt1,ypt1,c='b',s=80,zorder=20,marker='*')
    CSx2 = m.scatter(xpt2,ypt2,c='g',s=180,zorder=20,marker='*')  
    CSx3 = m.scatter(xpt3,ypt3,c='r',s=280,zorder=20,marker='*')     
if ( (plot_option == 22) or (plot_option == 23) ):
    x4, y4 = m(lons_iceage,lats_iceage)
    num_smooth_iterations = 1
    if num_smooth_iterations > 0: 
        smooth_amp = 0.75
        for ii in range(0,num_smooth_iterations):
            iceage_comp  = ndimage.gaussian_filter(iceage_comp,smooth_amp)
            iceage_comp_min  = ndimage.gaussian_filter(iceage_comp_min,smooth_amp)
            iceage_comp_max  = ndimage.gaussian_filter(iceage_comp_max,smooth_amp)
            iceage_comp_absmax  = ndimage.gaussian_filter(iceage_comp_absmax,smooth_amp)
            iceage_comp_absmin  = ndimage.gaussian_filter(iceage_comp_absmin,smooth_amp)
	    
    CS5 = m.contour(x4, y4, iceage_comp, [1.0,1.001], colors='c', linewidths=3.0)     
    #CS5a = m.contour(x4, y4, iceage_comp, [0.001,1.002], colors='c', linewidths=1.5)    
    CS5a = m.contour(x4, y4, iceage_comp_max, [1.0,1.001], colors='c', linewidths=1.5)    
    #CS5a = m.contour(x4, y4, iceage_comp_min, [1.0,1.001], colors='c', linewidths=1.5)    
    
    #CS5b = m.contour(x4, y4, iceage_comp_max, [1.0,1.001], colors='c', linewidths=3.0)
    
    #CS6a = m.contour(x4, y4, iceage_comp_absmin, [0.01,0.011], colors='m', linewidths=3.0)  
    #CS6a = m.contour(x4, y4, iceage_comp_absmax, [0.01,0.011], colors='m', linewidths=3.0)  
    
cbar.set_label(cbarlabel,size=20)

plt.savefig(imagedir + fimg_save, bbox_inches='tight')

#plt.colorbar()
if (showFig==True):
	plt.show()
else:
	print('Saving to: ', showFig)
	plt.savefig(showFig); plt.close()
