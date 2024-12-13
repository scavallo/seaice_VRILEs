# plot_nsidc_concentration_forpaper
#
# plots sea ice concentration with SLP overlaid for a specific time
#
# Steven Cavallo
# Decemeber 2024
################################
import os, sys, datetime, time
import numpy as np
import matplotlib as mpl
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy.ma as ma
import utilities_modules as um
import weather_modules as wm
import netCDF4
import pandas as pd

#import struct, numpy, gdal
from mstats import *

#Dimensions from https://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html
proj_latlon = [65. , 270.]
proj_latlon = [65. , 0.]
datezero = '2006082100'

cbar_fontsize = 20
label_fontsize = 20

fdir_out = '/data1/scavallo/data/seaice/concentration/daily/'

# Used only if concentration_source == 3:
fpath_seaice_concentration = '/data1/scavallo/data/seaice/concentration/daily/asi-n6250-20060821-v5.4.nc'

concentration_source = 3 # 1 for nsidc, 2 for era-i, 3 for Bremen
atmos_option = 2 # 1 for ERAI_sfc_' + yearatmos + '-01-01-to-' + yearatmos + '-12-31.nc
                 # 2 for era5_surface_2016073000.nc
seaice_version = 'v03r01'

fdir_in = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02202_V3/north/daily/'
fdir_work = '/home/scavallo/scripts/python_scripts/python_polar/'

fdir_atmos = '/data2/scavallo/ERA5/surface_fields/'
fname_atmos = 'era5_sfclevel_2006080100_2006083118.nc'

imagedir = '/home/scavallo/scripts/python_scripts/images/'
plot_option = 2 # 1 for concentration, 2 for extent, 3 for change in concentration
contour_option = 0 # 0 for slp in time range
                   # 1 for slp and wind minimum values in time range, 
                   # 2 for slp and wind at particular index (set ninds_before_datezero below)
		   # 3 for 2-m temperature
		   # 4 for 2-m temperature and slp
ninds_before_datezero = 0
num_smoothing_passes = 2
hr_adv = 24
ndays_back = 0
nhrs_back = hr_adv * ndays_back

if contour_option == 0:
    figname = 'conc_slp_timerange_' + datezero + '_' + str(ndays_back) + 'days.png'
elif contour_option == 1:
    figname = 'conc_slp_uvmax_timerange_' + datezero + '_' + str(ndays_back) + 'days.png'
elif contour_option == 2:
    figname = 'conc_slp_uvmax_onetime_' + datezero + '_' + str(ndays_back) + 'days.png'
elif contour_option == 3:
    figname = 'conc_t2m_' + datezero + '_' + str(ndays_back) + 'days.png'
elif contour_option == 4:
    figname = 'conc_t2m_slp_' + datezero + '_' + str(ndays_back) + 'days.png'

if concentration_source == 1:
    if 1 == 1:
        date1 = um.advance_time(datezero,-nhrs_back) 
        iternum = 0
        os.chdir(fdir_out)
        while date1 <= datezero:
            print(date1)

            yyyymmdd = date1[0:8]
            print(yyyymmdd)

            yearnow = np.int(yyyymmdd[0:4])
            if yearnow <= 1987:
                fname_in = 'seaice_conc_daily_nh_n07_' + yyyymmdd + '_' + seaice_version + '.nc'
                fpath_in = fdir_in + yyyymmdd[0:4] + '/' + fname_in
            elif ( (yearnow >=1988) & (yearnow <= 1991)):
                fname_in = 'seaice_conc_daily_nh_f08_' + yyyymmdd + '_' + seaice_version + '.nc'
                fpath_in = fdir_in + yyyymmdd[0:4] + '/' + fname_in
            elif ( (yearnow >=1992) & (yearnow <= 1995)):
                fname_in = 'seaice_conc_daily_nh_f11_' + yyyymmdd + '_' + seaice_version + '.nc'
                fpath_in = fdir_in + yyyymmdd[0:4] + '/' + fname_in
            elif ( (yearnow >=1996) & (yearnow <= 2007)): 
                fname_in = 'seaice_conc_daily_nh_f13_' + yyyymmdd + '_' + seaice_version + '.nc'
                fpath_in = fdir_in + yyyymmdd[0:4] + '/' + fname_in
            elif ( (yearnow >=2008) ):
                fname_in = 'seaice_conc_daily_nh_f17_' + yyyymmdd + '_' + seaice_version + '.nc'
                fpath_in = fdir_in + yyyymmdd[0:4] + '/' + fname_in

            fcheck = os.path.isfile(fdir_out + fname_in)
            print(fname_in)
            if fcheck == True :    
                print("File exists, skipping this time")
            else:	    
                cmdnow = 'wget -N ' + fpath_in
                try:
                    os.system(cmdnow)
                except:
                    print("File is not available")
                    pass


            date1 = um.advance_time(date1,hr_adv)	
            iternum += 1

    os.chdir(fdir_work)
    date1 = um.advance_time(datezero,-nhrs_back) 
    iternum = 0
    conc_arr = []
    while date1 <= datezero:
        print(date1)
        yyyymmdd = date1[0:8]
        yearnow = np.int(yyyymmdd[0:4])
        if yearnow < 1987:
            fpath_now = fdir_out +  '/seaice_conc_daily_nh_n07_' + yyyymmdd + '_v03r01.nc'
        elif ( (yearnow >=1987) & (yearnow <= 1991)):
            fpath_now = fdir_out +  '/seaice_conc_daily_nh_f08_' + yyyymmdd + '_v03r01.nc'
        elif ( (yearnow >=1992) & (yearnow <= 1995)):
            fpath_now = fdir_out + '/seaice_conc_daily_nh_f11_' + yyyymmdd + '_v03r01.nc'
        elif ( (yearnow >=1996) & (yearnow <= 2007)): 
            fpath_now = fdir_out + '/seaice_conc_daily_nh_f13_' + yyyymmdd + '_v03r01.nc'
        elif ( (yearnow >=2008) ):
            fpath_now = fdir_out + '/seaice_conc_daily_nh_f17_' + yyyymmdd + '_v03r01.nc'    
	
        date1 = um.advance_time(date1,hr_adv)
        try:
            f1 = netCDF4.Dataset(fpath_now,'r')
            lat_seaice = f1.variables['latitude'][:,:]
            lon_seaice = f1.variables['longitude'][:,:]
            linds = np.where(lon_seaice < 0.)
            lon_seaice[linds] = lon_seaice[linds] + 360.
            concnow = f1.variables['seaice_conc_cdr'][0,:,:].squeeze()	
            f1.close()
            conc_arr.append(concnow)  
        except:
            print("File not found.  Skipping this file")
            pass
    
    X2 = lon_seaice
    Y2 = lat_seaice    

    conc_arr = np.array(conc_arr).astype('float')
    dconc_now = conc_arr[1:,:,:] - conc_arr[0:-1,:,:]
    dconc = np.sum(dconc_now,0)


elif concentration_source == 2:
    
    yearatmos = datezero[0:4]
    
    f1= netCDF4.Dataset(fdir_atmos + fname_seaice,'r')
    lat_seaice = f1.variables['latitude'][:]
    lon_seaice = f1.variables['longitude'][:]
    concnow = f1.variables['ci'][:,:,:]
    f1.close()
    
    X2, Y2 = np.meshgrid(lon_seaice, lat_seaice)   
    del lat_seaice
    lat_seaice = Y2
    del lon_seaice
    lon_seaice = X2    

elif concentration_source == 3:
    f1 = netCDF4.Dataset(fpath_seaice_concentration,'r')
    conc_arr = f1.variables['z'][:]/100.
    lat_seaice = f1.variables['lat'][:]
    lon_seaice = f1.variables['lon'][:]
    f1.close()

    X2 = lon_seaice
    Y2 = lat_seaice   
    conc_arr = np.array(conc_arr).astype('float')


yymmdd_zero = datezero[0:8]
yyzero = datezero[0:4]
mmzero = datezero[4:6]
ddzero = datezero[6:8]
hhzero = datezero[8:10]
new_year_day = datetime.datetime(year=int(yyzero), month=int(mmzero), day=int(ddzero))
day_of_year = int(new_year_day.strftime('%j'))

if atmos_option == 1:
    sind = ((day_of_year-ndays_back)*4)-1
    eind = sind+(4*ndays_back)
    plotind = ((day_of_year)*4)-2-sind-ninds_before_datezero

    print(sind,eind, plotind)
elif atmos_option == 2:
    #datefirst = fname_atmos[19:29]
    datefirst = fname_atmos[14:24]
    print(datefirst)
    record_numnow = 0
    while (datefirst < datezero):
        print(datefirst)
        datefirst = um.advance_time(datefirst,6) 
        record_numnow += 1
    print(datefirst,record_numnow)    
    
    


if concentration_source == 2:    
    conc_arr = concnow[sind:eind,:,:].squeeze()    
    dconc_now = conc_arr[1:,:,:] - conc_arr[0:-1,:,:]
    zero_conc = np.zeros_like(conc_arr).astype('f')
    [c,d] = np.where(conc_arr[0,:,:] == 0)
    dconc = np.sum(dconc_now,0)
    dconc[c,d] = np.float('NaN')
    
yearatmos = datezero[0:4]
if atmos_option == 1:

    fname_atmos = 'ERAI_sfc_' + yearatmos + '-01-01-to-' + yearatmos + '-12-31.nc'
    #fname_atmos = 'ERAI_sfc_' + yearatmos + '-01-01-to-' + yearatmos + '-02-28.nc'
    print(fname_atmos)
    f2 = netCDF4.Dataset(fdir_atmos + fname_atmos,'r')
    atmoslat = f2.variables['latitude'][:]
    atmoslon = f2.variables['longitude'][:]
    u10 = f2.variables['u10'][sind:eind,:,:].squeeze()
    v10 = f2.variables['v10'][sind:eind,:,:].squeeze()
    slp = (f2.variables['msl'][sind:eind,:,:].squeeze())/100.
    f2.close
    numinds = eind - sind
if atmos_option == 2:
    #sind = 0
    #eind = 1

    #fname_atmos = 'erainterim_surface_' + yearatmos + mmzero + ddzero + hhzero + '.nc'
    #fname_atmos = 'ERAI_sfc_' + yearatmos + '-01-01-to-' + yearatmos + '-02-28.nc'
    print(fdir_atmos)
    print(fname_atmos)
    f2 = netCDF4.Dataset(fdir_atmos + fname_atmos,'r')
    atmoslat = f2.variables['latitude'][:]
    atmoslon = f2.variables['longitude'][:]
    u10 = f2.variables['u10'][record_numnow,:,:].squeeze()
    v10 = f2.variables['v10'][record_numnow,:,:].squeeze()
    slp = (f2.variables['msl'][record_numnow,:,:].squeeze())/100.
    t2m = f2.variables['t2m'][record_numnow,:,:].squeeze()
    f2.close
    numinds = 1
    

wind_mag = np.sqrt(u10**2.0 + v10**2.0)


if numinds > 1:
    if contour_option == 0:
        slp_plot = slp[-1,:,:].squeeze()
    if contour_option == 1:
	#slp_plot = np.nanmean(slp,0)
        slp_plot = slp[-1,:,:].squeeze()
        wind_plot = np.nanmean(wind_mag,0)
    elif contour_option == 2:
        slp_plot = slp[plotind,:,:].squeeze()
        wind_plot = wind_mag[plotind,:,:].squeeze()
    elif contour_option == 3:
        t2m_plot = t2m[plotind,:,:].squeeze()
        wind_plot = wind_mag[plotind,:,:].squeeze()
	
        dt2mdy,dt2mdx = wm.gradient_sphere(t2m_plot, atmoslat, atmoslon)
        t2mgrad = np.sqrt(dt2mdy**2.0 + dt2mdx**2.0)*1000.
	
	#t2m_plot = t2mgrad
	
        [a,b] = np.where(sfcp < 950.)
        t2m_plot[a,b] = float('NaN')
    elif contour_option == 4:
        t2m_plot = t2m[plotind,:,:].squeeze()
        wind_plot = wind_mag[plotind,:,:].squeeze()
        slp_plot = slp[plotind,:,:].squeeze()
        dt2mdy,dt2mdx = wm.gradient_sphere(t2m_plot, atmoslat, atmoslon)
        t2mgrad = np.sqrt(dt2mdy**2.0 + dt2mdx**2.0)*1000.
		
	
        [a,b] = np.where(sfcp < 950.)
        t2m_plot[a,b] = float('NaN')	
	
else:
    if contour_option == 0:
        slp_plot = slp
    if ( (contour_option == 1) or (contour_option == 2) ):
        slp_plot = slp
        wind_plot = wind_mag
    if (contour_option == 3):       
        t2m_plot = t2m
        wind_plot = wind_mag	

        dt2mdy,dt2mdx = wm.gradient_sphere(t2m_plot, atmoslat, atmoslon)
        t2mgrad = np.sqrt(dt2mdy**2.0 + dt2mdx**2.0)*1000.
	
	#t2m_plot = t2mgrad
	
        [a,b] = np.where(sfcp < 950.)
        t2m_plot[a,b] = float('NaN')
    if (contour_option == 4):       
        slp_plot = slp
        t2m_plot = t2m
        wind_plot = wind_mag	

        dt2mdy,dt2mdx = wm.gradient_sphere(t2m_plot, atmoslat, atmoslon)
        t2mgrad = np.sqrt(dt2mdy**2.0 + dt2mdx**2.0)*1000.
	
	#t2m_plot = t2mgrad
	
        [a,b] = np.where(sfcp < 950.)
        t2m_plot[a,b] = float('NaN')


	
X1, Y1 = np.meshgrid(atmoslon, atmoslat)
  
if plot_option == 1:
    cbar_min_ice = 0.1#0.15
    cbar_max_ice = 1.0
    cint_ice = 0.05#0.01
    cmap_opt_ice = plt.cm.Blues_r
    cbar_text = 'Minimum concentration'
    
    
    if 1 == 1 :
    	plotvar = np.nanmin(conc_arr,0)
    	plotvar = um.filter_numeric_nans(plotvar,cbar_min_ice,cbar_min_ice,'low')       
    	#plotvar = um.filter_numeric_nans(plotvar,cbar_min_ice,float('NaN'),'low')
    	plotvar = um.filter_numeric_nans(plotvar,cbar_max_ice,cbar_max_ice,'high')   
    else:
    	mstats(conc_arr)
    	plotvar = conc_arr[-1,:,:].squeeze()

    [a,b] = np.where(lat_seaice >= 86.5)
    plotvar[a,b] = 1.0;    
elif plot_option == 2:
    cbar_min_ice = 0.99    
    cmap_opt_ice = plt.cm.Blues_r
   
    plotvar_tmp = conc_arr
    plotvar_tmp = um.filter_numeric_nans(plotvar_tmp,cbar_min_ice,1.0,'high')    
    
    if concentration_source == 1:
        [a,b] = np.where(lat_seaice >= 87.5)
                
        plotvar_masks = np.copy(plotvar_tmp)
        plotvar_masks[:] = False
        plotvar_masks[a,b] = True

        plotvar = np.ma.array(plotvar_tmp,mask=plotvar_masks)
    else:
        plotvar = plotvar_tmp

    
elif plot_option == 3:
    cbar_min_ice = -0.3
    cbar_max_ice = 0
    cint_ice = 0.01        
    cmap_opt_ice = plt.cm.gist_heat
    dconc = um.filter_numeric_nans(dconc,cbar_min_ice,cbar_min_ice,'low')       
    dconc = um.filter_numeric_nans(dconc,cbar_max_ice,cbar_max_ice,'high')   
    cbar_text = 'Change in concentration'
    plotvar = dconc

if plot_option != 2:
    cflevs_ice =  np.arange(cbar_min_ice, cbar_max_ice + (cint_ice/2), cint_ice)  
    cflevs_ice_ticks = cflevs_ice[::10]
    cinds = np.where(np.abs(cflevs_ice_ticks) < cint_ice)
    cflevs_ice_ticks[cinds] = 0.0    
else:
    cflevs_ice = np.arange(0.9,1.002,0.005)
    cflevs_ice_ticks = cflevs_ice


cbar_min_slp = 924
cbar_max_slp = 1016
cint_slp = 4
cflevs_slp =  np.arange(cbar_min_slp, cbar_max_slp + (cint_slp/2), cint_slp)

if contour_option == 2:
    for ii in range(0,num_smoothing_passes):    
        wind_plot = ndimage.gaussian_filter(wind_plot,0.5)
    #mstats(wind_plot)
    minwindplot = np.percentile(wind_plot,75)

if contour_option == 3:
    for ii in range(0,num_smoothing_passes):    
        t2m_plot = ndimage.gaussian_filter(t2m_plot,0.5)

cbar_min_wind = 15
cbar_max_wind = 80
cint_wind = 5
cflevs_wind =  np.arange(cbar_min_wind, cbar_max_wind + (cint_wind/2), cint_wind)
#print(cflevs_wind)


cbar_min_t2m = 265
cbar_max_t2m = 300
cint_t2m = 2
cflevs_t2m =  np.arange(cbar_min_t2m, cbar_max_t2m + (cint_t2m/2), cint_t2m)

golden = (np.sqrt(5)+1.)/2.   
fig = plt.figure(figsize=(12.,12.), dpi=128)   # New figure
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])

m = Basemap(projection='npstere',boundinglat=proj_latlon[0],lon_0=proj_latlon[1],resolution='l')
m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
m.drawcountries(linewidth=1, color='#444444', zorder=5)
m.drawstates(linewidth=0.66, color='#444444', zorder=4)
m.fillcontinents(color='Wheat',lake_color='lightblue', zorder=1) 
m.drawmapboundary(fill_color='lightblue')
m.drawparallels(np.arange(-90, 90, 10))
m.drawmeridians(np.arange(0, 360, 30),labels=[True,True,True,True])
   
x1, y1 = m(X1, Y1)
x2, y2 = m(X2, Y2)

# Mask array above certain latitude where data are NaNs
#plotvar = np.ma.masked_where(lat_seaice>=86.5,plotvar)
print(cflevs_ice)
CS1 = m.contourf(x2,y2,plotvar,cmap=cmap_opt_ice,levels=cflevs_ice, extend='both',zorder=1)       	    
cb = plt.colorbar(CS1,ax=ax1,orientation="horizontal",pad=0.05)
cb.set_label(label='Concentration',size=16,weight='bold')
cb.ax.tick_params(labelsize=16)
if plot_option != 2:
    cbar = plt.colorbar(CS1, shrink=0.95, orientation='horizontal',ticks=cflevs_ice_ticks,extend='both',pad=0.05)	
    cbar.set_label(cbar_text,size=label_fontsize)
    cbar.ax.set_xticklabels(cflevs_ice_ticks,size=cbar_fontsize)
if contour_option == 0:
    CS2 = m.contour(x1,y1,slp_plot,levels=cflevs_slp, colors = '0.5', linewidths=2.5)   
    plt.clabel(CS2, cflevs_slp, fmt = '%i', inline=True, fontsize=10)
if ( (contour_option == 1) or (contour_option == 2) ):
    CS2 = m.contour(x1,y1,slp_plot,levels=cflevs_slp, colors = '0.5', linewidths=2.5)   
    plt.clabel(CS2, cflevs_slp, fmt = '%i', inline=True, fontsize=10)
    CS3 = m.contour(x1,y1,wind_plot,levels=cflevs_wind, colors = 'g', linewidths=2.75)   
    plt.clabel(CS3, cflevs_wind, fmt = '%i', inline=True, fontsize=10)   
if contour_option == 3:
    CS2 = m.contour(x1,y1,t2m_plot,levels=cflevs_t2m, colors = 'r', linewidths=1.5)   
    plt.clabel(CS2, cflevs_t2m, fmt = '%i', inline=True, fontsize=10)
if contour_option == 4:
    CS2 = m.contour(x1,y1,t2m_plot,levels=cflevs_t2m, colors = 'r', linewidths=1.5)   
    plt.clabel(CS2, cflevs_t2m, fmt = '%i', inline=True, fontsize=10)
    CS3 = m.contour(x1,y1,slp_plot,levels=cflevs_slp, colors = '0.5', linewidths=3)   
    plt.clabel(CS3, cflevs_slp, fmt = '%i', inline=True, fontsize=10)       
print(figname)
save_name = imagedir + figname	
plt.savefig(save_name, bbox_inches='tight')


