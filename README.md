# README for seaice_VRILEs

Code and data for analysis in `Sea ice loss in association with Arctic cyclones' by Cavallo, Frank, and Bitz (2024)

Python notebook (*.ipynb) scripts:

compute_deltaseaice_longterm_means
    Makes daily and monthly long-term climatologies of sea ice extent with errorbars

Run make_seaice_climatology_file
    Creates another long-term climatology file necessary for plots (but not in the calculations of VRILEs themselves)

make_seaice_vriles 
    Use this script to make the VRILE files

merge_lists
    Merges the sea ice loss files using the "mean removed" and "bwfilter" methods into a combined file called "allvriles"

filter_vriles
   Removes vRILES on the last day and first two days of each month due to the erroneous NSIDC data shift (See Meier and Stewart 2019)
 *** MUST be run AFTER merge_lists ***
    
plot_seaice_stats (Figure 3 in paper)
	Plots (1) Number of events vs. year since 1979, 
          (2) Change in SIE distribution decadal comparison
          (3) Accumulated sea ice loss from events
          (4) Percent sea ice extent reduction for warm season
          (5) Sept 1 sea ice extent vs. Sept 1 without VRILE loss
	      (6) Number of events vs. month, 
	      (7) Average change in extent for all events and fractional area loss
	      (8) Average change in extent for all events vs. climatological loss

print_VRILE_stats
    prints relevant VRILE stats in the format of the paper for quick translation

find_seaice_events
    Prints and plots the top sea ice events from an existing event file

plot_seaice_spectra_forpaper (Figure 2 in paper)
    Plots of spectra in sea ice extent, change in extent, red noise
	
plot_seaice_stats_wrt_longtermmean

Regular python scripts:
    plot_nsidc_concentration_forpaper.py plots sea ice concentration with ERA5 SLP for a specific time (Figure 1 in paper)
    
    make_era5_composite_fromevents_forpaper.py compiles the ERA5 data corresponding to VRILE event files (for Figure 4 in paper)
    
    plot_seaicemasks_events_forpaper.py uses VRILE dates, obtains sea ice concentrations, and calculates the masks and locations of VRILEs (for Figure 4 in paper)
    
    plot_era5_composites_forpaper.py plots VRILE locations/masks, SLP, and ice age composites (Figure 4 in paper)
    
    make_seaiceage_netcdffile.py gets and makes monthly sea ice age files (for Figure 4 in paper)
    
    utilities_modules.py contains numerous functions, such as advancing time
    
    weather_modules.py contains numerous functions related to weather calculations
    
    spectral_utilities.py contains various functions for computing spectra
    
    mstats.py prints statistics of a python variable
    
    seaIceMasks.py contains functions needed for plot_seaicemasks_events_forpaper.py

Sea ice data:
    NSIDC sea ice extent:
       [!https://noaadata.apps.nsidc.org/NOAA/G02135/north/daily/data/](https://noaadata.apps.nsidc.org/NOAA/G02135/north/daily/data/)

    NSIDC daily gridded concentrations:
       https://noaadata.apps.nsidc.org/NOAA/G02202_V4/north/daily/

    NSIDC sea ice age
       https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0611_seaice_age_v4/

    http://nsidc.org/data/G02202/versions/2#
    https://climatedataguide.ucar.edu/search/node/sea%20ice
    https://climatedataguide.ucar.edu/climate-data/sea-ice-concentration-data-overview-comparison-table-and-graphs
   
    NOAA/NSIDC: Concentration estimates are most reliable within the consolidated ice pack during cold, winter 
    conditions (errors ~5-10%). The estimates are least reliable close to the ice edge and during melt conditions, 
    where biases may be 20-30%. At any given location (grid cell) errors may be large (>50%) 
    (https://climatedataguide.ucar.edu/climate-data/sea-ice-concentration-noaansidc-climate-data-record)
