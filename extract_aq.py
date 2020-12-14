#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script extracts time-averaged (2014-2018) surface-level PM2.5 and NO2 
concentrations for administrative units in Europe and Brazil. 

:Authors:
    Gaige Hunter Kerr, <gaigekerr@gwu.edu>
"""
# # # # Relevant directories
DIR_ROOT = '/mnt/sahara/data1/COVID/'
DIR_AQ = DIR_ROOT+'AQ/'
DIR_SHAPE = DIR_ROOT+'geography/'
DIR_GPW = DIR_ROOT+'geography/GPW/'
DIR_OUT = DIR_ROOT+'code/dataprocessing/'

def extract_aq(vnuts): 
    """Extract 2014-2018 surface-level NO2 and PM25 estimates for 
    administrative units in Europe or Brazil. Note that this function reads
    air quality and population count datasets that have been regridded from 
    their native resolution (NO2 = 0.00833 x 0.00833 degree, PM2.5 = 0.01 x
    0.01 degree, GPW = 0.00833 x 0.00833 degree) to a uniform target grid with
    a resolution of 0.05 degres. To generate these regridded datasets, use 
    the following at the command line: 
        
    ncatted -O -a units,LAT,c,c,"degrees north" -a units,LON,c,c,"degrees east" INFILE.nc
    ncatted -a coordinates,VAR,m,c,"LON LAT" INFILE.nc
    cdo -remapbil,r7200x3600 INFILE.nc OUTFILE.nc
    
    Function calculates a simple arithmetic average for each unit as well as a 
    population-weighted average. Units of PM2.5 are ug/m3 and units of NO2 are
    ppbv. 
    
    Parameters
    ----------
    vnuts : int
        NUTS division; either 1, 2, or 3 

    Returns
    -------
    None    

    References
    ----------    
    Center for International Earth Science Information Network - CIESIN - 
    Columbia University. 2018. Gridded Population of the World, Version 4 
    (GPWv4): Population Density, Revision 11. Palisades, NY: NASA Socioeconomic 
    Data and Applications Center (SEDAC). https://doi.org/10.7927/H49C6VHW. 
    Accessed 26 November 2020. 
    
    Mohegh, Arash; Anenberg, Susan (2020): Global surface NO2 concentrations 
    1990-2019. figshare. Dataset. https://doi.org/10.6084/m9.figshare.12968114.v2.
    
    van Donkelaar, A., R.V Martin, M.Brauer, N. C. Hsu, R. A. Kahn, R. C Levy, 
    A. Lyapustin, A. M. Sayer, and D. M Winker, Global Estimates of Fine 
    Particulate Matter using a Combined Geophysical-Statistical Method with 
    Information from Satellites, Models, and Monitors, Environ. Sci. Technol, 
    doi: 10.1021/acs.est.5b05833, 2016.
    """
    import numpy as np
    import warnings
    import sys
    import pandas as pd
    import netCDF4 as nc
    import shapefile
    from shapely.geometry import shape, Point
    # Years of interest
    years = [2014, 2015, 2016, 2017, 2018]
    # Open population dataset 
    gpw = nc.Dataset(DIR_GPW+'gpw_v4_population_count_adjusted_to_2015'+\
        '_unwpp_country_totals_rev11_2020_regridded0.05.nc', 'r')
    gpw = gpw.variables['Band1'][:]
    gpw_full = gpw.filled(np.nan)
    # Open NO2 dataset for specified years 
    no2 = []
    for year in years: 
        no2_thisyear = nc.Dataset(DIR_AQ+'NO2/'+
            'anenberg_no2_%d_final_regridded0.05.nc'%(year), 'r')
        no2_thisyear = no2_thisyear['Band1'][:].data
        # Replace -999s with np.nan
        no2_thisyear[no2_thisyear==-999.] = np.nan
        no2.append(no2_thisyear)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')    
        no2_full = np.nanmean(np.stack(no2), axis=0)
    # Open PM2.5 dataset for specified years
    pm25 = []
    for year in years: 
        pm25_thisyear = nc.Dataset(DIR_AQ+'PM25/'+
            'ACAG_PM25_GWR_V4GL03_%s01_%s12_regridded0.05.nc'%(year,year), 'r')
        if year == years[0]:
            lat_full = pm25_thisyear['lat'][:].data
            lng_full = pm25_thisyear['lon'][:].data
        pm25.append(pm25_thisyear['PM25'][:].data)
    # Stack and average along time dimension
    pm25 = np.stack(pm25)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pm25 = np.nanmean(pm25, axis=0)
    # Rotate/flip PM25
    pm25_full = np.flipud(np.rot90(pm25))
    # # # # Select European domain for datasets to speed up looping 
    eu_south = np.abs(lat_full-35).argmin()
    eu_north = np.abs(lat_full-70).argmin()
    eu_east = np.abs(lng_full-50).argmin()
    eu_west = np.abs(lng_full-328).argmin()
    lat = lat_full[eu_south:eu_north]
    lng = lng_full[np.r_[:eu_east,eu_west:len(lng_full)]]
    pm25 = pm25_full[eu_south:eu_north,np.r_[:eu_east,eu_west:len(lng_full)]]
    no2 = no2_full[eu_south:eu_north,np.r_[:eu_east,eu_west:len(lng_full)]]
    gpw = gpw_full[eu_south:eu_north,np.r_[:eu_east,eu_west:len(lng_full)]]
    print('Air quality data read!')
    
    # # # # Read geometry '%vnuts
    r = shapefile.Reader(DIR_SHAPE+
        'NUTS_RG_10M_2016_4326_LEVL_%d/NUTS_RG_10M_2016_4326_LEVL_%d.shp'%(
        vnuts,vnuts))
    # Get shapes, records
    shapes = r.shapes()
    records = r.records()
    print('Shapefile read!')
    
    # Create empty pandas DataFrames for pollutants (weighted and unweighted). 
    # These dataframes will have a single column, and each row will correspond 
    # to a different administrative unit
    df_aq = pd.DataFrame(columns=['PM25','NO2'])
    df_aq_wt = pd.DataFrame(columns=['PM25','NO2'])
    # Variables for bar to indiciate progress iterating over shapes
    total = len(shapes)  # total number to reach
    bar_length = 30  # should be less than 100
    # Loop through shapes; each shapes corresponds to NUTS code
    for ishape in np.arange(0, len(shapes), 1):
        # Build a shapely polygon from shape
        polygon = shape(shapes[ishape]) 
        # Read a single record call the record() method with the record's index
        record = records[ishape]
        # Extract FID, NUTS of record 
        ifid = record['FID']
        # For each polygon, loop through model grid and check if grid cells
        # are in polygon
        i_inside, j_inside = [], []
        for i, ilat in enumerate(lat):
            for j, jlng in enumerate(lng): 
                point = Point(jlng, ilat)
                if polygon.contains(point) is True:
                    # Fill lists with indices of reanalysis in grid 
                    i_inside.append(i)
                    j_inside.append(j)
        # There are two cases where the above fails (i.e., doesn't yield 
        # points). The first is where the administrative unit is too small to 
        # intersect with any grid cells. In this case, find and average nearby
        # points. The second case is when the unit is outside the clipped 
        # region (i.e., oversees territories; e.g., Reunion Island for France).
        # In this case, set equal to NaN. The next if statement handles the 
        # second case
        if (len(i_inside)==0) and ((lng.min()<=polygon.centroid.x<=lng.max()) 
            and (lat.min()<=polygon.centroid.y<=lat.max())):
            lat_centroid = polygon.centroid.xy[1][0]
            lng_centroid = polygon.centroid.xy[0][0]
            lat_close = lat.flat[np.abs(lat-lat_centroid).argmin()]
            lng_close = lng.flat[np.abs(lng-lng_centroid).argmin()]            
            lat_close = np.where(lat==lat_close)[0][0]
            lng_close = np.where(lng==lng_close)[0][0]
            # Select and average the closest point and the nearest nine points
            #   (lat+1, lng-1)      (lat+1, lng)      (lat+1, lng+1)
            #   (lat, lng-1)         (lat, lng)       (lat, lng+1)
            #   (lat-1, lng-1)      (lat-1, lng)      (lat-1, lng+1)
            i_inside.extend([lat_close+1, lat_close+1, lat_close+1, lat_close, 
                lat_close, lat_close, lat_close-1, lat_close-1, lat_close-1])
            j_inside.extend([lng_close-1, lng_close, lng_close+1, lng_close-1, 
                lng_close, lng_close+1, lng_close-1, lng_close, lng_close+1])
        # # To check the geometries and make sure points are in the polygon
        # ax = plt.subplot2grid((1,1),(0,0), projection=ccrs.PlateCarree(
        #     central_longitude=0.))
        # ax.scatter(lng[j_inside], lat[i_inside])
        # ax.add_geometries([polygon], proj, edgecolor='k')
        if (len(i_inside)!=0) and (len(j_inside)!=0):
            pm25_unit = pm25[i_inside, j_inside]
            no2_unit = no2[i_inside, j_inside]
            gpw_unit = gpw[i_inside, j_inside]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                mask_pm25 = np.where((gpw_unit >= 0.0) & (np.isreal(pm25_unit) 
                    == True))[0]
                mask_no2 = np.where((gpw_unit >= 0.0) & (np.isreal(no2_unit) 
                    == True))[0]         
                # Total population in unit 
                sumpop = np.nansum(gpw_unit[mask_pm25])
                pm25_unit_wt = pm25_unit[mask_pm25]*gpw_unit[mask_pm25]/sumpop
                
                sumpop = np.nansum(gpw_unit[mask_no2])
                no2_unit_wt = no2_unit[mask_no2]*gpw_unit[mask_no2]/sumpop
                # Conduct population-weighted spatial average over region
                if sumpop > 0.0:
                    pm25_unit_wt = np.nansum(pm25_unit_wt)
                    no2_unit_wt = np.nansum(no2_unit_wt)
                else: 
                    pm25_unit_wt = np.nan
                    no2_unit_wt = np.nan
                # Conduct simple spatial average over region (i.e., no 
                # population weighting)
                pm25_unit = np.nanmean(pm25_unit)
                no2_unit = np.nanmean(no2_unit)
        else:
            pm25_unit = np.nan
            no2_unit = np.nan
            pm25_unit_wt = np.nan
            no2_unit_wt = np.nan
        # Add row corresponding to individual unit for unweighted case
        row_df = pd.DataFrame({'PM25':pm25_unit, 'NO2':no2_unit}, index=[ifid], 
            columns=['PM25','NO2'])
        df_aq = pd.concat([row_df, df_aq])
        # For weighted
        row_df = pd.DataFrame({'PM25':pm25_unit_wt, 'NO2':no2_unit_wt}, 
            index=[ifid], columns=['PM25','NO2'])
        df_aq_wt = pd.concat([row_df, df_aq_wt])        
        # Update progress bar
        percent = 100.0*ishape/total  
        sys.stdout.write('\r')
        sys.stdout.write("Completed: [{:{}}] {:>3}%".format('='*int(
            percent/(100.0/bar_length)), bar_length, int(percent)))
        sys.stdout.flush()
    # Write output files 
    df_aq.index.name = 'NUTS'
    df_aq.to_csv(DIR_OUT+'ERA5_NUTS%d_AQ.csv'%(vnuts), sep=',')
    df_aq_wt.index.name = 'NUTS'
    df_aq_wt.to_csv(DIR_OUT+'ERA5_NUTS%d_AQ_popwtd.csv'%(vnuts), sep=',')
    return 
    
for vnuts in [1,2,3]:
    extract_aq(vnuts)