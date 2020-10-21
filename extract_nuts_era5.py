#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module reads in daily average (or min/max) of a particular ERA5 variable
and conducts a simple spatial averaging over all grid cells within polygons
corresponding to NUTS classification or in Brazil. Timeseries of 
polygon-averages are output to a single CSV file

:Authors:
    Gaige Hunter Kerr, <gaigekerr@gwu.edu>
"""

def extract_nuts_era5(var, year, vnuts, focus_countries=None):
    """Function loops through polygons corresponding to NUTS' subdivisions 
    of EU countries and finds ERA5 grid cells in each subdivision. A simple 
    arithmetic average is conducted over grid cells in subdivision for 
    a variable/year of interest and a timeseries of that variable within the
    subdivision is written to an output .csv file
    
    Parameters
    ----------
    var : str
        ERA5 variable of interest; u10, v10, t2m, d2m, etc.
    year : int
        Year of interest
    vnuts : int
        NUTS division; either 1, 2, or 3 
    focus_countries : list, optional
        Two-letter codes referencing the countries of interest
        
    Returns
    -------
    None              
    """
    import numpy as np
    import warnings
    import sys
    import os, fnmatch
    import pandas as pd
    from netCDF4 import num2date, Dataset
    import shapefile
    from shapely.geometry import shape, Point
    
    # Relevant directories
    DIR_ROOT = '/mnt/sahara/data1/COVID/' 
    DIR_ERA = DIR_ROOT+'ERA5/daily/%d/'%year
    DIR_SHAPE = DIR_ROOT+'geography/NUTS_shapefiles/NUTS_RG_10M_2016_4326_LEVL_%s/'%vnuts
    DIR_OUT = DIR_ROOT+'code/dataprocessing/ERA5tables/'

    # Search ERA5 directory for all daily files of variable
    infiles = fnmatch.filter(os.listdir(DIR_ERA), 'ERA5_*_%s.nc'%var)
    infiles = [DIR_ERA+x for x in infiles]
    infiles.sort()
    
    # Open files and extract variable of interest and dimensional information
    dates = []
    varl = []
    for filen in infiles: 
        infile = Dataset(filen, 'r')
        # On first iteration extract lat/lon
        if filen == infiles[0]:
            lat = infile.variables['latitude'][:].data # degrees north 
            lng = infile.variables['longitude'][:].data # degrees east
        date = infile.variables['time']
        date = num2date(date[:], date.unit)
        dates.append(str(date[0]))
        varl.append(infile.variables[var][:])
    varl = np.array(varl)
    # Apply land-ocean mask; the documentation from ERA5 states tha t
    # grid boxes with a value of 0.5 and below can only be comprised of a 
    # water surface, so we set values less than this to 0 
    lsm = DIR_ROOT+'geography/era5-land_sea_mask.nc'
    lsm = Dataset(lsm, 'r')
    lsm = lsm.variables['lsm'][0].data
    mask = np.where(lsm<=0.5, np.nan, 0)
    # Add mask to variable 
    varl = varl+mask
    
    # Convert longitude from 0-360 deg to -180-180 deg. NUTS polygons appear
    # to be in these units
    lng = (lng+180)%360-180
    # A bit kludgey, but get rid of Southern Hemisphere and parts of North 
    # America to speed up looping through lats and lons to find intersections 
    # with shapefiles
    equator = np.where(lat==0.)[0][0]
    europe_east = np.where(lng==50.)[0][0]
    europe_west = np.where(lng==-70.)[0][0]
    lat = lat[:equator+1]
    lng_relevant = np.r_[0:europe_east,europe_west:len(lng)]
    lng = lng[lng_relevant]
    varl = varl[:,:,:equator+1,lng_relevant]    
    print('Data loaded!')
    
    # Create empty pandas DataFrame with columns corresponding to the dates
    # for each entry of varl. This DataFrame will be filled in the for loop 
    # below, and each row will corresponding to a different terroritial unit
    # in the EU 
    df = pd.DataFrame(columns=[x[:-9] for x in dates])
    # Read NUTS for appropriate division; note that the "RG" files should 
    # be used...the release-notes.txt file indicates that these are the 
    # RG: regions (multipolygons), which are appropriate for 
    r = shapefile.Reader(DIR_SHAPE+'NUTS_RG_10M_2016_4326_LEVL_%d.shp'%vnuts)
    # Get shapes, records
    shapes = r.shapes()
    records = r.records()
    print('Shapefile read!')

    # If it is desired to only create tables for particular countries, 
    # a list containing two-letter codes for focus countries should be passed
    # to the fuction as the optional argument "focus_countries"
    if focus_countries: 
        focus_countries_idx = []
        # Get NUTS-X codes for records in shapefile
        for i in np.arange(0, len(shapes), 1):
            record = records[i]
            country = record['FID']
            if country.startswith(tuple(focus_countries)):
                focus_countries_idx.append(i)
        shapes = map(shapes.__getitem__, focus_countries_idx)
        shapes = list(shapes)
        records = map(records.__getitem__, focus_countries_idx)
        records = list(records)
        
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
        # are in polygon (semi-slow and a little kludgey); see 
        # stackoverflow.com/questions/7861196/check-if-a-geopoint-with-
        # latitude-and-longitude-is-within-a-shapefile
        # for additional information
        i_inside, j_inside = [], []
        for i, ilat in enumerate(lat):
            for j, jlng in enumerate(lng): 
                point = Point(jlng, ilat)
                if polygon.contains(point) is True:
                    # Fill lists with indices of reanalysis in grid 
                    i_inside.append(i)
                    j_inside.append(j)
        # If the NUTS unit is too small to not intersect with the ERA5 grid
        # pick off and average the nearest 9 point
        if len(i_inside)==0:
            lat_centroid = polygon.centroid.xy[1][0]
            lng_centroid = polygon.centroid.xy[0][0]
            lat_close = lat.flat[np.abs(lat-lat_centroid).argmin()]
            lng_close = lng.flat[np.abs(lng-lng_centroid).argmin()]            
            lat_close = np.where(lat==lat_close)[0][0]
            lng_close = np.where(lng==lng_close)[0][0]
            # Select nearest six points
            # (lat_close+1, lng_close-1) (lat_close+1, lng_close)  (lat_close+1, lng_close+1)
            # (lat_close, lng_close-1)  (lat_close, lng_close)  (lat_close, lng_close+1)
            # (lat_close-1, lng_close-1) (lat_close-1, lng_close)  (lat_close-1, lng_close+1)
            i_inside.extend([lat_close+1, lat_close+1, lat_close+1, lat_close, 
                lat_close, lat_close, lat_close-1, lat_close-1, lat_close-1])
            j_inside.extend([lng_close-1, lng_close, lng_close+1, lng_close-1, 
                lng_close, lng_close+1, lng_close-1, lng_close, lng_close+1])
        # Select variable from reanalysis at grid cells 
        varl_nuts = varl[:, 0, i_inside, j_inside]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Conduct spatial averaging over region 
            varl_nuts = np.nanmean(varl_nuts, axis=1)
        # Add row corresponding to individual territory 
        row_df = pd.DataFrame([varl_nuts], index=[ifid], 
            columns=[x[:-9] for x in dates])
        df = pd.concat([row_df, df])
        # Update progress bar
        percent = 100.0*ishape/total  
        sys.stdout.write('\r')
        sys.stdout.write("Completed: [{:{}}] {:>3}%"
                         .format('='*int(percent/(100.0/bar_length)),
                                 bar_length, int(percent)))
        sys.stdout.flush()
    # Add index name and write, filename indicates NUTS division level, 
    # variable included in .csv file, and start/end dates of data
    df.index.name = 'NUTS'
    df.to_csv(DIR_OUT+'ERA5_NUTS%d_%s_%s_%s.csv'
        %(vnuts,var,pd.to_datetime(dates[0]).strftime('%Y%m%d'),
          pd.to_datetime(dates[-1]).strftime('%Y%m%d')), sep=',')
    return 

def extract_ibge_era5(var, year, adminlev):
    """Function loops through polygons corresponding to Brazilian ADMIN1-3 
    subdivisions find ERA5 grid cells in each subdivision. A simple average is
    thereafter conducted over these grid cells. If the polygon does not contain
    any ERA5 grid cells, a simple average value is determined with the nearest 
    nine grid cells. Output .csv file contains timeseries of variable of 
    interest for all administrative units in Brazil at a given level. 
    
    Parameters
    ----------
    var : str
        ERA5 variable of interest; u10, v10, t2m, d2m, etc.
    year : int
        Year of interest
    adminlev : int
        Administrative level; either 1, 2, or 3 
        
    Returns
    -------
    None              
    """
    import numpy as np
    import warnings
    import sys
    import os, fnmatch
    import pandas as pd
    from netCDF4 import num2date, Dataset
    import shapefile
    from shapely.geometry import Point, shape
        
    # As per Hamada's unified conventions, the second position should be
    # letter codes corresponding to the state/district. IBGE conventions 
    # for numerical representation of each current federation unit are found
    # at https://en.wikipedia.org/wiki/ISO_3166-2:BR
    code_to_abbrev = {11:'RO',12:'AC',13:'AM',14:'RR',15:'PA',16:'AP',
        17:'TO',21:'MA',22:'PI',23:'CE',24:'RN',25:'PB',26:'PE',27:'AL',
        28:'SE',29:'BA',31:'MG',32:'ES',33:'RJ',35:'SP',41:'PR',42:'SC',
        43:'RS',50:'MS',51:'MT',52:'GO',53:'DF'}
    
    # Relevant directories
    DIR_ROOT = '/mnt/sahara/data1/COVID/' 
    DIR_ERA = DIR_ROOT+'ERA5/daily/%d/'%year
    DIR_OUT = DIR_ROOT+'code/dataprocessing/ERA5tables/'
    DIR_SHAPE = DIR_ROOT+'geography/BR/'

    # Search ERA5 directory for all daily files of variable
    infiles = fnmatch.filter(os.listdir(DIR_ERA), 'ERA5_*_%s.nc'%var)
    infiles = [DIR_ERA+x for x in infiles]
    infiles.sort()
    
    # Open files and extract variable of interest and dimensional information
    dates = []
    varl = []
    for filen in infiles: 
        infile = Dataset(filen, 'r')
        # On first iteration extract lat/lon
        if filen == infiles[0]:
            lat = infile.variables['latitude'][:].data # degrees north 
            lng = infile.variables['longitude'][:].data # degrees east
        date = infile.variables['time']
        date = num2date(date[:], date.unit)
        dates.append(str(date[0]))
        varl.append(infile.variables[var][:])
    varl = np.array(varl)
    # Apply land-ocean mask; the documentation from ERA5 states tha t
    # grid boxes with a value of 0.5 and below can only be comprised of a 
    # water surface, so we set values less than this to 0 
    lsm = DIR_ROOT+'geography/era5-land_sea_mask.nc'
    lsm = Dataset(lsm, 'r')
    lsm = lsm.variables['lsm'][0].data
    mask = np.where(lsm<=0.5, np.nan, 0)
    # Add mask to variable 
    varl = varl+mask
    
    # Create empty pandas DataFrame with columns corresponding to the dates
    # for each entry of varl. This DataFrame will be filled in the for loop 
    # below, and each row will corresponding to a different administrative unit
    # in Brazil
    df = pd.DataFrame(columns=[x[:-9] for x in dates])
        
    # Convert longitude from 0-360 deg to -180-180 degrees
    lng = (lng+180)%360-180
    # Select Brazilian domain
    brazil_east = np.where(lng==-30.)[0][0]
    brazil_west = np.where(lng==-76.)[0][0]
    brazil_north = np.where(lat==7.)[0][0]
    brazil_south = np.where(lat==-34.)[0][0]
    lat = lat[brazil_north:brazil_south+1]
    lng = lng[brazil_west:brazil_east+1]
    varl = varl[:,:,brazil_north:brazil_south+1,brazil_west:brazil_east+1]    
    print('Data loaded!')
        
    # I couldn't get the damn geobr package to work so I downloaded state, 
    # local, and municipal files from https://www.ibge.gov.br/en/geosciences/
    # territorial-organization/territorial-organization/18890-municipal-mesh.html?=&t=downloads        
    # From this page go to municipio_2018 -> Brasil -> BR -> BR.zip
    # Admin 1 are states/federal district
    if adminlev==1:
        # Admin 1 are states/federal district, "unidades de federacao"
        r = shapefile.Reader(DIR_SHAPE+'BRUFE250GC_SIR.shp')
        # Get shapes, records
        shapes = r.shapes()
        records = r.records()
    # Admin 2 are meso regions
    if adminlev==2:
        r = shapefile.Reader(DIR_SHAPE+'BRMEE250GC_SIR.shp')
        shapes = r.shapes()
        records = r.records()
    # Admin 3 are micro regions
    if adminlev==3:
        r = shapefile.Reader(DIR_SHAPE+'BRMIE250GC_SIR.shp')
        shapes = r.shapes()
        records = r.records()
    print('Shapefile read!')    

    # Variables for bar to indiciate progress iterating over shapes
    total = len(shapes)  # total number to reach
    bar_length = 30  # should be less than 100
    # Loop through shapes; each shapes corresponds to NUTS code
    for ishape in np.arange(0, len(shapes), 1):
        # Build a shapely polygon from shape
        polygon = shape(shapes[ishape]) 
        # Read a single record call the record() method with the record's index
        record = records[ishape]
        # Build unified geospatial ID (from 
        # https://github.com/hsbadr/COVID-19)
        if adminlev==1:
            gid = ('BR'+record['CD_GEOCUF'])
        if adminlev==2:
            gid = ('BR'+record['CD_GEOCME'][:2]+
                record['CD_GEOCME'][2:])
        if adminlev==3:        
            gid = ('BR'+record['CD_GEOCMI'][:2]+
                record['CD_GEOCMI'][2:])
            # For each polygon, loop through model grid and check if grid cells
            # are in polygon (semi-slow and a little kludgey); see 
            # stackoverflow.com/questions/7861196/check-if-a-geopoint-with-
            # latitude-and-longitude-is-within-a-shapefile
            # for additional information
        i_inside, j_inside = [], []
        for i, ilat in enumerate(lat):
            for j, jlng in enumerate(lng): 
                point = Point(jlng, ilat)
                if polygon.contains(point) is True:
                    # Fill lists with indices of reanalysis in grid 
                    i_inside.append(i)
                    j_inside.append(j)
        # If the admin unit is too small to not intersect with the ERA5 grid
        # pick off and average the nearest 9 points
        if len(i_inside)==0:
            lat_centroid = polygon.centroid.xy[1][0]
            lng_centroid = polygon.centroid.xy[0][0]
            lat_close = lat.flat[np.abs(lat-lat_centroid).argmin()]
            lng_close = lng.flat[np.abs(lng-lng_centroid).argmin()]            
            lat_close = np.where(lat==lat_close)[0][0]
            lng_close = np.where(lng==lng_close)[0][0]
            # Select nearest nine points
            # (lat_close+1, lng_close-1) (lat_close+1, lng_close)  (lat_close+1, lng_close+1)
            # (lat_close, lng_close-1)  (lat_close, lng_close)  (lat_close, lng_close+1)
            # (lat_close-1, lng_close-1) (lat_close-1, lng_close)  (lat_close-1, lng_close+1)
            i_inside.extend([lat_close+1, lat_close+1, lat_close+1, lat_close, 
                lat_close, lat_close, lat_close-1, lat_close-1, lat_close-1])
            j_inside.extend([lng_close-1, lng_close, lng_close+1, lng_close-1, 
                lng_close, lng_close+1, lng_close-1, lng_close, lng_close+1])
        # # Check points in admin unit with the following
        # import matplotlib.pyplot as plt
        # plt.contourf(lng, lat, np.nanmean(varl,axis=0)[0])
        # plt.colorbar()
        # plt.title(gid)
        # plt.scatter(lng[j_inside],lat[i_inside])
        # plt.show()se, lng_close+1, lng_close-1, lng_close, lng_close+1])

        # Select variable from reanalysis at grid cells 
        varl_admin = varl[:, 0, i_inside, j_inside]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Conduct spatial averaging over region 
            varl_admin = np.nanmean(varl_admin, axis=1)
        # Add row corresponding to individual unit 
        row_df = pd.DataFrame([varl_admin], index=[gid], 
            columns=[x[:-9] for x in dates])
        df = pd.concat([row_df, df])
        # Update progress bar
        percent = 100.0*ishape/total  
        sys.stdout.write('\r')
        sys.stdout.write("Completed: [{:{}}] {:>3}%"
                         .format('='*int(percent/(100.0/bar_length)),
                                 bar_length, int(percent)))        
    # Add index name and write, filename indicates NUTS division level, 
    # variable included in .csv file, and start/end dates of data
    df.index.name = 'IBGE'
    df.to_csv(DIR_OUT+'ERA5_IBGE%d_%s_%s_%s.csv'
        %(adminlev,var,pd.to_datetime(dates[0]).strftime('%Y%m%d'),
          pd.to_datetime(dates[-1]).strftime('%Y%m%d')), sep=',')
    return

from datetime import datetime
# Extract all variables averaged over all NUTS units for 2020
era5vars = ['d2m', 'pev', 'sp', 'ssrd', 'swvl1', 'swvl2', 'swvl3', 'swvl4', 
    't2mmax', 't2mmin', 't2mavg', 'tp', 'u10', 'v10', 'slhf']
# Loop through years of interest
# for vnuts in [1,2,3]:
 #     for var in era5vars:
#         start = datetime.now()
#         print('Extracting %s for NUTS%d subdivisions!'%(var,vnuts))
#         extract_nuts_era5(var, 2020, vnuts)        
#         diff = datetime.now() - start
#         print('Finished in %d seconds!'%diff.seconds)
#         print('\n')

for adminlev in [1,2,3]:
    for var in era5vars:
        start = datetime.now()
        print('Extracting %s for Brazil Admin%d subdivisions!'%(var,adminlev))
        extract_ibge_era5(var, 2020, adminlev)     
        diff = datetime.now() - start
        print('\n')
        print('Finished in %d seconds!'%diff.seconds)
        print('\n')
        
        
        
        
