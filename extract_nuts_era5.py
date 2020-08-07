#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module reads in daily average (or min/max) of a particular ERA5 variable
and conducts a simple spatial averaging over all grid cells within polygons
corresponding to NUTS classification. Timeseries of polygon-averages are 
output to a single CSV file

:Authors:
    Gaige Hunter Kerr, <gaige.kerr@jhu.edu>
"""

def extract_nuts_era5(var, var_unified, year, vnuts, focus_countries=None):
    """Function loops through polygons corresponding to NUTS' subdivisions 
    of EU countries and finds ERA5 grid cells in each subdivision. A simple 
    arithmetic average is conducted over grid cells in subdivision for 
    a variable/year of interest and a timeseries of that variable within the
    subdivision is written to an output .csv file
    
    Parameters
    ----------
    var : str
        ERA5 variable of interest; u10, v10, t2m, d2m, etc.
    var_unified : str
        Unified variable name from https://docs.google.com/spreadsheets/d/
        1_fcdKu4c8c7zUSfF16071R-6bu2UpU_esLWKSA0YJJ8/edit?usp=sharing; should 
        match the var
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
    # Convert longitude from 0-360 deg to -180-180 deg. NUTS polygons appear
    # to be in these units
    lng = (lng+180)%360-180
    # A bit kludgey, but get rid of Southern Hemisphere and parts of North 
    # America to speed up looping through lats and lons to find intersections 
    # with shapefiles
    equator = np.where(lat==0.)[0][0]
    europe_east = np.where(lng==50.)[0][0]
    europe_west = np.where(lng==-25.)[0][0]
    lat = lat[:equator+1]
    lng_relevant = np.r_[0:europe_east,europe_west:len(lng)]
    lng = lng[lng_relevant]
    varl = varl[:,:,:equator+1,lng_relevant]    
    print('Data loaded!')
    
    # Create empty pandas DataFrame with columns corresponding to the dates
    # for each entry of varl. This DataFrame will be filled in the for loop 
    # below, and each row will corresponding to a different terroritial unit
    # in the EU 
    df = pd.DataFrame(columns=dates)
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
        # stackoverflow.com/questions/7861196/check-if-a-geopoint-with-latitude-and-longitude-is-within-a-shapefile
        # for additional information
        i_inside, j_inside = [], []
        for i, ilat in enumerate(lat):
            for j, jlng in enumerate(lng): 
                point = Point(jlng, ilat)
                if polygon.contains(point) is True:
                    # Fill lists with indices of reanalysis in grid 
                    i_inside.append(i)
                    j_inside.append(j)
        # Select variable from reanalysis at grid cells 
        varl_nuts = varl[:, 0, i_inside, j_inside]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Conduct spatial averaging over region 
            varl_nuts = np.nanmean(varl_nuts, axis=1)
        # Add row corresponding to individual territory 
        row_df = pd.DataFrame([varl_nuts], index=[ifid], columns=dates)
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
        %(vnuts,var_unified,pd.to_datetime(dates[0]).strftime('%Y%m%d'),
          pd.to_datetime(dates[-1]).strftime('%Y%m%d')), sep='\t')
    return 

from datetime import datetime
# Extract all variables averaged over all NUTS units for 2020
era5vars = ['d2m', 'pev', 'sp', 'ssrd', 'swvl1', 'swvl2', 'swvl3', 'swvl4', 
    't2mmax', 't2mmin', 't2mavg', 'tp', 'u10', 'v10', 'slhf']
era5vars_unified = ['DEW', 'PE', 'SP', 'SR', 'SM1', 'SM2', 'SM3', 'SM4', 
    'Tmax', 'Tmin', 'T', 'P', 'U', 'V', 'LH']
# Loop through years of interest
for vnuts in [1,2,3]:
    for var, var_unified in zip(era5vars, era5vars_unified):
        start = datetime.now()
        print('Extracting %s for NUTS%d subdivisions!'%(var,vnuts))
        extract_nuts_era5(var, var_unified, 2020, vnuts)        
        diff = datetime.now() - start
        print('Finished in %d seconds!'%diff.seconds)
        print('\n')
