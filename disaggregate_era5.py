
"""
"""
import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num
import netcdftime
from datetime import datetime
year= 2019
# ERA5 variables (n.b., each individual .nc file should contain hourly data 
# corresponding to individual variables in list VARS)
DDIR = '/mnt/sahara/data1/COVID/ERA5/'
# Variables as they appear in file names
FSTR = ['10m_u_component_of_wind',
    '10m_v_component_of_wind',
    '2m_dewpoint_temperature', 
    '2m_temperature',
    'evaporation',
    'potential_evaporation',
    'surface_pressure',
    'surface_solar_radiation_downwards',
    'total_precipitation',
    'volumetric_soil_water_layer_1']
# Variables as they appear as netCDF variable names. Note: the order must 
# exactly match order in FSTR
VAR = ['u10', 
    'v10',
    'd2m',
    't2m',
    'e',
    'pev',
    'sp',
    'ssrd',
    'tp',
    'swvl1']
# Loop through variables/files
for fstri in FSTR:
    var = VAR[np.where(np.array(FSTR)==fstri)[0][0]]
    infile = Dataset(DDIR+'era5-hourly-%s-%d.nc'%(fstri,year), 'r')
    tname = 'time'
    nctime = infile.variables[tname][:] # get values
    t_unit = infile.variables[tname].units # get unit (i.e., hours since...)
    try:
        t_cal = infile.variables[tname].calendar
    except AttributeError: # Attribute doesn't exist
        t_cal = u'gregorian' # or standard
    datevar = []
    datevar.append(netcdftime.num2date(nctime, units=t_unit, calendar=t_cal))
    datevar = np.array(datevar)[0]
    # Create date range that spans the time dimension
    dayrange = pd.date_range(datevar[0], datevar[-1])
    datevar = [x.date() for x in datevar]
    # Loop through days and find indices corresponding to each day
    for day in dayrange:
        day = day.date() 
        print(var, day)
        whereday = np.where(np.array(datevar)==day)[0]
        # Select variable of interest and hours in day 
        vararr = infile.variables[var]
        vararr = vararr[whereday]
        # Find daily mean of hourly values for all variables but temperature;
        # in the case of temperature, find the daily maximum and minimum 
        if var != 't2m': 
            vararr = np.nanmean(vararr,axis=0).data
            # From https://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one
            with Dataset(DDIR+'era5-hourly-%s-%d.nc'%(fstri,year), 'r') as \
                src, Dataset(DDIR+'daily/'+'ERA5_%s_%s.nc'%(
                day.strftime('%Y%m%d'),var), 'w') as dst:
                for name, dimension in src.dimensions.items():
                    if (name=='time') or (name=='expver'):
                        continue
                    dst.createDimension(name, len(dimension) if not 
                        dimension.isunlimited() else None)
                dst.createDimension('time', 1)
                for name, variable in src.variables.items():
                    if (name==var) or (name=='expver') or (name=='time'):
                        continue
                    x = dst.createVariable(name, variable.datatype, 
                        variable.dimensions)
                    dst.variables[name][:] = src.variables[name][:]
                dst.createVariable(var, np.float32, ('time', 'latitude', 'longitude'))        
                dst.variables[var][:] = vararr
                # Write timestamps to netCDF file
                time_unit_out = 'seconds since %s 00:00:00'%day.strftime('%Y-%m-%d')
                dateo = dst.createVariable('time', np.float64, ('time',))
                dateo[:] = date2num(datetime.combine(day, datetime.min.time()), 
                    units=time_unit_out)
                dateo.setncattr('unit',time_unit_out)
        else:
            # For daily maximum temperature 
            vararr_max = np.nanmax(vararr,axis=0).data
            with Dataset(DDIR+'era5-hourly-%s-%d.nc'%(fstri,year), 'r') as \
                src, Dataset(DDIR+'daily/'+'ERA5_%s_%smax.nc'%(
                day.strftime('%Y%m%d'),var), 'w') as dst:
                for name, dimension in src.dimensions.items():
                    if (name=='time') or (name=='expver'):
                        continue
                    dst.createDimension(name, len(dimension) if not 
                        dimension.isunlimited() else None)
                dst.createDimension('time', 1)
                for name, variable in src.variables.items():
                    if (name==var) or (name=='expver') or (name=='time'):
                        continue
                    x = dst.createVariable(name, variable.datatype, 
                        variable.dimensions)
                    dst.variables[name][:] = src.variables[name][:]
                dst.createVariable(var+'max', np.float32, 
                    ('time','latitude','longitude'))        
                dst.variables[var+'max'][:] = vararr_max
                time_unit_out = 'seconds since %s 00:00:00'%day.strftime('%Y-%m-%d')
                dateo = dst.createVariable('time', np.float64, ('time',))
                dateo[:] = date2num(datetime.combine(day, datetime.min.time()), 
                    units=time_unit_out)
                dateo.setncattr('unit',time_unit_out)                
            # For daily minimum temperature       
            vararr_min = np.nanmin(vararr,axis=0).data                      
            with Dataset(DDIR+'era5-hourly-%s-%d.nc'%(fstri,year), 'r') as \
                src, Dataset(DDIR+'daily/'+'ERA5_%s_%smin.nc'%(
                day.strftime('%Y%m%d'),var), 'w') as dst:
                for name, dimension in src.dimensions.items():
                    if (name=='time') or (name=='expver'):
                        continue
                    dst.createDimension(name, len(dimension) if not 
                        dimension.isunlimited() else None)
                dst.createDimension('time', 1)
                for name, variable in src.variables.items():
                    if (name==var) or (name=='expver') or (name=='time'):
                        continue
                    x = dst.createVariable(name, variable.datatype, 
                        variable.dimensions)
                    dst.variables[name][:] = src.variables[name][:]
                dst.createVariable(var+'min', np.float32, 
                    ('time','latitude','longitude'))        
                dst.variables[var+'min'][:] = vararr_min
                time_unit_out = 'seconds since %s 00:00:00'%day.strftime('%Y-%m-%d')
                dateo = dst.createVariable('time', np.float64, ('time',))
                dateo[:] = date2num(datetime.combine(day, datetime.min.time()), 
                    units=time_unit_out)
                dateo.setncattr('unit',time_unit_out)                
