#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
Kotal Genesis Potential for daily era5 data

- NCSU Large Scale and Tropical Dynamics
- A. Aiyyer

Uses:
     - metpy for units aware calculations

Updates:
Feb 15, 2023: Initial code version 1.0
'''


# ### Genesis Potential Index By Kotal et al.
# 
# - Kotal, S. D., Kundu, P. K., & Roy Bhowmik, S. K. (2009, August). Analysis
# of cyclogenesis parameter for developing and nondeveloping low-pressure
# systems over the Indian Sea. Natural Hazards, 50 (2), 398-402. doi:334
# https://doi.org/10.1007/s11069-009-9348-5
# 
# 
# 
# \begin{align}
#     KGPP &= \frac{\zeta_{850} *M*I}{S} &\text{  $\zeta_{850}>0, M>0, I>0$}\\
#          &= 0                       &\text{ $\zeta_{850}\leq0, M\leq0, I\leq0$}
# \end{align}
# 
# where, $\zeta_{850}$ is the relative vorticity at 850 hPa, RH is the mean relative humidity between 700 and 500 hPa, and I denotes the temperature difference between 850 and 500 hPa. S is the vertical wind shear between 200 and 850 hPa.  M = $\frac{|RH - 40|}{30}$ and I = $(T_{850} - T_{500})$. M is defined for $RH > 40$.
# 

# In[7]:


import numpy as np
import xarray as xr
import pandas as pd
from datetime import date
from numpy import absolute, exp, log

# Any import of metpy will activate the accessors
from metpy.units import units

#from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import mixing_ratio_from_specific_humidity, relative_humidity_from_dewpoint

from metpy.calc import absolute_vorticity


# In[101]:





# In[21]:


# daily era5
era5_sfc_dir = '/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/'
era5_pl_dir  = '/glade/collections/rda/data/ds633.0/e5.oper.an.pl/'


latS=0
latN=30
lonW=40
lonE=85


midLevel   =  500*units.hPa
lowLevel   =  850*units.hPa
upperLevel =  200*units.hPa

RH_Level1 = 700*units.hPa
RH_Level2 = 500*units.hPa


# In[105]:


# landsea mask

era5_invar_dir = '/glade/collections/rda/data/ds633.0/e5.oper.invariant/197901/'
fname = era5_invar_dir + 'e5.oper.invariant.128_172_lsm.ll025sc.1979010100_1979010100.nc'
ds  = xr.open_dataset(fname)   

lsm = ds.LSM[0,:,:].sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE))


# In[106]:





# In[22]:


# first set the years and months for the processing
# we shall allow for seclecting specific range of dates within each year

year_start = 2001
month_start = 5
day_start   = 1

year_end = 2010
month_end = 6
day_end   = 25

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]
# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]


# In[114]:


for a_date in dates_list:
    print( a_date.strftime('%Y%m%d') )
    b_date = a_date + pd.DateOffset(hours=23)
    
    # Layer Average rel humidity
    varId  = '157'
    varNam = 'r'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    infile = era5_pl_dir + fname    
    #print(infile)
    ds  = xr.open_dataset(infile)   
    relHum = ds.R.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=slice(RH_Level2,RH_Level1)).mean(["time","level"], keep_attrs=True).metpy.quantify()
    #print ('relHum min/max at midLevel hPa=', relHum.min().values, relHum.max().values)
    ds.close()

    
    
    
    # winds for shear and vorticity 
    varId  = '131'
    varNam = 'u'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025uv.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'    
    infile = era5_pl_dir + fname    
    #print(infile)
    ds  = xr.open_dataset(infile)   
    u = ds.U.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=lowLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('u min/max at low level hPa=', u.min().values, u.max().values)

    u2 = ds.U.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=upperLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('u2 min/max at upper level hPa=', u2.min().values, u2.max().values)
    ds.close()
  

    varId  = '132'
    varNam = 'v'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025uv.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    infile = era5_pl_dir + fname    
    #print(infile)
    ds  = xr.open_dataset(infile)   
    v = ds.V.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=lowLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('v min/max at low level hPa=', v.min().values, v.max().values)
    v2 = ds.V.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=upperLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('v2 min/max at upper hPa=', v2.min().values, v2.max().values)
    ds.close()
        
   
    Vshear = ( (u2 - u)**2 + (v2 - v)**2 )**0.5
    Vshear.attrs['long_name'] = 'wind shear'
    Vshear.attrs['units'] = 'm/s'    
    #print ('shear min/max at upper level hPa=', Vshear.min().values, Vshear.max().values)

     
    # temperature (3D field)
    varId  = '130'
    varNam = 't'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    infile = era5_pl_dir + fname    
    #print(infile)
    ds  = xr.open_dataset(infile)   
    T2 = ds.T.sel(level=midLevel, latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    T1 = ds.T.sel(level=lowLevel, latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    ds.close()
    
    
    # Low Level vorticity
    varId  = '138'
    varNam = 'vo'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    infile = era5_pl_dir + fname    
    #print(infile)
    ds  = xr.open_dataset(infile)   
    zeta = ds.VO.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=lowLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('relHum min/max at ', midLevel, ' hPa=', relHum.min().values, relHum.max().values)
    # convert to dimensionless 
    #print ('zeta min/max at lowLevel hPa=',  zeta.min().values, zeta.max().values)
    ds.close()


    
    
    # now calculate the KGPP
    
    # adjust relative humidity
    relHum = (relHum - 40*units.percent)/(30*units.percent)

    
    relHum = relHum.where(relHum>0.)
    zeta = zeta.where(zeta>0.)
    Tdiff = (T1-T2)
    Tdiff = Tdiff.where(Tdiff>0.)
    kgpp = (zeta*1.e5)*Tdiff*relHum/Vshear
    
    #kgpp = kgpp.where(lsm==0.0) # apply land mask

    dname = 'KGPP'

    ofile = '/glade/scratch/aiyyer/data/kgpp_gpi/kgpp' + a_date.strftime('%Y%m%d')+ '.nc'
    kgpp.to_dataset(name=dname).to_netcdf(ofile, "w", 
                                         encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel':1}},
                                         )


# In[112]:





# In[107]:





# In[ ]:




