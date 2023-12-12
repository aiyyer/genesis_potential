#!/usr/bin/env python
# coding: utf-8
'''
MPI and Genesis Potential for daily era5 data

- NCSU Large Scale and Tropical Dynamics
- A. Aiyyer
- parts of code adapted from  Wenchang Yang's TCI package


Uses:
     - metpy for units aware calculations
     - fortran wrapper for Emanuel's pcmin3 subroutine to calculate TC PI
       PI refers to potential intensity (min sea level pressure and max V)
       see: pcmin_2013.f90 

Updates:
Sep 23, 2021: Modified monthly data code to daily
'''
# In[ ]:


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

# usethe fortran subroutine via the pcmin wrapper
from pcmin import pcmin3



# In[ ]:


# daily era5
era5_sfc_dir = '/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/'
era5_pl_dir  = '/glade/collections/rda/data/ds633.0/e5.oper.an.pl/'


# In[ ]:


latS=0
latN=30
lonW=40
lonE=85


midLevel   =  700*units.hPa
lowLevel   =  850*units.hPa
upperLevel =  200*units.hPa

print('midLevel chosen for calculation of PI and Entropy deficit = ', midLevel)
print('lowLevel chosen for calculation of abs vort and shear = ', lowLevel)
print('upperLevel chosen for calculation of shear = ', upperLevel)
print(' ')
print(' ')


# In[ ]:


# first set the years and months for the processing
# we shall allow for seclecting specific range of dates within each year

year_start  = 1991
month_start = 5
day_start   = 5

year_end  = 2000
month_end = 6
day_end   = 20

date_series = [pd.date_range(date(i,month_start,day_start),date(i,month_end,day_end), freq ='D') for i in range(year_start,year_end+1)]
# date_series is a list of lists. Lets unpack it now
dates_list = [element for sublist in date_series for element in sublist]


# In[ ]:


dates_list[0].strftime("%Y%m%d")


# In[ ]:


# physical parameters
Lv = 2.555e6*units.joules/units.kg # J/kg
g = 9.81*units.meters/(units.second*units.second) #m/s**-2
c_p = 1005.7*units.joules/(units.kelvin*units.kg) # J/kg K   
Rd = 287*units.joules/(units.kelvin*units.kg) # J/kg K
Rv = 461 *units.joules/(units.kelvin*units.kg) # J/kg K   
epsilon = Rd/Rv
T0 = 273.15*units.K


# In[ ]:


def potential_intensity(sst, slp, plevels, temperature, mixingratio):
    '''xarray-wrapper of the FORTRAN module pcmin3_kflag.

    '''    
    
    nz = temperature.shape[0]
    ny = temperature.shape[1]
    nx = temperature.shape[2]
    
    
    print (nz,ny,nx)
    pmin, vmax, iflag = pcmin3(sst, slp, plevels, temperature, mixingratio, ny, nx, nz)
                
    return (pmin,vmax,iflag)


def saturated_vapor_pressure(T):
    '''calculated saturated water vapor pressure (in Pa) given temperature T (in Kelvin)'''    
    svp = (610.94*units.Pa)*exp(17.625*(T-T0)/(T-T0+243.04*units.kelvin))
    return svp

def mixing_ratio(p, e):
    '''calculate mixing ration given air pressure (p) and water vapor pressure (e)'''
    return epsilon*e/(p-e)


def moistEntropy(T, p, RH=None, q=None):
    '''calculate moist entropy given air temperature (T), pressure (p) and relative humidity (RH, 0-1) or specific humidity (q).
    The equation is: s = c_p*log(T) - Rd*log(p_d) + Lv*r_v/T - Rv*r_v*log(RH)'''
    if RH is None:
        assert q is not None, 'at least one of the two variables must be specified: relative humidity/specific humidity'
        r_v = mixing_ratio_by_q(q)
        e = vapor_pressure_by_mixing_ratio(p, r_v)
        RH = e/saturated_vapor_pressure(T)
    else: # RH is from input directly
        e = saturated_vapor_pressure(T) * RH
        r_v = mixing_ratio(p, e) 
        
    s = c_p*log(T/(1*units.kelvin)) - Rd*log((p-e)/(1*units.Pa) ) + Lv*r_v/T - Rv*r_v*log(RH)    
    s.attrs['long_name'] = 'moist entropy'
    #s.attrs['units'] = 'J/K/kg'
    return s


def entropyDeficit(sst, slp, Tb, RHb, p_m, Tm, RHm):
    '''calculate entropy deficity defined in Tang and Emanuel, 2012.
    sst: sea surface temperature (in Kelvin);
    slp: sea level pressure (in Pa);
    Tb: boundary layer air temperature;
    RHb: boundary layer relative humidity (0-1);
    p_m: middle troposphere pressure level (usually 6e4 Pa);
    Tm: middle troposphere air temperature;
    RHm: middle troposphere relative humidity (0-1).'''
 
    s_sst_star = moistEntropy(T=sst, p=slp, RH=1)
    s_b        = moistEntropy(T=Tb, p=slp, RH=RHb)
    s_m_star   = moistEntropy(T=Tm, p=p_m, RH=1)
    s_m        = moistEntropy(T=Tm, p=p_m, RH=RHm)
    
    # pipe/lambda messes with units interesting....
    chi = (s_m_star - s_m)/(s_sst_star - s_b)  #.pipe(lambda x: x.where(x>0)) # exclude values <= 0 
    
    #print(s_sst_star)
    chi.attrs['long_name'] = 'entropy deficit'
    return chi


# In[ ]:


for a_date in dates_list:
    print( a_date.strftime('%Y%m%d') )
    b_date = a_date + pd.DateOffset(hours=23)
       
            
        
    # All surface fields first 
    # these are stored in monthly files
    #--------------------------------------------------------------------------------------------------- 
    # sst (2D field)  Has data for entire month in one file
    varId  = '034'
    varNam = 'sstk'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.sfc.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m')  + '0100_' +  a_date.strftime('%Y%m') \
    + str(a_date.days_in_month) + '23.nc'
    
    infile = era5_sfc_dir + fname
    dsst  = xr.open_dataset(infile)   
    sst = dsst.SSTK.sel(time=slice(a_date,b_date), latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    print("Read SST: ", sst.min().values,  sst.max().values )
    dsst.close()
     
    #
    # Sea Level Pressure (2D field)
    varId  = '151'
    varNam = 'msl'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.sfc.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m')  + '0100_' +  a_date.strftime('%Y%m') \
    + str(a_date.days_in_month) + '23.nc'
        
    
    infile = era5_sfc_dir + fname    
    print(infile)
    dslp  = xr.open_dataset(infile)   
    slp = dslp.MSL.sel(time=slice(a_date,b_date), latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    print("Read MSLP: ", slp.min().values,  slp.max().values )
    dslp.close()
    #
 
    
    # 2m Temperature
    varId  = '167'
    varNam = '2t'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.sfc.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m')  + '0100_' +  a_date.strftime('%Y%m') \
    + str(a_date.days_in_month) + '23.nc'
     
    
    infile = era5_sfc_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    T2m = ds.VAR_2T.sel(time=slice(a_date,b_date), latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    print("Read 2T: ", T2m.min().values,  T2m.max().values )
    ds.close()
    #
    

    # 2m Dew Point
    varId  = '168'
    varNam = '2d'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.sfc.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m')  + '0100_' +  a_date.strftime('%Y%m') \
    + str(a_date.days_in_month) + '23.nc'
     
    
    infile = era5_sfc_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    DP2m = ds.VAR_2D.sel(time=slice(a_date,b_date), latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    print("Read 2DP: ", DP2m.min().values,  DP2m.max().values )
    ds.close()
    #
    
    RH2m = relative_humidity_from_dewpoint(T2m, DP2m)   
    print("2m Rel Hum = ", RH2m.min().values , RH2m.max().values)
    # RH2m is dimensionless values 0-1

    
     #--------------------------------------------------------------------------------------------------- 
   
    # Now data from pressure levels
    # ERA interim has Specific Hnmidity. We need mixing ratio
    # so, first read specific Humidity Specific Humidity q (3D field) 
    # 
    varId  = '133'
    varNam = 'q'

    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
    
    infile = era5_pl_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    q = ds.Q.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    ds.close()
   
    #era5 vertical leves are arranged from top to bottom
    q = q.reindex(level=q.level[::-1])
    
    print("Read Q: ", q.max().values,  q.min().values )
    # now use metpy function to calculate mixing ratio from specific humidity
    mixingRatio = mixing_ratio_from_specific_humidity(q) #kg/kg
    print ('mixing ratio min/max=', mixingRatio.min().values, mixingRatio.max().values)
    # output mixingRatio should be in kg/kg, i.e. dimensionless
    ds.close()
    
    
    # temperature (3D field)
    varId  = '130'
    varNam = 't'
 
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
     
    
    infile = era5_pl_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    temperature = ds.T.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE)).mean("time", keep_attrs=True).metpy.quantify()
    temperature = temperature.reindex(level=temperature.level[::-1])
    plevels = temperature.level.metpy.convert_units('hPa')
    ds.close()
    print ('T min/max=', temperature.min().values, temperature.max().values)
    
    
    # mid level rel humidity
    varId  = '157'
    varNam = 'r'
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025sc.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
       
    
    infile = era5_pl_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    relHum = ds.R.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=midLevel).mean("time", keep_attrs=True).metpy.quantify()
    #print ('relHum min/max at ', midLevel, ' hPa=', relHum.min().values, relHum.max().values)
    # convert to dimensionless 
    relHum = relHum/([100.]*units.percent)
    print ('relHum min/max at midLevel hPa=', relHum.min().values, relHum.max().values)


    
    # winds for shear and vorticity 
    varId  = '131'
    varNam = 'u'
    
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025uv.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
        
    
    infile = era5_pl_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    u = ds.U.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=lowLevel).mean("time", keep_attrs=True).metpy.quantify()
    print ('u min/max at low level hPa=', u.min().values, u.max().values)

    u2 = ds.U.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=upperLevel).mean("time", keep_attrs=True).metpy.quantify()
    print ('u2 min/max at upper level hPa=', u2.min().values, u2.max().values)
     

    varId  = '132'
    varNam = 'v'
  
    fname =   a_date.strftime('%Y%m') + '/e5.oper.an.pl.128_'+varId+'_'+varNam+'.ll025uv.' \
    +   a_date.strftime('%Y%m%d')  + '00_' +  a_date.strftime('%Y%m%d') + '23.nc'
      
    infile = era5_pl_dir + fname    
    print(infile)
    ds  = xr.open_dataset(infile)   
    v = ds.V.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=lowLevel).mean("time", keep_attrs=True).metpy.quantify()
    print ('v min/max at low level hPa=', v.min().values, v.max().values)
    v2 = ds.V.sel(latitude=slice(latN,latS), longitude=slice(lonW,lonE), level=upperLevel).mean("time", keep_attrs=True).metpy.quantify()
    print ('v2 min/max at upper hPa=', v2.min().values, v2.max().values)
               
        
    # Absolute vorticity at lowLevel (standard to take 850 hPa)
    absVor = absolute_vorticity(u,v,latitude=u.latitude)
    print("Abs Vor = ", absVor.min().values , absVor.max().values)
    
    Vshear = ( (u2 - u)**2 + (v2 - v)**2 )**0.5
    Vshear.attrs['long_name'] = 'wind shear'
    Vshear.attrs['units'] = 'm/s'    
    
    
    
    # Entropy deficit
    #p_m  = midLevel.to(units.pascal)
    #chi = entropyDeficit(sst=sst,slp=slp,Tb=T2m,RHb=RH2m,
    #        p_m=p_m,
    #        Tm=temperature.sel(level=midLevel),
    #        RHm=relHum)
            
    #print("Chi = ", chi.min().values , chi.max().values)
    
    P  = potential_intensity(sst, slp, plevels, temperature, mixingRatio)
  
    
    # Potential Intensity calculation
    vmax = xr.DataArray(name="vmax", dims=["latitude","longitude"], 
                    coords= dict(latitude=(["latitude"],sst.latitude.values),
                                 longitude=(["longitude"], sst.longitude.values)),
                    attrs=dict(long_name= 'maximum surface wind speed', units = 'm/s') )

    pmin = xr.DataArray(name="pmin", dims=["latitude","longitude"],
                    coords= dict(latitude=(["latitude"],   sst.latitude.values),
                                 longitude=(["longitude"], sst.longitude.values)),
                    attrs=dict(long_name= 'mininum central pressure', units = 'hPa') )

    GPI = xr.DataArray(name="GPI", dims=["latitude","longitude"],
                    coords= dict(latitude=(["latitude"],   sst.latitude.values),
                                 longitude=(["longitude"], sst.longitude.values)),
                    attrs=dict(long_name= 'Genesis Potential Index (Emanuel 2011)', units = '') )

    print(np.nanmin(P[0]), np.nanmax(P[1]))
    vmax = P[1]
    pmin = P[0]
    
    #----------------------------------------------------------------------------------------------------
    # GPI (Emanuel and Nolan 2004): |10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)
    GPI = (1e5 * absolute(absVor) )**(3/2) \
           * (relHum/0.5)**3 \
            * (vmax/70)**3 \
            * (1*units('m/s')+0.1*Vshear)**(-2)
    GPI.attrs['long_name'] = 'Genesis Potential Index (Emanuel & Nolan 2004)'
    #----------------------------------------------------------------------------------------------------

    # GPI (Emanuel 2010): |\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)
    #GPI = (absolute(1e5*absVor)**3) * (chi.where(chi>0)**(-4/3)) * ((vmax - 35).clip(min=0)**2) * ((25*units('m/s') + Vshear)**(-4))
    #----------------------------------------------------------------------------------------------------

    print("GPI = " , GPI.min().values, GPI.max().values, GPI.mean().values)

    dname = 'GPI'
    ofile = '/glade/scratch/aiyyer/data/era5_EN_GPI/EN' + a_date.strftime('%Y%m%d')+ '.nc'
    GPI.to_dataset(name=dname).to_netcdf(ofile, "w", 
                                         encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel':1}},
                                         )

    print('[saved]:', ofile)
    


# In[ ]:




