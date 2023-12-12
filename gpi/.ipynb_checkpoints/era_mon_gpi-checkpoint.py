''''
# MPI and Genesis Potential for monthly era interim data

- NCSU Large Scale and Tropical Dynamics
- A. Aiyyer
- parts of code adapted from  Wenchang Yang's TCI package


---------------------------------------------------------------------------
Uses:
     - metpy for units aware calculations
     - fortran wrapper for Emanuel's pcmin3 subroutine to calculate TC PI
       PI refers to potential intensity (min sea level pressure and max V)
       see: pcmin_2013.f90 

---------------------------------------------------------------------------
---------------------------------------------------------------------------
Updates:
Oct 20, 2021:  First version
Sep 11, 2022:  Removed code relating to graphics

---------------------------------------------------------------------------
---------------------------------------------------------------------------
---------------------------------------------------------------------------
'''



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


#----------------------------------------------------------------------------------------------------

# physical parameters
Lv = 2.555e6*units.joules/units.kg # J/kg
g = 9.81*units.meters/(units.second*units.second) #m/s**-2
c_p = 1005.7*units.joules/(units.kelvin*units.kg) # J/kg K   
Rd = 287*units.joules/(units.kelvin*units.kg) # J/kg K
Rv = 461 *units.joules/(units.kelvin*units.kg) # J/kg K   
epsilon = Rd/Rv
T0 = 273.15*units.K
#----------------------------------------------------------------------------------------------------



# first set the years and months for the processing
year_start = 1979
year_end   = 1979 #2018
date_series = [pd.date_range(date(i,6,1),date(i,10,1), freq ='MS') for i in range(year_start,year_end+1)]

dates = [element for sublist in date_series for element in sublist]
print ("First date in list = ", dates[0])
print ("Last  date in list = ", dates[-1])


#-------------------------------------------------------------------------------------------------------
def potential_intensity(sst, slp, plevels, temperature, mixingratio):
    '''xarray-wrapper of the FORTRAN module pcmin3_kflag.

    '''    
    
    nz = temperature.shape[0]
    ny = temperature.shape[1]
    nx = temperature.shape[2]
    pmin, vmax, iflag = pcmin3(sst, slp, plevels, temperature, mixingratio, ny, nx, nz)
                
   # pmin.attrs['long_name'] = 'mininum central pressure'
   # pmin.attrs['units'] = 'hPa'
   # vmax.attrs['long_name'] = 'maximum surface wind speed'
   # vmax.attrs['units'] = 'm/s'
   # iflag = iflag.astype('int32')
   # iflag.attrs['long_name'] = '1: OK; 0: no convergence; 2: CAPE routine failed.'
   
    
    pmin = xr.DataArray(data=pmin, name="pmin", dims=["latitude","longitude"],
    coords= dict(  latitude=(["latitude"],   sst.latitude.values),
                  longitude=(["longitude"], sst.longitude.values)),
    attrs=dict(long_name= 'mininum central pressure', units = 'hPa') )
    
    
    vmax = xr.DataArray(data=vmax, name="vmax", dims=["latitude","longitude"],
    coords= dict(  latitude=(["latitude"],   sst.latitude.values),
                  longitude=(["longitude"], sst.longitude.values)),
    attrs=dict(long_name= 'maximum surface wind speed', units = 'm/s') )
    
    iflag = xr.DataArray(data= iflag, name="iflag", dims=["latitude","longitude"],
    coords= dict(  latitude=(["latitude"],   sst.latitude.values),
                  longitude=(["longitude"], sst.longitude.values)),
    attrs=dict(long_name= '1: OK; 0: no convergence; 2: CAPE routine failed') )
    
    
    vmax.attrs['long_name'] = 'maximum surface wind speed'
    vmax.attrs['units'] = 'm/s'
   

    pmin=pmin.assign_coords({"time": sst.time.values})
    vmax=vmax.assign_coords({"time": sst.time.values})
    vmax=vmax.assign_coords({"time": sst.time.values})

    
    PI = xr.Dataset(dict(pmin=pmin, vmax=vmax, iflag=iflag))
    
    return PI

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

def plot_GPI(GPI):
    
    
    minlon = -170.
    maxlon =  -20.
    minlat = -4.
    maxlat =  30.

    minC = 0
    maxC = 10 #100
    intC = 1  #10
   
    levC = int((maxC-minC)/intC) + 1

# Generate figure (set its size (width, height) in inches)
    fig = plt.figure(figsize=(20, 4))

# Generate axes, using Cartopy, drawing coastlines, and adding features
    projection = ccrs.PlateCarree()

    ax = plt.axes(projection=projection)
    ax.coastlines(linewidths=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')

    # Import an NCL colormap
    newcmp = gvcmaps.BlAqGrYeOrRe

# Contourf-plot data

    contourLevels = np.arange(minC, maxC, intC)
    vmaxMap = GPI.plot.contourf(ax=ax,
                          transform=projection,
                          vmin=minC,
                          vmax=maxC,
                          levels=contourLevels,
                          cmap=newcmp,
                          add_colorbar=False)

# Add colorbar
    cbar = plt.colorbar(vmaxMap, ticks=np.arange(minC, maxC, intC))
#cbar.ax.set_yticklabels([str(i) for i in np.arange(minC, maxC, intC)])
#cbar.ax.set_yticklabels(str(conto

#=======================================================================================================
#=======================================================================================================
#=======================================================================================================


midLevel   =  600*units.hPa
lowLevel   =  850*units.hPa
upperLevel =  200*units.hPa

print('midLevel chosen for calculation of PI and Entropy deficit = ', midLevel)
print('lowLevel chosen for calculation of abs vort and shear = ', lowLevel)
print('upperLevel chosen for calculation of shear = ', upperLevel)
print(' ')
print(' ')


# loop over dates and read the required variables and calculate PI and GP
for d in dates:   
    # sst (2D field)
    fname = d.strftime('%Y%m%d') + '00.grib'
    infile = '/run/media/sagan/data100/data/era_monthly/sst/' + fname
    dsst  = xr.open_dataset(infile, engine="cfgrib", indexpath='')   
    sst = dsst.sst.sel(latitude=slice(30,-4), longitude=slice(180,360)).metpy.quantify()
    dsst.close()

    # decide whether a point has missing or non-missing SST
    # 
    oceanFlag = sst.pipe(lambda x: x*0==0)

    
    #
    # Sea Level Pressure (2D field)
    fname = d.strftime('%Y%m%d') + '00.mslp.grib'
    infile = '/run/media/sagan/data100//data/era_monthly/slp/' + fname
    dslp  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    slp = dslp.msl.sel(latitude=slice(30,-4), longitude=slice(180,360)).metpy.quantify()
    dslp.close()
    #

    # ERA interim has Specific Hunmidity. We need mixing ratio
    # so, first read specific Humidity
    # Specific Humidity q (3D field) 
    fname = d.strftime('%Y%m%d') + '00.pl.grib'
    infile = '/run/media/sagan/data100/data//era_monthly/plev/' + fname
    ds  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    q = ds.q.sel(latitude=slice(30,-4), longitude=slice(180,360)).metpy.quantify()

    # now use metpy function to calculate mixing ratio from specific humidity
    mixingRatio = mixing_ratio_from_specific_humidity(q) #kg/kg
    #print ('mixing ratio min/max=', mixingRatio.min().values, mixingRatio.max().values)
    # output mixingRatio should be in kg/kg, i.e. dimensionless
    
    #
    # temperature (3D field)
    temperature = ds.t.sel(latitude=slice(30,-4), longitude=slice(180,360)).metpy.quantify()
    plevels = ds.isobaricInhPa.metpy.convert_units('hPa')
    ds.close()

    #------------------------------------------------------------------------------
    
    print(sst.time.values)
    
    # Potential Intensity calculation
    #
    dim_z = 'isobaricInhPa'
    dim_y = 'latitude'
    dim_x = 'longitude'
    # Calculate the Potential Intensity
    PI = potential_intensity(sst, slp, plevels, temperature, mixingRatio) # dim_x, dim_y, dim_z)
    a = PI.sel(latitude=10, longitude=340, method="nearest")
    
    print("PI=", a.time.values.astype('datetime64[D]'), a.pmin.values, a.vmax.values,a.iflag.values)
    #------------------------------------------------------------------------------


    
    # now we calculate the genesis potential: Emanuel and Nolan 2004    

    # mid level rel humidity
    fname = d.strftime('%Y%m%d') + '00.pl.grib'
    infile = '/run/media/sagan/data100/data/era_monthly/plev/' + fname
    ds  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    relHum = ds.r.sel(latitude=slice(30,-4), longitude=slice(180,360), isobaricInhPa=midLevel).metpy.quantify()
    #print ('relHum min/max at ', midLevel, ' hPa=', relHum.min().values, relHum.max().values)
    # convert to dimensionless 
    relHum = relHum/([100.]*units.percent)
    #print ('relHum min/max at midLevel hPa=', relHum.min().values, relHum.max().values)

    
    # Aabsolute vorticity at lowLevel (standard to take 850 hPa)
    fname = d.strftime('%Y%m%d') + '00.uv.grib'
    infile = '/run/media/sagan/data100/data/era_monthly/plev/' + fname
    ds  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    u = ds.u.sel(latitude=slice(30,-4), longitude=slice(180,360), isobaricInhPa=lowLevel).drop('valid_time').metpy.quantify()
    v = ds.v.sel(latitude=slice(30,-4), longitude=slice(180,360), isobaricInhPa=lowLevel).drop('valid_time').metpy.quantify()
    # drop valid_time to avoid metpy warning about muliple time coords
    absVor = absolute_vorticity(u,v,latitude=u.latitude)
    #print("Abs Vor = ", absVor.min().values , absVor.max().values)
    

    #wind shear upper-lower winds (typically 200 - 850 hPa)
    u2 = ds.u.sel(latitude=slice(30,-4), longitude=slice(180,360), isobaricInhPa=upperLevel).drop('valid_time').metpy.quantify()
    v2 = ds.v.sel(latitude=slice(30,-4), longitude=slice(180,360), isobaricInhPa=upperLevel).drop('valid_time').metpy.quantify()
 
    Vshear = ( (u2 - u)**2 + (v2 - v)**2 )**0.5
    Vshear.attrs['long_name'] = 'wind shear'
    Vshear.attrs['units'] = 'm/s'    
    
    # 2m Temperature
    fname = 'T2m.' + d.strftime('%Y%m%d') + '00.grib'
    infile = '/run/media/sagan/data100/data/era_monthly/2mT/' + fname
    ds  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    T2m = ds.t2m.sel(latitude=slice(30,-4), longitude=slice(180,360)).drop('valid_time').metpy.quantify()
    #print(T2m) # needs to be 


    # 2m Dew Point
    fname = 'DP2m.' + d.strftime('%Y%m%d') + '00.grib'
    infile = '/run/media/sagan/data100/data/era_monthly/2mDP/' + fname
    ds  = xr.open_dataset(infile, engine="cfgrib", indexpath='')
    DP2m = ds.d2m.sel(latitude=slice(30,-4), longitude=slice(180,360)).drop('valid_time').metpy.quantify()
    RH2m = relative_humidity_from_dewpoint(T2m, DP2m)   
    #print("2m Rel Hum = ", RH2m.min().values , RH2m.max().values)
    # RH2m is dimensionless values 0-1

    # Entropy deficit
    p_m  = midLevel.to(units.pascal)
    chi = entropyDeficit(sst=sst,slp=slp,Tb=T2m,RHb=RH2m,
            p_m=p_m,
            Tm=temperature.sel(isobaricInhPa=midLevel),
            RHm=relHum)
            
    #print(chi)



    #----------------------------------------------------------------------------------------------------
    # GPI (Emanuel and Nolan 2004): |10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)
    #GPI = (1e5 * absolute(absVor) )**(3/2) \
    #        * (relHum/0.5)**3 \
    #        * (PI.vmax/70)**3 \
    #        * (1*units('m/s')+0.1*Vshear)**(-2)
    #GPI.attrs['long_name'] = 'Genesis Potential Index (Emanuel & Nolan 2004)'
    #----------------------------------------------------------------------------------------------------

    # GPI (Emanuel 2010): |\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)
    GPI = (absolute(1e5*absVor)**3) * (chi.where(chi>0)**(-4/3)) * ((PI.vmax - 35).clip(min=0)**2) * ((25*units('m/s') + Vshear)**(-4))
    GPI.attrs['long_name'] = 'Genesis Potential Index (Emanuel 2011)'


    
    # now discard GPI values that should be nan based on SST
    #
    #GPI = GPI.where(oceanFlag)
    
    
    print("GPI = " , GPI.min().values, GPI.max().values, GPI.mean().values)
    print(" ")
    #print(GPI)
    #plot_GPI(GPI)
    



    # uncomment below to write out to files

    dname = 'GPI'
    ofile = '/run/media/sagan/data100/data/era_monthly/gp2010/gpi'+  d.strftime('%Y%m%d') + '.nc'
    GPI.to_dataset(name=dname) \
            .to_netcdf(ofile, "w",
            encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel':1}},
                unlimited_dims='time')
    print('[saved]:', ofile)
    
