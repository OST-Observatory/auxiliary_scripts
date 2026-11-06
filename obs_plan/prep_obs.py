#! /usr/bin/python3
# -*- coding: utf-8 -*-


############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

###
#   Site
#
obs_site = [+52.409184, +12.973185, 39]


###
#   Obs time - iso format (YYYY-MM-DD HH:MM)
#
#   If left empty, we assume the current system time
#obs_time = '2021-05-03 00:00'
obs_time = ''


###
#   Object
#
#   Object coordinates - format: ra = '00h42m30s', dec = '+41d12m00s'
ra_obj  = ''
dec_obj = ''

#   Object name -> If only the object name is given, the object will be
#                  resolved by means of Sesame, Simbad or NED (mirrors when
#                  CDS does not answer)
name_obj = '?'


############################################################################
####                            Libraries                               ####
############################################################################

import sys

import numpy as np

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.time import Time
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, name_resolve
from astropy.coordinates import get_sun
from astropy.coordinates import get_body
from astropy.visualization import time_support

from astroplan import Observer
from astroplan import FixedTarget
from astroplan.plots import plot_sky

from ost_photometry.core.services import resolve_name


############################################################################
####                        Routines & definitions                      ####
############################################################################

# Colors status massages
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    
############################################################################
####                               Main                                 ####
############################################################################

###
#   Set object
#
if ra_obj != '' and dec_obj != '':
    obj_coord = SkyCoord(
        ra    = ra_obj, 
        dec   = dec_obj, 
        frame = 'icrs',
        )
        
elif name_obj != '':
    try:
        obj_coord = resolve_name(name_obj)
    except name_resolve.NameResolveError as err:
        print(
            bcolors.WARNING
            +f"   [Error] Object name not recognized: {err}"
            +bcolors.ENDC
            )
        sys.exit()
else:
    print(
        bcolors.WARNING
        +"   [Error] Neighter object coordinates nor object name given"
        +bcolors.ENDC
        )
    sys.exit()
    
    
###
#   Location
#
location=coord.EarthLocation(
    lat=obs_site[0], 
    lon=obs_site[1], 
    height=obs_site[2],
    )

#   Define the location by means of astroplan
ost = Observer(location=location, name="OST", timezone="Europe/Berlin")


###
#   Observation time
#
#   Check time & make Time object
if obs_time == '':
    time = Time.now()
else:
    time = Time(obs_time)

#   Sunset
sunset_tonight = ost.sun_set_time(time, which='nearest')

#   Sunrise
sunrise_tonight = ost.sun_rise_time(time, which='nearest')

#   Set approximate new end time, in case start > end 
#   -> This only occurs, if obs_time is set to a time during the day aka 
#      before Sunset
if sunset_tonight > sunrise_tonight:
    sunrise_tonight = sunrise_tonight + 1*u.d
    

###
#   Astroplan objects
#

#   Object as an astroplan object
fixed_target = FixedTarget(name=name_obj, coord=obj_coord)
if not ost.target_is_up(time, fixed_target):
    print(
        bcolors.WARNING
        +"   [Error] Object not visible"
        +bcolors.ENDC
        )
    sys.exit()


###
#   Airmass plots
#
#   astroplan's plot_airmass uses Axes.plot_date, which Matplotlib removed
#   (3.10): the airmass over the night is computed and plotted here
night = sunset_tonight - 1*u.hour + np.linspace(
    0, (sunrise_tonight - sunset_tonight).to_value(u.hour) + 2, 300
    )*u.hour
night_frame = coord.AltAz(obstime=night, location=location)
altitude = obj_coord.transform_to(night_frame).alt
sun_altitude = get_sun(night).transform_to(night_frame).alt
airmass = np.where(altitude > 0*u.deg, 1/np.sin(altitude), np.nan)
local_times = night.to_datetime(timezone=ost.timezone)

#   Make figure
fig = plt.figure(figsize=(11,7))
ax = fig.add_subplot()
ax.plot(local_times, airmass, label=name_obj)

#   Twilight (light) and night (dark)
ax.fill_between(
    local_times, 0, 1,
    where=(sun_altitude < 0*u.deg) & (sun_altitude >= -18*u.deg),
    transform=ax.get_xaxis_transform(), color='grey', alpha=0.2,
    )
ax.fill_between(
    local_times, 0, 1,
    where=sun_altitude < -18*u.deg,
    transform=ax.get_xaxis_transform(), color='grey', alpha=0.5,
    )

ax.set_ylim(3, 1)
ax.set_ylabel('Airmass')
ax.set_xlabel(f'Time ({ost.timezone})')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=ost.timezone))

#   Altitude axis
altitude_axis = ax.secondary_yaxis(
    'right',
    functions=(
        lambda x: np.degrees(np.arcsin(1/np.clip(x, 1, None))),
        lambda alt: 1/np.sin(np.radians(np.clip(alt, 1, 90))),
        ),
    )
altitude_axis.set_yticks([20, 30, 40, 50, 60, 90])
altitude_axis.set_ylabel('Altitude [deg]')

#   Position of the legend
plt.legend(loc=1, bbox_to_anchor=(1, 1), shadow=True) 

#   Save figure
t            = time.replicate(format='iso', copy=True)
t.out_subfmt = 'date'
plt.savefig(
    'airmass_'+name_obj.replace(' ','_')+'_'+str(t)+'.pdf'
    )
#plt.show()
plt.close()


###
#   Up time
#
#   Rise and set times
rise_time = ost.target_rise_time(time, fixed_target)
set_time  = ost.target_set_time(time, fixed_target)
    
if rise_time.jd != '--' and set_time.jd != '--':
    rise_time = rise_time + 5*u.minute
    set_time  = set_time  - 5*u.minute
else:
    rise_time = sunrise_tonight
    set_time  = sunset_tonight

start = np.max([sunset_tonight, set_time])
end   = np.min([sunrise_tonight, rise_time])
  

###
#   Sky chart
#   

#   Observation window
time_window = start + (end - start) * np.linspace(0, 1, 10)

#   Make figure
fig = plt.figure(figsize=(11,7))

#   Plot sky chart
plot_sky(
    fixed_target, 
    ost, 
    time_window, 
    #style_kwargs={'color':colors[i-jold]},
    ) 

#   Position of the legend
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), shadow=True) 

#   Save figure
t            = time.replicate(format='iso', copy=True)
t.out_subfmt = 'date'
plt.savefig(
    'sky-chart_'+name_obj.replace(' ','_')+'_'+str(t)+'.pdf'
    )
#plt.show()
plt.close()


###
#   Altitude vs. obs. time plot
#
    
#   Time between Sunrise and Sunset
dark_time = sunrise_tonight - sunset_tonight

#   Set spacing for the plot
delta_t = np.linspace(0, dark_time.jd*24+2, 100)*u.hour
times   = sunset_tonight-1*u.hour + delta_t

#   Calculate altitude and azimuth frame
altazframe = coord.AltAz(obstime=times, location=location)

#   Calculate position of the Sun and Moon
sunaltazs  = get_sun(times).transform_to(altazframe)
moonaltazs = get_body("moon", times, location).transform_to(altazframe)

#   Calculate position of the object
objaltazs  = obj_coord.transform_to(altazframe)  

#   Make figure
fig = plt.figure(figsize=(11,7))

#   Activate support for Time objects in the plots 
time_support(format='iso', scale='utc')

#   Plot Moon
plt.plot(times, moonaltazs.alt, color='r', label='Moon')  

#   Set format for the time label
plt.gcf().autofmt_xdate() 

#   Plot object
plt.scatter(
    times, 
    objaltazs.alt, 
    c=objaltazs.az.deg, 
    label=name_obj, 
    lw=0, 
    s=8,
    )  

#   Highlight twilight (gray)
plt.fill_between(
    times,
    0, 
    90, 
    sunaltazs.alt < -0*u.deg, 
    color='0.5', 
    zorder=0,
    )  

#   Highlight dark time (black)
plt.fill_between(
    times, 
    0, 
    90, 
    sunaltazs.alt < -18*u.deg, 
    color='k', 
    zorder=0,
    )  

#   Color bar
plt.colorbar().set_label('Azimuth [deg]')  

#   Position of the legend
plt.legend(loc='upper left', shadow=True)  

#   X range
plt.xlim(sunset_tonight-1*u.hour, sunrise_tonight+1*u.hour)

#   X ticks
nhours = int(dark_time.jd*24+2)
plt.xticks(sunset_tonight + np.arange(0, nhours, 2)*u.hour)  

#   Y range
plt.ylim(0, 90)  

#   Y label 
plt.ylabel('Altitude [deg]')  

#   Save figure
t            = time.replicate(format='iso', copy=True)
t.out_subfmt = 'date'
plt.savefig(
    name_obj.replace(' ','_')+'_'+str(t)+'.pdf'
    )
#plt.show()
plt.close()

