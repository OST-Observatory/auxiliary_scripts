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
#obs_time = '2021-03-30 00:00'
obs_time = ''

###
#   Path to the file with the cluster data
#
object_file = 'cluster.dat'


############################################################################
####                            Libraries                               ####
############################################################################

import sys

import numpy as np

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.table import Table
from astropy.table import vstack
from astropy.time import Time
import astropy.coordinates as coord
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon

from astroplan import Observer

############################################################################
####                               Main                                 ####
############################################################################

###
#   Read the data
#
cluster = Table.read(object_file, format='ascii')


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

#   Time between Sunrise and Sunset
dark_time = sunrise_tonight - sunset_tonight

#   Number of hours
nhours = int(dark_time.jd*24+2)


###
#   Coordinates of the objects
#
coord_cluster = coord.SkyCoord(
    ra=cluster['ra']*u.deg, 
    dec=cluster['dec']*u.deg, 
    frame='icrs',
    )


###
#   Add constellation
#

#   Find constellation
constellation = coord.get_constellation(coord_cluster, short_name=True)

#   Add 
cluster.add_column(constellation, name='constellation')


###
#   Sort cluster according to airmass - select suitable cluster
#

#   Create copyies of the cluster Table\
backup  = Table(cluster, copy=True)

#   Loop over all observing blocks (OBs) during the night -> 1OB = 1h
for i in range(0, nhours):
    time = sunset_tonight + i*u.hour
    
    #   Initialize system
    AlAz = coord.AltAz(
        location=location, 
        obstime=time
        )
    
    #   Transform coordinates to the Altitude-Azimuth system
    cluster_alt_az = coord_cluster.transform_to(AlAz)

    #   Add coordinates to the cluster Table
    cl_clea = Table(backup, copy=True)

    time.format = 'iso'
    cl_clea.add_columns(
        [time, i, cluster_alt_az.alt, cluster_alt_az.az],
        names=['time', 'slot', 'alt', 'az']
    )
    
    #   Calculate position of the Moon
    moon = get_moon(time).transform_to(AlAz)

    #   Calculate the separation between moon and target
    moon_separation = np.round(cluster_alt_az.separation(moon).degree, 1)
    
    #   Add moon separation to the cluster Table
    cl_clea.add_column(moon_separation, name='moon_d')
    
    #   Airmass
    airmass = np.round(cluster_alt_az.secz.value, 2)
    
    #   Add airmass to the cluster Table
    cl_clea.add_column(airmass, name='secz')
    
    #   Filter by altitude above horizon
    mask    = cl_clea['alt'] >= 1. 
    cl_alt  = cl_clea[mask]
    
    #   Filter by moon distance
    mask    = cl_alt['moon_d'] >= 90. 
    cl_alt  = cl_alt[mask]
    
    #   Filter by number of objects 
    mask    = cl_alt['n_stars'] >= 10. 
    cl_alt  = cl_alt[mask]

    
    
    #   secz = 1.0-1.5
    mask  = cl_alt['secz'] < 1.5
    if i == 0:
        secz1 = cl_alt[mask]
    else:
        secz1 = vstack([secz1, cl_alt[mask]])
    
    #   secz = 1.5-2.0
    mask  = cl_alt['secz'] >= 1.5
    dummy = cl_alt[mask]
    mask  = dummy['secz'] < 2.0
    if i == 0:
        secz2 = dummy[mask]
    else:
        secz2 = vstack([secz2, dummy[mask]])
        
    #   secz = 2.0-2.5
    mask  = cl_alt['secz'] >= 2.0
    dummy = cl_alt[mask]
    mask  = dummy['secz'] < 2.5
    if i == 0:
        secz3 = dummy[mask]
    else:
        secz3 = vstack([secz3, dummy[mask]])
        
    #   secz = 2.5-3.0
    mask  = cl_alt['secz'] >= 2.5
    dummy = cl_alt[mask]
    mask  = dummy['secz'] < 3.0
    if i == 0:
        secz4 = dummy[mask]
    else:
        secz4 = vstack([secz4, dummy[mask]])
        
    #   secz = 3.0-3.5
    mask  = cl_alt['secz'] >= 3.0
    dummy = cl_alt[mask]
    mask  = dummy['secz'] < 3.5
    if i == 0:
        secz5 = dummy[mask]
    else:
        secz5 = vstack([secz5, dummy[mask]])
        
    #   secz = 3.5-4.0
    mask  = cl_alt['secz'] >= 3.5
    dummy = cl_alt[mask]
    mask  = dummy['secz'] < 4.0
    if i == 0:
        secz6 = dummy[mask]
    else:
        secz6 = vstack([secz6, dummy[mask]])
        
    #   secz >= 4.0
    mask  = cl_alt['secz'] >= 4.0
    if i == 0:
        secz7 = cl_alt[mask]
    else:
        secz7 = vstack([secz7, cl_alt[mask]])
    
print('### Airmass < 1.5 ###')
print(secz1.group_by('slot'))
print('')

print('### 1.5 <= Airmass < 2.0 ###')
print(secz2.group_by('slot'))
print('')

print('### 2.0 <= Airmass < 2.5 ###')
print(secz3.group_by('slot'))
print('')
    
print('### 2.5 <= Airmass < 3.0 ###')
print(secz4.group_by('slot'))
print('')

print('### 3.5 <= Airmass < 3.5 ###')
print(secz5.group_by('slot'))
print('')

print('### 3.5 <= Airmass < 4.0 ###')
print(secz6.group_by('slot'))
print('')

print('### Airmass >= 4.0 ###')
print(secz7.group_by('slot'))
print('')
