# Given the current date/time and location, what is the APPROXIMATE 
# RA/Dec that corresponds with a certain altitude/azimuth? 
# Hardcode your latitude, longitude, and elevation in this
# script for ease of use while observing.
# 
# Limitations: does not take into account pressure, air temperature
# for refraction corrections. The output RA/Dec is APPROXIMATE
# for positioning a very wide-angle field of view, it is not meant
# for precision astrometry.
# 
# Use:
#   python altaz_to_radec.py ALT AZ
# Example:
#   python altaz_to_radec.py 41.25 359.42
#   Will tell the RA and Dec associated with altitude 41.0 deg, azimuth 0.0 deg,
#   which should be (approximately) Polaris from Salt Lake City.
#   Polaris:  RA = 02:31:48.7 hr, Dec = +89:15:51.3 deg

from astropy import units as u
import sys
import numpy as np
from astropy.time import Time
import datetime
from astropy import coordinates as coords
from astropy.coordinates import AltAz,SkyCoord

########
# HARDCODE YOUR CURRENT LATITUDE, LONGITUDE, AND ELEVATION HERE.
# Comment/uncomment sets of measurements.

# Westminster College: 40:44:16 N,111:52:16 W, elev 1344 m
locname='Westminster College'
lat=40.73778     # Degrees North
lon=-111.871111  # Degrees East
elev=1344        # meters

# Read in altitude and azimuth 
if len(sys.argv)<3:
    print('Usage:\npython altaz_to_radec.py ALT AZ \n')
    exit()
alt=np.float(sys.argv[1])
az=np.float(sys.argv[2])

########
# Get current UTC date/time from system. 
localnow=datetime.datetime.now()
now=datetime.datetime.utcnow()

print '============================================================='
print 'Local time: ',localnow.isoformat()
print 'UTC time:   ',now.isoformat()
print 'Confirm before trusting: MDT is UTC -06:00, MST is UTC -07:00'
print '============================================================='

########
# Convert Alt/Az to RA and Dec (J2000)
loc = coords.EarthLocation(lat = lat*u.deg,lon = lon*u.deg, height = elev*u.m)
t = Time(now,scale='utc') # Same result if us "t" or "now" for obstime, meh?
co = SkyCoord(alt = alt*u.deg, az = az*u.deg, obstime = now, frame = 'altaz',location=loc)

print 'Location, ',locname+': ',loc.lat,loc.lon,loc.height
print 'Altitude of ',alt,' degrees, azimuth of ',az,' degrees corresponds to (J2000):'
print co.fk5.to_string('hmsdms')
print '============================================================='
