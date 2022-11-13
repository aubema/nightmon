# Ephemerides for NightMon
# Maintain a database of important ephemierides parameters for
# operation and analysis of the NightMon data
# copyright Martin Aube 2022
#
# prior to use that script you need to install skyfield
# pip install skyfield
#
# ====================
import datetime as dt
import yaml
from pytz import timezone
from skyfield import almanac
from skyfield.api import N, W, E, wgs84, load
# read site coordinates
# Load Parameters
with open("input_params.in") as f:
    p = yaml.safe_load(f)
lat = p['LATITUDE']
lon = p['LONGITUDE']
elev = p['ELEVATION']
# Calculate twilights
# zone = timezone('GMT')
zone = timezone('EST')
now = zone.localize(dt.datetime.now())
print(now, zone)
midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
next_midnight = midnight + dt.timedelta(days=1)

ts = load.timescale()
t0 = ts.from_datetime(midnight)
t1 = ts.from_datetime(next_midnight)
eph = load('de421.bsp')
bluffton = wgs84.latlon(lat * N, lon * E)
f = almanac.dark_twilight_day(eph, bluffton)
times, events = almanac.find_discrete(t0, t1, f)

previous_e = f(t0).item()
for t, e in zip(times, events):
    tstr = str(t.astimezone(zone))[:16]
    if previous_e < e:
        print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')
    else:
        print(tstr, ' ', almanac.TWILIGHTS[previous_e], 'ends')
    previous_e = e
