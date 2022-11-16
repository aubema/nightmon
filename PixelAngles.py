# angle between pixels and the moon, sun and galactic altitude
#
# prior to use that script you need to install skyfield
# pip install skyfield
#
# ====================
import os

import matplotlib.pyplot as plt
import numpy as np
import yaml
from skyfield.api import E, N, load, wgs84
from skyfield.framelib import galactic_frame

# read site coordinates
# Load Parameters
home = os.path.expanduser("~")
with open(home + "/nightmon_config") as f:
    p = yaml.safe_load(f)

ts = load.timescale()
# time=ts.utc(2020, 1, 1, 10, 35, 7)
time = ts.now()

print(time)
print(time.utc_jpl())

# astronomical objects and ephemerides
eph = load("de421.bsp")
here = eph["earth"] + wgs84.latlon(
    p["LATITUDE"] * N, p["LONGITUDE"] * E, elevation_m=p["ELEVATION"]
)
here_now = here.at(time)
moon_position = here_now.observe(eph["moon"]).apparent()
sun_position = here_now.observe(eph["sun"]).apparent()
alts, azis, ds = sun_position.altaz()
altm, azim, dm = moon_position.altaz()

print(f"Moon altitude: {altm.degrees:.4f}")
print(f"Moon azimuth: {azim.degrees:.4f}")
print(f"Sun altitude: {alts.degrees:.4f}")
print(f"Sun azimuth: {azis.degrees:.4f}")
# half size of superpixel image (i.e. 1/4 of the original)
xzen = 1012
yzen = 758
xpol = 1012
ypol = 400
y, x = np.indices((1516, 2024))

# computing the distance to zenith in pixels
d = np.hypot(x - xzen, y - yzen)
z = p["Zconstant"] + p["Zslope"] * d + p["Zquad"] * d**2
# z = np.where(z > 90, np.nan, z)
azpol = np.arctan2(xpol - xzen, yzen - ypol) * 180 / np.pi
az = np.arctan2(-x + xzen, -y + yzen) * 180 / np.pi
az = az - azpol
az = np.where(az < 0, az + 360, az)

direction = here_now.from_altaz(alt_degrees=90 - z, az_degrees=az)
glat, glon, gdistance = direction.frame_latlon(galactic_frame)
galactic_lat = glat.degrees
theta_sun = direction.separation_from(sun_position).degrees
theta_moon = direction.separation_from(moon_position).degrees

# mask non sense zenith angles
mask = z > 90
az[mask] = np.nan
galactic_lat[mask] = np.nan
theta_sun[mask] = np.nan
theta_moon[mask] = np.nan
z[mask] = np.nan
# img [mask] = np.nan

np.save("SolarAngle", theta_sun)
np.save("MoonAngle", theta_moon)
np.save("GalacticLattitude", galactic_lat)
np.save("AzimuthAngle", az)
np.save("ZenithAngle", z)


plt.figure()
plt.imshow(np.round(theta_sun), cmap="rainbow")
plt.colorbar()
plt.title("Solar angle")

plt.figure()
plt.imshow(np.round(theta_moon), cmap="rainbow")
plt.colorbar()
plt.title("Moon angle")

plt.figure()
plt.imshow(np.round(galactic_lat), cmap="rainbow")
plt.colorbar()
plt.title("Galactic Latitude")

plt.figure()
plt.imshow(np.round(az), cmap="rainbow")
plt.colorbar()
plt.title("Azimuth")

#plt.show()
