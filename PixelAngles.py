# angle between pixels and the moon, sun and galactic altitude
#
# prior to use that script you need to install skyfield
# pip install skyfield
#
# ====================
import getopt
import sys
import matplotlib.pyplot as plt
import numpy as np
import yaml
from skyfield.api import E, N, load, wgs84
from skyfield.framelib import galactic_frame
from scipy import ndimage, misc
# read site coordinates
# Load Parameters
with open("input_params.in") as f:
    p = yaml.safe_load(f)

ts = load.timescale()
# time=ts.utc(2020, 1, 1, 10, 35, 7)
time = ts.now()


print(time)
print(time.utc_jpl())


def input(argv):
    Vfile = "undefined"
    Rfile = "undefined"
    try:
        opts, args = getopt.getopt(argv, "h:v:r:", ["help=", "vfile=", "rfile="])
    except getopt.GetoptError:
        pprint("test.py -v <Vfile> -r <Rfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("test.py -v <Vfile> -r <Rfile>")
            sys.exit()
        elif opt in ("-v", "--vfile"):
            Vfile = arg
        elif opt in ("-r", "--rfile"):
            Rfile = arg
    print("Johnson V file is ", Vfile)
    print("Johnson R file is ", Rfile)
    return Vfile , Rfile

Vfile, Rfile = input(sys.argv[1:])

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
#
#  LOAD IMAGES
#
#  CONVERT TO GREY
#
#
#

print(f"Moon altitude: {altm.degrees:.4f}")
print(f"Moon azimuth: {azim.degrees:.4f}")
print(f"Sun altitude: {alts.degrees:.4f}")
print(f"Sun azimuth: {azis.degrees:.4f}")

#for band, xzen, yzen, xpol, ypol, img in (
#    ("V", p["XzenithV"], p["YzenithV"], p["XpolarisV"], p["YpolarisV"], Vfile),
#    ("R", p["XzenithR"], p["YzenithR"], p["XpolarisR"], p["YpolarisR"], Rfile),):
#    print(f"Processing Johnson {band} camera...")
xzen=1012
yzen=758
xpol=1012
ypol=758
y, x = np.indices((1516, 2024))

# computing the distance to zenith in pixels
d = np.hypot(x - xzen, y - yzen)
z = p["Zconstant"] + p["Zslope"] * d + p["Zquad"] * d**2
#z = np.where(z > 90, np.nan, z)
azpol = np.arctan2(xpol - xzen,yzen - ypol)*180/np.pi
az = np.arctan2(-x + xzen,-y + yzen)*180/np.pi
az = (az - azpol)
az = np.where(az < 0, az + 360, az)

direction = here_now.from_altaz(alt_degrees=90 - z, az_degrees=az)
glat, glon, gdistance = direction.frame_latlon(galactic_frame)
galactic_lat = glat.degrees
theta_sun = direction.separation_from(sun_position).degrees
theta_moon = direction.separation_from(moon_position).degrees

# mask non sense zenith angles
az [z > 90] = np.nan
galactic_lat [z > 90] = np.nan
theta_sun [z > 90] = np.nan
theta_moon [z > 90] = np.nan
z [z > 90] = np.nan
#img [z > 90] = np.nan

plt.figure()
plt.imshow(np.round(theta_sun), cmap = 'rainbow')
plt.colorbar()
plt.title('Solar angle')

plt.figure()
plt.imshow(np.round(theta_moon), cmap = 'rainbow')
plt.colorbar()
plt.title('Moon angle')

plt.figure()
plt.imshow(np.round(galactic_lat), cmap = 'rainbow')
plt.colorbar()
plt.title('Galactic Latitude')

plt.figure()
plt.imshow(np.round(az), cmap = 'rainbow')
plt.colorbar()
plt.title('Azimuth')
y, x = np.indices((1516, 2024))
plt.show()
