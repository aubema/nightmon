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

for band, xzen, yzen, xpol, ypol, img in (
    ("V", p["XzenithV"], p["YzenithV"], p["XpolarisV"], p["YpolarisV"], Vfile),
    ("R", p["XzenithR"], p["YzenithR"], p["XpolarisR"], p["YpolarisR"], Rfile),):
    print(f"Processing Johnson {band} camera...")

    y, x = np.indices((1516, 2024))

    # computing the distance to zenith in pixels
    d = np.hypot(x - xzen, y - yzen)
    z = p["Zconstant"] + p["Zslope"] * d + p["Zquad"] * d**2
    #z = np.where(z > 90, np.nan, z)
    azpol = np.arctan2(xpol - xzen,yzen - ypol)*180/np.pi
    az = np.arctan2(-x + xzen,-y + yzen)*180/np.pi
    az = (az - azpol)
    az = np.where(az < 0, az + 360, az)

    # shift images
    z = ndimage.shift(z, [ (758-yzen) , (1012-xzen) ], mode='nearest')
    img = ndimage.shift(img, [ (758-yzen) , (1012-xzen) ], mode='nearest' )

    # rotate images
    z = ndimage.rotate(z, -azpol, reshape=False, mode='nearest')
    img = ndimage.rotate(img, -azpol, reshape=False, mode='nearest')

    # mask non sense zenith angles
    z [z > 90] = np.nan
    img [z > 90] = np.nan

    plt.figure()
    plt.imshow(np.round(img), cmap = 'rainbow')
    plt.colorbar()
    plt.title('Sky Image : Johnson ' + band)
plt.show()
