# angle between pixels and the moon, sun and galactic altitude
#
# prior to use that script you need to install skyfield
# pip install skyfield
#
# ====================
import getopt
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import yaml
import scipy
from skyfield.api import E, N, load, wgs84
from skyfield.framelib import galactic_frame
from scipy import ndimage, misc
from matplotlib.pyplot import imread

# read site coordinates
# Load Parameters
home = os.path.expanduser('~')
with open(home + "/nightmon_config") as f:

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
# TODO: open raw images Vfile and Rfile

# TODO: Convert raw to grayscale with coeffs RC GC and BC output array aray are Vgray and Rgray
RC = p["R2GRAYCOEF"]
GC = p["G2GRAYCOEF"]
BC = p["B2GRAYCOEF"]
Vgray =
Rgray =
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

for band, xzen, yzen, xpol, ypol, img, outf in (
    ("V", p["XzenithV"], p["YzenithV"], p["XpolarisV"], p["YpolarisV"], Vgray, "VImage.npy"),
    ("R", p["XzenithR"], p["YzenithR"], p["XpolarisR"], p["YpolarisR"], Rgray, "RImage.npy"),):
    print(f"Processing Johnson {band} camera...")
    #
    #
    #
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
    az = ndimage.shift(az, [ (758-yzen) , (1012-xzen) ], mode='nearest')
    z = ndimage.shift(z, [ (758-yzen) , (1012-xzen) ], mode='nearest')
    imag = ndimage.shift(imag, [ (758-yzen) , (1012-xzen) ], mode='nearest' )

    # rotate images
    az = ndimage.rotate(az, -azpol, reshape=False, mode='nearest')
    z = ndimage.rotate(z, -azpol, reshape=False, mode='nearest')
    imag = ndimage.rotate(imag, -azpol, reshape=False, mode='nearest')

    # mask non sense zenith angles
    imag [z > 90] = np.nan
    az [z > 90] = np.nan
    z [z > 90] = np.nan
    np.save(outf, imag)

    plt.figure()
    plt.imshow(imag, cmap = 'rainbow')
    plt.colorbar()
    plt.title('Sky Image : Johnson ' + band)

    plt.figure()
    plt.imshow(az, cmap = 'rainbow')
    plt.colorbar()
    plt.title('Azimuth angle' + band)

plt.figure()
plt.imshow(z, cmap = 'rainbow')
plt.colorbar()
plt.title('Zenith angle')


#plt.show()
