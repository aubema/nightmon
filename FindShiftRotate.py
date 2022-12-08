#!/usr/bin/python3
# angle between pixels and the moon, sun and galactic altitude
#
# prior to use that script you need to install skyfield
# pip install skyfield
#
# ====================
import getopt
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from datetime import date, timedelta
from astropy.convolution import Box2DKernel, Gaussian2DKernel, convolve
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from numpy import inf
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from scipy import ndimage, stats
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from skyfield import almanac
from skyfield.api import E, N, Star, load, wgs84
from skyfield.data import hipparcos
from skyfield.framelib import galactic_frame

from toolbox import open_raw, to_grayscale


# find indices of pixels of a list of alt az coordinates
def find_close_indices(array1, array2, value1, value2):
    # 1 = azimuth, 2 = elevation, 3 = radius, all in degrees
    minindex = [0, 0]
    for n in range(np.shape(value1)[0]):
        dist = np.sqrt(
            ((array1 - value1[n]) * np.cos(np.pi * (array2 + value2[n]) / 2 / 180)) ** 2
            + (array2 - value2[n]) ** 2
        )
        min = np.nanmin(dist)
        imin = np.where(dist == min)
        minindex = np.append(minindex, imin)
    return minindex


def airmass(elevat):
    zenith_rad = (90 - elevat) * np.pi / 180
    AM = (
        1.002432 * ((np.cos(zenith_rad)) ** 2)
        + 0.148386 * (np.cos(zenith_rad))
        + 0.0096467
    ) / (
        np.cos(zenith_rad) ** 3
        + 0.149864 * (np.cos(zenith_rad) ** 2)
        + 0.0102963 * (np.cos(zenith_rad))
        + 0.000303978
    )
    return AM


# reading command line parameters
def input(argv):
    Sfile = "undefined"
    # reading V and R images names from the input parameters
    # default extinction values will be replaces if provided to the command line arguments
    # RGO, D. K. (1985). Atmospheric Extinction at the Roque de los Muchachos Observatory, La Palma.
    # extinc_r = 0.0547
    # extinc_v = 0.102
    extinc = 0.102
    try:
        opts, args = getopt.getopt(
            argv,
            "h:s:c:",
            ["help=", "sfile=", "cam="],
        )
    except getopt.GetoptError:
        print("FindShiftRotate.py -s <Sfile> -c <Cam>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("FindShiftRotate.py -s <Sfile> -c <Cam>")
            sys.exit()
        elif opt in ("-s", "--sfile"):
            Sfile = arg
        elif opt in ("-c", "--cam"):
            Cam = arg
    print("Sky image file is :", Sfile)
    print("Camera is :", Cam)
    return Sfile, Cam


# ================================================
# MAIN
# default Parameters
user = "aubema"
sberr = 0
mflag = "False"
FWHM = 3
norm = ImageNormalize(stretch=SqrtStretch())

elmin = 5  # set the minimum elevation
# selecting magnitude limits for reference stars in the simbad database
# limits allow to remove saturated stars in the images
limiti = (
    3.7  # limitting stars magnitude NEED TO BE LOWER THAN 6 but 3.7 is a magic number
)
limits = 1.2
# load command line parameters
Sfile, Cam = input(sys.argv[1:])

# open image
print("Reading images...")
Simg = open_raw(Sfile)

# Load Parameters
# read NightMon config file
configpath = "/home/" + user + "/cameraorientation_config"
with open(cofigpath) as f:
    p = yaml.safe_load(f)
Site = p["Site"]
nightmonconfigpath = "/home/" + user + "/nightmon_config"
with open(nightmoncofigpath) as g:
    q = yaml.safe_load(g)


ts = load.timescale()
# set observer position
eph = load("de421.bsp")
here = eph["earth"] + wgs84.latlon(
    p["Latitude"] * N, p["Longitude"] * E, elevation_m=p["Altitude"]
)
# create the grayscale images
RC = 1
GC = 1
BC = 1
Sgray = to_grayscale(Simg, [RC, GC, BC], normalize=False)
ny = np.shape(Sgray)[0]
nx = np.shape(Sgray)[1]
# set the minimum elevation to the limit of the image if required
if elmin < 90 - 0.98 * (ny / 2 * q["Zslope"] + (ny / 2) ** 2 * q["Zquad"]):
    elmin = 90 - 0.98 * (ny / 2 * q["Zslope"] + (ny / 2) ** 2 * q["Zquad"])
print("Minimum elevation :", elmin)
# calculate maximum zenith angle
zemax = 90 - elmin
rnmax = (-q["Zslope"] + np.sqrt(q["Zslope"] ** 2 - 4 * q["Zquad"] * -zemax)) / (
    2 * q["Zquad"]
)

# creating elevation, azimith maps
print("Creating Azimuth and elevation maps...")
y, x = np.indices((ny, nx))
# computing the distance to zenith in pixels
d = np.hypot(x - nx / 2, y - ny / 2)
d[d < 0.5] = 0.5
z = p["Zslope"] * d + p["Zquad"] * d**2 + p["Zthird"] * d**3 + p["Zfourth"] * d**4
z[z < 0] = 0
# computing azimuth
az = np.arctan2(-x + nx / 2, -y + ny / 2) * 180 / np.pi
az = np.where(az < 0, az + 360, az)
# compute elevation angle
el = 90 - z

# initialize np arrays for backgrounds and stars fields
Sbkg = np.zeros([ny, nx])
Sstars = np.zeros([ny, nx])

# setting polaris position in the image (depend on the site, the time,
# the orientation of the camera and the camera model)
if Cam == "A":
    xpoli = p["XpolarisA"]
    ypoli = p["YpolarisA"]
elif Cam == "B":
    xpoli = p["XpolarisB"]
    ypoli = p["YpolarisB"]

# file names provided by NightMon are formatted in the following way
# YYYYY-MM-DD_hh-mm-ss_V.dng , the V letter represents the filter and can thus be
# replaced by R. Time and date are in UTC.
# time=ts.utc(2020, 9, 22, 4, 22, 53)
# time = ts.now()
# yesterday = time - timedelta(days = 1)
datearr = Sfile.replace("_", "-").split("-")
time = ts.utc(
    int(datearr[0]),
    int(datearr[1]),
    int(datearr[2]),
    int(datearr[3]),
    int(datearr[4]),
    int(datearr[5]),
)

here_now = here.at(time)

# creating star map with the SIMBAD stars database
simbadpath = "/home/" + user + "/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv"
ds = pd.read_csv(simbadpath, header=0, sep=";")

# locating Polaris
polaris = ds.loc[ds["identifier"] == "* alf UMi"]
radpolaris = polaris["coord1_ICRS,J2000/2000_"]
coordspolaris = radpolaris.to_numpy(dtype="str")
cpolaris = coordspolaris[0].split(" ")
etoilepol = Star(
    ra_hours=(float(cpolaris[0]), float(cpolaris[1]), float(cpolaris[2])),
    dec_degrees=(float(cpolaris[3]), float(cpolaris[4]), float(cpolaris[5])),
)
epol = here_now.observe(etoilepol).apparent()
alt_star, azi_star, distance = epol.altaz()
alt_pol = alt_star.degrees
azi_pol = azi_star.degrees
altpol = np.empty([2], dtype=float)
azipol = np.empty([2], dtype=float)
azipol[0] = azi_pol
altpol[0] = alt_pol
azipol[1] = azi_pol
altpol[1] = alt_pol
# position of Polaris on the image grid
print("Position Polaris on the image grid...")
polindex = find_close_indices(az, el, azipol, altpol)
polshape = int(np.shape(polindex)[0])
polshape = int(polshape / 2)
polindex = polindex.reshape(polshape, 2)
polindex = np.delete(polindex, (0), axis=0)
polshape = int(np.shape(polindex)[0])
print("Polaris located at : ", polindex[0, 1], polindex[0, 0])

stars_selected = ds[(ds["MagV"] < limiti) & (ds["MagV"] > limits)]
coordsradec = stars_selected["coord1_ICRS,J2000/2000_"]
coords = coordsradec.to_numpy(dtype="str")
identity = stars_selected["identifier"]
iden = identity.to_numpy(dtype="str")
magnitudev = stars_selected["MagV"]
magv = magnitudev.to_numpy(dtype="float")
magnituder = stars_selected["MagR"]
magr = magnituder.to_numpy(dtype="float")
# find pixel positions of the SIMBAD stars on the array grid corresponding to image size
print("Position simbad stars on the image grid...")
altstar = np.zeros(np.shape(coords)[0])
azistar = np.zeros(np.shape(coords)[0])
# create np array for stars
for i in range(np.shape(coords)[0]):
    posstar = coords[i].split(" ")
    rah = float(posstar[0])
    ram = float(posstar[1])
    ras = float(posstar[2])
    ded = float(posstar[3])
    dem = float(posstar[4])
    des = float(posstar[5])
    etoile = Star(ra_hours=(rah, ram, ras), dec_degrees=(ded, dem, des))
    etoi = here_now.observe(etoile).apparent()
    alt_star, azi_star, distance = etoi.altaz()
    alt_star = alt_star.degrees
    azi_star = azi_star.degrees
    azistar[i] = azi_star
    altstar[i] = alt_star
# keep stars above horizon limitz = 0.95 * ny/2 * p["Zslope"]
azistar = np.delete(azistar, np.where(altstar < elmin))
magv = np.delete(magv, np.where(altstar < elmin))
magr = np.delete(magr, np.where(altstar < elmin))
iden = np.delete(iden, np.where(altstar < elmin))
altstar = np.delete(altstar, np.where(altstar < elmin))
brightest = np.amin(magv)
# find position of simbad stars on the image array
index = find_close_indices(az, el, azistar, altstar)
ishape = int(np.shape(index)[0])
ishape = int(ishape / 2)
index = index.reshape(ishape, 2)
index = np.delete(index, (0), axis=0)
ishape = int(np.shape(index)[0])
# index first column = y and second = x
print("Number of SIMBAD reference stars above", elmin, "degrees :", ishape)

deltax = np.empty([3], dtype=float)
deltay = np.empty([3], dtype=float)
angle = np.empty([3], dtype=float)
imag = Sgray
shiftx = 0
shifty = 0
for i in range(3):
    shiftx = polindex[0, 1] - xpoli
    shifty = polindex[0, 0] - ypoli
    theta1 = np.arctan2(ypoli - ny / 2, xpoli - nx / 2) * 180 / np.pi
    theta2 = np.arctan2(polindex[0, 0] - ny / 2, polindex[0, 1] - nx / 2) * 180 / np.pi
    theta = theta1 - theta2
    angpol1 = np.arctan2(xpoli - nx / 2, ny / 2 - ypoli) * 180 / np.pi
    angpol = angpol1 + azi_pol
    imag = ndimage.shift(imag, [shifty, shiftx], mode="nearest")
    padX = [imag.shape[1] - polindex[0, 1], polindex[0, 1]]
    padY = [imag.shape[0] - polindex[0, 0], polindex[0, 0]]
    imgP = np.pad(imag, [padY, padX], "constant")
    imag = ndimage.rotate(imgP, theta, reshape=False, mode="nearest")
    imag = imag[padY[0] : -padY[1], padX[0] : -padX[1]]
    deltax[i] = shiftx
    deltay[i] = shifty
    angle[i] = theta
    # find background
    sigma_clip = SigmaClip(sigma=2.0)
    bkg_estimator = MedianBackground()
    # background image
    bkg = Background2D(
        imag,
        (9, 9),
        filter_size=(int(FWHM), int(FWHM)),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )
    imbkg = bkg.background
    # image without background
    imstars = imag - imbkg

    mean, median, std = sigma_clipped_stats(imstars, sigma=5.0)
    daofind = DAOStarFinder(fwhm=FWHM, threshold=5.0 * std)
    sources = daofind(imstars)

    # keep stars with elevation > elmin degrees
    for col in sources.colnames:
        sources[col].info.format = "%.8g"  # for consistent table output
    positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))
    Flux = sources["flux"]
    rn = np.hypot(positions[:, 0] - nx / 2, positions[:, 1] - ny / 2)
    sources = sources[rn < rnmax]
    positions = positions[rn < rnmax]
    Flux = Flux[rn < rnmax]
    maxflux = Flux.max()
    sources = sources[Flux > maxflux / 2.5 ** (limiti - brightest)]
    positions = positions[Flux > maxflux / 2.5 ** (limiti - brightest)]
    Flux = Flux[Flux > maxflux / 2.5 ** (limiti - brightest)]

    # find most probable Polaris
    dtopol = np.hypot(
        polindex[0, 1] - positions[:, 0], polindex[0, 0] - positions[:, 1]
    )
    mintopol = np.nanmin(dtopol)
    indtopol = np.where(dtopol == mintopol)
    xpoli = float(positions[indtopol, 0])
    ypoli = float(positions[indtopol, 1])
    # positions[:, 1] = ny - positions[:, 1]
    apertures = CircularAperture(positions, r=4.0)


print("Shift x : ", deltax)
print("Shift y : ", deltay)
print("Angle : ", angle)
print("Copy that in your nightmon_config file")

StarMatch = np.zeros([ishape, 9])
n = 0
nistars = ishape
# searching for correspondance between stars in simbad and found stars in image
dstar = np.zeros(np.shape(positions)[0])
for ns in range(ishape):
    # for npos in range(np.shape(positions)[0]-1):
    #     dstar[npos] = np.hypot(positions[npos,0]-index[ns,1],positions[npos,1]-index[ns,0])
    dstar = (
        np.hypot(positions[:, 0] - index[ns, 1], positions[:, 1] - index[ns, 0]) ** 2
        / Flux
    )
    dmin = np.amin(dstar)
    dmin_index = dstar.argmin()
    if dmin < 2000000:
        StarMatch[n, 0] = index[ns, 1]
        StarMatch[n, 1] = index[ns, 0]
        StarMatch[n, 2] = positions[dmin_index, 0]  # noeuds[0]
        StarMatch[n, 3] = positions[dmin_index, 1]  # noeuds[1]
        StarMatch[n, 4] = dmin
        StarMatch[n, 5] = (
            positions[dmin_index, 0] - index[ns, 1]
        )  # noeuds[0] - index[ns,1]
        StarMatch[n, 6] = (
            positions[dmin_index, 1] - index[ns, 0]
        )  # noeuds[1] - index[ns,0]
        StarMatch[n, 7] = magv[ns]
        StarMatch[n, 8] = Flux[dmin_index]
        n = n + 1
StarMatch[np.isnan(StarMatch)] = 0
StarMatch = np.delete(StarMatch, np.where(StarMatch == 0), axis=0)
avggap, mediangap, stdgap = sigma_clipped_stats(StarMatch[:, 4], sigma=2.0)
print("Average distance between nearest star :", avggap, "+/-", stdgap)
StarMatch = np.delete(
    StarMatch, np.where(StarMatch[:, 4] > avggap + 1 * stdgap), axis=0
)
StarMatch = np.delete(StarMatch, np.where(StarMatch[:, 8] == 0), axis=0)

title = "Stars correspondance"
plt.figure()
plt.plot(StarMatch[:, 2], StarMatch[:, 3], "or", markersize=2)
plt.plot(StarMatch[:, 0], StarMatch[:, 1], "ok", markersize=2)
plt.ylim(ny, 0)
plt.xlim(0, nx)
plt.title(title)
plt.legend(["Simbad reference stars", "Detected stars"])
plt.show()
