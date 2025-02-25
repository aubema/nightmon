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
from scipy.optimize import curve_fit
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
    Ifile = "undefined"
    # reading V and R images names from the input parameters
    # default extinction values will be replaces if provided to the command line arguments
    # RGO, D. K. (1985). Atmospheric Extinction at the Roque de los Muchachos Observatory, La Palma.
    # extinc_r = 0.0547
    # extinc_v = 0.102
    Extinc = 0
    Model = "RpiHQ"
    Cam = "A"
    Band = "JV"
    Calmet = "0"  # other option is fixed
    Zpoint = 0  # this value is only useful when selecting fixed calibration method otherwise ignored
    ExpF = 1. # this is the exposure correction factor compared to the exposure used for calibration. 
    Flip = 0
    # E.g. if the exposure is 12000 and the exposure used for calibration is 120000000, then exposure 
    # factor will be 12000/120000000 = 1/10000 = 0.0001
    try:
        opts, args = getopt.getopt(
            argv,
            "h:i:d:b:e:m:k:z:f:",
            [
                "help=",
                "ifile=",
                "dfile=",
                "extinc=",
                "band=",
                "model=",
                "calib=",
                "zerop=",
                "flip=",
            ],
        )
    except getopt.GetoptError:
        print(
            "ProcessNighMon.py -i <Ifile> -d <Dfile> -b <Band> -e <Extinc> -m <Model> -k <Calibration method> -z <Zeropoint> -f <1>"
        )
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(
                "ProcessNighMon.py -i <Ifile> -d <Dfile> -b <Band> -e <Extinc> -m <Model> -k <Calibration method> -z <Zeropoint> -f <1>"
            )
            sys.exit()
        elif opt in ("-i", "--ifile"):
            Ifile = arg
        elif opt in ("-d", "--dfile"):
            Dfile = arg
        elif opt in ("-b", "--band"):
            Band = arg
        elif opt in ("-e", "--extinc"):
            Extinc = arg
        elif opt in ("-m", "--model"):
            Model = arg
        elif opt in ("-k", "--calib"):
            Calmet = arg
        elif opt in ("-z", "--zerop"):
            Zpoint = arg
        elif opt in ("-f", "--flip"):
            Flip = arg
    print("Sky image file is :", Ifile)
    print("Dark frame file is :", Dfile)
    print("Band is :", Band)
    return Ifile, Dfile, Band, Extinc, Model, Calmet, Zpoint, Flip


def fit_func(x, a):
    # Curve fitting function
    return a * x


# ================================================
# MAIN
# default Parameters
user = "sand"
path = "/home/" + user + "/"
sberr = 0
calsb = 0
corcoef = 0
mflag = "False"
FWHM = 7
zeropoint = np.nan
Zpoint = 0
ZpointConfig = 0
max_cloud_cover = 4
Calmet = "0"
Calmetconfig = "0"
norm = ImageNormalize(stretch=SqrtStretch())

elmin = 15  # set the minimum elevation
# selecting magnitude limits for reference stars in the simbad database
# limits allow to remove saturated stars in the images
limiti = (
    4.7  # limitting stars magnitude NEED TO BE LOWER THAN 6 but 3.7 is a magic number
)
limits = -1
# load command line parameters
Ifile, Dfile, Band, Extinc, Model, Calmet, Zpoint, Flip = input(sys.argv[1:])
# Flip = 1 only for the bottom camera (camera C) in the NightMon-eco device
k = float(Extinc)
# TOTO : IL FAUT DEPLACER CECI ET LES EPHEMERIDES DANS LA BOUCLE DES FILTRES
# files names provided by NightMon are formatted in the following way
# YYYYY-MM-DD_hh-mm-ss_V.dng , the V letter represents the filter and can thus be
# replaced by R. Time and date are in UTC.
# time=ts.utc(2020, 9, 22, 4, 22, 53)
# time = ts.now()
# yesterday = time - timedelta(days = 1)
datearr = Ifile.replace("_", "-").replace(".", "-").split("-")
ts = load.timescale()
time = ts.utc(
    int(datearr[0]),
    int(datearr[1]),
    int(datearr[2]),
    int(datearr[3]),
    int(datearr[4]),
    int(datearr[5]),
)
Cam = datearr[6]

print("Camera is :", Cam)
Tint = int(datearr[7])
Gain = int(datearr[8])
# determine the R, G, B coefficients according to the band and the camera model
if Model == "A7S":
    if Band == "JV":
        RC = 0.65
        GC = 1.0
        BC = 0.4
    elif Band == "JR":
        RC = 1.0
        GC = 0.7
        BC = -1.0
    elif Band == "JB":
        RC = -0.65
        GC = 0.4
        BC = 1.0
    Zslope = 1.13163209e-01
    Zquad = 4.01469526e-06
    Zthird = -1.94853987e-08
    Zfourth = 4.63754847e-13
elif Model == "RpiHQ":
    if Band == "JV":
        RC = 1.0
        GC = 1.0
        BC = -0.2
    elif Band == "JR":
        RC = 1.0
        GC = 0.6
        BC = -0.8
    elif Band == "JB":
        RC = 0.1
        GC = -0.3
        BC = 1.0
    Zslope = 1.1362e-01 #1.13538e-01
    Zquad = 1.e-11 #1.0e-06
    Zthird = 0.
    Zfourth = 0.
elif Model == "RpiHQ-JFilters":
    if Band == "JV":
        RC = 1
        GC = 1
        BC = 0.6
    elif Band == "JR":
        RC = 1
        GC = 1
        BC = 1
    elif Band == "JB":
        RC = 1
        GC = 1
        BC = 1
    Zslope = 1.167e-01
    Zquad = 1.0e-09
    Zthird = 0
    Zfourth = 0
# open images
print("Process raw image: ", path + Ifile)
Simg = open_raw(path + Ifile)
print("Looking for image: ", path + Ifile)
Dimg = open_raw(Dfile)
# read site coordinates
# Load Parameters
configpath = "/home/" + user + "/nightmon_config"

with open(configpath) as f:

    p = yaml.safe_load(f)
Site = p["Site"]
Calmetconfig = p["Calmet"]
DefaultItime = p["Calib_itime"]
DefaultGain = p["Calib_gain"]
if Cam == "A":
    deltax = p["ShiftxA"]
    deltay = p["ShiftyA"]
    angle = p["AngleA"]
    if Band == "JB":
       ZpointConfig = p["ZpointA-JB"]
    elif Band == "JV":
       ZpointConfig = p["ZpointA-JV"]  
    elif Band == "JR": 
       ZpointConfig = p["ZpointA-JR"]      
elif Cam == "C":
    deltax = p["ShiftxC"]
    deltay = p["ShiftyC"]
    angle = p["AngleC"]
    if Band == "JB":
       ZpointConfig = p["ZpointC-JB"]
    elif Band == "JV":
       ZpointConfig = p["ZpointC-JV"]  
    elif Band == "JR": 
       ZpointConfig = p["ZpointC-JR"] 
if (
    Zpoint == 0 and ZpointConfig != 0
):  # if no zero point is provided in argument then use the value of the config file
    Zpoint = ZpointConfig
slp = 10 ** (-0.4 * float(Zpoint))
if (
    Calmet == "0" and Calmetconfig != "0"
):  # if no calib method is provided in argument then use the value of the config file
    Calmet = Calmetconfig
if Calmet == "0":
    Calmet = "stars"  # if no calib method is provided either in argument or in the configfile then use stars
# load sky brightness extraction points
pointspath = "/home/" + user + "/points_list"

pt = open(pointspath, "r")
pt.readline()  # skip one line
tpt = []
apt = []
ept = []
for line in pt:
    words = line.split()
    tpt.append(words[0])
    apt.append(float(words[1]))
    ept.append(float(words[2]))
pt.close()
num_pts = np.shape(tpt)[0]
print("Camera model is :", Model)
print("Calibration method :", Calmet)
print("Extinction is :", Extinc)
print("Number of points to extract :", num_pts)
# set observer position
eph = load("de421.bsp")
here = eph["earth"] + wgs84.latlon(
    p["Latitude"] * N, p["Longitude"] * E, elevation_m=p["Altitude"]
)
# create the grayscale images
Sgray = to_grayscale(Simg, [RC, GC, BC], normalize=False)
Dgray = to_grayscale(Dimg, [RC, GC, BC], normalize=False)
if Flip == "1":
    print("Flipping image vertically...")
    Sgray = np.flipud(Sgray)
    Dgray = np.flipud(Dgray)
ny = np.shape(Sgray)[0]
nx = np.shape(Sgray)[1]
# set the minimum elevation to the limit of the image if required
if elmin < 90 - 0.98 * (ny / 2 * Zslope + (ny / 2) ** 2 * Zquad):
    elmin = 90 - 0.98 * (ny / 2 * Zslope + (ny / 2) ** 2 * Zquad)
print("Minimum elevation :", elmin)
# calculate maximum zenith angle
zemax = 90 - elmin
rnmax = (-Zslope + np.sqrt(Zslope**2 - 4 * Zquad * -zemax)) / (2 * Zquad)

# remove darks
window = 7  #  1 deg ~= 10
kernel = Box2DKernel(width=window, mode="integrate")
dark = convolve(Dgray, kernel)
dark_norm = convolve(Sgray, kernel)
dmean = np.mean(dark[5:105, 5:105], axis=(0, 1))
imean = np.mean(Sgray[5:105, 5:105], axis=(0, 1))
dnorm = imean / dmean
dark = dark * dnorm
Sgray = Sgray - dark

# creating elevation, azimith maps
print("Creating Azimuth and elevation maps...")
y, x = np.indices((ny, nx))
# computing the distance to zenith in pixels
d = np.hypot(x - nx / 2, y - ny / 2)
d[d < 0.5] = 0.5
z = Zslope * d + Zquad * d**2 + Zthird * d**3 + Zfourth * d**4
z[z < 0] = 0
# computing azimuth
az = np.arctan2(-x + nx / 2, -y + ny / 2) * 180 / np.pi
# az = az - azpol
az = np.where(az < 0, az + 360, az)
# solid angle in sq arc second
sec2 = (
    (Zslope + 2 * Zquad * d + 3 * Zthird * d**2 + 4 * Zfourth * d**3)
    * 180
    * 3600**2
    * np.sin((z) * np.pi / 180)
) / (d * np.pi)
# compute elevation angle
el = 90 - z
ImgAirm = airmass(el)
ImgAirm[z >= 90] = np.nan
# generate flat image
if Model == "A7S":
    flat = (
        1.82888e-42 * d**16
        - 8.718595e-39 * d**15
        + 1.602995e-35 * d**14
        - 1.179508e-32 * d**13
        - 2.969461e-30 * d**12
        + 1.192837e-26 * d**11
        - 7.981811e-24 * d**10
        + 3.211854e-22 * d**9
        + 2.756806e-18 * d**8
        - 1.992123e-15 * d**7
        + 7.347057e-13 * d**6
        - 1.632845e-10 * d**5
        + 2.220187e-08 * d**4
        - 1.759099e-06 * d**3
        + 7.142221e-05 * d**2
        - 0.001388396 * d
        + 1.00525
    )
else:
    flat = (
        0.9999999999992863
        + 0.0003775455182209281 * d
        - 5.775598992847673e-07 * d**2
        - 4.419838608862352e-09 * d**3
        + 2.3664435720648865e-11 * d**4
        - 4.1361998468625023e-14 * d**5
        + 2.2896450343971182e-17 * d**6
    )
# flat = flat * sec2[ny / 2, nx / 2] / sec2
flat[z > 90] = 1
# plt.figure()
# plt.imshow(flat, cmap="magma")
# plt.colorbar()

# correct flat field
Sgray = Sgray / flat
# initialize np arrays for backgrounds and stars fields
Vbkg = np.zeros([ny, nx])
Rbkg = np.zeros([ny, nx])
Sstars = np.zeros([ny, nx])
imag = Sgray
print(f"Processing Johnson {Band} camera...")


basename = time.utc_strftime("%Y-%m-%d")
outname = path + "calibrated_" + Band + "_" + basename + "_sky.csv"
timestamp = time.utc_strftime("%Y-%m-%dT%H:%M:%S")
baseout = time.utc_strftime("%Y-%m-%d_%H-%M-%S")
calname = path + Band + "_calibration_stars_" + baseout + ".csv"
here_now = here.at(time)
# calculation sun and moon positions
moon_position = here_now.observe(eph["moon"]).apparent()
sun_position = here_now.observe(eph["sun"]).apparent()
alts, azis, ds = sun_position.altaz()
altm, azim, dm = moon_position.altaz()
altm = altm.degrees
azim = azim.degrees
alts = alts.degrees
azis = azis.degrees
moonphase = almanac.moon_phase(eph, time).degrees
# setting moon flag to 1 if its elevation is positive
if altm > 0:
    mflag = "Thrue"
print(f"Moon altitude: {altm:.4f}")
print(f"Moon azimuth: {azim:.4f}")
print(f"Sun altitude: {alts:.4f}")
print(f"Sun azimuth: {azis:.4f}")
# create the data file if it do not exists
if os.path.exists(outname) == False:
    o = open(outname, "w")
    first_line = "Loc_Name;Band;CCD_X_Position;CCD_Y_Position;Azimuth_Angle;Altitude_Angle;Airmass;Ext_Coef;Date;Moon;Clouds_Cover;Background_Sky_Brightness;err,Zero_Point,Cor_Coef,Sun_Angle,Moon_Angle,Galactic_Latitude,Moon_Phase \n"
    second_line = " ; ;(pixel);(pixel);(deg);(deg); ; ; ; ;(oktas);(mag/arcsec^2); ;(mag); ;(deg);(deg);(deg);(deg) \n"
    o.write(first_line)
    o.write(second_line)
    o.close()
# creating star map with the SIMBAD stars database
# TODO: restore the two following lines
# ds = pd.read_csv(
#    "/home/sand/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv", header=0, sep=";"
csvpath = "/home/" + user + "/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv"
ds = pd.read_csv(csvpath, header=0, sep=";")
# stars_selected=ds[ds['MagR'] < limit]

# try to find Polaris and shift and rotate image accordingly only when northern hemisphere
if float(p["Latitude"]) > 0:
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

    # position of Pole on the image grid
    print("Position Pole on the image grid...")

    altpole = np.empty([2], dtype=float)
    azipole = np.empty([2], dtype=float)
    azipole[0] = 0.
    altpole[0] = p["Latitude"]
    azipole[1] = 0.
    altpole[1] = p["Latitude"]
    poleindex = find_close_indices(az, el, azipole, altpole)
    poleshape = int(np.shape(poleindex)[0])
    poleshape = int(poleshape / 2)
    poleindex = poleindex.reshape(poleshape, 2)
    poleindex = np.delete(poleindex, (0), axis=0)
    poleshape = int(np.shape(poleindex)[0])
    print("Pole ideally located at : ", poleindex[0, 1], poleindex[0, 0])
    print(altpole,azipole)
    shiftx = 0
    shifty = 0
    for i in range(1):
        shiftx = float(deltax)
        shifty = float(deltay)
        theta = float(angle)
        imag = ndimage.shift(imag, [shifty, shiftx], mode="nearest")
        padX = [imag.shape[1] - poleindex[0, 1], poleindex[0, 1]]
        padY = [imag.shape[0] - poleindex[0, 0], poleindex[0, 0]]
        imgP = np.pad(imag, [padY, padX], "constant")
        imag = ndimage.rotate(imgP, theta, reshape=False, mode="nearest")
        imag = imag[padY[0] : -padY[1], padX[0] : -padX[1]]
# find background
sigma_clip = SigmaClip(sigma=3.0)
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

# determine cloud cover only for z < 40 deg to exclude near horizon obstacles
cld_deg = 30
cld_pix = cld_deg / Zslope
imstars_tmp = np.copy(imstars)

# plt.figure()
# plt.title("imstars_tmp")
# plt.imshow(imstars_tmp, cmap="magma")
# plt.colorbar()
# plt.figure()
# plt.title("imstars")
# plt.imshow(imstars, cmap="magma")
# plt.colorbar()
# plt.show()


# imstars_tmp = imstars #/ imbkg
imstars_tmp[z > cld_deg] = np.nan
starsstd = np.nanstd(imstars_tmp)
starsmean = np.nanmean(imstars_tmp)
threshold = starsmean + 7 * starsstd
stars_binary = np.full((ny, nx), 0.0)
stars_full = np.full((ny, nx), 1.0)
stars_full[z > cld_deg] = 0
stars_full_small = stars_full[
    int(ny / 2 - cld_pix) : int(ny / 2 + cld_pix),
    int(nx / 2 - cld_pix) : int(nx / 2 + cld_pix),
]
stars_binary[imstars_tmp >= threshold] = 1.0
stars_binary[z > cld_deg] = 0
stars_binary_small = stars_binary[
    int(ny / 2 - cld_pix) : int(ny / 2 + cld_pix),
    int(nx / 2 - cld_pix) : int(nx / 2 + cld_pix),
]


plt.figure()
plt.title("stars_binary_small")
plt.imshow(stars_binary_small, cmap="magma")
plt.colorbar()


stars_full[z > cld_deg] = 0
# set cloud detection window to about 5 deg (51)
window = 171  #  1 deg ~= 17
kernel = Box2DKernel(width=window, mode="integrate")
stars_count = convolve(stars_binary_small, kernel)
stars_count = stars_count * window * window
plt.figure()
plt.title("stars_count1")
plt.imshow(stars_count, cmap="magma")
plt.colorbar()


stars_count[stars_count <= 13] = 0
stars_count[stars_count > 13] = 1
stars_full_small = stars_full[
    int(ny / 2 - cld_pix) : int(ny / 2 + cld_pix),
    int(nx / 2 - cld_pix) : int(nx / 2 + cld_pix),
]
sec2_small = sec2[
    int(ny / 2 - cld_pix) : int(ny / 2 + cld_pix),
    int(nx / 2 - cld_pix) : int(nx / 2 + cld_pix),
]

plt.figure()
plt.title("stars_count")
plt.imshow(stars_count, cmap="magma")
plt.colorbar()
plt.figure()
plt.title("stars_full_small")
plt.imshow(stars_full_small, cmap="magma")
plt.colorbar()

# weighted with solid angle
cloud_cover = round(
    (1 - np.sum(stars_count * sec2_small) / np.sum(stars_full_small * sec2_small)) * 100
)
cloud_cover_okta = round(cloud_cover / 12.5)
print("Cloud cover (%, oktas) : ", cloud_cover, cloud_cover_okta)

if Calmet == "stars":
    if cloud_cover_okta > max_cloud_cover:
        print("Can't process the data for stars calibration methods under cloudy skies")
    else:
        stars_selected = ds[(ds["MagV"] < limiti) & (ds["MagV"] > limits)]
        coordsradec = stars_selected["coord1_ICRS,J2000/2000_"]
        coords = coordsradec.to_numpy(dtype="str")
        identity = stars_selected["identifier"]
        iden = identity.to_numpy(dtype="str")
        magnitudev = stars_selected["MagV"]
        magv = magnitudev.to_numpy(dtype="float")
        magnituder = stars_selected["MagR"]
        magr = magnituder.to_numpy(dtype="float")
        magnitudeb = stars_selected["MagB"]
        magb = magnitudeb.to_numpy(dtype="float")
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
            # keep stars above horizon limitz = 0.95 * ny/2 * Zslope
            azistar = np.delete(azistar, np.where(altstar < elmin))
            magv = np.delete(magv, np.where(altstar < elmin))
            magr = np.delete(magr, np.where(altstar < elmin))
            magb = np.delete(magb, np.where(altstar < elmin))
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
        # calculating airmass with A. T. Young, "AIR-MASS AND REFRACTION," Applied Optics, vol. 33,
        #    pp. 1108-1110, Feb 1994.
        AirM = airmass(altstar)
        Itime_cor = 1
        # Search for stars only if cloud_cover is lower than max_cloud_cover
        mean, median, std = sigma_clipped_stats(imstars, sigma=5.0)
        daofind = DAOStarFinder(fwhm=2, threshold=5.0 * std)
        sources = daofind(imstars)
        print(sources)
        # keep stars with elevation > elmin degrees
        for col in sources.colnames:
            sources[col].info.format = "%.8g"  # for consistent table output
        positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))

        positions = (np.rint(positions)).astype(int)
        xsa = positions[:, 0]
        ysa = positions[:, 1]
        Flux = np.zeros(np.shape(positions)[0])
        Back = np.zeros(np.shape(positions)[0])
        for nd in range(np.shape(positions)[0]):
            xsa = positions[nd, 0]
            ysa = positions[nd, 1]
            Flux[nd] = np.sum(imstars[ysa - 2 : ysa + 3, xsa - 2 : xsa + 3])
            Back[nd] = np.sum(imbkg[ysa - 2 : ysa + 3, xsa - 2 : xsa + 3])
        rn = np.hypot(positions[:, 0] - nx / 2, positions[:, 1] - ny / 2)
        sources = sources[rn < rnmax]
        positions = positions[rn < rnmax]

        Flux = Flux[rn < rnmax]
        Back = Back[rn < rnmax]
        maxflux = Flux.max()
        brightest = -1
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
        apertures = CircularAperture(positions, r=4.0)
        StarMatch = np.zeros([ishape, 10])
        StarName = np.empty(ishape, dtype="object")
        n = 0
        nistars = ishape
        # searching for correspondance between stars in simbad and found stars in image
        dstar = np.zeros(np.shape(positions)[0])
        for ns in range(ishape):
            dstar = np.hypot(
                positions[:, 0] - index[ns, 1], positions[:, 1] - index[ns, 0]
            )
            dweight = dstar**2 / Flux
            dweight_min = np.amin(dweight)
            dweight_min_index = dweight.argmin()
            dmin = np.amin(dstar)
            # try to match if distance is lower that 1 deg
            if dmin < 10:
                StarMatch[n, 0] = index[ns, 1]
                StarMatch[n, 1] = index[ns, 0]
                StarMatch[n, 2] = positions[dweight_min_index, 0]  # noeuds[0]
                StarMatch[n, 3] = positions[dweight_min_index, 1]  # noeuds[1]
                StarMatch[n, 4] = dmin
                StarMatch[n, 5] = (
                    positions[dweight_min_index, 0] - index[ns, 1]
                )  # noeuds[0] - index[ns,1]
                StarMatch[n, 6] = (
                    positions[dweight_min_index, 1] - index[ns, 0]
                )  # noeuds[1] - index[ns,0]
                if Band == "JV":
                    StarMatch[n, 7] = magv[ns]
                elif Band == "JR":
                    StarMatch[n, 7] = magr[ns]
                elif Band == "JB":
                    StarMatch[n, 7] = magb[ns]
                StarMatch[n, 8] = AirM[ns]
                StarMatch[n, 9] = Flux[dweight_min_index]
                StarName[n] = iden[ns]
                n = n + 1
        print("Number of matching stars : ", n, "/", ishape)
        StarMatch[np.isnan(StarMatch)] = 0
        StarMatch = np.delete(StarMatch, np.where(StarMatch == 0), axis=0)
        StarName = np.delete(StarName, np.where(StarMatch == 0), axis=0)
        avggap, mediangap, stdgap = sigma_clipped_stats(StarMatch[:, 4], sigma=2.0)
        print("Average gap between nearest star :", avggap, "+/-", stdgap)
        StarMatch = np.delete(
           StarMatch, np.where(StarMatch[:, 4] > avggap + 3 * stdgap), axis=0
        )
        StarName = np.delete(
            StarName, np.where(StarMatch[:, 4] > avggap + 3 * stdgap), axis=0
        )
        StarMatch = np.delete(StarMatch, np.where(StarMatch[:, 9] == 0), axis=0)
        StarName = np.delete(StarName, np.where(StarMatch[:, 9] == 0), axis=0)
        print("Number of matching stars : ", np.shape(StarMatch)[0], "/", ishape)
        # save stars match information for future Calibration
        c = open(calname, "w")
        first_line = "Ident,Band,Airmass,Ext_coef,Catalog_magnitude,Instrumental_flux,Rcoef,Gcoef,Bcoef \n"
        c.write(first_line)
        for nc in range(np.shape(StarMatch)[0]):
            # Ident , Band , Airmass , Ext_coef , Catalog_magnitude , Instrumental_flux"
            # to find the extinction after the night one can make a linear regression of the uncalibrated magnitude of a star
            # as a function of the airmass. Equation is m= m_0 + k * airmass
            # so that k is the slope for a given stars mesures at various airmass.
            cal_line = (
                StarName[nc]
                + ","
                + Band
                + ","
                + str("{:4.2f}".format(StarMatch[nc, 8]))
                + ","
                + str("{:5.3f}".format(k))
                + ","
                + str("{:5.3f}".format(StarMatch[nc, 7]))
                + ","
                + str("{:5.3f}".format(StarMatch[nc, 9]))
                + ","
                + str(RC)
                + ","
                + str(GC)
                + ","
                + str(BC)
                + " \n"
            )
            c.write(cal_line)
        c.close()

        if Band == "JV":
            lambm = 545
        elif Band == "JR":
            lambm = 700
        elif Band == "JB":
            lambm = 436

        # calculate uncalibrated mag, extinction, calibration factor
        uncalMag = -2.5 * np.log10(StarMatch[:, 9])
        # correct the extinction for the reference stars
        calMag = StarMatch[:, 7] + k * StarMatch[:, 8]

        ax = StarMatch[:, 9]
        ay = 10 ** (-0.4 * calMag)
        cr = np.corrcoef(ax, ay)
        corcoef = cr[0, 1]
        print("Correlation coefficient : ", corcoef)
        params = curve_fit(fit_func, ax, ay)
        slp = float(params[0])
        print("n ax=", np.shape(ax)[0])
        # filtering outliers (may be stars behind semi-transparent cloud or bad matching)
        deltacor = 1000
        ndelta = 0
        while deltacor > 0.001:
            residuals = ay - slp * ax
            res_mean, res_median, res_std = sigma_clipped_stats(residuals, sigma=3.0)
            residuals[residuals > res_mean + 2 * res_std] = -1000
            residuals[residuals < res_mean - 2 * res_std] = -1000
            axp = np.delete(ax, np.where(residuals == -1000))
            ayp = np.delete(ay, np.where(residuals == -1000))
            cr = np.corrcoef(axp, ayp)
            deltacor = corcoef
            corcoef = cr[0, 1]
            deltacor = abs(deltacor - corcoef)
            print("Correlation coefficient : ", corcoef)
            params = curve_fit(fit_func, axp, ayp)
            slp = float(params[0])
            gx = np.linspace(0, np.amax(ax), 100)
            gy = slp * gx
            if ndelta == 5:
                break
            ndelta += 1
            npts = np.shape(axp)[0]
            
            plt.figure()
            plt.plot(gx, gy, "r")
            plt.plot(ax, ay, "ob")
            plt.plot(axp, ayp, "or")
            plt.xlabel("Star's pixel values")
            plt.ylabel("10^(-0.4*CalMag)")
            #plt.show()
        if corcoef > 0.7 and npts > 9:
            calsb = 1

        file = path + Cam + "_" + Band + "_Stars_Match_" + baseout + ".png"
        title = Band + " Stars correspondance"
        plt.figure()
        plt.plot(StarMatch[:, 2], StarMatch[:, 3], "or", markersize=2)
        plt.plot(StarMatch[:, 0], StarMatch[:, 1], "ok", markersize=2)
        plt.ylim(ny, 0)
        plt.xlim(0, nx)
        plt.title(title)
        plt.legend(["Detected stars", "Simbad reference stars"])
        plt.savefig(file)

        title = Band + " calibration"
        file = path + Cam + "_"  + Band + "_calibration_" + baseout + ".png"
        plt.figure()
        plt.plot(gx, gy, "r")
        plt.plot(ax, ay, "ob")
        plt.plot(axp, ayp, "or")
        plt.xlabel("Star's pixel values")
        plt.ylabel("10^(-0.4*CalMag)")
        plt.title(title)
        plt.savefig(file)
elif Calmet == "fixed":
    calsb = 1
    Itime_cor = DefaultItime / Tint * DefaultGain / Gain
if calsb == 1:
    # print calibration slope slp
    print("Calibration slope (slp) :", slp)
    # replace zeros with a small non null value
    mask = imag <= 0
    imagtmp = np.copy(imag)
    imagtmp[mask] = imbkg[mask]
    imag = imagtmp
    imstars[imstars <= 0] = 0.0001
    clouds = imstars[imstars > 0.001]
    calMagTot = -2.5 * np.log10(imag * Itime_cor * slp)
    calMagBkg = -2.5 * np.log10(imbkg * Itime_cor * slp)
    calMagStr = -2.5 * np.log10(imstars * Itime_cor * slp)
    zeropoint = float(-2.5 * np.log10(slp))
    print("Zero point (mag) =", zeropoint)
    print("Mag(pixel) = Zeropoint +-2.5 * log10(R(pixel))")
    print(" where R the signal assigned to the pixel")

    # calibrated surface brightness in mag /sq arc sec
    calSbTot = calMagTot + 2.5 * np.log10(sec2)
    calSbBkg = calMagBkg + 2.5 * np.log10(sec2)
    calSbStr = calMagStr + 2.5 * np.log10(sec2)
    calSbTot[z > 90] = np.nan
    calSbBkg[z > 90] = np.nan
    calSbStr[z > 90] = np.nan

    norm1 = simple_norm(calSbBkg, "sqrt")
    title = Band + " background Surface Brightness"
    file = path + Cam + "_"  + Band + "_calSbBkg_" + baseout + ".png"
    plt.figure()
    plt.imshow(-calSbBkg, cmap="magma", vmin=-22, vmax=-12)
    plt.colorbar()
    plt.title(title)
    plt.savefig(file)

    norm1 = simple_norm(calSbTot, "sqrt")
    title = Band + " total Surface Brightness"
    file = path + Cam + "_"  + Band + "_calSbTot_" + baseout + ".png"
    plt.figure()
    plt.imshow(-calSbTot, cmap="magma", vmin=-22, vmax=-12)
    plt.colorbar()
    plt.title(title)
    plt.savefig(file)


# plt.figure()
# plt.title("stars bin")
# plt.imshow(stars_binary, cmap="inferno")
# plt.colorbar()
#
# plt.figure()
# plt.title("starcount")
# plt.imshow(stars_count, cmap="inferno")
# plt.colorbar()
# plt.figure()
# plt.title("full")
# plt.imshow(stars_full, cmap="inferno")
# plt.colorbar()

# Extract points and write data
index = find_close_indices(az, el, apt, ept)
pshape = int(np.shape(index)[0])
pshape = int(pshape / 2)
index = index.reshape(pshape, 2)
index = np.delete(index, (0), axis=0)
pshape = int(np.shape(index)[0])

for no in range(num_pts):
    posx = str(index[no, 1])
    posy = str(index[no, 0])
    direction = here_now.from_altaz(alt_degrees=ept[no], az_degrees=apt[no])
    glat, glon, gdistance = direction.frame_latlon(galactic_frame)
    galactic_lat = glat.degrees
    theta_sun = direction.separation_from(sun_position).degrees
    theta_moon = direction.separation_from(moon_position).degrees
    # calculate Airmass
    airmo = airmass(ept[no])
    # read V sky Brightness
    if calsb == 1:
        mago = calSbBkg[index[no, 0], index[no, 1]]
    else:
        mago = np.nan
    o = open(outname, "a")
    outputline = (
        Site
        + ";"
        + Band
        + ";"
        + posx
        + ";"
        + posy
        + ";"
        + str("{:6.2f}".format(apt[no]))
        + ";"
        + str("{:5.2f}".format(ept[no]))
        + ";"
        + str("{:4.2f}".format(airmo))
        + ";"
        + str("{:5.3f}".format(k))
        + ";"
        + timestamp
        + ","
        + mflag
        + ";"
        + str("{:d}".format(cloud_cover_okta))
        + ";"
        + str("{:6.3f}".format(mago))
        + ";"
        + str("{:6.3f}".format(sberr))
        + ";"
        + str("{:6.3f}".format(zeropoint))
        + ";"
        + str("{:6.3f}".format(corcoef))
        + ";"
        + str("{:6.2f}".format(theta_sun))
        + ";"
        + str("{:6.2f}".format(theta_moon))
        + ";"
        + str("{:6.2f}".format(galactic_lat))
        + ";"
        + str("{:6.2f}".format(moonphase))
        + "\n"
    )
    print(outputline)
    o.write(outputline)
    o.close()
# plt.show()
