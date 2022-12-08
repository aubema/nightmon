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
    Sfile = "undefined"
    # reading V and R images names from the input parameters
    # default extinction values will be replaces if provided to the command line arguments
    # RGO, D. K. (1985). Atmospheric Extinction at the Roque de los Muchachos Observatory, La Palma.
    # extinc_r = 0.0547
    # extinc_v = 0.102
    Extinc = 0.102
    try:
        opts, args = getopt.getopt(
            argv,
            "h:s:d:b:e:c:",
            ["help=", "sfile=", "dfile=", "extinc=", "band=", "cam="],
        )
    except getopt.GetoptError:
        print("ProcessNighMon.py -s <Sfile> -d <Dfile> -b <Band> -e <Extinc> -c <Cam>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(
                "ProcessNighMon.py -s <Sfile> -d <Dfile> -b <Band> -e <Extinc> -c <Cam>"
            )
            sys.exit()
        elif opt in ("-s", "--sfile"):
            Sfile = arg
        elif opt in ("-d", "--dfile"):
            Dfile = arg
        elif opt in ("-b", "--band"):
            Band = arg
        elif opt in ("-e", "--extinc"):
            Extinc = arg
        elif opt in ("-c", "--cam"):
            Cam = arg
    print("Sky image file is :", Sfile)
    print("Dark frame file is :", Dfile)
    print("Band is :", Band)
    print("Extinction is :", Extinc)
    print("Camera is :", Cam)
    return Sfile, Dfile, Band, Extinc, Cam


def fit_func(x, a):
    # Curve fitting function
    return a * x


# ================================================
# MAIN
# default Parameters
sberr = 0
mflag = "False"
FWHM = 7
norm = ImageNormalize(stretch=SqrtStretch())

elmin = 5  # set the minimum elevation
# selecting magnitude limits for reference stars in the simbad database
# limits allow to remove saturated stars in the images
limiti = (
    3.7  # limitting stars magnitude NEED TO BE LOWER THAN 6 but 3.7 is a magic number
)
limits = 1.2
# load command line parameters
Sfile, Dfile, Band, Extinc, Cam = input(sys.argv[1:])
# determine the R, G, B coefficients according to the band
if Band == "JV":
    RC = 0.1
    GC = 1
    BC = 0.1
elif Band == "JR":
    RC = 1
    GC = 0.1
    BC = 0.1
elif Band == "R":
    RC = 1
    GC = 0
    BC = 0
elif Band == "G":
    RC = 0
    GC = 1
    BC = 0
elif Band == "B":
    RC = 0
    GC = 0
    BC = 1
# open images
print("Reading images...")
Simg = open_raw(Sfile)
Dimg = open_raw(Dfile)
# read site coordinates
# Load Parameters
# home = os.path.expanduser("~")
# read NightMon config file
deltax = np.empty([3], dtype=float)
deltay = np.empty([3], dtype=float)
angle = np.empty([3], dtype=float)
# TODO: restore the following line
# with open("/home/sand/nightmon_config") as f:
with open("/home/aubema/nightmon_config") as f:

    p = yaml.safe_load(f)
Site = p["Site"]
if Cam == "A":
    deltax = p["ShiftxA"].split()
    deltay = p["ShiftyA"].split()
    angle = p["AngleA"].split()
elif Cam == "B":
    deltax = p["ShiftxB"].split()
    deltay = p["ShiftyB"].split()
    angle = p["AngleB"].split()

# load sky brightness extraction points
# TODO: restore the following line
# pt = open("/home/sand/points_list", "r")
pt = open("/home/aubema/points_list", "r")


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
print("Number of points to extract :", num_pts)
ts = load.timescale()
# set observer position
eph = load("de421.bsp")
here = eph["earth"] + wgs84.latlon(
    p["Latitude"] * N, p["Longitude"] * E, elevation_m=p["Altitude"]
)
# create the grayscale images
Sgray = to_grayscale(Simg, [RC, GC, BC], normalize=False)
Dgray = to_grayscale(Dimg, [RC, GC, BC], normalize=False)
ny = np.shape(Sgray)[0]
nx = np.shape(Sgray)[1]
# set the minimum elevation to the limit of the image if required
if elmin < 90 - 0.98 * (ny / 2 * p["Zslope"] + (ny / 2) ** 2 * p["Zquad"]):
    elmin = 90 - 0.98 * (ny / 2 * p["Zslope"] + (ny / 2) ** 2 * p["Zquad"])
print("Minimum elevation :", elmin)
# calculate maximum zenith angle
zemax = 90 - elmin
rnmax = (-p["Zslope"] + np.sqrt(p["Zslope"] ** 2 - 4 * p["Zquad"] * -zemax)) / (
    2 * p["Zquad"]
)

# remove darks
window = 7  #  1 deg clouds ~= 10
kernel = Box2DKernel(width=window, mode="integrate")
dark = convolve(Dgray, kernel)
dark_norm = convolve(Sgray, kernel)
dmean = np.mean(dark[5 : ny - 5, 5:105], axis=(0, 1))
imean = np.mean(Sgray[5 : ny - 5, 5:105], axis=(0, 1))
dnorm = imean / dmean
dark = dark * dnorm
# plt.figure()
# plt.imshow(Vgray, cmap="inferno", norm=norm)
# plt.colorbar()
# plt.title("img")
Sgray = Sgray - dark
# plt.figure()
# plt.imshow(Vgray, cmap="inferno", norm=norm)
# plt.colorbar()
# plt.title("img-dark")
# plt.show()

# **************************** supprimer
rm = 0.635
vv = np.array(
    [
        0.99903,
        0.9990698566039953,
        0.9981855774518988,
        0.9965508739460722,
        0.9943394574888771,
        0.9917250394826749,
        0.9888760750936606,
        0.9859422041998309,
        0.9829415915486667,
        0.9798735483421788,
        0.9767373857823781,
        0.9735324150712757,
        0.9702579474108823,
        0.966913294003209,
        0.9634977660502667,
        0.9600106747540663,
        0.9564513313166186,
        0.9528190469399347,
        0.9491131328260255,
        0.9453329001769019,
        0.9414776601945748,
        0.9375467240810552,
        0.9335394030383539,
        0.9294550082684819,
        0.9252928509734502,
        0.9210522423552696,
        0.9167324936159511,
        0.9123329159575057,
        0.9078528205819441,
        0.9032915186912773,
        0.8986483214875165,
        0.8939225401726723,
        0.8891134859487557,
        0.8842204700177777,
        0.8792428035817492,
        0.8741797978426812,
        0.8690307640025845,
        0.8637950132634701,
        0.8584718568273488,
        0.8530606058962317,
        0.8475605716721296,
        0.8419710653570536,
        0.8362913981530145,
        0.8305208812620231,
        0.8246588258860905,
        0.8187045432272276,
        0.8126573444874453,
        0.8065165408687546,
        0.8002814435731663,
        0.7939513638026915,
        0.7875250748828453,
        0.7810021458037676,
        0.7743799985558575,
        0.7676564974881454,
        0.7608295069496618,
        0.7538968912894375,
        0.7468565148565028,
        0.7397062419998887,
        0.7324439370686252,
        0.7250674644117431,
        0.717574688378273,
        0.7099634733172454,
        0.7022316835776909,
        0.69437718350864,
        0.6863978374591233,
        0.6782915097781714,
        0.6700560648148147,
        0.6616893669180839,
        0.6531892804370095,
        0.6445536697206221,
        0.6357803991179521,
        0.6268673329780303,
        0.6178123356498872,
        0.6086132714825532,
        0.599268004825059,
        0.589774400026435,
        0.580130321435712,
        0.5703336334019204,
        0.5603822002740906,
        0.5502738864012535,
        0.5400065561324394,
        0.5295780738166791,
        0.5189863038030028,
        0.5082291104404415,
        0.4973043580780254,
        0.486209911064785,
        0.4749436337497512,
        0.4635033904819544,
        0.4518870456104251,
        0.4400924634841938,
        0.4281175084522914,
        0.4159600448637481,
        0.4028579591928673,
        0.3841800074975651,
        0.3540685156254719,
        0.306664774044784,
        0.2361100732236971,
        0.1365457036304061,
        0.02552628311479194,
        0.02002000000000018,
    ]
)
dy = 0
dx = 0
yy = (np.linspace(1, ny, ny) - ny / 2.0 - dy) / ny
xx = (np.linspace(1, nx, nx) - nx / 2.0 - dx) / ny
[xxx, yyy] = np.meshgrid(xx, yy)
rrr = np.sqrt(xxx * xxx + yyy * yyy) / rm
rrr[rrr > 1.0] = 0.0
r = np.linspace(0, 1, 100)
flat = np.interp(rrr, r, vv)
# ****************************
# correct flat field
Sgray = Sgray / flat
# creating elevation, azimith maps
print("Creating Azimuth and elevation maps...")
y, x = np.indices((ny, nx))
# computing the distance to zenith in pixels
d = np.hypot(x - nx / 2, y - ny / 2)
d[d < 0.5] = 0.5
z = p["Zslope"] * d + p["Zquad"] * d**2
z[z < 0] = 0
# computing azimuth
az = np.arctan2(-x + nx / 2, -y + ny / 2) * 180 / np.pi
# az = az - azpol
az = np.where(az < 0, az + 360, az)
# solid angle in sq arc second
sec2 = ((p["Zslope"] + p["Zquad"]) * 180 * 3600**2 * np.sin((z) * np.pi / 180)) / (
    d * np.pi
)
# compute elevation angle
el = 90 - z
ImgAirm = airmass(el)
ImgAirm[z >= 90] = np.nan


# initialize np arrays for backgrounds and stars fields
Vbkg = np.zeros([ny, nx])
Rbkg = np.zeros([ny, nx])
Sstars = np.zeros([ny, nx])


print(f"Processing Johnson {Band} camera...")
# TOTO : IL FAUT DEPLACER CECI ET LES EPHEMERIDES DANS LA BOUCLE DES FILTRES
# files names provided by NightMon are formatted in the following way
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
basename = time.utc_strftime("%Y-%m-%d")
outname = "calibrated_" + Band + "_" + basename + "_sky.csv"
timestamp = time.utc_strftime("%Y-%m-%dT%H:%M:%S")
baseout = time.utc_strftime("%Y-%m-%d_%H-%M-%S")
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
# creating star map with the SIMBAD stars database
# TODO: restore the two following lines
# ds = pd.read_csv(
#    "/home/sand/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv", header=0, sep=";"
ds = pd.read_csv(
    "/home/aubema/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv", header=0, sep=";"
)
# stars_selected=ds[ds['MagR'] < limit]
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
# calculating airmass with A. T. Young, "AIR-MASS AND REFRACTION," Applied Optics, vol. 33,
#    pp. 1108-1110, Feb 1994.
AirM = airmass(altstar)
# create the data file if it do not exists
if os.path.exists(outname) == False:
    o = open(outname, "w")
    first_line = "# Loc_Name , Band , CCD_XY_position , ,  AzAlt_position , , Airmass , Ext_coef ,    \
    Date , Moon , Clouds ,   SkyBrightness , err , Zeropoint , Sun_Alt, Moon_Alt , Galactic_Lat , Moon_Phase \n"
    second_line = "# ,  , (pixel) , (pixel) , (deg) , (deg) ,  ,  ,  \
    , (%) , (mag/arcsec^2) ,  , mag , deg , deg ,  deg , deg , \n"
    o.write(first_line)
    o.write(second_line)
    o.close()
imag = Sgray
shiftx = 0
shifty = 0
for i in range(3):
    shiftx = float(deltax[i])
    shifty = float(deltay[i])
    theta = float(angle[i])
    imag = ndimage.shift(imag, [shifty, shiftx], mode="nearest")
    padX = [imag.shape[1] - polindex[0, 1], polindex[0, 1]]
    padY = [imag.shape[0] - polindex[0, 0], polindex[0, 0]]
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
# imstars[imstars < 0] = 0

mean, median, std = sigma_clipped_stats(imstars, sigma=3.0)
daofind = DAOStarFinder(fwhm=3, threshold=5.0 * std)
sources = daofind(imstars)

# keep stars with elevation > elmin degrees
for col in sources.colnames:
    sources[col].info.format = "%.8g"  # for consistent table output
positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))
Flux = sources["flux"]
Peak = sources["peak"]
rn = np.hypot(positions[:, 0] - nx / 2, positions[:, 1] - ny / 2)
sources = sources[rn < rnmax]
positions = positions[rn < rnmax]
Flux = Flux[rn < rnmax]
Peak = Peak[rn < rnmax]
maxflux = Flux.max()
sources = sources[Flux > maxflux / 2.5 ** (limiti - brightest)]
positions = positions[Flux > maxflux / 2.5 ** (limiti - brightest)]
Peak = Peak[Flux > maxflux / 2.5 ** (limiti - brightest)]
Flux = Flux[Flux > maxflux / 2.5 ** (limiti - brightest)]
# Flux = Flux[ rn < rnmax]
# find most probable Polaris
dtopol = np.hypot(polindex[0, 1] - positions[:, 0], polindex[0, 0] - positions[:, 1])
mintopol = np.nanmin(dtopol)
indtopol = np.where(dtopol == mintopol)
xpoli = float(positions[indtopol, 0])
ypoli = float(positions[indtopol, 1])
# positions[:, 1] = ny - positions[:, 1]
apertures = CircularAperture(positions, r=4.0)

StarMatch = np.zeros([ishape, 10])
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
        xs = round(StarMatch[n, 2])
        ys = round(StarMatch[n, 3])
        StarMatch[n, 4] = dmin
        StarMatch[n, 5] = (
            positions[dmin_index, 0] - index[ns, 1]
        )  # noeuds[0] - index[ns,1]
        StarMatch[n, 6] = (
            positions[dmin_index, 1] - index[ns, 0]
        )  # noeuds[1] - index[ns,0]
        if Band == "JV":
            StarMatch[n, 7] = magv[ns]
        elif Band == "JR":
            StarMatch[n, 7] = magr[ns]
        StarMatch[n, 8] = AirM[ns]
        StarMatch[n, 9] = np.sum(imstars[ys - 2 : ys + 2, xs - 2 : xs + 2])
        n = n + 1
StarMatch[np.isnan(StarMatch)] = 0
StarMatch = np.delete(StarMatch, np.where(StarMatch == 0), axis=0)
avggap, mediangap, stdgap = sigma_clipped_stats(StarMatch[:, 4], sigma=2.0)
print("Average distance between nearest star :", avggap, "+/-", stdgap)
StarMatch = np.delete(
    StarMatch, np.where(StarMatch[:, 4] > avggap + 1 * stdgap), axis=0
)
StarMatch = np.delete(StarMatch, np.where(StarMatch[:, 9] == 0), axis=0)
k = Extinc
if Band == "JV":
    lambm = 551
    cfactor = 1200
elif Band == "JR":
    lambm = 658
    cfactor = 400
print("Zenith atmospheric extinction (mag) :", k)

# calculate uncalibrated mag, extinction, calibration factor
uncalMag = -2.5 * np.log10(StarMatch[:, 9])
# correct the extinction for the reference stars
calMag = StarMatch[:, 7] + k * StarMatch[:, 8]

ax = StarMatch[:, 9]
ay = 10 ** (-0.4 * calMag)
params = curve_fit(fit_func, ax, ay)
slp = params[0]
gx = np.linspace(0, np.amax(ax), 100)
gy = slp * gx

# replace zeros with a small non null value
mask = imag <= 0
imagtmp = np.copy(imag)
imagtmp[mask] = imbkg[mask]
imag = imagtmp
calMagTot = -2.5 * np.log10(imag * slp)

calMagBkg = -2.5 * np.log10(imbkg * slp)
zeropoint = float(-2.5 * np.log10(slp))
print("Zero point (mag) =", zeropoint)
print("Mag(pixel) = Zeropoint +-2.5 * log10(R(pixel))")
print(" where R the signal assigned to the pixel")

# calibrated surface brightness in mag /sq arc sec
calSbTot = calMagTot + 2.5 * np.log10(sec2)
calSbBkg = calMagBkg + 2.5 * np.log10(sec2)
calSbTot[z > 90] = np.nan
calSbBkg[z > 90] = np.nan

title = Band + " calibration"
file = Band + "_calibration_" + baseout + ".png"
plt.figure()
plt.plot(gx, gy, "r")
plt.plot(ax, ay, "ob")
plt.xlabel("Star's pixel values")
plt.ylabel("10^(-0.4*CalMag)")
plt.title(title)
plt.savefig(file)

norm1 = simple_norm(calSbBkg, "sqrt")
title = Band + " background Surface Brightness"
file = Band + "_calSbBkg_" + baseout + ".png"
plt.figure()
plt.imshow(-calSbBkg, cmap="magma", vmin=-22, vmax=-16)
plt.colorbar()
plt.title(title)
plt.savefig(file)

norm1 = simple_norm(calSbTot, "sqrt")
title = Band + " total Surface Brightness"
file = Band + "_calSbTot_" + baseout + ".png"
plt.figure()
plt.imshow(-calSbTot, cmap="magma", vmin=-22, vmax=-16)
plt.colorbar()
plt.title(title)
plt.savefig(file)

file = Band + "_tars_Match_" + baseout + ".png"
title = Band + " Stars correspondance"
plt.figure()
plt.plot(StarMatch[:, 2], StarMatch[:, 3], "or", markersize=2)
plt.plot(StarMatch[:, 0], StarMatch[:, 1], "ok", markersize=2)
plt.ylim(ny, 0)
plt.xlim(0, nx)
plt.title(title)
plt.legend(["Detected stars", "Simbad reference stars"])
plt.savefig(file)

# determine cloud cover
threshold = np.amax(imstars) / cfactor
stars_binary = np.full((ny, nx), 0)
stars_full = np.full((ny, nx), 0)
stars_binary[imstars >= threshold] = 1.0
stars_full[imstars >= threshold / 100] = 1.0
stars_binary[z > 90] = 0
stars_full[z > 90] = 0
window = 19  #  1 deg clouds ~= 10
kernel = Box2DKernel(width=window, mode="integrate")
stars_count = convolve(stars_binary, kernel)
stars_count[z > 90] = 0.0
stars_full_count = convolve(stars_full, kernel)
stars_full_count[z > 90] = 0.0
stars_count[stars_count > 0] = 1
stars_full_count[stars_full_count > 0] = 1
# weighted with solid angle
cloud_cover = int(
    (1 - np.sum(stars_count * sec2) / np.sum(stars_full_count * sec2)) * 100
)
print("Cloud cover (%) : ", cloud_cover)

#
# plt.figure()
# plt.title("stars")
# plt.imshow(imstars, cmap="inferno" )
# plt.colorbar()

# Extract points and write data
index = find_close_indices(az, el, apt, ept)
pshape = int(np.shape(index)[0])
pshape = int(pshape / 2)
index = index.reshape(pshape, 2)
index = np.delete(index, (0), axis=0)
pshape = int(np.shape(index)[0])

for no in range(num_pts - 1):
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
    mago = calSbBkg[index[no, 0], index[no, 1]]
    o = open(outname, "a")
    outputline = (
        Site
        + " , "
        + Band
        + " , "
        + posx
        + " , "
        + posy
        + " , "
        + str("{:6.2f}".format(apt[no]))
        + " , "
        + str("{:5.2f}".format(ept[no]))
        + " , "
        + str("{:4.2f}".format(airmo))
        + " , "
        + str("{:5.3f}".format(k))
        + " , "
        + timestamp
        + " , "
        + mflag
        + " , "
        + str("{:5.1f}".format(cloud_cover))
        + " , "
        + str("{:6.3f}".format(mago))
        + " , "
        + str("{:6.3f}".format(sberr))
        + " , "
        + str("{:6.3f}".format(zeropoint))
        + " , "
        + str("{:6.2f}".format(theta_sun))
        + " , "
        + str("{:6.2f}".format(theta_moon))
        + " , "
        + str("{:6.2f}".format(galactic_lat))
        + " , "
        + str("{:6.2f}".format(moonphase))
        + "\n"
    )
    print(outputline)
    o.write(outputline)
    o.close()

plt.show()
