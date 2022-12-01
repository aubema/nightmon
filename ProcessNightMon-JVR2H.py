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


def closest_node(node, nodes):
    dist = np.sqrt(
        (nodes[cdist([node], nodes).argmin()][0] - node[0]) ** 2
        + (nodes[cdist([node], nodes).argmin()][1] - node[1]) ** 2
    )
    return nodes[cdist([node], nodes).argmin()], dist, cdist([node], nodes).argmin()


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
    Vfile = "undefined"
    Rfile = "undefined"
    Dfile = "undefined"
    # reading V and R images names from the input parameters
    # default extinction values will be replaces if provided to the command line arguments
    # RGO, D. K. (1985). Atmospheric Extinction at the Roque de los Muchachos Observatory, La Palma.
    extinc_r = 0.0547
    extinc_v = 0.102
    try:
        opts, args = getopt.getopt(
            argv,
            "h:v:r:d:ev:er:",
            ["help=", "vfile=", "dfile=", "rfile=", "extv=", "extr="],
        )
    except getopt.GetoptError:
        print("test.py -v <Vfile> -r <Rfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("test.py -v <Vfile> -r <Rfile>")
            sys.exit()
        elif opt in ("-v", "--vfile"):
            Vfile = arg
        elif opt in ("-r", "--rfile"):
            Rfile = arg
        elif opt in ("-d", "--dfile"):
            Dfile = arg
        elif opt in ("-ev", "--extv"):
            extinc_v = arg
        elif opt in ("-er", "--extr"):
            extinc_r = arg
    print("Johnson V file is :", Vfile)
    print("Johnson R file is :", Rfile)
    return Vfile, Rfile, Dfile, extinc_v, extinc_r


# ================================================
# MAIN
# default Parameters
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
Vfile, Rfile, Dfile, extinc_v, extinc_r = input(sys.argv[1:])
# open images
Vimg = open_raw(Vfile)
Rimg = open_raw(Rfile)
Dimg = open_raw(Dfile)
# read site coordinates
# Load Parameters
home = os.path.expanduser("~")
# read NightMon config file
with open(home + "/nightmon_config") as f:
    p = yaml.safe_load(f)
Site = p["Site"]
# reading coefficients to convert RGB into gray scale
RC = p["R2GrayCoef"]
GC = p["G2GrayCoef"]
BC = p["B2GrayCoef"]
# load sky brightness extraction points
pt = open(home + "/points_list", "r")
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
Vgray = to_grayscale(Vimg, [RC, GC, BC], normalize=False)
Rgray = to_grayscale(Rimg, [RC, GC, BC], normalize=False)
Dgray = to_grayscale(Dimg, [RC, GC, BC], normalize=False)
ny = np.shape(Vgray)[0]
nx = np.shape(Vgray)[1]
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
dark_norm = convolve(Vgray, kernel)
dmean = np.mean(dark[5 : ny - 5, 5:105], axis=(0, 1))
imean = np.mean(Vgray[5 : ny - 5, 5:105], axis=(0, 1))
dnorm = imean / dmean
dark = dark * dnorm
# plt.figure()
# plt.imshow(Vgray, cmap="inferno", norm=norm)
# plt.colorbar()
# plt.title("img")
Vgray = Vgray - dark
Rgray = Rgray - dark
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
Vgray = Vgray / flat
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
Vstars = np.zeros([ny, nx])
Rstars = np.zeros([ny, nx])

# ================================
# loop over the two bands
for band, xpoli, ypoli, imagi, imbkg, imstars in (
    (
        "V",
        p["XpolarisV"],
        p["YpolarisV"],
        Vgray,
        Vbkg,
        Vstars,
    ),
    (
        "R",
        p["XpolarisR"],
        p["YpolarisR"],
        Rgray,
        Rbkg,
        Rstars,
    ),
):
    print(f"Processing Johnson {band} camera...")
    # TOTO : IL FAUT DEPLACER CECI ET LES EPHEMERIDES DANS LA BOUCLE DES FILTRES
    # files names provided by NightMon are formatted in the following way
    # YYYYY-MM-DD_hh-mm-ss_V.dng , the V letter represents the filter and can thus be
    # replaced by R. Time and date are in UTC.
    # time=ts.utc(2020, 9, 22, 4, 22, 53)
    # time = ts.now()
    # yesterday = time - timedelta(days = 1)
    datearr = Vfile.replace("_", "-").split("-")
    time = ts.utc(
        int(datearr[0]),
        int(datearr[1]),
        int(datearr[2]),
        int(datearr[3]),
        int(datearr[4]),
        int(datearr[5]),
    )
    basename = time.utc_strftime("%Y-%m-%d")
    outname = "calibrated" + "_" + basename + "_sky.csv"
    timestamp = time.utc_strftime("%Y-%m-%dT%H:%M:%S")
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
    ds = pd.read_csv(
        home + "/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv", header=0, sep=";"
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
    imag = imagi
    shiftx = 0
    shifty = 0
    for i in range(2):
        shiftx = polindex[0, 1] - xpoli
        shifty = polindex[0, 0] - ypoli
        theta1 = np.arctan2(ypoli - ny / 2, xpoli - nx / 2) * 180 / np.pi
        theta2 = (
            np.arctan2(polindex[0, 0] - ny / 2, polindex[0, 1] - nx / 2) * 180 / np.pi
        )
        theta = theta1 - theta2
        angpol1 = np.arctan2(xpoli - nx / 2, ny / 2 - ypoli) * 180 / np.pi
        angpol = angpol1 + azi_pol
        imag = ndimage.shift(imag, [shifty, shiftx], mode="nearest")
        padX = [imag.shape[1] - polindex[0, 1], polindex[0, 1]]
        padY = [imag.shape[0] - polindex[0, 0], polindex[0, 0]]
        imgP = np.pad(imag, [padY, padX], "constant")
        imag = ndimage.rotate(imgP, theta, reshape=False, mode="nearest")
        imag = imag[padY[0] : -padY[1], padX[0] : -padX[1]]

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
        # Flux = Flux[ rn < rnmax]
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

    StarMatch = np.zeros([ishape, 10])
    n = 0
    nistars = ishape
    # searching for correspondance between stars in simbad and found stars in image
    dstar = np.zeros(np.shape(positions)[0])
    for ns in range(ishape):
        # for npos in range(np.shape(positions)[0]-1):
        #     dstar[npos] = np.hypot(positions[npos,0]-index[ns,1],positions[npos,1]-index[ns,0])
        dstar = (
            np.hypot(positions[:, 0] - index[ns, 1], positions[:, 1] - index[ns, 0])
            ** 2
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
            if band == "V":
                StarMatch[n, 7] = magv[ns]
            elif band == "R":
                StarMatch[n, 7] = magr[ns]
            StarMatch[n, 8] = AirM[ns]
            StarMatch[n, 9] = Flux[dmin_index]
            n = n + 1
    StarMatch[np.isnan(StarMatch)] = 0
    StarMatch = np.delete(StarMatch, np.where(StarMatch == 0), axis=0)
    avggap, mediangap, stdgap = sigma_clipped_stats(StarMatch[:, 4], sigma=2.0)
    print("Average distance between nearest star :", avggap, "+/-", stdgap)
    StarMatch = np.delete(
        StarMatch, np.where(StarMatch[:, 4] > avggap + 1 * stdgap), axis=0
    )
    StarMatch = np.delete(StarMatch, np.where(StarMatch[:, 9] == 0), axis=0)

    if band == "V":
        lambm = 551
        cfactor = 1200
        k = extinc_v
    elif band == "R":
        lambm = 658
        cfactor = 400
        k = extinc_r
    print("Zenith atmospheric extinction (mag) :", k)

    # calculate uncalibrated mag, extinction, calibration factor
    uncalMag = -2.5 * np.log10(StarMatch[:, 9])
    calMag = StarMatch[:, 7] + k * StarMatch[:, 8]
    deltam = uncalMag - calMag

    # fit de la magnitude hors atm en fonction de la mag non calibrée
    ci = np.polyfit(uncalMag, calMag, 1)
    slopei = ci[0]
    origini = ci[1]
    plt.figure()
    title = band + " Zero point"
    plt.plot(uncalMag, calMag, "ob")
    gx = np.linspace(np.amin(uncalMag), np.amax(uncalMag), 100)
    gy = ci[1] + gx * ci[0]
    plt.plot(gx, gy, "r")
    plt.title(title)
    plt.ylabel("Calibrated Magnitude")
    plt.xlabel("Uncalibrated Magnitude")
    file = band + "zeropoint_corr.png"
    plt.savefig(file)

    zeropoint = origini
    print("Zero point (mag) =", zeropoint)

    # filter outliers zscore = (x - mean)/std
    am = StarMatch[:, 8]
    co = np.polyfit(am, deltam, 1)
    slope = co[0]
    origin = co[1]
    fx = np.linspace(0, np.amax(am), 100)
    fy = co[1] + fx * co[0]
    title = band + " dDelta_Mag vs Air Mass"
    plt.figure()
    plt.plot(am, deltam, "or", markersize=2)
    plt.plot(fx, fy, "b")
    plt.title(title)
    plt.xlabel("Air Mass")
    plt.ylabel("Delta Magnitude")
    # atmospheric optical depth (aerosols + molecules)
    # tau = slope / (2.5 * np.log10(np.e))
    # extinction = slope
    #  molecular optical depth - Bodhaine et al. (1999)
    # tau_mol = (0.0021520*(1.0455996-341.29061*(lambm/1000.)**
    #     -2.-0.90230850*(lambm/1000.)**2.)/(1.+0.0027059889*(lambm/1000.)**
    #     -2.-85.968563*(lambm/1000.)**2.))
    # AOD = tau-tau_mol
    # print("Atmospheric optical depth =", tau)
    # print("Aerosol optical depth =",AOD," at ",lambm," nm")

    # calculate a calibrated magnitude
    # calMag = -2.5 * np.log10(pixel_value) - origin
    # calMagTot = -2.5 * np.log10(imag) - origini + k * AirM
    calMagBkg = -2.5 * np.log10(imbkg) + origini - k * ImgAirm
    # calibrated surface brightness in mag /sq arc sec
    # calSbTot = calMagTot + 2.5 * np.log10(sec2)
    calSbBkg = calMagBkg + 2.5 * np.log10(sec2)

    # plt.figure()
    # plt.imshow(ImgAirm, cmap="inferno", norm = norm )
    # plt.colorbar()
    # plt.title("airmass")

    norm1 = simple_norm(calSbBkg, "sqrt")
    title = band + " background Surface Brightness"
    file = band + "calSbBkg.png"
    plt.figure()
    plt.imshow(-calSbBkg, cmap="inferno")
    plt.colorbar()
    plt.savefig(file)
    plt.title(title)

    file = band + "Stars_Match.png"
    title = band + " Stars correspondance"
    plt.figure()
    plt.plot(StarMatch[:, 2], StarMatch[:, 3], "or", markersize=2)
    plt.plot(StarMatch[:, 0], StarMatch[:, 1], "ok", markersize=2)
    plt.ylim(ny, 0)
    plt.xlim(0, nx)
    plt.title(title)
    plt.legend(["Simbad reference stars", "Detected stars"])
    plt.savefig(file)

    # file = band  + "Sky.png"
    # title = band + " Stars image"
    # plt.figure()
    # plt.imshow(imstars, cmap="inferno" )
    # plt.colorbar()
    # plt.title(title)
    # plt.savefig('Sky.png')

    # determine cloud cover

    threshold = np.amax(imstars) / cfactor
    # stars_binary = np.full((ny,nx),np.nan)
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

    # plt.figure()
    # plt.imshow(stars_count, cmap="inferno" )
    # plt.colorbar()
    #
    # plt.figure()
    # plt.title("stars_full")
    # plt.imshow(stars_full_count, cmap="inferno" )
    # plt.colorbar()
    #
    # plt.figure()
    # plt.title("stars")
    # plt.imshow(imstars, cmap="inferno" )
    # plt.colorbar()

    if band == "V":
        np.save("BackgroundV.npy", calSbBkg)
    elif band == "R":
        np.save("BackgroundR.npy", calSbBkg)

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
            + band
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

# plt.show()
