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
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import yaml
from scipy import ndimage
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from skyfield.api import E, N, load, wgs84, Star
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from toolbox import open_raw, to_grayscale

# find indices of pixels within a radius of a list of alt az coordinates
def find_close_indices(array1, array2, value1, value2, value3):
    # 1 = azimuth, 2 = elevation, 3 = radius, all in degrees
    indices = np.zeros([1], dtype=int)
    minindex = np.zeros([1], dtype=int)
    for n in range(np.shape(value1)[0]):
        dist = np.sqrt(
            ((array1 - value1[n]) * np.cos(np.pi * (array2 - value2[n]) / 2 / 180)) ** 2
            + (array2 - value2[n]) ** 2
        )
        min = np.nanmin(dist)
        indi = np.where(dist < value3)
        mini = np.where(dist == min)
        indices = np.append(indices, indi)
        minindex = np.append(minindex, mini)
    indices = indices[indices > 0]
    minindex = minindex[minindex > 0]
    return indices, minindex


# default Parameters
mflag = 0
# read site coordinates
# Load Parameters
home = os.path.expanduser("~")
with open(home + "/nightmon_config") as f:
    p = yaml.safe_load(f)

# load sky brightness extraction points
with open(home + "/points_list") as pt:
    pts = yaml.safe_load(pt)
    pt.seek(0)
    num_pts = sum(1 for line in pt)
    num_pts = num_pts - 1
    print("Number of points to extract :", num_pts)

ts = load.timescale()
time = ts.utc(2020, 9, 22, 4, 22, 53)
# time = ts.now()
basename = time.utc_strftime("%Y%m%d")
timestamp = time.utc_strftime("%Y%m%d_%H%M%S")
# reading command line parameters
def input(argv):
    Vfile = "undefined"
    Rfile = "undefined"
    try:
        opts, args = getopt.getopt(argv, "h:v:r:", ["help=", "vfile=", "rfile="])
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
    print("Johnson V file is ", Vfile)
    print("Johnson R file is ", Rfile)
    return Vfile, Rfile


# astronomical objects and ephemerides
eph = load("de421.bsp")
here = eph["earth"] + wgs84.latlon(
    p["Latitude"] * N, p["Longitude"] * E, elevation_m=p["Altitude"]
)
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
# setting moon flag to 1 if its elevation is positive
if altm > 0:
    mflag = 1

print(f"Moon altitude: {altm:.4f}")
print(f"Moon azimuth: {azim:.4f}")
print(f"Sun altitude: {alts:.4f}")
print(f"Sun azimuth: {azis:.4f}")

Vfile, Rfile = input(sys.argv[1:])
Vimg = open_raw(Vfile)
Rimg = open_raw(Rfile)
# Vimg = 65535.*Vimg
# Rimg = 65535.*Vimg
RC = p["R2GrayCoef"]
GC = p["G2GrayCoef"]
BC = p["B2GrayCoef"]
Vgray = to_grayscale(Vimg, [RC, GC, BC], normalize=False)
Rgray = to_grayscale(Rimg, [RC, GC, BC], normalize=False)


ny = np.shape(Vgray)[0]
nx = np.shape(Vgray)[1]
Vbkg = np.zeros([ny, nx])
Rbkg = np.zeros([ny, nx])
Vstars = np.zeros([ny, nx])
Rstars = np.zeros([ny, nx])
for band, xzen, yzen, xpol, ypol, imag, imbkg, imstars, outf in (
    (
        "V",
        p["XzenithV"],
        p["YzenithV"],
        p["XpolarisV"],
        p["YpolarisV"],
        Vgray,
        Vbkg,
        Vstars,
        "VImage.npy",
    ),
    (
        "R",
        p["XzenithR"],
        p["YzenithR"],
        p["XpolarisR"],
        p["YpolarisR"],
        Rgray,
        Rbkg,
        Rstars,
        "RImage.npy",
    ),
):
    print(f"Processing Johnson {band} camera...")
    y, x = np.indices((ny, nx))
    # computing the distance to zenith in pixels
    d = np.hypot(x - xzen, y - yzen)
    z = p["Zconstant"] + p["Zslope"] * d + p["Zquad"] * d**2
    # calculate azimuth
    azpol = np.arctan2(xpol - xzen, yzen - ypol) * 180 / np.pi
    print("Rotation angle", azpol)
    az = np.arctan2(-x + xzen, -y + yzen) * 180 / np.pi
    az = az - azpol
    az = np.where(az < 0, az + 360, az)
    # shift images (must be replaced by a 3d rotation)
    az = ndimage.shift(az, [(ny / 2 - yzen), (nx / 2 - xzen)], mode="nearest")
    z = ndimage.shift(z, [(ny / 2 - yzen), (nx / 2 - xzen)], mode="nearest")
    imag = ndimage.shift(imag, [(ny / 2 - yzen), (nx / 2 - xzen)], mode="nearest")
    # rotate images
    az = ndimage.rotate(az, -azpol, reshape=False, mode="nearest")
    z = ndimage.rotate(z, -azpol, reshape=False, mode="nearest")
    imag = ndimage.rotate(imag, -azpol, reshape=False, mode="nearest")
    # find background
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    # background image
    bkg = Background2D(
        imag,
        (5, 5),
        filter_size=(3, 3),
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
    )
    imbkg = bkg.background
    # image without background
    imstars = imag - imbkg
    mean, median, std = sigma_clipped_stats(imstars, sigma=3.0)

    daofind = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
    sources = daofind(imstars)
    sources = sources[sources["peak"] > 0.5]

    # keep stars wit elevation > elmin degrees
    elmin = 30
    zemax = 90 - elmin
    sources = sources[
        np.sqrt(
            (sources["xcentroid"] - nx / 2) ** 2 + (sources["ycentroid"] - ny / 2) ** 2
        )
        < np.sqrt((xzen - xpol) ** 2 + (yzen - ypol) ** 2)
        * zemax
        / (90 - p["Latitude"])
    ]
    for col in sources.colnames:
        sources[col].info.format = "%.8g"  # for consistent table output
    print(sources)
    positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))
    # positions[:, 1] = ny - positions[:, 1]
    apertures = CircularAperture(positions, r=4.0)

    # saving sky image
    np.save(outf, imag)
    # mask non sense zenith angles
    # imag[z > 90] = np.nan
    # az[z > 90] = np.nan
    # imbkg[z > 90] = np.nan
    # imstars[z > 90] = np.nan
    # TODO : determine clear fraction - here we need to find a way to detect clouds
    cfrac = 0.0
    # compute elevation angle
    el = 90 - z
    direction = here_now.from_altaz(alt_degrees=90 - z, az_degrees=az)
    glat, glon, gdistance = direction.frame_latlon(galactic_frame)
    galactic_lat = glat.degrees
    theta_sun = direction.separation_from(sun_position).degrees
    theta_moon = direction.separation_from(moon_position).degrees
    norm = ImageNormalize(stretch=SqrtStretch())

    # plt.figure()
    # plt.imshow(imag, cmap="inferno", norm = norm )
    # plt.colorbar()
    # plt.title("V image")
    #
    # plt.figure()
    # plt.imshow(imbkg, cmap="inferno", norm = norm)
    # plt.colorbar()
    # plt.title("background")

    # plt.figure()
    # plt.imshow(imstars, cmap="inferno", norm = norm)
    # apertures.plot(color='white', lw=1.5, alpha=0.5)
    # plt.colorbar()
    # plt.title("stars")

    # plt.figure()
    # plt.imshow(imag, cmap="rainbow")
    # plt.colorbar()
    # plt.title("Sky Image : Johnson " + band)
    #
    # plt.figure()
    # plt.imshow(az, cmap="rainbow")
    # plt.colorbar()
    # plt.title("Azimuth angle " + band)
    #
    # plt.figure()
    # plt.imshow(z, cmap="rainbow")
    # plt.colorbar()
    # plt.title("Zenith angle " + band)
    #
    # plt.figure()
    # plt.imshow(el, cmap="rainbow")
    # plt.colorbar()
    # plt.title("Elevation angle " + band)

# creating star map
# with load.open(hipparcos.URL) as f:
#     df = hipparcos.load_dataframe(f)
#
# df = df[df['ra_degrees'].notnull()]
# df = df[df['magnitude'] <= 2.87] # 2.87 is about 170 highest v apparent luminosity
# bright_stars = Star.from_dataframe(df)
# astrometric = here_now.observe(bright_stars).apparent()
# ra, dec, distance = astrometric.radec()
# hip_num = df.index.values


# creating star map with the SIMBAD database
limit = 3.5  # limitting stars magnitude
ds = pd.read_csv(
    home + "/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv", header=0, sep=";"
)
stars_selected = ds[ds["MagV"] < limit]
stars_selected = ds[ds["MagR"] < limit]
coordsradec = stars_selected["coord1_ICRS,J2000/2000_"]
coords = coordsradec.to_numpy(dtype="str")
magnitudev = stars_selected["MagV"]
magv = magnitudev.to_numpy(dtype="float")
magnituder = stars_selected["MagR"]
magr = magnituder.to_numpy(dtype="float")
# create np array for stars
i = 0
altstar = np.zeros(np.shape(coords)[0])
azistar = np.zeros(np.shape(coords)[0])
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

# keep stars above horizon
azistar[altstar < 30] = np.nan
magv[altstar < 30] = np.nan
magr[altstar < 30] = np.nan
altstar[altstar < 30] = np.nan
# remove stars without R magnitude
azistar[np.isnan(magr)] = np.nan
altstar[np.isnan(magr)] = np.nan
magv[np.isnan(magr)] = np.nan

# remove nan
azistar = azistar[~np.isnan(azistar)]
altstar = altstar[~np.isnan(altstar)]
magv = magv[~np.isnan(magv)]
magr = magr[~np.isnan(magr)]


index, pindex = find_close_indices(az, el, azistar, altstar, 1)
pshape = int(np.shape(pindex)[0])
pshape = int(pshape / 2)
pindex = pindex.reshape(pshape, 2)
print(pindex)
# match database stars with detected stars


def closest_node(node, nodes):
    print(cdist([node], nodes).argmin(), nodes[cdist([node], nodes).argmin()], node)
    dist = np.sqrt(
        (nodes[cdist([node], nodes).argmin()][0] - node[0]) ** 2
        + (nodes[cdist([node], nodes).argmin()][1] - node[1]) ** 2
    )
    print("dist=", dist)
    return nodes[cdist([node], nodes).argmin()]


falseskyV = np.zeros([ny, nx], dtype=float)
falseskyR = np.zeros([ny, nx], dtype=float)
for ns in range(pshape):
    print(pindex[ns, 0], pindex[ns, 1])
    a_star = (pindex[ns, 0], pindex[ns, 1])
    closest_node(a_star, positions)

    falseskyV[pindex[ns, 0], pindex[ns, 1]] = magv[ns]
    falseskyR[pindex[ns, 0], pindex[ns, 1]] = magr[ns]

plt.figure()
plt.plot(positions[:, 0], positions[:, 1], "or")
plt.plot(pindex[:, 1], pindex[:, 0], "ob")

# falseskyV[falseskyV == 0] = np.nan
# falseskyR[falseskyR == 0] = np.nan
# falseskyV[el < 0] = np.nan
# falseskyR[el < 0] = np.nan
# plt.figure()
# plt.imshow(falseskyV, cmap="inferno", norm=norm)
# plt.colorbar()
# plt.title("Database stars in V band")
#
# plt.figure()
# plt.imshow(falseskyR, cmap="inferno", norm=norm)
# plt.colorbar()
# plt.title("Database stars in R band")

# Extract and write data
sberr = 0
# mflag = moon flag 1 moon in the sky 0 no moon
# cfrac = cloud fraction
pazi = np.zeros([1], dtype=float)
pele = np.zeros([1], dtype=float)
prad = np.zeros([1], dtype=float)
for no in range(num_pts):
    pos = "Pos" + str(no)
    posi = np.fromstring(pts[pos], dtype=float, sep=" ")
    pazi[0] = posi[0]
    pele[0] = posi[1]
    prad[0] = posi[2]
    if prad[0] < 0.4:
        prad[0] = 0.4
    index, pindex = find_close_indices(az, el, pazi, pele, prad)
    posx = int(pindex[1])
    posy = int(pindex[0])

    for band, imag, zerop in (
        (
            "V",
            Vgray,
            p["VZeroPoint"],
        ),
        (
            "R",
            Rgray,
            p["RZeroPoint"],
        ),
    ):
        outname = "SBJohnson" + "_" + band + basename + pos + ".dat"
        o = open(outname, "a")
        rad = np.mean(imag[index])
        # TODO : corriger la radiance pour le temps d'integration et le iso rad = rad / tint / iso
        sb = zerop * 10 ** (-0.4 * rad)
        sberr = 0.0
        outputline = (
            timestamp
            + " "
            + str(posx)
            + " "
            + str(posy)
            + " "
            + str("{:.6f}".format(pazi[0]))
            + " "
            + str("{:.6f}".format(pele[0]))
            + " "
            + str("{:.6f}".format(sb))
            + " "
            + str("{:.6f}".format(sberr))
            + " "
            + str(mflag)
            + " "
            + str("{:.6f}".format(cfrac))
            + " "
            + str("{:.6f}".format(np.mean(theta_sun[pindex[0], pindex[1]])))
            + " "
            + str("{:.6f}".format(np.mean(theta_moon[pindex[0], pindex[1]])))
            + " "
            + str("{:.6f}".format(np.mean(galactic_lat[pindex[0], pindex[1]])))
            + " "
            + str("{:.6f}".format(prad[0]))
            + "\n"
        )
        print(band, " : ", outputline)
        o.write(outputline)
        o.close()


# plt.figure()
# plt.imshow(np.round(theta_sun), cmap="rainbow")
# plt.colorbar()
# plt.title("Solar angle")

# plt.figure()
# plt.imshow(np.round(theta_moon), cmap="rainbow")
# plt.colorbar()
# plt.title("Moon angle")

# plt.figure()
# plt.imshow(np.round(galactic_lat), cmap="rainbow")
# plt.colorbar()
# plt.title("Galactic Latitude")
#
# plt.figure()
# plt.imshow(np.round(az), cmap="rainbow")
# plt.colorbar()
# plt.title("Azimuth")
#
# plt.figure()
# plt.imshow(np.round(z), cmap="rainbow")
# plt.colorbar()
# plt.title("Zenith angle")
#
plt.show()
