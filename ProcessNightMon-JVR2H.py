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
from numpy import inf
from scipy import ndimage, stats
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from skyfield.api import E, N, load, wgs84, Star
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
from skyfield import almanac
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import convolve
from astropy.convolution import  Box2DKernel
from astropy.convolution import Gaussian2DKernel
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from toolbox import open_raw, to_grayscale
from datetime import date, timedelta

# find indices of pixels of a list of alt az coordinates
def find_close_indices(array1, array2, value1, value2):
    # 1 = azimuth, 2 = elevation, 3 = radius, all in degrees
    minindex = np.zeros([1], dtype=int)
    for n in range(np.shape(value1)[0]):
        dist = np.sqrt(((array1-value1[n])*np.cos(np.pi*(array2-value2[n])/2/180))**2+(array2-value2[n])**2)
        min = np.nanmin(dist)
        mini = np.where( dist == min )
        minindex = np.append(minindex, mini)
    minindex = minindex[minindex > 0]
    return minindex

def closest_node(node, nodes):
    dist = np.sqrt((nodes[cdist([node], nodes).argmin()][0]-node[0])**2+(nodes[cdist([node], nodes).argmin()][1]-node[1])**2)
    return nodes[cdist([node], nodes).argmin()], dist, cdist([node], nodes).argmin()

def airmass(elevat):
    zenith_rad = (90 - elevat) * np.pi / 180
    AM = ( (1.002432*((np.cos(zenith_rad)) ** 2) +
            0.148386*(np.cos(zenith_rad)) + 0.0096467) /
            (np.cos(zenith_rad) ** 3 +
            0.149864*(np.cos(zenith_rad) ** 2) +
            0.0102963*(np.cos(zenith_rad)) + 0.000303978) )
    return AM

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
# ================================================
# MAIN
# default Parameters
sberr = 0
mflag = "False"
FWHM = 3
norm = ImageNormalize(stretch=SqrtStretch())
elmin = 10  # set the minimum elevation
# calculate maximum zenith angle
zemax = 90 - elmin
# read site coordinates
# Load Parameters
home = os.path.expanduser("~")
with open(home + "/nightmon_config") as f:
    p = yaml.safe_load(f)

# load sky brightness extraction points
pt = open(home + "/points_list","r")
pt.readline() # skip one line
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
time=ts.utc(2020, 9, 22, 4, 32, 53)
#time = ts.now()


basename = time.utc_strftime('%Y-%m-%d')
yesterday = time - timedelta(days = 1)

timestamp = time.utc_strftime('%Y-%m-%dT%H:%M:%S')
Site = p["Site"]
# =======================================
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
moonphase = almanac.moon_phase(eph, time).degrees


# setting moon flag to 1 if its elevation is positive
if (altm > 0):
    mflag = "Thrue"
print(f"Moon altitude: {altm:.4f}")
print(f"Moon azimuth: {azim:.4f}")
print(f"Sun altitude: {alts:.4f}")
print(f"Sun azimuth: {azis:.4f}")
# creating star map with the SIMBAD stars database
limit = 4  # limitting stars magnitude
ds = pd.read_csv(home + "/git/nightmon/data/simbad_lt_6Vmag_r1.8.csv",header=0,sep=";")
stars_selected = ds[ds['MagV'] < limit]
#stars_selected=ds[ds['MagR'] < limit]
polaris = ds.loc[ds["identifier"] == "* alf UMi"]
coordsradec = stars_selected["coord1_ICRS,J2000/2000_"]
radpolaris = polaris["coord1_ICRS,J2000/2000_"]
coords = coordsradec.to_numpy(dtype="str")
coordspolaris = radpolaris.to_numpy(dtype="str")
cpolaris = coordspolaris[0].split(" ")
etoilepol = Star(ra_hours=(float(cpolaris[0]), float(cpolaris[1]), float(cpolaris[2])),dec_degrees=(float(cpolaris[3]), float(cpolaris[4]), float(cpolaris[5])))
epol = here_now.observe(etoilepol).apparent()
alt_star, azi_star, distance = epol.altaz()
alt_pol=alt_star.degrees
azi_pol=azi_star.degrees
altpol = np.zeros(1)
azipol = np.zeros(1)
azipol[0] = azi_pol
altpol[0] = alt_pol
identity = stars_selected["identifier"]
iden = identity.to_numpy(dtype="str")
magnitudev = stars_selected["MagV"]
magv = magnitudev.to_numpy(dtype="float")
magnituder = stars_selected["MagR"]
magr = magnituder.to_numpy(dtype="float")

# create np array for stars
i=0
altstar = np.zeros(np.shape(coords)[0])
azistar = np.zeros(np.shape(coords)[0])
for i in range(np.shape(coords)[0]):
    posstar = coords[i].split(" ")
    rah=float(posstar[0])
    ram=float(posstar[1])
    ras=float(posstar[2])
    ded=float(posstar[3])
    dem=float(posstar[4])
    des=float(posstar[5])
    etoile = Star(ra_hours=(rah, ram, ras),dec_degrees=(ded, dem, des))
    etoi = here_now.observe(etoile).apparent()
    alt_star, azi_star, distance = etoi.altaz()
    alt_star=alt_star.degrees
    azi_star=azi_star.degrees
    azistar[i]=azi_star
    altstar[i]=alt_star
# keep stars above horizon
azistar[ altstar < elmin ] = np.nan
magv[ altstar < elmin ] = np.nan
magr[ altstar < elmin ] = np.nan
altstar[ altstar < elmin ] = np.nan
iden[ altstar < elmin ] = np.nan
# remove stars without R magnitude
azistar[ np.isnan(magr) ] = np.nan
altstar[ np.isnan(magr) ] = np.nan
magv[ np.isnan(magr) ] = np.nan
iden[ np.isnan(magr) ] = "unknown"
# remove nan
azistar = azistar[ ~np.isnan(azistar) ]
altstar = altstar[ ~np.isnan(altstar) ]
magv = magv[ ~np.isnan(magv) ]
magr = magr[ ~np.isnan(magr) ]
iden = np.delete(iden, np.where(iden == "unknown"))
# calculating airmass with A. T. Young, "AIR-MASS AND REFRACTION," Applied Optics, vol. 33,
#    pp. 1108-1110, Feb 1994.
AirM = airmass(altstar)
# =====================================
# sky image analysis
# reading V and R images names from the input parameters
Vfile, Rfile = input(sys.argv[1:])
Vimg = open_raw(Vfile)
Rimg = open_raw(Rfile)
# reading coefficients to convert RGB into gray scale
RC = p["R2GrayCoef"]
GC = p["G2GrayCoef"]
BC = p["B2GrayCoef"]
# create the grayscale images
Vgray = to_grayscale(Vimg, [RC, GC, BC], normalize=False)
Rgray = to_grayscale(Rimg, [RC, GC, BC], normalize=False)
ny = np.shape(Vgray)[0]
nx = np.shape(Vgray)[1]

# **************************** supprimer
rm = 0.635
vv = np.array([0.99903, 0.9990698566039953, 0.9981855774518988, 0.9965508739460722, 0.9943394574888771, 0.9917250394826749, 0.9888760750936606, 0.9859422041998309, 0.9829415915486667, 0.9798735483421788, 0.9767373857823781, 0.9735324150712757, 0.9702579474108823, 0.966913294003209, 0.9634977660502667, 0.9600106747540663, 0.9564513313166186, 0.9528190469399347, 0.9491131328260255, 0.9453329001769019, 0.9414776601945748, 0.9375467240810552, 0.9335394030383539, 0.9294550082684819, 0.9252928509734502, 0.9210522423552696, 0.9167324936159511, 0.9123329159575057, 0.9078528205819441, 0.9032915186912773, 0.8986483214875165, 0.8939225401726723, 0.8891134859487557, 0.8842204700177777, 0.8792428035817492, 0.8741797978426812, 0.8690307640025845, 0.8637950132634701, 0.8584718568273488, 0.8530606058962317, 0.8475605716721296, 0.8419710653570536, 0.8362913981530145, 0.8305208812620231, 0.8246588258860905, 0.8187045432272276, 0.8126573444874453, 0.8065165408687546, 0.8002814435731663, 0.7939513638026915, 0.7875250748828453, 0.7810021458037676, 0.7743799985558575, 0.7676564974881454, 0.7608295069496618, 0.7538968912894375, 0.7468565148565028, 0.7397062419998887, 0.7324439370686252, 0.7250674644117431, 0.717574688378273, 0.7099634733172454, 0.7022316835776909, 0.69437718350864, 0.6863978374591233, 0.6782915097781714, 0.6700560648148147, 0.6616893669180839, 0.6531892804370095, 0.6445536697206221, 0.6357803991179521, 0.6268673329780303, 0.6178123356498872, 0.6086132714825532, 0.599268004825059, 0.589774400026435, 0.580130321435712, 0.5703336334019204, 0.5603822002740906, 0.5502738864012535, 0.5400065561324394, 0.5295780738166791, 0.5189863038030028, 0.5082291104404415, 0.4973043580780254, 0.486209911064785, 0.4749436337497512, 0.4635033904819544, 0.4518870456104251, 0.4400924634841938, 0.4281175084522914, 0.4159600448637481, 0.4028579591928673, 0.3841800074975651, 0.3540685156254719, 0.306664774044784, 0.2361100732236971, 0.1365457036304061, 0.02552628311479194, 0.02002000000000018])
dy = 0
dx = 0
yy=(np.linspace(1,ny,ny)-ny/2.-dy)/ny
xx=(np.linspace(1,nx,nx)-nx/2.-dx)/ny
[xxx,yyy]=np.meshgrid(xx,yy)
rrr=np.sqrt(xxx*xxx+yyy*yyy)/rm
rrr[rrr>1.0]=0.0
r=np.linspace(0,1,100)
flat=np.interp(rrr,r,vv)
# ******************************

# plt.figure()
# plt.imshow(Vgray, cmap="rainbow")
# plt.colorbar()
# plt.title("Original")

Vgray = Vgray / flat

# plt.figure()
# plt.imshow(Vgray, cmap="rainbow")
# plt.colorbar()
# plt.title("Flat corrected ")
#
# plt.figure()
# plt.imshow(flat, cmap="rainbow")
# plt.colorbar()
# plt.title("flat ")

y, x = np.indices((ny, nx))
# computing the distance to zenith in pixels
d = np.hypot(x - nx/2, y - ny/2)
d[d < 0.5] = 0.5
#z = p["Zconstant"] + p["Zslope"] * d + p["Zquad"] * d**2
z =   p["Zslope"] * d + p["Zquad"] * d**2
z[z < 0] = 0
# computing azimuth
az = np.arctan2(-x + nx/2, -y + ny/2) * 180 / np.pi
#az = az - azpol
az = np.where(az < 0, az + 360, az)
# solid angle in sq arc second
sec2 = ((p["Zslope"]+p["Zquad"]) * 180 * 3600**2 * np.sin((z)*np.pi/180))/(d*np.pi)
# compute elevation angle
el = 90 - z





# plt.figure()
# plt.imshow(sec2, cmap="inferno")
# plt.colorbar()
# plt.title("Solid angle ")
# plt.figure()
# plt.imshow(z, cmap="inferno")
# plt.colorbar()
# plt.title("Zenith angle ")
# plt.show()
#
outname= "calibrated" + "_"  + basename + "_sky.csv"
# create the data file if it do not exists
if (os.path.exists(outname) == False):
    o = open(outname, "w")
    first_line =  "# Loc_Name , CCD_XY_position , ,  AzAlt_position , , Airmass ,      \
    Date , Moon , Clouds ,   Johnson_V , err , V_zeropoint , Johnson_R , err , R_zeropoint , Sun_Alt, Moon_Alt , Galactic_Lat , Moon_Phase \n"
    second_line = "# , (pixel) , (pixel) , (deg) , (deg) ,  ,  ,  \
    , (%) , (mag/arcsec^2) ,  , mag , (mag/arcsec^2) , , mag , deg , deg ,  deg , deg , \n"
    o.write(first_line)
    o.write(second_line)
    o.close()

# initialize np arrays for backgrounds and stars fields
Vbkg = np.zeros([ny,nx])
Rbkg = np.zeros([ny,nx])
Vstars = np.zeros([ny,nx])
Rstars = np.zeros([ny,nx])
# find pixel positions of the SIMBAD stars on the array grid corresponding to image size
index  = find_close_indices(az, el, azistar , altstar)
polindex = find_close_indices(az, el, azipol , altpol)
shape=int(np.shape(index)[0])
shape=int(shape/2)
index = index.reshape(shape,2)
# index first column = y and second = x
StarMatch = np.zeros([shape,11])
print("Number of SIMBAD reference stars : ", shape)
# ================================
# loop over the two bands
for band, xpoli, ypoli, imagi, imbkg, imstars   in (
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
    # iterative process to find correct zenith position xzen Yzenith
    xzen=nx/2
    yzen=nx/2
    xpol = xpoli
    ypol = ypoli
    gap = 100
    gaptx = 0
    gapty = 0
    angle = 0
    imag = imagi
    print("Searching zenith in image ",band," ...")
    while gap > 0.5 :
        # calculate polaris azimuth
        angpol =  np.arctan2(xpol - xzen, yzen - ypol) * 180 / np.pi
        angpol = angpol - azi_pol
        xpol = -np.hypot(xzen-xpol, yzen-ypol)*np.sin(angpol*np.pi/180.)+xzen
        ypol = -np.hypot(xzen-xpol, yzen-ypol)*np.cos(angpol*np.pi/180.)+yzen
        gapx = xpol - polindex[1]
        gapy = ypol - polindex[0]
        gap = np.hypot(gapx, gapy)
        xzen = xzen - gapx
        yzen = yzen - gapy
        xpol = xpol - gapx
        ypol = ypol - gapy
        gaptx = gaptx + gapx
        gapty = gapty + gapy
        angle = angle + angpol
    angpol =  np.arctan2(xpoli - xzen, yzen - ypoli) * 180 / np.pi
    angpol = angpol + azi_pol
    shiftx = xzen - nx /2
    shifty = yzen - ny /2
    print("Zenith position: ",xzen,yzen)
    imag = ndimage.shift(imag, [ shiftx,  shifty], mode="nearest")
    imag = ndimage.rotate(imag, angpol, reshape=False, mode="nearest")

    # find background
    sigma_clip = SigmaClip(sigma=2.)
    bkg_estimator = MedianBackground()
    # background image
    bkg = Background2D(imag, (9, 9), filter_size=(int(FWHM), int(FWHM)),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    imbkg = bkg.background
    # image without background
    imstars = imag - imbkg

    mean, median, std = sigma_clipped_stats(imstars, sigma=5.0)
    daofind = DAOStarFinder(fwhm=FWHM, threshold=5.*std)
    sources = daofind(imstars)
    print(sources)
    maxflux = sources["flux"].max()
    #sources = sources[sources['flux'] > maxflux / 5 ]
    # get the flux
    Flux = sources["flux"]
    # keep stars with elevation > elmin degrees
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    radius = np.hypot(positions[:,0]-nx/2,positions[:,1]-ny/2)
    sources = sources[radius > ny/3 ]
    #positions[:, 1] = ny - positions[:, 1]
    apertures = CircularAperture(positions, r=4.)


    ns=0
    n=0
    nistars = shape
    # searching for correspondance between stars in simbad and found stars in image
    dstar=np.zeros(np.shape(positions)[0])
    for ns in range(shape):
        for npos in range(np.shape(positions)[0]):
            dstar[npos] = np.hypot(positions[npos,0]-index[ns,1],positions[npos,1]-index[ns,0])
        dmin = np.amin(dstar)
        dmin_index = dstar.argmin()
        if ( dmin < 20 ):
            #noeuds, gap , ind_noeuds = closest_node(a_star, positions)  # FAUDRAIT QUE CETTE FONCTION RETOURNE LE FLUX
            #print("noeuds",noeuds)
            StarMatch[n,0] = index[ns,1]
            StarMatch[n,1] = index[ns,0]
            StarMatch[n,2] = positions[dmin_index,0] #noeuds[0]
            StarMatch[n,3] = positions[dmin_index,1] #noeuds[1]
            StarMatch[n,4] = dmin
            StarMatch[n,5] = positions[dmin_index,0] - index[ns,1] #noeuds[0] - index[ns,1]
            StarMatch[n,6] = positions[dmin_index,1] - index[ns,0] #noeuds[1] - index[ns,0]
            StarMatch[n,7] = magv[ns]
            StarMatch[n,8] = magr[ns]
            StarMatch[n,9] = AirM[ns]
            StarMatch[n,10] = Flux[dmin_index]
            n=n+1

    nfstars=np.shape(StarMatch)[0]
    Approx_cloud = 1 - nfstars / nistars
    avggap , mediangap, stdgap = sigma_clipped_stats(StarMatch[:,4], sigma=3.0)
    print("Average distance between nearest star: ", avggap)

    am = StarMatch[:,9]
    am = np.delete(am, np.where((StarMatch[:,4] > avggap + 2 * stdgap) & (StarMatch[:,4] < avggap - 2 * stdgap)))
    am = np.delete(am, np.where(StarMatch[:,10] == 0))
    mvtoa = StarMatch[:,7]
    mvtoa = np.delete(mvtoa, np.where((StarMatch[:,4] > avggap + 2 * stdgap) & (StarMatch[:,4] < avggap - 2 * stdgap)))
    mvtoa = np.delete(mvtoa, np.where(StarMatch[:,10] == 0))
    umv = StarMatch[:,10]
    umv = np.delete(umv, np.where((StarMatch[:,4] > avggap + 2 * stdgap) & (StarMatch[:,4] < avggap - 2 * stdgap)))
    umv = np.delete(umv, np.where(StarMatch[:,10] == 0))

    # calculate uncalibrated mag, extinction, calibration factor
    uncalMag = -2.5 * np.log10(umv)
    # print(uncalMag)
    # print(mvtoa)
    # print(Flux)
    deltam = uncalMag - mvtoa

    # filter outliers zscore = (x - mean)/std
    df = pd.DataFrame(zip(deltam, am))
    df = df[(np.abs(stats.zscore(df)) < 2.5).all(axis=1)]
    co = np.polyfit(df[1],df[0],1)

    slope = co[0]
    origin = co[1]
    zeropoint = -origin


    testcal = -2.5 * np.log10(umv) - origin
    plt.figure()
    plt.plot(mvtoa,testcal,"or",markersize=2)
    plt.title("verif calib")


    print("Slope = ",slope," Origin = ", origin)
    fx = np.linspace(0,np.amax(am),100)
    fy = co[1] + fx * co[0]
    # atmospheric optical depth (aerosols + molecules)
    tau = slope / (2.5 * np.log10(np.e))
    extinction = slope
    #  molecular optical depth - Bodhaine et al. (1999)
    if ( band == "V"):
        lambm = 551
    elif (band == "R"):
        lambm = 658
    tau_mol = (0.0021520*(1.0455996-341.29061*(lambm/1000.)**
        -2.-0.90230850*(lambm/1000.)**2.)/(1.+0.0027059889*(lambm/1000.)**
        -2.-85.968563*(lambm/1000.)**2.))
    AOD = tau-tau_mol
    print("Atmospheric optical depth =", tau)
    print("Aerosol optical depth =",AOD," at ",lambm," nm")
    print("Zenith atmospheric extinction (mag) =",extinction)
    print("Zero point (mag) =",zeropoint)

    # calculate a calibrated magnitude
    # calMag = -2.5 * np.log10(pixel_value) - origin
    calMagTot = -2.5 * np.log10(imag) - origin
    calMagBkg = -2.5 * np.log10(imbkg) - origin
    # calibrated surface brightness in mag /sq arc sec
    calSbTot = calMagTot + 2.5 * np.log10(sec2)
    calSbBkg = calMagBkg + 2.5 * np.log10(sec2)

    title = band + " background Surface Brightness"
    file =  band + "calSbBkg.png"
    plt.figure()
    plt.imshow(calSbBkg, cmap="inferno", norm = norm )
    plt.zlim(22, 16)
    plt.colorbar()
    plt.savefig(file)
    plt.title(title)

    file = band  + "zeropoint_corr.png"
    title = band + " delta_Mag vs Airmass"
    plt.figure()
    plt.plot(df[1],df[0],"or",markersize=2)
    plt.plot(fx,fy,'b')
    plt.title(title)
    plt.savefig(file)

    file = band  + "Stars_Match.png"
    title = band + " Stars correspondance"
    plt.figure()
    plt.plot(StarMatch[:,2],StarMatch[:,3],"or",markersize=2)
    plt.plot(StarMatch[:,0],StarMatch[:,1],"ok",markersize=2)
    plt.ylim(ny, 0)
    plt.xlim(0,nx)
    plt.title(title)
    plt.legend(["Simbad reference stars", "Detected stars"])
    plt.savefig(file)

    # file = band  + "Sky.png"
    # title = band + " Sky image"
    # plt.figure()
    # plt.imshow(imag, cmap="inferno", norm = norm )
    # plt.colorbar()
    # plt.title(title)
    # plt.savefig('Sky.png')

    # determine cloud cover
    threshold = 0.001
    stars_binary = np.full((ny,nx),np.nan)
    stars_full = np.full((ny,nx), 1.)
    stars_binary[imstars>=threshold] = 1.
    stars_binary[z > 90] = np.nan
    window = 11  #  1 deg clouds ~= 10
    kernel=Box2DKernel(width=window, mode='integrate')
    stars_count = convolve(stars_binary, kernel)
    stars_count[z > 90] = 0.
    stars_full[z > 90] = 0.
    stars_count=np.nan_to_num(stars_count,nan=0.0)
    stars_full=np.nan_to_num(stars_full,nan=0.0)
    # weighted with solid angle
    cloud_cover=1-np.sum(stars_count*sec2)/np.sum(stars_full*sec2)
    print("Cloud fraction : ",cloud_cover, Approx_cloud)

    if ( band == "V"):
        calSbBkgV = calSbBkg
        cloud_coverV = cloud_cover
        Approx_cloudV = Approx_cloud
        extinctionV = extinction
        AODV = AOD
        zeropointV = zeropoint
        os.system('mv zeropoint_corr.png zeropoint_corrV.png')
        errV = sberr
        # saving transformed sky image
        os.system('mv Sky.png SkyV.png')
        # saving sky background
        np.save("BackgroundV.npy", calSbBkgV)
        os.system('mv calSbBkg.png calSbBkgV.png')
        os.system('mv Stars_Match.png Stars_MatchV.png')

    elif ( band == "R" ):
        calSbBkgR = calSbBkg
        cloud_coverR = cloud_cover
        Approx_cloudR = Approx_cloud
        extinctionR = extinction
        AODR = AOD
        zeropointR = zeropoint
        os.system('mv zeropoint_corr.png zeropoint_corrR.png')
        errR = sberr
        # saving transformed sky image
        os.system('mv Sky.png SkyR.png')
        # saving sky background
        np.save("BackgroundR.npy", calSbBkgR)
        os.system('mv calSbBkg.png calSbBkgR.png')
        os.system('mv Stars_Match.png Stars_MatchR.png')

# Extract points and write data
index  = find_close_indices(az, el, apt , ept)
shape=int(np.shape(index)[0])
shape=int(shape/2)
index = index.reshape(shape,2)

for no in range(num_pts):
    posx=str(index[no,1])
    posy=str(index[no,0])
    direction = here_now.from_altaz(alt_degrees=ept[no], az_degrees=apt[no])
    glat, glon, gdistance = direction.frame_latlon(galactic_frame)
    galactic_lat = glat.degrees
    theta_sun = direction.separation_from(sun_position).degrees
    theta_moon = direction.separation_from(moon_position).degrees
    # calculate Airmass
    airmo = airmass(ept[no])
    # read V sky Brightness
    Vmago = calSbBkgV[index[no,0],index[no,1]]
    # read R sky Brightness
    Rmago = calSbBkgR[index[no,0],index[no,1]]
    cloud_covero = (cloud_coverV + cloud_coverR) / 2
    Approx_cloudo = (Approx_cloudV + Approx_cloudR) / 2
    # check reliability of the cloud cover values


    o = open(outname, "a")
    outputline = Site + " , " + posx + " , " + posy + " , " + \
    str("{:6.2f}".format(apt[no])) + " , " + \
    str("{:5.2f}".format(ept[no])) + " , " + \
    str("{:3.1f}".format(airmo))  + " , " + timestamp + " , " + mflag +  " , " + \
    str("{:5.1f}".format(cloud_covero*100)) + " , " + \
    str("{:6.3f}".format(Vmago)) + " , " + \
    str("{:6.3f}".format(errV)) + " , "  + \
    str("{:6.3f}".format(zeropointV)) + " , "  + \
    str("{:6.3f}".format(Rmago)) + " , " + \
    str("{:6.3f}".format(errR)) + " , " + \
    str("{:6.3f}".format(zeropointV)) + " , "  + \
    str("{:6.2f}".format(theta_sun)) + " , " + \
    str("{:6.2f}".format(theta_moon)) + " , " + \
    str("{:6.2f}".format(galactic_lat)) + " , " + \
    str("{:6.2f}".format(moonphase)) + "\n"
    print(band," : ",outputline)
    o.write(outputline)
    o.close()

plt.show()
