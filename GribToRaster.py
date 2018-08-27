'''
Version 0.1
Created on February 1, 2018
@author: Shane Motley

@Purpose: 1) Pull data from grib file into numpy array
          2) Convert data from a given coordinate system into lat/lng.
          3) Given a polygon (e.g. river basin), determine which points in the raster file are in the polygon
          4) Create a numpy mask and set points outside the polygon to zero and perform necessary calculations.
          5) Create any plots needed.

@Usage: User can call program and get the forecast for the following calls:
        -l: for a single point (lat / lon). EX: -l 34.5 -104.3
        -hr: hourly (will interpolate forecast hours to generate hourly values)
        -help: help menu


@VersionHistory:
    v0.1: 2/01/18 Beta version
'''

from osgeo import osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from pyproj import Proj
import os, sys
import numpy as np
import scipy.ndimage
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from datetime import datetime
from PIL import Image
import time
import copy
import pytz
import fiona
import argparse
import logging
import GRIB_Download as gd

class Grib():
    '''
    class object for holding the grib data
    '''

    def __init__(self):
        self.model = ""
        self.modelRun = ""
        self.date = np.array([])
        self.variable = ""
        self.level = ""
        self.acc = "" #0-3 hr accumulation or 1hr forecast
        self.gribAll = ""
        self.units = ""
        self.basin = ""
        self.displayunits = ""
        self.mask = np.array([])
        self.data = np.array([])
        self.data4326 = np.array([])
        self.basinTotal = np.array([])
        self.hours = np.array([])

#you will want to remove "type", this is just for testing.
def handle_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('model',
                    nargs="*",
                    default='namNest',
                    help = 'Avaliable Models: nam, namNest, gfs, hrrr')

    parser.add_argument('modelRun',
                        nargs="*",
                        default='06',
                        help = 'Time (in UTC) of model run. For nam, namNest, and gfs models, the times to chose '
                               'are 00, 06, 12, and 18')

    parser.add_argument('date',
                    nargs="*",
                    default=time.strftime("%Y%m%d"),
                    help = 'Date of model run in YYYYMMDD format')  # This will be today's date in yyyymmdd format)

    parser.add_argument('-v','--variable',
                        dest='variable', default='APCP03',
                        required=False,
                        help='The grib file variable to examine (EX: -v APCP03)')

    parser.add_argument('-b', '--basin',
                        dest='basin', default='French_Meadows',
                        required=False,
                        help='The basin to calculate Total SWE. Options inlude:\n'
                             'Hell_Hole, French_Meadows, or MFP')

    parser.add_argument('-u', '--displayunits',
                        dest='displayunits', default='US',
                        required=False,
                        help='Show Units In US or SI')

    parser.add_argument('-l', '--level',
                        dest='level', default='0-SFC',
                        required=False,
                        help='Variable level EX: -l 0-SFC ')

    parser.add_argument('-m', '--map',
                        dest='map', default=True, required=False,
                        help='use this option if you want to output a matplotLib map')

    parser.add_argument('-mt', '--mbtiles',
                        dest='mbtiles', default=True, required=False,
                        help='use this option if you want to output tiles')

    parser.add_argument('-dt', '--tiles4326',
                        dest='tiles4326', default=True, required=False,
                        help='use this option if you want to output tiles')

    parser.add_argument('-p', '--plot',
                        dest='plot', default=True, required=False,
                        help='use this option if you want to output a plot')

    parser.add_argument('-pt', '--pointvalue',
                        dest='point_value', default=None, required=False,
                        help='Value of raster at a lng/lat. Useage: -pt -120.399 39.1339')

    #parser.add_argument('-hr', '--hours',
    #                    dest='hourBlock', default=24, action='store', type=int,
    #                    required=False, nargs=2,
    #                    help='Integer: 12 or 24. Display Forecast as Sunday, Monday or Sunday, Sunday Night, Monday etc')

    #parser.add_argument('--startDay',
    #                    dest='startDay', default=0, action='store', type=int,
    #                    required=False,
    #                    help='Integer: Number of days to advance day1. Default is 0; today'
    #                    )


    args = parser.parse_args()

    return args

def main():
    global inputArgs, grib, logger, dir_path
    dir_path = os.path.dirname(os.path.realpath(__file__))

    logger = logging.getLogger()
    inputArgs = handle_args(sys.argv)
    grib = Grib()
    grib.model = inputArgs.model
    grib.level = inputArgs.level
    grib.variable = inputArgs.variable
    grib.displayunits = inputArgs.displayunits
    grib.modelRun = inputArgs.modelRun

    ##########Debugging############
    # Hard Coding For Debugging
    #inputArgs.date = '20180222'
    inputArgs.mbtiles = True
    inputArgs.tiles4326 = True
    grib.model = 'namNest'
    grib.displayunits = "US"
    inputArgs.point_value = (-120.4388, 39.08109)  # French Meadows User defined: Used to get value of the raster at a given lon/lat
    #inputArgs.point_value = (-121.0547, 38.921513) # PCWA User defined: Used to get value of the raster at a given lon/lat
    #inputArgs.point_value = (38.727543,-121.181441) #9380 Oak Leaf
    ########End Debugging##########

    grib.basin = inputArgs.basin

    ####################################################################################################################
    ###################################### PARAMETER SELECTION #########################################################
    ####################################################################################################################
    grib.acc = False      #Accumulation Parameter. When True, each raster will be a sum of previous fcst hours.
    #Precipitation Parameters
    #grib.variable, grib.level, grib.acc = "APCP",'0-SFC',True   # Accumulated Precipitation [kg/m^2] (Note: kg/m^2 = mm)
    #grib.variable, grib.level, grib.acc = "BGRUN", '0-SFC',True # Baseflow-Groundwater Runoff [kg/m^2]
    #grib.variable, grib.level, grib.acc = "SSRUN", '0-SFC',True # Storm Surface Runoff [kg/m^2]
    #grib.variable, grib.level, grib.acc = "EVP", '0-SFC',True   # Evaporation [kg/m^2]
    #grib.variable, grib.level, grib.acc = "WATR", '0-SFC',True  # GFS Parameter: Water Runoff [kg/m^2]

    #Snowfall / Depth Parameters
    #grib.variable, grib.level, grib.acc = "SNOM",'0-SFC',True   # Snow Melt [kg/m^2]
    #grib.variable, grib.level = "WEASD",'0-SFC'                 # Water Equivalent of Accumulated Snow Depth [kg/m^2]
    #grib.variable, grib.level = "SNOD",  '0-SFC'                # Snow Depth [m]
    #grib.variable, grib.level = "CSNOW", '0-SFC'                # Categorical Snow [-]

    #Radar / Satellite Simulation
    #grib.variable, grib.level = "REFD",  '1000-HTGL'            # Simulated Reflectivity [dBz]
    #grib.variable, grib.level = "MAXREF",'1000-HTGL'            # Hourly Maximum of Simulated Reflectivity at 1 km AGL [dB]
    #grib.variable, grib.level = "BRTMP", '0-NTAT'               # Brightness Temperature [K]
    #grib.variable, grib.level = "TMP", '0-CTL'                  # Cloud Top Temperature [K]

    #Convection
    #grib.variable, grib.level = "LTNG",  '0-SFC'                # Lightning [non-dim]
    #grib.variable, grib.level = "MXUPHL", '1000-HTGL'           # Hourly Maximum of Updraft Helicity over Layer 2km to 5 km AGL [m^2/s^2]
    #grib.variable, grib.level = "PWAT",  '0-EATM'               # Precipitable Water [kg/m^2]
    #grib.variable, grib.level = "CAPE",  '0-SFC'                # Convective Available Potential Energy [J/kg]

    #Surface Temperature / Dewpoint / Wind
    #grib.variable, grib.level = "TMAX",  '2-HTGL'               # Maximum Temperature [K]
    #grib.variable, grib.level = "TMIN",  '2-HTGL'               # Maximum Temperature [K]
    grib.variable, grib.level = "TMP",   '0-SFC'                 # Temperature [K]
    #grib.variable, grib.level = "DPT",   '2-HTGL'               # Dew Point Temperature [K]
    #grib.variable, grib.level = "GUST",  '0-SFC'                # Wind Speed (Gust) [m/s]

    #Solar Radiation
    #grib.variable, grib.level = "CSDSF", '0-SFC'                # Clear Sky Downward Solar Flux [W/m^2]
    #grib.variable, grib.level = "DSWRF", '0-SFC'                # Downward Short-Wave Radiation Flux [W/m^2]
    ####################################################################################################################
    ###################################### END PARAMETER SELECTION #####################################################
    ####################################################################################################################

    # Bounding box will clip the raster to focus in on a region of interest (e.g. CA) This makes the raster MUCH smaller
    # and easier to work with. See gdal.Open -> gdal.Translate below for where this is acutally used.
    # grib.bbox = [-125.0, 50.0, -115.0, 30.0]        # [upper left lon, upper left lat, lower left lon, lower left lat]
    # grib.bbox = [-123.75, 40.9798, -112.5, 31.9522]  # tile 5/5/19
    grib.bbox = [-134.12207, 52.628703, -60.9284436, 21.146938]  # namNest Bounds
    #grib.bbox = [-130.12207, 52.628703, -80.9284436, 21.146938]  # random Bounds
    # grib.bbox = [225.87793, 52.628703, 299.0715564, 21.146938]  # namNest Bounds for gfs

    grib.hours = model_hours()                      #Gets the forecast hours avalible for a given model.

    grib = gribData(grib)                           #Fills the grib object with data (including metadata).
    projInfo = grib.gribAll.GetProjection()         #Projection of grib file
    geoinformation = grib.gribAll.GetGeoTransform() #[0] x min, [1] x resolution, [3] ymax, [5] yres

    xres = geoinformation[1]
    yres = geoinformation[5]
    xmin = geoinformation[0]
    xmax = geoinformation[0] + (xres * grib.gribAll.RasterXSize)
    ymin = geoinformation[3] + (yres * grib.gribAll.RasterYSize)
    ymax = geoinformation[3]

    spatialRef = osr.SpatialReference()             #Get spatial reference for converting this data into lon/lat (or x/y if it's already in lon/lat)
    spatialRef.ImportFromWkt(projInfo)
    spatialRefProj = spatialRef.ExportToProj4()

    # create a grid of xy (or lat/lng) coordinates in the original projection
    xy_source = np.mgrid[xmin:xmax:xres, ymax:ymin:yres]
    xx, yy = xy_source

    # A numpy grid of all the x/y into lat/lng
    # This will convert your projection to lat/lng (it's this simple).
    lons, lats = Proj(spatialRefProj)(xx, yy, inverse=True)

    mxlat, minlat = np.max(lats), np.min(lats)
    mxlon, minlon = np.max(lons), np.min(lons)
    test = makePNG()
    # Find the center point of each grid box.
    # This says move over 1/2 a grid box in the +x direction and move down (since yres is -) in the
    # y direction. Also, the +yres (remember, yres is -) is saying the starting point of this array will
    # trim off the y direction by one row (since it's shifted off the grid)
    xy_source_centerPt = np.mgrid[xmin + (xres / 2):xmax:xres, ymax + (yres / 2):ymin:yres]
    xxC, yyC = xy_source_centerPt

    grib.mask = createMask(xxC, yyC, spatialRefProj)
    if isinstance(grib.mask, ndarray):
        grib.basinTotal = calculateBasin(grib, xres, yres)
        if grib.basin == 'Hell_Hole':  # Part of this basin is SMUD's teritory, so remove 92% of water in this basin
            grib.basin = 'Hell_Hole_SMUD'  # This is just to get the correct directory structure
            grib.mask = createMask(xxC, yyC, spatialRefProj) #Overwrite to get new mask.
            smudBasinTotal = calculateBasin(grib, xres, yres)
            print(
                "Extracting 92% of the SWE values from SMUD Basin...\n" + "Current Basin Total: " + str(
                    grib.basinTotal))
            grib.basinTotal = np.asarray(grib.basinTotal) - (0.92 * np.asarray(smudBasinTotal))
            print("Smud Total: " + str(smudBasinTotal) + "\n New Total: " + str(grib.basinTotal))
            grib.basin = 'Hell_Hole'  # reset back
        print(grib.basinTotal)
    else:
        grib.basinTotal = None
        Warning("Basin Shapefile not found. A basin total can not be created, continuing...")


    if inputArgs.map == True:
        fig = plt.figure(figsize=(8, 8))
        #m = Basemap(llcrnrlon=-125.25, llcrnrlat=36,
        #            urcrnrlon=-119.25, urcrnrlat=41.75,
        #            resolution='i', projection='lcc',
        #            lat_1=38.5, lat_2=38.5, lat_0=38.5, lon_0=262.5)

        m = Basemap(llcrnrlon=-122.8, llcrnrlat=37.3,
                    urcrnrlon=-119.0, urcrnrlat=40.3)

        # plot the data (first layer)
        #lons = scipy.ndimage.zoom(lons, 3)
        #lats = scipy.ndimage.zoom(lats, 3)
        #m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=2000, verbose=True)
        value = 0
        for hr in grib.hours:  # start at hour 3 and go to the last hour in increments of 3
            if inputArgs.point_value != None:
                userLon, userLat = inputArgs.point_value
                value = value_at_point(xx,yy, userLon, userLat, hr, spatialRefProj,m) + value
                # plot on map
                m.plot(userLon, userLat, 'bo')
                plt.annotate(str(round(value * 0.03937,2))+'"', m(userLon-0.15, userLat-0.12))
            makeMap(lons, lats, hr, m)
    if inputArgs.plot == True:
        makePlot()
    return

def gribData(gribObj):
    baseFolder = os.path.join(dir_path,'grib_files')
    grib_dir = os.path.join(baseFolder, inputArgs.date, gribObj.model)
    downloadReset = True  # This prevents an infinite loop from occurring.
    grib_var = []
    grib_var4326 = []
    grib_date = []
    i=0
    if grib.model == 'gfs':
        acc_hrs = 6 #The gfs does precip accumulation every six hours
    else:
        acc_hrs = 3 #All other models do precip accumulation every six hours
    # Remove hour 0 from the list since we don't want the analysis hour.
    del grib.hours[0]
    while i < len(gribObj.hours):
        hr = gribObj.hours[i]
        i += 1
        #Put Hours in correct format
        if hr < 10:
            hr = '0' + str(hr)
            if grib.model == 'gfs': #GFS reports hours as ### (e.g. 003, 006, 009 etc).
                hr = '0' + str(hr)
        else:
            if grib.model == 'gfs' and hr < 100:
                hr = '0' + str(hr)

        #Get file name (if this is an accumulation parameter, get every 3rd hour)
        try:
            if gribObj.acc == False:
                filename = os.path.abspath(
                    grib_dir + '/' + os.fsdecode(inputArgs.date + '_' + inputArgs.modelRun + '_' + str(hr) + '.grib'))
            else: # for precip, we want every third hour (the total precip parameter is 3-6 hr precip)
                if int(hr) % acc_hrs == 0 and int(hr) != 0:
                    # Open the dataset
                    filename = os.path.abspath(grib_dir + '/' + os.fsdecode(
                        inputArgs.date + '_' + inputArgs.modelRun + '_' + str(hr) + '.grib'))
                else:
                    continue

            gribObj.gribAll = gdal.Open(filename, GA_ReadOnly)
            # Notes: This next section is important because:
            #       1) We are quickly reducing the size of raster using the "projWin=" parameter.
            #       2) We are transforming the SNODAS grid from latlon (EPSG:4326) to XY (EPSG:3857)
            # We MUST transform this into XY coordinates because we are making calculations on the rasters based off of
            # meters (not decimal degrees). For example, once we use gdal.Warp, it will reproject it into xy coordinates, which
            # we can then use to find the x,y resolution of each grid box in meters. Now, we will know that for every raster
            # grid cell, we can calculate any type of Volume calculation we need, for example 100 mm of rainfall
            # over a 1,000 x 1,000 m grid will give us ~8 acre feet.

            # Clip the raster to a given bounding box (makes the raster much easier to work with)
            projInfo = grib.gribAll.GetProjection()
            try:
                spatialRef = osr.SpatialReference()  # Get spatial reference for converting this data into lon/lat (or x/y if it's already in lon/lat)
                spatialRef.ImportFromWkt(projInfo)

                cord_sys = spatialRef.GetAttrValue("UNIT|AUTHORITY", 1)
            except:
                print("NO COORDINATE SYSTEM FOUND! Transforming to XY")
                cord_sys = '4326'
            if cord_sys != '3875':
                # If we want tiles, we have to specify a bounding box. If we want a bounding box, the ONLY way to make
                # one is to put the coordinate system in EPSG4326
                if inputArgs.tiles4326 == True:
                    # Before we can apply a bounding box, we have to prject this into lat/lng coordinate system
                    # Note: GFS is already in a lat/lng coordinate system, so we don't need to do anything.
                    #if grib.model != 'gfs':
                        #gribObj.gribAll = gdal.Warp('/vsimem/temp.dat', gribObj.gribAll,dstSRS='EPSG:4326')
                    gribObj.gribAll = gdal.Translate('/vsimem/temp.dat', gribObj.gribAll, projWin=grib.bbox, projWinSRS='EPSG:4326')
                    gribObj.gribAll = gdal.Warp('/vsimem/temp.dat', gribObj.gribAll, dstSRS='EPSG:4326')

                    if gribObj.gribAll != None:
                        print("Successfully Reprojected Coordinate System into Lat/Lng Coords")
                        band_id, gribObj.units, gdate = find_band_number(gribObj, hr)
                        band = gribObj.gribAll.GetRasterBand(band_id)

                        data = BandReadAsArray(band)
                        data[data == 9999.0] = 0
                        grib_var4326.append(data)
                        grib_date.append(gdate)

                        #projInfo = grib.gribAll.GetProjection()  # Projection of grib file
                        #geoinformation = grib.gribAll.GetGeoTransform()  # [0] x min, [1] x resolution, [3] ymax, [5] yres
                        #spatialRef = osr.SpatialReference()  # Get spatial reference for converting this data into lon/lat (or x/y if it's already in lon/lat)
                        #spatialRef.ImportFromWkt(projInfo)
                        #spatialRefProj = spatialRef.ExportToProj4()
                        #xres = geoinformation[1]
                        #yres = geoinformation[5]
                        #xmin = geoinformation[0]
                        #xmax = geoinformation[0] + (xres * grib.gribAll.RasterXSize)
                        #ymin = geoinformation[3] + (yres * grib.gribAll.RasterYSize)
                        #ymax = geoinformation[3]
                        #print("Bounds \n (max lat / max lon = " + str(ymax) + "/" +
                        #      str(xmax) + "\n" + "min lat / min lon: " + str(ymin) +"/" + str(xmin))

                    else:
                        print("Failed to Reproject into Lat/Lng")


                # Note: Even if we put it into EPSG4326 for a bounding box, we still have to put it back into 3857 to make the png AND to calculate Acre Feet.
                gribObj.gribAll = gdal.Warp('/vsimem/temp.dat', gribObj.gribAll,
                                           dstSRS='EPSG:3857')  # If you wanted to put it into x/y coords

                print("Successfully Reprojected Coordinate System From Lat/Lng to X/Y")

            #hr = int(filename.rsplit('.grib', 1)[0][-2:])
            band_id, gribObj.units, gdate = find_band_number(gribObj, hr)
            band = gribObj.gribAll.GetRasterBand(band_id)

            data = BandReadAsArray(band)
            # data[data == 9999.0] = np.nan #if you resample (smooth) and there are ANY nan's in the array, it will ALL become nan's
            data[data == 9999.0] = 0
            grib_var.append(data)
            grib_date.append(gdate)
            downloadReset = True  # Grib file contained needed info, no download needed, reset to True.
        except Exception as e:
            print (e)
            logging.error("ERROR: File " + inputArgs.date+'_'+inputArgs.modelRun+'_'+str(hr)+'.grib' + " (with parameter " +
                  grib.variable +") Not Found in " + grib_dir)
            print("ERROR: File " + inputArgs.date+'_'+inputArgs.modelRun+'_'+str(hr)+'.grib' + " (with parameter " +
                  grib.variable +") Not Found in " + grib_dir)
            if downloadReset == True:
                i -= 1 #Can cause infinte loop if a reset switch is not used.
                downloadReset = False
                gd.main(grib, grib.modelRun)

    gribObj.data = np.array(grib_var)
    gribObj.data4326 = np.array(grib_var4326)
    gribObj.date = np.array(grib_date)
    #test = makePNG()
    return gribObj

def find_band_number(gribObj, hr):
    '''
    Finds the band number inside the GRIB file, given the variable and the level names
    '''
    for i in range(1,gribObj.gribAll.RasterCount + 1):
        band = gribObj.gribAll.GetRasterBand(i)
        metadata = band.GetMetadata()
        targetVar = gribObj.variable
        band_level = metadata['GRIB_SHORT_NAME']
        band_variable = metadata['GRIB_ELEMENT']
        time_stamp = int(metadata['GRIB_VALID_TIME'][:10])

            #time_stamp = int(metadata['GRIB_VALID_TIME'][:12])
            #gfs has strange naming conventions where things like precip will be named APCP03, APCP06.
            # This just says, if the variable we're looking for matches at least
            #some part of the band_variable, then it's the one we're looking for
        if gribObj.variable == 'APCP' and gribObj.variable in band_variable:
            targetVar = 'APCP03'
        grib_date=(datetime.fromtimestamp(time_stamp, pytz.timezone('US/Pacific')))
        print("band_level", band_level, "for: " + grib_date.strftime("%a %b %d %H:%M") + '\n', "band_variable", band_variable)
        if (targetVar == band_variable) and (gribObj.level == band_level):
            units = metadata['GRIB_UNIT']
            return i, units, grib_date


def model_hours():
    """
    The temporal resolution changes for each weather model. This function will return a list containing the valid
    hours for a given model.
    Args:
        None; however, grib.model must be defined.
    Returns:
        <int list> valid hours for a given model

    """
    if grib.model not in ['nam','namNest','gfs','hrrr']:
        logger.exception(e)
        raise ValueError("Check model name, " + grib.model + " is not a valid Model. " )

    maxhr = 2  # start with a value, which will be changed based on what model we're using
    if grib.model == 'nam':
        maxhr = 84
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))       # The +1 is since the first hr is 0 (so we need an array len of 85)
        firsthrs = maxhr[0:36]              # the nam has hourly data up to hour 36, so get 0-35
        secondhrs = maxhr[36::3]            # after hour 36, we need every third hour (36,39,42, etc)
        grib.hours = list(firsthrs) + list(
            secondhrs)                      # merge two lists together to get (0-36, then 36-84 every 3rd hr).
    if grib.model == 'namNest':
        maxhr = 60
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))       # The +1 is since the first hr is 0 (so we need an array len of 61)
        firsthrs = maxhr[0::1]              # the nam has hourly data up to hour 60, so increment up to hour 60
        grib.hours = list(firsthrs)         # merge two lists together to get (0-36, then 36-84 every 3rd hr).
    if grib.model == 'gfs':
        maxhr = 180
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))       # The +1 is since the first hr is 0 (so we need an array len of 85)
        firsthrs = maxhr[0:36]              # the nam has hourly data up to hour 36, so get 0-35
        secondhrs = maxhr[36::3]            # after hour 36, we need every third hour (36,39,42, etc)
        grib.hours = list(firsthrs) + list(
            secondhrs)                      # merge two lists together to get (0-36, then 36-84 every 3rd hr).
        # hrs = maxhr[0::3] #the gfs has hourly data up to hour 36, so get 0-35
    if grib.model == 'hrrr':
        maxhr = 18
        # The HRRR goes out to 36 hours for the following model runs: 00, 06, 12, and 18;
        if int(grib.modelRun) % 6 == 0:
            maxhr = 36
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))        # The +1 is since the first hr is 0 (so we need an array len of 15)
        grib.hours = maxhr[0::1]
    return grib.hours

def createMask(xxC,yyC,spatialRefProj):
    """
    This will create a mask for a given basin, which will be defined by grib.basin
    :param xxC: <np array> The center points of the x grid
    :param yyC: <np array> The center points in the y grid
    :param spatialRefProj: <gdal obj> Current Projection
    Note:
        Make sure the .geojson basin files exist. The files were created using the USGS streamstats online app:
        https://water.usgs.gov/osw/streamstats/
    :return: <np array> masked array.
    """
    # Get the lat / lon of each grid box's center point.
    lons_centerPt, lats_centerPt = Proj(spatialRefProj)(xxC, yyC, inverse=True)

    # given a bunch of lat/lons in the geotiff, we can get the polygon points.
    try:
        sf = fiona.open(os.path.join(dir_path,'Shapefiles',grib.basin,grib.basin+'.geojson'))
    except:
        logger.error('No basin shapefile found in ' +
                     dir_path + '/Shapefiles/' + grib.basin + '/' + grib.basin + '.geojson' +
                     '\n A mask was not created.')
        return
    geoms = [feature["geometry"] for feature in sf]
    poly = Polygon(geoms[0]['coordinates'][0])  # For a simple square, this would be 4 points, but it can be thousands of points for other objects.
    polyMinLon, polyMinLat, polyMaxLon, polyMaxLat = sf.bounds  # Get the bounds to speed up the process of generating a mask.

    # create a 1D numpy mask filled with zeros that is the exact same size as the lat/lon array from our projection.
    mask = np.zeros(lons_centerPt.flatten().shape)

    # Debugging: FOR CENTER POINT OF GRID BOX
    # These are test variables to see where the center points are (plot with basemaps to prove to yourself they're in the right spot).
    #debugCenterX = np.zeros(lats_centerPt.flatten().shape)
    #debugCenterY = np.zeros(lons_centerPt.flatten().shape)

    # Create Mask by checking whether points in the raster are within the bounds of the polygon. Instead of checking
    # every single point in the raster, just focus on points within the max/min bounds of the polygon (it's slow as hell
    # if you don't do that).
    i = 0  # counter
    for xp, yp in zip(lons_centerPt.flatten(), lats_centerPt.flatten()):
        if ((polyMinLon <= xp <= polyMaxLon) and (polyMinLat <= yp <= polyMaxLat)):
            mask[i] = (Point(xp, yp).within(poly))
            # Debugging FOR CENTER POINT OF GRID BOX: If you want to visualize the center point
            #       of each grid box that's found in the polygon,
            #       include this below and then you can put a dot (via m.plot)
            #if (Point(xp, yp).within(poly)):
                #debugCenterX[i] = xp
                #debugCenterY[i] = yp
        i += 1

    # mask = ([Point(xp, yp).within(poly) for xp, yp in zip(lons.flatten(), lats.flatten()) if
    #         ((polyMinLon <= xp <= polyMaxLat) and (polyMinLat <= yp <= polyMaxLat))])
    mask = np.reshape(mask, (xxC.shape))
    return mask

def calculateBasin(gribObj, xres, yres):
    """
    For a given basin, convert the total rainfall in each cell to acre feet and add them all up.
    :param gribObj:
    :param xres:
    :param yres:
    :return:
    """
    #maskArea = np.sum(mask) * xres * yres  # mask area in m^2 (total number of pixels * area of one pixel).
    basinTotal = 0
    if gribObj.units == '[kg/(m^2)]':
        value_mm = gribObj.data
        value_m = value_mm * 0.001  # mm to m.
        value_inches = value_mm * 0.03937  # mm to inches.
        #value_Cubed = totalPrecip_mm * abs(xres) * abs(yres)  # change each grid box to total precip in cubic meters
        value_AF = (value_m * abs(xres) * abs(yres))/ 1233.48  # 1 AF = 1233.48 m^3
        basinTotal = [np.sum(gribObj.mask * value_AF[i].T) for i in range(0, len(gribObj.data))]
    return np.sum(basinTotal)

def value_at_point(xx,yy, userLon, userLat, hr, spatialRefProj,m):
    # Calling a Proj class instance with the arguments lon, lat will convert lon/lat (in degrees) to x/y
    # native map projection coordinates (in meters).
    # If optional keyword 'inverse' is True (default is False), the inverse transformation
    # from x/y to lon/lat is performed.
    userX, userY = Proj(spatialRefProj)(userLon, userLat)

    # To find the value at a given lat/lon
    # 1) First, take the lat / lon and convert it to X/Y for the given projection.
    #    The reason you do this is because the numpy grid will contain a 2D grid for both the x values and y values
    #    but the x value at point x1,y1 will be same at the point x1,y2 (BUT THE LAT/LON would both be different)
    #    Since you're searching for a unique value (a specific grid box) you want to find the X value you're looking
    #    for and the Y value you're looking for. Once you have your X/Y grid box, you can find the value for that grid box.
    xIdx = (np.abs(xx - userX)).T.argmin()
    yIdx = (np.abs(yy - userY)).argmin()

    raster_val = grib.data[grib.hours.index(hr)].T[xIdx][yIdx]
    all_values = grib.data[0:].T[xIdx][yIdx]

    return raster_val

def makeMap(lons,lats,hr,m):

    if grib.acc == True:
        raster = sum(grib.data[0:grib.hours.index(hr)], axis=0) #cumulative
    else:
        raster = grib.data[grib.hours.index(hr)]  # 1 hr forecast (not cumulative)

    if grib.displayunits == 'US' and grib.units == '[kg/(m^2)]':
        raster = raster * 0.03937
        #grib.units = 'inches'

    if grib.displayunits == 'US' and grib.units == '[m]':
        raster = raster * 39.37

    # YOU CAN NOT HAVE NAN VAULES IN nparray BEFORE DOING scipy.ndimage.zoom
    #raster = scipy.ndimage.zoom(raster, 3)


    #YOU CAN NOT PUT NAN VAULES IN BEFORE DOING scipy.ndimage.zoom
    raster[raster < 0.01] = np.nan #this will prevent values of zero from being plotted.

    # For some rasters (like solar) they will be all zero at certain hours. Therefore, you can't calculate a
    # percentile off of a zero value and we will just set the max and min values for the color bar to arbitrary numbers
    minVal, maxVal = 0, 7
    if np.nanmax(raster) > 1:
        #maxVal = int(np.nanpercentile(raster, 99, interpolation='linear'))
        minVal = int(np.nanpercentile(raster, 1, interpolation='linear'))

    im = m.pcolormesh(lons, lats, raster.T, cmap=plt.cm.jet, vmin=minVal,
                      vmax=maxVal)  # return 50th percentile, e.g median., latlon=True)
    cb = m.colorbar(mappable=im, location='right', label='Total Precip (in)')

    #map = folium.Map([37, -120], zoom_start=3)

    #folium.plugins.ImageOverlay(
    #    image=raster,
    #    bounds=[[21.146938, -134.12207], [52.628703, -60.9284436]],  # ([SW lat, SW lng], [NE lat, NE lng])
    #    colormap = plt.cm.jet,
    #    opacity = 0.5,
    #).add_to(map)

    #map.add_child(folium.LatLngPopup())
    #map.save('index.html')

    subTitle = ''
    if isinstance(grib.basinTotal, ndarray):
        subTitle = "Total AF Change from today " + grib.date[grib.hours.index(hr)].strftime("%a %b %d %H:%M") + " : "  + \
                    str(int((grib.basinTotal[grib.hours.index(hr)]-grib.basinTotal[0])))

    #im1 = m.pcolormesh(lons, lats, mask, cmap=plt.cm.jet, vmin=0.1, vmax=4, latlon=True, alpha=0.5)
    #cb = plt.colorbar( orientation='vertical', fraction=0.05, shrink=0.7)

    ########################################################################################################################
    #Debugging: Test to prove a given lat/lng pair is accessing the correct precip grid box:
    #totalPrecip[1,1] = 4 #increase the variable by some arbitrary amount so it stands out.
    #xpts, ypts = m(lons_centerPt[1,1],lats_centerPt[1,1]) #This should be in the dead center of grid[1,1]
    #m.plot(xpts,ypts, 'bo')

    #Test to see a if the point at [x,y] is in the upper right corner of the cell (it better be!)
    #xpts, ypts = m(lons[1,1],lats[1,1]) #should be in upper right corner of cell
    #m.plot(xpts,ypts, 'bo')

    #Test to see the location of center points of each grid in polygon
    #debug_Xpoly_center_pts, debug_Ypoly_center_pts = m(debug_centerLonsInPoly,debug_centerLatsInPoly)
    #m.plot(debug_Xpoly_center_pts,debug_Ypoly_center_pts, 'bo')

    #Test to see all points
    #xtest, ytest = m(lons,lats)
    #m.plot(xtest,ytest, 'bo')
    ########################################################################################################################

    #plot shapefile
    try:
        m.readshapefile(dir_path + '/Shapefiles/' + grib.basin + '/' + grib.basin + '_EPSG4326', grib.basin + '_EPSG4326')
    except:
        logging.error("Basin Shapefile was not found in: " +
                      dir_path + '/Shapefiles/' + grib.basin + '/' + grib.basin + '_EPSG4326', grib.basin + '_EPSG4326')
        Warning("No Basin Shape File Found...Basins will remain absent from the map")

    output_dir = os.path.join(dir_path,'images',grib.variable)
    if grib.basin != None:
        output_dir = os.path.join(dir_path, 'images', grib.basin, inputArgs.date, grib.variable)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # annotate
    m.drawcountries()
    m.drawstates()
    #m.drawrivers()
    m.drawcounties(color='darkred')
    #m.drawcoastlines(linewidth=.5)
    #plt.title(grib.basin+' '+grib.variable+' '+grib.level+' ('+grib.units+') '+
    #          grib.date[grib.hours.index(hr)].strftime("%a %b %d %H:%M") + '\n' + subTitle)

    plt.title('NAM Model: Total Precipitation (in) through midnight Sunday')

    plt.savefig(os.path.join(output_dir,grib.model+str(hr)+'_'+grib.variable+'_'+grib.level+'.png'),dpi=775)
    print("File " + grib.model+str(hr)+'_'+grib.variable+'_'+grib.level+'.png Saved')
    plt.show()
    plt.close()



def makePlot():
    plt.title('Total AF within '+grib.basin+' basin')
    plt.bar(np.arange(len(grib.basinTotal)),grib.basinTotal, color='blue')
    label = (grib.basinTotal[x] for x in range(0,len(grib.basinTotal)))
    for i, v in enumerate(grib.basinTotal):
        plt.text(i, v + .5, str(int(v)))
    plt.xlabel('Forecast Hour')
    plt.ylabel('Total Acre Feet')
    plt.savefig('BAR_GRAPH_TOTAL_ACRE_FEET.png', dpi=775)
    plt.show()

def makePNG():
    for hr in grib.hours:
        base_folder = 'C:/xampp/htdocs/demos/build/tiles'
        output_dir = os.path.join(base_folder, inputArgs.date, grib.model, str(hr))

        #IF YOU USE GDAL.TRANSLATE with projWin, the projection MUST be in lat/lng (EPSG 4326)!
        #If you try to run GDAL.TRANSLATE in any other projection other than 4326 it will not work.

        #tileProjection = gdal.Warp('/vsimem/temp.dat', grib.gribAll,
        #                            dstSRS='EPSG:4326')  # If you wanted to put it into x/y coord
        #tileProjection = gdal.Translate('temp.dat', tileProjection, projWin=grib.bbox, bandList=[band_id])
        #tileProjection = gdal.Warp('/vsimem/temp.dat', tileProjection,
        #                          dstSRS='EPSG:3857')  # If you wanted to put it into x/y coord

        if grib.acc == True:
            #gfs only does accumulation every 3 hours.
            if grib.model == 'gfs':
                raster = sum(grib.data4326[0:grib.hours.index(hr/3)], axis=0)  # cumulative
            else:
                raster = sum(grib.data4326[0:grib.hours.index(hr)], axis=0) #cumulative
        else:
            raster = grib.data4326[grib.hours.index(hr)]  # 1 hr forecast (not cumulative)

        if grib.displayunits == 'US' and grib.units == '[kg/(m^2)]':
            raster = raster * 0.03937
            grib.units = 'inches'
            raster = raster * 20 # To get a nicer color profile

        if grib.displayunits == 'US' and grib.units == '[m]':
            raster = raster * 39.37

        if grib.displayunits == 'US' and grib.units == '[C]':
            raster = ((raster * 9 / 5) + 32 + 100)
            raster[raster==132] = np.nan

        if grib.displayunits == 'US' and grib.units == '[dB]':
            raster = raster * 4

        # YOU CAN NOT HAVE NAN VAULES IN nparray BEFORE DOING scipy.ndimage.zoom
        #raster = scipy.ndimage.zoom(raster, 3)


        #YOU CAN NOT PUT NAN VAULES IN BEFORE DOING scipy.ndimage.zoom
        raster[raster < 0.01] = np.nan #this will prevent values of zero from being plotted.

        im = Image.fromarray((raster).astype('uint8'))
        #im = Image.fromarray(((raster * 9 / 5) + 32 + 150).astype('uint8'))
        # basewidth = 100
        # wpercent = (basewidth / float(im.size[0]))
        # hsize = int((float(im.size[1]) * float(wpercent)))
        # im = im.resize((basewidth, hsize), Image.ANTIALIAS)
        tile_dir = os.path.join(output_dir, '5','5')
        if not os.path.exists(tile_dir):
            os.makedirs(tile_dir)
        im.save(tile_dir + '/radar.png', "PNG")
    return


if __name__ == "__main__":
    main()