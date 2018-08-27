'''
Version 1.3
Created on December 4, 2014
@author: Shane Motley

@Purpose: Downloads portion of grib data files from the NOMADS server and places
         into the folder YYYYMMDD/Model/

@Usage: User can scroll down and change the downloaded parameters by uncommenting
        various fields. User can also change the model to download (e.g. nam, gfs).

@VersionHistory:
    v1.0:  12/4/2014: Added a model parameter at main() --> main(model) to
           either take a value of None if the program is being run
           directly, or to take a user's model if the program is being
           called on by another python program such as makeGribs.py
           This makes gribDownload.py an independent program if needed.
           v1.1:  6/26/2015: Converted to python 3
                  v1.2:  7/31/2015: Changed folder structure for EC2
                  v1.3:  8/24/2015: Added support for GFS and HRRR

'''
import urllib.request, urllib.parse, urllib.error
import urllib.parse
import dateutil.parser
from datetime import datetime, timedelta
import tarfile
import gzip
import requests
import os
import time


def main(gribObj, hour):
    ######################################################################################
    ############################USER ADJUSTABLE VARIABLES#################################
    forecast_date = time.strftime("%Y%m%d")  # This will be today's date in yyyymmdd format
    #forecast_date = '20180222'
    forecast_time = gribObj.modelRun  # What time the forecast is (00, 06, 12, 18)
    if gribObj == None:
        model = 'namNest'
    else:
        model = gribObj.model
    if hour == None:
        forecast_time = '06'  # What time the forecast is (00, 06, 12, 18)
    else:
        forecast_time = hour
    ######################################################################################
    ######################################################################################

    baseFolder = os.path.join(os.path.dirname(os.path.realpath(__file__)),'grib_files')
    output_dir = os.path.join(baseFolder,forecast_date,model)
    maxhr = 2  # start with a value, which will be changed based on what model we're using
    if model == 'nam':
        maxhr = 84
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))  # The +1 is since the first hr is 0 (so we need an array len of 85)
        firsthrs = maxhr[0:36]  # the nam has hourly data up to hour 36, so get 0-35
        secondhrs = maxhr[36::3]  # after hour 36, we need every third hour (36,39,42, etc)
        hrs = list(firsthrs) + list(secondhrs)  # merge two lists together to get (0-36, then 36-84 every 3rd hr).
    if model == 'namNest':
        maxhr = 60
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))  # The +1 is since the first hr is 0 (so we need an array len of 85)
        firsthrs = maxhr[0::1]  # the nam has hourly data up to hour 60, so get 0-59
        hrs = list(firsthrs)  # merge two lists together to get (0-36, then 36-84 every 3rd hr).
    if model == 'gfs':
        maxhr = 180
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))  # The +1 is since the first hr is 0 (so we need an array len of 85)
        firsthrs = maxhr[0:36]  # the nam has hourly data up to hour 36, so get 0-35
        secondhrs = maxhr[36::3]  # after hour 36, we need every third hour (36,39,42, etc)
        hrs = list(firsthrs) + list(secondhrs)  # merge two lists together to get (0-36, then 36-84 every 3rd hr).
        # hrs = maxhr[0::3] #the gfs has hourly data up to hour 36, so get 0-35
    if model == 'hrrr':
        maxhr = 16
        # Get a range of all the hours we'll need to loop through
        maxhr = range(0, (maxhr + 1))  # The +1 is since the first hr is 0 (so we need an array len of 15)
        hrs = maxhr[0::1]

    # USE THE FIELDS BELOW FOR FULL DOMAIN
    top_lat=90 #Top of bounding box (North)
    bottom_lat=-90 #Bottom of bounding box (South)
    left_lon=0 #Left of bounding box (West)
    right_lon=360 #Right of bounding box (East)
    if gribObj != None:
        left_lon,top_lat,right_lon,bottom_lat = gribObj.bbox

    for hr in hrs:
        if model == 'nam':
            if hr < 10:
                hr = '0' + str(hr)
            griburl = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_nam.pl?'
            griburl = griburl + 'file=nam.t' + str(forecast_time) + 'z.awphys' + str(hr) + '.tm00.grib2'
        if model == 'namNest':
            if hr < 10:
                hr = '0' + str(hr)
            griburl = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl?'
            griburl = griburl + 'file=nam.t' + str(forecast_time) + 'z.conusnest.hiresf' + str(hr) + '.tm00.grib2'
        if model == 'gfs':
            if hr < 100:
                hr = '0' + str(hr)
            if int(hr) < 10:
                # Note: If hr < 10, it's also less than 100, so this will add 00 to the hr.
                hr = '0' + str(hr)
            griburl = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?'
            griburl = griburl + 'file=gfs.t' + str(forecast_time) + 'z.pgrb2.0p25.f'
            griburl = griburl + str(hr)
        if model == 'hrrr':
            if hr < 10:
                hr = '0' + str(hr)
            griburl = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?'
            griburl = griburl + 'file=hrrr.t' + str(forecast_time) + 'z.wrfsfcf'
            griburl = griburl + str(hr) + '.grib2'

        # Select atmospheric levels
        griburl = griburl + '&lev_surface=on'  # surface mb level
        #griburl=griburl+'&lev_2_m_above_ground=on'  #2 m above ground
        #griburl=griburl+'&lev_10_m_above_ground=on'  #10 m above ground
        #griburl = griburl + '&lev_80_m_above_ground=on'  # 80 m above ground
        #griburl = griburl + '&lev_1000_m_above_ground=on'  # 80 m above ground
        #griburl = griburl + '&lev_30-0_mb_above_ground=on'  # 0 to 30mb above ground (approx 0 to 280 m)

        # griburl=griburl+'&lev_1000_mb=on'  #1000 mb level
        # griburl=griburl+'&lev_975_mb=on'   #975 mb level
        # griburl=griburl+'&lev_950_mb=on'   #950 mb level
        # griburl=griburl+'&lev_925_mb=on'   #925 mb level
        # griburl=griburl+'&lev_900_mb=on'   #900 mb level
        # griburl=griburl+'&lev_850_mb=on'   #850 mb level
        # griburl=griburl+'&lev_800_mb=on'   #800 mb level
        # griburl=griburl+'&lev_750_mb=on'   #750 mb level
        # griburl=griburl+'&lev_700_mb=on'   #700 mb level
        # griburl=griburl+'&lev_650_mb=on'   #650 mb level
        # griburl=griburl+'&lev_600_mb=on'   #600 mb level
        # griburl=griburl+'&lev_550_mb=on'   #550 mb level
        # griburl=griburl+'&lev_500_mb=on'   #500 mb level
        # griburl=griburl+'&lev_450_mb=on'   #450 mb level
        # griburl=griburl+'&lev_400_mb=on'   #400 mb level
        # griburl=griburl+'&lev_350_mb=on'   #350 mb level
        # griburl=griburl+'&lev_300_mb=on'   #300 mb level
        # griburl=griburl+'&lev_250_mb=on'   #250 mb level
        # griburl=griburl+'&lev_200_mb=on'   #200 mb level
        # griburl=griburl+'&lev_150_mb=on'   #150 mb level
        # griburl=griburl+'&lev_100_mb=on'   #100 mb level
        # griburl=griburl+'&lev_70_mb=on'    #70 mb level
        # griburl=griburl+'&lev_30_mb=on'    #30 mb level
        # griburl=griburl+'&lev_20_mb=on'    #20 mb level
        # griburl=griburl+'&lev_10_mb=on'    #10 mb level

        #These variables will be added to the download (passed from GribToRaster.py)
        if gribObj != None:
            if gribObj.level == '0-SFC':
                griburl = griburl + '&lev_surface=on'  # surface mb level
            if gribObj.level == '1000-HTGL':
                griburl = griburl + '&lev_1000_m_above_ground=on'  # 1000 m above ground
            if gribObj.level == '2-HTGL':
                griburl = griburl + '&lev_2_m_above_ground=on'  # 2 m above ground
            if gribObj.level == '0-CTL':
                griburl = griburl + '&lev_cloud_top=on'  # Cloud Top
            if gribObj.level == '0-NTAT':
                griburl = griburl + '&lev_top_of_atmosphere=on'  # 2 m above ground

        # Select variables

        # griburl=griburl+'&var_HGT=on'  #Height (geopotential m)
        #griburl = griburl + '&var_RH=on'  # Relative humidity (%)
        griburl = griburl + '&var_TMP=on'  # Temperature (K)
        #griburl = griburl + '&var_UGRD=on'  # East-West component of wind (m/s)
        #griburl = griburl + '&var_VGRD=on'  # North-South component of wind (m/s)
        #griburl = griburl + '&var_DSWRF=on'  # Downward Shortwave Radiation Flux (W/m^2)
        #griburl=griburl+'&var_VVEL=on' #Vertical Windspeed (Pa/s)
        #griburl = griburl + '&var_WATR=on'  # GFS Water Runoff(kg/m^3 which is mm)
        griburl = griburl + '&var_APCP=on'  # Total Precipitation (kg/m^3 which is mm)
        #griburl = griburl + '&var_SSRUN=on'  # Storm Surface Runoff (kg/m^3 which is mm)
        #griburl = griburl + '&var_SNOM=on'  # Snow melt (kg/m^3 which is mm)
        #griburl = griburl + '&var_REFD=on'  # Reflectivity (dB)
        griburl = griburl + '&var_WEASD=on'  # Water Equivalent of Accumulated Snow Depth (kg/m^3 which is mm)
        griburl = griburl + '&var_SNOD=on'  # Snow Depth (kg/m^3 which is mm)

        if gribObj != None:
            griburl = griburl + '&var_' + gribObj.variable + '=on'

        # Select bounding box
        #griburl = griburl + '&subregion=' #TURN THIS LINE ON IF YOU ARE ACTUALLY USING A SUBREGION!!
        griburl = griburl + '&leftlon=' + str(left_lon) #SEE LINE ABOVE
        griburl = griburl + '&rightlon=' + str(right_lon)
        griburl = griburl + '&toplat=' + str(top_lat)
        griburl = griburl + '&bottomlat=' + str(bottom_lat)

        # Select date and time
        if model == 'nam':
            griburl = griburl + '&dir=%2Fnam.' + forecast_date

        if model == 'namNest':
            griburl = griburl + '&dir=%2Fnam.' + forecast_date

        if model == 'gfs':
            griburl = griburl + '&dir=%2Fgfs.' + forecast_date + forecast_time

        # WE HAVE TO TREAT THE HRRR AS A SPECIAL CASE AND OVERWRITE THE EXISTING URL
        # CREATED ABVOE, SO WE START FROM SCRATCH.
        if model == 'hrrr':
            griburl = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?'
            griburl = griburl + 'file=hrrr.t' + str(forecast_time) + 'z.wrfsfcf'
            griburl = griburl + str(hr) + '.grib2'
            griburl = griburl + '&lev_10_m_above_ground=on'  # 10 m above ground
            griburl = griburl + '&lev_80_m_above_ground=on'  # 80 m above ground
            griburl = griburl + '&var_UGRD=on'  # East-West component of wind (m/s)
            griburl = griburl + '&var_VGRD=on'  # North-South component of wind (m/s)
            # Select bounding box
            griburl = griburl + '&subregion='
            griburl = griburl + '&leftlon=' + str(left_lon)
            griburl = griburl + '&rightlon=' + str(right_lon)
            griburl = griburl + '&toplat=' + str(top_lat)
            griburl = griburl + '&bottomlat=' + str(bottom_lat)
            griburl = griburl + '&dir=%2Fhrrr.' + forecast_date

        print(griburl)
        print('Downloading GRIB file for date ' + forecast_date + ' time ' + forecast_time + ', forecasting ' + str(
            hr) + ' hours ahead...')
        # URL Constructed, get parameters.
        connect_timeout = 15.0  # Connection will timeout after xx seconds
        read_timeout = 15.0  # Download will time out after xx seconds
        # Use a `Session` instance to customize how `requests` handles making HTTP requests.
        session = requests.Session()
        session.mount("http://", requests.adapters.HTTPAdapter(max_retries=10))
        session.mount("https://", requests.adapters.HTTPAdapter(max_retries=102))
        try:
            webfile = session.post(griburl, timeout=(connect_timeout, read_timeout))
        except requests.exceptions.HTTPError as e:
            print ("You got an HTTPError:"), e.message
        # webfile=urllib.request.urlopen(griburl,'POST'.encode(encoding='utf-8'))
        print("Download complete.  Saving...")

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        local_filename = forecast_date + '_' + forecast_time + '_' + str(hr) + '.grib'
        local_filename = os.path.join(output_dir, local_filename)
        try:
            local_file = open(local_filename,
                              'wb')  # 'wb' is the same as 'w' in UNIX, Linux, MacOS, but it is the safest notation.
        except IOError:
            print ("Can Not Open", local_file)
        else:
            local_file.write(webfile.content)
            print('Requested grib data written to file: ' + local_filename)
            local_file.close()
            webfile.close()

# If you run the program from here (i.e. it's not called from
# another program), then the name is "__main__" and we will
# pass a model name of "None"
if __name__ == "__main__":
    main(None, None)
