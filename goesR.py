from osgeo import osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import os, sys
from datetime import datetime, timedelta
import pytz
import netCDF4 as nc
from subprocess import call
from mapbox import Uploader
import uuid
import botocore
import boto3
import re

gdaladdoFile = 'C:/Users/smotley/AppData/Local/Continuum/anaconda3/Library/bin/gdaladdo.exe'
mb_access_token = os.environ['MAPBOX_TOKEN']
uploader = Uploader(access_token=mb_access_token)


def main():
    hours_of_data = 1
    last_file = False #Do you want just the latest file, or all files since
    bbox = [-134.12207, 52.628703, -60.9284436, 21.146938]  # namNest Bounds
    gribData(bbox, last_file, hours_of_data)


def gribData(bbox, grab_last_file, hours_of_data):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sat_dir = os.path.join(dir_path,'satellite')
    #product = 'ABI-L1b-RadF' #Full Disk
    product = 'ABI-L1b-RadC' #CONUS
    band_num = 2 #ABI Band 1 can effectively be thought of as the blue channel of visible light seen by the satellite. Used on its own it can be effective in identifying aerosols (dust, haze, smoke) suspended in the atmosphere, in contrast to clouds and ground features. These can be identified by locating areas of moderate to low brightness with a milky or featureless appearance on ABI band 1 versus the 'Red' band (ABI 2) where these features effectively disappear. This is because of the difference in scattering of red light verses blue light.
                 #ABI Band 2 can effectively be thought of as the red channel of visible light seen by the satellite. Used on its own it is effective in identifying clouds in contrast to ground and large bodies of water. At 0.64 μm, it is effectively the same as the traditional 'Visible' channel we are familiar with from GOES-13 but with better optics onboard GOES-16 (GOES East), this product is now transmitted at a much higher spatial resolution. Additionally, because of improved transmission equipment and strategies, the temporal resolution is also improved. Briefly, this means the data now comes in more frequently with higher clarity.

    now = datetime.utcnow()
    past = now - timedelta(hours=hours_of_data)

    day_of_year = now.timetuple().tm_yday
    yyyy = now.year
    hh = now.hour

    day_of_year2 = past.timetuple().tm_yday
    yyyy2 = past.year
    hh2 = past.hour

    s3 = boto3.resource('s3')
    bucket = s3.Bucket('noaa-goes16')

    # get keys for restricted files
    prefix = '{}/{:04d}/{:d}/{:d}/{}{:02d}{}'.format(product,yyyy,day_of_year,hh,'OR_'+product+'-M3C',band_num,'_G16_')
    prefix2 = '{}/{:04d}/{:d}/{:d}/{}{:02d}{}'.format(product, yyyy2, day_of_year2, hh2, 'OR_' + product + '-M3C', band_num,
                                                     '_G16_')
    objs = bucket.objects.filter(Prefix=prefix)
    objs2 = bucket.objects.filter(Prefix=prefix2)

    file_paths = [o.key for o in objs if o.last_modified > (datetime.utcnow() - timedelta(hours=hours_of_data)).replace(tzinfo=pytz.UTC)]
    file_paths2 = [o.key for o in objs2 if o.last_modified > (datetime.utcnow() - timedelta(hours=hours_of_data)).replace(tzinfo=pytz.UTC)]
    file_paths = file_paths+file_paths2
    file_paths.sort()

    if grab_last_file == True:
        # Must place file_paths in [] since we only have one file and we need to pass a list.
        upload_latest_files([file_paths[-1]], s3, bbox)
    else:
        upload_latest_files(file_paths, s3, bbox)

    print("PROCESS COMPLETE")
    return


def upload_latest_files(latest_files, s3, bbox):
    visHr = 0
    for file in reversed(latest_files):
        myDownloadDir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'satellite','band2')

        outfile_name = file.split('/')[-1] #Save the file as ##+fileName where hh is the number of the
        file_creation_date = (re.search("_c(.*).nc", outfile_name)).group(1) #Find the string between _c and .nc -> doyHHMM
        try:
            #Check to see if file is already downloaded. If not, download the file.
            if os.path.isfile(os.path.join(myDownloadDir,outfile_name)) == False:
                print("Downloading " + file)
                s3.Bucket('noaa-goes16').download_file(file, os.path.join(myDownloadDir,outfile_name))
                print("Saved " + file + " as " + outfile_name)
                filename = os.path.join(myDownloadDir,outfile_name)
                tilename = os.path.join(myDownloadDir,file_creation_date+'.mbtiles')
                uploadname = 'vis'+file_creation_date #We are putting 'vis' in the name so that we can search for 'vis' and delete these tiles when needed.

                print("Opening " + filename + " to Warp file into EPSG: 4326")
                g16nc = gdal.Open('NETCDF:"' + filename + '":Rad', GA_ReadOnly)
                ncAll = gdal.Warp('/vsimem/temp.dat', g16nc, dstNodata=-999.0,
                                  dstSRS='EPSG:4326')
                # ncAll = gdal.Translate('/vsimem/temps.dat', ncAll, projWin=bbox) #Note, the virtual memory file must have a dif name.
                band = ncAll.GetRasterBand(1)
                #data = BandReadAsArray(band)

                print("Creating tileset: " + tilename)
                gdal.Translate(tilename, ncAll, format='mbtiles', projWin=bbox, scaleParams=[[0, 3000, 0, 255]])

                print("Calling gdaladdo: " + tilename)
                call([gdaladdoFile, '-r', 'nearest', '-b', '1', tilename, '2 4 8'])
                visHr += 1
                with open(tilename, 'rb') as src:
                    print("UPLOADING TILE "+tilename)
                    upload_resp = uploader.upload(src, uploadname)
                    # print(uploader.status('smotley.OR2').json())
                    print("DONE UPLOADING TILE "+tilename)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                print("The file " + file + " does not exist.")
            else:
                raise
    cleanUp(True)

def cleanUp(delete_old_tiles):
    curDir = os.path.dirname(os.path.realpath(__file__))
    satDir = os.path.join(curDir,'satellite','band2')
    for filename in os.listdir(satDir):
        file_date = os.path.getmtime(os.path.join(satDir,filename))
        if file_date < (datetime.utcnow() - timedelta(hours=24)).timestamp():
           os.remove(os.path.join(satDir,filename))
    if delete_old_tiles == True:
        tiles_on_server = uploader.list().json()
        for tile in tiles_on_server:
            cdate = datetime.strptime(tile['created'], "%Y-%m-%dT%H:%M:%S.%fZ")
            tname = tile['name']
            if 'vis' in tname and cdate < datetime.utcnow() - timedelta(hours=5):
                uploader.delete(tname)




if __name__ == "__main__":
    main()