'''
Version 1.3
Created on December 4, 2014
@author: Shane Motley

@Purpose: Downloads portion of grib data files from the NOMADS server and places
         into the folder YYYYMMDD/Model/

@Usage: User can scroll down and change the downloaded parameters by uncommenting
        various fields. User can also change the model to download (e.g. nam, gfs).

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


def main(model, analysis_date):
    ######################################################################################
    ############################USER ADJUSTABLE VARIABLES#################################
    if analysis_date == None:
        analysis_date = time.strftime("%Y%m%d")  # This will be today's date in yyyymmdd format
    #analysis_date = '20100401'
    code = 1034  # 1034 is SWE, 1036 is snow depth
    if model == None:
        model = 'SNODAS'
    ######################################################################################
    ######################################################################################

    outdir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'grib_files')
    dt = analysis_date
    snodas_url_str = None
    snodas_outdir = os.path.join(outdir, dt, 'snodas')
    if not os.path.exists(snodas_outdir):
        os.makedirs(snodas_outdir)

    dt = dateutil.parser.parse(str(dt), fuzzy=True)
    # Note: unmasked products (beyond CONUS) are only available from 2010-present
    if dt >= datetime(2003, 9, 30) and dt < datetime(2010, 1, 1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/%Y/%m_%b/SNODAS_%Y%m%d.tar'
        tar_subfn_str_fmt = 'us_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    elif dt >= datetime(2010, 1, 1):
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/%Y/%m_%b/SNODAS_unmasked_%Y%m%d.tar'
        tar_subfn_str_fmt = './zz_ssmv1%itS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    else:
        print("No SNODAS data available for input date")

    if snodas_url_str is not None:
        snodas_url = dt.strftime(snodas_url_str)

        ################DOWNLOAD FILE##############
        fn = os.path.split(snodas_url)[-1]
        if snodas_outdir is not None:
            fn = os.path.join(snodas_outdir, fn)
        if not os.path.exists(fn):
            try:
                from urllib.request import urlretrieve
            except ImportError:
                from urllib import urlretrieve
            print("Retrieving: %s" % snodas_url)
            # Add progress bar
            urlretrieve(snodas_url, fn)
        snodas_tar_fn=fn
        ############DOWNLOAD COMPLETE#############

        print("Unpacking SNODAS dataset")
        tar = tarfile.open(snodas_tar_fn)
        # gunzip to extract both dat and Hdr files, tar.gz
        for ext in ('dat', 'Hdr'):
            tar_subfn_str = tar_subfn_str_fmt % (code, ext)
            tar_subfn_gz = dt.strftime(tar_subfn_str)
            tar_subfn = os.path.splitext(tar_subfn_gz)[0]
            print(tar_subfn)
            if snodas_outdir is not None:
                tar_subfn = os.path.join(snodas_outdir, tar_subfn)
            if not os.path.exists(tar_subfn):
                # Should be able to do this without writing intermediate gz to disk
                tar.extract(tar_subfn_gz)
                with gzip.open(tar_subfn_gz, 'rb') as f:
                    outf = open(tar_subfn, 'wb')
                    outf.write(f.read())
                    outf.close()
                os.remove(tar_subfn_gz)

        # Need to delete 'Created by module comment' and 'Last modified by module comment' line from Hdr, can contain too many characters
        snodas_fn = tar_subfn
        f = open(snodas_fn)
        output = []
        for line in f:
            if len(line) < 200:
                output.append(line)
        f.close()
        f = open(snodas_fn, 'w')
        f.writelines(output)
        f.close()
        print("Done Saving SNODAS File...")

# If you run the program from here (i.e. it's not called from
# another program), then the name is "__main__" and we will
# pass a model name of "None"
if __name__ == "__main__":
    main(None, None)
