#!/usr/bin/env python

# gdalogr_catalogue.py
# Purpose: Catalog all raster & vector datasources/layers found in a directory tree
# Usage: python gdalogr_catalogue.py <search path>
# Sends hierarchical JSON to stdout
# Requires GDAL/OGR libraries: http://pypi.python.org/pypi/GDAL

# Author: Tyler Mitchell, Matthew Perry Jan-2008

# CHANGELOG
# JAN-06 - Initial release of ogr_catalog5.py
# 27-DEC-07 - fixed try statement so it fails more gracefully
# 5-JAN-08 - converted to using XML output, pipe to a text file
# 7-JAN-08 - started refactoring, raster output complete
# 8-JAN-08 - added rudimentary vector output
# 16-JAN-08 - renamed outputs entities, added higher level elements/summary stats
# 12-FEB-08 - Added bad hack to output INSERT statements if you add 2nd argument in command "SQL".  e.g. python gdalogr_catalogue.py ../ SQL | grep INSERT - hacked for Markus :)
# 25-MAR-12 - Rework arg handling a bit, made indentation consistent, 
# 09-JAN-13 - use a pure-python data structure and convert directly to json 

'''
TODO
# DONE - higher level attributes about process: num of files, dirs search, timestamp
# DONE - filesize, user/owner, moddate/timestamp for entries
- extent values to GML or basic WKT bbox
# DONE - decide on checksum process for determining changes
- decide on process -> datasource linking (timestamp?) for top level relations
'''
import logging
import os, sys
import hashlib
import json
from optparse import OptionParser, OptionGroup
from string import strip
from time import asctime
from datetime import datetime
from pyproj import Proj
from osgeo import gdal
from osgeo import osr
from osgeo import ogr

gdal.UseExceptions()
ogr.UseExceptions()

logging.basicConfig( stream=sys.stderr, level=logging.WARNING )
log = logging.getLogger("catalog")

def spinning_cursor():
    cursor='/-\|'
    i = 0
    while 1:
        yield cursor[i]
        i = (i + 1) % len(cursor)

def getCatalog(startpath):
    skiplist = ['.svn','.shx','.dbf', '.prj', '.aux.xml', '.e00', '.adf']
    gdal.PushErrorHandler()
    pathwalker = os.walk(startpath)
    counterraster = 0
    countervds = 0
    cursor = spinning_cursor()

    dirlist, filelist = [], []

    catalog = { 
        'raster_data': [],
        'vector_data': [],
    }

    for eachpath in pathwalker:
        startdir = eachpath[0]

        alldirs = eachpath[1]
        dirlist += alldirs

        allfiles = eachpath[2]
        filelist += allfiles

        for eachfile in allfiles + alldirs:
            currentfile = os.path.join(startdir, eachfile)
            raster, vector = None, None
            starttime = datetime.now()
            if (not skipfile(currentfile,skiplist)):
                raster, vector = tryopends(currentfile)
            if raster:
                try:
                    raster = processraster(raster, counterraster, currentfile)
                    raster['time'] = (datetime.now() - starttime).total_seconds()
                    catalog['raster_data'].append(raster)
                    counterraster += 1
                except NotGeographic:
                    pass
            if vector:
                if skipfile(vector.GetName(), skiplist):
                    raise Exception("This should not happen")
                try:
                    vector = processvds(vector, countervds, currentfile)
                    vector['time'] = (datetime.now() - starttime).total_seconds()
                    catalog['vector_data'].append(vector)
                    countervds += 1
                except NotGeographic:
                    pass

        sys.stderr.write("\rFound %d vector and %d rasters datasets.   %s" % (countervds, 
            counterraster, cursor.next()) )
        sys.stderr.flush()
    
    sys.stderr.write("\n")
    catalog['meta'] = processmeta(startpath, skiplist, dirlist, filelist)
    return catalog

class NotGeographic(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        log.info(message)

def skipfile(filepath, skiplist):
    skipstatus = None
    for skipitem in skiplist:
        if filepath.find(skipitem) > 0: 
            skipstatus = True
            return True
        else:
            skipstatus = False
    return skipstatus
  
def tryopends(filepath):
    dsogr, dsgdal = False, False
    # GetFeatureCount and Extent are Way too slow and called many times TODO
    skip_drivers = ['AVCBin',] 

    try:
        dsgdal = gdal.OpenShared(filepath)
    except RuntimeError:
        pass
    except gdal.GDALError:
        pass

    try:
        dsogr = ogr.OpenShared(filepath)
    except RuntimeError:
        pass
    except ogr.OGRError:
        pass

    if dsogr and dsogr.GetDriver().GetName() in skip_drivers:
        dsogr = False

    return dsgdal, dsogr

def extentToLatLon(extent, proj):
    if proj == '' or proj is None:
        return None
    p1 = Proj(proj, preserve_units=True) 
    if p1.is_latlong():
        return extent
    x1,y1 = p1(extent[0], extent[1], inverse=True)
    x2,y2 = p1(extent[2], extent[3], inverse=True)
    return x1, y1, x2, y2

def processraster(raster, counterraster, currentpath):
    rastername = raster.GetDescription()
    log.debug(rastername)
    bandcount = raster.RasterCount
    geotrans = strip(str(raster.GetGeoTransform()),"()")
    geotrans = [float(strip(x)) for x in geotrans.split(",")]
    if geotrans == [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]:
        raise NotGeographic("No geotrans found; not a geographic dataset")
    driver = raster.GetDriver().LongName
    rasterx = raster.RasterXSize
    rastery = raster.RasterYSize
    wkt = raster.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    proj = srs.ExportToProj4()

    # llx, lly, urx, ury
    extent = (geotrans[0], 
              geotrans[3] + (geotrans[5] * rastery), 
              geotrans[0] + (geotrans[1] * rasterx), 
              geotrans[3])
    latlon_extent = extentToLatLon(extent, proj)

    resultsbands = []
    for bandnum in range(1, bandcount + 1): # 1-indexed
        band = raster.GetRasterBand(bandnum)
        minval, maxval = band.ComputeRasterMinMax(1) # approx_ok=1
        overviews = band.GetOverviewCount()
        resultseachband = {
                'bandId': bandnum, 
                'min': minval,
                'max': maxval, 
                'overviews': overviews
        }
        resultsbands.append(resultseachband)

    resultsraster = { 
            'bands': resultsbands, 
            'rasterId': counterraster, 
            'name': rastername, 
            'bandcount': bandcount, 
            'geotrans': geotrans, 
            'driver': driver, 
            'rasterX': rasterx, 
            'rasterY': rastery, 
            'projection': wkt,
            'proj': proj,
            'latlonextent': latlon_extent,
            'extent': extent,
            'filestats': getFileStats(currentpath)
    }

    return resultsraster
  
def processmeta(startpath, skiplist, dirlist, filelist):
    processValues = {
            'SearchPath':startpath,
            'LaunchPath':os.getcwd(),
            'UserHome':os.getenv("HOME"),
            'IgnoredString':" ".join(map(str, skiplist)),
            'DirCount':int(len(dirlist)),
            'FileCount':int(len(filelist)),
            'Timestamp':asctime()
    }
    return processValues

def processvds(vector, countervds, currentpath):
    vdsname = vector.GetName()
    log.debug(vdsname)
    vdsformat = vector.GetDriver().GetName()
    vdslayercount = vector.GetLayerCount()
    filestats = getFileStats(currentpath)
    if filestats['fileType'] == 'Directory' and vdsformat == 'ESRI Shapefile':
        raise NotGeographic("Just a directory of shapefiles; nothing to see here")

    resultslayers = []
    for layernum in range(vdslayercount): #process all layers
        layer = vector.GetLayer(layernum)
        spatialref = layer.GetSpatialRef()
        if spatialref:
            layerproj = spatialref.ExportToProj4()
        else:
            layerproj = None
        layername = layer.GetName()
        layerfcount = str(layer.GetFeatureCount())
        layerextentraw = strip(str(layer.GetExtent()),"()")
        le = [float(x) for x in layerextentraw.split(',')]
        # reorder to llx, lly, urx, ury
        layerextent = (le[0], le[2], le[1], le[3])
        layerftype = featureTypeName(layer.GetLayerDefn().GetGeomType())

        if layerproj:
            latlon_extent = extentToLatLon(layerextent, layerproj)
        else:
            latlon_extent = None

        resultseachlayer = {
            'layerId': layernum, 
            'name': layername, 
            'proj': layerproj, 
            'featuretype': layerftype, 
            'featurecount': layerfcount, 
            'extent': layerextent,
            'latlonextent': latlon_extent
        }
        resultslayers.append(resultseachlayer)

    resultsvector = { 
            'datasourceId': str(countervds), 
            'name': vdsname, 
            'format': vdsformat, 
            'layercount': str(vdslayercount), 
            'filestats': filestats,
            'layers': resultslayers
    }

    return resultsvector

def featureTypeName( type ):
    if type == ogr.wkbUnknown:
        return 'Unknown'
    elif type == ogr.wkbPoint:
        return 'Point'
    elif type == ogr.wkbLineString:
        return 'LineString'
    elif type == ogr.wkbPolygon:
        return 'Polygon'
    elif type == ogr.wkbMultiPoint:
        return 'MultiPoint'
    elif type == ogr.wkbMultiLineString:
        return 'MultiLineString'
    elif type == ogr.wkbMultiPolygon:
        return 'MultiPolygon'
    elif type == ogr.wkbGeometryCollection:
        return 'GeometryCollection'
    elif type == ogr.wkbNone:
        return 'None'
    elif type == ogr.wkbLinearRing:
        return 'LinearRing'
    else:
        return 'Unknown'


def sqlOutput(tableName, valueDict):
     sqlStatement = "INSERT INTO %s %s VALUES %s;" % (tableName, tuple((valueDict.keys())),tuple(valueDict.values()))
     print sqlStatement

def sqlCreateTables():
    processColumns = "SearchPath VARCHAR, LaunchPath VARCHAR, UserHome  VARCHAR, IgnoredString VARCHAR, DirCount  INTEGER, FileCount INTEGER, Timestamp VARCHAR"
    tables = ('process',) #,'dataset','layer','raster','band')

    for table in tables:
        sqlStatement = "CREATE TABLE %s (%s);" % (table, processColumns)
        print sqlStatement

def getFileStats(filepath):
    mode, ino, dev, nlink, user_id, group_id, file_size, \
            time_accessed, time_modified, time_created = os.stat(filepath)

    if os.path.isfile(filepath):
        file_type = "File"
    else: 
        file_type = "Directory"
    try:
        import pwd # not available on all platforms
        userinfo = pwd.getpwuid(user_id)
    except (ImportError, KeyError):
        user_name = "N/A"
        user_full_name = "N/A"
    else:
        user_name = userinfo[0]
        user_full_name = userinfo[4]
    full_path = os.path.abspath(filepath)
    md5_key = (full_path, user_name, file_size, time_modified, time_created)
    md5_digest = getMd5HexDigest(md5_key)
    resultsFileStats = {
            'fullPath': full_path, 
            'userId': user_id, 
            'groupId': group_id, 
            'fileSize': file_size, 
            'timeAccessed': str(time_accessed), 
            'timeModified': str(time_modified), 
            'timeCreated': str(time_created), 
            'fileType': file_type, 
            'userName': user_name, 
            'userFullName': user_full_name, 
            'uniqueDigest': md5_digest
    }
    return resultsFileStats

def getMd5HexDigest(encodeString):
    m = hashlib.md5()
    m.update(str(encodeString))
    return m.hexdigest()

if __name__ == '__main__':
    parser = OptionParser(usage="gdalogr_catalog.py /path/to/search -f output.json")
    parser.add_option("-f","--file", action="store", type="string", dest="outfile", 
            help="Output file (if specified, no stdout)" )

    (options, args) = parser.parse_args()

    startpath = None
    if len(args) == 1:
        startpath = args[0]

    if not startpath or not os.path.exists(startpath):
        parser.error("Please supply an valid search directory")

    catalog = getCatalog(startpath)
    if options.outfile:
        with open(options.outfile,'w') as out:
            out.write(json.dumps(catalog, indent=2))
            log.info("%s written" % options.outfile)
    else:
        print json.dumps(catalog, indent=2)
