#!/usr/bin/env python

# gdalogr_catalogue.py
# Purpose: Catalog all raster & vector datasources/layers found in a directory tree
# Usage: python gdalogr_catalogue.py <search path>
# Sends hierarchical XML to stdout
# Requires GDAL/OGR libraries: http://pypi.python.org/pypi/GDAL
# Requires Python > 2.3 for itertools.tee function.

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

'''
TODO
# DONE - higher level attributes about process: num of files, dirs search, timestamp
# DONE - filesize, user/owner, moddate/timestamp for entries
- extent values to GML or basic WKT bbox
# DONE - decide on checksum process for determining changes
- decide on process -> datasource linking (timestamp?) for top level relations
'''
from osgeo import gdal
from osgeo import osr
from osgeo import ogr

import logging
import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser, OptionGroup
from string import strip
from time import asctime
from ElementTree_pretty import prettify

logging.basicConfig( stream=sys.stderr, level=logging.DEBUG )
log = logging.getLogger("catalog")

def startup(startpath):
    skiplist = ['.svn','.shx','.dbf', '.prj', '.aux.xml', '.e00', '.adf']
    gdal.PushErrorHandler()
    pathwalker = os.walk(startpath)
    counterraster = 0
    countervds = 0

    dirlist, filelist = [], []
    for eachpath in pathwalker:
        startdir = eachpath[0]

        dirlist += eachpath[1]
        filelist += eachpath[2]

        alldirs = eachpath[1]
        allfiles = eachpath[2]

        for eachdir in alldirs:
            currentdir = os.path.join(startdir,eachdir)
            raster,vector = None, None
            if (not skipfile(currentdir,skiplist)):
                raster,vector = tryopends(currentdir)
            if raster:
                try:
                    resultsraster,resultsFileStats = processraster(raster,counterraster,currentdir)
                    xmlraster = outputraster(resultsraster, counterraster, countervds, resultsFileStats, xmlroot)
                    counterraster += 1
                except NotGeographic:
                    pass
            if vector:
                try:
                    resultsvds,resultsFileStats = processvds(vector,countervds,currentdir)
                    xmlvector = outputvector(resultsvds, counterraster, countervds, resultsFileStats, xmlroot)
                    countervds += 1
                except NotGeographic:
                    pass
        for eachfile in allfiles:
            currentfile = "/".join([startdir, eachfile])
            raster, vector = None, None
            if (not skipfile(currentfile,skiplist)):
                raster, vector = tryopends(currentfile)
            if raster:
                try:
                    resultsraster,resultsFileStats = processraster(raster, counterraster, currentfile)
                    xmlraster = outputraster(resultsraster, counterraster, countervds, resultsFileStats, xmlroot)
                    counterraster += 1
                except NotGeographic:
                    pass
            if vector:
                if not skipfile(vector.GetName(), skiplist):
                    try:
                        resultsvds,resultsFileStats = processvds(vector, countervds, currentfile)
                        xmlvector = outputvector(resultsvds,counterraster,countervds,resultsFileStats,xmlroot)
                        countervds += 1
                    except NotGeographic:
                        pass

    xmlcatalog = appendXML(xmlroot, "CatalogueProcess")
    appendXML(xmlcatalog, "SearchPath", startpath)
    appendXML(xmlcatalog, "LaunchPath", os.getcwd())
    appendXML(xmlcatalog, "UserHome", os.getenv("HOME"))
    appendXML(xmlcatalog, "IgnoredStrings", str(skiplist))
    appendXML(xmlcatalog, "DirCount", str(len(dirlist)))
    appendXML(xmlcatalog, "FileCount", str(len(filelist)))
    appendXML(xmlcatalog, "Timestamp", asctime())
    
    if options.printSql: 
        processValues = {'SearchPath':startpath,'LaunchPath':os.getcwd(),'UserHome':os.getenv("HOME"),'IgnoredString':" ".join(map(str, skiplist)),'DirCount':int(len(dirlist)),'FileCount':int(len(filelist)),'Timestamp':asctime()}
        print sqlOutput('process',processValues)

class NotGeographic(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        log.warn(message)

def startXML():
    xmlroot = ET.Element("DataCatalogue")
    return xmlroot

def appendXML(elementroot, subelement, subelstring=None):
    newelement = ET.SubElement(elementroot, subelement)
    newelement.text = subelstring
    return newelement

def writeXML(xmlroot, outfile):
    xmltree = ET.ElementTree(xmlroot)
    xmltree.write(outfile)

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
    skip_drivers = ['AVCBin',] # GetFeatureCount and Extent are Way too slow and called many times TODO
    try:
        dsgdal = gdal.OpenShared(filepath)
    except gdal.GDALError:
        pass

    try:
        dsogr = ogr.OpenShared(filepath)
    except ogr.OGRError:
        pass

    if dsogr and dsogr.GetDriver().GetName() in skip_drivers:
        dsogr = False

    return dsgdal, dsogr

def extentToLatLon(extent, proj):
    from pyproj import Proj
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
    extent = (geotrans[0], geotrans[3] + (geotrans[5] * rastery), geotrans[0] + ( geotrans[1] * rasterx ), geotrans[3], )
    latlon_extent = extentToLatLon(extent, proj)

    resultsbands = {}
    resultsFileStats = fileStats(currentpath)
    for bandnum in range(bandcount):
        band = raster.GetRasterBand(bandnum+1)
        min, max = band.ComputeRasterMinMax(1) #approx_ok=1
        overviews = band.GetOverviewCount()
        resultseachband = {'bandId': str(bandnum+1), 'min': str(min),'max': str(max), 'overviews': str(overviews)}
        resultseachbandShort = {'bandId': bandnum+1, 'min': min,'max': max, 'overviews': str(overviews)}
        resultsbands[str(bandnum+1)] = resultseachband
        if options.printSql: 
            print sqlOutput('band',resultseachbandShort)
    resultsraster = { 
            'bands': resultsbands, 
            'rasterId': str(counterraster), 
            'name': rastername, 
            'bandcount': str(bandcount), 
            'geotrans': str(geotrans), 
            'driver': str(driver), 
            'rasterX': str(rasterx), 
            'rasterY': str(rastery), 
            'projection': wkt,
            'proj': proj,
            'latlonextent': strip(str(latlon_extent),"()"),
            'extent': strip(str(extent), "()")
    }
    resultsrasterShort =  {
            'rasterId':counterraster, 
            'name': rastername, 
            'bandcount': bandcount, 
            'geotrans': str(geotrans), 
            'driver': driver, 
            'rasterX': rasterx, 
            'rasterY': rastery, 
            'projection': wkt
    }
    if options.printSql: 
        print sqlOutput('raster',resultsrasterShort)

    return resultsraster, resultsFileStats
  
def outputraster(resultsraster, counterraster, countervds, resultsFileStats, xmlroot):
    xmlraster = appendXML(xmlroot, "RasterData")
    statfileStats = outputFileStats(resultsFileStats, xmlraster)
    for rasteritem, rastervalue in resultsraster.iteritems(): # for each raster attribute
        if rasteritem <> 'bands':
            appendXML(xmlraster, rasteritem, rastervalue)
        if rasteritem == 'bands':
            for banditem, bandvalue in rastervalue.iteritems(): # for each band
                xmlband = appendXML(xmlraster, "RasterBand")
                for banditemdetails, bandvaluedetails in bandvalue.iteritems():
                    appendXML(xmlband, banditemdetails, bandvaluedetails)
    return True

def processvds(vector, countervds,currentpath):
    vdsname = vector.GetName()
    log.debug(vdsname)
    vdsformat = vector.GetDriver().GetName()
    vdslayercount = vector.GetLayerCount()
    resultslayers = {}
    resultsFileStats = fileStats(currentpath)
    if resultsFileStats['fileType'] == 'Directory' and vdsformat == 'ESRI Shapefile':
        raise NotGeographic("Just a directory of shapefiles; nothing to see here")

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

        # the following throws all the attributes into dictionaries of attributes, 
        # some of which are other dictionaries
        # resultseachlayer = 1 layer attributes
        # resultslayers = dict. of all layers and their attributes
        # resultsvds = datasource attributes
        # resultsvector = dict of datasource attributes, plus a dict of all layers
        # Note all get saved as strings, which isn't what you'd want for SQL output
        resultseachlayer = {
            'layerId': str(layernum+1), 
            'name': layername, 
            'proj': layerproj, 
            'featuretype': str(layerftype), 
            'featurecount': str(layerfcount), 
            'extent': strip(str(layerextent),"()"),
            'latlonextent': strip(str(latlon_extent),"()"),
        }
        resultslayers[str(layernum+1)] = resultseachlayer
        sqlstringvlay = "INSERT INTO layer %s VALUES %s;" % (
                ('layerId','datasourceId','name','featurecount','extent'), 
                (layernum+1,countervds,layername,int(layerfcount),layerextentraw))
        if options.printSql: 
            print sqlOutput('layer',resultseachlayer)
        #if (layerftype <> 'UNKNOWN'):
        #    Mapping(vector,layerextentraw,layername,layerftype) # mapping test
    resultsvds = { 
            'datasourceId': str(countervds), 
            'name': vdsname, 
            'format': vdsformat, 
            'layercount': str(vdslayercount) 
    }
    sqlstringvds = "INSERT INTO datasource %s VALUES %s;" % (
            ('datasourceId','name','format','layercount'), 
            (countervds, vdsname, vdsformat, int(vdslayercount)))
    resultsvector = { 'resultsvds': resultsvds, 'resultslayers': resultslayers } 
    if options.printSql: 
        print sqlOutput('dataset',resultsvds)

    return resultsvector,resultsFileStats

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

def outputvector(resultsvector, counterraster, countervds, resultsFileStats,xmlroot):
    xmlvector = appendXML(xmlroot, "VectorData")
    statfileStats = outputFileStats(resultsFileStats, xmlvector)
    for vectoritem, vectorvalue in resultsvector.iteritems(): # resultsvector includes two dictionaries
        if vectoritem == 'resultslayers':
            for layeritem, layervalue in vectorvalue.iteritems(): # vectorvalue contains a dictionary of the layers
                xmlvectorlayer = appendXML(xmlvector, "VectorLayer")
                for layeritemdetails, layervaluedetails in layervalue.iteritems(): # layervalue contains layer attributes
                    appendXML(xmlvectorlayer, layeritemdetails, layervaluedetails)
        else:
            for vectordsitem, vectordsvalue in vectorvalue.iteritems(): # vectorvalue contains datasource attributes
                appendXML(xmlvector, vectordsitem, vectordsvalue)
    return True

def sqlOutput(tableName, valueDict):
     sqlStatement = "INSERT INTO %s %s VALUES %s;" % (tableName, tuple((valueDict.keys())),tuple(valueDict.values()))
     print sqlStatement

def sqlCreateTables():
    processColumns = "SearchPath VARCHAR, LaunchPath VARCHAR, UserHome  VARCHAR, IgnoredString VARCHAR, DirCount  INTEGER, FileCount INTEGER, Timestamp VARCHAR"
    tables = ('process',) #,'dataset','layer','raster','band')

    for table in tables:
        sqlStatement = "CREATE TABLE %s (%s);" % (table, processColumns)
        print sqlStatement

### TODO functions below...

def xmlDtdOutput():
    import zipfiles
    # output dtd that corresponds to the xml, or is it schema?

def checkZip(currentfile):
    import zipfiles
    # check if it can read zips

def openZip(currentfile):
    import zipfiles
    # extract files and catalogue them

def fileStats(filepath):
    mode, ino, dev, nlink, user_id, group_id, file_size, time_accessed, time_modified, time_created = os.stat(filepath)
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
    resultsFileStats = {'fullPath': str(full_path), 'userId': str(user_id), 'groupId': str(group_id), 'fileSize': str(file_size), 'timeAccessed': str(time_accessed), 'timeModified': str(time_modified), 'timeCreated': str(time_created), 'fileType': file_type, 'userName': user_name, 'userFullName': user_full_name, 'uniqueDigest': md5_digest}
    return resultsFileStats

def outputFileStats(resultsFileStats, xmlroot):
    xmlfilestats = appendXML(xmlroot, "FileStats")
    for statitem, statvalue in resultsFileStats.iteritems():
        appendXML(xmlfilestats, statitem, statvalue)
    return True

def outputXml(root,newelement):
    SubElement(root,newelement)
    return 
  
def getMd5HexDigest(encodeString):
    import md5
    m = md5.new()
    m.update(str(encodeString))
    return m.hexdigest()

class Mapping:
    def __init__(self,datasource,extent,layername,layerftype):
        import mapscript
        from time import time
        tmap = mapscript.mapObj()
        print "checkpoint 1"
        map.setSize(400,400)
        #ext = rectObj(-180,-90,180,90)
        ext = mapscript.rectObj(extent[0],extent[2],extent[1],extent[3]) # some trouble with some bad extents in my test data
        map.extent = ext
        map.units = mapscript.MS_DD # should be programmatically set
        lay = mapscript.layerObj(map)
        lay.name = "Autolayer"
        lay.units = mapscript.MS_DD
        if (layerftype == 'RASTER'):
            lay.data = datasource.GetDescription()
        else:
            lay.data = datasource.GetName()
        print lay.data
        lay.status = mapscript.MS_DEFAULT
        cls = mapscript.classObj(lay)
        sty = mapscript.styleObj()
        col = mapscript.colorObj(0,0,0)
        #symPoint = mapscript.symbolObj
        map.setSymbolSet("symbols/symbols_basic.sym")
        if (layerftype == 'POINT'): 
            lay.type = mapscript.MS_LAYER_POINT
            sty.setSymbolByName = "achteck"
            sty.width = 100
            sty.color = col
        elif (layerftype == 'LINE'): 
            lay.type = mapscript.MS_LAYER_LINE
            sty.setSymbolByName = "punkt"
            sty.width = 5
            sty.color = col
        elif (layerftype == 'POLYGON'): 
            lay.type = mapscript.MS_LAYER_POLYGON
            sty.setSymbolByName = "circle"
            sty.width = 10
            sty.outlinecolor = col
        elif (layerftype == 'RASTER'): 
            lay.type = mapscript.MS_LAYER_RASTER
            sty.setSymbolByName = "squares"
            sty.size = 10
            sty.color = col
        #sty.setSymbolByName(map,symname)
        #sty.size = symsize
        cls.insertStyle(sty)
        try:
            img = map.draw()
            img.save(str(time()) + "auto.gif")
            map.save(str(time()) + "auto.map")
        except MapServerError:
            return None
        # add layer
        # assign datasource to layer
        # add basic styling
        # apply styling to layer
        # open output image
        # write, close, cleanup


if __name__ == '__main__':
    parser = OptionParser(usage="gdalogr_catalog.py [options] -d /path/to/search")
    parser.add_option("-d","--dir", action="store", type="string", dest="directory", 
            help="Top level folder to start scanning from")
    parser.add_option("-f","--file", action="store", type="string", dest="outfile", 
            help="Output file (if specified, no stdout)" )

    group = OptionGroup(parser, "Hack Options", "May not function without advanced knowledge")
    group.add_option("-s","--sql", action="store_true", dest="printSql", 
            help="Output results in SQL INSERT statements instead of XML")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    startpath = options.directory
    if not startpath and len(args) >= 1:
        startpath = args[0]

    if not startpath or not os.path.exists(startpath):
        parser.error("Please supply an valid search directory")

    xmlroot = startXML()
    startup(startpath)
    if options.outfile:
        #writeXML(xmlroot, options.outfile)
        with open(options.outfile,'w') as out:
            out.write(prettify(xmlroot))
            log.info("%s written" % options.outfile)
    else:
        print prettify(xmlroot)
