#!/usr/bin/python
# -*- coding: utf-8 -*-

# ******************************************************************************
#  $Id: gdal2tiles.py 27349 2014-05-16 18:58:51Z rouault $
#
# Project:  Google Summer of Code 2007, 2008 (http://code.google.com/soc/)
# Support:  BRGM (http://www.brgm.fr)
# Purpose:  Convert a raster into TMS (Tile Map Service) tiles in a directory.
#           - generate Google Earth metadata (KML SuperOverlay)
#           - generate simple HTML viewer based on Google Maps and OpenLayers
#           - support of global tiles (Spherical Mercator) for compatibility
#               with interactive web maps a la Google Maps
# Author:   Klokan Petr Pridal, klokan at klokan dot cz
# Web:      http://www.klokan.cz/projects/gdal2tiles/
# GUI:      http://www.maptiler.org/
#
###############################################################################
# Copyright (c) 2008, Klokan Petr Pridal
# Copyright (c) 2010-2013, Even Rouault <even dot rouault at mines-paris dot org>
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
# ******************************************************************************

import sys

try:
    from osgeo import gdal
    from osgeo import osr
except:
    import gdal
    print 'You are using "old gen" bindings. gdal2tiles needs "new gen" bindings.'
    sys.exit(1)

import os
import math

try:
    from PIL import Image
    import numpy
    import osgeo.gdal_array as gdalarray
except:

    # 'antialias' resampling is not available

    pass

__version__ = '$Id: gdal2tiles.py 27349 2014-05-16 18:58:51Z rouault $'

resampling_list = (
    'average',
    'near',
    'bilinear',
    'cubic',
    'cubicspline',
    'lanczos',
    'antialias',
    )

# =============================================================================
# =============================================================================
# =============================================================================

__doc__globalmaptiles = \
    """
globalmaptiles.py

Global Map Tiles as defined in Tile Map Service (TMS) Profiles
==============================================================

Functions necessary for generation of global tiles used on the web.

More info at:

http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
http://wiki.osgeo.org/wiki/WMS_Tiling_Client_Recommendation
http://msdn.microsoft.com/en-us/library/bb259689.aspx
http://code.google.com/apis/maps/documentation/overlays.html#Google_Maps_Coordinates

Created by Klokan Petr Pridal on 2008-07-03.
Google Summer of Code 2008, project GDAL2Tiles for OSGEO.

In case you use this class in your product, translate it to another language
or find it usefull for your project please let me know.
My email: klokan at klokan dot cz.
I would like to know where it was used.

Class is available under the open-source GDAL license (www.gdal.org).
"""

import math

MAXZOOMLEVEL = 32

# ---------------------
# TODO: Finish Zoomify implemtentation!!!

class Zoomify(object):

    """
    Tiles compatible with the Zoomify viewer
    ----------------------------------------
    """

    def __init__(
        self,
        width,
        height,
        tilesize=256,
        tileformat='jpg',
        ):
        """Initialization of the Zoomify tile tree"""

        self.tilesize = tilesize
        self.tileformat = tileformat
        imagesize = (width, height)
        tiles = (math.ceil(width / tilesize), math.ceil(height
                 / tilesize))

        # Size (in tiles) for each tier of pyramid.

        self.tierSizeInTiles = []
        self.tierSizeInTiles.push(tiles)

        # Image size in pixels for each pyramid tierself

        self.tierImageSize = []
        self.tierImageSize.append(imagesize)

        while imagesize[0] > tilesize or imageSize[1] > tilesize:
            imagesize = (math.floor(imagesize[0] / 2),
                         math.floor(imagesize[1] / 2))
            tiles = (math.ceil(imagesize[0] / tilesize),
                     math.ceil(imagesize[1] / tilesize))
            self.tierSizeInTiles.append(tiles)
            self.tierImageSize.append(imagesize)

        self.tierSizeInTiles.reverse()
        self.tierImageSize.reverse()

        # Depth of the Zoomify pyramid, number of tiers (zoom levels)

        self.numberOfTiers = len(self.tierSizeInTiles)

        # Number of tiles up to the given tier of pyramid.

        self.tileCountUpToTier = []
        self.tileCountUpToTier[0] = 0
        for i in range(1, self.numberOfTiers + 1):
            self.tileCountUpToTier.append(self.tierSizeInTiles[i
                    - 1][0] * self.tierSizeInTiles[i - 1][1]
                    + self.tileCountUpToTier[i - 1])

    def tilefilename(
        self,
        x,
        y,
        z,
        ):
        """Returns filename for tile with given coordinates"""

        tileIndex = x + y * self.tierSizeInTiles[z][0] \
            + self.tileCountUpToTier[z]
        return os.path.join('TileGroup%.0f' % math.floor(tileIndex
                            / 256), '%s-%s-%s.%s' % (z, x, y,
                            self.tileformat))


# =============================================================================
# =============================================================================
# =============================================================================

class GDAL2Tiles(object):

    # -------------------------------------------------------------------------

    def process(self):
        """The main processing function, runs all the main steps of processing"""

        # Opening and preprocessing of the input file

        self.open_input()

        # Generation of main metadata files and HTML viewers

        self.generate_metadata()

        # Generation of the lowest tiles

        self.generate_base_tiles()

        # Generation of the overview tiles (higher in the pyramid)

        self.generate_overview_tiles()

    # -------------------------------------------------------------------------

    def error(self, msg, details=''):
        """Print an error message and stop the processing"""

        if details:
            self.parser.error(msg + '''

''' + details)
        else:
            self.parser.error(msg)

    # -------------------------------------------------------------------------

    def progressbar(self, complete=0.0):
        """Print progressbar for float value 0..1"""

        gdal.TermProgress_nocb(complete)

    # -------------------------------------------------------------------------

    def stop(self):
        """Stop the rendering immediately"""

        self.stopped = True

    # -------------------------------------------------------------------------

    def __init__(self, arguments):
        """Constructor function - initialization"""

        self.stopped = False
        self.input = None
        self.output = None

        # Tile format

        self.tilesize = 256
        self.tiledriver = 'PNG'
        self.tileext = 'png'

        # Should we read bigger window of the input raster and scale it down?
        # Note: Modified leter by open_input()
        # Not for 'near' resampling
        # Not for Wavelet based drivers (JPEG2000, ECW, MrSID)

        self.scaledquery = True

        # How big should be query window be for scaling down
        # Later on reset according the chosen resampling algorightm

        self.querysize = 4 * self.tilesize

        # Should we use Read on the input file for generating overview tiles?
        # Note: Modified later by open_input()
        # Otherwise the overview tiles are generated from existing underlying tiles

        self.overviewquery = False

        # RUN THE ARGUMENT PARSER:

        self.optparse_init()
        (self.options, self.args) = \
            self.parser.parse_args(args=arguments)
        if not self.args:
            self.error('No input file specified')

        # POSTPROCESSING OF PARSED ARGUMENTS:

        # Workaround for old versions of GDAL

        try:
            if self.options.verbose and self.options.resampling \
                == 'near' or gdal.TermProgress_nocb:
                pass
        except:
            self.error('This version of GDAL is not supported. Please upgrade to 1.6+.'
                       )

            # ,"You can try run crippled version of gdal2tiles with parameters: -v -r 'near'")

        # Is output directory the last argument?

        # Test output directory, if it doesn't exist

        if os.path.isdir(self.args[-1]) or len(self.args) > 1 \
            and not os.path.exists(self.args[-1]):
            self.output = self.args[-1]
            self.args = self.args[:-1]

        # More files on the input not directly supported yet

        if len(self.args) > 1:
            self.error('Processing of several input files is not supported.'
                       ,
                       """Please first use a tool like gdal_vrtmerge.py or gdal_merge.py on the files:
gdal_vrtmerge.py -o merged.vrt %s"""
                        % ' '.join(self.args))

            # TODO: Call functions from gdal_vrtmerge.py directly

        self.input = self.args[0]

        # Default values for not given options

        if not self.output:

            # Directory with input filename without extension in actual directory

            self.output = \
                os.path.splitext(os.path.basename(self.input))[0]

        if not self.options.title:
            self.options.title = os.path.basename(self.input)

        if self.options.url and not self.options.url.endswith('/'):
            self.options.url += '/'
        if self.options.url:
            self.options.url += os.path.basename(self.output) + '/'

        # Supported options

        self.resampling = None

        if self.options.resampling == 'average':
            try:
                if gdal.RegenerateOverview:
                    pass
            except:
                self.error("'average' resampling algorithm is not available."
                           ,
                           "Please use -r 'near' argument or upgrade to newer version of GDAL."
                           )
        elif self.options.resampling == 'antialias':

            try:
                if numpy:
                    pass
            except:
                self.error("'antialias' resampling algorithm is not available."
                           ,
                           'Install PIL (Python Imaging Library) and numpy.'
                           )
        elif self.options.resampling == 'near':

            self.resampling = gdal.GRA_NearestNeighbour
            self.querysize = self.tilesize
        elif self.options.resampling == 'bilinear':

            self.resampling = gdal.GRA_Bilinear
            self.querysize = self.tilesize * 2
        elif self.options.resampling == 'cubic':

            self.resampling = gdal.GRA_Cubic
        elif self.options.resampling == 'cubicspline':

            self.resampling = gdal.GRA_CubicSpline
        elif self.options.resampling == 'lanczos':

            self.resampling = gdal.GRA_Lanczos

        # User specified zoom levels

        self.tminz = None
        self.tmaxz = None
        if self.options.zoom:
            minmax = self.options.zoom.split('-', 1)
            minmax.extend([''])
            (min, max) = minmax[:2]
            self.tminz = int(min)
            if max:
                self.tmaxz = int(max)
            else:
                self.tmaxz = int(min)

        # Output the results

        if self.options.verbose:
            print ('Options:', self.options)
            print ('Input:', self.input)
            print ('Output:', self.output)
            print 'Cache: %s MB' % (gdal.GetCacheMax() / 1024 / 1024)
            print ''

    # -------------------------------------------------------------------------

    def optparse_init(self):
        """Prepare the option parser for input (argv)"""

        from optparse import OptionParser, OptionGroup
        usage = 'Usage: %prog [options] input_file(s) [output]'
        p = OptionParser(usage, version='%prog ' + __version__)
        p.add_option(
            '-r',
            '--resampling',
            dest='resampling',
            type='choice',
            choices=resampling_list,
            help="Resampling method (%s) - default 'average'"
                % ','.join(resampling_list),
            )
        p.add_option('-z', '--zoom', dest='zoom',
                     help="Zoom levels to render (format:'2-5' or '10')."
                     )
        p.add_option('-e', '--resume', dest='resume',
                     action='store_true',
                     help='Resume mode. Generate only missing files.')
        p.add_option('-a', '--srcnodata', dest='srcnodata',
                     metavar='NODATA',
                     help='NODATA transparency value to assign to the input data'
                     )
        p.add_option('-v', '--verbose', action='store_true',
                     dest='verbose',
                     help='Print status messages to stdout')

        # TODO: MapFile + TileIndexes per zoom level for efficient MapServer WMS
            # g = OptionGroup(p, "WMS MapServer metadata", "Options for generated mapfile and tileindexes for MapServer")
            # g.add_option("-i", "--tileindex", dest='wms', action="store_true"
            #                 help="Generate tileindex and mapfile for MapServer (WMS)")
            # p.add_option_group(g)

        p.set_defaults(
            verbose=False,
            url='',
            webviewer='all',
            copyright='',
            resampling='average',
            resume=False,
            )

        self.parser = p

    # -------------------------------------------------------------------------

    def open_input(self):
        """Initialization of the input raster, reprojection if necessary"""

        gdal.AllRegister()

        # Initialize necessary GDAL drivers

        self.out_drv = gdal.GetDriverByName(self.tiledriver)
        self.mem_drv = gdal.GetDriverByName('MEM')

        if not self.out_drv:
            raise Exception("The '%s' driver was not found, is it available in this GDAL build?"
                            , self.tiledriver)
        if not self.mem_drv:
            raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?"
                            )

        # Open the input file

        if self.input:
            self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
        else:
            raise Exception('No input file was specified')

        if self.options.verbose:
            print ('Input file:', '( %sP x %sL - %s bands)'
                   % (self.in_ds.RasterXSize, self.in_ds.RasterYSize,
                   self.in_ds.RasterCount))

        if not self.in_ds:

            # Note: GDAL prints the ERROR message too

            self.error("It is not possible to open the input file '%s'."
                        % self.input)

        # Read metadata from the input file

        if self.in_ds.RasterCount == 0:
            self.error("Input file '%s' has no raster band"
                       % self.input)

        if self.in_ds.GetRasterBand(1).GetRasterColorTable():

            # TODO: Process directly paletted dataset by generating VRT in memory

            self.error('Please convert this file to RGB/RGBA and run gdal2tiles on the result.'
                       ,
                       """From paletted file you can create RGBA file (temp.vrt) by:
gdal_translate -of vrt -expand rgba %s temp.vrt
then run:
gdal2tiles temp.vrt"""
                       % self.input)

        # Get NODATA value

        self.in_nodata = []
        for i in range(1, self.in_ds.RasterCount + 1):
            if self.in_ds.GetRasterBand(i).GetNoDataValue() != None:
                self.in_nodata.append(self.in_ds.GetRasterBand(i).GetNoDataValue())
        if self.options.srcnodata:
            nds = list(map(float, self.options.srcnodata.split(',')))
            if len(nds) < self.in_ds.RasterCount:
                self.in_nodata = (nds
                                  * self.in_ds.RasterCount)[:self.in_ds.RasterCount]
            else:
                self.in_nodata = nds

        if self.options.verbose:
            print 'NODATA: %s' % self.in_nodata

        #
        # Here we should have RGBA input dataset opened in self.in_ds
        #

        if self.options.verbose:
            print ('Preprocessed file:', '( %sP x %sL - %s bands)'
                   % (self.in_ds.RasterXSize, self.in_ds.RasterYSize,
                   self.in_ds.RasterCount))

        # Are the reference systems the same? Reproject if necessary.

        self.out_ds = None

        if not self.out_ds:
            self.out_ds = self.in_ds

        #
        # Here we should have a raster (out_ds) in the correct Spatial Reference system
        #

        # Get alpha band (either directly or from NODATA value)

        self.alphaband = self.out_ds.GetRasterBand(1).GetMaskBand()
        if self.alphaband.GetMaskFlags() & gdal.GMF_ALPHA \
            or self.out_ds.RasterCount == 4 or self.out_ds.RasterCount \
            == 2:

            # TODO: Better test for alpha band in the dataset

            self.dataBandsCount = self.out_ds.RasterCount - 1
        else:
            self.dataBandsCount = self.out_ds.RasterCount

        # Read the georeference

        self.out_gt = self.out_ds.GetGeoTransform()

        # originX, originY = self.out_gt[0], self.out_gt[3]
        # pixelSize = self.out_gt[1] # = self.out_gt[5]

        # Test the size of the pixel

        # MAPTILER - COMMENTED
        # if self.out_gt[1] != (-1 * self.out_gt[5]) and self.options.profile != 'raster':
            # TODO: Process corectly coordinates with are have swichted Y axis (display in OpenLayers too)
            # self.error("Size of the pixel in the output differ for X and Y axes.")

        # Report error in case rotation/skew is in geotransform (possible only in 'raster' profile)

        if (self.out_gt[2], self.out_gt[4]) != (0, 0):
            self.error('Georeference of the raster contains rotation or skew. Such raster is not supported. Please use gdalwarp first.'
                       )

            # TODO: Do the warping in this case automaticaly

        #
        # Here we expect: pixel is square, no rotation on the raster
        #

        # Output Bounds - coordinates in the output SRS

        self.ominx = self.out_gt[0]
        self.omaxx = self.out_gt[0] + self.out_ds.RasterXSize \
            * self.out_gt[1]
        self.omaxy = self.out_gt[3]
        self.ominy = self.out_gt[3] - self.out_ds.RasterYSize \
            * self.out_gt[1]

        # Note: maybe round(x, 14) to avoid the gdal_translate behaviour, when 0 becomes -1e-15

        if self.options.verbose:
            print ('Bounds (output srs):', round(self.ominx, 13),
                   self.ominy, self.omaxx, self.omaxy)

        #
        # Calculating ranges for tiles in different zoom levels
        #

        log2 = lambda x: math.log10(x) / math.log10(2)  # log2 (base 2 logarithm)

        self.nativezoom = \
            int(max(math.ceil(log2(self.out_ds.RasterXSize
                / float(self.tilesize))),
                math.ceil(log2(self.out_ds.RasterYSize
                / float(self.tilesize)))))

        if self.tmaxz < self.nativezoom:
            self.tmaxz = self.nativezoom

        if self.options.verbose:
            print ('Native zoom of the raster:', self.nativezoom)

        # Get the minimal zoom level (whole raster in one tile)

        if self.tminz == None:
            self.tminz = 0

        # Get the maximal zoom level (native resolution of the raster)

        if self.tmaxz == None:
            self.tmaxz = self.nativezoom

        # Generate table with min max tile coordinates for all zoomlevels

        self.tminmax = list(range(0, self.tmaxz + 1))
        self.tsize = list(range(0, self.tmaxz + 1))
        for tz in range(0, self.tmaxz + 1):
            tsize = 2.0 ** (self.nativezoom - tz) * self.tilesize
            (tminx, tminy) = (0, 0)
            tmaxx = int(math.ceil(self.out_ds.RasterXSize / tsize)) \
                - 1
            tmaxy = int(math.ceil(self.out_ds.RasterYSize / tsize)) \
                - 1
            self.tsize[tz] = math.ceil(tsize)
            self.tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

        # Function which generates SWNE in LatLong for given tile

        self.tileswne = lambda x, y, z: (0, 0, 0, 0)

    # -------------------------------------------------------------------------

    def generate_metadata(self):
        """Generation of main metadata files and HTML viewers (metadata related to particular tiles are generated during the tile processing)."""

        if not os.path.exists(self.output):
            os.makedirs(self.output)

        (west, south) = (self.ominx, self.ominy)
        (east, north) = (self.omaxx, self.omaxy)

        self.swne = (south, west, north, east)

        # Generate openlayers.html

        if self.options.webviewer in ('all', 'openlayers'):
            if not self.options.resume \
                or not os.path.exists(os.path.join(self.output,
                    'openlayers.html')):
                f = open(os.path.join(self.output, 'openlayers.html'
                         ), 'w')
                f.write(self.generate_openlayers())
                f.close()

        # Generate tilemapresource.xml.

        if not self.options.resume \
            or not os.path.exists(os.path.join(self.output,
                                  'tilemapresource.xml')):
            f = open(os.path.join(self.output, 'tilemapresource.xml'),
                     'w')
            f.write(self.generate_tilemapresource())
            f.close()

    # -------------------------------------------------------------------------

    def generate_base_tiles(self):
        """Generation of the base tiles (the lowest in the pyramid) directly from the input raster"""

        print 'Generating Base Tiles:'

        if self.options.verbose:

            # mx, my = self.out_gt[0], self.out_gt[3] # OriginX, OriginY
            # px, py = self.mercator.MetersToPixels( mx, my, self.tmaxz)
            # print "Pixel coordinates:", px, py, (mx, my)

            print ''
            print 'Tiles generated from the max zoom level:'
            print '----------------------------------------'
            print ''

        # Set the bounds

        (tminx, tminy, tmaxx, tmaxy) = self.tminmax[self.tmaxz]

        # Just the center tile
        # tminx = tminx+ (tmaxx - tminx)/2
        # tminy = tminy+ (tmaxy - tminy)/2
        # tmaxx = tminx
        # tmaxy = tminy

        ds = self.out_ds
        tilebands = self.dataBandsCount + 1
        querysize = self.querysize

        if self.options.verbose:
            print ('dataBandsCount: ', self.dataBandsCount)
            print ('tilebands: ', tilebands)

        # print tminx, tminy, tmaxx, tmaxy

        tcount = (1 + abs(tmaxx - tminx)) * (1 + abs(tmaxy - tminy))

        # print tcount

        ti = 0

        tz = self.tmaxz
        yrange = range(tminy, tmaxy + 1)

        for ty in yrange:
            for tx in range(tminx, tmaxx + 1):

                if self.stopped:
                    break
                ti += 1
                tilefilename = os.path.join(self.output, str(tz),
                        str(tx), '%s.%s' % (ty, self.tileext))
                if self.options.verbose:
                    print (ti, '/', tcount, tilefilename)  # , "( TileMapService: z / x / y )"

                if self.options.resume and os.path.exists(tilefilename):
                    if self.options.verbose:
                        print 'Tile generation skiped because of --resume'
                    else:
                        self.progressbar(ti / float(tcount))
                    continue

                # Create directories for the tile

                if not os.path.exists(os.path.dirname(tilefilename)):
                    os.makedirs(os.path.dirname(tilefilename))

                # print "\tgdalwarp -ts 256 256 -te %s %s %s %s %s %s_%s_%s.tif" % ( b[0], b[1], b[2], b[3], "tiles.vrt", tz, tx, ty)

                # Don't scale up by nearest neighbour, better change the querysize
                # to the native resolution (and return smaller query tile) for scaling

                tsize = int(self.tsize[tz])  # tilesize in raster coordinates for actual zoom
                xsize = self.out_ds.RasterXSize  # size of the raster in pixels
                ysize = self.out_ds.RasterYSize
                if tz >= self.nativezoom:
                    querysize = self.tilesize  # int(2**(self.nativezoom-tz) * self.tilesize)

                rx = tx * tsize
                rxsize = 0
                if tx == tmaxx:
                    rxsize = xsize % tsize
                if rxsize == 0:
                    rxsize = tsize

                rysize = 0
                if ty == tmaxy:
                    rysize = ysize % tsize
                if rysize == 0:
                    rysize = tsize
                ry = ty * tsize

                (wx, wy) = (0, 0)
                (wxsize, wysize) = (int(rxsize / float(tsize)
                        * self.tilesize), int(rysize / float(tsize)
                        * self.tilesize))

                if self.options.verbose:
                    print ('\tReadRaster Extent: ', (rx, ry, rxsize,
                           rysize), (wx, wy, wxsize, wysize))

                # Query is in 'nearest neighbour' but can be bigger in then the tilesize
                # We scale down the query to the tilesize by supplied algorithm.

                # Tile dataset in memory

                dstile = self.mem_drv.Create('', self.tilesize,
                        self.tilesize, tilebands)
                data = ds.ReadRaster(
                    rx,
                    ry,
                    rxsize,
                    rysize,
                    wxsize,
                    wysize,
                    band_list=list(range(1, self.dataBandsCount + 1)),
                    )
                alpha = self.alphaband.ReadRaster(
                    rx,
                    ry,
                    rxsize,
                    rysize,
                    wxsize,
                    wysize,
                    )

                if self.tilesize == querysize:

                    # Use the ReadRaster result directly in tiles ('nearest neighbour' query)

                    dstile.WriteRaster(
                        wx,
                        wy,
                        wxsize,
                        wysize,
                        data,
                        band_list=list(range(1, self.dataBandsCount
                                + 1)),
                        )
                    dstile.WriteRaster(
                        wx,
                        wy,
                        wxsize,
                        wysize,
                        alpha,
                        band_list=[tilebands],
                        )
                else:

                    # Note: For source drivers based on WaveLet compression (JPEG2000, ECW, MrSID)
                    # the ReadRaster function returns high-quality raster (not ugly nearest neighbour)
                    # TODO: Use directly 'near' for WaveLet files
                    # Big ReadRaster query in memory scaled to the tilesize - all but 'near' algo

                    dsquery = self.mem_drv.Create('', querysize,
                            querysize, tilebands)

                    # TODO: fill the null value in case a tile without alpha is produced (now only png tiles are supported)
                    # for i in range(1, tilebands+1):
                    #   dsquery.GetRasterBand(1).Fill(tilenodata)

                    dsquery.WriteRaster(
                        wx,
                        wy,
                        wxsize,
                        wysize,
                        data,
                        band_list=list(range(1, self.dataBandsCount
                                + 1)),
                        )
                    dsquery.WriteRaster(
                        wx,
                        wy,
                        wxsize,
                        wysize,
                        alpha,
                        band_list=[tilebands],
                        )

                    self.scale_query_to_tile(dsquery, dstile,
                            tilefilename)
                    del dsquery

                del data

                if self.options.resampling != 'antialias':

                    # Write a copy of tile to png/jpg

                    self.out_drv.CreateCopy(tilefilename, dstile,
                            strict=0)

                del dstile

                if not self.options.verbose:
                    self.progressbar(ti / float(tcount))

    # -------------------------------------------------------------------------

    def generate_overview_tiles(self):
        """Generation of the overview tiles (higher in the pyramid) based on existing tiles"""

        print 'Generating Overview Tiles:'

        tilebands = self.dataBandsCount + 1

        # Usage of existing tiles: from 4 underlying tiles generate one as overview.

        tcount = 0
        for tz in range(self.tmaxz - 1, self.tminz - 1, -1):
            (tminx, tminy, tmaxx, tmaxy) = self.tminmax[tz]
            tcount += (1 + abs(tmaxx - tminx)) * (1 + abs(tmaxy
                    - tminy))

        ti = 0

        # querysize = tilesize * 2

        for tz in range(self.tmaxz - 1, self.tminz - 1, -1):
            (tminx, tminy, tmaxx, tmaxy) = self.tminmax[tz]
            yrange = range(tminy, tmaxy + 1)
            for ty in yrange:
                for tx in range(tminx, tmaxx + 1):

                    if self.stopped:
                        break

                    ti += 1
                    tilefilename = os.path.join(self.output, str(tz),
                            str(tx), '%s.%s' % (ty, self.tileext))

                    if self.options.verbose:
                        print (ti, '/', tcount, tilefilename)  # , "( TileMapService: z / x / y )"

                    if self.options.resume \
                        and os.path.exists(tilefilename):
                        if self.options.verbose:
                            print 'Tile generation skiped because of --resume'
                        else:
                            self.progressbar(ti / float(tcount))
                        continue

                    # Create directories for the tile

                    if not os.path.exists(os.path.dirname(tilefilename)):
                        os.makedirs(os.path.dirname(tilefilename))

                    dsquery = self.mem_drv.Create('', 2
                            * self.tilesize, 2 * self.tilesize,
                            tilebands)

                    # TODO: fill the null value
                    # for i in range(1, tilebands+1):
                    #   dsquery.GetRasterBand(1).Fill(tilenodata)

                    dstile = self.mem_drv.Create('', self.tilesize,
                            self.tilesize, tilebands)

                    # TODO: Implement more clever walking on the tiles with cache functionality
                    # probably walk should start with reading of four tiles from top left corner
                    # Hilbert curve

                    children = []

                    # Read the tiles and write them to query window

                    for y in range(2 * ty, 2 * ty + 2):
                        for x in range(2 * tx, 2 * tx + 2):
                            (minx, miny, maxx, maxy) = self.tminmax[tz
                                    + 1]
                            if x >= minx and x <= maxx and y >= miny \
                                and y <= maxy:
                                dsquerytile = \
                                    gdal.Open(os.path.join(self.output,
                                        str(tz + 1), str(x), '%s.%s'
                                        % (y, self.tileext)),
                                        gdal.GA_ReadOnly)

                                if ty:
                                    tileposy = y % (2 * ty) \
* self.tilesize
                                elif ty == 0 and y == 1:
                                    tileposy = self.tilesize
                                else:
                                    tileposy = 0

                                if tx:
                                    tileposx = x % (2 * tx) \
    * self.tilesize
                                elif tx == 0 and x == 1:
                                    tileposx = self.tilesize
                                else:
                                    tileposx = 0
                                dsquery.WriteRaster(
                                    tileposx,
                                    tileposy,
                                    self.tilesize,
                                    self.tilesize,
                                    dsquerytile.ReadRaster(0, 0,
        self.tilesize, self.tilesize),
                                    band_list=list(range(1, tilebands
        + 1)),
                                    )
                                children.append([x, y, tz + 1])

                    self.scale_query_to_tile(dsquery, dstile,
                            tilefilename)

                    # Write a copy of tile to png/jpg

                    if self.options.resampling != 'antialias':

                        # Write a copy of tile to png/jpg

                        self.out_drv.CreateCopy(tilefilename, dstile,
                                strict=0)

                    if self.options.verbose:
                        print (
                            '\tbuild from zoom',
                            tz + 1,
                            ' tiles:',
                            (2 * tx, 2 * ty),
                            (2 * tx + 1, 2 * ty),
                            (2 * tx, 2 * ty + 1),
                            (2 * tx + 1, 2 * ty + 1),
                            )

                    if not self.options.verbose:
                        self.progressbar(ti / float(tcount))

    # -------------------------------------------------------------------------

    def scale_query_to_tile(
        self,
        dsquery,
        dstile,
        tilefilename='',
        ):
        """Scales down query dataset to the tile dataset"""

        querysize = dsquery.RasterXSize
        tilesize = dstile.RasterXSize
        tilebands = dstile.RasterCount

        if self.options.resampling == 'average':

            # Function: gdal.RegenerateOverview()

            for i in range(1, tilebands + 1):

                # Black border around NODATA
                # if i != 4:
                #   dsquery.GetRasterBand(i).SetNoDataValue(0)

                res = gdal.RegenerateOverview(dsquery.GetRasterBand(i),
                        dstile.GetRasterBand(i), 'average')
                if res != 0:
                    self.error('RegenerateOverview() failed on %s, error %d'
                                % (tilefilename, res))
        elif self.options.resampling == 'antialias':

            # Scaling by PIL (Python Imaging Library) - improved Lanczos

            array = numpy.zeros((querysize, querysize, tilebands),
                                numpy.uint8)
            for i in range(tilebands):
                array[:, :, i] = \
                    gdalarray.BandReadAsArray(dsquery.GetRasterBand(i
                        + 1), 0, 0, querysize, querysize)
            im = Image.fromarray(array, 'RGBA')  # Always four bands
            im1 = im.resize((tilesize, tilesize), Image.ANTIALIAS)
            if os.path.exists(tilefilename):
                im0 = Image.open(tilefilename)
                im1 = Image.composite(im1, im0, im1)
            im1.save(tilefilename, self.tiledriver)
        else:

            # Other algorithms are implemented by gdal.ReprojectImage().

            dsquery.SetGeoTransform((
                0.0,
                tilesize / float(querysize),
                0.0,
                0.0,
                0.0,
                tilesize / float(querysize),
                ))
            dstile.SetGeoTransform((
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                ))

            res = gdal.ReprojectImage(dsquery, dstile, None, None,
                    self.resampling)
            if res != 0:
                self.error('ReprojectImage() failed on %s, error %d'
                           % (tilefilename, res))

    # -------------------------------------------------------------------------

    def generate_tilemapresource(self):
        """
        Template for tilemapresource.xml. Returns filled string. Expected variables:
          title, north, south, east, west, isepsg4326, projection, publishurl,
          zoompixels, tilesize, tileformat
        """

        args = {}
        args['title'] = self.options.title
        (args['south'], args['west'], args['north'], args['east']) = \
            self.swne
        args['tilesize'] = self.tilesize
        args['tileformat'] = self.tileext
        args['publishurl'] = self.options.url

        s = \
            """<?xml version="1.0" encoding="utf-8"?>
    <TileMap version="1.0.0" tilemapservice="http://tms.osgeo.org/1.0.0">
      <Title>%(title)s</Title>
      <Abstract></Abstract>
      <SRS></SRS>
      <BoundingBox minx="%(west).14f" miny="%(south).14f" maxx="%(east).14f" maxy="%(north).14f"/>
      <Origin x="%(west).14f" y="%(south).14f"/>
      <TileFormat width="%(tilesize)d" height="%(tilesize)d" mime-type="image/%(tileformat)s" extension="%(tileformat)s"/>
      <TileSets profile="raster">
""" \
            % args
        for z in range(self.tminz, self.tmaxz + 1):
            s += \
                """        <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" \
                % (args['publishurl'], z, 2 ** (self.nativezoom
                   - z) * self.out_gt[1], z)
        s += """      </TileSets>
    </TileMap>
    """
        return s

# =============================================================================
# =============================================================================
# =============================================================================

if __name__ == '__main__':
    argv = gdal.GeneralCmdLineProcessor(sys.argv)
    if argv:
        gdal2tiles = GDAL2Tiles(argv[1:])
        gdal2tiles.process()
