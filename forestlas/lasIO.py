import sys, os, struct, array, math, datetime, tempfile, shutil, subprocess
import multiprocessing
import numpy as np

### from forestlas.lasStructure import *
from lasStructure import *

class lasIO:

    """
    This allows reading and writing of .las files (currently supports 1.1 and 1.2).
    .las files can be read as one file(s), gridded or a plot of data extracted. 
    Output types include .las, a numpy array or as .znr file where only the height
    and the "number of returns" metadata remain. 
    """

    def __init__(self, path, out=False, znr=False, verbose=False, search=".",
                 tmpDir=False, keepTemp=False):

        """
        Functions creates holder for .las file and setups required dictionaries
        etc.
        
        Parameters
        ----------
        path: File path or list
            file path to tile .las file or directory containing .las files,
        also excepts a list of file paths.
        
        out: Path to directory or path to save .las to, Default None
            If None then saves output to os.getcwd()
        
        znr: Boolean. Default True
            Save output as a .znr file where only the height and "number of
        returns metadata is retained.
        
        search: String. Default "."
            Can be used to filter .las files if a directory is supplied for
        "path"
        
        tmpDir: File path to temporary directory. Default False
            Multiprocessing is not available with this module but is by others 
        that call lasIO.  This specifies the temporary directory.
        
        keepTemp: Boolean. Default True
            When False all temporary files are kept.  It is important to run 
            removeTemp() to clean up before exiting
        
        Returns
        -------
        out: self
         
        """

        self.verbose = verbose
        self.keepTemp = keepTemp

        ### parse args and create file structure
        # is path a directory or file
        if isinstance(path, list) :
            self.tileList = path
            self.dir = os.path.split(path[0])[0]
            self.tile = None
        elif os.path.isfile(path) :
            self.dir = os.path.split(path)[0]
            if self.dir == '': self.dir = os.getcwd()
            self.tile = os.path.split(path)[1]
            self.tileList = [os.path.join(self.dir, self.tile)]
        elif os.path.isdir(path) :
            self.dir = path
            self.tile = None
            tileList = os.listdir(self.dir)
            self.tileList = [tile for tile in tileList if (tile.endswith(".las")
                                                           or tile.endswith(".laz")
                                                           or tile.endswith(".znr"))
                                                           and tile.find(search) != -1]
        else:
            raise NameError ("path not recognised")

        if len(self.tileList) == 0:
            raise IOError("There are no .las, .laz or .znr tiles in {}".format(path))

        # create temporary directory in %temp%
        if tmpDir:
            self.tempDirectory = tmpDir
            if os.path.isdir(self.tempDirectory) is False:
                os.makedirs(self.tempDirectory)
        else:
            tempDirectoryName = "lidar.processing." + str(np.random.randint(0, 9999999)) + ".tmp"
            self.tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
            os.makedirs(self.tempDirectory)

        #tempDirectoryName = "lidar.processing." + str(np.random.randint(0, 9999999)) + ".tmp"
        #self.tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        #os.makedirs(self.tempDirectory)

        if self.verbose:    print "temporary directory at: {}".format(self.tempDirectory)

        # create output directory
        if not out:
            self.savePath = self.dir
            self.saveLAS = None
        elif os.path.isdir(out) :
            self.savePath = out
            self.saveLAS = None
        else:
            self.savePath = os.path.split(out)[0]
            self.saveLAS = os.path.split(out)[1]

        # znr
        if znr:
            self.znr = True
        else:   self.znr = False

        # global variables    
        self.numberTiles = 1
        for i, las in enumerate(self.tileList):
            las = os.path.join(self.dir, las)
            if i == 0:
                self.globalHeader = parseHeader(las)
                # self.globalHeader["numptrecords"] -= 1
                if self.globalHeader['filesig'] == "ZNRF":
                    self.ptFrmt, self.dt = znrStruct()
                    self.znr = True
                else:
                    self.ptFrmt, self.dt = getPtFrmt(self.globalHeader, verbose=True)
                self.vlr = getVLR(self.globalHeader["headersize"], las)
                self.h = parseHeader(las)
            else:
                lasHeader = parseHeader(os.path.join(path, las))
                self.globalHeader["numptrecords"] += lasHeader["numptrecords"]
                if lasHeader["xmax"] > self.globalHeader["xmax"]:
                    self.globalHeader["xmax"] = lasHeader["xmax"]
                if lasHeader["xmin"] < self.globalHeader["xmin"]:
                    self.globalHeader["xmin"] = lasHeader["xmin"]
                if lasHeader["ymax"] > self.globalHeader["ymax"]:
                    self.globalHeader["ymax"] = lasHeader["ymax"]
                if lasHeader["ymin"] < self.globalHeader["ymin"]:
                    self.globalHeader["ymin"] = lasHeader["ymin"]
                self.numberTiles += 1

        self.x_centre = np.mean([self.globalHeader["xmax"], self.globalHeader["xmin"]])
        self.y_centre = np.mean([self.globalHeader["ymax"], self.globalHeader["ymin"]])

        if self.verbose:    print "number of tiles to process: {}".format(self.numberTiles)
        self.area = (self.globalHeader["xmax"] - self.globalHeader["xmin"]) * (self.globalHeader["ymax"] - self.globalHeader["ymin"])
        # self.pointDensity = self.globalHeader["numptrecords"] / self.area

        # admin!
        self.counter = self.globalHeader["numptrecords"] // 20
        if self.globalHeader["numptrecords"] > 100000000:
            self.counter = self.globalHeader["numptrecords"] // 100
        self.badPoints = 0
        self.filesOpen = 0

    def laz2las(self, tile):

        """ 
        If input file is a .laz file then calls a LAStools to 
        decompress to .las. Requires LAStools to be installed and in 
        the PATH to be in the system environment variables
        """

        print "converting to .las: {}".format(os.path.abspath(tile))
        os.chdir(self.dir)
        tileName = os.path.splitext(os.path.split(tile)[1])[0] + ".las"
        tempFile = os.path.join(self.tempDirectory, tileName)
        cmd = ["las2las", "-i", tile, "-o", tempFile]
        subprocess.call(cmd)
        return tempFile

    def tiling(self, resolution, takeSample=False, buffer=0):

        """ 
        Used to tile .las file
        
        Parameters
        ----------
        resolution: int
            Tile resolution
        
        takeSample: int, Default False
            Samples data at specified factor
            
        Returns
        -------
        out: self
        
        """

        if self.globalHeader["filesig"] == "ZNRF":
            raise NameError("Can not tile a .znr file as there are no x y coordinates")

        self.totPoints = 0
        self.globalHeader["guid1"] = resolution

        self.xmax = (float(np.ceil(self.globalHeader["xmax"])) // resolution) * resolution # to make nice neat boxes!
        self.xmin = (float(np.floor(self.globalHeader["xmin"])) // resolution) * resolution
        self.ymax = (float(np.ceil(self.globalHeader["ymax"])) // resolution) * resolution
        self.ymin = (float(np.floor(self.globalHeader["ymin"])) // resolution) * resolution

        if self.tile == None:
            dirName = str(resolution) + "m_TILES"
        else:
            dirName = os.path.splitext(self.tile)[0] + "_" + str(resolution) + "m_TILES"
        self.tilePath = os.path.join(self.savePath, dirName)
        if os.path.isdir(self.tilePath) :
            shutil.rmtree(self.tilePath)
        os.makedirs(self.tilePath)

        # tile dictionary and other variables
        self.xy_dictionary = {}
        for i, x in enumerate(np.arange(self.xmin, self.xmax + resolution, resolution)):
            for j, y in enumerate(np.arange(self.ymin, self.ymax + resolution, resolution)):
                self.xy_dictionary[(x, y)] = {"xmin":x - buffer,
                                              "ymin":y - buffer,
                                              "zmin":999,
                                              "xmax":x + resolution + buffer,
                                              "ymax":y + resolution + buffer,
                                              "zmax": -999,
                                              "num_rtn": {1:0, 2:0, 3:0, 4:0, 5:0},
                                              "i":0,
                                              "tempFile": os.path.abspath(os.path.join(self.tempDirectory, str(x) + "_" + str(y) + ".temp")),
                                              "outFile": os.path.abspath(os.path.join(self.tilePath, str(x) + "_" + str(y) + ".las")),
                                              "isOpen": False,
                                              "lastTouched": None
                                             }

        keys = np.array(self.xy_dictionary.keys(), dtype=[('x', int), ('y', int)])

        if self.verbose:    print "number of plots: {}".format(len(self.xy_dictionary))

        for tile in self.tileList:

            tile = os.path.join(self.dir, tile)
            self.h = parseHeader(tile)

            if self.h["filesig"] == "LASF":
                if self.h["guid2"] == 1:
                    tile = self.laz2las(tile)

            self.h["num_rtn"] = {}

            if takeSample:
                sample = self.generateSample(takeSample)
                if self.verbose:
                    print "random sample produced: {}".format(len(sample))
            else:
                sample = range(self.h['numptrecords'])

            with open(os.path.join(self.dir, tile), 'rb') as fh:

                fh.seek(self.h["offset"])
                numPoints = 0

                for i in sample: # loops through all points

                    if i%self.counter == 0 and self.verbose:
                        print "{}% | {} of {} points selected | {}".format(np.round((np.float(self.totPoints)/self.globalHeader['numptrecords'])*100), numPoints, self.globalHeader['numptrecords'], datetime.datetime.now())

                    fh.seek(self.h["offset"] + (numPoints * self.h['pointreclen'])) # searches to beginning of point

                    point_dictionary = extract_return(self, fh)
                    X, Y = point_dictionary['x'], point_dictionary['y']
                    KEYS = keys[(X >= keys['x'] - buffer) & (X < keys['x'] + resolution + buffer) &
                                (Y >= keys['y'] - buffer) & (Y < keys['y'] + resolution + buffer)]

                    for key in KEYS:

                        self.xy_dictionary[tuple(key)] = self.writePoint(self.xy_dictionary[tuple(key)], point_dictionary)

                        if point_dictionary["rtn_num"] < 6 and point_dictionary["rtn_num"] > 0:
                            self.xy_dictionary[tuple(key)]["num_rtn"][point_dictionary["rtn_num"]] += 1

                    numPoints += 1
                    self.totPoints += 1

                    # closes file handles that have been unused for a while
                    if self.filesOpen == 40:
                        self.closeFiles()
                        self.filesOpen = 20

            # deletes .las tiles that were converted from .laz 
            if self.tempDirectory in tile:
                os.unlink(tile)

        for key in self.xy_dictionary.keys():
            if self.xy_dictionary[key]["isOpen"] is not False:
                self.xy_dictionary[key]["isOpen"].close()

        if self.verbose :    print "number of bad points = {}".format(self.badPoints)
        if self.znr :
            self.ptFrmt, self.dt = znrStruct()
            self.globalHeader["pointreclen"] = 5
        else:
            if resolution >= 1: self.h['guid1'] = resolution

        self.outTileCount = 0

        return self

    def all(self, takeSample=False):

        """ 
        Returns complete .las file
        
        Parameters
        ----------
        
        takeSample: int, Default False
            Samples data at specified factor
            
        Returns
        -------
        out: self
        
        """

        self.totPoints = 0
        self.globalHeader["guid1"] = 0

        x, y = self.globalHeader["xmin"], self.globalHeader["ymin"]
        if self.saveLAS == None:
            self.saveLAS = str(int(x)) + "_" + str(int(y)) + "_OUT.las"
        self.xy_dictionary = {}
        self.xy_dictionary['all'] = { "xmin":x, \
                                      "ymin":y, \
                                      "zmin":999, \
                                      "xmax":self.globalHeader["xmax"], \
                                      "ymax":self.globalHeader["ymax"], \
                                      "zmax": -999, \
                                      "num_rtn": {1:0, 2:0, 3:0, 4:0, 5:0}, \
                                      "i":0,    \
                                      "outFile": os.path.abspath(os.path.join(self.savePath, self.saveLAS)), \
                                      "tempFile": os.path.abspath(os.path.join(self.tempDirectory, str(x) + "_" + str(y) + ".temp")), \
                                      "isOpen": False, \
                                      "lastTouched": None
                                      }

        for tile in self.tileList:

            tile = os.path.join(self.dir, tile)
            self.h = parseHeader(tile)

            if self.h["filesig"] == "LASF":
                if self.h["guid2"] == 1:
                    tile = self.laz2las(tile)

            self.h["num_rtn"] = {}

            if takeSample:
                sample = self.generateSample(takeSample)
                if self.verbose :
                    print "random sample produced: {}".format(len(sample))
            else:
                sample = range(self.h['numptrecords'])

            with open(os.path.join(self.dir, tile), 'rb') as fh:

                fh.seek(self.h["offset"])
                numPoints = 0

                for i in sample: # loops through all points
                #for i in range(10): # number of rows

                    try:
                        if i%self.counter == 0 and self.verbose:
                            self.printProgress()
                    except: pass

                    fh.seek(self.h["offset"] + (numPoints * self.h['pointreclen'])) # searches to beginning of point

                    try:
                        point_dictionary = extract_return(self, fh)

                        if self.znr is False:
                            self.xy_dictionary["all"] = self.writePoint(self.xy_dictionary["all"], point_dictionary)
                            # populates number of returns
                            if point_dictionary["rtn_num"] < 6 and point_dictionary["rtn_num"] > 0:
                                self.xy_dictionary["all"]["num_rtn"][point_dictionary["rtn_num"]] += 1

                        else:
                            self.xy_dictionary["all"] = self.writeZNR(self.xy_dictionary["all"], point_dictionary)

                        numPoints += 1
                        self.totPoints += 1
                    except:
                        self.badPoints += 1
                        continue

        for key in self.xy_dictionary.keys():
            if self.xy_dictionary[key]["isOpen"] is not False:
                self.xy_dictionary[key]["isOpen"].close()

        if self.verbose :    print "number of bad points: {}".format(self.badPoints)
        if self.znr :
            self.ptFrmt, self.dt = znrStruct()
            self.globalHeader["pointreclen"] = 5

        return self

    def plot(self, centre_x, centre_y, extent=24, round=False):
        """
        
        Returns a plotwise array of points from the tile defined with 
        lasIO with plot centre at centre_x and centre_y and area equal 
        to (radius*2)**2 if round=False and pi*radius**2 if round=True.
        
        The radius defaults to 10 metres. Plots are square by default 
        but can be circular with round=True.
        
        Returns self which returns a numpy array if asArray is cploted, 
        or can saved as .las, .txt or .xyz by cploting exportLAS, 
        exportTXT or exportXYZ respectively.
        
        Paramaters
        ----------
        centre_x, centre_y : int or float
            Cartesian coordinates of plot centre (in the same coordinate 
        system as data)

        extent : int, float or tuple with length of 2
            Diameter of round plot or extent of square plot.  Will except
            a tuple of two ints or floats.
            
        round : Boolean, defualt False
            If False a square plot is returned, if True a round plot is
        returned
        
        Returns
        -------
        out : self
        
        """

        if isinstance(extent, tuple):
            extent_x = extent[0] / 2.
            extent_y = extent[1] / 2.
        else:
            extent_x, extent_y = extent / 2., extent / 2.

        xmin = centre_x - extent_x
        xmax = centre_x + extent_x
        ymin = centre_y - extent_y
        ymax = centre_y + extent_y

        self.totPoints = 0
        self.globalHeader["guid1"] = 0

        x, y = self.globalHeader["xmin"], self.globalHeader["ymin"]
        if self.saveLAS == None:
            self.saveLAS = str(int(x)) + "_" + str(int(y)) + "_OUT.las"
        self.xy_dictionary = {}
        self.xy_dictionary['plot'] = {"xmin":xmin,
                                      "ymin":ymin,
                                      "zmin":999,
                                      "xmax":xmax,
                                      "ymax":ymax,
                                      "zmax": -999,
                                      "num_rtn": {1:0, 2:0, 3:0, 4:0, 5:0},
                                      "i":0,
                                      "outFile": os.path.abspath(os.path.join(self.savePath, self.saveLAS)),
                                      "tempFile": os.path.abspath(os.path.join(self.tempDirectory, str(x) + "_" + str(y) + ".PLOT.temp")),
                                      "isOpen": False,
                                      "lastTouched": None
                                      }

        for tile in self.tileList:

            tile = os.path.join(self.dir, tile)
            self.h = parseHeader(tile)
            self.tile = tile

            if self.h["xmax"] < xmin or self.h["xmin"] > xmax or self.h["ymax"] < ymin or self.h["ymin"] > ymax:
                #print "skipping {}: out of bounds".format(tile)
                continue

            if self.h["filesig"] == "LASF":
                if self.h["guid2"] == 1:
                    tile = self.laz2las(tile)

            self.h["num_rtn"] = {}

            with open(tile, 'rb') as fh:

                fh.seek(self.h["offset"])

                for i in range(self.h['numptrecords']): # loops through plot points

                    if i%self.counter == 0 and self.verbose :
                        self.printProgress()

                    fh.seek(self.h["offset"] + (i * self.h['pointreclen'])) # searches to beginning of point

                    try:

                        # test x point first ...
                        fh.seek(self.h['offset'] + (i * self.h['pointreclen']))
                        x = fh.read(4)
                        x = struct.unpack('=' + 'L', x)[0]
                        x = (x * self.h['xscale'] ) + self.h['xoffset']
                        if x < xmin or x > xmax:
                            continue

                        # test y point next ...
                        fh.seek(self.h['offset'] + (i * self.h['pointreclen'] + 4))
                        y = fh.read(4)
                        y = struct.unpack('=' + 'L', y)[0]
                        y = (y * self.h['yscale']) + self.h['yoffset']
                        if y < ymin or y > ymax:
                            continue

                        # extract round plot
                        if round and round_plot((x, y), centre_x, centre_y, xmax, xmin, ymax, ymin, extent_x) == 0:
                            continue

                        fh.seek(self.h["offset"] + (i * self.h['pointreclen'])) # searches to beginning of point

                        if not self.znr:
                            point_dictionary = extract_return(self, fh)

                            self.xy_dictionary["plot"] = self.writePoint(self.xy_dictionary["plot"], point_dictionary)

                            # populates number of returns
                            if point_dictionary["rtn_num"] < 6 and point_dictionary["rtn_num"] > 0:
                                self.xy_dictionary["plot"]["num_rtn"][point_dictionary["rtn_num"]] += 1

                        else:
                            point_dictionary = extract_return(self, fh)
                            self.xy_dictionary["plot"] = self.writeZNR(self.xy_dictionary["plot"], point_dictionary)

                        self.totPoints += 1

                    except:
                        self.badPoints += 1
                        continue

        if self.xy_dictionary["plot"]["isOpen"]:
            self.xy_dictionary["plot"]["isOpen"].close()

        if self.verbose :    print "number of bad points = {}".format(self.badPoints)
        if self.znr :
            self.ptFrmt, self.dt = znrStruct()
            self.globalHeader["pointreclen"] = 5

        return self

    def fromGrid(self, csv, resolution=None):
        """

        Returns a plotwise array of points from the tile defined with
        lasIO with plot centre at centre_x and centre_y and area equal
        to (radius*2)**2 if round=False and pi*radius**2 if round=True.

        The radius defaults to 10 metres. Plots are square by default
        but can be circular with round=True.

        Returns self which returns a numpy array if asArray is cploted,
        or can saved as .las, .txt or .xyz by cploting exportLAS,
        exportTXT or exportXYZ respectively.

        Paramaters
        ----------
        centre_x, centre_y : int or float
            Cartesian coordinates of plot centre (in the same coordinate
        system as data)

        extent : int, float or tuple with length of 2
            Diameter of round plot or extent of square plot.  Will except
            a tuple of two ints or floats.

        round : Boolean, defualt False
            If False a square plot is returned, if True a round plot is
        returned

        Returns
        -------
        out : self

        """

        self.totPoints = 0
        self.globalHeader["guid1"] = 0

        grid = np.loadtxt(csv, skiprows=1, delimiter=',', dtype=([('x', np.float), ('y', np.float)]))
        grid['x'] = grid['x'].astype(int)
        grid['y'] = grid['y'].astype(int)

        if not resolution:
            resolution = grid['x'][1] - grid['x'][0]

        xmin = grid['x'].min() - (resolution / 2.)
        xmax = grid['x'].max() + (resolution / 2.)
        ymin = grid['y'].min() - (resolution / 2.)
        ymax = grid['y'].max() + (resolution / 2.)

        if self.verbose:
            print 'grid resolution:', resolution
            print 'aiming to produce {} tiles'.format(len(grid))

        self.xy_dictionary = {}
        for x, y in grid:
            self.xy_dictionary[(x, y)] = {"xmin":x - (resolution / 2.),
                                          "ymin":y - (resolution / 2.),
                                          "zmin":999,
                                          "xmax":x + (resolution / 2.),
                                          "ymax":y + (resolution / 2.),
                                          "zmax": -999,
                                          "num_rtn": {1:0, 2:0, 3:0, 4:0, 5:0},
                                          "i":0,
                                          "outFile": os.path.abspath(os.path.join(self.savePath, "{}_{}.PLOT.las".format(x, y))),
                                          "tempFile": os.path.abspath(os.path.join(self.tempDirectory, "{}_{}.PLOT.temp".format(x, y))),
                                          "isOpen": False,
                                          "lastTouched": None
                                          }

        keys = np.array(self.xy_dictionary.keys(), dtype=[('x', int), ('y', int)])

        for tile in self.tileList:

            tile = os.path.join(self.dir, tile)
            self.h = parseHeader(tile)
            self.tile = tile

            if self.h["xmax"] < xmin or self.h["xmin"] > xmax or self.h["ymax"] < ymin or self.h["ymin"] > ymax:
                continue

            if self.h["filesig"] == "LASF":
                if self.h["guid2"] == 1:
                    tile = self.laz2las(tile)

            self.h["num_rtn"] = {}

            with open(tile, 'rb') as fh:

                fh.seek(self.h["offset"])

                for i in range(self.h['numptrecords']): # loops through plot points

                    if i%self.counter == 0 and self.verbose :
                        self.printProgress()

                    fh.seek(self.h["offset"] + (i * self.h['pointreclen'])) # searches to beginning of point

                    # test x point first ...
                    fh.seek(self.h['offset'] + (i * self.h['pointreclen']))
                    x = fh.read(4)
                    x = struct.unpack('=' + 'L', x)[0]
                    x = (x * self.h['xscale'] ) + self.h['xoffset']
                    if not xmin < x < xmax:
                        continue

                    # test y point next ...
                    fh.seek(self.h['offset'] + (i * self.h['pointreclen'] + 4))
                    y = fh.read(4)
                    y = struct.unpack('=' + 'L', y)[0]
                    y = (y * self.h['yscale']) + self.h['yoffset']
                    if not ymin < y < ymax:
                        continue

                    # extract round plot
                    inGrid = False
                    idx = [(grid['x'] > self.h["xmin"]) | (grid['x'] <= self.h["xmax"]) |
                           (grid['y'] > self.h["ymin"]) | (grid['y'] <= self.h["ymax"])]
                    for row in grid[idx]:
                        if ((row[0] - (resolution / 2.)) < x < (row[0] + (resolution / 2.)) and
                            (row[1] - (resolution / 2.)) < y < (row[1] + (resolution / 2.))):
                            inGrid = True
                            break

                    if not inGrid: continue

                    fh.seek(self.h["offset"] + (i * self.h['pointreclen'])) # searches to beginning of point

                    # try:
                    point_dictionary = extract_return(self, fh)
                    X, Y = point_dictionary['x'], point_dictionary['y']
                    KEYS = keys[(X >= keys['x'] - (resolution / 2.)) & (X <= keys['x'] + (resolution / 2.)) &
                                (Y >= keys['y'] - (resolution / 2.)) & (Y <= keys['y'] + (resolution / 2.))]

                    for key in KEYS:

                        self.xy_dictionary[tuple(key)] = self.writePoint(self.xy_dictionary[tuple(key)], point_dictionary)

                        if point_dictionary["rtn_num"] < 6 and point_dictionary["rtn_num"] > 0:
                            self.xy_dictionary[tuple(key)]["num_rtn"][point_dictionary["rtn_num"]] += 1

                    self.totPoints += 1

                    # except:
                    #     self.badPoints += 1
                    #     continue

                    # closes file handles that have been unused for a while
                    if self.filesOpen == 40:
                        self.closeFiles()
                        self.filesOpen = 20

            # deletes .las tiles that were converted from .laz
            if self.tempDirectory in tile:
                os.unlink(tile)

        for key in self.xy_dictionary.keys():
            if self.xy_dictionary[key]["isOpen"] is not False:
                self.xy_dictionary[key]["isOpen"].close()

        if self.verbose :    print "number of bad points = {}".format(self.badPoints)
        if self.znr :
            self.ptFrmt, self.dt = znrStruct()
            self.globalHeader["pointreclen"] = 5
        else:
            if resolution >= 1: self.h['guid1'] = resolution

        self.outTileCount = 0

        return self

    def writePoint(self, tile, d):

        # opens file for writing points to

        if not type(tile["isOpen"]) == 'file':
            tile["isOpen"] = open(tile["tempFile"], "ab")
            self.filesOpen += 1

        # writes point
        for i in self.ptFrmt:
            if i[0] == 'return_grp':
                byte = ((d['scan_edge'] & 1) << 7) | ((d['scan_dir'] & 1) << 6) | ((d['rtn_tot'] & 7) << 3) | (d['rtn_num'] & 7)
            elif i[0] == 'x':
                byte = (d['x'] - self.h['xoffset']) / self.h['xscale']
            elif i[0] == 'y':
                byte = (d['y'] - self.h['yoffset']) / self.h['yscale']
            elif i[0] == 'z':
                byte = (d['z'] - self.h['zoffset']) / self.h['zscale']
            else:
                byte = d[i[0]]

            tile["isOpen"].write(struct.pack('=' + i[2], byte))

        # updates header information
        if d['z'] > tile["zmax"]:
            tile["zmax"] = d['z']
        if d['z'] < tile["zmin"]:
            tile["zmin"] = d['z']
        tile["i"] += 1

        tile["lastTouched"] = datetime.datetime.now()

        return tile

    def writeZNR(self, tile, d):

        # opens file for writing points to 
        if tile["isOpen"] == False:
            tile["isOpen"] = open(tile["tempFile"], "ab")
            self.filesOpen += 1

        # writes point
        for i in znrStruct()[0]: # retrieves struct from lasStructure
            if i[0] == 'z':
                byte = (d['z'] / .01 )
            else:
                byte = d[i[0]]
            tile["isOpen"].write(struct.pack('=' + i[2], byte))

        # updates header information
        if d['z'] > tile["zmax"]:
            tile["zmax"] = d['z']
        if d['z'] < tile["zmin"]:
            tile["zmin"] = d['z']
        tile["i"] += 1

        tile["lastTouched"] = datetime.datetime.now()

        return tile

    def asDic(self):

        """
        Returns array as a dictionary where the keys are
        the central coordinates and the values are a list of
        tuples (height, number of returns). This is useful
        for plotting or calculating continuous variables across
        a plot.
        """

        arr = {}
        for i, key in enumerate(self.xy_dictionary.keys()):
            tile = self.xy_dictionary[key]
            if tile["i"] == 0:
                del tile
                continue
            a = np.zeros(tile["i"], dtype = self.dt)
            with open(tile["tempFile"], "rb") as fh:
                for j, line in enumerate(range(tile["i"])):
                    fh.seek(j * self.globalHeader['pointreclen'])
                    d = extract_return(self, fh)
                    for field in d:
                        a[j][field] = d[field]
            arr[(tile["xmin"], tile["ymin"])] = a
            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

        return arr

    def asArray(self):

        """
        Returns all points in a single array or as a list of
        arrays if there is more than one tile.  If there is 
        more than one tile then it may be preferable to use
        asDic command.
        """

        arr = []
        for i, key in enumerate(self.xy_dictionary.keys()):
            tile = self.xy_dictionary[key]
            if tile["i"] == 0:
                del tile
                continue
            if self.znr:
                a = np.zeros(tile["i"], dtype = znrStruct()[1])
            else:
                a = np.zeros(tile["i"], dtype = self.dt)
            with open(tile["tempFile"], "rb") as fh:
                for j, line in enumerate(range(tile["i"])):
                    fh.seek(j * self.globalHeader['pointreclen'])
                    d = extract_return(self, fh)
                    for field in d:
                        a[j][field] = d[field]
            arr.append(a)
            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

        if len(arr) == 1: # if list arr has only one array...
            arr = arr[0]  # ...return only the array

        return arr

    def exportLAS(self, out=False):

        """
        Writes tile(s) to .las format. Requires znr=False or
        will throw an exception
        """

        if self.znr :    raise NameError("Processed as ZNR, can not export to .LAS")
        if len(self.xy_dictionary) > 10:
            nTileCounter = len(self.xy_dictionary) // 10
        else:   nTileCounter = len(self.xy_dictionary)

        tileCount = 0

        for i, key in enumerate(self.xy_dictionary.keys()):
            if self.xy_dictionary[key]["i"] == 0:
                del self.xy_dictionary[key]
                continue
            else:
                tileCount += 1

            if i%nTileCounter == 0 and i > 0 and self.verbose:
                print "{:.0f}% | {} of {} tiles exported | {}".format( np.float(i) / len(self.xy_dictionary) * 100, i, len(self.xy_dictionary), datetime.datetime.now())

            tile = self.xy_dictionary[key]
            self.h['gensoftware'] = 'CRC207 LiDAR analysis software  '
            self.h['sysid'] = 'CRC207 LiDAR analysis software  '
            self.h['xmin'] = tile['xmin']
            self.h['xmax'] = tile['xmax']
            self.h['ymin'] = tile['ymin']
            self.h['ymax'] = tile['ymax']
            self.h['zmin'] = tile['zmin']
            self.h['zmax'] = tile['zmax']
            self.h['numptbyreturn'] = tuple([tile["num_rtn"][i] for i in range(1, 6)]) # sorting out the rtn_num tuple
            self.h['numptrecords'] = tile["i"]
            self.h['guid2'] = 0
            if self.h['numvlrecords'] > 0:
                self.h['offset'] = 313
                self.h['numvlrecords'] = 1

            if out:
                if out.endswith('.las') and len(self.xy_dictionary) > 1:
                    outFile = out[:-4] + '.' + str(self.outTileCount) + '.las'
                elif out.endswith('.las') and os.path.isdir(os.path.split(out)[0]):
                    outFile = out
                elif os.path.isdir(out):
                    outFile = os.path.join(out, os.path.split(tile["outFile"])[1])
                else:
                    raise IOError('out path not recognised')
            else:
                outFile = tile["outFile"]
            self.savePath = os.path.split(outFile)[0]

            with open(outFile, "wb") as outOpen:
                for j in headerstruct():
                    if j[2] == 'c':
                        outOpen.write(self.h[j[0]])
                    elif j[3] > 1:
                        outOpen.write(struct.pack('=' + str(j[3]) + j[2], *self.h[j[0]]))
                    else:
                        outOpen.write(struct.pack('=' + j[2] , self.h[j[0]]))

                ## write VLR
                if self.h['numvlrecords'] > 0:
                    # keeps only the first VLR e.g. the projection data
                    outOpen.write(self.vlr[:86])

                ## write points
                outOpen.seek(self.h['offset'])
                with open(tile["tempFile"], "rb") as o:
                    points = o.read()
                outOpen.write(points)

            tile["isOpen"] = "written_to"

            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if tileCount > 0:
            print "100% | {} of {} tiles exported | {}".format(len(self.xy_dictionary),
                                                               len(self.xy_dictionary),
                                                               datetime.datetime.now())

            if len(self.xy_dictionary) == 1:
                print ".las file written to {}".format(outFile)
            else:
                print ".las file(s) written to {}".format(os.path.split(outFile)[0])
        else:
            print "! no tiles to export !"

        if not self.keepTemp:  shutil.rmtree(self.tempDirectory)

        return self

    def np2ZNR(self, arr):

        """ 
        Can be used to export a numpy array in .znr format
        """

        if self.znr == False:   raise NameError ("Data processed as .LAS and therefore wcan not be saved as .ZNR")

        self.ptFrmt, self.dt = znrStruct()

        if len(self.xy_dictionary) > 10:
            nTileCounter = len(self.xy_dictionary) // 10
        else:   nTileCounter = len(self.xy_dictionary)

        for i, key in enumerate(self.xy_dictionary.keys()):
            if self.xy_dictionary[key]["i"] == 0:
                del self.xy_dictionary[key]
                continue

            if i%nTileCounter == 0:
                print "{}% | {} of {} tiles exported | {}".format(np.round((i/np.float(len(self.xy_dictionary)))*100), i, len(self.xy_dictionary), datetime.datetime.now())

            tile = self.xy_dictionary[key]
            self.h['filesig'] = "ZNRF"
            self.h['gensoftware'] = 'CRC207 LiDAR analysis software  '
            self.h['sysid'] = 'CRC207 LiDAR analysis software  '
            self.h['xmin'] = tile['xmin']
            self.h['xmax'] = tile['xmax']
            self.h['ymin'] = tile['ymin']
            self.h['ymax'] = tile['ymax']
            self.h['zmin'] = tile['zmin']
            self.h['zmax'] = tile['zmax']
            self.h['numptrecords'] = tile["i"]
            self.h['resolution'] = self.globalHeader['guid1']
            self.h['pointreclen'] = 5
            self.h['headersize'] = 92
            self.h['offset'] = 178
            self.h['zscale'] = 0.01
            self.h['zoffset'] = 0

            outFile = os.path.splitext(tile["outFile"])[0] + ".znr"
            with open(outFile, "wb") as out:
                for i in znrHeaderStruct():
                    if i[2] == 'c':
                        out.write(self.h[i[0]])
                    else:
                        out.write(struct.pack('=' + i[2] , self.h[i[0]]))

                ## write VLR
                out.write(self.vlr)

                ## write points
                for line in arr:
                    for i in znrStruct()[0]: # retrieves struct from lasStructure
                        if i[0] == 'z':
                            byte = line[0] / .01
                        else:
                            byte = line[1]
                        out.write(struct.pack('=' + i[2], byte))

            tile["isOpen"] = "written_to"
            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if self.verbose :    print ".znr file written to {}".format(os.path.split(outFile)[0])

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

    def np2LAS(self, arr, out=False):

        """ 
        Can be used to export a numpy array in .las format
        """

        if self.znr:   raise NameError ("Data processed as .ZNR and therefore wcan not be saved as .LAS")

        self.h['gensoftware'] = 'CRC207 LiDAR analysis software  '
        self.h['sysid'] = 'CRC207 LiDAR analysis software  '
        self.h['xmin'] = arr['x'].min()
        self.h['xmax'] = arr['x'].max()
        self.h['ymin'] = arr['y'].min()
        self.h['ymax'] = arr['y'].max()
        self.h['zmin'] = arr['z'].min()
        self.h['zmax'] = arr['z'].max()
        ### sorting out the rtn_num tuple
        rtn = np.zeros(5).astype(int)
        for row in arr:
            rtn_num = row['rtn_num'] - 1
            if rtn_num < 5:
                rtn[rtn_num] += 1
        self.h['numptbyreturn'] = tuple(rtn)
        self.h['numptrecords'] = len(arr)
        self.h['guid2'] = 0

        if out:
            if out.endswith('.las'):
                outFile = out
            elif os.path.isdir(out):
                outFile = os.path.join(out, os.path.split(tile["outFile"])[1])
            else:
                raise IOError('out path not recognised')
        else:
            x = str(np.mean(arr['x']).astype(int))
            y = str(np.mean(arr['y']).astype(int))
            outFile = os.path.join(self.savePath, "{}_{}.las".format(x, y))

        with open(outFile, "wb") as out:
            for j in headerstruct():
                if j[2] == 'c':
                    out.write(self.h[j[0]])
                elif j[3] > 1:
                    out.write(struct.pack('=' + str(j[3]) + j[2], *self.h[j[0]]))
                else:
                    out.write(struct.pack('=' + j[2] , self.h[j[0]]))

            ## write VLR
            if self.h['numvlrecords'] > 0:
                out.write(self.vlr)

            ## write points
            out.seek(self.h['offset'])
            for d in arr:
                for i in self.ptFrmt:
                    if i[0] == 'return_grp':
                        byte = ((d['scan_edge'] & 1) << 7) | ((d['scan_dir'] & 1) << 6) | ((d['rtn_tot'] & 7) << 3) | (d['rtn_num'] & 7)
                    elif i[0] == 'x':
                        byte = (d['x'] - self.h['xoffset']) / self.h['xscale']
                    elif i[0] == 'y':
                        byte = (d['y'] - self.h['yoffset']) / self.h['yscale']
                    elif i[0] == 'z':
                        byte = (d['z'] - self.h['zoffset']) / self.h['zscale']
                    else:
                        byte = d[i[0]]
                    out.write(struct.pack('=' + i[2], byte))

        if self.verbose: print ".las file written to {}".format(outFile)

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

        return outFile

    def exportZNR(self):

        if self.znr == False:   raise NameError ("Data processed as .LAS and therefore cannot be saved as .ZNR")

        self.ptFrmt, self.dt = znrStruct()
        self.globalHeader["pointreclen"] = 5
        if len(self.xy_dictionary) > 10:
            nTileCounter = len(self.xy_dictionary) // 10
        else:   nTileCounter = len(self.xy_dictionary)

        for i, key in enumerate(self.xy_dictionary.keys()):
            if self.xy_dictionary[key]["i"] == 0:
                del self.xy_dictionary[key]
                continue

            print len(self.xy_dictionary)
            if i%nTileCounter == 0:
                print "{}% | {} of {} tiles exported | {}".format(np.round((i/np.float(len(self.xy_dictionary)))*100), i, len(self.xy_dictionary), datetime.datetime.now())

            tile = self.xy_dictionary[key]
            self.h['filesig'] = "ZNRF"
            self.h['gensoftware'] = 'CRC207 LiDAR analysis software  '
            self.h['sysid'] = 'CRC207 LiDAR analysis software  '
            self.h['xmin'] = tile['xmin']
            self.h['xmax'] = tile['xmax']
            self.h['ymin'] = tile['ymin']
            self.h['ymax'] = tile['ymax']
            self.h['zmin'] = tile['zmin']
            self.h['zmax'] = tile['zmax']
            self.h['numptrecords'] = tile["i"]
            self.h['resolution'] = self.globalHeader['guid1']
            self.h['pointreclen'] = 5
            self.h['headersize'] = 92
            self.h['offset'] = 178
            self.h['zscale'] = 0.01
            self.h['zoffset'] = 0

            outFile = os.path.splitext(tile["outFile"])[0] + ".znr"
            with open(outFile, "wb") as out:
                for i in znrHeaderStruct():
                    if i[2] == 'c':
                        out.write(self.h[i[0]])
                    else:
                        out.write(struct.pack('=' + i[2] , self.h[i[0]]))

                ## write VLR
                out.write(self.vlr)

                ## write points
                with open(tile["tempFile"], "rb") as o:
                    points = o.read()
                out.write(points)

            tile["isOpen"] = "written_to"
            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if len(self.xy_dictionary) == 0:
            print "No tiles written"
        elif len(self.xy_dictionary) == 1:
            print ".znr file written to {}".format(tile["outFile"])
        else:
            if self.verbose :    print ".znr files written to {}".format(os.path.split(outFile)[0])
            if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

    def LAS2txt(self, enum=False, outFile=False):

        for i, key in enumerate(self.xy_dictionary.keys()):
            tile = self.xy_dictionary[key]
            if tile["i"] == 0:
                del tile
                continue
            if enum:
                savePath = os.path.join(os.path.split(tile["outFile"])[0], "{}.txt".format(i))
            elif outFile:
                if os.path.isdir(outFile):
                    savePath = os.path.join(outFile, "{}.txt".format(i))
                else:
                    savePath = outFile
            else: savePath = os.path.splitext(tile["outFile"])[0] + ".txt"
            a = np.zeros(tile["i"], dtype = self.dt)
            with open(tile["tempFile"], "rb") as fh:
                for j, line in enumerate(range(tile["i"])):
                    fh.seek(j * self.globalHeader['pointreclen'])
                    d = extract_return(self, fh)
                    for field in d:
                        a[j][field] = d[field]
            np.savetxt(savePath, np.transpose([a['x'], \
                                               a['y'], \
                                               a['z'], \
                                               a['class'], \
                                               a['i'], \
                                               a['scan_ang'], \
                                               a['rtn_num'], \
                                               a['rtn_tot']]), \
                                               fmt='%.1f', delimiter=',')
            if self.verbose: print '.txt saved to:', savePath
            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

    def exportXYZ(self):

        for i, key in enumerate(self.xy_dictionary.keys()):
            tile = self.xy_dictionary[key]
            if tile["i"] == 0:
                del tile
                continue
            self.savePath = os.path.splitext(tile["outFile"])[0] + ".txt"
            a = np.zeros(tile["i"], dtype = self.dt)
            with open(tile["tempFile"], "rb") as fh:
                for j, line in enumerate(range(tile["i"])):
                    fh.seek(j * self.globalHeader['pointreclen'])
                    d = extract_return(self, fh)
                    for field in d:
                        a[j][field] = d[field]
            np.savetxt(self.savePath, np.transpose([a['x'], a['y'], a['z']]), fmt='%.1f', delimiter=',')

            if self.keepTemp is False:  os.unlink(tile["tempFile"])

        if self.keepTemp is False:  shutil.rmtree(self.tempDirectory)

        if self.verbose: print 'XYZ saved to:', self.savePath

        return self

    def closeFiles(self):

        l = {}
        for key in self.xy_dictionary.keys():
            if self.xy_dictionary[key]["lastTouched"] is not None:
                l[key] = self.xy_dictionary[key]["lastTouched"]
        for i, (key, val) in enumerate(zip(l.keys(), sorted(l.values()))):
            if self.xy_dictionary[key]["isOpen"] != False:
                self.xy_dictionary[key]["isOpen"].close()
                self.xy_dictionary[key]["isOpen"] = False
                if i == 20:
                    return

    def printProgress(self):

        print "{}% | {} of {} points processed | {}".format(np.round((np.float(self.totPoints)/self.globalHeader['numptrecords'])*100), self.totPoints, self.globalHeader['numptrecords'], datetime.datetime.now())

    def generateSample(self, sample):

        return np.random.choice(self.h['numptrecords'],
                                size=int(self.h['numptrecords'] / sample),
                                replace=False)

    def removeTemp(self):

        if os.path.isdir(self.tempDirectory):
            for file in os.listdir(self.tempDirectory):
                os.unlink(os.path.join(self.tempDirectory, file))

        shutil.rmtree(self.tempDirectory)
        if self.verbose: print "{} has been deleted".format(self.tempDirectory)


def parseHeader(filename):

    """
    returns header information as a dictionary.
    """

    with open(filename,'rb') as fh:
        header = {'infile':filename}

        if fh.read(4) == "ZNRF":
            headerStructType = znrHeaderStruct()
        else:   headerStructType = headerstruct()

        fh.seek(0)

        for i in headerStructType:
            if i[2] == 'c':
                value = fh.read(i[1])
            elif i[3] > 1:
                value = struct.unpack( '=' + str(i[3]) + i[2] , fh.read(i[1]) )
            else:
                value = struct.unpack( '=' + i[2] , fh.read(i[1]) )[0]
            header[i[0]] = value

    if headerStructType == headerstruct():
        if header["pointformat"] > 127: # it is a .laz file
            header["pointformat"] -= 128
            header["numvlrecords"] -= 1
            header["offset"] -= 100
            header["guid2"] = 1

    return header


def extract_return(self, fh):

    if self.znr:
        ptFrmt = znrStruct()[0]
    else: ptFrmt = self.ptFrmt

    point_dictionary = {} # dictionary for storing values in for each point

    for ent in ptFrmt:
        byte = fh.read(ent[1])
        val = struct.unpack('=' + ent[2] , byte)[0]
        if ent[0] == 'x':
            val = (val * self.h['xscale'] ) + self.h['xoffset']
        if ent[0] == 'y':
            val = (val * self.h['yscale'] ) + self.h['yoffset']
        if ent[0] == 'z':
            val = (val * self.h['zscale'] ) + self.h['zoffset']
        if ent[0] == 'return_grp':
            point_dictionary['rtn_num'] = val & 7
            #if point_dictionary['rtn_num'] == 0:
            #    raise Exception
            point_dictionary['rtn_tot'] = (val >> 3) & 7
            point_dictionary['scan_dir'] = (val >> 6) & 1
            point_dictionary['scan_edge'] = (val >> 7)
            continue # required so that 'return_grp' is not added to dictionary

        point_dictionary[ent[0]] = val

    if point_dictionary["z"] > 1000000:
        raise NameError("z very high: {}".format(point_dictionary["z"]))

    return point_dictionary


def round_plot(point, plot_x, plot_y, xmax, xmin, ymax, ymin, r):

    x, y = point

    #NW
    if plot_x <= x and plot_y >= y:
        adj = x - plot_x
        opp = plot_y - y
        hyp = math.sqrt(math.pow(adj, 2) + math.pow(opp, 2))

    #SW
    elif plot_x <= x and plot_y <= y:
        adj = x - plot_x
        opp = y - plot_y
        hyp = math.sqrt(math.pow(adj, 2) + math.pow(opp, 2))

    #SE
    elif plot_x >= x and plot_y <= y:
        adj = plot_x - x
        opp = y - plot_y
        hyp = math.sqrt(math.pow(adj, 2) + math.pow(opp, 2))

    #NE
    elif plot_x >= x and plot_y >= y:
        adj = plot_x - x
        opp = plot_y - y
        hyp = math.sqrt(math.pow(adj, 2) + math.pow(opp, 2))

    if hyp > r: return 0


def getPtFrmt(globalHeader, verbose=False):

    # return structure
    if "txt2las" in globalHeader["gensoftware"]:
        ptFrmt, dt = point_fmtLTstruct()
    elif globalHeader["pointformat"] == 0:
        ptFrmt, dt = point_fmt0struct()
    elif globalHeader["pointformat"] == 1:
        ptFrmt, dt = point_fmt1struct()
    elif globalHeader["pointformat"] == 3:
        ptFrmt, dt = point_fmt3struct()

    return ptFrmt, dt


def getVLR(headerSize, las):
        fh = open(os.path.join(las))
        fh.seek(headerSize)
        vlr = fh.read(86)
        fh.close()
        return vlr


if __name__ == '__main__':

    # import glob
    #
    # A = []
    # for las in glob.glob('/Users/phil/ALS/RiverHealth/Las/originalLas/*.las'):
    #     name = os.path.split(las)[1]
    #     x, y = int(name[1:4]), int(name[5:9])
    # #     plt.scatter(x, y)
    #     if y > 5950: A.append(las)

    # T = lasIO('/Users/phil/ALS/WC/spl/tile_20/WC1_10m_TILES', verbose=True).tiling(5).exportLAS('/Users/phil/ALS/WC/spl/tile_20/WC1_5m_TILES')
    T = lasIO('/Users/phil/ALS/WC/spl/tile_20/ForestLAS_tutorial/LAS/large_tile/WC45_SUB.las', verbose=True, keepTemp=True).tiling(20).exportLAS()
    # T = lasIO('/Users/phil/ALS/WC/spl/tile_20/WC1_test', verbose=True).all().exportLAS('/Users/phil/ALS/WC/spl/tile_20/WC1_test/all.las')
    # T = lasIO('/Users/phil/ALS/RiverHealth/Las/originalLas', verbose=True).fromGrid('/Users/phil/Dropbox/regular_points_3_30_test.csv').exportLAS('/Users/phil/ALS/WC/spl/tile_20/WC1_grid/works')

    # lasIO('/Users/phil/ALS/WC/spl/tile_20/WC1.las', verbose=True, keepTemp=True).all().exportLAS()