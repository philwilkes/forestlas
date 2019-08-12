import os
import struct
import math
import datetime
import tempfile
import shutil
import subprocess
import multiprocessing
import random
import numpy as np
from lasStructure import *

class lasIO:

    """
    This allows reading and writing of .las files (currently supports 1.1 and 1.2).
    .las files can be read as one file(s), gridded or a plot of data extracted. 
    Output types include .las, a numpy array or as .znr file where only the height
    and the "number of returns" metadata remain. 
    """

    def __init__(self, path, out=False, verbose=False, search=".",
                 tmpDir=False, keepTemp=False, number_of_processes=1,
                 copy=False, create_temp=True):

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

        search: String. Default "."
            Can be used to filter .las files if a directory is supplied for
            "path"
        
        tmpDir: File path to temporary directory. Default False
            Multiprocessing is not available with this module but is by others 
            that call lasIO.  This specifies the temporary directory.
        
        keepTemp: Boolean. Default True
            When False all temporary files are kept.  It is important to run 
            removeTemp() to clean up before exiting

        number_of_processes: Int, Default 1
            If processing a number of tiles or using a regular grid, Python's
            multiprocessing canbe envoked.  This variable sets the number of
            cores that are utilised.

        copy: Boolean, Default False
            Wether to copy files to the temp directory before processing.  This
            can speed thinks up if reading a number of large files simultaneously
            from a remote drive.
        
        Returns
        -------
        out: self
         
        """

        self.verbose = verbose
        self.keepTemp = keepTemp
        self.number_of_processes = number_of_processes

        ### parse args and create file structure
        # is path a directory or file
        if isinstance(path, list) :
            self.tile_list = path
            self.dir = os.path.split(path[0])[0]
            self.tile = None
        elif os.path.isfile(path) :
            self.dir = os.path.split(path)[0]
            if self.dir == '': self.dir = os.getcwd()
            self.tile = os.path.split(path)[1]
            self.tile_list = [os.path.join(self.dir, self.tile)]
        elif os.path.isdir(path) :
            self.dir = path
            self.tile = None
            tile_list = os.listdir(self.dir)
            self.tile_list = [tile for tile in tile_list if (tile.endswith(".las")
                                                           or tile.endswith(".laz"))
                                                           and tile.find(search) != -1]
        else:
            raise NameError ("path not recognised")

        if len(self.tile_list) == 0:
            raise IOError("There are no .las or .laz tiles in {}".format(path))

        # create temporary directory in %temp%
        if create_temp:
            if tmpDir:
                self.tempDirectory = tmpDir
                if os.path.isdir(self.tempDirectory) is False:
                    os.makedirs(self.tempDirectory)
            else:
                tempDirectoryName = "lidar.processing." + str(np.random.randint(0, 9999999)) + ".tmp"
                self.tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
                os.makedirs(self.tempDirectory)

        if self.verbose: print "temporary directory at: {}".format(self.tempDirectory)

        # create output directory
        if not out:
            self.savePath = self.dir
            self.saveLAS = None
        elif os.path.isdir(out):
            self.savePath = out
            self.saveLAS = None
        else:
            self.savePath = os.path.split(out)[0]
            self.saveLAS = os.path.split(out)[1]

        # global variables    
        self.numberTiles = 1
        for i, las in enumerate(self.tile_list):
            las = os.path.join(self.dir, las)
            if i == 0:
                self.globalHeader = parseHeader(las)
                # self.globalHeader["numptrecords"] = lasHeader["numptrecords"]
                if self.globalHeader['filesig'] == "ZNRF":
                    raise Exception('ZNR is deprecated, use an older version of lasIO')
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

        if self.verbose:    print "number of tiles to process:", self.numberTiles
        if self.numberTiles > 1 and self.verbose:    print "processing tiles from:", self.dir
        # self.area = (self.globalHeader["xmax"] - self.globalHeader["xmin"]) * (self.globalHeader["ymax"] - self.globalHeader["ymin"])
        # self.pointDensity = self.globalHeader["numptrecords"] / self.area

        # admin!
        self.counter = self.globalHeader["numptrecords"] // 20
        if self.globalHeader["numptrecords"] > 1e7:
            self.counter = self.globalHeader["numptrecords"] // 100
        self.badPoints = 0
        self.copy = copy
        if self.copy and self.verbose:
            print 'tiles are copied to temp'
        self.resolution = None

    def grid(self, csv, resolution=None, take_sample=False):

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

        self.take_sample = take_sample

        self.Global_totPoints = 0
        self.Global_files_open = 0
        self.globalHeader["guid1"] = 0

        self.grid = np.loadtxt(csv, skiprows=1, delimiter=',', dtype=([('x', np.float), ('y', np.float)]))
        self.grid['x'] = self.grid['x'].astype(int)
        self.grid['y'] = self.grid['y'].astype(int)

        if not resolution:
            self.resolution = self.grid['x'][1] - self.grid['x'][0]
        else:
            self.resolution = resolution

        self.xmin = self.grid['x'].min() - (self.resolution / 2.)
        self.xmax = self.grid['x'].max() + (self.resolution / 2.)
        self.ymin = self.grid['y'].min() - (self.resolution / 2.)
        self.ymax = self.grid['y'].max() + (self.resolution / 2.)

        if self.tile == None:
            dirName = str(self.resolution) + "m_GRID"
        else:
            dirName = os.path.splitext(self.tile)[0] + "_" + str(self.resolution) + "m_GRID"
        self.savePath = os.path.join(self.savePath, dirName)
        if os.path.isdir(self.savePath):
            shutil.rmtree(self.savePath)

        if self.verbose:
            print 'grid resolution:', self.resolution
            print 'aiming to produce {} tiles'.format(len(self.grid))

        self.xy_dictionary = dict()
        for x, y in self.grid:
            self.xy_dictionary[(x, y)] = {"xmin":x - (self.resolution / 2.),
                                          "ymin":y - (self.resolution / 2.),
                                          "zmin":999,
                                          "xmax":x + (self.resolution / 2.),
                                          "ymax":y + (self.resolution / 2.),
                                          "zmax": -999,
                                          "num_rtn": {1:0, 2:0, 3:0, 4:0, 5:0},
                                          "i":0,
                                          "outFile": os.path.abspath(os.path.join(self.savePath, "{}_{}.PLOT.las".format(x, y))),
                                          "tempFile": os.path.abspath(os.path.join(self.tempDirectory, "{}_{}.PLOT.temp".format(x, y))),
                                          "isOpen": False,
                                          "lastTouched": None

                                          }

        self.keys = np.array(self.xy_dictionary.keys(), dtype=[('x', int), ('y', int)])

        for tile in self.tile_list:

            if len(os.path.split(tile)) > 1:
                    tile = os.path.join(self.dir, tile)

            if self.copy:
                shutil.copyfile(tile, os.path.join(self.tempDirectory, os.path.split(tile)[1]))
                tile = os.path.join(self.tempDirectory, os.path.split(tile)[1])
            h = parseHeader(tile)

            if h["filesig"] == "LASF":
                if h["guid2"] == 1:
                    tile = self.laz2las(tile)

            if self.take_sample:
                sample = self.generateSample(self.take_sample)
                if self.verbose:
                    print "random sample produced: {}".format(len(sample))
            else:
                sample = range(h['numptrecords'])

            grid = self.grid[(self.grid['x'] >= h['xmin']) & (self.grid['x'] <= h['xmax']) &
                             (self.grid['y'] >= h['ymin']) & (self.grid['y'] <= h['ymax'])]

            with open(os.path.join(self.dir, tile), 'rb') as fh:

                fh.seek(h["offset"])

                for i in sample: # loops through all points

                    if self.Global_totPoints%self.counter == 0 and self.verbose:
                        print "{:.0f}% | {} of {} new tiles created | {} | {}".format((self.Global_totPoints / float(self.globalHeader['numptrecords'])) * 100,
                                                                                 len(os.listdir(self.tempDirectory)), len(self.grid),
                                                                                 os.path.split(tile)[1],
                                                                                 datetime.datetime.now())
                    self.Global_totPoints += 1

                    try:
                        # test x point first ...
                        fh.seek(h['offset'] + (i * h['pointreclen']))
                        x = struct.unpack('=' + 'L', fh.read(4))[0]
                        x = (x * h['xscale'] ) + h['xoffset']
                        if not self.xmin < x < self.xmax:
                            continue

                        # test y point next ...
                        fh.seek(h['offset'] + (i * h['pointreclen'] + 4))
                        y = struct.unpack('=' + 'L', fh.read(4))[0]
                        y = (y * h['yscale']) + h['yoffset']
                        if not self.ymin < y < self.ymax:
                            continue

                        # extract round plot
                        inGrid = False
                        idx = [(grid['x'] > h["xmin"]) | (grid['x'] <= h["xmax"]) |
                               (grid['y'] > h["ymin"]) | (grid['y'] <= h["ymax"])]
                        for row in grid[idx]:
                            if ((row[0] - (self.resolution / 2.)) < x < (row[0] + (self.resolution / 2.)) and
                                (row[1] - (self.resolution / 2.)) < y < (row[1] + (self.resolution / 2.))):
                                inGrid = True
                                break

                        if not inGrid: continue

                        fh.seek(h["offset"] + (i * h['pointreclen'])) # searches to beginning of point

                        point_dictionary = extract_return(self, fh)
                        X, Y = point_dictionary['x'], point_dictionary['y']
                        KEYS = self.keys[(X >= self.keys['x'] - (self.resolution / 2.)) & (X < self.keys['x'] + (self.resolution / 2.)) &
                                         (Y >= self.keys['y'] - (self.resolution / 2.)) & (Y < self.keys['y'] + (self.resolution / 2.))]

                        for key in KEYS:
                            self.xy_dictionary[tuple(key)] = self.writePoint(self.xy_dictionary[tuple(key)], point_dictionary, h)

                    except:
                        self.badPoints += 1

                # deletes .las tiles that were converted from .laz or copies
                if self.tempDirectory in tile: os.unlink(tile)

        if self.resolution >= 1: h['guid1'] = self.resolution

        return self

    def writePoint(self, tile, d, h):

        pointString = ''
        for i in self.ptFrmt:
            if i[0] == 'return_grp':
                byte = ((d['scan_edge'] & 1) << 7) | ((d['scan_dir'] & 1) << 6) | \
                       ((d['rtn_tot'] & 7) << 3) | (d['rtn_num'] & 7)
            elif i[0] == 'x':
                byte = (d['x'] - h['xoffset']) / h['xscale']
            elif i[0] == 'y':
                byte = (d['y'] - h['yoffset']) / h['yscale']
            elif i[0] == 'z':
                byte = (d['z'] - h['zoffset']) / h['zscale']
            else:
                byte = d[i[0]]
            pointString += struct.pack('=' + i[2], byte)

        with open(tile["tempFile"], "ab") as o:
            o.write(pointString)

        # updates header information
        if d['z'] > tile["zmax"]:
            tile["zmax"] = d['z']
        if d['z'] < tile["zmin"]:
            tile["zmin"] = d['z']
        if 0 < d["rtn_num"] < 6:
            tile["num_rtn"][d["rtn_num"]] += 1
        tile["i"] += 1

        return tile

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
            a = np.zeros(tile["i"], dtype = self.dt)
            with open(tile["tempFile"], "rb") as fh:
                for j, line in enumerate(range(tile["i"])):
                    fh.seek(j * self.globalHeader['pointreclen'])
                    d = extract_return(self, fh)
                    for field in d:
                        a[j][field] = d[field]
            arr.append(a)
            if self.keepTemp is False: os.unlink(tile["tempFile"])

        if self.keepTemp is False: shutil.rmtree(self.tempDirectory)

        if len(arr) == 1: # if list arr has only one array...
            arr = arr[0]  # ...return only the array

        return arr

    def exportLAS(self, out=False):

        """
        Writes tile(s) to .las format.
        """

        if len(self.xy_dictionary) > 10:
            nTileCounter = len(self.xy_dictionary) // 10
        else:   nTileCounter = len(self.xy_dictionary)

        self.xy_dictionary = {key:values for key, values in self.xy_dictionary.items() if values['i'] > 0}

        if len(self.xy_dictionary) > 0:
            for i, key in enumerate(self.xy_dictionary.keys()):

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
                if self.resolution != None:
                    self.h['guid1'] = self.resolution

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
                    if not os.path.isdir(self.savePath):
                        os.makedirs(self.savePath)
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
            raise Exception('ZNR is deprecated - use an older version of lasIO')
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

    ptFrmt = self.ptFrmt

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
#
#
# if __name__ == '__main__':
#
#     import shutil
#     path = '/Users/phil/ALS/WC/spl/tile_20/ForestLAS_tutorial/LAS/large_tile'
#     shutil.rmtree(os.path.join(path, 'WC45_SUB_20m_TILES'))
#     os.makedirs(os.path.join(path, 'WC45_SUB_20m_TILES'))
#     start = datetime.datetime.now()
#     # lasIO(os.path.join(path, 'WC45_SUB_5m_SUBSET'), verbose=True, number_of_processes=4).all().exportLAS()
#     lasIO(os.path.join(path, 'WC45_SUB.las'), verbose=True, number_of_processes=1).grid(os.path.join(path, 'coords_10.csv'), resolution=10).exportLAS()
#     # lasIO(os.path.join(path, 'WC45_SUB_5m_SUBSET'), verbose=True, number_of_processes=8).tiling(1000).exportLAS(os.path.join(path, 'WC45_SUB_20m_TILES'))
#     print datetime.datetime.now() - start