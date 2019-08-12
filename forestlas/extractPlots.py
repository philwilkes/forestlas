"""
Extract points for a given trap size around a given coordinate
"""

import os, sys, subprocess, csv, shutil, tempfile
import MySQLdb as mdb
import numpy as np
from lasIO import *
import multiprocessing
import traceback

class extract:

    def __init__(self, verbose=False, processes=4):

        self.verbose = verbose
        self.cwd = os.getcwd()
        self.maxProcess = int(processes)

    def query(self, points_raw, src="point", out=False, extent=24, 
              round=True, flightlines=False, search=".", tmpDir=True,
              namePrefix=''):
        
        """
        Is it a points query or an extent?
        """
        try:
            self.fl = flightlines
            self.round = round
            self.extent = extent
            self.search = "%" + search + "%"
            self.tmpDir = tmpDir
            self.namePrefix = namePrefix
            
            if not out:
                parentDir, outDir = os.path.split(points_raw)
                if len(parentDir) == 0: parentDir = os.getcwd()
                self.outDir = os.path.join(parentDir, os.path.splitext(outDir)[0]) 
                self.outLAS = False
            elif os.path.isdir(out):
                self.outDir = out
                self.outLAS = False
            elif out.endswith('.las'):
                self.outDir, self.outLAS = os.path.split(out)
            else:
                raise IOError('output destination not recognised')
        
            if os.path.isdir(self.outDir) is False:
                os.makedirs(self.outDir)
        
            if self.verbose == True:
                print "src: {}".format(src)
                print "extent: {}".format(self.extent)
                print "round: {}".format(round)
                print "flightlines selected: {}".format(self.fl)
                
            if src == "point":
                self.read_points(points_raw)
            else:
                plots = read_extent(points_raw)
        
        except Exception as err:
            print traceback.format_exc()
        finally:
            os.chdir(self.cwd)
 
    def read_points(self, points_raw):

        self.tiles_to_process = []
    
        with open(points_raw) as pr:
            
            for i in pr.read().split('\n'):
                if i.find('#') != -1:
                    continue
                elif i != "":
                    i = str.split(i, ',')
                    point = str(i[0])
                    x = float(i[1])
                    y = float(i[2])
                else:    continue
               
                #names point if not 
                if point == None:
                    point = str(x) + '_' + str(y)
                
                #print "processing point: {} x: {} y: {}".format(point, x, y)
                
                # calculates search self.extent
                xmin = float(x)-(self.extent/2.)
                ymin = float(y)-(self.extent/2.)
                xmax = float(x)+(self.extent/2.)
                ymax = float(y)+(self.extent/2.)

                if xmin > ymin:
                    xmin, ymin, xmax, ymax = ymin, xmin, ymax, xmax
                
                self.tiles_to_process.append({"xmin":xmin, "ymin":ymin, 
                                              "xmax":xmax, "ymax":ymax, 
                                              "point":point})

        if self.maxProcess == 0:
            self.select_tiles(xmin, ymin, xmax, ymax, point)
        else: self.mp()

    def mp(self):

        listI = 0

        for i in range((len(self.tiles_to_process) / self.maxProcess) + 1):
        
            jobs = []
    
            try:
                if (self.maxProcess * listI) + self.maxProcess < len(self.tiles_to_process):
                    processingList = self.tiles_to_process[self.maxProcess * listI: (self.maxProcess * listI) + self.maxProcess]
                else:   processingList = self.tiles_to_process[self.maxProcess * listI:]

                
            
                for j, t in enumerate(processingList): # limits number of images run at once
                    p = multiprocessing.Process(target=self.select_tiles, args=(t["xmin"], 
                                                                                t["ymin"], 
                                                                                t["xmax"], 
                                                                                t["ymax"], 
                                                                                t["point"], ))
                    jobs.append(p)
                    p.start()
    
                for proc in jobs:
                    proc.join()
                
            except:
                print self.tiles_to_process[0:1]
                print self.maxProcess, type(self.maxProcess)
                print listI, type(listI)
                raise NameError
    
            listI += 1
    
    def select_tiles(self, xmin, ymin, xmax, ymax, point):

        ## connects to db
        try:
            con = mdb.connect('127.0.0.1', 'seo', 'lidar', 'lidar') # server, user, pw, db
            cur = con.cursor()
        except Exception, err:
            raise Exception(err)

        ## selects either raw flightlines or processed .las files
        if self.fl == True:
            fl = '"%classed"'
        else:
            fl = '"%_height"'

        # queries db for available als tiles
        query = ('select concat(path, "/", tilename, ".",  format) from tiles where \
            (xmin > %(xmin)s and xmin < %(xmax)s and ymin > %(ymin)s and ymin < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax > %(xmin)s and xmax < %(xmax)s and ymax > %(ymin)s and ymax < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmin > %(xmin)s and xmin < %(xmax)s and ymax > %(ymin)s and ymax < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax > %(xmin)s and xmax < %(xmax)s and ymin > %(ymin)s and ymin < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax > %(xmax)s and xmin < %(xmin)s and ymax > %(ymin)s and ymax < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax > %(xmax)s and xmin < %(xmin)s and ymin > %(ymin)s and ymin < %(ymax)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax < %(xmax)s and xmax > %(xmin)s and ymax > %(ymax)s and ymin < %(ymin)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmin < %(xmax)s and xmin > %(xmin)s and ymax > %(ymax)s and ymin < %(ymin)s and (path like "%(name)s" or tilename like "%(name)s")) or \
            (xmax > %(xmax)s and xmin < %(xmin)s and ymax > %(ymax)s and ymin < %(ymin)s and (path like "%(name)s" or tilename like "%(name)s"))\
            ' % {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax, 'name':self.search})

#         if self.verbose is True:
#             print query

        # creates list of tiles from query
        cur.execute(query)
        ans = cur.fetchall()
        self.tile_list = [a[0] for a in ans]
        cur.close() # closes cursor
        con.close() # closes db connection

        if len(self.tile_list) == 0 and self.verbose == True:
            print query

        ## self.extent of tile
        print 'processing %s tile(s) for point %s' % (len(self.tile_list), point)
        if self.verbose:
            for tile in self.tile_list:
                print tile

        if len(self.tile_list) > 0:
            if not self.outLAS:
                self.outLAS = os.path.join(self.outDir, 
                                           '{}_'.format(self.namePrefix) + point + '.las')
            self.clipLAS(xmin, ymin, xmax, ymax)
        else:
            print "outside of bounds: {}".format(point)
            with open(os.path.join(self.outDir, "log.csv"), "a") as log:
                log.write(",".join([point, "Out of bounds", str(np.mean([xmax, xmin])), 
                                    str(np.mean([ymax, ymin])), "\n"]))
            
    def clipLAS(self, xmin, ymin, xmax, ymax):
    
        x = np.mean([xmax, xmin])
        y = np.mean([ymax, ymin])
        
        if not self.tmpDir or self.maxProcess == 0:
            tempDirectoryName = "lidar.processing." + str(np.random.randint(0, 9999999)) + ".tmp"
        else:    
            pid = multiprocessing.current_process()._identity[0]
            tempDirectoryName = "lidar.processing." + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + ".tmp"
        self.tmpDir = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        
        lasIO(self.tile_list, 
              out=self.outLAS, 
              verbose=self.verbose, 
              znr=False, 
              search=self.search,
              tmpDir=self.tmpDir,
              keepTemp=False).plot(x, y,
                                   extent=self.extent, 
                                   round=self.round).exportLAS()

        if self.verbose:
            print ".las exported to: {}".format(self.outLAS)
            
        if not os.path.isfile(self.outLAS):
            with open(os.path.join(self.outDir, "log.csv"), "a") as log:
                log.write(",".join([point, ",".join(self.tile_list), str(x), str(y), "\n"]))
            
if __name__=='__main__':
    extract(verbose=False).query(sys.argv[1])
