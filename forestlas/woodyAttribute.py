import sys, os, glob, multiprocessing, datetime, shutil, tempfile
import numpy as np
import matplotlib.pyplot as plt
from osgeo import osr, gdal

from forestlas.lasIO import *
from forestlas.canopyComplexity import *

class woodyAttribute:
    
    """
    Input is a directory containing tiled .las files and generates
    a raster of desired metrics, the defult is canopy height

    TO DO:
    1: Get projection information from .las
    """
    
    def __init__(self, verbose=False, maxProcess=4):
        self.verbose = verbose
        self.maxProcess = maxProcess
        
    def canopyHeight(self, directory, resolution=None, height=95):
    
        self.metric = 'height' 
        self.height = height
        
        self.processSetup(directory, resolution)
        return self
        
    def canopyComplexity(self, directory, alpha=.3, N=20, resolution=None):
              
        self.metric = 'complexity'
        self.alpha = alpha
        self.N = N

        self.processSetup(directory, resolution)
        return self
        
    def Pgap(self, directory, N=20, resolution=None, ):
              
        self.metric = 'Pgap'   
        self.N = N

        self.processSetup(directory, resolution)
        return self
        
    def fhd(self, directory, resolution=None, ):
              
        self.metric = 'fhd'   

        self.processSetup(directory, resolution)
        return self
        
    def Cv(self, directory, resolution=None, ):
              
        self.metric = 'Cv'   

        self.processSetup(directory, resolution)
        return self
        
    def fracCover(self, directory, resolution=None, cellSize=1, threshold=.3):
    
        self.metric = 'FC'
        self.cellSize = cellSize
        self.threshold = threshold 

        self.processSetup(directory, resolution)
        return self

    def processSetup(self, directory, resolution):            

        self.dir = directory
        self.resolution = resolution

        self.LASlist = sorted(glob.glob(os.path.join(self.dir, '*.znr')))
        if len(self.LASlist) == 0:
            self.LASlist = sorted(glob.glob(os.path.join(self.dir, '*.las')))
        self.counter = len(self.LASlist) // 100

        if len(self.LASlist) == 0: raise Exception('No .las files in {}'.format(self.dir))
        
        self.createGrid()
        self.mp(self.LASlist)
        self.v, self.X, self.Y = self.populateGrid(self.plot_dictionary)
        return self
    
    def calculateComplexity(self, las):
        
        ### Generate plot profile ###
        
        # calculate plot centre
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        model = CanopyComplexity(mp=tempDirectory).fromLAS(las).CHP(method='model')
        
        if model.z.max() < 2:
            self.plot_dictionary[x, y] = 0
        
        else:
            results = np.zeros(self.N)

            for i in range(self.N):
                z, zw = model.simulateCloud
                chp = CanopyComplexity().fromSample(z, zw).CHP(alpha=self.alpha)
                results[i] = chp.layerCount
             
            if len(results) > 0:
                self.plot_dictionary[x, y] = results.mean()
            else:    self.plot_dictionary[x, y] = 0
            
    def calculatePgap(self, las):
        
        ### Generate plot profile ###
        
        # calculate plot centre
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        
        pgap = CanopyComplexity(mp=tempDirectory).fromLAS(las, threshold = 1).Pgap()

        self.plot_dictionary[x, y] = pgap.min()
        
    def calculateFC(self, las):
        LAS = lasIO.lasIO(las)
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        tiles = LAS.tiling(self.cellSize).asDic()
        o = np.zeros(len(tiles))
        for i, tile in enumerate(tiles.keys()):
            T = tiles[tile]
            if len(T) > 0 and T['z'].max() > 1.7:
                pgap = CanopyComplexity().fromSample(T['z'], T['rtn_tot'], threshold=0).Pgap()
                o[i] = pgap.pgap[np.where(pgap.z > 1.7)][0]
            else:
                o[i] = 0
        self.plot_dictionary[x, y] = len(o[o > self.threshold]) / float(len(o))
    
    def calculateHeight(self, las):
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        las = lasIO.lasIO(las, tmpDir=tempDirectory, keepTemp=False).all().asArray()
        self.plot_dictionary[x, y] = np.percentile(las['z'], self.height)
        
    def calculateCv(self, las):
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        las = lasIO(las, tmpDir=tempDirectory, keepTemp=False).all().asArray()
        self.plot_dictionary[x, y] = las['z'].std() / las['z'].mean()
        
    def calculateFHD(self, las):
        x, y = self.tileCentre(parseHeader(las)) # calculate plot centre
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        chp = CanopyComplexity(mp=tempDirectory, keepTemp=False).fromLAS(las)
        fhd = []
        for h in np.linspace(2, chp.z.max() // 2 * 2, chp.z.max() // 2):
            pi = sum(chp.zw[(chp.z > h) & (chp.z < h + 2)]) / chp.zw.sum()
            if pi == 0:
                continue
            fhdh = pi * np.log(pi)
            if np.isnan(fhdh):    fhdh = 0
            fhd.append(fhdh)
        self.plot_dictionary[x, y] = -sum(fhd)
            
    def mp(self, listItems):
        
        ''' Carries out the multiprocessing grunt work '''
        
        manager = multiprocessing.Manager()
        self.plot_dictionary = manager.dict()
        self.global_x = manager.list()
        self.global_y = manager.list()
        
        listI = 0            
        
        for i in range((len(listItems) / self.maxProcess) + 1):
            
            if self.verbose:
                print '{:.2f}% | processing job {} of {} | {}'.format((float(i) / ((len(listItems) / self.maxProcess) + 1)) * 100, \
                                                                      i, (len(listItems) / self.maxProcess) + 1, \
                                                                      datetime.datetime.now())
            jobs = []
            
            if (self.maxProcess * listI) + self.maxProcess < len(listItems):
                processingList = listItems[self.maxProcess * listI: (self.maxProcess * listI) + self.maxProcess]
            else:   processingList = listItems[self.maxProcess * listI:]
            
            for j, las in enumerate(processingList): # limits number of lass run at once
                
                p = False
                processMap = {'height':self.calculateHeight, 'complexity':self.calculateComplexity,
                              'fhd':self.calculateFHD, 'Cv':self.calculateCv, 'Pgap':self.calculatePgap,
                              'FC':self.calculateFC}

                p = multiprocessing.Process(target=processMap[self.metric], args=(las, ))
                
                if p:
                    jobs.append(p)
                    p.start()
            
            for proc in jobs:
                proc.join()
            
            listI += 1
            
    def createGrid(self):
                    
        # generate parmaeters from .las header
        header = parseHeader(self.LASlist[0]) # .las header to dictionary
    
        # search for predefined resolution if not there calculate 
        # tile resolution from data
        if self.resolution != None:
            pass
        elif 'guid1' in header.keys() and header['guid1'] > 0:
            self.resolution = header['guid1']
        elif 'resolution' in header.keys():
            self.resolution = header['resolution']
        else:
            self.resolution = header['xmax'] - header['xmin'] # a wild guess!

        if self.verbose == True:    print 'grid resolution: {}'.format(self.resolution)

        # grabself.vlr info
        # will use this at a later date to grab projection info
        self.vlr = getVLR(header['headersize'], self.LASlist[0])

        return self

    def populateGrid(self, plot_dictionary):
        
        self.v, self.X, self.Y = np.meshgrid(0., np.unique(self.global_x), np.unique(self.global_y))
       
        for key in plot_dictionary.keys():
            idx = [(self.X == key[0]) & (self.Y == key[1])]
            self.v[idx] = np.float(plot_dictionary[key])
        
        self.v = np.rot90(self.v.reshape(np.shape(self.v)[0], np.shape(self.v)[2]))

        print self.X
        return self.v, self.X, self.Y
    
    def tileCentre(self, header):

        x_range = header['xmax'] - header['xmin']
        x_centre = (header['xmax'] - (x_range / 2.))
        x_min = x_centre - (self.resolution / 2.)
        x = x_min // self.resolution
        self.global_x.append(x)
        
        y_range = header['ymax'] - header['ymin']
        y_centre = (header['ymax'] - (y_range / 2.))
        y_max = y_centre + (self.resolution / 2.) 
        y = y_max // self.resolution
        self.global_y.append(y)
        
        return x, y

    def exportTiff(self, saveTo=False):
    
        print 'writing to tiff'

        if saveTo:
            if os.path.isfile(saveTo):
                ans = raw_input('Image already exists, overwrite? (Y|N): ').lower()
                if ans == 'y':
                    shutil.rmtree(path)
                else:
                    raise NameError('Change save filepath')
            elif os.path.isdir(saveTo):
                saveTo = os.path.join(saveTo, '{}_{:.2f}.tif'.format(self.metric, self.resolution))
        else:
            saveTo = os.path.join(self.dir, '{}_{:.2f}.tif'.format(self.metric, self.resolution))
             
        driver =  gdal.GetDriverByName('GTiff')
        xDim = int((self.X.max() - self.X.min()) + 1)
        yDim = int((self.Y.max() - self.Y.min()) + 1)
        dataset = driver.Create(saveTo, xDim, yDim, 1, gdal.GDT_Float32)

        # set projection
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(28355)
        dataset.SetProjection(srs.ExportToWkt())
        
        # set transformation
        dataset.SetGeoTransform([self.X.min() * self.resolution, self.resolution, 0, self.Y.max() * self.resolution, 0, 0 - self.resolution])
        
        # write raster and close
        dataset.GetRasterBand(1).WriteArray(self.v.astype(float))
        dataset = None # closes .tif

        print '{} tif saved at {}'.format(self.metric, saveTo)

    def asArray(self):
        return self.v

if __name__ == '__main__':

    G = woodyAttribute(verbose=True, maxProcess=4).canopyHeight('/Users/phil/ALS/WC/spl/tile_20/WC1_5m_TILES').exportTiff(saveTo='/Users/phil/ALS/WC/spl/tile_20/WC1_5m_TILES/height_5.00.old.tif')