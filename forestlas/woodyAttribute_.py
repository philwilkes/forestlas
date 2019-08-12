import os
import glob
import multiprocessing
import datetime
import tempfile
import numpy as np

from forestlas.lasIO import *

class canopyStructure:
    
    """
    Input is a directory containing tiled .las files and generates
    a raster of desired metrics, the defult is canopy height

    TO DO:
    1: Get projection information from .las
    """
    
    def __init__(self, las_files, alpha=.3, N=10, resolution=None,
                 layers=False, canopyDensity=False, underDensity=False,
                 height=False, pgap=False, baseHeight=False,
                 canopyPgap=False, all=False, point_only=False,
                 verbose=False, number_of_processes=4, threshold=2, points_per_voxel=2):

        self.verbose = verbose
        self.number_of_processes = number_of_processes

        if not layers \
       and not canopyDensity \
       and not underDensity \
       and not height \
       and not pgap \
       and not baseHeight \
       and not canopyPgap \
       and not all \
       and not point_only:
            raise Exception('no method chosen')

        if all:
            layers = True
            canopyDensity = True
            underDensity = True
            height = True
            pgap = True
            baseHeight = True
            canopyPgap = True

        self.alpha = alpha
        self.N = N
        self.threshold = threshold
        self.points_per_voxel = points_per_voxel

        self.metrics = {'layers':False,
                        'canopyDensity':False,
                        'underDensity':False,
                        'height': False,
                        'baseHeight': False,
                        'pgap': False,
                        'canopyPgap':False}

        if layers:
            self.metrics['layers'] = True
        if canopyDensity:
            self.metrics['canopyDensity'] = True
        if underDensity:
            self.metrics['underDensity'] = True
        if height:
            self.metrics['height'] = True
        if pgap:
            self.metrics['pgap'] = True
        if baseHeight:
            self.metrics['baseHeight'] = True
        if canopyPgap:
            self.metrics['canopyPgap'] = True
        if point_only:
            self.metrics['point_only'] = True

        if isinstance(las_files, list):
            self.dir = os.path.split(las_files[0])[0]
            if not os.path.isdir(self.dir):
                if os.path.isfile(os.path.join(os.getcwd(), las_files[0])):
                    self.dir = os.getcwd()
                    self.LASlist = las_files
                else:
                    raise Exception('las files not in cwd and no directory supplied (suggest using glob)')
            else:
                self.LASlist = las_files
        elif os.path.isdir(las_files):
            self.dir = las_files
            self.LASlist = sorted(glob.glob(os.path.join(self.dir, '*.las')))
        else:
            raise Exception('input needs to be a directory or list of .las files')

        self.resolution = resolution
        self.counter = len(self.LASlist) // 100

        if len(self.LASlist) == 0: raise Exception('No .las files in {}'.format(self.dir))
        if self.verbose: print 'number of tiles to process:', len(self.LASlist)

        self.createGrid()
        self.mp(self.LASlist)
        self.populateGrid()

    def base(chp):
        return chp.threshold if len(chp.crownBase) == 0 else chp.crownBase[-1]

    def calculateStructure(self, las, x, y, lasF):

        from forestlas.canopyComplexity import CanopyComplexity

        ### Generate plot profile ###

        try:
            if las['z'].max() < self.threshold or len(las[las['z'] >= self.threshold]) < self.points_per_voxel:
                for metric in self.metrics.keys():
                    if metric != 'height': self.metrics[metric][x, y] = np.nan

            else:
                model = CanopyComplexity().fromLAS(las).CHP(method='model')
                results = {metric:np.zeros(self.N) for metric in self.metrics.keys()}
                for i in range(self.N):
                    sim = model.simulateCloud()
                    if sim['z'].max() < self.threshold or len(sim[sim['z'] >= self.threshold]) < self.points_per_voxel:
                        results['layers'][i] = 0
                        results['canopyDensity'][i] = 0
                        results['underDensity'][i] = 1
                        results['baseHeight'][i] = np.nan
                        results['pgap'] = 1
                        results['canopyPgap'] = np.nan
                    else:
                        chp = CanopyComplexity().fromLAS(sim, top_threshold=100).CHP(method='sample', alpha=self.alpha)
                        for metric in results.keys():
                            if metric == 'layers':
                                results[metric][i] = chp.layerCount
                            if metric == 'canopyDensity':
                                results[metric][i] = chp.fd[chp.zx >= self.base(chp)].sum() / float(chp.fd.sum())
                            if metric == 'underDensity':
                                results[metric][i] = chp.fd[chp.zx < self.base(chp)].sum() / float(chp.fd.sum())
                            if metric == 'baseHeight':
                                results[metric][i] = np.nan if len(chp.crownBase) == 0 else chp.crownBase[-1]
                            if metric == 'pgap':
                                results[metric][i] = chp.pgap.min()
                            if metric == 'canopyPgap':
                                results[metric][i] = chp.pgap[chp.zx == (chp.threshold if len(chp.crownBase) == 0
                                                                                       else chp.crownBase[-1])][0]

                for metric in self.metrics.keys():
                    if metric != 'height': self.metrics[metric][x, y] = results[metric].mean()

        except Exception as err:
            print '!!!', lasF, err, '!!!'

        return self

    def mp(self, listItems):

        ''' Carries out the multiprocessing grunt work '''

        manager = multiprocessing.Manager()
        for metric, do in self.metrics.items():
            if self.metrics[metric]:
                self.metrics[metric] = manager.dict()
            else:
                del self.metrics[metric]

        # self.plot_dictionary = manager.dict()
        self.global_x = manager.list()
        self.global_y = manager.list()

        listI = 0

        for i in range((len(listItems) / self.number_of_processes) + 1):

            if self.counter > 0 and i%self.counter == 0 and self.verbose:
                print '{:.2f}% | processing job {} of {} | {}'.format((float(i) / ((len(listItems) / self.number_of_processes) + 1)) * 100, \
                                                                      i, (len(listItems) / self.number_of_processes) + 1, \
                                                                      datetime.datetime.now())
            jobs = []

            if (self.number_of_processes * listI) + self.number_of_processes < len(listItems):
                processingList = listItems[self.number_of_processes * listI: (self.number_of_processes * listI) + self.number_of_processes]
            else:   processingList = listItems[self.number_of_processes * listI:]

            for j, las in enumerate(processingList): # limits number of lass run at once

                p = multiprocessing.Process(target=self.readLAS, args=(las, ))

                if p:
                    jobs.append(p)
                    p.start()

            for proc in jobs:
                proc.join()

            listI += 1

    def readLAS(self, lasF):

        # read las file and send to different processes
        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = 'lidar.processing.' + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + '.tmp'
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)

        las = lasIO(lasF, tmpDir=tempDirectory)
        X, Y = self.tileCentre(las)

        if 'point_only' in self.metrics.keys():
            self.metrics['point_only'][X, Y] = 1

        else:
            las = las.all().asArray()

            if 'height' in self.metrics.keys():
                self.metrics['height'][X, Y] = np.percentile(las['z'], 95)

            if len(set(self.metrics.keys()).intersection(['layers', 'canopyDensity',
                                                          'underDensity', 'pgap',
                                                          'baseHeight'])) > 0:
                self.calculateStructure(las, X, Y, lasF)

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

    def populateGrid(self):

        self.metricGrids = {key:None for key in self.metrics.keys()}

        for metric in self.metrics.keys():

            # array, self.X, self.Y = np.meshgrid(np.nan, np.unique(self.global_x), np.unique(self.global_y))
            array, self.X, self.Y = np.meshgrid(np.nan,
                                                np.arange(min(self.global_x), max(self.global_x) + 1, self.resolution),
                                                np.arange(min(self.global_y), max(self.global_y) + 1, self.resolution))

            for key in self.metrics[metric].keys():
                idx = [(self.X == key[0]) & (self.Y == key[1])]
                array[idx] = np.float(self.metrics[metric][key])

            self.metricGrids[metric] = np.rot90(array.reshape(np.shape(array)[0], np.shape(array)[2]))

        return self

    def tileCentre(self, las):


        x = np.floor(las.x_centre)
        self.global_x.append(x)

        y = np.round(las.y_centre)
        self.global_y.append(y)

        return x, y

    def exportTiff(self, saveTo=False):

        from osgeo import osr, gdal

        for metric in self.metrics.keys():
            if saveTo:
                if os.path.isfile(saveTo):
                    ans = raw_input('Image already exists, overwrite? (Y|N): ').lower()
                    if ans == 'y':
                        os.unlink(saveTo)
                    else:
                        raise NameError('Change save filepath')
                    savePath = saveTo
                elif os.path.isdir(saveTo):
                    savePath = os.path.join(saveTo, '{}_{:.2f}.tif'.format(metric, self.resolution))
                elif saveTo.endswith('.tif'):
                    dir = os.path.split(saveTo)[0]
                    if os.path.isdir(dir):
                        savePath = saveTo
                    else:
                        raise Exception('{} is not a directory'.format(dir))
            else:
                savePath = os.path.join(self.dir, '{}_{:.2f}.tif'.format(metric, self.resolution))

            driver =  gdal.GetDriverByName('GTiff')
            xDim = len(np.unique(self.X))
            yDim = len(np.unique(self.Y))
            dataset = driver.Create(savePath, xDim, yDim, 1, gdal.GDT_Float32)

            # set projection
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(28355)
            dataset.SetProjection(srs.ExportToWkt())

            dataset.SetGeoTransform([self.X.min() - float(self.resolution) / 2.,
                                     self.resolution,
                                     0,
                                     self.Y.max() + float(self.resolution) / 2.,
                                     0,
                                     0 - self.resolution])

            # set no data value
            dataset.GetRasterBand(1).SetNoDataValue(-1)

            # write raster and close
            dataset.GetRasterBand(1).WriteArray(self.metricGrids[metric])
            del dataset # closes .tif

            print '{} tif saved at {}'.format(metric, savePath)

    def asArray(self):
        return self.metricGrids

if __name__ == '__main__':

    import glob
    # path = '/Users/phil/ALS/WC/spl/tile_20/ForestLAS_tutorial/LAS/large_tile'
    path = '/Users/phil/ALS/WC/spl/tile_20/ForestLAS_tutorial/LAS/large_tile/WC45_SUB_20m_TILES'
    L = glob.glob(os.path.join(path, '*.las'))
    # L = ['/Users/phil/Google_Drive/RiverHealth/DATA/AREA_A/10/309231.0_6004106.0.PLOT.las']
    # canopyStructure(L, point_only=True, resolution=10, verbose=True)
    start = datetime.datetime.now()
    canopyStructure(L, height=True, verbose=True, number_of_processes=4).exportTiff()
    print datetime.datetime.now() - start
    # G = canopyStructure(glob.glob('/Users/phil/Google_Drive/RiverHealth/DATA/AREA_C/100/*.las'), height=True, verbose=True, number_of_processes=4, resolution=100)#