# Wilkes et al. in review. Using discrete-return ALS for the assessment 
# of vertical canopy structure across diverse forest types
# Methods in Ecology and Evolution

from scipy.stats import *
import numpy as np
from scipy.interpolate import UnivariateSpline
import glob
import os
import tempfile
import multiprocessing
import numpy.lib.recfunctions as rfn

from forestlas.lasIO import *

np.seterr(invalid="ignore")

class CanopyComplexity:
    
    def __init__(self, verbose=False, mp=False):
        
        self.verbose = verbose
        
        if mp is not False:
            self.tempDir = mp
        else:
            self.tempDir = False
            
    def fromLAS(self, las, threshold=2.0, top_threshold=99, z_scale=False):
        
        self.threshold = threshold

        if isinstance(las, str):
            lasRAW = lasIO(las, tmpDir=self.tempDir, keepTemp=False).all().asArray() # reads las file
        
            if lasRAW['rtn_tot'].max() == 0:
                lasRAW['rtn_tot'][:] = 1
            
            self.las = lasRAW[lasRAW['rtn_tot'] > 0] # inds lines where rtn_tot == 0 and removes them
            self.z_scale = parseHeader(las)['zscale']

        elif isinstance(las, np.ndarray):

            self.las = las
            if not z_scale:
                self.z_scale = .01

        else:
            raise Exception('input needs to path to LAS file or a LAS file array')

        self.z = self.las['z']
        self.zw = self.las['rtn_tot']
        self.lenLAS = float(len(self.las))

        rtn_weight = np.around(1. / self.las['rtn_tot'], decimals=2)
        self.total = np.sum(rtn_weight) # sum of weighted returns
        
        # removes ground and outliers
        idx = [(self.z > threshold) & (self.z < np.percentile(self.z, top_threshold))]
        self.z = self.z[idx] # ...from z
        self.zw = rtn_weight[idx] # ... and from weighted heights
        self._create_bins()

        return self

    def _create_bins(self):

        # create bin array
        #### required because of numpy float issue ###
        factor = 1 / self.z_scale
        z_min, z_max = self.threshold * factor, int((self.z.max() + (self.z_scale*1000)) * factor)
        ##############################################
        self.bins = z_max-(z_min-1) # number of bins
        self.zxOrig = np.linspace(z_min, z_max, self.bins) / factor # "x-axis"

        return self

    def Pgap(self, frequency=False):
    
        if len(self.z) == 0:
            self.pgap = np.array([1])
        else:
            #calculatePgap
            if not frequency:
                self.pgap = np.zeros(self.bins) # array for populating with pgap
                for (i, height_below_toc) in enumerate(self.zxOrig):
                    idx = [self.z >= height_below_toc] # index all rows >= height z
                    weight_of_returns_above_z = sum(self.zw[idx]) # calculates sum of weighted returns above height z
                    self.pgap[i] = 1. - (weight_of_returns_above_z / self.total) # populates pgap with proportion of weight
            else:
                #calculateFreqPgap
                self.pgap = np.zeros(self.bins) # array for populating with pgap
                for (i, height_below_toc) in enumerate(sorted(self.zxOrig)):
                    idx = [self.z >= height_below_toc] # index all rows >= height z
                    num_of_returns_above_z = len(self.z[idx]) # calculates len of z array above height z
                    self.pgap[i] = 1. - (num_of_returns_above_z / self.lenLAS) # populates pgap with proportion of weight
        
        return self
    
    def CHP(self, method="sample", alpha=.3, frequency=False, normalise=False, noise=0.05):

        """
        method is "model" for a log transformed and normalised CHP and "sample" for anything else!
        frequncy uses retrun frequency instead of weight
        """

        if len(self.z) == 0:
            self.pgap = np.array([1])
            self.fd = np.array([])
            self.sd = np.array([])
            self.spline = None
            self.ps = np.array([])
            self.zx = np.array([1])
            self.layerLocation = np.array([])
            self.layerCount = 0
            self.crownBase = np.array([])
        
        else:
            if method == 'sample':
                normalise = False
                self.alpha = alpha
                log_transform = False
            elif method == "model":
                normalise = True
                self.alpha = 0
                log_transform = True
            elif method == 'log_transform':
                normalise = False
                self.alpha = alpha
                log_transform = True
            else:
                raise Exception('method not recognised')

            self.Pgap()
            
            #smooth_pgap
            self.spline = UnivariateSpline(self.zxOrig, self.pgap, s=self.alpha)
            self.ps = self.spline(self.zxOrig)
        
            # clips ps and zx vectors to maximum height
            self.ps = self.ps[np.where(self.zxOrig < self.z.max() * 1)]
            self.pgap = self.pgap[np.where(self.zxOrig < self.z.max() * 1)]
            self.zx = self.zxOrig[np.where(self.zxOrig < self.z.max() * 1)]
            
            # log transformation ... or not
            if log_transform:
                self.cc = -np.log(self.ps)
            else:
                self.cc = 1 - self.ps

            # first_derivative
            #  self.fd1 = np.hstack([np.diff(self.cc[::-1]) / self.z_scale, 0])[::-1]
            self.fd = -(np.gradient(self.cc) / self.z_scale)

            # removes any lingering negative values
            if self.fd.min() < 0:
                self.fd = self.fd + (-self.fd.min())

            # normalise so sum of vector == 1, required for probability
            if normalise:
                self.fd = self.fd * ((1-self.pgap.min())/np.sum(self.fd))
        
            #second_derivative
            self.sd =  np.gradient(self.fd) / self.z_scale 
        
            #number_of_modes
            signs = np.diff(self.sd/abs(self.sd)) # finds zero crossings
            idx = [signs == -2] 
            potentialLayerLocation = self.zx[idx] # and their height
            layerAmplitude = self.fd[idx] # and the signal amplitude
            maxAmplitude = self.fd.max() # and the maximum amplitude for the CHP
            self.layerLocation = [layer for i, layer in enumerate(potentialLayerLocation) if layerAmplitude[i] > maxAmplitude * noise] # and filters noise
            self.layerCount = len(self.layerLocation)

            idx = [signs == 2] 
            self.crownBase = self.zx[idx] # and their height
        
        return self

    def simulateCloud(self):

        """
        This simulates a height vector and then assigns heights a weight
        """

        self.zAll = np.hstack([0, self.zx]) # add ground to height bins
        
        # height weighting vector
        self.heightWeight = np.hstack([self.pgap.min(), self.fd]) # weights for each height bin and add ground weight
        #self.heightWeight = np.ma.masked_array(self.heightWeight, np.isnan(self.heightWeight))
        
        # return weighting vector
        self.returnWeight = {} # dictionary to store self.returnWeight
        for h in np.unique(np.floor(self.zAll)): # rounds xs down to 1 metre bins
            # selects returns with heights in bin h
            self.returnWeight[h] = {} # creates dictionary within self.returnWeight to store count of NoR values
            idx = [(self.z > h) & (self.z <= h + 1)]
            NoR = self.zw[idx] # total number of returns for returns in range h to h+1
            for rtn in np.unique(NoR):
                self.returnWeight[h][rtn] = len(NoR[NoR == rtn]) # counts number of returns by rtn_tot
            sumNoR = np.sum(self.returnWeight[h].values()) # counts number of returns in bin
            for rtn in self.returnWeight[h].keys(): # return values in height bin
                # self.returnWeight[h][rtn] = self.returnWeight[h][rtn] / np.float(sumNoR) # calculates weight
                self.returnWeight[h][rtn] /= np.float(sumNoR)

        # Simulated height
        self.simHeightAll = np.random.choice(self.zAll, max(self.lenLAS, 100), p=self.heightWeight)
        #self.simHeightAll = np.random.choice(self.zAll, 100, p=self.heightWeight)

        # Simulated weights
        self.simRtnWeight = np.zeros(len(self.simHeightAll))  # array to store weights
        for i, z in enumerate(self.simHeightAll):
            hgt_bin = np.floor(z)
            for offset in [0, 1, -1, 2, -2, 3, -3]: # if weighting bin is empty selects a neigbour
                try:
                    rtn_num = self.returnWeight[hgt_bin + offset].keys() # return numbers
                    rtn_weights = self.returnWeight[hgt_bin + offset].values() # probability of return number
                    """ return number randomly chosen (weighted by probability of return) """
                    self.simRtnWeight[i] = np.random.choice(rtn_num, 1, p=rtn_weights) # assignment of NoR value (divided by 1)
                    break
                except:
                    self.simRtnWeight[i] = 1
        
        return rfn.merge_arrays([np.array(self.simHeightAll, dtype=[('z', np.float)]),
                                 np.array(self.simRtnWeight, dtype=[('rtn_tot', np.float)])])
        
class bootstrapComplexity:
    
    import multiprocessing
    
    def __init__(self, las, verbose=False, processes=1, N=100):

        self.N = N
        self.verbose = verbose
        
        if type(las) is str:
            if os.path.isdir(las):
                self.l = glob.glob(os.path.join(las, "*.znr"))
                if len(self.l) == 0: self.l = glob.glob(os.path.join(las, "*.las"))
            elif os.path.isfile(las): self.l = [las]
        elif type(las) is list:
            self.l = las
        else:
            raise IOError("No .las or .znr files in {}".format(las))
            
        self.mp(processes)
        
    def chp(self, las):
    
        if self.verbose: print 'processing:', las

        pid = multiprocessing.current_process()._identity[0]
        tempDirectoryName = "lidar.processing." + str(np.random.mtrand.RandomState(pid).randint(0, 9999999)) + ".tmp"
        tempDirectory = os.path.join(tempfile.gettempdir(), tempDirectoryName)
        
        self.bsCHP = np.zeros(self.N)
        
        self.chp = CanopyComplexity(mp=tempDirectory).fromLAS(las).canopyHeightProfile("model")

        for i in range(self.N):
    
            z, zw = self.chp.simulateCloud

            if z.max() < 2 or len(z) < 2 or 0 in zw:
                self.bsCHP[i] = 0
            else:
                sample = CanopyComplexity().fromSample(z, zw).canopyHeightProfile()
                self.bsCHP[i] = sample.layerCount
        
        if type(las) is str: plotName = os.path.split(os.path.splitext(las)[0])[1]
        else: plotName = 'array'
        self.chp_dictionary[plotName] = self.bsCHP
        return self
  
    def mp(self, maxProcess):

        listI = 0
        manager = multiprocessing.Manager()
        self.chp_dictionary = manager.dict()

        for i in range((len(self.l) / maxProcess) + 1):
        
            jobs = []
    
            if (maxProcess * listI) + maxProcess < len(self.l):
                processingList = self.l[maxProcess * listI: (maxProcess * listI) + maxProcess]
            else:   processingList = self.l[maxProcess * listI:]
    
            for j, las in enumerate(processingList): # limits number of images run at once
                p = multiprocessing.Process(target=self.chp, args=(las, ))
                jobs.append(p)
                p.start()
    
            for proc in jobs:
                proc.join()
    
            listI += 1
        
        self.bsCHPmutiple = dict(self.chp_dictionary)

if __name__ == '__main__':

    path = '/Users/phil/ALS/WC/spl/tile_20/WC1_5m_TILES/383475.0_5828910.0.las'
    las = lasIO(path).all().asArray()
    # plot = 'PE2744N2556'
    # las_path = os.path.join(path, plot + '.las')
    # las = CanopyComplexity().fromSample(las['z'], las['rtn_tot']).CHP('model')
    las = CanopyComplexity().fromLAS(las).CHP('model')
    chp = CanopyComplexity().fromLAS(las.simulateCloud()).CHP()
    print chp.zw

    
