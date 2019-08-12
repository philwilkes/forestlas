import numpy as np
import os
from forestlas.lasIO import *

# Wilkes et al. 2015. Understanding the effects of ALS pulse density 
# for metric retrieval across diverse forest types. PE&RS

class thinCloud:
    
    def __init__(self, las, ppsm):
    
        if isinstance(las, str):
            las = lasIO(las).all().asArray()
        
        self.ppsm = ppsm
        self.las = las
        self.xy = 1. / np.sqrt(ppsm)
        
        self.xmin = np.round(self.las['x'].min())
        self.xmax = np.round(self.las['x'].max())
        self.ymin = np.round(self.las['y'].min())
        self.ymax = np.round(self.las['y'].max())
        
        
    def thinned(self, Ex=0, Ey=0, fsw=1):
    
        """
        fsw = first return search window, it is best to set this to the number of
        iterations but sometimes this is to small.
        
        RETURNS self.thinnedCloud POINT CLOUD
        """
        
        fsw = (self.xy / 2 ) / fsw

        self.thinnedCloud = np.array([], dtype=self.las.dtype)
        first = self.las[self.las['rtn_num'] == 1]
        other = self.las[self.las['rtn_num'] != 1]
        self.sw = (self.las['z'].max() - self.las['z'].min()) * np.tan(np.deg2rad(5))

        if self.sw == 0: self.sw = 1

        for ii, i in enumerate(np.linspace(self.xmin, self.xmax, 
                                          (self.xmax - self.xmin) / self.xy + 1) + Ex):
            for jj, j in enumerate(np.linspace(self.ymin, self.ymax, 
                                              (self.ymax - self.ymin) / self.xy + 1) + Ey):

                f = first[(first['x'] > (i - fsw)) & 
                          (first['x'] <= (i + fsw)) &
                          (first['y'] > (j - fsw)) & 
                          (first['y'] <= (j + fsw))]
                o = other[(other['x'] > (i - self.sw)) & (other['x'] <= (i + self.sw)) &
                          (other['y'] > (j - self.sw)) & (other['y'] <= (j + self.sw))]
                
                if len(f) > 0:
                    fx = f['x'] - i
                    fy = f['y'] - j
                    distance_m = np.hypot(fx, fy)
                    idx = [distance_m == distance_m.min()]
                    f = f[idx][0]
                    rtn_tot = f["rtn_tot"]
                    self.thinnedCloud = np.hstack([self.thinnedCloud, f])
                    if rtn_tot > 1:
                        for rtn_num in range(2, rtn_tot + 1):
                            rn = o[o["rtn_num"] == rtn_num]
                            if len(rn) > 0:
                                rx = rn['x'] - f['x']
                                ry = rn['y'] - f['y']
                                distance_m = np.hypot(rx, ry)
                                idx = [distance_m == distance_m.min()]
                                rn = rn[idx][0]
                                self.thinnedCloud = np.hstack([self.thinnedCloud, rn])
        
        return self.thinnedCloud
