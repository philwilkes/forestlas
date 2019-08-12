import numpy as np
import matplotlib.pyplot as plt
import math

class crownDim:
    
    def __init__(self, verbose=False, plot=False):
        
        self.verbose = verbose
        self.plot = plot
        
        if plot is not False:
            self.a = plot
            
    def genus_params(self, allometric, dbh, height=None, genus=None, species=None, height_class=None, position="none"):
        
        self.genus = genus
        self.position = position
        
        if allometric == 'scanlan':
            ### Scanlan 2010 ###
            self.height = 15.868 + (0.411 * dbh) - 0.003 * (dbh - 53.346)**2
            self.base = 13.622 + (0.255 * dbh) - 0.002 * (dbh - 53.346)**2
            self.radius = 1.966 + (0.115 * dbh)
        
        elif allometric == 'rushworth':        
            if species == "polyanthemos" and position == 1:
                self.height = 6.73589841952 + 0.199401565137 * dbh
                self.radius = 2.44438391898 + 0.0968236918828 * dbh
                self.base = 5.64863434095 + 0.0508669633267 * dbh
            elif species == "polyanthemos" and position == 2:
                self.height = 8.03659852077 + 0.0980603043708 * dbh
                self.radius = 1.53665512596 + 0.153832682204 * dbh
                self.base = 4.77895559438 + -0.0501664079081 * dbh
            elif species == "melliodora" and position == 1:
                self.height = 13.1074336283 + 0.181415929204 * dbh
                self.radius = -0.0258262268705 + 0.261383748994 * dbh
                self.base = -3.84514883347 + 0.421962992759 * dbh
            elif species == "melliodora" and position == 2:
                self.height = 4.62598991008 + 0.592333644974 * dbh
                self.radius = -2.06257925378 + 0.404778261046 * dbh
                self.base = 4.13866554657 + 0.0782895004268 * dbh
            elif species == "tricarpa" and position == 1:
                self.height = 11.2973396649 + 0.118405746654 * dbh
                self.radius = 1.90859934333 + 0.11935877297 * dbh
                self.base = 7.15742290658 + -0.0309056611952 * dbh
            elif species == "tricarpa" and position == 2:
                self.height = 5.05200102811 + 0.277243417286 * dbh
                self.radius = -0.306595001588 + 0.223036436419 * dbh
                self.base =  5.20539453592 + -0.00489743984257 * dbh
            elif species == "microcarpa": 
                self.height = 6.14528846053 + 0.541168944217 * dbh
                self.radius = 0.937658231167 + 0.180962968401 * dbh
                self.base = 1.10902140512 + 0.3315224801 * dbh
            elif species == "macrorrhyncha":
                self.height = 3.56849545646 + 0.489722695707 * dbh
                self.radius = 0.192442248096 + 0.16049402261 * dbh
                self.base = 4.13546317166 + 0.10044532714 * dbh
            elif position == 1:
                self.height = 3.05768707007 + 3.75103882462 * np.log(dbh)
                self.radius = -8.08923000772 + 4.14669504717 * np.log(dbh)
                self.base = 7.26517208273 + -0.0301284879315 * dbh            
            elif position == 2:
                self.height = 0.481343823925 + 3.74241504451 * np.log(dbh)
                self.radius = -6.32450180289 + 3.57228666209 * np.log(dbh)
                self.base = 5.34980523765 + -0.0219543192422 * dbh
            else:
                print species
            
            #self.base = 0.616759536617 + 3.07542949969 * np.log(self.height / 2.42)            

            ### Rushworth field work
            #self.height = 5.59 * np.log(0.48 * dbh)
            #self.base = 4.91 * np.log(0.22 * dbh)
            #self.radius = 3.28 * np.log(0.15 * dbh)
            
        elif allometric == 'field':
            if "eucalypt" in genus.lower():
                if height_class == 0 and position <= 1:  # low height and dominant
                    self.height = 1.86064881159 + 4.33064331125 * np.log(dbh) # n:306, rmse: 1.45, r2: 0.68
                    self.radius = -4.40385183196 + 2.99180690117 * np.log(dbh) # n:196, rmse: 0.51, r2: 0.84
                    self.base = 3.77967694134 + 0.293670137086 * self.height # n:169, rmse: 1.01, r2: 0.22
                if height_class == 1 and position <= 1:  # med height and dominant
                    self.height = -6.3614102534 + 11.5476586081 * np.log(dbh) # n:285, rmse: 3.98, r2: 0.40
                    self.radius = -30.9357885227 + 9.7432187604 * np.log(dbh)# n:272, rmse: 1.40, r2: 0.75
                    self.base = 11.0629490818 + 0.410104416904 * self.height # n:233, rmse: 3.08, r2: 0.39
                # if height_class == 2 and position == 1:  # tall height and dominant (all trees)
                   # self.height = -14.8638124803 + 16.0836040014 * np.log(dbh) # n:22, rmse: 5.78, r2: 0.50
                   # self.radius = -30.9357885227 + 9.7432187604 * np.log(dbh) # combined with height class >0 as n = 2
                   # self.base = 11.0629490818 + 0.410104416904 * self.height # ditto   
                # if height_class == 2 and position == 1:  # tall height and dominant (regnans only)
                    # self.height = -51.6987300663 + 22.503790377 * np.log(dbh) # n:22, rmse: 5.78, r2: 0.50
                    # self.radius = -17.7715368458 + 7.04971631466 * np.log(dbh) # combined with height class >0 as n = 2
                    # self.base = 9.83446096282 + 0.285250808259 * self.height # ditto
                if height_class == 2 and position == 1:  # tall height and dominant (regnans only)
                    self.height = height # n:22, rmse: 5.78, r2: 0.50
                    self.radius = -17.7715368458 + 7.04971631466 * np.log(dbh) # combined with height class >0 as n = 2
                    self.base = 0.6 * height # ditto
                if height_class == 0 and position == 2:  # low height and suppressed
                    self.height = -2.20552357108 + 4.98820033612 * np.log(dbh) # n:246, rmse: 1.77, r2: 0.39
                    self.radius = -6.02277337577 + 3.42232555079 * np.log(dbh) # n:132, rmse: 0.76, r2: 0.48
                    self.base = 0.84309744471 + 0.39574379931 * self.height # n:132, rmse: 0.99, r2: 0.35
                if height_class > 0 and position == 2:  # med height and suppressed
                    self.height = -10.3155403958 + 9.37056330951 * np.log(dbh) # n:23, rmse: 2.51, r2: 0.70
                    self.radius = -6.02277337577 + 3.42232555079 * np.log(dbh) # combined with height class = 0 as n = 2
                    self.base = 0.84309744471 + 0.39574379931 * self.height # ditto
            else:  # non-eucalypt
                if height_class == 0:  # tall height and dominant
                    self.height = -7.25044287906 + 6.92940716229 * np.log(dbh) # n:97, rmse: 2.07, r2: 0.55
                    self.radius = -7.42177845525 + 3.955374589 * np.log(dbh) # n:77, rmse: 0.63, r2: 0.75
                    self.base = 1.99951178234 + 0.459496058702 * self.height # n:77, rmse: 1.00, r2: 0.69                   
                if height_class == 1:  # low height and suppressed
                    self.height = -11.4257385832 + 8.43514568771 * np.log(dbh) # n:86, rmse: 1.58, r2: 0.73
                    self.radius = -4.94874693431 + 3.05793566579 * np.log(dbh) # n:11, rmse: 0.56, r2: 0.69
                    self.base = 2.02235653003 + 0.485452812202 * self.height # n:11, rmse: 0.77, r2: 0.54
                if height_class == 2:  # med height and suppressed
                    self.height = -13.944895656 + 9.01486400668 * np.log(dbh) # n:75, rmse: 1.53, r2: 0.82
                    self.radius = -10.1429095331 + 5.06708994335 * np.log(dbh) # n:32, rmse: 0.54, r2: 0.72
                    self.base = 2.37924762553 + 0.421325043146 * self.height # n:32, rmse: 0.49, r2: 0.75
        
            self.height = height
            
            if height > 50 or species == "regnans":
                self.base = (1-0.29) * height
            elif height > 10:
                self.base = 0.55 * height
            else: self.base = 0.7 * height
            #print dbh, self.height, self.radius, self.base, height_class, position
        
        else:   raise NameError("allometric relationship keyword not understood")
        
        #if self.base > self.height:
        #    self.height, self.base = self.base, self.height
        
        if self.verbose:
            print "dbh: {}, h: {}, r: {}, b: {}, a:{}".format(dbh, self.height, self.radius, self.base, xxx)

        return self
    
    def ellipsoid_volume(self):
        
        crown_depth = self.height - self.base
        crown_centre = self.height - (crown_depth / 2)
        a = crown_depth / 2 # semimajor self.ais
        b = self.radius / 2 # semiminor self.ais
        
        L = {}
                
        if self.plot is not False:
            for t in range(0, 360):
                x = a*math.sin(math.radians(t))
                y = b*math.cos(math.radians(t))
                self.a.scatter(y, x + crown_centre, s=1, c="g", edgecolor="none")
                #if "ucalypt" in self.genus:
                #    x = (a * .75) * math.sin(math.radians(t))
                #    y = (b * .75) * math.cos(math.radians(t))
                #    self.a.scatter(y, x + (a / 1.375) + self.base, s=1, c="r", edgecolor="none")
        
        eccentricity = np.sqrt(b**2./a**2.)
        hyp = a * eccentricity
        bins = np.linspace(-np.floor(a), np.floor(a), np.floor(a) * 2 + 1)
        
        for i, opp in enumerate(bins):
            angle = math.asin((opp * eccentricity) / hyp)
            adj = math.cos(angle) * hyp
            area = np.pi * adj**2
            # calculates "hollow" eucalypt crowns
            if "ucalypt" in self.genus:
                j = float(i) / len(bins)
                if j <= 0.75:
                    factor = self.exp(j, 0.26229889, -1.52996682)
                    area = area * factor
                # calculates suppressed trees
                if self.position == 2:
                    area *= 0.5
            if "othofagus" in self.genus:
                area *= 1.2
            opp += crown_centre
            L[opp] = area

        return L
    
    def exp(self, x, a, b):
        return a * np.exp(-b * x)