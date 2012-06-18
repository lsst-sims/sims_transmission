"""
modtranAtmos.py
--
Contains the class Atmosphere that can produce simulated
realistic atmosphere parameters for MODTRAN
List of parameters that can be output :
o3                    - Ozone
h2o                   - Water vapor
vis0, visamp, visaz   - Aerosol parameters
ihaze                 - Aerosol driver
iseas                 - Season
ivulc                 - Volcanic aerosols
icld                  - Clouds
ivsa                  - Vertical structure algorithm
zaer11, zaer12        - Lower aerosol layer
scale1                - Scale for lower aerosol layer
zaer21, zaer22        - Higher aerosol layer
scale2                - Scale for higher aerosol layer
"""
import numpy as np
from numpy.random import uniform
from scipy.interpolate import splrep splev

import MJDtools

class Atmosphere(object):
    """Base class for retrieving a set of simulated MODTRAN
    atmosphere parameters between two dates in mjd."""
    
    # fixed meteo sequence step to 1/100 of a day = 14.4 min
    _mjd_step = 0.01
    _south_hem = ['Summer', 'Spring', 'Winter', 'Fall']
    _north_hem = ['Winter', 'Fall', 'Summer', 'Spring']
    
    def __init__(self, mjds, mjde, npoints=None):
        """Initialize class"""
        self.mjds = float(mjds)
        self.mjde = float(mjde)
        self.npoints = npoints

        self._initMjdArray()
        
    def _initMjdArray(self):
        """Create mjd array"""
        mjds = round(self.mjds, 2)
        mjde = round(self.mjde, 2) + _mjd_step
        if self.npoint:
            self.mjd_step = (mjde - mjds) / self.npoints
        else:
            self.mjd_step = _mjd_step
        self.mjd_arr = np.arange(mjds, mjde, self.mjd_step)
        self.mjds = mjds
        self.mjde = mjde

    def _init_o3(self):
        '''
        self._o3_spline = ...
        '''
        pass
    

    def _init_h2o(self):
        '''
        self._h2o_spline = ...
        '''
        pass

    def _init_aer(self):
        self.vis0 = np.zeros_like(self.mjd_arr)
        self.visaz = np.zeros_like(self.mjd_arr)
        self.visamp = np.zeros_like(self.mjd_arr)
        pass

    def _init_vulc(self, n):
        """Create n periods of fresh volcanic aerosol with
        random durations (0-50 days) """
        if (self.mjde-self.mjds) < n*50:
            raise ValueError('Too many volcanic periods')
        ivulc  = np.zeros_like(self.mjd_arr, dtype='int')
        tstart = uniform(self.mjds, self.mjde - 60, n)
        tend = tstart + uniform(0, 50)
        for i in xrange(n):
            ivulc += np.where((self.mjd_arr >= tstart[i]) &
                             (self.mjd_arr < tend[i]), 2, 0)
        ivulc[np.where(ivulc>2)] = 2
        self.ivulc_arr = ivulc

    def init_main_parameters(self, nvulc = 5):
        self._init_o3()
        self._init_h2o()
        self._init_aer()
        self._init_vulc(nvulc)
        
    def ozone(self):
        pass

    def vapor(self):
        pass

    def model(self, idx):
        """Seasonal model"""
        models = [2,6,3,6]
        return models[self.iseas(idx)]

    def ihaze(self):
        """Aerosols drivers
           0 : none
           1 : rural
           6 : tropospheric
        """
        return 1

    def iseas(self, idx, addName=None):
        """Season at Cerro Pachon site
           0 : Summer
           1 : Spring
           2 : Winter
           3 : Fall
        """
        iseas = MJDtools.getSeason(self.mjd_arr[idx])
        if addName:
            return iseas, self._south_hem[iseas]
        else:
            return iseas

    def ivulc(self, idx):
        """ """
        return ivulc_arr[idx]

    def ivsa(self):
        """Switch on vertical structure algorith"""
        return 1

    def zaer11(self):
        """Bottom lower aerosol layer"""
        return 0.0
    
    def zaer12(self):
        """Top lower aerosol layer"""
        return 3.0
            
    def scale1(self, idx):
        if not _initialized_aer:
            self.init_aer()
        return np.fix(-7.08 * np.log(self.vis0[idx]) + 29.)

    def zaer21(self):
        """Bottom higher aerosol layer"""
        return 0.0
    
    def zaer22(self):
        """Top higher aerosol layer"""
        return 3.0
        
    def scale2(self):
        return 1.0
