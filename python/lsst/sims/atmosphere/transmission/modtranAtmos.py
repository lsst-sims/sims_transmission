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
import numpy.fft as npf
import numpy.random as npr
import scipy.interpolate as intp
import os

import MJDtools
import modtranTools

# fixed meteo sequence step to 1/100 of a day = 14.4 min
_mjd_step = 0.01
_default_seed = 10
_south_hem = ['Summer', 'Spring', 'Winter', 'Fall']
_north_hem = ['Winter', 'Fall', 'Summer', 'Spring']

class Atmosphere(object):
    """Base class for retrieving a set of simulated MODTRAN
    atmosphere parameters between two dates in mjd.
    """
    def __init__(self, mjds, mjde, npoints=None, seed=_default_seed):
        """Initialize class"""
        self.mjds = float(mjds)
        self.mjde = float(mjde)
        self.npoints = npoints
        self.seed = seed

        self._initMjdArray()
        
    def _initMjdArray(self):
        """Create mjd array"""
        mjds = round(self.mjds, 2)
        mjde = round(self.mjde, 2) + _mjd_step
        if self.npoints:
            self.mjd_step = (mjde - mjds) / self.npoints
        else:
            self.mjd_step = _mjd_step
        self.mjd_arr = np.arange(mjds, mjde, self.mjd_step)
        self.mjds = mjds
        self.mjde = mjde

    def _init_o3(self):
        """Initialize the atmospheric ozone sequence"""
        self._o3_arr = np.zeros_like(self.mjd_arr)
        data = self._simulate_o3()
        ids = MJDtools.MJDtoindex(self.mjds)
        data_scaled = data[ids:]
        mjd_scaled = np.arange(data_scaled.shape[0]) + int(self.mjds)
        data_spl = intp.UnivariateSpline(mjd_scaled, data_scaled)
        self._o3_arr = data_spl(self.mjd_arr)

    def _simulate_o3(self):
        """From satellite datasets, retrieve the main temporal
        variability and simulate 8 years of daily data
        """
        transDir = '/Users/alexandreboucaud/work/LSST/calib/LSST_svn/transmission/trunk/python/lsst/sims/atmosphere/transmission'
        #transDir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        o3avg_file = os.path.join(transDir, 'datafiles/ozone8Yavg.npy')
        #o3_file = os.path.join(transDir, 'datafiles/ozone8Ymasked.dat')
        o3flat_file = os.path.join(transDir, 'datafiles/ozone8Yflatmasked.dat')
        # Get a smooth seasonal variation
        o3avg = np.load(o3avg_file)
        o3smooth = modtranTools.movingAverage(o3avg - np.median(o3avg), 60)
        o3seasSpl = intp.InterpolatedUnivariateSpline(np.arange(365), o3smooth) 

        # Data corrected from seasonal var.
        flat_data = np.ma.load(o3flat_file)
        mean_data = np.mean(flat_data)
        data = flat_data - mean_data

        # FFT
        ndays = int(len(data))
        fft = npf.fft(data)
        # Random
        mu, sig = 0, np.abs(fft)
        npr.seed(self.seed)
        Famp = npr.normal(mu,sig)
        npr.seed(self.seed)
        phi = npr.uniform(0, 2*np.pi, ndays)
        Fphase = np.cos(phi) + 1j*np.sin(phi)
        random_data = npf.ifft(Famp * Fphase)
        season = o3seasSpl(np.arange(ndays)%365)
        return np.real(random_data) + mean_data + season
        
    def _init_h2o(self):
        """Initialize the atmospheric water vapor sequence
        """
        self._h2o_arr = np.zeros_like(self.mjd_arr)

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
        tstart = npr.uniform(self.mjds, self.mjde - 60, n)
        tend = tstart + npr.uniform(0, 50)
        for i in xrange(n):
            ivulc += np.where((self.mjd_arr >= tstart[i]) &
                             (self.mjd_arr < tend[i]), 2, 0)
        ivulc[np.where(ivulc>2)] = 2
        self._ivulc_arr = ivulc

    def init_main_parameters(self, nvulc = 5):
        """Initialize the main atmospheric parameters"""
        self._init_o3()
        self._init_h2o()
        self._init_aer()
        self._init_vulc(nvulc)
        
    def ozone(self, idx):
        """Atmosperic ozone"""
        return self._o3_arr[idx]

    def vapor(self, idx):
        """Atmospheric water vapor"""
        return self._h2o_arr[idx]
        
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
        """Volcanic aerosols"""
        return self._ivulc_arr[idx]

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


if __name__=="__main__":
    atm = Atmosphere(5538, 5539.5)
    atm.init_main_parameters(nvulc = 0)
    
