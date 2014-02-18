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
import numpy
import scipy
import os

import lsst.sims.atmosphere.transmission.MJDtools as MJDtools
import lsst.sims.atmosphere.transmission.modtranTools as modtranTools
import lsst.sims.atmosphere.transmission.aerSim as asl


# Maximum number of simulated days
# Currently limited by water vapor with 7 years
_max_sim_days = 365. * 7
_default_seed = 10
_south_hem = ['Summer', 'Spring', 'Winter', 'Fall']
_north_hem = ['Winter', 'Fall', 'Summer', 'Spring']


class Atmosphere(object):
    """Base class for retrieving a set of simulated MODTRAN
    atmosphere parameters between two dates in mjd.
    """
    def __init__(self, seed=_default_seed):
        
        """Initialize Atmosphere class"""
        # Starting/ending date        
        # Seed for the randoms
        self.seed = seed
        # Ozone spline
        self._o3_spl = None
        # Water vapor spline
        self._h2o_spl = None
        # Aerosols parameter splines
        self._aer_p0_spl = None
        self._aer_p1_spl = None
        self._aer_p2_spl = None
        self.started_year = 0 

    # def _init_mjdArray(self):
    #     """Create mjd array to serve as base for splines"""
    #     dmjd = float(self.mjde - self.mjds)
    #     mjd_step = dmjd / self.npoints
    #     if mjd_step < _mjd_minstep:
    #         mjd_step = _mjd_minstep
    #     nmjd = dmjd / float(mjd_step)
    #     if nmjd < _mjd_nmin:
    #         nmjd = _mjd_nmin
    #     self._mjd_arr = numpy.linspace(self.mjds, self.mjde, nmjd)

    def _init_o3(self):
        """Initialize the atmospheric ozone sequence"""
        # self._o3_arr = numpy.zeros_like(self._mjd_arr)
        #ids = self.MJDtoindex(self.mjds, seas=True)
        #ide = self.MJDtoindex(self.mjde, seas=True)
        #print "index min max", ids, ide
        data = self._simulate_o3()
        # ids = MJDtools.MJDtoindex(self._mjd_arr[0], seas=True)
        # ide = MJDtools.MJDtoindex(self._mjd_arr[-1], seas=True)
        #data_scaled = data[ids:ide]
        mjd_scaled = numpy.arange(data.size)
        # MODTRAN SCALING
        data *= self._modtran_ozone_scalefactor(mjd_scaled)
        # data_spl = scipy.interpolate.UnivariateSpline(mjd_scaled, data_scaled)
        # self._o3_arr = data_spl(self._mjd_arr)
        self._o3_spl = scipy.interpolate.UnivariateSpline(
            mjd_scaled, data)

        

    def _simulate_o3(self):
        """From satellite datasets, retrieve the main temporal variability of
        ozone and simulate 8 years of daily data"""
        # transDir = '/Users/alexandreboucaud/work/LSST/calib/LSST_svn/\
        #     transmission/trunk/python/lsst/sims/atmosphere/transmission'
        transDir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        o3avg_file = os.path.join(transDir, 'data/ozone8Yavg.npy')
        #o3_file = os.path.join(transDir, 'data/ozone8Ymasked.dat')
        o3flat_file = os.path.join(transDir, 'data/ozone8Yflatmasked.dat')
        # Get a smooth seasonal variation
        o3avg = numpy.load(o3avg_file)
        o3smooth = modtranTools.movingAverage(o3avg - numpy.median(o3avg), 60)
        o3seasSpl = scipy.interpolate.InterpolatedUnivariateSpline(
            numpy.arange(365), o3smooth)
        # Data corrected from seasonal var.
        flat_data = numpy.ma.load(o3flat_file)
        mean_data = numpy.mean(flat_data)
        data = flat_data - mean_data
        # FFT
        ndays = int(len(data))
        fft = numpy.fft.fft(data)
        # Random
        mu, sig = 0, numpy.abs(fft)
        numpy.random.seed(self.seed)
        Famp = numpy.random.normal(mu, sig)
        numpy.random.seed(self.seed)
        phi = numpy.random.uniform(0, 2 * numpy.pi, ndays)
        Fphase = numpy.cos(phi) + 1j*numpy.sin(phi)
        random_data = numpy.fft.ifft(Famp * Fphase)
        season = o3seasSpl(numpy.arange(ndays) % 365)
        res = numpy.real(random_data) + mean_data + season
        return res

    def _modtran_ozone_scalefactor(self, mjd_arr):
        """Ozone scaling factor depends on MODTRAN seasonal model
        which depends on the season.
        Hence we can relate the ozone scaling factor straight to
        the season
        Finally it is normalized to 275 DU"""
        o3_sc_factor = numpy.array([0.8292, 0.8005, 0.7299, 0.8005])
        seas_idx = numpy.asarray([MJDtools.getSeason(mjd) for mjd in mjd_arr])
        sc_out = o3_sc_factor[seas_idx]
        return sc_out / 275.

    def _init_h2o(self):
        """Initialize the atmospheric water vapor sequence"""
        # of data pts a day => factor scl
        scl = 2
        #ids = self.MJDtoindex(self.mjds, seas=None)
        #ide = self.MJDtoindex(self.mjde, seas=None)
        #data_scaled = self._simulate_h2o(scl*ids, scl*ide)
        data_scaled = self._simulate_h2o()
        mjd_scaled = numpy.arange(data_scaled.size)/float(scl)
        data_scaled *= self._modtran_h2o_scalefactor(mjd_scaled)
        # data_spl = scipy.interpolate.UnivariateSpline(mjd_scaled, data_scaled)
        # self._h2o_arr = data_spl(self.mjd_arr)
        self._h2o_spl = scipy.interpolate.UnivariateSpline(
            mjd_scaled, data_scaled)

    def _simulate_h2o(self):
        """From satellite datasets, retrieve the main temporal variability of
        water vapor and simulate 7 years of data with 2 data points per day """
        # The computation is made in log
        # transDir = '/Users/alexandreboucaud/work/LSST/calib/LSST_svn/\
        #     transmission/trunk/python/lsst/sims/atmosphere/transmission'
        transDir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        wvavg_file = os.path.join(transDir, 'data/wv7Yavgmasked.dat')
        wvavg_stdfile = os.path.join(transDir, 'data/wv7Ystdmasked.dat')
        wvflat_file = os.path.join(transDir, 'data/wv7Yflatmasked.dat')
        # Get a smooth seasonal variation wv_mean(t)
        wvavg_ma = numpy.load(wvavg_file)
        wvsmooth = modtranTools.movingAverage(wvavg_ma, 60)
        wvseasSpl = scipy.interpolate.InterpolatedUnivariateSpline(
            numpy.arange(2 * 365), wvsmooth)
        # Get a smooth standard deviation sigma_wv(t)
        wvstd_ma = numpy.load(wvavg_stdfile)
        wvstdsmooth = modtranTools.movingAverage(wvstd_ma, 60)
        wvstdSpl = scipy.interpolate.InterpolatedUnivariateSpline(
            numpy.arange(2 * 365), wvstdsmooth)
        # Data corrected from seasonal variation and amplitude
        data = numpy.ma.load(wvflat_file)
        # FFT
        ndays = int(len(data))
        fft = numpy.fft.fft(data)
        # Random
        mu, sig = 0, numpy.abs(fft)
        numpy.random.seed(self.seed)
        Famp = numpy.random.normal(mu, sig*2.0)
        numpy.random.seed(self.seed)
        phi = numpy.random.uniform(0, 2*numpy.pi, ndays)
        Fphase = numpy.cos(phi) + 1j*numpy.sin(phi)
        random_data = numpy.fft.ifft(Famp * Fphase)
        season_range = wvseasSpl(numpy.arange(ndays) % (2 * 365))
        std_range = wvstdSpl(numpy.arange(ndays) % (2 * 365))
        res = numpy.exp(numpy.real(random_data) * std_range + season_range)
        return res

    def _modtran_h2o_scalefactor(self, mjd_arr):
        """Water vapor scaling factor depends on MODTRAN seasonal model
        which depends on the season. Hence we can relate it straight to
        the season."""
        h2o_sc_factor = numpy.array([2.5396, 1.4164, 0.8534, 1.4164])
        seas_idx = numpy.asarray([MJDtools.getSeason(mjd) for mjd in mjd_arr])
        sc_out = h2o_sc_factor[seas_idx]
        return sc_out / 449.23

    def _init_aer(self):
        """
        for each polynomial coefficient fit with spline ie spline scipy object
        """
        data = self._simulate_aer()
        print data.shape
        print data[99:102]
        data_p0, data_p1, data_p2, stdev = data.T
        mjd_scaled = numpy.arange(data_p0.size)
        print "_init_aer:  Start ", data_p0.size
        print mjd_scaled
        
        self._aer_p0_spl = scipy.interpolate.UnivariateSpline(
            mjd_scaled, data_p0, s=0)
        self._aer_p1_spl = scipy.interpolate.UnivariateSpline(
            mjd_scaled, data_p1, s=0)
        self._aer_p2_spl = scipy.interpolate.UnivariateSpline(
            mjd_scaled, data_p2, s=0)
        print "_init_aer:  End"
        
        
    def _simulate_aer(self):
        """Simulate the aerosol optical depth as a function of wavelength
        and return 2nd degree polynomial fitting parameters
        (cf. aerSim.py for more info)"""
        print "_simulate_aer: start"
        data = asl.getAerParameters(self.seed, airmass='0')
        print "_simulate_aer: end"
        arrdata = numpy.array(data)
        return arrdata

    def init_main_parameters(self, nvulc=0):
        """Initialize the main atmospheric parameters"""
        # if not self._initialized_array:
            # self._init_mjdArray()
        self._init_o3()
        self._init_h2o()
        self._init_aer()
        # self._init_vulc(nvulc)

    def ozone(self, mjd):
        """Return atmosperic ozone for a given mjd"""
        # return self._o3_arr[idx]
        return self._o3_spl(mjd)

    def wvapor(self, mjd):
        """Return atmospheric water vapor for a given mjd"""
        # return self._h2o_arr[idx]
        return self._h2o_spl(mjd)

    def aerosols(self, mjd):
        """Return atmospheric aerosols parameters for a given mjd"""
        ret = numpy.array([self._aer_p0_spl(mjd),
                self._aer_p1_spl(mjd),
                self._aer_p2_spl(mjd)])
        return ret.ravel()

    def model(self, mjd):
        """Seasonal model"""
        models = [2, 6, 3, 6]
        # return models[self.iseas(idx)]
        return models[MJDtools.getSeason(mjd)]

    def ihaze(self):
        """Aerosols drivers
           0 : none
           1 : rural
           6 : tropospheric
        """
        return 0

    def _init_vulc(self, n):
        """Create n periods of fresh volcanic aerosol with
        random durations (0-50 days) """
        if (self.mjde-self.mjds) < n*50:
            raise ValueError('Too many volcanic periods')
        ivulc = numpy.zeros_like(self.mjd_arr, dtype='int')
        tstart = numpy.random.uniform(self.mjds, self.mjde - 60, n)
        tend = tstart + numpy.random.uniform(0, 50)
        for i in xrange(n):
            ivulc += numpy.where((self.mjd_arr >= tstart[i]) &
                                (self.mjd_arr < tend[i]), 2, 0)
        ivulc[numpy.where(ivulc > 2)] = 2
        self._ivulc_arr = ivulc

    def ivulc(self, idx):
        """Volcanic aerosols"""
        return self._ivulc_arr[idx]

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
        if not self._initialized_aer:
            self.init_aer()
        return numpy.fix(-7.08 * numpy.log(self.vis0[idx]) + 29.)

    def zaer21(self):
        """Bottom higher aerosol layer"""
        return 0.0

    def zaer22(self):
        """Top higher aerosol layer"""
        return 3.0

    def scale2(self):
        return 1.0


if __name__ == "__main__":
    print "by main "
    atm = Atmosphere(12)
    atm.init_main_parameters(nvulc=0)
