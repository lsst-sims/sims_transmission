#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
auxtelSequence.py

Class used to generate sequences of atmospheric parameters for auxiliary telescope.

Input data:
----------
A list of dictionaries sorted by date, with 4 specific keys:
    'ID'    the sequence id
    'MJD'   the date in mjd format
    'ZANG'  the zenith angle of this pointing
    'TGRAY' the gray extinction factor

Calling method:
--------------
>>> atmoseq = AtmosphereSequence(auxtel_pointing)
>>> atmoseq.generateParameters(seed)

=>  generates parameters for the pointings with the provided seed

>>> seq.getAtmTrans('tmp', save=True)

=> computes the atmospheric extinction in $MODTRAN_DATADIR/tmp/
=> outputfiles in the form tmp.*
=> saves the final transmission data into $MODTRAN_DATADIR/tmp/tmp_final.plt

Description :
-----------
The AtmosphereSequence
- inputs a live opsim dictionary / list of dictionaries OR a filename,
- verifies the data contains at least the following keys :
'obsHistID', 'expMJD', 'fieldRA', 'fieldDEC',
- initializes a pointing sequence by defining the date boundaries,
- calls the Atmosphere class in the modtranAtmos.py and generates atmospheric
parameters for the given pointings,
- writes those parameters in a file / database (using the pickle module),
- loads parameters from database on demand,
- runs MODTRAN without the aerosols for each selected pointing (one at a time)
and returns output,
- computes the aerosol extinction spectrum for the same pointings,
- combines both spectra to output the final atmosphere transmission either
live or saved into a file.

"""

import numpy as np
import os
import copy

from lsst.sims.atmosphere.transmission.modtranAtmos import Atmosphere
from lsst.sims.atmosphere.transmission.modtranCardsNoAer import ModtranCards
import lsst.sims.atmosphere.transmission.modtranTools as modtranTools
import lsst.sims.atmosphere.transmission.MJDtools as MJDtools
import lsst.sims.atmosphere.transmission.aerSim as asl

S_Verbose = 1

RAD2DEC = 180. / np.pi
DEG2RAD = np.pi / 180.

_opsim_keys = ['obsHistID', 'expMJD', 'fieldRA', 'fieldDEC']
#_modtran_keys = ['O3', 'H2O', 'IHAZE', 'ISEAS', 'IVULC', 'ICLD', 'IVSA', 'VIS',
#               'ZAER11', 'ZAER12', 'SCALE1', 'ZAER21', 'ZAER22', 'SCALE2']
_default_seed = 10


class AuxTelSequence(object):
    """This class is used to generate sequences of atmospheric parameters,
    matching long-term weather trends, and accounting for the airmass/pointing
    location of each observation. These atmospheric parameters are the
    basic inputs for MODTRAN, so can be used to generate an atmospheric
    transmission curve.
    The class generates these parameters for a series of opsim pointings.
    """
    def __init__(self, seq_list):
        """Instantiate an AuxTelSequence."""
        # List of SORTED sequence dictionaries
        self.visits = seq_list
        # List of MODTRAN dictionaries
        self.modtran_visits = []
        # List of aerosol parameters
        self.aerosol_visits = []
        # MODTRAN wavelengths array
        self.modtran_wl = []
        # Transmittance array
        self.transmittance = []
        # aero transmittance
        self.lAeroTrans = []

    def changeList(self, newlist):
        """Input a new sequence list."""
        self.visits = newlist

    def initPointingSequence(self):
        """Initialize input parameters for the atmospheric simulation."""
        # Number of visits
        self.npoints = len(self.visits)
        # Starting date
        self.mjds = self.visits[0]['MJD']
        started_year =  MJDtools.fromMJD(self.mjds)[0]        
        if (self.mjds - MJDtools.toMJD(started_year,12, 21)) < 0 :
            started_year -= 1
        self.offsetMJD = MJDtools.toMJD(started_year, 12, 21)


    def generateParameters(self, seed=_default_seed):
        """Generate the atmospheric parameters over time.
        Returns a list of dictionaries containing the modtran
        information for each opsim visit (in the same order).
        """
        self.initPointingSequence()
        # Instantiate the Atmosphere class
        self.atmos = Atmosphere(seed)
        # Generate main atmosphere parameters sequence
        self.atmos.init_main_parameters()
        # Associate a value of these parameters for each pointing
        for visit_dict in self.visits:
            # Get ID and date
            vis_id = visit_dict['ID']
            mjd = visit_dict['MJD']
            z_angle = visit_dict['ZANG']
            # Get atmosphere parameters
            modtran_dict = self.fillModtranDictionary(mjd-self.offsetMJD, vis_id, z_angle)
            self.modtran_visits.append(modtran_dict)
            # used modulo to do interpolation  no extrapolation !
            # may be not coherent with aerosol saison ? (JMC)    
            aerosolOnly = self.atmos.aerosols(mjd % asl.N_DAYS)
            self.aerosol_visits.append(aerosolOnly + (z_angle,))
            if S_Verbose==1:
                print "generateParameters mjd",mjd
                print aerosolOnly, z_angle, self.aerosol_visits[-1]
        # Init transmission array
        self.initTransmissionArray(len(self.modtran_visits))
        print 'generateParameters FIN'

    def fillModtranDictionary(self, mjd, obsid, z_angle):
        """Return a dictionary filled with all Modtran parameters
        for a given visit."""
        if not self.atmos:
            raise ValueError('Atmosphere class not called')
        mdict = {}
        mdict['ID'] = obsid
        mdict['ZANGLE'] = z_angle
        mdict['MODEL'] = self.atmos.model(mjd)
        mdict['O3'] = self.atmos.ozone(mjd)
        mdict['H2O'] = self.atmos.wvapor(mjd)
        # put aerosols computation outsite of MODTRAN
        mdict['IHAZE'] = self.atmos.ihaze()
        return mdict

    def initTransmissionArray(self, nruns):
        """Initialize class attribute for atmospheric transmission"""
        if len(self.modtran_wl) == 0:
            self.initModtranWavelengths()

        self.transmittance = np.zeros((nruns, len(self.modtran_wl)))

    def initModtranWavelengths(self):
        """Load in memory the sorted tabulated wavelengths in nanometers at which
        MODTRAN outputs data."""
        main_dir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        modtranwfile = os.path.join(main_dir, 'data/modtranwl.txt')
        self.modtran_wl = np.loadtxt(modtranwfile)

    def getAtmTrans(self, outfile='tmp', modversion=4, save=True):
        """Go through the process of getting an atmospheric transmission"""
        # Set the name of MODTRAN output file
        self.outfilename = outfile
        # Init the modtran card class
        self.base_modcard = self.initModtran(modversion)
        # Loop over the MODTRAN runs
        for run in xrange(len(self.visits)):
            self.runModtran(run)
            # Read MODTRAN output
            modtrans = self.getModtranExtinction()
            # Compute aerosols transmission
            aertrans = self.getAerTransmittance(run)
            self.lAeroTrans.append(aertrans)
            print "getAerTransmittance(%d) mean %f, std %f"%(run, aertrans.mean(), aertrans.std()) 
            # Retrieve gray extinction
            t_gray = self.visits[run]['TGRAY']
            # Multiply both to get full transmission
            self.transmittance[run] = t_gray * modtrans * aertrans
            print t_gray,  modtrans, aertrans
        if save:
            self.saveTrans()
        # Data saved into a file named '%s_final.plt' % self.outfilename

    def initModtran(self, modversion):
        """Call the modtranCards class in order to run MODTRAN for
        selected visits"""
        #modcard.setDefaults()
        atmTransDir = os.getenv('ATMOSPHERE_TRANSMISSION_DIR')
        if not atmTransDir:
            raise Exception('If ATMOSPHERE_TRANSMISSION_DIR env not set, \
                    must specify template filename.')
        if modversion == 4:
            templatefile = os.path.join(atmTransDir, 'data/Cardtemplate.dat_michel')
            formatfile = os.path.join(atmTransDir, 'data/FormatParameters.dat_michel')
        elif modversion == 5:
            templatefile = os.path.join(atmTransDir, 'data/Cardtemplate.dat')
            formatfile = os.path.join(atmTransDir, 'data/FormatParameters.dat')
        else:
            raise ValueError('MODTRAN version not supported')
        # Call for the MODTRAN card class
        modcard = ModtranCards()
        modcard.readCardTemplate(templatefile)
        modcard.readParameterFormats(formatfile)

        return modcard

    def runModtran(self, run):
        """Call the modtranCards class in order to run MODTRAN for
        selected visits"""
        # Create a copy of the base modtran card
        modcard = copy.deepcopy(self.base_modcard)
        # Write the cards to the disk
        modcard.writeModtranCards(self.modtran_visits[run], self.outfilename)
        modcard.runModtran()

    def getModtranExtinction(self):
        """Read MODTRAN atmosphere extinction output file '.plt' and store the
        runs into a 2D array"""
        print  self.modtran_wl
        if len(self.modtran_wl)==0:
            self.initModtranWavelengths()

        modtranDataDir = os.getenv('MODTRAN_DATADIR')
        # MODTRAN transmission outputfile
        outputfile = '{0}/{0}.plt'.format(self.outfilename)
        outputpath = os.path.join(modtranDataDir, outputfile)
        # Initialize array for transmittance
        trans = np.zeros(len(self.modtran_wl))
        # print "outputpath",outputpath
        with open(outputpath, 'r') as outf:
            # File starts with data - no header
            # I use negative indices because MODTRAN prints the wavelengths
            # in decreasing order
            idx = len(self.modtran_wl) -1
            for line in outf:
                lElt = line.split()
                # starstwith not available in python 2.6.4
                if len(lElt) == 0 : continue
                if lElt[0][0] == '$': continue
                values = line.strip().split()
                if idx < 0:
                    print abs(idx) , len(self.modtran_wl)
                    print values[0]
                    raise ValueError("Too many values to unpack from MODTRAN \
                        outputfile.")
                trans[idx] = float(values[1])
                idx -= 1
        return trans
        # MODTRAN transmittance stored



    def getAerTransmittance(self, run):
        """
        Compute the atmospheric transmittance due to aerosols
        call aerSim module to compute transmittance
        """
        if self.modtran_wl == []:
            self.initModtranWavelengths()
        # Compute Vandermonde matrix
        # degree is hardcoded but only 2nd degree polynomial is used
        wlFit = self.modtran_wl
        if S_Verbose >=1: print "wlFit: ",wlFit
        #vdm_wl = np.vander(wlFit, 3)
        # Get polynomial roots as well as zenith angle
        p0, p1, p2, z_ang = self.aerosol_visits[run]
        if S_Verbose >=1: print "aerosol_visits[%d]:%s"%(run, str(self.aerosol_visits[run]))
        pfit = np.array([p0, p1, p2])
        # Retrieve airmass from zenith angle
        airmass = modtranTools.zenith2airmass(z_ang, site='lsst', unit='rad')
        if S_Verbose >=1: print "airmass: ",airmass
        aero = asl.getAerTransmittance(wlFit, pfit, airmass)
        #np.save("aero%d"%run, aero )        
        return aero 


    def saveTrans(self):
        """Save full extinction into an ascii file

        The structure is:
            wavelength transm[run1] transm[run2] transm[run3] etc.
        """
        if S_Verbose: print "array transmittance:",self.transmittance
        modtranDataDir = os.getenv('MODTRAN_DATADIR')
        outputfile = '{0}/{0}_final.plt'.format(self.outfilename)
        outputfileAero = '{0}/{0}_aero.plt'.format(self.outfilename)
        outputpath = os.path.join(modtranDataDir, outputfile)
        with open(outputpath, 'w') as transmf:
            transmf.write('# FINAL ATMOSPHERE TRANSMISSION\n')
            for val in xrange(len(self.modtran_wl)):
                data = '\t'.join('{0:f}'.format(self.transmittance[run][val])
                                for run in xrange(len(self.visits)))
                line = '{0}\t{1}\n'.format(self.modtran_wl[val], data)
                transmf.write(line)
        with open(outputfileAero, 'w') as transmf:
            transmf.write('# AEROSOL TRANSMISSION\n')
            for val in xrange(len(self.modtran_wl)):
                data = '\t'.join('{0:f}'.format(self.lAeroTrans[run][val])
                                for run in xrange(len(self.visits)))
                line = '{0}\t{1}\n'.format(self.modtran_wl[val], data)
                transmf.write(line)
        print "{0} MODTRAN computed spectra written in {1}".format(
            len(self.visits), outputpath)
