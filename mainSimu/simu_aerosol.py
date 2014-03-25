'''
Created on 25 mars 2014

@author: colley
'''


import lsst.sims.atmosphere.transmission.aerSim as aSol


def doFileAeroParameter(pFileOut, seed=1):
    oSim = aSol.SimuAeroStruct()
    aSol.fitAeroParameters(oSim, seed, '0', pFileOut)
    
    
#
# MAIN
#
S_PathFileSimu = "../testSImuAero.txt"

doFileAeroParameter(S_PathFileSimu)
