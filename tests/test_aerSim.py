'''
Created on 25 mars 2014

@author: colley
'''
import lsst.sims.atmosphere.transmission.aerSim as aSol
import numpy as np
import matplotlib.pyplot as pl


def test_getAerParametersFromFile(pFile):
    aCoef = aSol.getAerParametersFromFile(pFile)
    aSol.S_PlotLevel = 4
    wl = np.linspace(300, 900, 1000)
    aSol.getAerTransmittance(wl, aCoef[10, :3], 1.0)
    aSol.getAerTransmittance(wl, aCoef[11, :3], 1.0)
    aSol.getAerTransmittance(wl, aCoef[12, :3], 1.0)
    
    
#
# TEST
#
S_PathFileSimu = "../testSImuAero.txt"

test_getAerParametersFromFile(S_PathFileSimu) 

pl.show()

###########
# TESTING #
###########

# nval = 30
# imin = 650
# imax = imin + nval
# outdict = main()
# am = '0'
# numesh = np.arange(NU_MIN, NU_MAX, NU_STEP)
# lsstwl = 1.e7 / numesh  # convert frequencies (cm-1)into wvlength (nm)
# lwv = np.log(lsstwl)

# vcur = np.vander(lwv, P_DEG + 1)

# testarr = np.zeros((nval, len(aerstr)))
# for iwl, wl in enumerate(aerstr):
#     # testarr[:, iwl] = outdict[wl][am][imin:imax]
#     testarr[:, iwl] = outdict[wl][am]

# Construct the Vandermonde matrix for standard deviation
# vx = np.vander(laer_wl, P_DEG + 1)
# aer_wl = np.exp(laer_wl)

# col = ['b', 'r', 'm', 'g', 'c']
# xplot = np.arange(N_DAYS - 1)
# xplot2 = np.arange(0, N_DAYS - 1, 0.5)

# f, ax = plt.subplots(3, 2, sharex=True)
# pf = np.zeros((nval, 3))
# for day in xrange(nval):
#     lvect = testarr[day]
#     # Fit
#     pfit = np.polyfit(laer_wl, lvect, P_DEG)
#     pf[day] = pfit

# p0 = pf[:, 0]
# p1 = pf[:, 1]
# p2 = pf[:, 2]
# p0sp = splrep(xplot, p0)
# p1sp = splrep(xplot, p1)
# p2sp = splrep(xplot, p2)
# p02 = splev(xplot2, p0sp)
# p12 = splev(xplot2, p1sp)
# p22 = splev(xplot2, p2sp)

# ax[0, 0].plot(p0, 'r,')
# ax[1, 0].plot(p1, 'b,')
# ax[2, 0].plot(p2, 'g,')
# ax[0, 1].plot((p02-p0)/p0, 'r,')
# ax[1, 1].plot((p12-p1)/p1, 'b,')
# ax[2, 1].plot((p22-p2)/p2, 'g,')

# plt.show()
