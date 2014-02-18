# -*- coding: utf-8 -*-
'''
Created on 11 févr. 2014

@author: colley
'''
import numpy as np
import matplotlib.pyplot as plt

#np.random.seed(10)

np.random.seed(620)


import matplotlib.pyplot as pl
from lsst.sims.atmosphere.transmission.modtranAtmos import *



def test_aerosolData():
    atmos = Atmosphere(10)
    atmos.init_main_parameters()
    coef = atmos.aerosols(100)
    print coef
    wl = numpy.linspace(300, 1000, 1000)
    asl.getAerTransmittance(wl, coef, 1.0)
    wl = np.linspace(300, 1000, 512)
    lday = [100,101, 200, 300, 400]
    lTrans = [asl.getAerTransmittance(wl, atmos.aerosols(_d) , 1.0) for _d in lday]
    print lTrans[0]
    plt.figure()
    plt.title("aerosol transmission from spline")
    [plt.plot(wl, tr,'*') for tr in lTrans]
    strDay = ["day "+str(_d) for _d in lday]
    plt.grid()
    plt.xlabel("nm")
    plt.legend(strDay, loc="best")


def test_aerosolData2():
    atmos = Atmosphere(10)
    atmos.init_main_parameters()
    coef = atmos.aerosols(100)
    print coef
    wl = numpy.linspace(300, 1000, 1000)
    asl.getAerTransmittance(wl, coef, 1.0)
    wl = np.linspace(300, 1000, 512)
    lday = np.array([100 ,101, 200, 300, 400], dtype=np.float32)
    lday += 0.5
    print lday
    lTrans = [asl.getAerTransmittance(wl, atmos.aerosols(_d) , 1.0) for _d in lday]
    print lTrans[0]
    plt.figure()
    plt.title("aerosol transmission from spline")
    [plt.plot(wl, tr) for tr in lTrans]
    strDay = ["day "+str(_d) for _d in lday]
    plt.grid()
    plt.xlabel("nm")
    plt.legend(strDay, loc="best")

def test_getAerTransmittance():
    coef = numpy.array([  1.65058424e+00,  -2.15406291e+01 ,  6.63529808e+01])
    coef = numpy.array([  1.5,  -2.15406291e+01 ,  6.63529808e+01])
    wl = numpy.linspace(300, 1000, 1000)
    asl.getAerTransmittance(wl, coef, 1.0)

#
# MAIN
#


test_aerosolData2()
#test_getAerTransmittance()
pl.show()
