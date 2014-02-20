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


#
# TESTS
#

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
# ETUDES
#

def evolutionAerosol():
    """
    create plot image for movie
    """
    atmos = Atmosphere(10)
    atmos.init_main_parameters()
    coef = atmos.aerosols(100)
    print coef
    wl = numpy.linspace(300, 1000, 1000)
    asl.getAerTransmittance(wl, coef, 1.0)
    wl = np.linspace(300, 1000, 512)    
    pasDay = 120
    nbPlot = pasDay*24
    lday01 = np.linspace(0, pasDay, nbPlot, True)
    for d1,idx in zip(lday01, range(nbPlot)):
        print idx,d1
        plt.figure(1)
        plt.title("aerosol transmission")        
        [plt.plot(wl,asl.getAerTransmittance(wl, atmos.aerosols(d1+_li*pasDay) , 1.0)) for _li in range(3)]
        strDay = ["day %.2f"%(d1+_li*pasDay) for _li in range(3)]
        plt.ylim(ymin=0.8, ymax=1.0)
        plt.grid()
        plt.xlabel("nm")
        plt.legend(strDay, loc=4)
        plt.savefig("/home/colley/temp/lsst/movie/aero%04d.png"%idx)
        plt.close()


#
# MAIN
#

evolutionAerosol()
#test_aerosolData2()
#test_getAerTransmittance()
pl.show()
