'''
Created on 25 mars 2014

@author: colley
'''
import lsst.sims.atmosphere.transmission.aerSim as aSol
import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as spi
import scipy.signal as sps

aSol.S_PlotLevel = 2


def smooth(sig, nWindow):
    """
    lissage d'un signal par un fenetre gaussienne de nWindow echantillons
    """
    w = sps.gaussian(nWindow,nWindow/10)
    s = w.sum()
    a = sps.fftconvolve(sig, w/s)
    print len(sig)
    print len(a)
    return a[nWindow/2-1:-nWindow/2]

#
# TEST
#

def test_getAerParametersFromFile(pFile):
    aCoef = aSol.getAerParametersFromFile(pFile)
    aSol.S_PlotLevel = 4
    wl = np.linspace(300, 900, 1000)
    aSol.getAerTransmittance(wl, aCoef[10, :3], 1.0)
    aSol.getAerTransmittance(wl, aCoef[11, :3], 1.0)
    aSol.getAerTransmittance(wl, aCoef[12, :3], 1.0)
    

def test_plotRawData():
    oSimu = aSol.SimuAeroStruct()
    oSimu.plotRawData()


#
# CHECK
#

def getOptDeep(tps, lInter):
    return np.array([oInter(tps) for oInter in lInter])


def interpolSelectRawData(pRaw):
    idxSel = np.where(pRaw[2] < 3)[0]
    #idxSel = np.where(np.logical_and(pRaw[2] > 1.5, pRaw[2] < 2))[0]
    rData = pRaw[:,idxSel]    
    lInter = [spi.interp1d(rData[0], optDeep) for optDeep in rData[3:8]]
    wlNM = np.array([380, 440, 500, 675, 870])
    pl.figure()
    print rData.shape
    pl.plot(wlNM, np.exp(-getOptDeep(rData[0,2000], lInter)))
    pl.plot(wlNM, np.exp(-getOptDeep(rData[0,2500], lInter)))
    pl.plot(wlNM, np.exp(-getOptDeep(rData[0,3000], lInter)))
    pl.ylim(ymin=0.8, ymax=1.0)
    pl.xlabel("nm")
    pl.grid()
    return rData, lInter

    
    
def viewRawData():
    oSimu = aSol.SimuAeroStruct()
    #idx = np.where(366*2<oSimu.RawData[0] and oSimu.RawData[0 > 366)[0]
    #rData = oSimu.RawData[:,idx]
    rData = oSimu.RawData
    pl.figure()
    pl.plot(rData[0], rData[4])
    pl.grid()
    pl.figure()
    pl.plot(rData[0], rData[4]*100)
    pl.plot(rData[0], rData[2],".")
    pl.grid()
    # selection entre 1 et 1.5 airmass
    idxSel = np.where(oSimu.RawData[2] < 1.5)[0]
    rData = oSimu.RawData[:,idxSel]
    pl.figure()
    pl.plot(rData[0], rData[4]*100)
    pl.plot(rData[0], rData[2],".")
    pl.grid()
    # fit spline
    #tck = spi.splrep(rData[0], rData[4]*100, k=1)
    #tck = spi.splrep(rData[0], rData[4]*100, s=8000)    
    xfin = np.linspace(rData[0][0], rData[0][-1], 500*len(rData[0]))
    deltaXfin_mn = 24*60*(rData[0][-1] - rData[0][0]) / len(xfin)
    print deltaXfin_mn
    #pl.figure()
    #pl.plot(xfin, spi.splev(xfin, tck))
    #pl.plot(rData[0], rData[4]*100,'.')
    #pl.plot(tck[0], spi.splev(tck[0], tck),'*')
    #pl.grid()
    # linear 
    interLin = spi.interp1d(rData[0], rData[4]*100)
    sigGrid = interLin(xfin)
    pl.figure()
    pl.title('MaunaLoa_1997-2009_final.dat: optical depth at 380nm')
    #pl.plot(xfin, sigGrid)
    pl.plot(rData[0], rData[4])
    sifSmooth = smooth(sigGrid, int(3*60/deltaXfin_mn))
    print sifSmooth.shape
    #pl.plot(xfin, sifSmooth)
    pl.xlabel('day')
    pl.grid()
    # 
    wlNM = np.array([380, 440, 500, 675, 870])
    pl.figure()
    print rData.shape
    pl.plot(wlNM, np.exp(-rData[3:8,200]))
    pl.plot(wlNM, np.exp(-rData[3:8,250]))
    pl.plot(wlNM, np.exp(-rData[3:8,300]))
    pl.ylim(ymin=0.8, ymax=1.0)
    pl.xlabel("nm")
    pl.grid()
    interpolSelectRawData(oSimu.RawData)
    

def multiplotRawData():
    oSimu = aSol.SimuAeroStruct()
    data, lInter = interpolSelectRawData(oSimu.RawData)
    beginDay = 4050
    beginDay = 12
    pasDay = 120
    nbPlot = pasDay*24
    wl = np.array([380, 440, 500, 675, 870])
    lday01 = np.linspace(0, pasDay, nbPlot, True)+beginDay
    for d1,idx in zip(lday01, range(nbPlot)):
        print idx,d1
        pl.figure()
        pl.title("aerosol trans. %s"% oSimu.getFileRawData())        
        [pl.plot(wl,np.exp(-getOptDeep(d1+_li*pasDay, lInter))) for _li in range(3)]
        strDay = ["day %.2f"%(d1+_li*pasDay) for _li in range(3)]
        pl.ylim(ymin=0.8, ymax=1.0)
        pl.grid()
        pl.xlabel("nm")
        pl.legend(strDay, loc=4)
        pl.savefig("/home/colley/temp/lsst/movie/aeroRaw%04d.png"%idx)
        pl.close()
        
        
def multiplotRawData2():
    """
    try to improve time execution , not validate
    """
    oSimu = aSol.SimuAeroStruct()
    data, lInter = interpolSelectRawData(oSimu.RawData)
    beginDay = 4050
    pasDay = 2
    nbPlot = pasDay*24
    wl = np.array([380, 440, 500, 675, 870])
    lday01 = np.linspace(0, pasDay, nbPlot, True)+beginDay
    aTps = np.outer(np.ones(3),lday01 )
    aTps  += np.array([0,pasDay, 2* pasDay]).reshape(3,1)
    aTps = aTps.ravel()
    aTrans = np.exp(-np.array([oInter(aTps) for oInter in lInter])).T
    for d1, idx in zip(lday01, range(nbPlot)):
        print idx,d1
        pl.figure()
        pl.title("aerosol transmission")        
        [pl.plot(wl,aTrans[idx+_li]) for _li in range(3)]
        strDay = ["day %.2f"%(d1+_li*pasDay) for _li in range(3)]
        pl.ylim(ymin=0.8, ymax=1.0)
        pl.grid()
        pl.xlabel("nm")
        pl.legend(strDay, loc=4)
        pl.savefig("/home/colley/temp/lsst/movie/aeroRaw%04d.png"%idx)
        pl.close()
 
#
# CLASS AeronetInterpol
#  

def test_AeronetInterpol_plotXXX():
    obj = aSol.AeronetInterpol()
    obj.plotHistoAM()
    obj.plotAM()
    obj.plotFMF()


def test_AeronetInterpol_interpolSplineInTime():
    obj = aSol.AeronetInterpol()
    xfin=np.linspace(300, 1100, 1000)
    obj.interpolSplineInTime(23, xfin)
    obj.interpolSplineInTime(50, xfin)
    obj.interpolSplineInTime(170, xfin)
    obj.interpolSplineInTime(220, xfin)
    obj.interpolSplineInTime(300, xfin)
    obj.interpolSplineInTime(340, xfin)
    obj.interpolSplineInTime(360, xfin)


def AeronetInterpol_studyFit():
    obj = aSol.AeronetInterpol()
    tps = 99.9
    aidxTps = np.where(np.logical_and(obj.RawData[0]> tps, obj.RawData[0]<tps+2))[0]
    print aidxTps    
    idx = aidxTps[0]
    print tps 
    print "time : ",obj.RawData[0,idx-2: idx+2]
    print obj.RawData[0,idx-1: idx+1]
    obj.plotTransmission(idx)
    xObs = obj.Wl
    yObs = np.exp(-obj.RawData[3:8, idx])
    
    #tck = spi.splrep(xObs, yObs, k=2, t=np.array([550])) 
    tck = spi.splrep(xObs, yObs, k=2, task=-1,t=np.array([580])) 
    tck = spi.splrep(xObs, yObs, k=2,t=np.array([600])) 
    xfin=np.linspace(300, 900, 1000)
    pl.plot(xfin, spi.splev(xfin, tck))
    print tck       


def test_AeronetInterpol_interpolONeil():
    obj = aSol.AeronetInterpol()
    obj.interpolONeil(1110)
    obj.interpolONeil(2110)
    obj.interpolONeil(3110)
    obj.interpolONeil(4110)


def test_AeronetInterpol_interpolONeilInTime():
    obj = aSol.AeronetInterpol()
    xfin = np.linspace(300, 1100, 1000)
    obj.interpolONeilInTime(23, xfin)
    obj.interpolONeilInTime(50, xfin)
    obj.interpolONeilInTime(170, xfin)
    obj.interpolONeilInTime(220, xfin)
    obj.interpolONeilInTime(300, xfin)
    obj.interpolONeilInTime(340, xfin)
    obj.interpolONeilInTime(360, xfin)


def AeronetInterpol_multiplotONeil():
    oSim = aSol.AeronetInterpol()
    beginDay = 230
    endDay = 130 +beginDay
    pasDay = 1./25
    nbPlot = int((endDay - beginDay)/pasDay)
    lday01 = np.linspace(beginDay, endDay, nbPlot, True)
    xfin = np.linspace(300, 1100, 1000)
    for d1,idx in zip(lday01, range(nbPlot)):
        print idx,d1
        pl.figure()
        pl.title("aerosol trans. %s fit O'Neil"% oSim.getFileRawData())        
        oSim.interpolONeilInTimeCache(d1, xfin)
        pl.xlabel("nm")
        pl.savefig("/home/colley/temp/lsst/movie/aeroONeil%04d.png"%idx)
        pl.close()
    

def check_correl():
    oSim = aSol.AeronetInterpol()
    print "aer380"
    print np.corrcoef(oSim.RawData[[2,3]])
    print "aer440"
    print np.corrcoef(oSim.RawData[[2,4]])
    print "aer500"
    print np.corrcoef(oSim.RawData[[2,5]])
    print "aer675"
    print np.corrcoef(oSim.RawData[[2,6]])
    print "aer870"
    print np.corrcoef(oSim.RawData[[2,7]])
    
    print "aer380 cor ?"
    copyData = oSim.RawData[[2,3]].copy()
    copyData[1] = copyData[1]*copyData[0]
    print np.corrcoef(copyData)
    
    print np.corrcoef(oSim.RawData[[1,3]])
    print np.corrcoef(oSim.RawData[[0,1]])
    
    print np.corrcoef(oSim.RawData[3:8,:])
     
#
# MAIN
#

S_PathFileSimu = "../testSImuAero.txt"

#test_getAerParametersFromFile(S_PathFileSimu)
#test_plotRawData()
#interpolRawData()
#viewRawData()
#multiplotRawData()

#
# CLASS AeronetInterpol

#test_AeronetInterpol_plotXXX()
#AeronetInterpol_studyFit()
#test_AeronetInterpol_interpolONeil()
#test_AeronetInterpol_interpolSplineInTime()
#test_AeronetInterpol_interpolONeilInTime()
#AeronetInterpol_multiplotONeil()
check_correl()


pl.show()
