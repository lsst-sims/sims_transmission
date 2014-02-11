# -*- coding: utf-8 -*-
'''
Created on 7 févr. 2014

@author: colley
'''
import lsst.sims.atmosphere.transmission.auxtelSequence as ats
import numpy as np

seed = 10

np.random.seed(seed)


def basicSimu():
    """
    """
    lPos = [{'ID':1,'MJD': 50.11,'ZANG':0.1, 'TGRAY':1.0} , 
            {'ID':2,'MJD': 62,'ZANG':0.1, 'TGRAY':1.0},
             {'ID':2,'MJD': 75,'ZANG':0.1, 'TGRAY':1.0}]
    lPos = [{'ID':1,'MJD': 56647.11,'ZANG':0.0, 'TGRAY':1.0},
            {'ID':2,'MJD': 56647.11,'ZANG': np.deg2rad(10), 'TGRAY':1.0},
            {'ID':3,'MJD': 56647.11,'ZANG': np.deg2rad(20), 'TGRAY':1.0},
            {'ID':4,'MJD': 56647.11,'ZANG': np.deg2rad(40), 'TGRAY':1.0},
            {'ID':5,'MJD': 56647.11,'ZANG': np.deg2rad(60), 'TGRAY':1.0}]
    lPos = [{'ID':1,'MJD': 56647.11,'ZANG':0.0, 'TGRAY':1.0},
            {'ID':2,'MJD': 56647.11,'ZANG': 10, 'TGRAY':1.0},
            {'ID':3,'MJD': 56647.11,'ZANG': 20, 'TGRAY':1.0},
            {'ID':4,'MJD': 56647.11,'ZANG': 40, 'TGRAY':1.0},
            {'ID':5,'MJD': 56647.11,'ZANG': 60, 'TGRAY':1.0}]
    lPos = [{'ID':1, 'MJD': 56727.11, 'ZANG':0.0, 'TGRAY':1.0},
            {'ID':1, 'MJD': 56747.51, 'ZANG':0.0, 'TGRAY':1.0},
            {'ID':1, 'MJD': 56748.51, 'ZANG':0.0, 'TGRAY':1.0}]
    # create object simu
    oSim = ats.AuxTelSequence(lPos)
    oSim.generateParameters(seed)
    oSim.getAtmTrans('tmp', save=True)




#
# MAIN
#

basicSimu()