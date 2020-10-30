"""
Created on Oct 30 16:34 2020
    Engineering models for vertical axis wind turbine.
@author: Ming Huang
"""

import math
import numpy as np


def GaussianVAWT(Ct, x, y, z, zh, D, H, k_starD = None, k_starH = None, Ia = None):
    """    A kinematic wake model transposed from BPA model for HAWT, proposed by Abkar(2018).

    :param Ct: Thrust coefficient
    :type Ct: float/array
    :param x: non-dimentional downstream distance
    :type x: float/array
    :param y: non-dimentional horizontal coordinate
    :type y: float/array
    :param z: non-dimentional vertical coordinate
    :type z: float/array
    :param zh: Center height
    :type zh: float/array
    :param D: Rotor diameter
    :type D: float/array
    :param H: Rotor Height
    :type H: float/array
    :param k_starD: [description]
    :type k_starD: float/array
    :param k_starH: [description]
    :type k_starH: float/array
    :param Ia: [description]
    :type Ia: float/array
    
    """
    if k_starD is None:
        if Ia <= 0.134 and Ia >= 0.069:
            k_star = 0.3837*Ia + 0.003678;
            k_starD,k_starH = k_star,k_star
        
    Dp = np.sqrt((1+np.sqrt(1-Ct))/(2*np.sqrt(1-Ct)))
    x0 = 1/k_starD*(np.sqrt(Ct/8)-0.2*Dp/D)
    sigmaD = 0.2*Dp*D + k_starD*x
    sigmaH = 0.2*Dp*H + k_starH*x
    EXP_VALUE = -0.5*(np.square((z-zh)/sigmaH) + np.square(y/sigmaD))
    Uxc_over_U0 = np.sqrt(abs(1-(Ct/(2.*math.pi*(sigmaD)*(sigmaH)/(math.pi*H*D/4)))))
    DU_U0 = (1 - Uxc_over_U0)*np.exp(EXP_VALUE)
    return DU_U0,sigmaD,sigmaH
    
def TophatVAWT(Ct,D,H,kD,kH,x,y=None,z=None,zh=0):
    """    A kinematic wake model transposed from Jensen model for HAWT, proposed by Abkar(2018).
    :param Ct: Thrust coefficient
    :type Ct: float/array
    :param x: non-dimentional downstream distance
    :type x: float/array
    :param y: non-dimentional horizontal coordinate
    :type y: float/array
    :param z: non-dimentional vertical coordinate
    :type z: float/array
    :param zh: Center height
    :type zh: float/array
    :param D: Rotor diameter
    :type D: float/array
    :param H: Rotor Height
    :type H: float/array
    :param k_starD: [description]
    :type k_starD: float/array
    :param k_starH: [description]
    :type k_starH: float/array
    """
    a = (1- np.sqrt(1 - Ct))/2
    Hw = (1 + 2*kH*x/H)*H; Dw =  (1 + 2*kD*x/D)*D
    Htop = zh+ 0.5*Hw; Hbot = zh - 0.5*Hw
    Dl = -0.5*Dw; Dr = 0.5*Dw
    Ur_over_U0 = 1 - 2*a/(1 + 2*kH*x/H)/(1 + 2*kD*x/D)
    if y != None and z != None:
        Ur_over_U0[z >= Htop] = 1; Ur_over_U0[z <= Hbot] = 1
        Ur_over_U0[y <= Dl] = 1; Ur_over_U0[y >= Dr] = 1

    return Ur_over_U0, Hw, Dw

