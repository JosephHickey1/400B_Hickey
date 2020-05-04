#!/usr/bin/env python
# coding: utf-8

import numpy as np
import astropy.units as u

def Read1(filename):
    
    #Put all the data into a sorted ordered array
    data = np.genfromtxt(filename,dtype=None,names=True)
    time = 1
    N = 1
    return time, N, data

