#!/usr/bin/env python
# coding: utf-8

import numpy as np
import astropy.units as u

def Read(filename):
    
#Open file in same directory
    file = open(filename,'r')
#Readout first line and split it into component word and value. Assign units to value
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
#Readout second line and split it and assign units.
    line2 = file.readline()
    label, value = line2.split()
    N = float(value)
#Close file
    file.close()
#Put all the other data into a sorted ordered array
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    return time, N, data

