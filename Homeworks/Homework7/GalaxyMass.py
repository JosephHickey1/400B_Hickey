#!/usr/bin/env python
# coding: utf-8

# Joseph Hickey  
# ASTR400B Homework 3  
# 2/5/2020  
# Create a function which returns the total mass of a component of a galaxy




import numpy as np
import astropy.units as u
from ReadFile import Read
import matplotlib.pyplot as plt





def ComponentsMass(file,particletype):
    #file takes the name of the file containing the mass data fo a given galaxy
    #particletype takes 1:Halo, 2:Disk, or 3:Bulge particles from the file
    
#Call the Read function passing it the file of interest
    time, N, data = Read(file)
#Get the array positions for all particles of a given type
    index = np.where(data['type'] == particletype)
#Get Masses of particles
    m = data['m'][index]
#Sum the new masses list to find total mass of given particle type
    Mass = sum(m) * (10**10 * u.Msun)
    Mass = np.around(Mass.to(10**12 * u.Msun),3)
    
    return Mass


