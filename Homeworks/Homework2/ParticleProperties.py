#!/usr/bin/env python
# coding: utf-8

# Joseph Hickey  
# ASTR400B Homework 2  
# 1/29/2020  
# Write code which returns the properties of a requested particle  

# In[26]:


import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(file,particletype,particleno):
#Call the Read function passing it the file of interest
    time, N, data = Read(file)
#Get the array positions for all particles of a given type
    index = np.where(data['type'] == particletype)
#Get information on specified particle and convert it to floats with units
    x,y,z = data['x'][index], data['y'][index], data['z'][index]
    vx,vy,vz = data['vx'][index], data['vy'][index], data['vz'][index]
    m = data['m'][index]
    xpos,ypos,zpos = float(x[particleno])*u.kpc, float(y[particleno])*u.kpc, float(z[particleno])*u.kpc
    xvel,yvel,zvel = float(vx[particleno])*u.km/u.s, float(vy[particleno])*u.km/u.s, float(vz[particleno])*u.km/u.s
#Calculate magnitude of position and velocity of particle rounded to 3 places
    mass = m[particleno]*(10**10)*u.Msun
    distance = np.around(((xpos**2 + ypos**2 + zpos**2)**0.5),3)
    velocity = np.around(((xvel**2 + yvel**2 + zvel**2)**0.5),3)
    
    return distance, velocity, mass
    
    


# In[29]:


#Call previously defined function
d,v,m = ParticleInfo('MW_000.txt',3,99)
print(d,v,m)
#Convert distance in kiloparsecs to lightyears
d.to(u.lyr)

