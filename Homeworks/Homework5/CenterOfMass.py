#!/usr/bin/env python
# coding: utf-8

# Joseph Hickey  
# 2/12/2020  
# Astr400B Homework 4  
# Create a class which calculates the center of mass position and velocity of a given galaxy  
# 

# In[1]:


import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read


# In[2]:


#Define the class
class CenterOfMass:
    
#Initialize the class
    def __init__(self,filename,particletype):
        
        #Self refers to an aspect of the class that becomes a local variable
        self.time, self.total, self.data = Read(filename)
        
        #Index of particle type given
        self.index = np.where(self.data['type'] == particletype)
        
        #Get info on all particles of given type and put it in new arrays
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        
    
#Function returns center of mass in three axes
    def COMdefine(self,a,b,c,m):
        #a, b, c, and m are arrays of positions or velocities and mass
        Acom = np.sum(a * m)/np.sum(m)
        
        Bcom = np.sum(b * m)/np.sum(m)
        
        Ccom = np.sum(c * m)/np.sum(m)
        
        return Acom, Bcom, Ccom
        
    
#Returns position and velocity of COM
    def COM_P(self,delta):
        #Takes particletype and a tolerance of change as it iterates
        
        #First calculation from COMdefine
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)
        
        #Get initial values for the iteration by setting COMdefine zero to be local zero
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
        #Maximum distance to a particle from the guessed COM/2 to recalculate
        RMAX = max(RNEW)/2.
        Change = 1000.0
        
        #Iterate until the change is within the tolerance
        while(Change > delta):
            #Get all particles within RMAX
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
            
            #Get the new COM with COMdefine
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)
            
            #Check if the difference is within tolerance
            Change = np.abs(RCOM - RCOM2)
            RMAX /= 2.
            
            #Reset the separations and rmax
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
            
            #Set initial values to new values for next iteration
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2
            
            VEC = [XCOM,YCOM,ZCOM]
            
        
        #Set units and return the components
        VEC *= u.kpc
        VEC = np.around(VEC,2)
        
        return VEC
    
    
    def COM_V(self,COMX,COMY,COMZ):
        #Use the refined COM position to set your origin point and calculate the COM velocity
        RVMAX = 15.0 * u.kpc
        
        #Get locations of all particles w.r.t. COM
        xV = self.x - COMX.value
        yV = self.y - COMY.value
        zV = self.z - COMZ.value
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        #Get index for all particles within RVMAX
        indexV = np.where(RV < RVMAX.value)
        
        #Get velocity components for those particles
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
        
        #Compute COM velocity using COMdefine
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)
        
        VVEC = [VXCOM, VYCOM, VZCOM]
        VVEC = VVEC * u.km / u.s
        VVEC = np.around(VVEC,2)
        
        return VVEC
    
    
