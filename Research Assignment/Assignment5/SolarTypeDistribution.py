#!/usr/bin/env python
# coding: utf-8

# ## Evolution of the spatial distribution of stars in bound orbits at the solar radius in M31

# With this code I will find the spatial distribution of stars within M31 matching the kinematic properties of the sun as the MW-M31 merger commences. The final environments of the stars will be appropriately histogrammed to provide a probability distribution.

# The redistribution of material during the merger is interesting because it can have an effect on the distribution of metal rich dust post-merger. This metal rich dust is the environment in which stars form and, if it migrates, there can be interesting features in the merger product such as a flocculant distrbution of star formation.
# The tidal forcing which causes the redistribution will change not only the position distribution but a velocity distribution as well. This should in theory lead to a more dispersion supported merger product rather than an entirely rotation supported galaxy.

# In[2]:


# Import useful libraries and developed modules.
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from ReadFile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentsMass


# In[164]:


# This class provides a mask for the data read in by ReadFile.py which can be used at every time step
# to track the movements of particles starting close to a given radius which are bound to their host
class Distribution:
    
    # Read in files and assign global variables for use elsewhere
    def __init__(self,galaxy,rd=5,rb=1,rh=52):
        # galaxy is the shorthand name for the target galaxy: MW, M31, M33
        # rd, rb, rh are the scale radii of the disk, bulge, and halo respectively
        # for the target galaxy. Defaults are set for M31 as that is the assignment
        
        # Create unitless constant for calculation of potentials later
        self.G = const.G.to(u.kpc**3 / u.Msun / u.Gyr**2).value
        
        # Convert input name to a filename to read from
        galaxyname = str(galaxy) + "_000.txt"
        # Read in data making the array global
        T,N,self.Gal = Read(galaxyname)
        
        # Call the CenterOfMass class and use it to find zero points in the reference frame later
        GalCOM = CenterOfMass(galaxyname,2)
        # Gets the position of the COM using disk particles
        self.Galpos = GalCOM.COM_P(0.1,4).value
        # Gets the velocity of the COM using it's position
        self.Galvel = GalCOM.COM_V(self.Galpos[0],self.Galpos[1],self.Galpos[2]).value
        
        # Scale radii and the masses of each particle type are read in or set manually
        self.rdisk = rd
        self.Mdisk = ComponentsMass(galaxyname,2).value * 1e12
        self.rbulge = rb
        self.Mbulge = ComponentsMass(galaxyname,3).value * 1e12
        self.rhalo = rh
        self.Mhalo = ComponentsMass(galaxyname,1).value * 1e12
        
    
    
    # This function calculates the potential of a particle due to the parts of  a
    # galaxy that can be modeled as a Hernquist sphere (i.e. the bulge and halo)
    def HernquistPotential(self,M,r_a,x,y,z):
        # M is the mass of either the halo or bulge of the galaxy
        # r_a is the scale radius of the bulge or halo
        # x,y,z are the current position of the particle in the galaxy's frame
        
        # r is the magnitude of the separation
        r = np.sqrt(x**2 + y**2 + z**2)
        
        # Calculate potential      >>>>       Hernquist 1990, ApJ...356..359H
        potential = self.G * M / (r + r_a)
        
        # Returns the scalar potential
        return potential
        
    
    # This function calculates the potential of a particle due to the extended disk of a galaxy
    # The disk cannot be modeled as a Herquist sphere but rather as a torus
    def MiyamotoNagaiPotential(self,M,rd,x,y,z):
        # M is the mass of the disk
        # rd is the scale radius of the disk
        # x,y,z are the current position of the particle in the galaxy's frame
        
        # R is the distance the COM and the particle on an x-y projection
        R = np.sqrt(x**2 + y**2)
        
        # B is a corrected z position to make up for the lack of spherical symmetry of the disk
        B = rd + np.sqrt(z**2 + (rd / 5.)**2)
        
        #  Calculate the potential      >>>>       Miyamoto & Nagai 1975 PASJ...27..533M
        potential = self.G * M / ((R**2 + B**2)**0.5)
        
        # Returns the scalar potentiaL
        return potential
        
    
    # This function sums the acceleration due to each component of M31
    def TotalPotential(self,pos):
        # pos is a vector containing the current particle position
        
        # Call earlier functions with global variables and current position
        bulgepot = self.HernquistPotential(self.Mbulge,self.rbulge,pos[0],pos[1],pos[2])
        diskpot = self.MiyamotoNagaiPotential(self.Mdisk,self.rdisk,pos[0],pos[1],pos[2])
        halopot = self.HernquistPotential(self.Mhalo,self.rhalo,pos[0],pos[1],pos[2])
        
        # Returns a total acceleration vector as a numpy array
        return bulgepot + diskpot + halopot
        
    
    # Translate particle positions/velocities to COM positions/velocities
    # Probably could have been done in __init__ function
    def FrameCorrection(self):
        # Takes no inputs
        
        # Subtract a value from each component of the position/velocity data
        position = np.array([self.Gal['x']-self.Galpos[0], self.Gal['y']-self.Galpos[1], self.Gal['z']-self.Galpos[2]])
        velocity = np.array([self.Gal['vx']-self.Galvel[0], self.Gal['vy']-self.Galvel[1], self.Gal['vz']-self.Galvel[2]])
        
        # Returns a pair of same sized 2-D arrays
        return position,velocity
        
    # Determines the escape velocity of a particle at a given radius
    # If velocity is below this, it is bound. 
    # Could be refined to take only near circular orbits
    def VBound(self,pos):
        # pos is a 2-D array of particle relative position information
        
        # Vesc = sqrt(2*potential) and the potential is position dependent
        
        # Initialize an empty array to store data
        Vesc = np.zeros(len(self.Gal['x']))
        
        # Loop over each particle to return their escape speeds
        for i in range(len(self.Gal['x'])):
            Vesc[i] = np.sqrt(2*self.TotalPotential(pos[:,i]))
        
        # Reurns an array with a length equal to the input array
        return Vesc
    
    
    # Creates a mask for the initial data that was read in which identifies the particles 
    # with Sun-like orbital characteristics which should be usable in future snapshots
    def Mask(self,radius):
        # Takes in a radius to show stars within 10% of
        
        # Gets the positions and velocities of all particles in the COM reference frame
        pos, vel = self.FrameCorrection()
        
        # Uses coordinates to find magnitude of distance and velocity
        dist = np.sqrt(sum(pos**2))
        speed = np.sqrt(sum(vel**2))
        
        # Sets 10% bounds on the radius to narrow the sample
        bounds = [radius*1.1,radius*0.9]
        
        # Creates an index where each entry is a disk particle within the range of radii and is bound to the system
        index = np.where((self.Gal['type'] == 2) & (dist <= bounds[0]) & (dist >= bounds[1]) & (np.absolute(speed) < self.VBound(pos)))
        
        # Returns an index array
        return index
        
    
    # Creates a 2-D histogram of the data
    def Histogram(self,coord1,coord2,radius=8):
        # coord1 & coord2 are the axes you want to project information from (x,y,z,vx,vy,vz)
        # radius is the distance from the COM of the galaxy at which you want to find particles
        
        # Make sure inputs are in a string format
        coord1 = str(coord1)
        coord2 = str(coord2)
        
        # Apply mask to global data
        stars = self.Gal[Mask(radius)]
        
        # Create a plot of the two data axes. 
        plt.hist2d(stars[coord1],stars[coord2],bins=(len(stars[coord1])/40))
        
