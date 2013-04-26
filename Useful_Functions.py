from Background_OH import *
from Constants import *
from math import *
import numpy as np
import scipy.integrate as integrate
    
# Cosmological dynamics (zdot functions etc.)
def get_H(H0,O_l,O_m,O_r,O_k, z): # Get H(z)
  E = lambda z: O_l + O_m * (1+z)**3 + O_r * (1+z)**4 + O_k * (1+z)**2
	return H0 * np.sqrt(E(z)) 

def get_zdot(H0,O_l,O_m,O_r,O_k, z):  # Get zdot
    zdot_s = H0*(1+z) - get_H(H0,O_l,O_m,O_r,O_k, z) # /s
    return zdot_s * (365*24*60**2) # in /yr 

def gaussian(x,A,mu,sigma):
   return A * exp(-0.5*((x-mu)/sigma)**2)
   
# Round to specified number (base). eg. round(470,50) = 450
def myround(x, base):
   return int(base * round(float(x)/base))
   
def getNphotons(E_tot,mu): # Takes energy in ergs, wavelength in A
   E = h*c/(mu*10**(-10.)) # Energy of one photon in J
   return (E_tot*10**(-7.))/E # N photons
   
def getEphotons(N,mu): # Takes wavelength in A; returns E in ergs
   E = h*c/(mu*10**(-10.)) # Energy of one photon in J
   return N * E * 10**(7.) # Energy of all photons in ergs
   
def ListToDict(Keys,Values): # Maps list of values, list of keys -> dict
   assert (len(Values) == len(Keys)),"Length of List (" + str(len(Values)) + \
   ") should equal length of Keys (" + str(len(Keys)) + ")"
   
   dictionary = dict(zip(Keys, Values))
   return dictionary

def ConvertToDict(Keys,List,SortKey): # Maps list of lists -> list of dicts
   '''Keys: The keys to assign to each element of the list
   List: The list of data to be converted to a dictionary
   SortKey: The key used to order the list of dictionaries (increasing)'''
   for i,elem in enumerate(List):
      dict = ListToDict(Keys,elem)
      List[i] = dict
   return sorted(List, key=lambda k: k[SortKey]) 
   
def get_wavelength(line, z):
    return (1+z)*{
        'H_alpha': 6562.8, # Rest wavelengths / A
        'OII': 3727.1,
        'OIII': 4958.9,
    }[line]
    
def Get_obs_z_dot(mu,wavelength_shift,z): 
   obs_z_dot = (wavelength_shift / mu) / wait_time * (1+z) # Why is this 1+z factor necessary!!!?
   return obs_z_dot  
   
def get_bin(mu):    
   base = int(bin_wavelength_range/2.) # what to round to
   bin_no = int((myround(mu,base)-lower_wavelength_range)/base) - 1
   return bin_no
   
# Takes mu, sigma in A, flux in erg/s/cm2. 
# Returns amplitude in ergs/s/A/m2
def get_amplitude(mu, sigma, flux): # Brute force amplitude of gaussian (in ergs/s/A)

   amplitude_tolerance = 10.**-5

   A1,A2 = 0,10**(-3) # Initial amplitudes to brute force
   integral = integrate.quad(lambda x: gaussian(x,1,mu,sigma),mu-10*sigma,mu+10*sigma)[0]
   while 1:
	  A = A1 + (A2 - A1)/2.
	  # Integral of gaussian [-inf,inf] = A/2 (1+erf(x-mu/sigma))
	  try_flux = A * integral
	  error = try_flux - flux
  
	  if(abs(error/flux) < amplitude_tolerance): # How tight does this need to be?
		 return A * 10**4. # Return amplitude in ergs/s/A/m2 (10^4 takes cm2 to m2)
	  if(error > 0):
		 A2 -= (A2-A1)/2.
	  if(error < 0):
		 A1 += (A2-A1)/2. # Existing python optimisation routines aren't tight enough
