from Process import Correlate
from Useful_Functions import *
from Constants import *
from Graph import *
from pylab import *
from math import *
import numpy as np
import scipy as sp
import random
import os,sys
from Sky_Background import Get_Sky_Background
import time
import cPickle as pickle
import scipy.integrate as integrate
import scipy.signal as ss

''' OBSERVE
Returns number of photons / SRE. 

Input Parameters:

   ELG parameters: mu, sigma, amplitude
   The bin the ELG is in: bin_no
   Telescope diameter: tel_diam (default 8m)
   Resolution: resolution (default 20000)

   Switches:

   Noise: on/off
   Sky background: on/off

Returns:

  SRE_wavelengths (list): Wavelength of each spectral resolution element (SRE)
	SRE_energies (list): The signal energy in each SRE
	SRE_photons (list): The number of photons collected in each SRE

'''
def ObserveELG(mu,sigma,amplitude,
                bin_wavelength_range=100.,SRE_width=0.5,
	            tel_diam=8.,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
                n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
                noise='on',background='on',background_subtraction='on',plot='on'):

   tel_area = pi * (tel_diam/2.0)**2.
   fibre_area = pi*(fibre_diam/2.)**2 # sq. arcseconds
   
   N_SRE = int(bin_wavelength_range / SRE_width)
   wavelength_start = myround(mu, int(bin_wavelength_range/2.)) - int(bin_wavelength_range/2.)

   def add_noise(signal_energy,SRE_wavelength,t_obs): # Observing time only used for read noise
 
      predicted_photons = getNphotons(signal_energy,SRE_wavelength)
      
      # Add shot noise
      try: # Fails if too many photons
         actual_photons = np.random.poisson(predicted_photons) # Number of photons actually hitting the sensor over whole time period 
      except:
         actual_photons = int(predicted_photons) 
      
      # Add read noise at each exposure. See http://bit.ly/Svgdcl
      stdev_electrons = sqrt((readout_error**2. + dark_current)*n_pix)
      # Add readout noise in quadrature at for multiple exposures
      exposures = int(t_obs/t_int)
      total_stdev_electrons = sqrt(exposures * stdev_electrons**2.) # Normal stdev add in quadrature
      
      actual_photons += ((random.gauss(0,total_stdev_electrons))) # abs (bit dodgy)?
#       if actual_photons < 0: # Noise can go negative
#          actual_photons = 0

      # Convert back to energy
      noise = (actual_photons-predicted_photons)/predicted_photons
      return getEphotons(int(actual_photons),SRE_wavelength) # Returns energy in ergs

   labels = (wavelength_start,SRE_width) # Information to pass to graph
   _SRE_wavelengths, _SRE_energies, _SRE_photons = [],[],[]

   ######### PROCESS SIGNAL IN EACH SRE AND ADD NOISE#########
   for SRE_no in range(0,N_SRE):
      
      # SRE_wavelength is half way in the bin
      SRE_wavelength = wavelength_start + (SRE_no + 0.5)*SRE_width
      # Multiply by SRE_width to get rid of /A
      if background == 'on':
         # Sky background 
         Background_Photons = Get_Sky_Background(SRE_wavelength) # photons/s/sqm/A/sqarcsec
         Background_Photons = Background_Photons * t_obs * tel_area * SRE_width * fibre_area # photons
         sky_background =  Background_Photons * getEphotons(1,SRE_wavelength) # erg
      elif background == 'off':
         sky_background = 0. # erg
      # Note multiplication by SRE_width. Higher resolution means fewer photons per bin
      ELR_signal_energy = gaussian(SRE_wavelength,amplitude,mu,sigma) * t_obs * SRE_width * tel_area # Signal energy in ergs
      # Incorporate throughput as only a fraction of the light makes it through to the detector
      raw_signal_energy = throughput * (ELR_signal_energy + sky_background) # Ideal signal before processing
      #SRE_energy = raw_signal_energy + noise * noise_suppressor # Combine to make noisy processed signal
      if noise == 'on':
         SRE_energy = add_noise(raw_signal_energy,SRE_wavelength,t_obs) # In ergs
      elif noise == 'off':
         SRE_energy = raw_signal_energy
      
      # Used for spectra with high background levels
      if background_subtraction == 'on':
         ref_background = throughput * sky_background
         #print SRE_energy,ref_background
         SRE_energy = SRE_energy - ref_background
         # if SRE_energy < 0: SRE_energy = 0; # Can go negative
      
      # Update SRE
      _SRE_wavelengths.append(SRE_wavelength) # In A
      _SRE_energies.append(SRE_energy)
      _SRE_photons.append(int(getNphotons(SRE_energy,SRE_wavelength)))
   
   # End for loop
   
   if plot == 'on': # Plot the spectra to have a look
      if random.random() < 1. : # Plot ~ 1 in 10 spectra
         SRE_nos = list(xrange(N_SRE)) # List of integers from [0,N_SRE)
         Plot_Bar_Chart(_SRE_wavelengths, _SRE_photons, labels)
    
   return _SRE_wavelengths,_SRE_energies, _SRE_photons


''' RV_HARPS
Returns a predicted velocity error (exoplanet hunting method) from an input spectrum and an SRE width

Input parameters:

	Spectra_Wavelengths (list)
	Spectra_Photons (list)
	SRE_width (float)

Returns:

	Quality factor Q
	Expected velocity error dV (cm/s)

'''

def Correlate(lam,Spectra_Phot_t, Spectra_Phot_t0, SRE_width=0.5):
      
   def quadratic_Method(X,k):
      y1 = X[k-1]
      y2 = X[k]
      y3 = X[k+1]
      d = (y3-y1)/(2*(2*y2 - y1 - y3))
      return k + d

   # Flipped convolution is a correlation
   corr = ss.fftconvolve(np.asarray(Spectra_Phot_t), np.asarray(Spectra_Phot_t0[::-1]))
   
   # Get peak of correlation to the nearest SRE
   peak = np.argmax(corr)
   
   # Get sub-pixel peak using quadratic interpolation
   px = quadratic_Method(corr, peak)

   # Get wavelength_shift from peak
   SRE_shift = px - (len(Spectra_Phot_t0) - 1)
   dlam = SRE_shift * SRE_width
   
   # dV / V = c dlam / lam
   dV = c * dlam / lam # m/s
   
   return dV * 100 # cm/s

def RV_HARPS(Spectra_Wavelengths,Spectra_Phot,SRE_width=0.5,det_noise=4.): # In photoelectrons
   
   Weights = []
   
   i = 0
   for l_i,A0_i in zip(Spectra_Wavelengths,Spectra_Phot):
   
      if i != (len(Spectra_Phot) - 1):
         dA0_dl_i =  float(Spectra_Phot[i+1]-Spectra_Phot[i])/SRE_width # Photoelectrons / A
      else: dA0_dl_i = 0 # Cant find derivative off end of the plot
      
      #print dA0_dl_i ; sleep(0.1)
      
      num = l_i**2. * dA0_dl_i**2. # photoelectrons**2
      denom = A0_i + det_noise**2. # photoelectrons**2
      
      W_i = num/denom
      Weights.append(W_i)
      
      i = i + 1
    
   Q = sqrt(np.sum(Weights)) / sqrt(np.sum(Spectra_Phot)) # Quality factor
    
   dV_RMS =  c * 1. / sqrt(np.sum(Weights))   # m/s
   
   return Q,dV_RMS * 100 # cm/s
  
def Unit_Test_ObserveELG(F):
	
	mu,sigma = 5000.,5.
	amplitude = get_amplitude(mu,sigma,F) # erg/s/m2/A
	print 'amp',amplitude

# Check integral if need be	
# 	import scipy.integrate as integrate
# 	integral = integrate.quad(lambda x: gaussian(x,amplitude,mu,sigma),mu-10*sigma,mu+10*sigma)[0]
# 	print integral
	
	_wavelengths,_energies,_photons =  ObserveELG(mu,sigma,amplitude,
					bin_wavelength_range=500.,SRE_width=0.5,
					tel_diam=8.,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
					n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
					noise='on',background='on',background_subtraction='on',plot='on')
					
	for i in range(0,len(_wavelengths)):
		print _wavelengths[i], _energies[i], _photons[i]

def RV_errors():
	_logF = np.linspace(-18.,-10.,20)
	_F = np.power(10,_logF)
	tel_diam = 8.
	
	for F in _F:
	
		mu,sigma = 5000.,5.
		amplitude = get_amplitude(mu,sigma,F) # erg/s/m2/A
	
		# Create one noisless spectrum for exoplanet HARPS technique
		_wavelengthsNnoise1,_energiesNnoise1,_photonsNnoise1 =  ObserveELG(mu,sigma,amplitude,
						bin_wavelength_range=100.,SRE_width=0.5,
						tel_diam=tel_diam,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
						n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
						noise='off',background='on',background_subtraction='on',plot='off')
			
		# Create two noisy spectra for cross-correlation	
		_dV_Corr = []
		for i in range(0,25):
			#print F,i
			_wavelengths1,_energies1,_photons1 =  ObserveELG(mu,sigma,amplitude,
							bin_wavelength_range=100.,SRE_width=0.5,
							tel_diam=tel_diam,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
							n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
							noise='on',background='on',background_subtraction='on',plot='off')
							
			_wavelengths2,_energies2,_photons2 =  ObserveELG(mu,sigma,amplitude,
							bin_wavelength_range=100.,SRE_width=0.5,
							tel_diam=tel_diam,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
							n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
							noise='on',background='on',background_subtraction='on',plot='off')
	
			_dV_Corr.append([abs(Correlate(mu,_photons2, _photons1, SRE_width=0.5))])
						
		dV_HARPS = RV_HARPS(_wavelengthsNnoise1,_photonsNnoise1,SRE_width=0.5,det_noise=4.)[1]
		dV_CORR = np.mean(_dV_Corr)
		#print F, dV_HARPS, dV_CORR
		print SRE_width, dV_HARPS


def megagal_errors():
	F = 10**(-15.) # 10**-17 or 10**-15
	tel_diam = 8.
	
	_logN_bin = np.linspace(0,5,100.)
	N_bin = np.power(10,_logN_bin)
	
	for N in N_bin:
		mu,sigma = 5000.,5.
		amplitude = get_amplitude(mu,sigma,F) # erg/s/m2/A
	
		# Create one noisless spectrum for exoplanet HARPS technique
		_wavelengthsNnoise1,_energiesNnoise1,_photonsNnoise1 =  ObserveELG(mu,sigma,amplitude,
						bin_wavelength_range=100.,SRE_width=0.5,
						tel_diam=tel_diam,t_obs = 1000.*60**2.,t_int = 3.*60**2.,
						n_pix = 4,throughput = 0.1,dark_current = 0,readout_error = 4,fibre_diam=3.,
						noise='off',background='on',background_subtraction='on',plot='off')
						
		# Make megagalaxy
		_photonsNnoise1 = [elem * N for elem in _photonsNnoise1]
						
		dV_HARPS = RV_HARPS(_wavelengthsNnoise1,_photonsNnoise1,SRE_width=0.5,det_noise=4.)[1]
		print N, dV_HARPS

# Monte carlo is independent to the rest
def Monte_Carlo():

	_zLF = pickle.load(open('zLF_OII_Urquhart_24dF_BRIGHT.p', "rb" ) ) 
	
	flux_cutoff = 10**(-15.25)
	L_cutoff = 10**42.25
	_zLF = [elem for elem in _zLF if elem[2] > flux_cutoff] # Apply flux cutoff
	_zLF = [elem for elem in _zLF if elem[1] > L_cutoff]
	
	# Telescope parameters
	throughput,D,t_int = 0.25,8.,1000
	
	_zLF_dV = []
	for i,gal in enumerate(_zLF):
		z,L,F = gal[0],gal[1],gal[2]
		
		
		dV_5 = 4500 * (throughput/0.1)**(-0.5) * (F/10**(-15.))**(-0.5) \
		 * (D/8.)**(-1.) * (t_int/1000)**(-0.5) # cm/s
		
		# Get sigma
		L_star = 10.**(0.54*z + 41.56)
		sigma = 12.5 * (L/L_star)**(-0.5)
		# Scale by sigma
		A,B,C = 1.23,-12.62,1084.5
		x = sigma/5.
		dV = dV_5 * x*(A*x**2. + B*x + C) / 1073.11
		
		_zLF_dV.append([z,L,F,dV])
		print z,L,F,sigma,dV
	
	pickle.dump(_zLF_dV, open('zLF_OII_Urquhart_24dF_10-16_FINAL.p', "wb" ))
	#pickle.load(open('zLF_OII_Urquhart_24dF_10-16_FINAL.p', "rb" ) )
	
	# Monte Carlo data
	zstart,zfin = 0.5,1.5
	zstep = 0.25
	Nbins = int((zfin-zstart)/zstep)
	
	for i in range(0,Nbins):
		zmin,zmax = zstart + i * zstep, zstart + (i+1) * zstep
		
		bin_dV = [elem[3] for elem in _zLF_dV if zmin < elem[0] < zmax] # Bin up data	
		
		frac = [elem**(-2.) for elem in bin_dV]
		sigma_ov = ( math.sqrt (np.sum(frac)))**(-1.)
		
		gal_in_bin = len(bin_dV)
		print zmin,zmax,sigma_ov,gal_in_bin

###
#Unit_Test_ObserveELG(10**-15.)

# Two routines - RV errors of megagal errors
#megagal_errors()
Monte_Carlo()
#RV_errors()
