import numpy as np
import LF_Data
import sys
from cosmolopy.luminosityfunction import schechterL
from cosmolopy.distance import luminosity_distance,comoving_volume
import scipy
import random
import math
import time
import cPickle as pickle
from pylab import Rectangle
import matplotlib.pyplot as plt
import matplotlib

def DensityPlot(x,y,z):
  # Create a density plot
	import numpy as np
	import matplotlib.pyplot as plt
	import scipy.interpolate
	
	x,y,z = np.array(x),np.array(y),np.array(z)
	
	for i in range(0,len(x)):
	   print x[i],y[i],z[i]
	
	# Set up a regular grid of interpolation points
	xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
	#xi, yi = np.meshgrid(xi, yi)
	
	# Interpolate
	#rbf = scipy.interpolate.Rbf(x, y, z, function = 'linear')
	densityfit = scipy.interpolate.interp2d(x, y, z, kind = 'linear')
	zi = densityfit(xi, yi)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	# Plot labels
	for i in range(0,len(z)):
	   ax.annotate(
		  str(z[i]), 
		  xy = (x[i],y[i])
		  , xytext = (0, 0),
		  textcoords = 'offset points', ha = 'center', va = 'center', color='white',size=13, fontname='Times New Roman')
	
	
	# set your ticks manually
	#ax.xaxis.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5])
	ax.xaxis.set_ticks([0.5,0.75,1.0,1.25,1.50])
	
	plt.xlabel('redshift',fontname='Times New Roman',size=16)
	plt.ylabel('log L [erg/s]', fontname='Times New Roman',size=16)
	
	plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
			   extent=[0.5, 1.5, 38, 45.], aspect='auto',cmap = 'jet')
	#plt.scatter(x, y, c=z)
	
	matplotlib.rcParams.update({'font.size': 13, 'font': 'Times New Roman'})
	
	ax.grid(color='#B0B0B0',linestyle='-',linewidth=1.)
	
	_F = [10**(-18.),10**(-17.),10**(-16.),10**(-15.),10**(-14.)]
	for F in _F:
		_z = np.linspace(0.001,1.5,100)
		_L = []
		for z in _z:
			DL = luminosity_distance(z, **cosmo) # In MPc
			DL = DL * (10**6. * pc) * 10**2. # in cm
			L = math.log10(4 * math.pi * DL**2. * F)
			_L.append(L)
		plt.plot(_z,_L,'w--')
          
	#plt.xlim(0.0,0.5)
	plt.xlim(0.5,1.5)
	plt.ylim(38,44)
	plt.savefig('LvsZ.png')


pc = 3.08567758*10**(16.) #m

# Telescope FOV
line = 'O_II'
FOV = 3. # deg
Tel_Area = math.pi * (FOV/2.)**2.

### Gridding parameters
no_zsteps, no_Lsteps = 4,6
zstart, zstep = 0.5, 0.25
Lstart, Lmult = 10**38., 10.
###

O_l,O_m,O_r,O_k,h = 0.7,0.3,0.0,0.0,0.7 # Cosmological parameters
cosmo = {'omega_k_0': O_k, 'omega_M_0' : O_m, 'omega_lambda_0' : O_l, 'h' : h}

Dist = []
for i in range(0,no_zsteps):
   for j in range(0,no_Lsteps):
      # Define params
      zmin = zstart + i*zstep
      zmax = zstart + (i+1)*zstep
      z = zmin + 0.5*(zmax-zmin) # Effective z to draw LF from
      
      Lmin = Lstart * Lmult**(j) 
      Lmax = Lstart * Lmult**(j+1) 

      # Compute comoving volume in Mpc-3
      dV = ( comoving_volume(zmax,**cosmo) - comoving_volume(zmin,**cosmo) ) * (Tel_Area / 41253.) # Mpc3

      # Get Phi in Mpc-3
      alpha,logLStar,logPhiStar = LF_Data.retrieve_LF(z,line)
      LStar,PhiStar = 10.**(logLStar),10.**(logPhiStar)
      
      n = scipy.integrate.quad(lambda L: schechterL(L,PhiStar, alpha, LStar),Lmin, Lmax)[0]  

      # Number = number density * volume
      N = n * dV

      # Compute flux
      DL = luminosity_distance(z, **cosmo) # In MPc
      DL = DL * (10**6. * pc) * 10**2. # in cm
      F = Lmin / (4 * math.pi * DL**2.) # ergs/s/cm2

      Dist.append([zmin,zmax,math.log10(Lmin),math.log10(Lmax), int(N), F])

# Plot L vs. z
x = [elem[0] + 0.5 * (elem[1] - elem[0]) for elem in Dist] # av z
y = [elem[2] + 0.5 * (elem[3] - elem[2]) for elem in Dist] # av log10
z = [elem[4] for elem in Dist] 

#DensityPlot(x,y,z)

# Make zLF distributions
total = np.sum(z)
_zLF = []
for box in Dist:
    
	def getF(z,L):
		DL = luminosity_distance(z, **cosmo) # In MPc
		DL = DL * (10**6. * pc) * 10**2. # in cm
		return L / (4 * math.pi * DL**2.) # ergs/s/cm2

	N = box[4]
	zmin,zmax = box[0], box[1]
	logLmin, LogLmax = box[2],box[3]
	
	for gal in range(0,N):
		z = random.uniform(zmin,zmax)
		
		log10L = random.uniform(logLmin, LogLmax)
		L = 10.**log10L 
		
		F = getF(z,L)
		
		_zLF.append([z,L,F])
		print z,L,F
		print len(_zLF)/float(total)*100,'% complete'


pickle.dump(_zLF, open('zLF_OII_Urquhart.p', "wb" ))






      
