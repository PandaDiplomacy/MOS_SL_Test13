# Compare with the following paper:
# Observations of the OH airflow emission, 1993, Yamashita
# Astronomical society of the pacific (value of 590 photons/s/m/arcsec2/micron at 1.65microns) 

# Flux data taken from: 
# UBVRI Data taken from E-ELT background analysis: http://www.eso.org/sci/facilities/eelt/science/drm/tech_data/background/
# JH Data Taken from Simons, Tokunaga: http://irtfweb.ifa.hawaii.edu/IRrefdata/paper_I.pdf

from scipy.interpolate import UnivariateSpline
import numpy as np
import math
import matplotlib.pyplot as plt

band_colours = ('DarkMagenta','DarkBlue','yellow','DarkRed','crimson','black','green','orange','pink')

# Reference Figure from Yamashita
# H = 18.1 at 1.6 microns => F = 7.33E-18 (factor of 4 lower than E-ELT data)

# Data from (the over-optimistic) OH suppression paper http://arxiv.org/pdf/0801.3870v1.pdf
# U = [3597, 625, 1.06*10**(-17.), 1.06*10**(-17.)]
# B = [4377, 890, 7.35*10**(-18.), 7.35*10**(-18.)]
# V = [5488, 830, 7.83*10**(-18.), 7.83*10**(-18.)]
# R = [6515, 1443, 8.62*10**(-18.), 3.76*10**(-18.)]
# I = [7981, 1499, 1.07*10**(-17.), 2.04*10**(-18.)] 
# J = [12483, 1600, 1.23*10**(-16.), 4.45*10**(-19.)]  
# H = [16313, 2900, 7.73*10**(-16.), 1.82*10**(-19.)] 

# Data taken from the E-ELT background data (Includes OH suppression?)
U = [3597, 625, 1.06*10**(-17.), 190.]
B = [4377, 890, 7.35*10**(-18.), 150.]
V = [5488, 830, 7.83*10**(-18.), 210.]
R = [6515, 1443, 3.76*10**(-18.), 340.]
I = [7981, 1499, 2.04*10**(-18.), 500.] 
J = [12483, 1600, 1.94*10**(-17.), 1200.] 
H = [16313, 2900, 3.16*10**(-17.), 2300.] 
K = [22000, 2500, 2.13*10**(-17.), 2300.]

Background_Data = [U]+[B]+[V]+[R]+[I]+[J]+[H]+[K]

# Manipulate the data into m_wavelengths, sky flux, and sky flux (OH suppressed)
m_wavelengths = [elem[0] for elem in Background_Data]
FWHM = [elem[1] for elem in Background_Data]
F_bgd = [elem[2] for elem in Background_Data]
N_bgd = [elem[3] for elem in Background_Data] # Photons/s/micron/m2/arcsec2

# Convert to log Y axis for a nicer plot
F_bgd = [math.log10(flux)for flux in F_bgd]

# Fit the data using a splined interpolation (no smoothing)
sky_fit = UnivariateSpline(m_wavelengths,F_bgd,s=0)
sky_fit_phot = UnivariateSpline(m_wavelengths,N_bgd,s=0)

# Generate 100 points of fitted data for no suppression and OH suppression fits
xs = np.linspace(min(m_wavelengths),max(m_wavelengths),1000)
ys = sky_fit(xs)

def Get_Sky_Background(wavelength): # Returns background in ergs/s/cm2/A/arcsec2

   # Sanity check
   if wavelength > 16000:
      #print "NO BACKGROUND DATA"
      return 0
      
   log_flux = sky_fit(wavelength) 
   
   bkgd_phot = sky_fit_phot(wavelength) # Photons/s/micron/m2/arcsec2
   return bkgd_phot * 10**-4. # Photons/s/A/m2/arcsec2
   #return 10**(log_flux) * 10**4. # ergs/s/m2/A/arcsec2 (10**4. converts cm2 to m2)

def Plot_Background(m_wavelengths,F_bgd,xs,ys):

   plt.scatter(m_wavelengths,F_bgd)
   plt.plot(xs,ys,'-',label='Integrated sky brightness')

   plt.title("Sky background as a function of wavelength")
   plt.xlabel("Wavelength / A"), plt.ylabel("log flux (ergs/s/cm2/A/arcsec2)") 
   x_min,x_max = (m_wavelengths[0] - FWHM[0]), (m_wavelengths[-1] + FWHM[-1])
   plt.xlim(x_min,x_max)
   plt.annotate

   # Shade in filters (optional)
   
   previous = 0
   for i in range(0,len(band_colours)-1):
      if m_wavelengths[i]-FWHM[i] < previous:
         min,max = previous, m_wavelengths[i]+FWHM[i]
      else: 
         min,max = m_wavelengths[i]-FWHM[i],m_wavelengths[i]+FWHM[i]
      plt.axvspan(min,max,facecolor=band_colours[i],alpha=0.35)
      previous = m_wavelengths[i]+FWHM[i]
   
   
   labels=['U','B','V','R','I','J','H','K']
   
   for label, x, y in zip(labels, m_wavelengths,F_bgd):
      plt.annotate(
         label, 
         xy = (x, y), xytext = (10, 12),
         textcoords = 'offset points', ha = 'right', va = 'bottom',
         bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white', alpha = 0.5))

   plt.show()
   
###
def main():
   Plot_Background(m_wavelengths,F_bgd,xs,ys)
   for i,j in zip(xs,ys):
      print i,j

if __name__ == '__main__': 
   main()
