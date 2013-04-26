import matplotlib.pyplot as plt
from Constants import *
from math import *
import numpy as np
import random
import time

def Plot_Bar_Chart(x_array, y_array, labels=""):
   if labels != "":
      wavelength_start, wavelength_step = (round(labels[0],2)), (round(labels[1],2))
   
      description = "SRE 0 = " + str(wavelength_start) + "A, SRE " + str(len(x_array)) + " = " + \
                 str(wavelength_start + len(x_array)*wavelength_step) + "A, SRE_step = " + str(wavelength_step) + "A"


      wavelength_ticks = [wavelength_start + (SRE_no * wavelength_step) for SRE_no in x_array]
   
   pos = np.arange(len(y_array))
   width = 1.0
   plt.ylabel('Flux /ergs ')
   plt.xlabel('Pixel #')
   plt.ylim(-0.5*max(y_array),1.1*max(y_array))
   plt.step(x_array,(y_array))
   plt.show()
