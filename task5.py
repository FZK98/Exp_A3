# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 10:12:57 2021

@author: maria
"""
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#load data
wl, n, k=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm

wl_integer = np.arange(np.min(wl),np.max(wl),1) #every integer wavelength from 330-2500
n_integer = np.interp(wl_integer,wl,n) #refractive index for each integer wavelength

#define the sellmeier refractive index function
def sellmeier_fit(wavelength, a1,a2, a3,b1, b2,b3):
	wavelength=wavelength/1000 #convert wavelength to micrometres to make starting estimates work
	tempFraction1 = a1*(wavelength**2)/(wavelength**2-b1)
	tempFraction2=a2*(wavelength**2)/(wavelength**2-b2)
	tempFraction3=a3*(wavelength**2)/(wavelength**2-b3)
	return np.sqrt(1+tempFraction1+tempFraction2+tempFraction3)

#fit the sellmeier function to the data with initial guesses of coefficients
popt, _ = sp.optimize.curve_fit(sellmeier_fit, wl, n, [1.04, 0.23, 1.01, 0.006, 0.02, 103.5]) #starting estimates from wikipedia

#make list of refractive indices for each wl using sellmeier equation with fitted coefficients
sellmeiery = []
for i in wl:
	tempN=sellmeier_fit(i, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
	sellmeiery.append(tempN)
	
#plot results
plt.figure()
plt.plot(wl,n, 'x', label='data')
plt.plot(wl_integer, n_integer, 'o', label = "interpolated")
plt.plot(wl, sellmeiery, label='sellmeier')
plt.ylabel("refractive index n")
plt.xlabel("wavelength (nm)")
plt.legend()
plt.grid()
a_BK7 = popt[0:3]
b_BK7 = popt[3:6]
print("the sellmeier coefficeints are A: ["+str(a_BK7)+"], B:["+str(b_BK7)+"].")

#plot difference in results
plt.figure()
plt.plot(wl, (sellmeiery-n))
plt.ylabel("difference in refractive index")
plt.xlabel("wavelength (nm)")
plt.grid()