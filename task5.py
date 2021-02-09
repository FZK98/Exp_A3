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
wl, n, k=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) 

#define the sellmeier refractive index function
def sellmeier(wavelength, a1,a2, a3,b1, b2,b3):
	tempFraction1 = a1*(wavelength**2)/(wavelength**2-b1**2)
	tempFraction2=a2*(wavelength**2)/(wavelength**2-b2**2)
	tempFraction3=a3*(wavelength**2)/(wavelength**2-b3**2)
	return np.sqrt(1+tempFraction1+tempFraction2+tempFraction3)

#fit the sellmeier function to the data with initial guesses of coefficients
popt, _ = sp.optimize.curve_fit(sellmeier, wl, n, [0.39, 0.45, 0.39,116.6, 116.6, -116]) 

#make list of refractive indices for each wl using sellmeier equation with fitted coefficients
sellmeiery = []
for i in wl:
	tempN=sellmeier(i, popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
	sellmeiery.append(tempN)
	
#plot results
plt.figure()
plt.plot(wl,n, label='data')
plt.plot(wl, sellmeiery, label='sellmeier 3 terms')
plt.legend()
plt.grid()
#a_BK7 = popt[0:2]
#b_BK7 = popt[3:5]
#print("the sellmeier coefficeints are A: ["+str(a_BK7)+"], B:["+str(b_BK7)+"].")

#plot difference in results
#plt.figure()
#plt.plot(wl, (sellmeiery-n))
#plt.grid()