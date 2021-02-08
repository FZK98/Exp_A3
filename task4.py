# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 17:10:19 2021

@author: maria
"""
import numpy as np
import matplotlib.pyplot as plt
# =============================================================================
# import data from https://www.filmetrics.com/refractive-index-database/BK7/Float-Glass
# sellmeier coeffs from https://web.archive.org/web/20151011033820/http://www.lacroixoptical.com/sites/default/files/content/LaCroix%20Dynamic%20Material%20Selection%20Data%20Tool%20vJanuary%202015.xlsm
# =============================================================================
#wl=wavelength in nm, n=refractive index, k=wavevector
wl, n, k=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) 
b_BK7=[1.03961212, 0.231792344, 1.01046945]
c_BK7=[0.006000699, 0.020017914, 103.560653] #in micrometres^2

# =============================================================================
# fuction to find interpolated refractive indices
# =============================================================================
wl_integer = np.arange(np.min(wl),np.max(wl),1) #every integer wavelength from 330-2500
n_integer = np.interp(wl_integer,wl,n) #refractive index for each integer wavelength
def interpolatedN(wavelength):
	loc=np.where(wavelength==wl_integer)
	return n_integer[loc]

#plot the (interpolated) refractive index
plt.figure()
plt.plot(wl_integer, n_integer, label="interpolated")
#plt.plot(wl, n, 'x',label="raw data")
plt.grid()
plt.ylabel("refractive index n")
plt.xlabel("wavelength (nm)")
# =============================================================================
# sellmeier method
# =============================================================================
#a and b are empirically-found Sellmeier constants
def sellmeier(wavelength, b=b_BK7, c=c_BK7):
	wavelength = wavelength/1000 #convert wavelength to micrometres
	tempFraction = [] #the fraction part of the function to be summed over
	for i in range(len(b)):
		tempFraction.append((b[i]*wavelength**2)/(wavelength**2-c[i]**2))
	refractiveIndex = np.sqrt(1 + np.sum(tempFraction))
	return refractiveIndex

#make a list of refractive indeices for each integer wavelength using sellmeier method
sellmeierVals = []
for i in wl_integer:
	sellmeierVals.append(sellmeier(i))
#plot the sellmeier reflective index on the same figure as above
plt.plot(wl_integer, sellmeierVals, label="sellmeier")
plt.legend()

# =============================================================================
#  compare the two methods
# =============================================================================
def wlCompare(wavelength, b=b_BK7, c=c_BK7):
	interpTemp=interpolatedN(wavelength)
	sellmeierTemp = sellmeier(wavelength, b, c)
	return abs(interpTemp-sellmeierTemp)

#plot the difference between the two methods
methodDiff = []
for i in wl_integer:
	methodDiff.append(wlCompare(i))
plt.figure()
plt.plot(wl_integer, methodDiff)
plt.grid()
plt.ylabel("difference in refractive index")
plt.xlabel("wavelength (nm)")