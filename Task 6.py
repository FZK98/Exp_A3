# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:29:03 2021

@author: User
"""

import numpy as np 
import matplotlib.pyplot as plt 

Gold_wl , Gold_n = np.loadtxt("Real_Gold_Data.txt", delimiter='\t',skiprows=1, unpack=True)
Gold_wl , Gold_k = np.loadtxt("Im_Gold_Data.txt", delimiter='\t',skiprows=1, unpack=True)

wl = Gold_wl 
n = Gold_n
k = Gold_k 


# =============================================================================
# fuction to find interpolated refractive indices
# =============================================================================
wl = wl*1000 # convert wavelengths into nm 
wl_integer = np.arange(np.min(wl),np.max(wl),1) #every integer wavelength from 330-2500
n_integer = np.interp(wl_integer, wl, n) #refractive index for each integer wavelength
k_integer = np.interp(wl_integer, wl, k)
# def interpolatedN(wavelength):
# 	loc=np.where(wavelength==wl_integer)
# 	return n_integer[loc], k_integer[loc]

#plot the (interpolated) refractive index
plt.figure()
plt.plot(wl_integer, n_integer, label="interpolated n")
plt.plot(wl_integer, k_integer, label="interpolated k")
plt.legend()
#plt.plot(wl, n, 'x',label="raw data")
plt.grid()
plt.ylabel("refractive index n")
plt.xlabel("wavelength (nm)")

# =============================================================================
# find the complex refractive index at some input wavelength
# =============================================================================
#input some wavelength into the function, returns the real and complex values of the refractive index
def interpolated_ncomp(wavelength):
	loc=np.where(wavelength==wl_integer)
	return n_integer[loc], k_integer[loc]

