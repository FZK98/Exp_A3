# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:09:59 2021

@author: maria
"""
import numpy as np
import scipy as sp
#user input variables
angleIncidence = 0
wavelength = 0
#load necessary data
BK7_guess = np.loadtxt("BK7_guess.txt")
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold , n_gold = np.loadtxt("Real_Gold_Data.txt", delimiter='\t',skiprows=1, unpack=True) #gold wavelength in um
_ , k_gold = np.loadtxt("Im_Gold_Data.txt", delimiter='\t',skiprows=1, unpack=True)

#the sellmeier equation for given wavelength and coefficients - returns n for transparent medium
def sellmeier_fit(wl, a1,a2, a3,b1, b2,b3):
	wl=wl/1000 #convert wavelength to micrometres to make starting estimates work
	tempFraction1 = a1*(wl**2)/(wl**2-b1)
	tempFraction2=a2*(wl**2)/(wl**2-b2)
	tempFraction3=a3*(wl**2)/(wl**2-b3)
	return np.sqrt(1+tempFraction1+tempFraction2+tempFraction3)

#fits the sellmeier equation to given data to find coefficients for transparent material 
def sellmeier_coefficients(wl_data, n_data, initialGuess):
	popt, _ = sp.optimize.curve_fit(sellmeier_fit, wl_data, n_data, initialGuess)
	return(popt) #return fitted Sellemeier coefficients
	
#interpolation functions for a given wavelength - returns complex n for non-transparent medium
def interpolate_n(wl, wl_data, n_data, k_data):
	wl_data = wl_data*1000 # convert wavelengths into nm 
	wl_integer = np.arange(np.min(wl_data),np.max(wl_data),1) #every integer wavelength from 330-2500
	n_integer = np.interp(wl_integer, wl_data, n_data) #refractive index for each integer wavelength
	k_integer = np.interp(wl_integer, wl_data, k_data)
	loc=np.where(wavelength==wl_integer)
	return n_integer[loc], k_integer[loc]

#function to return n and k for any wavelength, given the data and materialtype
def get_refractiveIndex(wl, wl_data, n_data, k_data, initialGuess, transparent=True):
	if transparent == True:
		coeffs_temp = sellmeier_coefficients(wl_data, n_data, initialGuess)
		print(coeffs_temp)
		return(sellmeier_fit(wl, coeffs_temp[0],coeffs_temp[1],coeffs_temp[2],coeffs_temp[3],coeffs_temp[4],coeffs_temp[5]))
	else:
		return(interpolate_n(wl, wl_data, n_data, k_data))
