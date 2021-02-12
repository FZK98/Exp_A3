# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:28:53 2021

@author: maria
"""
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
#load necessary data
BK7_guess = np.loadtxt("BK7_guess.txt")
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold, n_gold, k_gold = np.loadtxt('Au.txt', delimiter = '\t', skiprows=1, unpack=True) #gold wavelength in nm
#can maunally set the parameters (no real users)
theta_i = np.pi/2  # in radians - this is what numpy handles 
user_wl = 400 #nm
polarization = "s" #or "p"
N = 2 #this is the number of layers, not including the substrate
materials = ["air", "BK7"] #which materials constitute the layers
d = [10, 4] #list of the thickness of each layer (decide units)

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
def interpolate_n(wl, wl_data, n_data):
	wl_integer = np.arange(np.min(wl_data),np.max(wl_data),1) #every integer wavelength from 330-2500
	n_integer = np.interp(wl_integer, wl_data, n_data) #refractive index for each integer wavelength
	loc=np.where(wl==wl_integer)
	return n_integer[loc]
def interpolate_k(wl, wl_data, k_data):
	wl_integer = np.arange(np.min(wl_data),np.max(wl_data),1) #every integer wavelength from 330-2500
	k_integer = np.interp(wl_integer, wl_data, k_data)
	loc=np.where(wl==wl_integer)
	return  k_integer[loc]

#function to return n and k for any wavelength, given the data and materialtype
def get_refractiveIndex_n(wl, wl_data, n_data, initialGuess, transparent=True):
	if transparent == True:
		coeffs_temp = sellmeier_coefficients(wl_data, n_data, initialGuess)
		sellmeier_n = sellmeier_fit(wl, coeffs_temp[0],coeffs_temp[1],coeffs_temp[2],coeffs_temp[3],coeffs_temp[4],coeffs_temp[5])
		return sellmeier_n
	else:
		opaque_n = interpolate_n(wl, wl_data, n_data)
		return opaque_n
def get_refractiveIndex_k(wl, wl_data, k_data):
	return float(interpolate_k(wl, wl_data, k_data))

#refractive index for each layer		
refractive_index_n = []
refractive_index_k = []
for i in range(N):
	if materials[i] == "air":
		refractive_index_n.append(1)
		refractive_index_k.append(0)
	elif materials[i] == "BK7":
		refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_BK7, n_BK7, BK7_guess, True))
		refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_BK7, k_BK7))
	elif materials[i] == "gold":
		refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_gold, n_gold, 0, False))
		refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_gold, k_gold))
#	elif materials[i] == "MgF2":
#		refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_MgF2, n_MgF2, MgF2_guess, True))
#		refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_MgF2, k_MgF2))
	else:
		refractive_index_n.append("unknown material")
		refractive_index_k.append("unknown material")
		
#======= Angles for each layer
def layer_angles(n):
    angle_list = [theta_i] #will become n+1 items long
    for i in range(len(n)-1):
        sin_tj = n[i]*np.sin(angle_list[i])/ n[i+1]
        theta_j = np.arcsin(sin_tj)
        angle_list.append(theta_j)
    return (angle_list) #returns a list of angle in each layer plus initial incident angle

#======= z wavenumber component for each layer
def wavenumber_zi (wl, n, kappa, angle_list):
    k0 = 2*np.pi/wl
    im_kz = [0]
    re_kz = [k0]
    for i in range(len(n)-1):
        j = i + 1 
        new_im_kz = k0 * np.cos(angle_list[j]) * kappa[j]
        new_re_kz = k0 * np.cos(angle_list[j]) * n[j]
        im_kz.append(new_im_kz)
        re_kz.append(new_re_kz)
    return(re_kz, im_kz) #returns the list of the wavevector(z) in each layer 
	
angle_list = layer_angles(refractive_index_n) #this will be the list of angles in each layer found by the prev function 
wavenumber_list = wavenumber_zi(user_wl, refractive_index_n, refractive_index_k, angle_list)
# =============================================================================
# we now have the following lists of data:
	# real refractive index: refractive_index_n
	# imaginary refractive index: refractive_index_k
	# angle for each layer: angle_list
	# kz for each layer: wavenumber_list
# now to build the matrices
# =============================================================================
def P_elements(re_kz, im_kz, d):
	P_00_element = np.exp(complex(0,1)*re_kz*d)*np.exp(-im_kz*d) 
	P_11_element = np.exp(complex(0,-1)*re_kz*d)*np.exp(im_kz*d)
	return(P_00_element, P_11_element)
    
#make the matrix, calling on each calculated element, will need to call on this N times
def Make_P_Matrix(P_00_element, P_11_element):
    P = np.array([[P_00_element, 0], 
                  [0, P_11_element]])
    return (P)

p_matrix_dictionary={}
for i in range(len(d)):
	p_00_temp, p_11_temp = P_elements(wavenumber_list[i][0], wavenumber_list[i][1], d[i])
	p_temp = Make_P_Matrix(p_00_temp, p_11_temp)
	mat_name = "P"+str(i)
	p_matrix_dictionary[mat_name] = p_temp