# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 00:18:00 2021

@author: User
"""
import numpy as np 

#======= input parameters==========
#need angle, and wavelength and polarisation to model the interaction

# #allows the user to input variable values
# theta_i = input("Enter the angle of incidence: ")
# user_wl = input("Enter the wavelength of light: ")
# polarization = input("is your light s or p polarized? (enter s/p)")


#can maunally set the parameters (no real users)
theta_i = np.pi/2  # in radians - this is what numpy handles 
user_wl = 400 #nm
polarization = "s" #or "p"
N = 10 #this is the number of layers, not including the substrate
d = [] #list of the thickness of each layer (decide units)


#======= Angles for each layer
# n = [] #this will be a list with the (real) refractive index of each layer, will be n items long
def layer_angles(n):
    angle_list = [user_wl] #will become n+1 items long
    for i in range(len(n)):
        sin_tj = n[i]*np.sin(angle_list[i])/ n[i+1]
        theta_j = np.arcsin(sin_tj)
        angle_list.append(theta_j)
    return (angle_list) #returns a list of angle in each layer plus initial incident angle



#======= z wavenumber component for each layer

# kappa = [] # this will be a list of the imaginary components of the refractive index in each layer, will be n items long    
# angle_list = [] #this will be the list of angles in each layer found by the prev function 
def wavenumer_zi (wl, n, kappa, angle_list):
    k0 = 2*np.pi/wl
    im_kz = [0]
    re_kz = [k0]
    for i in range (len(n)):
        j = i + 1 
        new_im_kz = k0 * np.cos(angle_list[j]) * kappa[j]
        new_re_kz = k0 * np.cos(angle_list[j]) * n[j]
        im_kz.append(new_im_kz)
        re_kz.append(new_re_kz)
    return(re_kz, im_kz) #returns the list of the wavevector(z) in each layer 



#====== determine the P matrix : 
def P_elements(re_kz, im_kz, d):
    P_00_elements = []
    P_11_elements = []
    for i in range(len(d)):
        P_00_elements[i] = np.exp(complex(0,1)*re_kz[i]*d[i])*np.exp(-im_kz[i]*d[i]) #this will return a list of elements
        P_11_elements[i] = np.exp(complex(0,-1)*re_kz[i]*d[i])*np.exp(im_kz[i]*d[i])
    return(P_00_elements, P_11_elements)
    
#make the matrix, calling on each calculated element, will need to call on this N times
def Make_P_Matrix(P_00_element, P_11_element):
    P = np.array([[P_00_element, 0], 
                  [0, P_11_element]])
    return (P)

#====== determine the T matrix:
    
#first need to find the values of each matrix element. 

def XSP(re_kz, im_kz):
    
    

XSM 

XPP

XPM 

        











    