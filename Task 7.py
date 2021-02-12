# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 00:18:00 2021

@author: User
"""
import numpy as np 

# ============================================================================
# input parameters
# =============================================================================

#need angle, and wavelength and polarisation to model the interaction

# #allows the user to input variable values
# theta_i = input("Enter the angle of incidence: ")
# user_wl = input("Enter the wavelength of light: ")
# polarization = input("is your light s or p polarized? (enter s/p)")


#can maunally set the parameters (no real users)
theta_i = np.pi/2  # in radians - this is what numpy handles 
user_wl = 400 #nm
polarization = "s" #or "p"
if polarization == "s":
    polarization = True 
else:
    polarization = False 

N = 10 #this is the number of layers, not including the substrate
d = [] #list of the thickness of each layer (decide units)



# ============================================================================
# Angles for each layer
# =============================================================================

# nr = [] #this will be a list with the (real) refractive index of each layer, will be n items long
def layer_angles(nr):
    angle_list = [user_wl] #will become n+1 items long
    for i in range(len(nr)):
        sin_tj = nr[i]*np.sin(angle_list[i])/ nr[i+1]
        theta_j = np.arcsin(sin_tj)
        angle_list.append(theta_j)
    return (angle_list) #returns a list of angle in each layer plus initial incident angle




# ============================================================================
# z wavenumber component for each layer
# =============================================================================

# kappa = [] # this will be a list of the imaginary components of the refractive index in each layer, will be n items long    
# angle_list = [] #this will be the list of angles in each layer found by the prev function 
def wavenumer_zi (wl, nr, kappa, angle_list):
    k0 = 2*np.pi/wl
    im_kz = [0]
    re_kz = [k0]
    for i in range (len(nr)):
        j = i + 1 
        new_im_kz = k0 * np.cos(angle_list[j]) * kappa[j]
        new_re_kz = k0 * np.cos(angle_list[j]) * nr[j]
        im_kz.append(new_im_kz)
        re_kz.append(new_re_kz)
    return(re_kz, im_kz) #returns the list of the wavevector(z) in each layer 


# ============================================================================
# determine the P matrix : 
# =============================================================================
 
def P_elements(re_kz, im_kz, d):
    P_00_elements = []
    P_11_elements = []
    for i in range(len(d)):
        P_00_elements[i] = np.exp(complex(0,1)*re_kz[i]*d[i])*np.exp(-im_kz[i]*d[i]) #this will return a list of elements
        P_11_elements[i] = np.exp(complex(0,-1)*re_kz[i]*d[i])*np.exp(im_kz[i]*d[i])
    return(P_00_elements, P_11_elements)
    
#make the matrix, calling on each calculated element, will need to call on this N-1 times
def Make_P_Matrix(P_00_element, P_11_element):
    P = np.array([[P_00_element, 0], 
                  [0, P_11_element]])
    return (P)

def All_P_Matrices(P_00_elements, P_11_elements): #there are N-1 elements 
    Plist = np.empty_like((2,2))
    for i in range(len(P_00_elements)):
        Pi = np.array([[P_00_elements[i], 0], 
                  [0, P_11_elements[i]]])
        Plist = np.append(Plist, Pi)
    P_arrays = np.reshape(Plist,(len(P_00_elements),2,2))

    return (P_arrays)

# ============================================================================
# determine the T matrix:
# =============================================================================

#first need to find the values of each matrix element. 
def XS_real(re_kz): #input is a list of real refractive indexes for each layer 
    for i in range(len(re_kz)-1): #goes from the 0th vacuum layer to the N-1th (non-substrate) layer (the nth layer is the substrate and it doesnt have a bottom face)
        ki = re_kz[i]
        kj = re_kz[i+1]
        Term1 = np.sqrt(kj/ki)
        Term2 = np.sqrt(ki/kj)
        XSP = 0.5*(Term1 + Term2)
        XSM = 0.5*(Term1 - Term2)
    return(XSP, XSM)

def XS(re_kz, im_kz): #input is a list refractive indexes for each layer (0 and substrate inclusive)
    for i in range(len(re_kz)-1): 
        ki = complex(re_kz[i], im_kz[i])
        kj = complex(re_kz[i+1], im_kz[i+1])
        Term1 = np.sqrt(kj/ki)
        Term2 = np.sqrt(ki/kj)
        XSP = 0.5*(Term1 + Term2)
        XSM = 0.5*(Term1 - Term2)
    return(XSP, XSM)

def XP_real(re_kz, nr):#re_kz has N+1 elements, nr has N+1 elements
    for i in range (len(re_kz)):
        ki = re_kz[i] 
        kj = re_kz[i+1]
        ni = nr[i]
        nj = nr[i+1]
        Term1 = (nj/ni)*np.sqrt(ki/kj)
        Term2 = (ni/nj)*np.sqrt(kj/ki)
        XPP = 0.5*(Term1 + Term2)
        XPM = 0.5*(Term1 - Term2)
    return(XPP, XPM)

    
def XP(re_kz, im_kz, nr, kappa):
    for i in range (len(re_kz)):
        ki = complex(re_kz[i], im_kz[i])
        kj = complex(re_kz[i+1], im_kz[i+1])
        ni = complex(nr[i], kappa[i])
        nj = complex(nr[i+1], kappa[i+1])
        Term1 = (nj/ni)*np.sqrt(ki/kj)
        Term2 = (ni/nj)*np.sqrt(kj/ki)
        XPP = 0.5*(Term1 + Term2)
        XPM = 0.5*(Term1 - Term2)
    return(XPP, XPM)

def Make_T(XP, XM): #will need to be called on N times 
    T = np.array([[XP,XM], 
                  [XM, XP]])
    return(T)
    
def All_T_Matrices(XP, XM):
    Tlist = np.empty_like((2, 2))
    for i in range(len(XP)):
        Ti = np.array([[XP, XM],
                       [XM,XP]])
        Tlist = np.append(Tlist, Ti)
    T_arrays = np.reshape(Tlist, (len(XP),2,2))
    return(T_arrays) #this is the 3d array - a list of 2x2 matrices


# ============================================================================
# Determine the M Matrix 
# =============================================================================
# def Make_M(N, polarization): #N is the integer number of layers we have including the substrate
    
            
 











    
