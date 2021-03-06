# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 03:09:48 2021

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import minimize

#load necessary data
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_MgF2, n_MgF2, k_MgF2=np.loadtxt('MgF2.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold, n_gold, k_gold = np.loadtxt('Au.txt', delimiter = '\t', skiprows=1, unpack=True) #gold wavelength in nm

def transfer_matrix(theta_i, user_wl, polarisation,materials,d):
    #theta_i = 0
    #user_wl = 700 #nm
    polarization = "s" #or "p"
    if polarization == "s":
        polarization = True 
    else:
        polarization = False 
    #materials = ["vacuum","MgF2","BK7"] #which materials constitute the layers
    #d = [0,10,0] #list of the thickness of each layer (decide units). zero for vacuum and substrate for length consistency
    
#interpolation functions for a given wavelength - returns complex n for non-transparent medium
    def interpolate_n(wl, wl_data, n_data):
        wl_integer = np.arange(np.rint(np.min(wl_data)),np.rint(np.max(wl_data)),1) #every integer wavelength from 330-2500
        n_integer = np.interp(wl_integer, wl_data, n_data) #refractive index for each integer wavelength
        loc=np.where(wl==wl_integer)
        return n_integer[loc]
    def interpolate_k(wl, wl_data, k_data):
        wl_integer = np.arange(np.rint(np.min(wl_data)),np.rint(np.max(wl_data)),1) #every integer wavelength from 330-2500
        k_integer = np.interp(wl_integer, wl_data, k_data)
        loc=np.where(wl==wl_integer)
        return  k_integer[loc]

#function to return n and k for any wavelength, given the data and materialtype
    def get_refractiveIndex_n(wl, wl_data, n_data):
        opaque_n = interpolate_n(wl, wl_data, n_data)
        return opaque_n
    def get_refractiveIndex_k(wl, wl_data, k_data):
        return float(interpolate_k(wl, wl_data, k_data))

#refractive index for each layer        
    refractive_index_n = []
    refractive_index_k = []
    for i in range(len(materials)): #vacuum needs to be included for these
        if materials[i] == "air":
            refractive_index_n.append(1)
            refractive_index_k.append(0)
        elif materials[i] == "BK7":
            refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_BK7, n_BK7))
            refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_BK7, k_BK7))
        elif materials[i] == "gold":
            refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_gold, n_gold))
            refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_gold, k_gold))
        elif materials[i] == "MgF2":
            refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_MgF2, n_MgF2))
            refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_MgF2, k_MgF2))
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
        re_kz = [k0*np.cos(theta_i)]
        for i in range(len(n)-1):
            j = i + 1 
            new_im_kz = k0 * np.cos(angle_list[j]) * kappa[j]
            new_re_kz = k0 * np.cos(angle_list[j]) * n[j]
            im_kz.append(new_im_kz)
            re_kz.append(new_re_kz)
        return(re_kz, im_kz) #returns the list of the wavevector(z) in each layer 
        
#======= construct lists of angles and wavenumbers
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
    
#======= functions to build P matrices
    def P_elements(re_kz, im_kz, d):
        P_00_elements = []
        P_11_elements = []
        for i in range(1,len(d)-1): #don't need a P matrix for the substrate
            P_00_elements.append(np.exp(complex(0,1)*re_kz[i]*d[i])*np.exp(-im_kz[i]*d[i])) #this will return a list of elements
            P_11_elements.append(np.exp(complex(0,-1)*re_kz[i]*d[i])*np.exp(im_kz[i]*d[i]))
        return(P_00_elements, P_11_elements)
       
    def All_P_Matrices(P_00_elements, P_11_elements): #there are N-1 elements 
        Plist = []
        for i in range(len(P_00_elements)):
            Pi = np.array([[complex(P_00_elements[i]), 0], 
                  [0, complex(P_11_elements[i])]])
            Plist.append(Pi)
        return Plist
#======= build P matrices
    p_elements_list = P_elements(wavenumber_list[0], wavenumber_list[1], d)
    p_list = All_P_Matrices(p_elements_list[0], p_elements_list[1])

#======= functions to build T matrices
    def XS(re_kz, im_kz): #input is a list refractive indexes for each layer (0 and substrate inclusive)
        XSP=[]
        XSM=[]
        for i in range(len(re_kz)-1): 
            ki = complex(re_kz[i], im_kz[i])
            kj = complex(re_kz[i+1], im_kz[i+1])
            Term1 = np.sqrt(kj/ki)
            Term2 = np.sqrt(ki/kj)
            XSP.append(0.5*(Term1 + Term2))
            XSM.append(0.5*(Term1 - Term2))
        return(XSP, XSM)
    
    def XP(re_kz, im_kz, nr, kappa):
        XPP=[]
        XPM=[]
        for i in range (len(re_kz)-1):
            ki = complex(re_kz[i], im_kz[i])
            kj = complex(re_kz[i+1], im_kz[i+1])
            ni = complex(nr[i], kappa[i])
            nj = complex(nr[i+1], kappa[i+1])
            Term1 = (nj/ni)*np.sqrt(ki/kj)
            Term2 = (ni/nj)*np.sqrt(kj/ki)
            XPP.append(0.5*(Term1 + Term2))
            XPM.append(0.5*(Term1 - Term2))
        return(XPP, XPM)

    def Make_T(XP, XM): #will need to be called on N times 
        T = np.array([[XP,XM], 
                      [XM, XP]])
        return(T)

    def All_T_Matrices(XP, XM):
        Tlist = []
        for i in range(len(XP)):
            Ti = np.array([[XP[i], XM[i]],
                       [XM[i],XP[i]]])
            Tlist.append(Ti)   
        return(Tlist) #this is the 3d array - a list of 2x2 matrices

#======= build T matrices
    if polarization == True:
        list_XS_complex=XS(wavenumber_list[0], wavenumber_list[1])
        list_T_complex = All_T_Matrices(list_XS_complex[0], list_XS_complex[1])
    else:
        list_XP_complex = XP(wavenumber_list[0],wavenumber_list[1], refractive_index_n, refractive_index_k)
        list_T_complex = All_T_Matrices(list_XP_complex[0], list_XP_complex[1])

#======= functions to build M matrix
    def Make_M(re_kz, im_kz, d, polarization): 
        Plist=p_list
        Tlist=list_T_complex 
    #initializes M, therefore there are an equal amount of T and P matrices to be added to M 
        M = Tlist[0] 
        for i in range (len(Tlist)-1): #iterates through the process N-1 times
            M = np.matmul(Plist[i], M) #add the P of the ith layer 
            M = np.matmul(Tlist[i+1], M) #adds the T of the ith + 1 layer 
        return (M) #this is the fully calculated M matrix
    
#======= build M matrix
    matrix_M = Make_M(wavenumber_list[0], wavenumber_list[1], d, polarization)
#======= function to return r and t
    def rt_solver(M):
        r = - M[1,0] / M[1, 1]
        t =   M[0, 0] + (M[0, 1] * r )
        return(r, t)
#======= return r and t    
    answer_r_t=rt_solver(matrix_M)
    return(answer_r_t)
    
# =============================================================================
# now to run the transfer matrix for different wavelength and different layer thcikness to find minimised reflectivity    
# =============================================================================
plain_glass = transfer_matrix(0, 600, "s",materials = ["air","gold","BK7"],d=[0,0,0])
plain_glass_R=abs(plain_glass[0])**2
plain_glass_T=abs(plain_glass[1])**2

wl_test=np.arange(330, 2470, 1)
wl_test_R = []
wl_test_T = []
wl_test_A = []
d_test = np.arange(0, 100, 1)
d_test_R = []
d_test_T = []
d_test_T_increasing=[] #interpolation function needs x values to be increasing
d_test_A1 = []
d_test_A2 = []
d_test_A3 = []
d_test_A4 = []
d_test_A5 = []
d_test_A6 = []
for i in wl_test:
    rtemp_wl, ttemp_wl = transfer_matrix(0, i, "s",materials = ["air","gold","BK7"],d=[0,50,0])
    wl_test_R.append(abs(rtemp_wl)**2)
    wl_test_T.append(abs(ttemp_wl)**2)
    wl_test_A.append(1 - (abs(rtemp_wl)**2) - (abs(ttemp_wl)**2) )
for i in d_test:
	#set the wavelengths to whatever needed
    rtemp_d1, ttemp_d1 = transfer_matrix(0, 400, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    rtemp_d2, ttemp_d2 = transfer_matrix(0, 500, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    rtemp_d3, ttemp_d3 = transfer_matrix(0, 800, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    rtemp_d4, ttemp_d4 = transfer_matrix(0, 1000, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    rtemp_d5, ttemp_d5 = transfer_matrix(0, 2000, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    rtemp_d6, ttemp_d6 = transfer_matrix(0, 2400, "s", materials = ["air","gold","BK7"],d=[0,i,0])
    d_test_R.append(abs(rtemp_d1)**2)
    d_test_T.append(abs(ttemp_d1)**2)
    d_test_A1.append(1 - (abs(rtemp_d1)**2) - (abs(ttemp_d1)**2) )
    d_test_A2.append(1 - (abs(rtemp_d2)**2) - (abs(ttemp_d2)**2) )
    d_test_A3.append(1 - (abs(rtemp_d3)**2) - (abs(ttemp_d3)**2) )
    d_test_A4.append(1 - (abs(rtemp_d4)**2) - (abs(ttemp_d4)**2) )
    d_test_A5.append(1 - (abs(rtemp_d5)**2) - (abs(ttemp_d5)**2) )
    d_test_A6.append(1 - (abs(rtemp_d6)**2) - (abs(ttemp_d6)**2) )

# =============================================================================
# find depth of gold for 10^3 attenuation
# =============================================================================
for i in range(1,len(d_test_T)):
    d_test_T_increasing.append(d_test_T[-i])
d_test_T_increasing.append(d_test_T[0])
three_orders_attenuation=np.interp(plain_glass_T/1000, d_test_T_increasing, d_test) #this is the index for the reverse list
three_orders_attenuation=len(d_test)-three_orders_attenuation #corrected back value, only works if line 202 has spacing of 1
print("depth of gold needed to attenuate by 3 orders of magnitude is: "+str(three_orders_attenuation)+"nm.")

# =============================================================================
# find the wavelength and depth of layer that provide minimum R and maximum T
# =============================================================================
def wl_minimise_R(wl_test, wl_test_R):
	minR=np.min(wl_test_R)
	loc=np.where(wl_test_R==minR)
	return(wl_test[loc])
def d_minimise_R(d_test, d_test_R):
    minR=np.min(d_test_R)
    loc=np.where(d_test_R==minR)
    return(d_test[loc])

def wl_maximize_T(wl_test, wl_test_T):
	maxT=np.max(wl_test_T)
	loc=np.where(wl_test_T==maxT)
	return(wl_test[loc])
def d_maximize_T(d_test, d_test_T):
    maxT=np.max(d_test_T)
    loc=np.where(d_test_T==maxT)
    return(d_test[loc])

def wl_maximize_A(wl_test, wl_test_A):
	maxA=np.max(wl_test_A)
	loc=np.where(wl_test_A==maxA)
	return(wl_test[loc])
def d_maximize_A(d_test, d_test_A):
    maxA=np.max(d_test_A)
    loc=np.where(d_test_A==maxA)
    return(d_test[loc])


# =============================================================================
# plot reflectivity against wavelength and thickness
# =============================================================================
plt.figure()
plt.plot(wl_test, wl_test_R, label="optimal wavelength = "+str(wl_minimise_R(wl_test,wl_test_R)))
plt.title("wavelength against R")
#plt.legend()
plt.grid()
plt.xlabel("wavelength [nm]")
plt.ylabel("Reflectance R")
plt.figure()
plt.plot(d_test, d_test_R, label="optimal depth = "+str(d_minimise_R(d_test,d_test_R)))
#plt.legend()
plt.grid()
plt.xlabel("layer depth [nm]")
plt.ylabel("Reflectance R")
plt.title("depth of layer index against R")
#
## =============================================================================
## plot tranmission against wavelength and thickness
## =============================================================================
plt.figure()
plt.plot(wl_test, wl_test_T, label="optimal wavelength = "+str(wl_maximize_T(wl_test,wl_test_T)))
plt.title("wavelength against T")
#plt.legend()
plt.grid()
plt.xlabel("wavelength [nm]")
plt.ylabel("Transmittance T")
plt.figure()
plt.plot(d_test, d_test_T, label="optimal depth = "+str(d_maximize_T(d_test,d_test_T)))
#plt.legend()
plt.grid()
plt.xlabel("layer depth [nm]")
plt.ylabel("Transmittance T")
plt.title("depth of layer index against T")

# =============================================================================
# plot absobance against wavelength and thickness
# =============================================================================
plt.figure()
plt.plot(wl_test, wl_test_A, label="optimal wavelength = "+str(wl_maximize_A(wl_test,wl_test_A)))
plt.title("wavelength against A")
#plt.legend()
plt.grid()
plt.xlabel("wavelength [nm]")
plt.ylabel("Absorbance A")
plt.figure()
plt.plot(d_test, d_test_A1, label="400nm optimal depth = "+str(d_maximize_A(d_test,d_test_A1)))
plt.plot(d_test, d_test_A2, label="500nm optimal depth = "+str(d_maximize_A(d_test,d_test_A2)))
plt.plot(d_test, d_test_A3, label="600nm optimal depth = "+str(d_maximize_A(d_test,d_test_A3)))
plt.plot(d_test, d_test_A4, label="700nm optimal depth = "+str(d_maximize_A(d_test,d_test_A4)))
plt.plot(d_test, d_test_A5, label="800nm optimal depth = "+str(d_maximize_A(d_test,d_test_A5)))
#plt.plot(d_test, d_test_A6, label="2400nm optimal depth = "+str(d_maximize_A(d_test,d_test_A6)))
plt.legend()
plt.grid()
plt.xlabel("layer depth [nm]")
plt.ylabel("Absorbance A")
plt.title("depth of layer index against A")
