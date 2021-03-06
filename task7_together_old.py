# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:28:53 2021

@author: maria
"""
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
#load necessary data
#BK7_guess = np.loadtxt("BK7_guess.txt")
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_MgF2, n_MgF2, k_MgF2=np.loadtxt('MgF2.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold, n_gold, k_gold = np.loadtxt('Au.txt', delimiter = '\t', skiprows=1, unpack=True) #gold wavelength in nm
theta_i = 0#np.pi/2  # in radians - this is what numpy handles 
user_wl = 700 #nm 
polarization = "s" #or "p"
if polarization == "s":
	polarization = True 
else:
    polarization = False 
materials = ["air","MgF2","BK7"]
d = [0,0.0006,0]
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
#refractive index for each layer        
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

def P_elements(re_kz, im_kz, d):
    P_00_elements = []
    P_11_elements = []
    for i in range(1,len(d)-1): #don't need a P matrix for the substrate
        P_00_elements.append(np.exp(complex(0,1)*re_kz[i]*d[i])*np.exp(-im_kz[i]*d[i])) #this will return a list of elements
        P_11_elements.append(np.exp(complex(0,-1)*re_kz[i]*d[i])*np.exp(im_kz[i]*d[i]))
    return(P_00_elements, P_11_elements)
    
#make the matrix, calling on each calculated element, will need to call on this N-1 times
def Make_P_Matrix(P_00_element, P_11_element):
    P = np.array([[P_00_element, 0], 
                  [0, P_11_element]])
    return (P)
p_elements_list = P_elements(refractive_index_n, refractive_index_k, d)

#outputs all the P matrices in a 3D array 
#def All_P_Matrices(P_00_elements, P_11_elements): #there are N-1 elements 
#    Plist = np.empty_like((2,2))
#    for i in range(len(P_00_elements)):
#        Pi = np.array([[P_00_elements[i], 0], 
#                  [0, P_11_elements[i]]])
#        Plist = np.append(Plist, Pi)
#    P_arrays = np.reshape(Plist,(len(P_00_elements),2,2))
#    return (P_arrays)
def All_P_Matrices(P_00_elements, P_11_elements): #there are N-1 elements 
    Plist = []
    for i in range(len(P_00_elements)):
        Pi = np.array([[complex(P_00_elements[i]), 0], 
                  [0, complex(P_11_elements[i])]])
        Plist.append(Pi)
    return Plist

#making all the P matrices for each layer
#p_matrix_dictionary={} #stores P matrix for each layer
#p_matrix_list = []
#for i in range(len(p_elements_list[0])):
#    p_00_temp, p_11_temp = complex(p_elements_list[0][i]), complex(p_elements_list[1][i])
#    p_temp = Make_P_Matrix(p_00_temp, p_11_temp)
#    print(p_temp)
#    mat_name = "P"+str(i)
#    p_matrix_dictionary[mat_name] = p_temp
#    p_matrix_list.append(p_temp)
    
#first need to find the values of each matrix element. 
p_list = All_P_Matrices(p_elements_list[0], p_elements_list[1])

def XS_real(re_kz): #input is a list of real refractive indexes for each layer
    XSP=[]
    XSM=[]
    for i in range(len(re_kz)-1): #goes from the 0th vacuum layer to the N-1th (non-substrate) layer (the nth layer is the substrate and it doesnt have a bottom face)
        ki = float(re_kz[i])
        kj = float(re_kz[i+1])
        Term1 = np.sqrt(kj/ki)
        Term2 = np.sqrt(ki/kj)
        XSP.append(0.5*(Term1 + Term2))
        XSM.append(0.5*(Term1 - Term2))
    return(XSP, XSM)

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

def XP_real(re_kz, nr):#re_kz has N+1 elements, nr has N+1 elements
    XPP=[]
    XPM=[]
    for i in range (len(re_kz)-1):
        ki = float(re_kz[i]) 
        kj = float(re_kz[i+1])
        ni = nr[i]
        nj = nr[i+1]
        Term1 = (nj/ni)*np.sqrt(ki/kj)
        Term2 = (ni/nj)*np.sqrt(kj/ki)
        XPP.append(0.5*(Term1 + Term2))
        XPM.append(0.5*(Term1 - Term2))
    return(XPP, XPM)

    
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
    #T_arrays = np.reshape(Tlist, (len(XP),2,2))
    return(Tlist) #this is the 3d array - a list of 2x2 matrices
if polarization == True:
    list_XS_real = XS_real(wavenumber_list[0])
    list_XS_complex=XS(wavenumber_list[0], wavenumber_list[1])
    list_T_real = All_T_Matrices(list_XS_real[0], list_XS_real[1])
    list_T_complex = All_T_Matrices(list_XS_complex[0], list_XS_complex[1])
else:
    list_XP_real = XP_real(wavenumber_list[0], refractive_index_n)
    list_XP_complex = XP(wavenumber_list[0],wavenumber_list[1], refractive_index_n, refractive_index_k)
    list_T_real = All_T_Matrices(list_XP_real[0], list_XP_real[1])
    list_T_complex = All_T_Matrices(list_XP_complex[0], list_XP_complex[1])
	
def Make_M(re_kz, im_kz, d, polarization): 
    Plist=p_list
    Tlist=list_T_complex #can change to real
    #initializes M, therefore there are an equal amount of T and P matrices to be added to M 
    M = Tlist[0] 
    for i in range (len(Tlist)-1): #iterates through the process N-1 times
        M = np.matmul(Plist[i], M) #add the P of the ith layer 
        M = np.matmul(Tlist[i+1], M) #adds the T of the ith + 1 layer 
    return (M) #this is the fully calculated M matrix
matrix_M = Make_M(wavenumber_list[0], wavenumber_list[1], d, polarization)    
def rt_solver(M):
    r = - M[1,0] / M[1, 1]
    t =   M[0, 0] + (M[0, 1] * r )
    return(r, t)
answer_r_t=rt_solver(matrix_M)          
print(abs(answer_r_t[0])**2)  