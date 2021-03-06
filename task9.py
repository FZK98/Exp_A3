# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:56:15 2021

@author: maria
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:19:09 2021

@author: maria
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:28:53 2021

@author: maria
"""
import numpy as np
import matplotlib.pyplot as plt
#load necessary data
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_MgF2, n_MgF2, k_MgF2=np.loadtxt('MgF2.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold, n_gold, k_gold = np.loadtxt('Au.txt', delimiter = '\t', skiprows=1, unpack=True) #gold wavelength in nm

def transfer_matrix(theta_i, user_wl, polarization,materials,d, n_ideal=0):
    
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
            if n_ideal==0:
                refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_MgF2, n_MgF2))
                refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_MgF2, k_MgF2))
            else:
                refractive_index_n.append(n_ideal)
                refractive_index_k.append(0)
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
    def wavenumber_zi (wl, n, kappa, angle_list, theta_i):
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
    wavenumber_list = wavenumber_zi(user_wl, refractive_index_n, refractive_index_k, angle_list, theta_i)
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
    if polarization == "s":
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
plain_glass = transfer_matrix(0, 400, "s", ["air","BK7"],[0,0],0)
plain_glass_R=abs(plain_glass[0])**2

wl_test=np.arange(330, 900, 1)
wl_test_R = []
n_test=np.arange(0.5,2.5,0.05)
n_test_R=[]
d_test = np.arange(5, 500, 1)
d_test_R = []
theta_test = np.linspace(0,np.pi/2, 50)
theta_test_R_s=[]
theta_test_R_p=[]
theta_test_R_s1=[]
theta_test_R_p1=[]
theta_test_R_s2=[]
theta_test_R_p2=[]
theta_test_R_s3=[]
theta_test_R_p3=[]
theta_test_R_s4=[]
theta_test_R_p4=[]
for i in n_test:
    rtemp_r, ttemp_r = transfer_matrix(0, 400, "s",materials = ["air","MgF2","BK7"],d=[0,400/(4*i),0],  n_ideal=i)
    n_test_R.append(abs(rtemp_r)**2)
for i in wl_test:
    rtemp_wl, ttemp_wl = transfer_matrix(0, i, "s",materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    wl_test_R.append(abs(rtemp_wl)**2)
for i in d_test:
    rtemp_d, ttemp_d = transfer_matrix(0, 400, "s", materials = ["air","MgF2","BK7"],d=[0,i,0])
    d_test_R.append(abs(rtemp_d)**2)
for i in theta_test:
    rtemp_theta_s, ttemp_theta_s = transfer_matrix(i, 350, "s", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_s.append(abs(rtemp_theta_s)**2)
    rtemp_theta_p, ttemp_theta_p = transfer_matrix(i, 350, "p", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_p.append(abs(rtemp_theta_p)**2)
    rtemp_theta_s, ttemp_theta_s = transfer_matrix(i, 450, "s", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_s1.append(abs(rtemp_theta_s)**2)
    rtemp_theta_p, ttemp_theta_p = transfer_matrix(i, 450, "p", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_p1.append(abs(rtemp_theta_p)**2)
    rtemp_theta_s, ttemp_theta_s = transfer_matrix(i, 600, "s", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_s2.append(abs(rtemp_theta_s)**2)
    rtemp_theta_p, ttemp_theta_p = transfer_matrix(i, 600, "p", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_p2.append(abs(rtemp_theta_p)**2)
    rtemp_theta_s, ttemp_theta_s = transfer_matrix(i, 700, "s", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_s3.append(abs(rtemp_theta_s)**2)
    rtemp_theta_p, ttemp_theta_p = transfer_matrix(i, 700, "p", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_p3.append(abs(rtemp_theta_p)**2)
    rtemp_theta_s, ttemp_theta_s = transfer_matrix(i, 800, "s", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_s4.append(abs(rtemp_theta_s)**2)
    rtemp_theta_p, ttemp_theta_p = transfer_matrix(i, 800, "p", materials = ["air","MgF2","BK7"],d=[0,400/(4*1.38),0])
    theta_test_R_p4.append(abs(rtemp_theta_p)**2)

# =============================================================================
# find the wavelength and depth of layer that provide minimum
# =============================================================================
def wl_minimise(wl_test, wl_test_R):
    minR=np.min(wl_test_R)
    loc=np.where(wl_test_R==minR)
    return(wl_test[loc])
def d_minimise(d_test, d_test_R):
    minR=np.min(d_test_R)
    loc=np.where(d_test_R==minR)
    return(d_test[loc])
def n_minimise(n_test, n_test_R):
    minR=np.min(n_test_R)
    loc=np.where(n_test_R==minR)
    return(n_test[loc])
def theta_minimise(theta_test, test_R):
    minR=np.min(test_R)
    loc=np.where(test_R==minR)
    return(theta_test[loc])
    
    
# =============================================================================
# plot reflectivity against parameter
# =============================================================================
#wavelength
plt.figure()
plt.plot(wl_test, wl_test_R, label="optimal wavelength = "+str(wl_minimise(wl_test,wl_test_R)))
plt.title("wavelength against R")
plt.hlines(plain_glass_R,np.min(wl_test), np.max(wl_test), label="glass with no layer")
plt.legend()
plt.grid()
plt.xlabel("wavelength [nm]")
plt.ylabel("Reflectance R")
#depth of layer
plt.figure()
plt.plot(d_test, d_test_R, label="optimal depth = "+str(d_minimise(d_test,d_test_R)))
plt.hlines(plain_glass_R,np.min(d_test), np.max(d_test), label="glass with no layer")
plt.legend()
plt.grid()
plt.xlabel("layer depth [nm]")
plt.ylabel("Reflectance R")
plt.title("depth of layer index against R")
#incident angle
plt.figure()
plt.plot(theta_test, theta_test_R_p, label="350nm brewster angle = "+str(theta_minimise(theta_test,theta_test_R_p))+" rad")
plt.plot(theta_test, theta_test_R_p1, label="450 nm brewster angle = "+str(theta_minimise(theta_test,theta_test_R_p1))+" rad")
plt.plot(theta_test, theta_test_R_p2, label="600nm brewster angle = "+str(theta_minimise(theta_test,theta_test_R_p2))+" rad")
plt.plot(theta_test, theta_test_R_p3, label="700 nm brewster angle = "+str(theta_minimise(theta_test,theta_test_R_p3))+" rad")
plt.plot(theta_test, theta_test_R_p4, label="800 nm brewster angle = "+str(theta_minimise(theta_test,theta_test_R_p4))+" rad")
plt.xticks(ticks=np.arange(0, 5*np.pi/8, np.pi/8), labels=["0", "$\pi$/8", "$\pi$/4", "3$\pi$/8", "$\pi$/2"])
plt.legend()
plt.grid()
plt.xlabel("angle of incidence [rad]")
plt.ylabel("Reflectance R")
plt.title("angle of incidence against R, p polarisation")
plt.figure()
plt.plot(theta_test, theta_test_R_s,label="350 nm")
plt.plot(theta_test, theta_test_R_s1,label="450 nm")
plt.plot(theta_test, theta_test_R_s2,label="600 nm")
plt.plot(theta_test, theta_test_R_s3,label="700 nm")
plt.plot(theta_test, theta_test_R_s4,label="800 nm")
plt.xticks(ticks=np.arange(0, 5*np.pi/8, np.pi/8), labels=["0", "$\pi$/8", "$\pi$/4", "3$\pi$/8", "$\pi$/2"])
plt.legend()
plt.grid()
plt.xlabel("angle of incidence [rad]")
plt.ylabel("Reflectance R")
plt.title("angle of incidence against R, s polarization")
#refractive index
plt.figure()
plt.plot(n_test, n_test_R,label="optimal n = "+str(n_minimise(n_test,n_test_R)))
plt.hlines(plain_glass_R,0,2.5, label="glass with no layer")
plt.legend()
plt.grid()
plt.xlabel("refractive index")
plt.ylabel("Reflectance R")
