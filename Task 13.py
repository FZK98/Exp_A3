# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:07:24 2021

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import minimize

#load necessary data
wl_BK7, n_BK7, k_BK7=np.loadtxt('BK7.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_MgF2, n_MgF2, k_MgF2=np.loadtxt('MgF2.txt', delimiter = '\t', skiprows=1, unpack=True ) #BK7 wavelength in nm
wl_gold, n_gold, k_gold = np.loadtxt('Au.txt', delimiter = '\t', skiprows=1, unpack=True) #gold wavelength in nm
#wl_Ta2O5, n_Ta2O5, k_Ta2O5 = np.loadtxt(open("Ta2O5.csv"), delimiter=",", skiprows=1, unpack=True)
#wl_Ta2O5*=1000
wl_Ta2O5, n_Ta2O5, k_Ta2O5 = np.loadtxt(open("Ta2O5_2.csv"), delimiter=",", skiprows=1, unpack=True)
wl_Ta2O5*=1000

def transfer_matrix(theta_i, user_wl, polarization,materials,d):
    #theta_i = 0
    #user_wl = 700 #nm
   # polarization = "s" #or "p"
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
        elif materials[i] == "Ta2O5":
            refractive_index_n.append(get_refractiveIndex_n(user_wl, wl_Ta2O5, n_Ta2O5))
            refractive_index_k.append(get_refractiveIndex_k(user_wl, wl_Ta2O5, k_Ta2O5))
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

 #interpolate refractive index for analysis
	
def interpolate_n(wl, wl_data, n_data):
    wl_integer = np.arange(np.rint(np.min(wl_data)),np.rint(np.max(wl_data)),1) #every integer wavelength from 330-2500
    n_integer = np.interp(wl_integer, wl_data, n_data) #refractive index for each integer wavelength
    loc=np.where(wl==wl_integer)
    return n_integer[loc]

# =============================================================================
# N against R for given wavelength
## =============================================================================
number_test=np.arange(0,50,1)
number_test_R=[]
for i in number_test:
    mg_d_test = 633/(4*float(interpolate_n(633, wl_MgF2, n_MgF2)))
    ta2o5_d_test = 633/(4*float(interpolate_n(633, wl_Ta2O5, n_Ta2O5)))
    materials_chosen=["Ta2O5","MgF2"]*i
    depths_chosen=[ta2o5_d_test,mg_d_test]*i
    materials_N=["air"]
    depths_N=[0]
    for j in range(2*i):
        materials_N.append(materials_chosen[j])
        depths_N.append(depths_chosen[j])
    materials_N.append("BK7")    
    depths_N.append(0)
    #print(materials_N)
    #print(depths_N)
    number_test_R_temp, number_test_T_temp = transfer_matrix(0,633,"s",materials_N,depths_N)
    number_test_R.append(abs(number_test_R_temp)**2)
plt.figure()
plt.plot(number_test, number_test_R)    
plt.grid()
plt.xlabel("number of periods N")
plt.ylabel("R")


ninetynine=np.round(np.interp(0.9999, number_test_R, number_test))
print("need "+str(ninetynine)+" layers to get 99.99% reflectance")

# =============================================================================
# Varying incident Angle against R for a for a given N
## =============================================================================

#make a dielectric stack with N periods at a given wavelength and angle 
wl_bragg = 633
def stack(theta, wl, polarisation, N):
    mg_d = wl_bragg/(4*float(interpolate_n(wl, wl_MgF2, n_MgF2)))
    ta2o5_d= wl_bragg/(4*float(interpolate_n(wl, wl_Ta2O5, n_Ta2O5)))
    materials_chosen=["Ta2O5","MgF2"]*N
    depths_chosen=[ta2o5_d,mg_d]*N
    materials_N=["air"]
    depths_N=[0]
    for j in range(2*N):
        materials_N.append(materials_chosen[j])
        depths_N.append(depths_chosen[j])
    materials_N.append("BK7")    
    depths_N.append(0)
    number_test_R_temp, number_test_T_temp = transfer_matrix(theta,wl, polarisation ,materials_N,depths_N)
    number_test_R = (abs(number_test_R_temp)**2)
    return(number_test_R)

incident_angles = np.linspace(0 , np.pi/2 , 90)
angle_test_R_S = []
angle_test_R_P = []
for i in incident_angles :
    angle_test_R_S.append(stack(i, 633,"s", 12))
    angle_test_R_P.append(stack(i, 633,"p", 12))
    
    


plt.figure()
plt.plot(incident_angles, angle_test_R_S)    
plt.grid()
plt.xlabel("Angle of incidence")
plt.ylabel("R")
plt.title("Angle of incidence against Reflectance (s polarised)")

plt.figure()
plt.plot(incident_angles, angle_test_R_P)    
plt.grid()
plt.xlabel("Angle of incidence")
plt.ylabel("R")
plt.title("Angle of incidence against Reflectance (p polarised)")

# =============================================================================
# R against wavelength for different angles
# =============================================================================
wl_bragg =633
def stack_fixed_d(theta, wl, polarisation, N):
    mg_d = wl_bragg/(4*float(interpolate_n(wl, wl_MgF2, n_MgF2)))
    ta2o5_d= wl_bragg/(4*float(interpolate_n(wl, wl_Ta2O5, n_Ta2O5)))
    materials_chosen=["Ta2O5","MgF2"]*N
    depths_chosen=[ta2o5_d,mg_d]*N
    materials_N=["air"]
    depths_N=[0]
    for j in range(2*N):
        materials_N.append(materials_chosen[j])
        depths_N.append(depths_chosen[j])
    materials_N.append("BK7")    
    depths_N.append(0)
    number_test_R_temp, number_test_T_temp = transfer_matrix(theta,wl, polarisation ,materials_N,depths_N)
    number_test_R = (abs(number_test_R_temp)**2)
    return(number_test_R)

wl_test = np.arange(400,900,1)
wl_test_R1s =[]
wl_test_R2s =[]
wl_test_R3s =[]
wl_test_R4s =[]
wl_test_R5s = []

wl_test_R1p =[]
wl_test_R2p =[]
wl_test_R3p =[]
wl_test_R4p =[]
wl_test_R5p =[]

number_periods=8

for i in wl_test:
	wl_test_R1s.append(stack_fixed_d(0,i,"s",number_periods))
	wl_test_R2s.append(stack_fixed_d(np.pi/16,i,"s",number_periods))
	wl_test_R3s.append(stack_fixed_d(np.pi/8,i,"s",number_periods))
	wl_test_R4s.append(stack_fixed_d(3*np.pi/16,i,"s",number_periods))
	wl_test_R5s.append(stack_fixed_d(np.pi/4,i,"s",number_periods))
	
	
	wl_test_R1p.append(stack_fixed_d(0,i,"p",number_periods))
	wl_test_R2p.append(stack_fixed_d(np.pi/16,i,"p",number_periods))
	wl_test_R3p.append(stack_fixed_d(np.pi/8,i,"p",number_periods))
	wl_test_R4p.append(stack_fixed_d(3*np.pi/16,i,"p",number_periods))
	wl_test_R5p.append(stack_fixed_d(np.pi/4, i, "p", number_periods))
	
plt.figure()
plt.plot(wl_test, wl_test_R1s, label="0")
#plt.plot(wl_test, wl_test_R2s, label="$\pi$/16")
plt.plot(wl_test, wl_test_R3s, label="$\pi$/8")
#plt.plot(wl_test, wl_test_R4s, label="3$\pi$/8")
plt.plot(wl_test, wl_test_R5s, label="$\pi$/4")
plt.legend()
plt.grid()
plt.title("s polarization")
plt.xlabel("wavelength(nm)")
plt.ylabel("Reflectance R")

plt.figure()
plt.plot(wl_test, wl_test_R1p, label="0")
#plt.plot(wl_test, wl_test_R2p, label="$\pi$/16")
plt.plot(wl_test, wl_test_R3p, label="$\pi$/8")
#plt.plot(wl_test, wl_test_R4p, label="3$\pi$/8")
plt.plot(wl_test, wl_test_R5p, label="$\pi$/4")
plt.legend()
plt.grid()
plt.title("p polarization")
plt.xlabel("wavelength(nm)")
plt.ylabel("Reflectance R")
# =============================================================================
# contour plots for wavelength and angle against R
# =============================================================================
wl_test=np.arange(400,900,1)
theta_test=np.arange(0,np.pi/2,0.1)
contours = np.zeros((len(theta_test),len(wl_test)))
contourp = np.zeros((len(theta_test),len(wl_test)))
for i in range(len(wl_test)):
	for j in range(len(theta_test)):
		contours[j,i]=stack(theta_test[j],wl_test[i],"s",6)
		contourp[j,i]=stack(theta_test[j],wl_test[i],"p",6)
plt.figure()
plt.contourf(contours,levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.colorbar()
plt.xticks(ticks=[0,100,200,300,400,500],labels=["400","500","600","700","800","900"])
plt.yticks(ticks=[0,3,6,9,12,15],labels=["0","$\pi$/10","$\pi$/5","$3\pi$/10","$2\pi$/5","$\pi$/2"])
plt.xlabel("wavelength(nm)")
plt.ylabel("Angle of incidence")
plt.title("s polarization")

plt.figure()
plt.contourf(contourp,levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
plt.colorbar()
plt.xticks(ticks=[0,100,200,300,400,500],labels=["400","500","600","700","800","900"])
plt.yticks(ticks=[0,3,6,9,12,15],labels=["0","$\pi$/10","$\pi$/5","$3\pi$/10","$2\pi$/5","$\pi$/2"])
plt.xlabel("wavelength(nm)")
plt.ylabel("Angle of incidence")
plt.title("p polarization")
