import numpy as np
import math
import random
import h5py
import os
import os, sys
from pylab import *
from math import acos, sin, cos, sqrt, pi, exp, log, floor
from abc import ABCMeta, abstractmethod
from inspect import signature

from astropy import constants as const
from astropy import units as u
import astropy.units as u
from astropy.utils import isiterable
from astropy.cosmology import FlatLambdaCDM

from scipy import optimize as op
from scipy import signal as ss
from scipy import interpolate
from time import time
import random

M_min = 10.5
M_max = 16.5

d2 = np.loadtxt("/home2/guillermo/Halos/prueba/Halos_tree_DOC_PID_Vmax_all_Mass.txt")# Leemos el fichero de halos
f_12 = d2[:,0]
f2 = int(len(f_12))
c_12 = d2[0,:]
c = int(len(c_12))
mvir_l2 = []


#Formato h5
filename = "/home/aknebe/Projects/UNITSIM/SAMs/SAGE/ELGs/UNITSIM1/UNITSIM1_model_z1.321_ELGs.h5" # Leemos el fichero de galaxias
f = h5py.File(filename, 'r')
c = list(f.keys())
data0 = np.array(f['HostHaloID'])
data1 = np.array(f['MainHaloID'])
data2 = np.array(f['MainMhalo'])
data3 = np.array(f['Mhalo'])
data4 = np.array(f['Mstar'])
data5 = np.array(f['Xpos'])
data6 = np.array(f['Xvel'])
data7 = np.array(f['Ypos'])
data8 = np.array(f['Yvel'])
data9 = np.array(f['Zpos'])
data10 = np.array(f['Zvel'])
data11 = np.array(f['logFHalpha'])
data12 = np.array(f['logFHalpha_att'])
datos = int(len(data0))
print(datos)

sel = np.where((data12 > log10(1.32502*10**(-16))))[0] # Cutting of objects that meet flux condition for HOD pessimistic model
#sel = np.where((data12 > log10(2*10**(-16))))[0]
#sel = np.where((data12 > log10(1.1645*10**(-16))))[0]

data0F = data0[sel]
data1F = data1[sel]
data5F = data5[sel]
data6F = data6[sel]
data7F = data7[sel]
data8F = data8[sel]
data9F = data9[sel]
data10F = data10[sel]
data2F = log10(data2[sel])

sel3 = np.where((data0F != data1F) )[0] # We search for satellite galaxies that meet the flow condition
SM = data1F[sel3] # MainHalo array of all satellite that meet the flux condition
SM_0 = data0F[sel3]
SM_x = data5F[sel3]
SM_y = data7F[sel3]
SM_z = data9F[sel3]
SM_vx = data6F[sel3]
SM_vy = data8F[sel3]
SM_vz = data10F[sel3]
SM_M = data2F[sel3]

sel1 = np.where((data0F == data1F) )[0] # We search all the central galaxies
C = data0F[sel1] # HostHalo array of all central galaxies
C1 = data1F[sel1]
C_x = data5F[sel1]
C_y = data7F[sel1]
C_z = data9F[sel1]
C_vx = data6F[sel1]
C_vy = data8F[sel1]
C_vz = data10F[sel1]
C_M = data2F[sel1]
nsm = len(SM)

Nsat = float(len(SM_M))
print("Numero de galaxias Satellites segun SAGE")
print((Nsat))

Ncent = float(len(C_M))
print("Numero de galaxias Centrales segun SAGE")
print((Ncent))

mvir_l_central = []
mvir_l_sate = []
mvir_Tot_l = []


selwc = np.isin( SM , C  )
Gswc = SM[selwc]
Nswc = float(len(Gswc))
Mswc = SM_M[selwc]

selwoc = np.isin( SM , C , invert = True )
Gswoc = SM[selwoc]
Nswoc = float(len(Gswoc))
Mswoc = SM_M[selwoc]

print("Numero de Galaxias satellites con central segun SAGE")
print(Nswc)

print("Numero de Galaxias satellites sin central segun SAGE")
print(Nswoc)


Masas_Halos = []
for i in range( 0 , f2  , 1):# Creamos un array de Halos en el que cada Halo se correspode con un Halo Host que consta de las corresponientes que queremos guardar: ID , Rvir , Rs , Rsk , x , y , z
      b2 = d2[i,:]
      ID_Halo = int(b2[1])
      #mvir2 = float(b2[9])
      Rvir = float(b2[10])
      x_h = float(b2[3])
      y_h = float(b2[4])
      z_h = float(b2[5])
      Rs = float(b2[11])
      Rsk = float(b2[14])
      #mvirr2 = math.log10(float(mvir2))
      PID = float(b2[13])
      if PID == -1.:
        mvir2 = float(b2[17])
        mvirr2 = log10((mvir2))
        #c = np.array([ID_Halo,mvirr2,Rvir])
        c = np.array([ID_Halo,mvirr2,Rvir,x_h,y_h,z_h,Rs,Rsk])
        mvir_l2.append(c)
        Masas_Halos.append(mvirr2)

m_Halos = np.array ( mvir_l2 )
print("numero de Halos totales")
N_halos = float(len(m_Halos))
print( N_halos )
mvirr22 = np.array (Masas_Halos)


N_bin = 70.
l_bin = (M_max - M_min)/(N_bin)# Definimos una longitud de los bines en los que vamos a calcular el HOD

#Mass = np.arange(10.3, 16.5 + l_bin , l_bin) #Generamos el array que definen los correspondientes bienes de masa a estudiar
#print("Bines de Masa")
#print(Mass)
#print("longitud del bin")
#print(l_bin)
#Mass_Bin = Mass[19:21]
#print("Bin de Masa elegido")
#print(Mass_Bin)

#for i in range( 0 , 1 , N_bin):

N_C = []
N_S = []
N_swc = []
N_swoc = []
N_Halos = []

N_cm_v = []
N_sm_v = []
N_swcm_v = []
N_swocm_v = []

M_minimos = []
M_maximos = []

K1 = []
K2 = []

ii = M_min
while ii < M_max:
    
    minf = ii    
    msup = ii + l_bin

    soldc_Bin = np.where((minf <=  C_M)  &  ( C_M < msup))[0]
    M_C_Bin = C_M[soldc_Bin]
    N_C_Bin = float(len(M_C_Bin))
    #print("Numero de galaxias Centrales segun SAGE dentro del bin de masas")
    #print(N_C_Bin)

    solds_Bin = np.where((minf <= SM_M )  &  ( SM_M  < msup))[0]
    M_S_Bin = SM_M [solds_Bin]
    N_S_Bin = float(len(M_S_Bin))
    #print("Numero de galaxias Satellites segun SAGE dentro del bin de masas")
    #print(N_S_Bin)

    soldswc_Bin = np.where((minf <= Mswc )  &  ( Mswc < msup))[0]
    M_Swc_Bin = Mswc[soldswc_Bin]
    N_Swc_Bin = float(len(M_Swc_Bin))
    #print("Numero de galaxias satellites con central dentro del bin de masas")
    #print(N_Swc_Bin)

    soldswoc_Bin = np.where((minf <= Mswoc )  &  ( Mswoc < msup))[0]
    M_Swoc_Bin = Mswoc[soldswoc_Bin]
    N_Swoc_Bin = float(len(M_Swoc_Bin))
    #print("Numero de galaxias satellites sin central dentro del bin de masas")
    #print(N_Swoc_Bin)
    
    soldHalos_Bin = np.where((minf <= mvirr22 ) &  ( mvirr22 < msup))[0]
    M_Halos_Bin = mvirr22[soldHalos_Bin]
    N_Halos_Bin = float(len(M_Halos_Bin))
    #N_cm = N_C_Bin/N_Halos_Bin
    #N_sm = N_S_Bin/N_Halos_Bin
    #N_swcm = N_Swc_Bin/N_Halos_Bin
    #N_swocm = N_Swoc_Bin/N_Halos_Bin
    
    ii += l_bin
    if N_C_Bin != 0.0 or N_S_Bin != 0.0:
        N_cm = N_C_Bin/N_Halos_Bin
        N_sm = N_S_Bin/N_Halos_Bin
        N_swcm = N_Swc_Bin/N_Halos_Bin
        N_swocm = N_Swoc_Bin/N_Halos_Bin     
        
        M_minimos.append(float(minf))
        M_maximos.append(float(msup))


        N_C.append(int(N_C_Bin))        
        N_S.append(int(N_S_Bin))
        N_swc.append(int(N_Swc_Bin))
        N_swoc.append(int(N_Swoc_Bin))
        N_Halos.append(int(N_Halos_Bin))
        N_cm_v.append(float(N_cm))
        N_sm_v.append(float(N_sm))
        N_swcm_v.append(float(N_swcm))
        N_swocm_v.append(float(N_swocm))

        if N_Swc_Bin != 0.0 and N_Swoc_Bin != 0.0:
            k1 = (N_Swc_Bin)*(1./N_C_Bin)*(1./N_S_Bin)*N_Halos_Bin
            k2 = (1 - k1 * (N_C_Bin/N_Halos_Bin))/(1 - (N_C_Bin/N_Halos_Bin))
            #k2 = 1 - k1
            K1.append(float(k1))
            K2.append(float(k2))

        elif N_Swc_Bin == 0.0 and N_Swoc_Bin != 0.0:
            k1 = 0.0
            k2 = (N_Swoc_Bin)*(1./(N_Halos_Bin-N_C_Bin))*(1./N_S_Bin)*N_Halos_Bin
            K1.append(k1)
            K2.append(k2)
        
        elif N_Swc_Bin != 0.0 and N_Swoc_Bin == 0.0:
            k1 = (N_Swc_Bin)*(1./N_C_Bin)*(1./N_S_Bin)*N_Halos_Bin
            k2 = 0.0
            K1.append(k1)
            K2.append(k2)
        
        else:
            k1 = 0.0
            k2 = 0.0
            K1.append(k1)
            K2.append(k2)

"""    
Halos_Bin = []
for i in range( 0, int(len(m_Halos)) , 1 ):# Halos
    Lh_i = m_Halos[i]
    M_hi = Lh_i[1]
    if  Mass_Bin[0] <= M_hi < Mass_Bin[1]:
    	Halos_Bin.append((Lh_i))

    
N_halos_Bin = float(len(Halos_Bin))
print("Numero de Halos dentro del bin de masas")
print(N_halos_Bin)
"""
print("Comparación")
print(sum(N_C))
print(sum(N_S))
print(sum(N_swc))
print(sum(N_swoc))
print(sum(N_Halos))

N_C_global = float((sum(N_C))/(sum(N_Halos)))
N_S_global = float((sum(N_S))/(sum(N_Halos)))
N_swc_global = float((sum(N_swc))/(sum(N_Halos)))
N_swoc_global = float((sum(N_swoc))/(sum(N_Halos)))

k1_global_1 = sum(N_swc)/(sum(np.array(N_S) * np.array( N_cm_v )))
k2_global_1 = sum(N_swoc)/(sum(np.array(N_S) * (1 - np.array(N_cm_v))))

k1_global = (sum(N_swc)*sum(N_Halos))/(sum(N_S)*sum(N_C))
k2_global = (1 - k1_global*(sum(N_C)/sum(N_Halos)))/(1 - (sum(N_cm_v)/sum(N_Halos)))


"""
print("Probabilidad de Centrales")
N_C_global = float((sum(N_C))/((N_halos)))
print(N_C_global)
N_S_global = float((sum(N_S))/((N_halos)))
N_swc_global = float((sum(N_swc))/((N_halos)))
N_swoc_global = float((sum(N_swoc))/((N_halos)))
"""
print("Parametros")
print(k1_global)
print(k2_global)
print(k1_global_1)
print(k2_global_1)
"""
k1_global = (N_swc_global)*(1./N_C_global)*(1./N_S_global)
k2_global = (1 -k1_global * N_C_global)/(1 - N_C_global)
print(k1_global)
print(k2_global)

k1_global_1 = (N_swc_global)*(1./N_S_global)
k2_global_2 = 1 - k1_global_1
print(k1_global_1)
print(k2_global_2)
"""

outfile = open( "Fichero_prueba_Correct_70_bin_1.txt", 'w')


for i in range(0,int(len(M_minimos)),1):
        ii = int(i)
        o = M_minimos[i]
        t = M_maximos[i]
        r = K1[i]
        u = K2[i]
        v = N_cm_v[i]
        l = N_sm_v[i]
        s = N_swcm_v[i]
        g = N_swocm_v[i]
        a = N_C[i]
        b = N_S[i]
        c = N_swc[i]
        d = N_swoc[i]
        e = N_Halos[i]
        
        if ii == 0 :
            outfile.write('# Mmim Mmax k1 k2 Vm_Cen Vm_Sat Vm_swc Vm_swoc N_Cen N_Sat N_swc N_swoc N_Halos K1_global K2_global K1_global_1 K2_global_1 \n')
            outfile.write(' %f %f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %d %d %d %d %d %10.4e %10.4e %10.4e %10.4e \n' %(o,t,r,u,v,l,s,g,a,b,c,d,e,k1_global,k2_global,k1_global_1,k2_global_1))
        
        else:
            outfile.write(' %f %f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %d %d %d %d %d \n' %(o,t,r,u,v,l,s,g,a,b,c,d,e))
outfile.close()







"""
out = np.transpose(np.array([ N_halos , Ncent , Nsat , Nswc , Nswoc , N_halos_Bin , N_C_Bin , N_S_Bin , N_Swc_Bin , N_Swoc_Bin , N_cm , N_sm  , N_swcm , N_swocm , k1 , k2]))
np.savetxt("Galaxy_Bin_HoD_Model_Mass200c_40_K_const_Comprobacion.dat", out , header = "N_halos_Tot N_cen_Tot N_Sat_Tot N_Swc_Tot N_Swoc_Tot N_Halos_Bin N_Cen_Bin N_Sat_Bin N_Swc_Bin N_Swoc_Bin Vm_Cen Vm_Sat Vm_Swc Vm_Swoc K1 K2 ", fmt = "%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e ")
"""
"""
HRvir = []
HMass = []
HID = []
Hx = []
Hy = []
Hz = []
HRs = []
HRsk = []

for e in range( 0 , len(Halos_Bin) , 1):
    a = Halos_Bin[e]
    Rvh = a[2]
    Masah = a[1]
    IDh = a[0]
    xxh = a[3]
    yyh = a[4]
    zzh = a[5]
    Rsh = a[6]
    Rskh = a[7]
    HRvir.append(Rvh)
    HMass.append(Masah)
    HID.append(IDh)
    Hx.append(xxh)
    Hy.append(yyh)
    Hz.append(zzh)
    HRs.append(Rsh)
    HRsk.append(Rskh)


outfile = open( "HoD_s_UNITSIM_Model_M_halo_main_3_Rs_PDF_PID_Mod_N_M200c_40_Bin_Comprovacion.txt", 'w')

for i in range(0,int(len(HRvir)),1):
        
        o = HRvir[i]
        t = HMass[i]
        r = HID[i]
        u = Hx[i]
        v = Hy[i]
        l = Hz[i]
        s = HRs[i]
        g = HRsk[i]

        outfile.write(' %f %f %f %f %f %f %f %f \n' %(float(r),float(t),float(o),float(u),float(v),float(l),float(s),float(g)))
outfile.close()
"""
