#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 20:11:03 2026

@author: morgan
"""


import numpy as np
import matplotlib.pyplot as plt
#from tqdm import tqdm
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from tqdm import tqdm 
import os

def rdm1(psi, n_max):
    C = psi.reshape(n_max, n_max)      # C[i, j]
    rho1 = C @ C.conj().T              # rho1[i, i'] = sum_j C[i,j] C*[i',j]
    return rho1

def rdm2(psi, n_max):
    C = psi.reshape(n_max, n_max)
    rho2 = C.conj().T @ C              # rho2[j, j'] = sum_i C*[i,j] C[i,j']
    return rho2

def purity(rho):
    return np.real(np.trace(rho @ rho))



def entropie(rho):
    S=0
    for i in range(len(rho)):
        S=S+rho[i,i]*np.log(rho[i,i])
    return -S


def propagation(H,psi0,dt):
    A1=(np.identity(len(H))+1j*dt/2*H)
    A2=np.dot((np.identity(len(H))-1j*dt/2*H),psi0)
    return np.linalg.solve(A1,A2)



def Energiemoyenne(H,rho):
    return np.trace(np.dot(H,rho))


def H_b(nmax,L):
    H=np.zeros((nmax,nmax))
    for i in range(nmax):
        H[i,i]=(i+1)**2*np.pi**2/(2*L**2)
    return H 


Ll = 1 # left box size
Lr = np.pi/np.e# right box size
C1 = -800
n_max =70
dim = n_max**2  # number of basis functions per particle
dt=10
Tmin=0
Tmax=1000
T=np.arange(Tmin,Tmax,dt)


H0=np.zeros((dim,dim))
##terme cinetique
for i in range(n_max):
    for j in range(n_max):
        idx1 = i*n_max + j
        kinetic = (i+1)**2 * np.pi**2 / (2 * Ll**2)+(j+1)**2 * np.pi**2 / (2 * Lr**2)
        H0[idx1, idx1] += kinetic
dossier = "/Users/morgan/Documents/postdoc/thermodynamique_quantique/papier_long"
fichier_V = os.path.join(dossier, "V_nmax{}.dat")
f = fichier_V.format(n_max)
if os.path.getsize(f) != dim * dim * 8:
    raise ValueError(f"{f} : taille incompatible avec nmax = {n_max}")
V = np.memmap(f, dtype=np.float64, mode="r", shape=(dim, dim))
V = np.array(V)   

# V=np.loadtxt("potentiel_nmax"+str(n_max)+".txt")    
H=H0+C1*V



psi=np.zeros(len(H),complex)
Purete=np.zeros((len(T),2))
Ener=np.zeros((len(T),2))
Entro=np.zeros((len(T),2))

# psi[3*(n_max)+29]=1


###simulation9
# psi =np.zeros(dim, dtype=complex)
# psi[48*n_max + 19] = 1.0

###simulation11
psi_0= np.zeros(dim, dtype=complex)
j0=0
psi_0[5*n_max + (22)] = 1.0
psi=psi_0
rho1=rdm1(psi,n_max)
rho2=rdm2(psi,n_max)
HG=H_b(n_max,Ll)
HD=H_b(n_max,Lr)


Purete[0,0]=np.trace(np.dot(rho1,rho1))
Purete[0,1]=np.trace(np.dot(rho2,rho2))

Ener[0,0]= Energiemoyenne(HG,rho1)
Ener[0,1]= Energiemoyenne(HD,rho2)

Entro[0,0]=entropie(rho1)
Entro[0,1]=entropie(rho2)

# import sys
# sys.exit()

# ============================================================
# PREPARATION CRANK-NICOLSON UNE SEULE FOIS
# ============================================================

N = H.shape[0]

Id = np.eye(N, dtype=complex)

# A_plus = Id + 1j * dt/2.0 * H
# A_minus = Id - 1j * dt/2.0 * H

# Factorisation LU UNE SEULE FOIS
# lu, piv = lu_factor(A_plus)

# def propagation_fast(psi):
#     rhs = A_minus @ psi
#     return lu_solve((lu, piv), rhs)

D,P=np.linalg.eig(H)
# D=np.diag(D)

P_inv=P.conj().T
psi_eig = P.conj().T @ psi_0

rhoo1=np.zeros(len(rho1))
rhoo2=np.zeros(len(rho2))
cpt=0
t_moy = 0.2 * T[-1]          #
t_debut = T[-1] - t_moy

for i in tqdm(range(len(T)-1)):

    psi = P @ (np.exp(-1j * D * T[i]) * psi_eig)

    # psi = propagation_fast(psi)

    rho1 = rdm1(psi, n_max)
    rho2 = rdm2(psi, n_max)

    Purete[i+1, 0] = np.real(np.trace(rho1 @ rho1))
    Purete[i+1, 1] = np.real(np.trace(rho2 @ rho2))

    Ener[i+1, 0] = np.real(Energiemoyenne(HG, rho1))
    Ener[i+1, 1] = np.real(Energiemoyenne(HD, rho2))

    Entro[i+1, 0] = np.real(entropie(rho1))
    Entro[i+1, 1] = np.real(entropie(rho2))

    if T[i] >= t_debut:
        rhoo1 = rhoo1 + rho1
        rhoo2 = rhoo2 + rho2
        cpt = cpt + 1

rhoo1 = rhoo1 / cpt
rhoo2 = rhoo2 / cpt


A=np.zeros((len(T),6))
A[:,0]=T
A[:,1]=Purete[:,0]
A[:,2]=Ener[:,0]
A[:,3]=Ener[:,1]
A[:,4]=Entro[:,0]
A[:,5]=Entro[:,1]


plt.figure(2)
plt.plot(T,Purete[:,0],color="k",label="purity")   

plt.figure(3)
plt.plot(T,Ener[:,0],'k',label="Energy (left)")
plt.plot(T,Ener[:,1],'b',label="Energy (right)")
plt.xlabel("Time $t$")
plt.ylabel("Energy $E$")

plt.legend()



plt.figure(4)
Energie1 = np.linalg.eigh(HG)[0]   # eigh: HG hermitien -> valeurs propres réelles triées
Energie2 = np.linalg.eigh(HD)[0]
Energie1=np.sort(Energie1)
Energie2=np.sort(Energie2)