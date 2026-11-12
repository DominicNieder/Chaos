#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import os


# ============================================================
# PARAMETERS
# ============================================================

Ll = 1.0
Lr = np.pi / np.e
C1 = 1.0

n_max = 100
Nx = 20 * n_max
Ny = 20 * n_max

eps = 10**-15

# Nombre de paires (i,k) calculées simultanément
block_size = 50


# ============================================================
# GRID
# ============================================================

x_l = np.linspace(-Ll, 0.0, Nx)
x_r = np.linspace(0.0, Lr, Ny)

dx = x_l[1] - x_l[0]
dy = x_r[1] - x_r[0]


# ============================================================
# BASIS FUNCTIONS
# ============================================================

n = np.arange(1, n_max + 1)[:, None]

psiL = np.sqrt(2.0 / Ll) * np.sin(
    n * np.pi / Ll * (x_l[None, :] + Ll)
)

psiR = np.sqrt(2.0 / Lr) * np.sin(
    n * np.pi / Lr * x_r[None, :]
)


# ============================================================
# COULOMB POTENTIAL
# ============================================================

Vpot = C1 / np.sqrt(
    (x_r[None, :] - x_l[:, None])**2 + eps
)

# quadrature directement incorporée
Vpot *= dx * dy


# ============================================================
# UNIQUE INDEX PAIRS
#
# only:
#
#       i <= k
#       j <= l
# ============================================================

pairs = [
    (i, k)
    for i in range(n_max)
    for k in range(i, n_max)
]

Npairs = len(pairs)

print("n_max =", n_max)
print("Nombre de paires complètes :", n_max**2)
print("Nombre de paires indépendantes :", Npairs)


# ============================================================
# PRODUCTS OF BASIS FUNCTIONS
# ============================================================

print("Construction des produits symétriques...")

pairL = np.empty((Npairs, Nx), dtype=np.float64)
pairR = np.empty((Npairs, Ny), dtype=np.float64)

for p, (i, k) in enumerate(pairs):

    pairL[p] = psiL[i] * psiL[k]
    pairR[p] = psiR[i] * psiR[k]


# ============================================================
# COMPUTE INDEPENDENT MATRIX ELEMENTS
#
# W[p,q] corresponds to
#
# p = (i,k)
# q = (j,l)
#
# W[p,q] = V_{ij,kl}
# ============================================================

print("Calcul des éléments indépendants du potentiel...")

W = np.empty(
    (Npairs, Npairs),
    dtype=np.float64
)


for start in tqdm(range(0, Npairs, block_size)):

    stop = min(start + block_size, Npairs)

    # intégration en x
    #
    # (block,Nx) @ (Nx,Ny)
    #
    # -> (block,Ny)

    tmp = pairL[start:stop] @ Vpot

    # intégration en y
    #
    # (block,Ny) @ (Ny,Npairs)
    #
    # -> (block,Npairs)

    W[start:stop] = tmp @ pairR.T


# ============================================================
# FULL V MATRIX
# ============================================================

dim = n_max**2

filename = f"V_nmax{n_max}.dat"

print()
print("Dimension de V :", dim, "x", dim)

size_GB = dim**2 * 8 / 1024**3

print(f"Taille théorique : {size_GB:.2f} Go")
print("Fichier :", filename)


# ============================================================
# MEMMAP
#
# V est réellement une matrice dim x dim,
# mais stockée directement sur disque.
#
# Elle se comporte presque comme un ndarray.
# ============================================================

V = np.memmap(
    filename,
    dtype=np.float64,
    mode="w+",
    shape=(dim, dim)
)


# ============================================================
# RECONSTRUCTION USING SYMMETRIES
# ============================================================

print("Reconstruction de la matrice complète...")


for p, (i, k) in enumerate(tqdm(pairs)):

    for q, (j, l) in enumerate(pairs):

        value = W[p, q]

        # ----------------------------------------------------
        # Base product index:
        #
        # |i,j> -> i*n + j
        # |k,l> -> k*n + l
        #
        # Main element:
        #
        # V_{ij,kl}
        # ----------------------------------------------------

        a = i * n_max + j
        b = k * n_max + l

        V[a, b] = value
        V[b, a] = value                # Hermiticity


        # ----------------------------------------------------
        # symmetry i <-> k
        #
        # V_{kj,il}
        # ----------------------------------------------------

        a2 = k * n_max + j
        b2 = i * n_max + l

        V[a2, b2] = value
        V[b2, a2] = value


        # ----------------------------------------------------
        # symmetry j <-> l
        #
        # V_{il,kj}
        # ----------------------------------------------------

        a3 = i * n_max + l
        b3 = k * n_max + j

        V[a3, b3] = value
        V[b3, a3] = value


        # ----------------------------------------------------
        # both symmetries
        #
        # V_{kl,ij}
        #
        # already related by Hermiticity, but explicit
        # assignment causes no problem.
        # ----------------------------------------------------

        a4 = k * n_max + l
        b4 = i * n_max + j

        V[a4, b4] = value
        V[b4, a4] = value


# force writing to disk
V.flush()


print()
print("Matrice V construite.")


# # ============================================================
# # KINETIC ENERGY
# # ============================================================

# E_L = (
#     np.arange(1, n_max + 1, dtype=float)**2
#     * np.pi**2
#     / (2.0 * Ll**2)
# )

# E_R = (
#     np.arange(1, n_max + 1, dtype=float)**2
#     * np.pi**2
#     / (2.0 * Lr**2)
# )

# Ekin = (
#     E_L[:, None]
#     +
#     E_R[None, :]
# ).ravel()


# ============================================================
# OPTIONAL: ADD KINETIC ENERGY DIRECTLY TO V

# If you want H = T + V
# ============================================================

H_filename = f"H_nmax{n_max}.dat"

H = np.memmap(
    H_filename,
    dtype=np.float64,
    mode="w+",
    shape=(dim, dim)
)


# Copy potential
#
# better row by row to avoid huge RAM allocation

print("Construction de H = T + V...")

copy_block = 1000

for start in tqdm(range(0, dim, copy_block)):

    stop = min(start + copy_block, dim)

    H[start:stop] = V[start:stop]


# kinetic part is diagonal

# diag = np.arange(dim)

# H[diag, diag] += Ekin

# H.flush()


print()
print("Hamiltonien complet construit.")
print("Fichier :", H_filename)

# n_max = 100
# dim = n_max**2

# H3 = np.memmap(
#     "H_nmax50.dat",
#     dtype=np.float64,
#     mode="r",              # lecture seule
#     shape=(dim, dim)
# )