#!/usr/bin/env python3

import numpy as np , pylab, tkinter
import math
import matplotlib.pyplot as plt
import mpmath as mp

mub=5.78901e-5          # Bohr magneton (eV/T)
hbar=0.658212           # Planck's constant (eV.fs/rad)
g=2.0                   # Lande factor (adim)
gyro=g*mub/hbar         # gyromag ratio (rad/fs/T)
alpha=0.01              # damping coefficient
pi=math.pi

Bnrm=10.0               # mag. field (T)
Bext = np.array([0.0, 0.0, 1.0])
Sn = 2.0                # spin norm (in # of muB)
S = np.array([1.0, 0.0, 0.0])

N=500000                # number of timesteps
dt=0.1                  # timestep (fs)

# Rodrigues rotation formula
def rotation_matrix(axis, theta):
  """
  Return the rotation matrix associated with counterclockwise
  rotation about the given axis by theta radians
  """
  axis = np.asarray(axis)
  a = math.cos(theta / 2.0)
  b, c, d = -axis * math.sin(theta / 2.0)
  aa, bb, cc, dd = a * a, b * b, c * c, d * d
  bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
  return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
      [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
      [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

# calc. precession field
def calc_rot_vector(Fi,Sp):
  rot = gyro*Sn*Bnrm*(Fi-alpha*np.cross(Fi,Sp))
  return rot
  
# np.set_printoptions(precision=4)
for t in range (0,N):
  wf = calc_rot_vector(Bext,S)
  theta=dt*np.linalg.norm(wf)
  axis=wf/np.linalg.norm(wf)
  S = np.dot(rotation_matrix(axis, theta), S)
  en = -hbar*gyro*Sn*Bnrm*np.dot(S,Bext)
  # print res. in ps for comparison with LAMMPS
  print(t*dt/1000.0,S[0],S[1],S[2],en)

