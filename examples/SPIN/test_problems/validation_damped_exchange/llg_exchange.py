#!/usr/bin/env python3

import numpy as np , pylab, tkinter
import math
import matplotlib.pyplot as plt
import mpmath as mp

hbar=0.658212           # Planck's constant (eV.fs/rad)
J0=0.05                 # per-neighbor exchange interaction (eV)
S1 = np.array([1.0, 0.0, 0.0])
S2 = np.array([0.0, 1.0, 0.0])
alpha=0.01              # damping coefficient
pi=math.pi

N=30000                 # number of timesteps
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

# calculating precession field of spin Sr
def calc_rot_vector(Sr,Sf):
  rot = (J0/hbar)*(Sf-alpha*np.cross(Sf,Sr))/(1.0+alpha**2)
  return rot
  
# second-order ST decomposition as implemented in LAMMPS 
for t in range (0,N):
  # advance s1 by dt/4
  wf1 = calc_rot_vector(S1,S2)
  theta=dt*np.linalg.norm(wf1)*0.25
  axis=wf1/np.linalg.norm(wf1)
  S1 = np.dot(rotation_matrix(axis, theta), S1)
  # advance s2 by dt/2
  wf2 = calc_rot_vector(S2,S1)
  theta=dt*np.linalg.norm(wf2)*0.5
  axis=wf2/np.linalg.norm(wf2)
  S2 = np.dot(rotation_matrix(axis, theta), S2)
  # advance s1 by dt/2
  wf1 = calc_rot_vector(S1,S2)
  theta=dt*np.linalg.norm(wf1)*0.5
  axis=wf1/np.linalg.norm(wf1)
  S1 = np.dot(rotation_matrix(axis, theta), S1)
  # advance s2 by dt/2
  wf2 = calc_rot_vector(S2,S1)
  theta=dt*np.linalg.norm(wf2)*0.5
  axis=wf2/np.linalg.norm(wf2)
  S2 = np.dot(rotation_matrix(axis, theta), S2)
  # advance s1 by dt/4
  wf1 = calc_rot_vector(S1,S2)
  theta=dt*np.linalg.norm(wf1)*0.25
  axis=wf1/np.linalg.norm(wf1)
  S1 = np.dot(rotation_matrix(axis, theta), S1)
  # calc. average magnetization
  Sm = (S1+S2)*0.5
  # calc. energy
  # en = -hbar*(np.dot(S1,wf1)+np.dot(S2,wf2))
  en = -2.0*J0*(np.dot(S1,S2))
  # print res. in ps for comparison with LAMMPS
  print(t*dt/1000.0,Sm[0],Sm[1],Sm[2],en)
  # print(t*dt/1000.0,S1[0],S2[0],S1[1],S2[1],S1[2],S2[2],en)
