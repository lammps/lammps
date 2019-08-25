"""Contains simple algorithms.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Functions:
   matrix_exp: Computes the exponential of a square matrix via a Taylor series.
   stab_cholesky: A numerically stable version of the Cholesky decomposition.
   h2abc: Takes the representation of the system box in terms of an upper
      triangular matrix of column vectors, and returns the representation in
      terms of the lattice vector lengths and the angles between them
      in radians.
   h2abc_deg: Takes the representation of the system box in terms of an upper
      triangular matrix of column vectors, and returns the representation in
      terms of the lattice vector lengths and the angles between them in
      degrees.
   abc2h: Takes the representation of the system box in terms of the lattice
      vector lengths and the angles between them, and returns the
      representation in terms of an upper triangular lattice vector matrix.
   invert_ut3x3: Inverts a 3*3 upper triangular matrix.
   det_ut3x3(h): Finds the determinant of a 3*3 upper triangular matrix.
   eigensystem_ut3x3: Finds the eigenvector matrix and eigenvalues of a 3*3
      upper triangular matrix
   exp_ut3x3: Computes the exponential of a 3*3 upper triangular matrix.
   root_herm: Computes the square root of a positive-definite hermitian
      matrix.
   logsumlog: Routine to accumulate the logarithm of a sum
"""

__all__ = ['matrix_exp', 'stab_cholesky', 'h2abc', 'h2abc_deg', 'abc2h',
           'invert_ut3x3', 'det_ut3x3', 'eigensystem_ut3x3', 'exp_ut3x3',
            'root_herm', 'logsumlog' ]

import numpy as np
import math
from ipi.utils.messages import verbosity, warning

def logsumlog(lasa, lbsb):
   """Computes log(|A+B|) and sign(A+B) given log(|A|), log(|B|),
   sign(A), sign(B).

   Args:
      lasa: (log(|A|), sign(A)) as a tuple
      lbsb: (log(|B|), sign(B)) as a tuple

   Returns:
      (log(|A+B|), sign(A+B)) as a tuple
   """

   (la,sa) = lasa
   (lb,sb) = lbsb

   if (la > lb):
      sr = sa
      lr = la + np.log(1.0 + sb*np.exp(lb-la))
   else:
      sr = sb
      lr = lb + np.log(1.0 + sa*np.exp(la-lb))

   return (lr,sr)

def matrix_exp(M, ntaylor=15, nsquare=15):
   """Computes the exponential of a square matrix via a Taylor series.

   Calculates the matrix exponential by first calculating exp(M/(2**nsquare)),
   then squaring the result the appropriate number of times.

   Args:
      M: Matrix to be exponentiated.
      ntaylor: Optional integer giving the number of terms in the Taylor series.
         Defaults to 15.
      nsquare: Optional integer giving how many times the original matrix will
         be halved. Defaults to 15.

   Returns:
      The matrix exponential of M.
   """

   n = M.shape[1]
   tc = np.zeros(ntaylor+1)
   tc[0] = 1.0
   for i in range(ntaylor):
      tc[i+1] = tc[i]/(i+1)

   SM = np.copy(M)/2.0**nsquare

   EM = np.identity(n,float)*tc[ntaylor]
   for i in range(ntaylor-1,-1,-1):
      EM = np.dot(SM,EM)
      EM += np.identity(n)*tc[i]

   for i in range(nsquare):
      EM = np.dot(EM,EM)
   return EM

def stab_cholesky(M):
   """ A numerically stable version of the Cholesky decomposition.

   Used in the GLE implementation. Since many of the matrices used in this
   algorithm have very large and very small numbers in at once, to handle a
   wide range of frequencies, a naive algorithm can end up having to calculate
   the square root of a negative number, which breaks the algorithm. This is
   due to numerical precision errors turning a very tiny positive eigenvalue
   into a tiny negative value.

   Instead of this, an LDU decomposition is used, and any small negative numbers
   in the diagonal D matrix are assumed to be due to numerical precision errors,
   and so are replaced with zero.

   Args:
      M: The matrix to be decomposed.
   """

   n = M.shape[1]
   D = np.zeros(n,float)
   L = np.zeros(M.shape,float)
   for i in range(n):
      L[i,i] = 1.
      for j in range(i):
         L[i,j] = M[i,j]
         for k in range(j):
            L[i,j] -= L[i,k]*L[j,k]*D[k]
         if (not D[j] == 0.0):
            L[i,j] = L[i,j]/D[j]
      D[i] = M[i,i]
      for k in range(i):
         D[i] -= L[i,k]*L[i,k]*D[k]

   S = np.zeros(M.shape,float)
   for i in range(n):
      if (D[i]>0):
         D[i] = math.sqrt(D[i])
      else:
         warning("Zeroing negative element in stab-cholesky decomposition: " + str(D[i]), verbosity.low)
         D[i] = 0
      for j in range(i+1):
         S[i,j] += L[i,j]*D[j]
   return S

def h2abc(h):
   """Returns a description of the cell in terms of the length of the
      lattice vectors and the angles between them in radians.

   Args:
      h: Cell matrix in upper triangular column vector form.

   Returns:
      A list containing the lattice vector lengths and the angles between them.
   """

   a = float(h[0,0])
   b = math.sqrt(h[0,1]**2 + h[1,1]**2)
   c = math.sqrt(h[0,2]**2 + h[1,2]**2 + h[2,2]**2)
   gamma = math.acos(h[0,1]/b)
   beta = math.acos(h[0,2]/c)
   alpha = math.acos(np.dot(h[:,1], h[:,2])/(b*c))

   return a, b, c, alpha, beta, gamma

def h2abc_deg(h):
   """Returns a description of the cell in terms of the length of the
      lattice vectors and the angles between them in degrees.

   Args:
      h: Cell matrix in upper triangular column vector form.

   Returns:
      A list containing the lattice vector lengths and the angles between them
      in degrees.
   """

   (a, b, c, alpha, beta, gamma) = h2abc(h)
   return a, b, c, alpha*180/math.pi, beta*180/math.pi, gamma*180/math.pi

def abc2h(a, b, c, alpha, beta, gamma):
   """Returns a lattice vector matrix given a description in terms of the
   lattice vector lengths and the angles in between.

   Args:
      a: First cell vector length.
      b: Second cell vector length.
      c: Third cell vector length.
      alpha: Angle between sides b and c in radians.
      beta: Angle between sides a and c in radians.
      gamma: Angle between sides b and a in radians.

   Returns:
      An array giving the lattice vector matrix in upper triangular form.
   """

   h = np.zeros((3,3) ,float)
   h[0,0] = a
   h[0,1] = b*math.cos(gamma)
   h[0,2] = c*math.cos(beta)
   h[1,1] = b*math.sin(gamma)
   h[1,2] = (b*c*math.cos(alpha) - h[0,1]*h[0,2])/h[1,1]
   h[2,2] = math.sqrt(c**2 - h[0,2]**2 - h[1,2]**2)
   return h

def invert_ut3x3(h):
   """Inverts a 3*3 upper triangular matrix.

   Args:
      h: An upper triangular 3*3 matrix.

   Returns:
      The inverse matrix of h.
   """

   ih = np.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]*ih[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]*ih[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[0,0] - ih[0,0]*h[0,2]*ih[2,2]
   return ih

def eigensystem_ut3x3(p):
   """Finds the eigenvector matrix of a 3*3 upper-triangular matrix.

   Args:
      p: An upper triangular 3*3 matrix.

   Returns:
      An array giving the 3 eigenvalues of p, and the eigenvector matrix of p.
   """

   eigp = np.zeros((3,3), float)
   eigvals = np.zeros(3, float)

   for i in range(3):
      eigp[i,i] = 1
   eigp[0,1] = -p[0,1]/(p[0,0] - p[1,1])
   eigp[1,2] = -p[1,2]/(p[1,1] - p[2,2])
   eigp[0,2] = -(p[0,1]*p[1,2] - p[0,2]*p[1,1] + p[0,2]*p[2,2])/((p[0,0] - p[2,2])*(p[2,2] - p[1,1]))

   for i in range(3):
      eigvals[i] = p[i,i]
   return eigp, eigvals

def det_ut3x3(h):
   """Calculates the determinant of a 3*3 upper triangular matrix.

   Note that the volume of the system box when the lattice vector matrix is
   expressed as a 3*3 upper triangular matrix is given by the determinant of
   this matrix.

   Args:
      h: An upper triangular 3*3 matrix.

   Returns:
      The determinant of h.
   """

   return h[0,0]*h[1,1]*h[2,2]

MINSERIES=1e-8
def exp_ut3x3(h):
   """Computes the matrix exponential for a 3x3 upper-triangular matrix.

   Note that for 3*3 upper triangular matrices this is the best method, as
   it is stable. This is terms which become unstable as the
   denominator tends to zero are calculated via a Taylor series in this limit.

   Args:
      h: An upper triangular 3*3 matrix.

   Returns:
      The matrix exponential of h.
   """
   eh = np.zeros((3,3), float)
   e00 = math.exp(h[0,0])
   e11 = math.exp(h[1,1])
   e22 = math.exp(h[2,2])
   eh[0,0] = e00
   eh[1,1] = e11
   eh[2,2] = e22

   if (abs((h[0,0] - h[1,1])/h[0,0])>MINSERIES):
      r01 = (e00 - e11)/(h[0,0] - h[1,1])
   else:
      r01 = e00*(1 + (h[0,0] - h[1,1])*(0.5 + (h[0,0] - h[1,1])/6.0))
   if (abs((h[1,1] - h[2,2])/h[1,1])>MINSERIES):
      r12 = (e11 - e22)/(h[1,1] - h[2,2])
   else:
      r12 = e11*(1 + (h[1,1] - h[2,2])*(0.5 + (h[1,1] - h[2,2])/6.0))
   if (abs((h[2,2] - h[0,0])/h[2,2])>MINSERIES):
      r02 = (e22 - e00)/(h[2,2] - h[0,0])
   else:
      r02 = e22*(1 + (h[2,2] - h[0,0])*(0.5 + (h[2,2] - h[0,0])/6.0))

   eh[0,1] = h[0,1]*r01
   eh[1,2] = h[1,2]*r12

   eh[0,2] = h[0,2]*r02
   if (abs((h[2,2] - h[0,0])/h[2,2])>MINSERIES):
      eh[0,2] += h[0,1]*h[0,2]*(r01 - r12)/(h[0,0] - h[2,2])
   elif (abs((h[1,1] - h[0,0])/h[1,1])>MINSERIES):
      eh[0,2] += h[0,1]*h[0,2]*(r12 - r02)/(h[1,1] - h[0,0])
   elif (abs((h[1,1]-h[2,2])/h[1,1])>MINSERIES):
      eh[0,2] += h[0,1]*h[0,2]*(r02 - r01)/(h[2,2] - h[1,1])
   else:
      eh[0,2] += h[0,1]*h[0,2]*e00/24.0*(12.0 + 4*(h[1,1] + h[2,2] - 2*h[0,0]) + (h[1,1] - h[0,0])*(h[2,2] - h[0,0]))

   return eh

def root_herm(A):
   """Gives the square root of a hermitian matrix with real eigenvalues.

   Args:
      A: A Hermitian matrix.

   Returns:
      A matrix such that itself matrix multiplied by its transpose gives the
      original matrix.
   """

   if not (abs(A.T - A) < 1e-10).all():
      raise ValueError("Non-Hermitian matrix passed to root_herm function")
   eigvals, eigvecs = np.linalg.eigh(A)
   ndgrs = len(eigvals)
   diag = np.zeros((ndgrs,ndgrs))
   for i in range(ndgrs):
      if eigvals[i] >= 0:
         diag[i,i] = math.sqrt(eigvals[i])
      else:
         warning("Zeroing negative element in matrix square root: " + str(eigvals[i]), verbosity.low)
         diag[i,i] = 0
   return np.dot(eigvecs, np.dot(diag, eigvecs.T))

