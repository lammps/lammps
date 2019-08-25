"""Contains functions for doing the inverse and forward normal mode transforms.

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


Classes:
   nm_trans: Uses matrix multiplication to do normal mode transformations.
   nm_rescale: Uses matrix multiplication to do ring polymer contraction
      or expansion.
   nm_fft: Uses fast-Fourier transforms to do normal modes transformations.

Functions:
   mk_nm_matrix: Makes a matrix to transform between the normal mode and bead
      representations.
   mk_rs_matrix: Makes a matrix to transform between one number of beads and
      another. Higher normal modes in the case of an expansion are set to zero.
"""

__all__ = ['nm_trans', 'nm_rescale', 'nm_fft']

import numpy as np
from ipi.utils.messages import verbosity, info

def mk_nm_matrix(nbeads):
   """Gets the matrix that transforms from the bead representation
   to the normal mode representation.

   If we return from this function a matrix C, then we transform between the
   bead and normal mode representation using q_nm = C . q_b, q_b = C.T . q_nm

   Args:
      nbeads: The number of beads.
   """

   b2nm = np.zeros((nbeads,nbeads))
   b2nm[0,:] = np.sqrt(1.0)
   for j in range(nbeads):
      for i in range(1, nbeads/2+1):
         b2nm[i,j] = np.sqrt(2.0)*np.cos(2*np.pi*j*i/float(nbeads))
      for i in range(nbeads/2+1, nbeads):
         b2nm[i,j] = np.sqrt(2.0)*np.sin(2*np.pi*j*i/float(nbeads))
   if (nbeads%2) == 0:
      b2nm[nbeads/2,0:nbeads:2] = 1.0
      b2nm[nbeads/2,1:nbeads:2] = -1.0
   return b2nm/np.sqrt(nbeads)

def mk_rs_matrix(nb1, nb2):
   """Gets the matrix that transforms a path with nb1 beads into one with
   nb2 beads.

   If we return from this function a matrix T, then we transform between the
   system with nb1 bead and the system of nb2 beads using q_2 = T . q_1

   Args:
      nb1: The initial number of beads.
      nb2: The final number of beads.
   """

   if (nb1 == nb2):
      return np.identity(nb1,float)
   elif (nb1 > nb2):
      b1_nm = mk_nm_matrix(nb1)
      nm_b2 = mk_nm_matrix(nb2).T

      #builds the "reduction" matrix that picks the normal modes we want to keep
      b1_b2 = np.zeros((nb2, nb1), float)
      b1_b2[0,0] = 1.0
      for i in range(1, nb2/2+1):
         b1_b2[i,i] = 1.0
         b1_b2[nb2-i, nb1-i] = 1.0
      if (nb2 % 2 == 0):
         #if we are contracting down to an even number of beads, then we have to
         #pick just one of the last degenerate modes to match onto the single
         #stiffest mode in the new path
         b1_b2[nb2/2, nb1-nb2/2] = 0.0

      rs_b1_b2 = np.dot(nm_b2, np.dot(b1_b2, b1_nm))
      return rs_b1_b2*np.sqrt(float(nb2)/float(nb1))
   else:
      return mk_rs_matrix(nb2, nb1).T*(float(nb2)/float(nb1))


class nm_trans:
   """Helper class to perform beads <--> normal modes transformation.

   Attributes:
      _b2nm: The matrix to transform between the bead and normal mode
         representations.
      _nm2b: The matrix to transform between the normal mode and bead
         representations.
   """

   def __init__(self, nbeads):
      """Initializes nm_trans.

      Args:
         nbeads: The number of beads.
      """

      self._b2nm = mk_nm_matrix(nbeads)
      self._nm2b = self._b2nm.T

   def b2nm(self, q):
      """Transforms a matrix to the normal mode representation.

      Args:
         q: A matrix with nbeads rows, in the bead representation.
      """

      return np.dot(self._b2nm,q)

   def nm2b(self, q):
      """Transforms a matrix to the bead representation.

      Args:
         q: A matrix with nbeads rows, in the normal mode representation.
      """

      return np.dot(self._nm2b,q)


class nm_rescale:
   """Helper class to rescale a ring polymer between different number of beads.

   Attributes:
      _b1tob2: The matrix to transform between a ring polymer with 'nbeads1'
         beads and another with 'nbeads2' beads.
      _b2tob1: The matrix to transform between a ring polymer with 'nbeads2'
         beads and another with 'nbeads1' beads.
   """

   def __init__(self, nbeads1, nbeads2):
      """Initializes nm_rescale.

      Args:
         nbeads1: The initial number of beads.
         nbeads2: The rescaled number of beads.
      """

      self._b1tob2 = mk_rs_matrix(nbeads1,nbeads2)
      self._b2tob1 = self._b1tob2.T*(float(nbeads1)/float(nbeads2))

   def b1tob2(self, q):
      """Transforms a matrix from one value of beads to another.

      Args:
         q: A matrix with nbeads1 rows, in the bead representation.
      """

      return np.dot(self._b1tob2,q)

   def b2tob1(self, q):
      """Transforms a matrix from one value of beads to another.

      Args:
         q: A matrix with nbeads2 rows, in the bead representation.
      """

      return np.dot(self._b2tob1,q)



class nm_fft:
   """Helper class to perform beads <--> normal modes transformation
      using Fast Fourier transforms.

   Attributes:
      fft: The fast-Fourier transform function to transform between the
         bead and normal mode representations.
      ifft: The inverse fast-Fourier transform function to transform
         between the normal mode and bead representations.
      qdummy: A matrix to hold a copy of the bead positions to transform
         them to the normal mode representation.
      qnmdummy: A matrix to hold a copy of the normal modes to transform
         them to the bead representation.
      nbeads: The number of beads.
      natoms: The number of atoms.
   """

   def __init__(self, nbeads, natoms):
      """Initializes nm_trans.

      Args:
         nbeads: The number of beads.
         natoms: The number of atoms.
      """

      self.nbeads = nbeads
      self.natoms = natoms
      try:
         import pyfftw
         info("Import of PyFFTW successful", verbosity.medium)
         self.qdummy = pyfftw.n_byte_align_empty((nbeads, 3*natoms), 16, 'float32')
         self.qnmdummy = pyfftw.n_byte_align_empty((nbeads//2+1, 3*natoms), 16, 'complex64')
         self.fft = pyfftw.FFTW(self.qdummy, self.qnmdummy, axes=(0,), direction='FFTW_FORWARD')
         self.ifft = pyfftw.FFTW(self.qnmdummy, self.qdummy, axes=(0,), direction='FFTW_BACKWARD')
      except ImportError: #Uses standard numpy fft library if nothing better
                          #is available
         info("Import of PyFFTW unsuccessful, using NumPy library instead", verbosity.medium)
         self.qdummy = np.zeros((nbeads,3*natoms), dtype='float32')
         self.qnmdummy = np.zeros((nbeads//2+1,3*natoms), dtype='complex64')
         def dummy_fft(self):
            self.qnmdummy = np.fft.rfft(self.qdummy, axis=0)
         def dummy_ifft(self):
            self.qdummy = np.fft.irfft(self.qnmdummy, n=self.nbeads, axis=0)
         self.fft = lambda: dummy_fft(self)
         self.ifft = lambda: dummy_ifft(self)

   def b2nm(self, q):
      """Transforms a matrix to the normal mode representation.

      Args:
         q: A matrix with nbeads rows and 3*natoms columns,
            in the bead representation.
      """

      if self.nbeads == 1:
         return q
      self.qdummy[:] = q
      self.fft()
      if self.nbeads == 2:
         return self.qnmdummy.real/np.sqrt(self.nbeads)

      nmodes = self.nbeads/2

      self.qnmdummy /= np.sqrt(self.nbeads)
      qnm = np.zeros(q.shape)
      qnm[0,:] = self.qnmdummy[0,:].real

      if self.nbeads % 2 == 0:
         self.qnmdummy[1:-1,:] *= np.sqrt(2)
         (qnm[1:nmodes,:], qnm[self.nbeads:nmodes:-1,:]) = (self.qnmdummy[1:-1,:].real, self.qnmdummy[1:-1,:].imag)
         qnm[nmodes,:] = self.qnmdummy[nmodes,:].real
      else:
         self.qnmdummy[1:,:] *= np.sqrt(2)
         (qnm[1:nmodes+1,:], qnm[self.nbeads:nmodes:-1,:]) = (self.qnmdummy[1:,:].real, self.qnmdummy[1:,:].imag)

      return qnm

   def nm2b(self, qnm):
      """Transforms a matrix to the bead representation.

      Args:
         qnm: A matrix with nbeads rows and 3*natoms columns,
            in the normal mode representation.
      """

      if self.nbeads == 1:
         return qnm
      if self.nbeads == 2:
         self.qnmdummy[:] = qnm
         self.ifft()
         return self.qdummy*np.sqrt(self.nbeads)

      nmodes = self.nbeads/2
      odd = self.nbeads - 2*nmodes  # 0 if even, 1 if odd

      qnm_complex = np.zeros((nmodes+1, len(qnm[0,:])), complex)
      qnm_complex[0,:] = qnm[0,:]
      if not odd:
         (qnm_complex[1:-1,:].real, qnm_complex[1:-1,:].imag) = (qnm[1:nmodes,:], qnm[self.nbeads:nmodes:-1,:])
         qnm_complex[1:-1,:] /= np.sqrt(2)
         qnm_complex[nmodes,:] = qnm[nmodes,:]
      else:
         (qnm_complex[1:,:].real, qnm_complex[1:,:].imag) = (qnm[1:nmodes+1,:], qnm[self.nbeads:nmodes:-1,:])
         qnm_complex[1:,:] /= np.sqrt(2)

      self.qnmdummy[:] = qnm_complex
      self.ifft()
      return self.qdummy*np.sqrt(self.nbeads)
