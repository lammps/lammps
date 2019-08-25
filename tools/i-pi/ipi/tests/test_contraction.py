"""Tests ring polymer contraction.

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
"""

import sys
sys.path.append("../")
sys.path.append("../../")

from ipi.utils import nmtransform
import numpy as np
from numpy.testing import assert_almost_equal as assert_equals

def check_up_and_down_scaling(n, q):
    """Check if q expanding and then contracting a ring polymer is a no-op.

    Args:
       n: The number of beads in the scaled ring polymer.
       q: The original position array.
    """

    rescale = nmtransform.nm_rescale(q.shape[0], n)
    print "Initial position of the beads:"
    print q, q.shape, (q.shape[0], n)

    # rescale up to the n beads
    beads_n = rescale.b1tob2(q)
    print "Upscaled to %d beads:"%n
    print beads_n, beads_n.shape

    beads_final = rescale.b2tob1(beads_n)
    print "Final position of the beads:"
    print beads_final

    assert_equals(q, beads_final)
    return beads_n

def check_rpc_consistency(n, q):
   """Check if q expanding and then contracting a ring polymer is a no-op.

   Args:
      n: The number of beads in the scaled ring polymer.
      q: The original position array.
   """

   rescale1 = nmtransform.nm_rescale(q.shape[0], n)
   rescale2 = nmtransform.nm_rescale(n,q.shape[0])

   beads_n=rescale1.b1tob2(q)
   beads_1=rescale1.b2tob1(beads_n)
   beads_2=rescale2.b1tob2(beads_n)

   assert_equals(beads_1, beads_2)

def check_centroid_pos(n, q):
   """Check if expanding and then contracting a ring polymer
   maintains the centroid.

   Args:
      n: The number of beads in the scaled ring polymer.
      q: The original position array.
   """

   beads_big = check_up_and_down_scaling(n, q)
   rescale_big = nmtransform.mk_rs_matrix(n, 1)
   rescale_q = nmtransform.mk_rs_matrix(q.shape[0], 1)

   centroid_big = np.dot(rescale_big, beads_big)
   centroid_q = np.dot(rescale_q, q)

   assert_equals(centroid_q, centroid_big)

numbers_to_check = range(10, 56, 9)
def test_1_to_n():
   """One bead tests."""

   for n in numbers_to_check:
      q = np.array([[0.0,0.0,0.0, 1.0,0.0,0.0]])
      yield check_up_and_down_scaling, n, q
      yield check_rpc_consistency, n, q
      yield check_centroid_pos, n, q

def test_2_to_n():
   """Two bead tests."""

   for n in numbers_to_check:
      q = np.array([[0.0,0.0,0.0, 1.0,0.0,0.0],
                    [0.0,0.1,0.0, 1.0,0.1,0.0]])
      yield check_up_and_down_scaling, n, q
      yield check_rpc_consistency, n, q
      yield check_centroid_pos, n, q

def test_3_to_n():
   """Three bead tests."""

   for n in numbers_to_check:
      q = np.array([[0.0, 0.0,0.0, 1.0, 0.0,0.0],
                    [0.0, 0.1,0.0, 1.0, 0.1,0.0],
                    [0.0,-0.1,0.0, 1.0,-0.1,0.0]])
      yield check_up_and_down_scaling, n, q
      yield check_rpc_consistency, n, q
      yield check_centroid_pos, n, q

def test_4_to_n():
   """Four bead tests."""

   for n in numbers_to_check:
      q = np.array([[0.0, 0.0,0.0, 1.0, 0.0,0.0],
                    [0.0, 0.1,0.0, 1.0, 0.1,0.0],
                    [0.0, 0.2,0.0, 1.0, 0.2,0.0],
                    [0.0,-0.1,0.0, 1.0,-0.1,0.0]])
      yield check_up_and_down_scaling, n, q
      yield check_rpc_consistency, n, q
      yield check_centroid_pos, n, q
