# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

################################################################################
# LAMMPS data structures
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

class NeighList:
    """This is a wrapper class that exposes the contents of a neighbor list.

    It can be used like a regular Python list. Each element is a tuple of:

    * the atom local index
    * its number of neighbors
    * and a pointer to an c_int array containing local atom indices of its
      neighbors

    Internally it uses the lower-level LAMMPS C-library interface.

    :param lmp: reference to instance of :py:class:`lammps`
    :type  lmp: lammps
    :param idx: neighbor list index
    :type  idx: int
    """
    def __init__(self, lmp, idx):
        self.lmp = lmp
        self.idx = idx

    def __str__(self):
        return "Neighbor List ({} atoms)".format(self.size)

    def __repr__(self):
        return self.__str__()

    @property
    def size(self):
        """
        :return: number of elements in neighbor list
        """
        return self.lmp.get_neighlist_size(self.idx)

    def get(self, element):
        """
        :return: tuple with atom local index, numpy array of neighbor local atom indices
        :rtype:  (int, int, ctypes.POINTER(c_int))
        """
        iatom, numneigh, neighbors = self.lmp.get_neighlist_element_neighbors(self.idx, element)
        return iatom, numneigh, neighbors

    # the methods below implement the iterator interface, so NeighList can be used like a regular Python list

    def __getitem__(self, element):
        return self.get(element)

    def __len__(self):
        return self.size

    def __iter__(self):
        inum = self.size

        for ii in range(inum):
            yield self.get(ii)
