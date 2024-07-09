/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#ifndef LMP_STL_READER_H
#define LMP_STL_READER_H

#include "pointers.h"

namespace LAMMPS_NS {

class STLReader : protected Pointers {
 public:
  STLReader(class LAMMPS *);
  ~STLReader();
  int read_file(const char *, double **);

 private:
  int ntris,maxtris;
  double **tris;
};

}    // namespace LAMMPS_NS

#endif
