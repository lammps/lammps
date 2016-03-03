/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_RANDOM_MT_H
#define LMP_RANDOM_MT_H

#include "pointers.h"

namespace LAMMPS_NS {

class RanMT : protected Pointers {
 public:
  RanMT(class LAMMPS *, int);
  virtual ~RanMT();
  double uniform();
  double gaussian();
  double gaussian_z();

 private:
  // some internal constants
  enum {MT_S=7,MT_U=11,MT_T=15,MT_L=18,MT_R=31,MT_M=397,MT_N=624};
  uint32_t _m[MT_N];            // state of RNG
  int _idx;                     // twist index
  int _save;                    // for gaussian distributed RNG
  int _nlayers;                 // number of layers for ziggurat method
  double *_xlayers, *_ylayers;  // layer distribution for ziggurat method
  double _second;               // stored RNG for gaussian distributed RNG

  uint32_t _randomize();        // generate new 32-bit integer
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Marsaglia random # generator

The initial seed for this random number generator must be a positive
integer less than or equal to 900 million.

*/
