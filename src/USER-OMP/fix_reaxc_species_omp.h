/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(reax/c/species/omp,FixReaxCSpeciesOMP)

#else

#ifndef LMP_FIX_REAXC_SPECIES_OMP_H
#define LMP_FIX_REAXC_SPECIES_OMP_H

#include "pair_reaxc_omp.h"
#include "fix_reaxc_species.h"

#define BUFLEN 1000

namespace LAMMPS_NS {
  
  class FixReaxCSpeciesOMP : public FixReaxCSpecies {
    
  public:
    FixReaxCSpeciesOMP(class LAMMPS *, int, char **);
    ~FixReaxCSpeciesOMP(){};
    virtual void init();
    
  private:
    class PairReaxCOMP *reaxc;
    
  };
}

#endif
#endif
