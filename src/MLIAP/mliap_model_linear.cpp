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

#include "mliap_model_linear.h"
#include "pair_mliap.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::MLIAPModelLinear(LAMMPS* lmp, char* coefffilename, PairMLIAP* pairmliap_in) : 
  MLIAPModel(lmp, coefffilename, pairmliap_in)
{
  nonlinearflag = 0;
  ndescriptors = ncoeffall - 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelLinear::~MLIAPModelLinear(){}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors for each atom dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelLinear::gradient(NeighList* list, double **descriptors, double **beta, int eflag)
{
  int i;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = pairmliap->map[itype];
    double* coeffi = coeffelem[ielem];

    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation

    if (eflag) {

      // energy of atom I

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E_i = beta.B_i

      for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*descriptors[ii][icoeff];
      
      pairmliap->e_tally(i,etmp);
    }
  }
}

