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

FixStyle(reax/c/species,FixReaxCSpecies)

#else

#ifndef LMP_FIX_REAXC_SPECIES_H
#define LMP_FIX_REAXC_SPECIES_H

#include "stdio.h"
#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixReaxCSpecies : public Fix {
 public:
  FixReaxCSpecies(class LAMMPS *, int, char **);
  ~FixReaxCSpecies();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();
  double compute_vector(int);

 private:
  class PairReaxC *reaxc;
  class NeighList *list;

  int me, nprocs, nmax, nlocal, ntypes, ntotal;
  int nrepeat, irepeat, repeat, nfreq, posfreq;
  int Nmoltype, vector_nmole, vector_nspec;
  int *Name, *MolName, *NMol, *nd, *MolType, *molmap;
  double *clusterID;

  double bg_cut;
  double *tmpq, **tmpx;
  double **BOCut, **abo;

  FILE *fp, *pos;
  int eleflag, posflag, padflag, multipos;
  int singlepos_opened, multipos_opened;
  char *ele, **eletype, *posspec, *filepos;

  void Output_ReaxC_Bonds(bigint, FILE *);
  void GatherBondOrder(struct _reax_list*);
  void FindMolecule();
  void SortMolecule(int &);
  void FindSpecies(int, int &);
  void WriteFormulae(int, int);
  void OpenPos();
  void WritePos(int, int);

  int CheckExistence(int, int);
  int nint(const double &);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  double memory_usage();

  bigint nvalid, nextvalid();
  struct _reax_list *lists;

};
}

#endif
#endif
