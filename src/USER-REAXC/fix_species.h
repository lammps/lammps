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

FixStyle(species,FixSpecies)

#else

#ifndef LMP_FIX_SPECIES_H
#define LMP_FIX_SPECIES_H

#include "fix.h"
#include "pointers.h"

#include "pair_reax_c.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"
// #include "pair_foo.h"

#define MAXSPECBOND 12
#define BUFLEN 1000

namespace LAMMPS_NS {

class FixSpecies : public Fix {
 public:
  FixSpecies(class LAMMPS *, int, char **);
  ~FixSpecies();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void end_of_step();

 private:
  int eleflag, posflag, multipos, padflag, setupflag;
  int me, nprocs, nmax, nlocal, nghost, nall, ntypes, ntotal;
  int nrepeat, nfreq, Nmoltype;
  int *Name, *MolName, *NMol, *nd, *MolType, *molmap, *clusterID;

  double bg_cut, **BOCut;

  char *ele, **eletype;
  char **tmparg;
  char *posspec, *filepos;

  FILE *fp, *pos;

  void Output_ReaxC_Bonds(bigint, FILE *);
  void create_compute();
  void create_fix();
  void FindMolecule();
  void SortMolecule(int &);
  void FindSpecies(int, int &);
  void WriteFormulas(int, int);
  int CheckExistence(int, int);

  int nint(const double &);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  double memory_usage();

  bigint nvalid;

  class NeighList *list;
  class FixAveAtom *f_SPECBOND;
  class PairReaxC *reaxc;
  // class PairFoo *foo;

};
}

#endif
#endif
