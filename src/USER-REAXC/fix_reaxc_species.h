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
#include "pair_reax_c.h"
#include "reaxc_defs.h"

#define MaxMolTypes 80000
#define BUFLEN 1000

namespace LAMMPS_NS {

class FixReaxCSpecies : public Fix {
 public:
  FixReaxCSpecies(class LAMMPS *, int, char **);
  ~FixReaxCSpecies();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();


 private:
  int me, nprocs, Nmoltype, nmax, ntypes;
  int nrepeat, nfreq, irepeat, fixnvery;
  int *numcount, *nd, *tag, *type;
  int *MolID, *NMol, *MolType, *Mol, *Name, *MolName;
  int multifile, padflag, bondflag;
  double **BO, **BOCut, *MolMass, *masstype;
  double *sbo, *nlp, *avq;
  char *ele, *filename, **eletype;
  FILE *fp, *mass, *bond, *read;

  // converting BO to BL
  int **BLcount;
  double **BLsum;
  FILE *blfp;

  void OpenFile();
  void ReadDict(FILE *);
  void OutputReaxCBonds(bigint, FILE *);
  void PassBuffer(reax_system*, double *, int &);
  void RecvBuffer(reax_system*, double *, int, int);
  void FindBond(reax_system*, int);
  void FindMolecule(reax_system*, int, int &);
  void FindSpecies(reax_system*, int, int , int, int &, int &);
  int CheckExistence(int, int);
  void FindFormulas(int, int, int);
  void PlotSpectra(reax_system*, int, int, int, int);
  void WriteBond(reax_system*, int, int, double *, int);
  void allocate();
/*
  void ClusterID(reax_system*);
  int pack_comm(int, int*, double*, int, int*);
  void unpack_comm(int, int, double*);
*/
  bigint nvalid, nextvalid();
  int nint(const double &);

  reax_system *system;
  class PairReaxC *reaxc;
  class NeighList *list;
};
}

#endif
#endif
