// clang-format off
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

#ifdef FIX_CLASS
// clang-format off
FixStyle(reaxff/species,FixReaxFFSpecies);
FixStyle(reax/c/species,FixReaxFFSpecies);
// clang-format on
#else

#ifndef LMP_FIX_REAXC_SPECIES_H
#define LMP_FIX_REAXC_SPECIES_H

#include "fix.h"

#define BUFLEN 1000

namespace LAMMPS_NS {

typedef struct {
  double x, y, z;
} AtomCoord;

class FixReaxFFSpecies : public Fix {
 public:
  FixReaxFFSpecies(class LAMMPS *, int, char **);
  ~FixReaxFFSpecies() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void post_integrate() override;
  double compute_vector(int) override;

 protected:
  int nmax, nlocal, ntypes, ntotal;
  int nrepeat, nfreq, posfreq, compressed, ndelspec;
  int Nmoltype, vector_nmole, vector_nspec;
  int *Name, *MolName, *NMol, *nd, *MolType, *molmap, *mark;
  int *Mol2Spec;
  double *clusterID;
  AtomCoord *x0;

  double bg_cut;
  double **BOCut;

  std::vector<std::string> del_species;

  FILE *fp, *pos, *fdel;
  int eleflag, posflag, multipos, padflag, setupflag;
  int delflag, specieslistflag, masslimitflag;
  int delete_Nlimit, delete_Nlimit_varid;
  std::string delete_Nlimit_varname;
  int delete_Nsteps, *delete_Tcount;
  double massmin, massmax;
  int singlepos_opened, multipos_opened, del_opened;
  char *ele, **eletype, *filepos, *filedel;

  void Output_ReaxFF_Bonds(bigint, FILE *);
  AtomCoord chAnchor(AtomCoord, AtomCoord);
  virtual void FindMolecule();
  void SortMolecule(int &);
  void FindSpecies(int, int &);
  void WriteFormulas(int, int);
  void DeleteSpecies(int, int);
  int CheckExistence(int, int);

  int nint(const double &);
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void OpenPos();
  void WritePos(int, int);
  double memory_usage() override;

  bigint nvalid;

  class NeighList *list;
  class FixAveAtom *f_SPECBOND;
  class PairReaxFF *reaxff;
};
}    // namespace LAMMPS_NS

#endif
#endif
