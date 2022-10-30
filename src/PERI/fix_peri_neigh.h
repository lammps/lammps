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
FixStyle(PERI_NEIGH,FixPeriNeigh);
// clang-format on
#else

#ifndef LMP_FIX_PERI_NEIGH_H
#define LMP_FIX_PERI_NEIGH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPeriNeigh : public Fix {
  friend class PairPeri;
  friend class PairPeriPMB;
  friend class PairPeriPMBOMP;
  friend class PairPeriLPS;
  friend class PairPeriVES;
  friend class PairPeriEPS;
  friend class PairPeriLPSOMP;
  friend class ComputeDamageAtom;
  friend class ComputePlasticityAtom;

 public:
  FixPeriNeigh(class LAMMPS *, int, char **);
  ~FixPeriNeigh() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void min_setup(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  int first;                            // flag for first time initialization
  int maxpartner;                       // max # of peridynamic neighs for any atom
  int *npartner;                        // # of neighbors for each atom
  tagint **partner;                     // neighs for each atom, stored as global IDs
  double **deviatorextention;           // Deviatoric extension
  double **deviatorBackextention;       // Deviatoric back extension
  double **deviatorPlasticextension;    // Deviatoric plastic extension
  double *lambdaValue;
  double **r0;                       // initial distance to partners
  double **r1;                       // instanteneous distance to partners
  double *thetaValue;                // dilatation
  double *vinter;                    // sum of vfrac for bonded neighbors
  double *wvolume;                   // weighted volume of particle
  int isPMB, isLPS, isVES, isEPS;    // which flavor of PD

  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
