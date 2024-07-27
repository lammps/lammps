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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(write_psf,WritePsf);
// clang-format on
#else

#ifndef LMP_WRITE_PSF_H
#define LMP_WRITE_PSF_H

#include "command.h"

namespace LAMMPS_NS {

class WritePsf : public Command {
 public:
  WritePsf(class LAMMPS *);
  void command(int, char **) override;
  void write(const std::string &);

 private:
  int me, nprocs;
  int lmapflag;
  int igroup, groupbit;    // group that WritePsf is performed on
  FILE *fp;
  bigint nbonds_local, nbonds;
  bigint nangles_local, nangles;
  bigint ndihedrals_local, ndihedrals;
  bigint nimpropers_local, nimpropers;
  int **atom_iarray_psf;

  void header();
  void atoms();
  void bonds();
  void angles();
  void dihedrals();
  void impropers();
  int count();

  int size_one;    // # of quantities for one atom
  int nme;         // # of atoms in this dump from me

  bigint ntotal;            // total # of per-atom lines in snapshot
  int reorderflag;          // 1 if OK to reorder instead of sort
  bigint ntotal_reorder;    // # of atoms that must be in snapshot
  int nme_reorder;          // # of atoms I must own in snapshot
  tagint idlo;              // lowest ID I own when reordering

  int maxbuf;     // size of buf
  double *buf;    // memory for atom quantities
  int maxids;     // size of ids
  int maxsort;    // size of bufsort, idsort, index
  int maxproc;    // size of proclist
  tagint *ids;    // list of atom IDs, if sorting on IDs
  double *bufsort;
  tagint *idsort;
  int *index, *proclist;

  void sort();
#if defined(LMP_QSORT)
  static int idcompare(const void *, const void *);
#else
  static int idcompare(const int, const int, void *);
#endif

  class Irregular *irregular;


};


}    // namespace LAMMPS_NS

#endif
#endif
