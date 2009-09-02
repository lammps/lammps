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

#ifndef DUMP_CFG_H
#define DUMP_CFG_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpCFG : public DumpCustom {
 public:
  DumpCFG(class LAMMPS *, int, char **);
  ~DumpCFG();
  void init();
  void write_header(int);
  void write_data(int, double *);
  int modify_param2(int, char **);

 private:
  int ntypes;                // # of atom types
  char **typenames;	     // array of element names for each type
  char **auxname;            // name strings of auxiliary properties
  int nchosen;               // # of lines to be written on a writing proc
  int nlines;                // # of lines transferred from buf to rbuf
  int *recvcounts;           // # of data (tag info) received from all procs
  int *displs;               // displacement of data (tag info) within buffer
  int *tags;                 // atom IDs of all the chosen atoms collected
  double **rbuf;             // buf of data lines for data lines rearrangement

  void create_sorted_tags();
  static int tag_compare(const void *, const void *);
};

}

#endif
