/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(read_dump,ReadDump)

#else

#ifndef LMP_READ_DUMP_H
#define LMP_READ_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadDump : protected Pointers {
 public:
  ReadDump(class LAMMPS *);
  ~ReadDump() {}
  void command(int, char **);

private:
  int me,nprocs;
  FILE *fp;
  int dimension;
  int triclinic;

  bigint nstep;            // timestep to find in dump file
  int boxflag;             // use dump file box params
  int replaceflag,addflag; // flags for processing dump snapshot atoms
  int trimflag,purgeflag;
  int scaledflag;          // user setting for coordinate scaling
  int scaled;              // actual setting for coordinate scaling
  int format;              // style of dump file
  int compressed;          // flag for dump file compression

  int nfield;              // # of fields to extract from dump file
  int *fieldtype;          // type of each field = X,VY,IZ,etc
  char **fieldlabel;       // user specified label for field
  double **fields;         // per-atom field values

  double box[3][3];         // dump file box parameters
  double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;
  double xprd,yprd,zprd;

  bigint nsnapatoms;        // # of atoms in dump file shapshot

  int npurge,nreplace,ntrim,nadd;     // stats on processed atoms
  int addproc;                        // proc that should add next atom
  int yindex,zindex;                  // field index for Y,Z coords

  int *uflag;               // set to 1 if snapshot atom matches owned atom
  int *ucflag,*ucflag_all;  // set to 1 if snapshot chunk atom was processed

  class ReadDumpNative *reader;     // class that reads native dump file

  void process_atoms(int);
  void delete_atoms(int *);

  double xfield(int, int);
  double yfield(int, int);
  double zfield(int, int);

  void open(char *);
};

}

#endif
#endif
