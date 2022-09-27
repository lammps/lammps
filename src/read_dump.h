/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(read_dump,ReadDump);
// clang-format on
#else

#ifndef LMP_READ_DUMP_H
#define LMP_READ_DUMP_H

#include "command.h"

namespace LAMMPS_NS {

class ReadDump : public Command {
 public:
  ReadDump(class LAMMPS *);
  ~ReadDump() override;
  void command(int, char **) override;

  void store_files(int, char **);
  void setup_reader(int, char **);
  bigint seek(bigint, int);
  void header(int);
  bigint next(bigint, bigint, int, int);
  void atoms();
  int fields_and_keywords(int, char **);

 private:
  int me, nprocs;

  char **files;       // list of input dump files to process
  int nfile;          // # of dump files to process (each may be parallel)
  int currentfile;    // current open file (0 to nfile-1)

  MPI_Comm clustercomm;              // comm for proc cluster that reads/shares a file
  int me_cluster, nprocs_cluster;    // proc ID and count for my read cluster

  int multiproc;          // 0 = each dump file is a single file
                          // 1 = each dump file is parallel (multiple files)
  int multiproc_nfile;    // number of parallel files in one dump file

  int nreader;       // # of parallel dump files read by my cluster
  int firstfile;     // index of 1st dump file my cluster reads
                     //   (0 to multiproc_nfile-1)
  int filereader;    // 1 if this proc reads from a dump file(s)
  int parallel;      // 1 if parallel reading (e.g. via ADIOS2)

  int dimension;    // same as in Domain
  int triclinic;

  int boxflag;                 // overwrite simulation box with dump file box params
  int timestepflag;            // overwrite simulation timestep with dump file timestep
  int replaceflag, addflag;    // flags for processing dump snapshot atoms
  int trimflag, purgeflag;
  int scaleflag;        // user 0/1 if dump file coords are unscaled/scaled
  int wrapflag;         // user 0/1 if dump file coords are unwrapped/wrapped
  char *readerstyle;    // style of dump files to read

  int nnew;             // # of dump file atoms this proc owns
  int nfield;           // # of fields to extract from dump file
  int *fieldtype;       // type of each field = X,VY,IZ,etc
  char **fieldlabel;    // user specified label for field
  double **fields;      // per-atom field values
  int maxnew;           // allocation size of fields array
  double **buf;         // read buffer

  int scaled;     // 0/1 if dump file coords are unscaled/scaled
  int wrapped;    // 0/1 if dump file coords are unwrapped/wrapped

  double box[3][3];                                   // dump file box parameters
  double xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;    // dump snapshot box params
  double xprd, yprd, zprd;

  bigint *nsnapatoms;    // # of atoms in one snapshot from
                         //   one (parallel) dump file
                         // nreader-length vector b/c a reader proc
                         //   may read from multiple parallel dump files

  int npurge, nreplace, ntrim, nadd;    // stats on processed atoms
  int yindex, zindex;                   // field index for Y,Z coords

  class Reader **readers;    // class that reads a dump file
                             // nreader-length list of readers if proc reads
                             //   from multiple parallel dump files

  void read_atoms();
  void process_atoms();
  void migrate_old_atoms();
  void migrate_new_atoms();
  void migrate_atoms_by_coords();

  void setup_multiproc();
  int whichtype(char *);

  double xfield(int, int);
  double yfield(int, int);
  double zfield(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
