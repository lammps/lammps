/* -*- c++ -*- ----------------------------------------------------------
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

#include <mpi.h>
#include "pointers.h"

namespace LAMMPS_NS {

class ReadDump : protected Pointers {
 public:
  ReadDump(class LAMMPS *);
  ~ReadDump();
  void command(int, char **);

  void store_files(int, char **);
  void setup_reader(int, char **);
  bigint seek(bigint, int);
  void header(int);
  bigint next(bigint, bigint, int, int);
  void atoms();
  int fields_and_keywords(int, char **);

private:
  int me,nprocs;

  char **files;            // list of input dump files to process
  int nfile;               // # of dump files to process (each may be parallel)
  int currentfile;         // current open file (0 to nfile-1)

  MPI_Comm clustercomm;    // comm for proc cluster that reads/shares a file
  int me_cluster,nprocs_cluster;  // proc ID and count for my read cluster

  int multiproc;           // 0 = each dump file is a single file
                           // 1 = each dump file is parallel (multiple files)
  int multiproc_nfile;     // number of parallel files in one dump file

  int nreader;             // # of parallel dump files read by my cluster
  int firstfile;           // index of 1st dump file my cluster reads
                           //   (0 to multiproc_nfile-1)
  int filereader;          // 1 if this proc reads from a dump file(s)

  int dimension;           // same as in Domain
  int triclinic;

  int boxflag;             // overwrite simulation with dump file box params
  int replaceflag,addflag; // flags for processing dump snapshot atoms
  int trimflag,purgeflag;
  int scaleflag;           // user 0/1 if dump file coords are unscaled/scaled
  int wrapflag;            // user 0/1 if dump file coords are unwrapped/wrapped
  char *readerstyle;       // style of dump files to read

  int nnew;                // # of dump file atoms this proc owns
  int nfield;              // # of fields to extract from dump file
  int *fieldtype;          // type of each field = X,VY,IZ,etc
  char **fieldlabel;       // user specified label for field
  double **fields;         // per-atom field values
  int maxnew;              // allocation size of fields array
  double **buf;            // read buffer

  int scaled;              // 0/1 if dump file coords are unscaled/scaled
  int wrapped;             // 0/1 if dump file coords are unwrapped/wrapped

  double box[3][3];         // dump file box parameters
  double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;  // dump snapshot box params
  double xprd,yprd,zprd;

  bigint *nsnapatoms;       // # of atoms in one snapshot from
                            //   one (parallel) dump file
                            // nreader-length vector b/c a reader proc
                            //   may read from multiple parallel dump files

  int npurge,nreplace,ntrim,nadd;     // stats on processed atoms
  int yindex,zindex;                  // field index for Y,Z coords

  class Reader **readers;   // class that reads a dump file
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

}

#endif
#endif

/* ERROR/WARNING messages:

E: Read_dump command before simulation box is defined

The read_dump command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump file does not contain requested snapshot

Self-explanatory.

E: Unrecognized dump reader style

The choice of dump reader style via the format keyword is unknown.

E: No box information in dump, must use 'box no'

UNDOCUMENTED

E: Read_dump triclinic status does not match simulation

Both the dump snapshot and the current LAMMPS simulation must
be using either an orthogonal or triclinic box.

E: Read_dump field not found in dump file

Self-explanatory.

E: Read_dump xyz fields do not have consistent scaling/wrapping

Self-explanatory.

E: All read_dump x,y,z fields must be specified for scaled, triclinic coords

For triclinic boxes and scaled coordinates you must specify all 3 of
the x,y,z fields, else LAMMPS cannot reconstruct the unscaled
coordinates.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: Read dump of atom property that isn't allocated

Self-explanatory.

E: Duplicate fields in read_dump command

Self-explanatory.

E: If read_dump purges it cannot replace or trim

These operations are not compatible.  See the read_dump doc
page for details.

E: Read_dump cannot use 'add keep' without atom IDs

UNDOCUMENTED

E: Cannot add atoms if dump file does not store atom type

UNDOCUMENTED

U: No box information in dump. You have to use 'box no'

Self-explanatory.

*/
