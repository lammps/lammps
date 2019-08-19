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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "reader_xyz.h"
#include <cstdlib>
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024        // max line length in dump file

enum{ID,TYPE,X,Y,Z};

/* ---------------------------------------------------------------------- */

ReaderXYZ::ReaderXYZ(LAMMPS *lmp) : Reader(lmp)
{
  line = new char[MAXLINE];
  fieldindex = NULL;
  nstep = 0;
}

/* ---------------------------------------------------------------------- */

ReaderXYZ::~ReaderXYZ()
{
  delete [] line;
  memory->destroy(fieldindex);
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderXYZ::read_time(bigint &ntimestep)
{
  char *eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) return 1;

  // first line has to have the number of atoms
  // truncate the string to the first whitespace,
  //   so force->bnumeric() does not hiccup

  for (int i=0; (i < MAXLINE) && (eof[i] != '\0'); ++i) {
    if (eof[i] == '\n' || eof[i] == '\r' || eof[i] == ' ' || eof[i] == '\t') {
      eof[i] = '\0';
      break;
    }
  }
  natoms = force->bnumeric(FLERR,line);
  if (natoms < 1)
    error->one(FLERR,"Dump file is incorrectly formatted");

  // skip over comment/title line

  read_lines(1);

  // fake time step numbers

  ntimestep = nstep;

  // count this frame

  ++nstep;
  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderXYZ::skip()
{
  // invoke read_lines() in chunks no larger than MAXSMALLINT

  int nchunk;
  bigint nremain = natoms;
  while (nremain) {
    nchunk = MIN(nremain,MAXSMALLINT);
    read_lines(nchunk);
    nremain -= nchunk;
  }
}

/* ----------------------------------------------------------------------
   read remaining header info:
     return natoms
     box bounds, triclinic (inferred), fieldflag (1 if any fields not found),
     xyz flags = from input scaleflag & wrapflag
   if fieldflag set:
     match Nfield fields to per-atom column labels
     allocate and set fieldindex = which column each field maps to
     fieldtype = X,VX,IZ etc
     fieldlabel = user-specified label or NULL if use fieldtype default
   xyz flag = scaledflag if has fieldlabel name, else set by x,xs,xu,xsu
   only called by proc 0
------------------------------------------------------------------------- */

bigint ReaderXYZ::read_header(double /*box*/[3][3], int &boxinfo, int &/*triclinic*/,
                              int fieldinfo, int nfield,
                              int *fieldtype, char **/*fieldlabel*/,
                              int scaleflag, int wrapflag, int &fieldflag,
                              int &xflag, int &yflag, int &zflag)
{

  nid = 0;

  // signal that we have no box info at all

  boxinfo = 0;

  // if no field info requested, just return

  if (!fieldinfo) return natoms;

  memory->create(fieldindex,nfield,"read_dump:fieldindex");

  // for xyz we know nothing about the style of coordinates,
  // so caller has to set the proper flags

  xflag = 2*scaleflag + wrapflag + 1;
  yflag = 2*scaleflag + wrapflag + 1;
  zflag = 2*scaleflag + wrapflag + 1;

  // copy fieldtype list for supported fields

  fieldflag = 0;
  for (int i = 0; i < nfield; i++) {
    if ( (fieldtype[i] == X) ||
         (fieldtype[i] == Y) ||
         (fieldtype[i] == Z) ||
         (fieldtype[i] == ID) ||
         (fieldtype[i] == TYPE) ) {
      fieldindex[i] = fieldtype[i];
    } else {
      fieldflag = 1;
    }
  }

  return natoms;
}

/* ----------------------------------------------------------------------
   read N atom lines from dump file
   stores appropriate values in fields array
   return 0 if success, 1 if error
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderXYZ::read_atoms(int n, int nfield, double **fields)
{
  int i,m,rv;
  char *eof;
  int mytype;
  double myx, myy, myz;

  for (i = 0; i < n; i++) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");

    ++nid;
    rv = sscanf(line,"%*s%lg%lg%lg", &myx, &myy, &myz);
    if (rv != 3)
      error->one(FLERR,"Dump file is incorrectly formatted");

    // XXX: we could insert an element2type translation here
    // XXX: for now we flag unrecognized types as type 0,
    // XXX: which should trigger an error, if LAMMPS uses it.
    mytype = atoi(line);

    for (m = 0; m < nfield; m++) {
      switch (fieldindex[m]) {
      case X:
        fields[i][m] = myx;
        break;
      case Y:
        fields[i][m] = myy;
        break;
      case Z:
        fields[i][m] = myz;
        break;
      case ID:
        fields[i][m] = nid;
        break;
      case TYPE:
        fields[i][m] = mytype;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   read N lines from dump file
   only last one is saved in line
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderXYZ::read_lines(int n)
{
  char *eof = NULL;
  if (n <= 0) return;
  for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");
}
