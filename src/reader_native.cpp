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

#include "reader_native.h"
#include <cstring>
#include <cstdlib>
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024        // max line length in dump file

// also in read_dump.cpp

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,Q,IX,IY,IZ,FX,FY,FZ};
enum{UNSET,NOSCALE_NOWRAP,NOSCALE_WRAP,SCALE_NOWRAP,SCALE_WRAP};

/* ---------------------------------------------------------------------- */

ReaderNative::ReaderNative(LAMMPS *lmp) : Reader(lmp)
{
  line = new char[MAXLINE];
  words = NULL;
  fieldindex = NULL;
}

/* ---------------------------------------------------------------------- */

ReaderNative::~ReaderNative()
{
  delete [] line;
  delete [] words;
  memory->destroy(fieldindex);
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderNative::read_time(bigint &ntimestep)
{
  char *eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) return 1;

  // skip over unit and time information, if present.

  if (utils::strmatch(line,"^\\s*ITEM: UNITS\\s*$"))
    read_lines(2);

  if (utils::strmatch(line,"^\\s*ITEM: TIME\\s*$"))
    read_lines(2);

  if (!utils::strmatch(line,"^\\s*ITEM: TIMESTEP\\s*$"))
    error->one(FLERR,"Dump file is incorrectly formatted");

  read_lines(1);
  int rv = sscanf(line,BIGINT_FORMAT,&ntimestep);
  if (rv != 1)
    error->one(FLERR,"Dump file is incorrectly formatted");

  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::skip()
{
  read_lines(2);
  bigint natoms;
  int rv = sscanf(line,BIGINT_FORMAT,&natoms);
  if (rv != 1)
    error->one(FLERR,"Dump file is incorrectly formatted");

  read_lines(5);

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
     xyz flags = UNSET (not a requested field), SCALE/WRAP as in enum
   if fieldflag set:
     match Nfield fields to per-atom column labels
     allocate and set fieldindex = which column each field maps to
     fieldtype = X,VX,IZ etc
     fieldlabel = user-specified label or NULL if use fieldtype default
   xyz flags = scaleflag+wrapflag if has fieldlabel name,
     else set by x,xs,xu,xsu
   only called by proc 0
------------------------------------------------------------------------- */

bigint ReaderNative::read_header(double box[3][3], int &boxinfo, int &triclinic,
                                 int fieldinfo, int nfield,
                                 int *fieldtype, char **fieldlabel,
                                 int scaleflag, int wrapflag, int &fieldflag,
                                 int &xflag, int &yflag, int &zflag)
{
  bigint natoms;
  int rv;

  read_lines(2);
  rv = sscanf(line,BIGINT_FORMAT,&natoms);
  if (rv != 1)
    error->one(FLERR,"Dump file is incorrectly formatted");

  boxinfo = 1;
  triclinic = 0;
  box[0][2] = box[1][2] = box[2][2] = 0.0;
  read_lines(1);
  if (line[strlen("ITEM: BOX BOUNDS ")] == 'x') triclinic = 1;

  read_lines(1);
  if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[0][0],&box[0][1]);
  else rv = 3 - sscanf(line,"%lg %lg %lg",&box[0][0],&box[0][1],&box[0][2]);
  if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

  read_lines(1);
  if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[1][0],&box[1][1]);
  else rv = 3 - sscanf(line,"%lg %lg %lg",&box[1][0],&box[1][1],&box[1][2]);
  if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

  read_lines(1);
  if (!triclinic) rv = 2 - sscanf(line,"%lg %lg",&box[2][0],&box[2][1]);
  else rv = 3 - sscanf(line,"%lg %lg %lg",&box[2][0],&box[2][1],&box[2][2]);
  if (rv != 0) error->one(FLERR,"Dump file is incorrectly formatted");

  read_lines(1);

  // if no field info requested, just return

  if (!fieldinfo) return natoms;

  // exatract column labels and match to requested fields

  char *labelline = &line[strlen("ITEM: ATOMS ")];

  nwords = atom->count_words(labelline);
  char **labels = new char*[nwords];
  labels[0] = strtok(labelline," \t\n\r\f");
  if (labels[0] == NULL) {
    delete[] labels;
    return 1;
  }
  for (int m = 1; m < nwords; m++) {
    labels[m] = strtok(NULL," \t\n\r\f");
    if (labels[m] == NULL) {
      delete[] labels;
      return 1;
    }
  }

  // match each field with a column of per-atom data
  // if fieldlabel set, match with explicit column
  // else infer one or more column matches from fieldtype
  // xyz flag set by scaleflag + wrapflag (if fieldlabel set) or column label

  memory->create(fieldindex,nfield,"read_dump:fieldindex");

  int s_index,u_index,su_index;
  xflag = UNSET;
  yflag = UNSET;
  zflag = UNSET;

  for (int i = 0; i < nfield; i++) {
    if (fieldlabel[i]) {
      fieldindex[i] = find_label(fieldlabel[i],nwords,labels);
      if (fieldtype[i] == X) xflag = 2*scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Y) yflag = 2*scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Z) zflag = 2*scaleflag + wrapflag + 1;
    }

    else if (fieldtype[i] == ID)
      fieldindex[i] = find_label("id",nwords,labels);
    else if (fieldtype[i] == TYPE)
      fieldindex[i] = find_label("type",nwords,labels);

    else if (fieldtype[i] == X) {
      fieldindex[i] = find_label("x",nwords,labels);
      xflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("xs",nwords,labels);
        u_index = find_label("xu",nwords,labels);
        su_index = find_label("xsu",nwords,labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          xflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          xflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          xflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Y) {
      fieldindex[i] = find_label("y",nwords,labels);
      yflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("ys",nwords,labels);
        u_index = find_label("yu",nwords,labels);
        su_index = find_label("ysu",nwords,labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          yflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          yflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          yflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Z) {
      fieldindex[i] = find_label("z",nwords,labels);
      zflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("zs",nwords,labels);
        u_index = find_label("zu",nwords,labels);
        su_index = find_label("zsu",nwords,labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          zflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          zflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          zflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == VX)
      fieldindex[i] = find_label("vx",nwords,labels);
    else if (fieldtype[i] == VY)
      fieldindex[i] = find_label("vy",nwords,labels);
    else if (fieldtype[i] == VZ)
      fieldindex[i] = find_label("vz",nwords,labels);

    else if (fieldtype[i] == FX)
      fieldindex[i] = find_label("fx",nwords,labels);
    else if (fieldtype[i] == FY)
      fieldindex[i] = find_label("fy",nwords,labels);
    else if (fieldtype[i] == FZ)
      fieldindex[i] = find_label("fz",nwords,labels);

    else if (fieldtype[i] == Q)
      fieldindex[i] = find_label("q",nwords,labels);

    else if (fieldtype[i] == IX)
      fieldindex[i] = find_label("ix",nwords,labels);
    else if (fieldtype[i] == IY)
      fieldindex[i] = find_label("iy",nwords,labels);
    else if (fieldtype[i] == IZ)
      fieldindex[i] = find_label("iz",nwords,labels);
  }

  delete [] labels;

  // set fieldflag = -1 if any unfound fields

  fieldflag = 0;
  for (int i = 0; i < nfield; i++)
    if (fieldindex[i] < 0) fieldflag = -1;

  // create internal vector of word ptrs for future parsing of per-atom lines

  words = new char*[nwords];

  return natoms;
}

/* ----------------------------------------------------------------------
   read N atom lines from dump file
   stores appropriate values in fields array
   return 0 if success, 1 if error
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::read_atoms(int n, int nfield, double **fields)
{
  int i,m;
  char *eof;

  for (i = 0; i < n; i++) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");

    // tokenize the line

    words[0] = strtok(line," \t\n\r\f");
    for (m = 1; m < nwords; m++)
      words[m] = strtok(NULL," \t\n\r\f");

    // convert selected fields to floats

    for (m = 0; m < nfield; m++)
      fields[i][m] = atof(words[fieldindex[m]]);
  }
}

/* ----------------------------------------------------------------------
   match label to any of N labels
   return index of match or -1 if no match
------------------------------------------------------------------------- */

int ReaderNative::find_label(const char *label, int n, char **labels)
{
  for (int i = 0; i < n; i++)
    if (strcmp(label,labels[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   read N lines from dump file
   only last one is saved in line
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderNative::read_lines(int n)
{
  char *eof = NULL;
  if (n <= 0) return;
  for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");
}
