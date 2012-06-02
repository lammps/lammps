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

#include "string.h"
#include "stdlib.h"
#include "read_dump_native.h"
#include "atom.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024        // max line length in dump file

enum{ID,TYPE,X,Y,Z,VX,VY,VZ,IX,IY,IZ};
enum{UNSET,UNSCALED,SCALED};

/* ---------------------------------------------------------------------- */

ReadDumpNative::ReadDumpNative(LAMMPS *lmp) : Pointers(lmp)
{
  dimension = domain->dimension;
  triclinic = domain->triclinic;

  line = new char[MAXLINE];
  words = NULL;
  fieldindex = NULL;
}

/* ---------------------------------------------------------------------- */

ReadDumpNative::~ReadDumpNative()
{
  delete [] line;
  delete [] words;
  delete [] fieldindex;
}

/* ---------------------------------------------------------------------- */

void ReadDumpNative::init(FILE *fpcaller)
{
  fp = fpcaller;
}

/* ----------------------------------------------------------------------
   proc 0 scans dump file until reaching snapshot with timestamp = Nstep
   extract natoms and box bounds from snapshot
   set fieldindex for specified fields and overall scaled setting
   error check on current vs new box and fields
   NOTE: error checking should maybe be moved to parent
------------------------------------------------------------------------- */

void ReadDumpNative::scan(bigint nstep, 
			  int nfield_caller, int *fieldtype, 
			  char **fieldlabel, int scaledflag,
			  bigint &natoms, double box[3][3], int &scaled)
{
  int nchunk,triclinic_snap,s_index,u_index,su_index;
  bigint ntimestep,nremain;
  char *bc,*names;

  nfield = nfield_caller;

  read_lines(1);

  while (1) {
    if (strstr(line,"ITEM: TIMESTEP") != line)
      error->one(FLERR,"Incorrectly formatted dump file");
    read_lines(1);
    sscanf(line,BIGINT_FORMAT,&ntimestep);

    if (ntimestep > nstep)
      error->one(FLERR,"Dump file does not contain requested snapshot");

    read_lines(2);
    sscanf(line,BIGINT_FORMAT,&natoms);

    // skip snapshot
    // invoke read_lines() in chunks no larger than MAXSMALLINT

    if (ntimestep < nstep) {
      read_lines(5);

      nremain = natoms;
      while (nremain) {
	nchunk = MIN(nremain,MAXSMALLINT);
	read_lines(nchunk);
	nremain -= nchunk;
      }

      read_lines(1);
    } else break;
  }

  // found correct snapshot
  // read box size and boundary conditions

  triclinic_snap = 0;
  box[0][2] = box[1][2] = box[2][2] = 0.0;
  read_lines(1);
  bc = &line[strlen("ITEM: BOX BOUNDS ")];
  if (bc[0] == 'x') {
    triclinic_snap = 1;
    bc = &bc[9];
  }
  
  char boundstr[9];
  domain->boundary_string(boundstr);
  if (strstr(bc,boundstr) != bc) 
    error->warning(FLERR,"Read_dump boundary flags do not match simulation");
  
  read_lines(1);
  if (!triclinic_snap) sscanf(line,"%lg %lg",&box[0][0],&box[0][1]);
  else sscanf(line,"%lg %lg %lg",&box[0][0],&box[0][1],&box[0][2]);
  read_lines(1);
  if (!triclinic_snap) sscanf(line,"%lg %lg",&box[1][0],&box[1][1]);
  else sscanf(line,"%lg %lg %lg",&box[1][0],&box[1][1],&box[1][2]);
  read_lines(1);
  if (!triclinic_snap) sscanf(line,"%lg %lg",&box[2][0],&box[2][1]);
  else sscanf(line,"%lg %lg %lg",&box[2][0],&box[2][1],&box[2][2]);
  
  // read ITEM: ATOMS line
  // labels = column labels
  
  read_lines(1);
  names = &line[strlen("ITEM: ATOMS ")];
  
  nwords = atom->count_words(names);
  char **labels = new char*[nwords];
  labels[0] = strtok(names," \t\n\r\f");
  if (labels[0] == NULL) 
    error->one(FLERR,"Incorrect atom format in dump file");
  for (int m = 1; m < nwords; m++) {
    labels[m] = strtok(NULL," \t\n\r\f");
    if (labels[m] == NULL) 
      error->one(FLERR,"Incorrect atom format in dump file");
  }
  
  // match each field with column
  
  fieldindex = new int[nfield];
  int xflag = UNSET;
  int yflag = UNSET;
  int zflag = UNSET;
  
  for (int i = 0; i < nfield; i++) {
    if (fieldlabel[i]) {
      fieldindex[i] = find_label(fieldlabel[i],nwords,labels);
      if (fieldtype[i] == X) xflag = scaledflag;
      else if (fieldtype[i] == Y) yflag = scaledflag;
      else if (fieldtype[i] == Z) zflag = scaledflag;
    }
    
    else if (fieldtype[i] == ID)
      fieldindex[i] = find_label("id",nwords,labels);
    else if (fieldtype[i] == TYPE)
      fieldindex[i] = find_label("type",nwords,labels);
    
    else if (fieldtype[i] == X) {
      fieldindex[i] = find_label("x",nwords,labels);
      xflag = UNSCALED;
      if (fieldindex[i] < 0) {
	fieldindex[i] = nwords;
	s_index = find_label("xs",nwords,labels);
	u_index = find_label("xu",nwords,labels);
	su_index = find_label("xsu",nwords,labels);
	if (s_index >= 0 && s_index < fieldindex[i]) {
	  fieldindex[i] = s_index;
	  xflag = SCALED;
	}
	if (u_index >= 0 && u_index < fieldindex[i]) {
	  fieldindex[i] = u_index;
	  xflag = UNSCALED;
	}
	if (su_index >= 0 && su_index < fieldindex[i]) {
	  fieldindex[i] = su_index;
	  xflag = SCALED;
	}
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;
      
    } else if (fieldtype[i] == Y) {
      fieldindex[i] = find_label("y",nwords,labels);
      yflag = UNSCALED;
      if (fieldindex[i] < 0) {
	fieldindex[i] = nwords;
	s_index = find_label("ys",nwords,labels);
	u_index = find_label("yu",nwords,labels);
	su_index = find_label("ysu",nwords,labels);
	if (s_index >= 0 && s_index < fieldindex[i]) {
	  fieldindex[i] = s_index;
	  yflag = SCALED;
	}
	if (u_index >= 0 && u_index < fieldindex[i]) {
	  fieldindex[i] = u_index;
	  yflag = UNSCALED;
	}
	if (su_index >= 0 && su_index < fieldindex[i]) {
	  fieldindex[i] = su_index;
	  yflag = SCALED;
	}
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;
      
    } else if (fieldtype[i] == Z) {
      fieldindex[i] = find_label("z",nwords,labels);
      zflag = UNSCALED;
      if (fieldindex[i] < 0) {
	fieldindex[i] = nwords;
	s_index = find_label("zs",nwords,labels);
	u_index = find_label("zu",nwords,labels);
	su_index = find_label("zsu",nwords,labels);
	if (s_index >= 0 && s_index < fieldindex[i]) {
	  fieldindex[i] = s_index;
	  zflag = SCALED;
	}
	if (u_index >= 0 && u_index < fieldindex[i]) {
	  fieldindex[i] = u_index;
	  zflag = UNSCALED;
	}
	if (su_index >= 0 && su_index < fieldindex[i]) {
	  fieldindex[i] = su_index;
	  zflag = SCALED;
	}
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;
      
    } else if (fieldtype[i] == VX)
      fieldindex[i] = find_label("vx",nwords,labels);
    else if (fieldtype[i] == VY)
      fieldindex[i] = find_label("vy",nwords,labels);
    else if (fieldtype[i] == VZ)
      fieldindex[i] = find_label("vz",nwords,labels);
    else if (fieldtype[i] == IX)
      fieldindex[i] = find_label("ix",nwords,labels);
    else if (fieldtype[i] == IY)
      fieldindex[i] = find_label("iy",nwords,labels);
    else if (fieldtype[i] == IZ)
      fieldindex[i] = find_label("iz",nwords,labels);
  }
  
  // error checks
  
  if ((triclinic_snap && !triclinic) || 
      (!triclinic_snap && triclinic))
    error->one(FLERR,"Read_dump triclinic setting does not match simulation");
  
  for (int i = 0; i < nfield; i++)
    if (fieldindex[i] < 0) 
      error->one(FLERR,"Read_dump field not found in dump file");
  
  // set overall scaling of coordinates
  // error if x,y,z scaling is not the same
  
  scaled = MAX(xflag,yflag);
  scaled = MAX(zflag,scaled);
  if ((xflag != UNSET && xflag != scaled) ||
      (yflag != UNSET && yflag != scaled) ||
      (zflag != UNSET && zflag != scaled))
    error->one(FLERR,"Read_dump x,y,z fields do not have consistent scaling");
  
  // scaled, triclinic coords require all 3 x,y,z fields to perform unscaling
  
  if (scaled == SCALED && triclinic) {
    int flag = 0;
    if (xflag != scaled) flag = 1;
    if (yflag != scaled) flag = 1;
    if (dimension == 3 && zflag != scaled) flag = 1;
    if (flag)
      error->one(FLERR,"All read_dump x,y,z fields must be specified for "
		 "scaled, triclinic coords");
  }
  
  delete [] labels;
  
  // create vector of word ptrs for future parsing of per-atom lines
  
  words = new char*[nwords];
}

/* ----------------------------------------------------------------------
   proc 0 reads N atom lines from dump file
   stores appropriate values in fields array
------------------------------------------------------------------------- */

void ReadDumpNative::read(int n, double **fields)
{
  int i,m;
  char *eof;

  for (int i = 0; i < n; i++) {
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

int ReadDumpNative::find_label(const char *label, int n, char **labels)
{
  for (int i = 0; i < n; i++)
    if (strcmp(label,labels[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   only last one is saved in line
------------------------------------------------------------------------- */

void ReadDumpNative::read_lines(int n)
{
  char *eof;
  for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");
}
