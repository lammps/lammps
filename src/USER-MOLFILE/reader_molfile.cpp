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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "reader_molfile.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "molfile_interface.h"

using namespace LAMMPS_NS;
typedef MolfileInterface MFI;

enum{ID,TYPE,X,Y,Z};
#define SMALL 1.0e-6

/* ---------------------------------------------------------------------- */

ReaderMolfile::ReaderMolfile(LAMMPS *lmp) : Reader(lmp)
{
  mf = NULL;
  coords = NULL;
  types = NULL;
  fieldindex = NULL;
  nstep = 0;
  me = comm->me;
}

/* ---------------------------------------------------------------------- */

ReaderMolfile::~ReaderMolfile()
{
  if (me == 0) {
    memory->destroy(fieldindex);

    memory->destroy(types);
    memory->destroy(coords);
    if (mf) delete mf;
  }
}

/* ----------------------------------------------------------------------
   pass on settings to find and load the proper plugin
   called by all processors.
------------------------------------------------------------------------- */
void ReaderMolfile::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal molfile reader command");

  if (me == 0) {
    mf = new MolfileInterface(arg[0],MFI::M_READ);

    const char *path = (const char *) ".";
    if (narg > 1)
      path=arg[1];

    if (mf->find_plugin(path)!= MFI::E_MATCH)
      error->one(FLERR,"No suitable molfile plugin found");

    if (screen)
      fprintf(screen,"Dump reader uses molfile plugin: %s\n",
              mf->get_plugin_name());

    if (logfile)
      fprintf(logfile,"Dump reader uses molfile plugin: %s\n",
              mf->get_plugin_name());
  }
}

/* ----------------------------------------------------------------------
   try to open given file through plugin interface
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderMolfile::open_file(const char *file)
{
  int rv;
  char str[1024];

  // close open file, if needed.
  if (mf->is_open()) mf->close();

  rv = mf->open(file,&natoms);

  if (rv != MFI::E_NONE) {
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  if (natoms < 1) {
    sprintf(str,"No atoms in file %s",file);
    error->one(FLERR,str);
  }

  memory->create(types,natoms,"reader:types");
  memory->create(coords,3*natoms,"reader:coords");

  // initialize system properties, if available
  if (mf->has_props()) {
    mf->structure();
    mf->property(MFI::P_TYPE,types);

  } else {
    for (int i=0; i < natoms; ++i)
      types[i] = 1;
  }
}

/* ----------------------------------------------------------------------
   close current file
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderMolfile::close_file()
{
  mf->close();
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   only called by proc 0
------------------------------------------------------------------------- */

int ReaderMolfile::read_time(bigint &ntimestep)
{
  int rv;

  // try to read in the time step (coordinates and cell only)
  rv = mf->timestep(coords, NULL, cell, NULL);
  if (rv != 0) return 1;

  // we fake time step numbers.
  ntimestep = nstep;
  nstep++;

  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot from timestamp onward
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderMolfile::skip()
{
  // since we can only signal EOF to the caller in ::read_time(), we
  // have to read the entire timestep always there and this is a NOP.
  ;
}

/* ----------------------------------------------------------------------
   read remaining header info:
     return natoms
     box bounds, triclinic (inferred), fieldflag (1 if any fields not found),
     xyz flag = UNSET (not a requested field), SCALED, UNSCALED
   if fieldflag set:
     match Nfield fields to per-atom column labels
     allocate and set fieldindex = which column each field maps to
     fieldtype = X,VX,IZ etc
     fieldlabel = user-specified label or NULL if use fieldtype default
   xyz flag = scaledflag if has fieldlabel name, else set by x,xs,xu,xsu
   only called by proc 0
------------------------------------------------------------------------- */

bigint ReaderMolfile::read_header(double box[3][3], int &triclinic,
                                  int fieldinfo, int nfield,
                                  int *fieldtype, char **fieldlabel,
                                  int scaledflag, int &fieldflag,
                                  int &xflag, int &yflag, int &zflag)
{
  nid = 0;

  // signal that we have no box info at all so far.
  triclinic = -1;

  // heuristics to determine if we have boxinfo and whether
  // we have an orthogonal box.
  // XXX: add some tolerances for rounding and some sanity checks.
  if (fabs(static_cast<double>(cell[0]*cell[1]*cell[2])) > SMALL) {
    if ((cell[3] != 90.0f) || (cell[4] != 90.0f) ||
        (cell[5] != 90.0f)) triclinic = 1;
    else triclinic = 0;
  }

  // if no field info requested, just return
  if (!fieldinfo) return natoms;

  memory->create(fieldindex,nfield,"read_dump:fieldindex");

  // we know nothing about the scaling style of coordinates, 
  // so the caller has to set the proper flag.
  xflag = scaledflag;
  yflag = scaledflag;
  zflag = scaledflag;

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

void ReaderMolfile::read_atoms(int n, int nfield, double **fields)
{
  int i,m, mytype;

  for (i = 0; i < n; i++) {

    // FIXME: we need to find a way to signal parser errors here. XXX
    mytype = 1; // atoi(property())
    if (mytype < 1) mytype = 1;
    ++nid;

    for (m = 0; m < nfield; m++) {
      switch (fieldindex[m]) {
      case X:
        fields[i][m] = static_cast<double>(coords[3*nid-3]);
        break;
      case Y:
        fields[i][m] = static_cast<double>(coords[3*nid-2]);
        break;
      case Z:
        fields[i][m] = static_cast<double>(coords[3*nid-1]);
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
