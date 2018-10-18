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

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "reader_molfile.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "molfile_interface.h"
#include "math_const.h"

using namespace LAMMPS_NS;
typedef MolfileInterface MFI;
using namespace MathConst;

enum{ID,TYPE,X,Y,Z,VX,VY,VZ};
#define SMALL 1.0e-6

// true if the difference between two floats is "small".
// cannot use fabsf() since it is not fully portable.
static bool is_smalldiff(const float &val1, const float &val2)
{
  return (fabs(static_cast<double>(val1-val2)) < SMALL);
}

/* ---------------------------------------------------------------------- */

ReaderMolfile::ReaderMolfile(LAMMPS *lmp) : Reader(lmp)
{
  mf = NULL;
  coords = NULL;
  vels = NULL;
  types = NULL;
  fieldindex = NULL;
  nstep = 0;
  needvels = 0;
  me = comm->me;
}

/* ---------------------------------------------------------------------- */

ReaderMolfile::~ReaderMolfile()
{
  if (me == 0) {
    memory->destroy(fieldindex);

    memory->destroy(types);
    memory->destroy(coords);
    memory->destroy(vels);
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
    snprintf(str,1024,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  if (natoms < 1) {
    snprintf(str,1024,"No atoms in file %s",file);
    error->one(FLERR,str);
  }

  memory->create(types,natoms,"reader:types");
  memory->create(coords,3*natoms,"reader:coords");
  if (mf->has_vels())
    memory->create(vels,3*natoms,"reader:vels");

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

  // try to read in the time step (coordinates, velocities and cell only)
  rv = mf->timestep(coords, vels, cell, NULL);
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
     xyz flags = from input scaleflag & wrapflag
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
                                  int scaleflag, int wrapflag, int &fieldflag,
                                  int &xflag, int &yflag, int &zflag)
{
  nid = 0;

  // signal that we have no box info at all so far.
  triclinic = -1;

  // heuristics to determine if we have boxinfo (first if)
  // and whether we have an orthogonal box (second if)
  if (!is_smalldiff(cell[0]*cell[1]*cell[2], 0.0f)) {
    if (is_smalldiff(cell[3],90.0f) && is_smalldiff(cell[4],90.0f) &&
        is_smalldiff(cell[5],90.0f)) {
      triclinic = 0;
      // we have no information about the absolute location
      // of the box, so we assume that the origin is in the middle.
      // also we cannot tell periodicity. we assume, yes.
      box[0][0] = -0.5*static_cast<double>(cell[0]);
      box[0][1] =  0.5*static_cast<double>(cell[0]);
      box[0][2] =  0.0;
      box[1][0] = -0.5*static_cast<double>(cell[1]);
      box[1][1] =  0.5*static_cast<double>(cell[1]);
      box[1][2] =  0.0;
      box[2][0] = -0.5*static_cast<double>(cell[2]);
      box[2][1] =  0.5*static_cast<double>(cell[2]);
      box[2][2] =  0.0;
    } else {

      triclinic = 1;

      const double la = static_cast<double>(cell[0]);
      const double lb = static_cast<double>(cell[1]);
      const double lc = static_cast<double>(cell[2]);
      const double alpha = static_cast<double>(cell[3]);
      const double beta  = static_cast<double>(cell[4]);
      const double gamma = static_cast<double>(cell[5]);

      const double lx = la;
      const double xy = lb * cos(gamma/90.0*MY_PI2);
      const double xz = lc * cos(beta/90.0*MY_PI2);
      const double ly = sqrt(lb*lb - xy*xy);
      const double yz = (fabs(ly) > SMALL) ?
        (lb*lc*cos(alpha/90.0*MY_PI2) - xy*xz) / ly : 0.0;
      const double lz = sqrt(lc*lc - xz*xz - yz*yz);

      /* go from box length to boundary */
      double xbnd;

      xbnd = 0.0;
      xbnd = (xy < xbnd) ? xy : xbnd;
      xbnd = (xz < xbnd) ? xz : xbnd;
      xbnd = (xy+xz < xbnd) ? (xy + xz) : xbnd;
      box[0][0] = -0.5*lx + xbnd;

      xbnd = 0.0;
      xbnd = (xy > xbnd) ? xy : xbnd;
      xbnd = (xz > xbnd) ? xz : xbnd;
      xbnd = (xy+xz > xbnd) ? (xy + xz) : xbnd;
      box[0][1] =  0.5*lx+xbnd;
      box[0][2] =  xy;

      xbnd = 0.0;
      xbnd = (yz < xbnd) ? yz : xbnd;
      box[1][0] = -0.5*ly+xbnd;

      xbnd = 0.0;
      xbnd = (yz > xbnd) ? yz : xbnd;
      box[1][1] =  0.5*ly+xbnd;
      box[1][2] =  xz;

      box[2][0] = -0.5*lz;
      box[2][1] =  0.5*lz;
      box[2][2] =  yz;
    }
  }

  // if no field info requested, just return
  if (!fieldinfo) return natoms;

  memory->create(fieldindex,nfield,"read_dump:fieldindex");

  // we know nothing about the style of coordinates,
  // so caller has to set the proper flags

  xflag = 2*scaleflag + wrapflag + 1;
  yflag = 2*scaleflag + wrapflag + 1;
  zflag = 2*scaleflag + wrapflag + 1;

  // copy fieldtype list for supported fields

  fieldflag = 0;
  needvels = 0;
  for (int i = 0; i < nfield; i++) {
    if ( (fieldtype[i] == X) ||
         (fieldtype[i] == Y) ||
         (fieldtype[i] == Z) ||
         (fieldtype[i] == ID) ||
         (fieldtype[i] == TYPE) ) {
      fieldindex[i] = fieldtype[i];
    } else if ( (fieldtype[i] == VX) ||
                (fieldtype[i] == VY) ||
                (fieldtype[i] == VZ) ) {
      fieldindex[i] = fieldtype[i];
      needvels = 1;
    } else {
      fieldflag = 1;
    }
  }

  if ((needvels > 0) && (!mf->has_vels()))
    error->one(FLERR,"Molfile plugin does not support reading velocities");

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
  int i,m,mytype;
  char buf[16];

  for (i = 0; i < n; i++) {
    ++nid;

    if (mf->property(MFI::P_TYPE,nid-1,buf) != MFI::P_NONE) {
      mytype = atoi(buf);
    } else mytype = 0;

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
      case VX:
        fields[i][m] = static_cast<double>(vels[3*nid-3]);
        break;
      case VY:
        fields[i][m] = static_cast<double>(vels[3*nid-2]);
        break;
      case VZ:
        fields[i][m] = static_cast<double>(vels[3*nid-1]);
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
