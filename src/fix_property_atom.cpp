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

#include "stdlib.h"
#include "string.h"
#include "fix_property_atom.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{MOLECULE,INTEGER,DOUBLE};

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix property/atom command");

  restart_peratom = 1;

  int iarg = 3;
  nvalue = narg-iarg;
  style = new int[nvalue];
  index = new int[nvalue];

  molecule_flag = 0;

  nvalue = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mol") == 0) {
      if (atom->molecule_flag)
        error->all(FLERR,"Fix property/atom mol when atom_style "
                   "already has molecule attribute");
      if (molecule_flag) 
        error->all(FLERR,"Fix property/atom cannot specify mol twice");
      style[nvalue] = MOLECULE;
      atom->molecule_flag = molecule_flag = 1;
      nvalue++;
    } else if (strstr(arg[iarg],"i_") == arg[iarg]) {
      style[nvalue] = INTEGER;
      int tmp;
      index[nvalue] = atom->find_custom(&arg[iarg][2],tmp);
      if (index[nvalue] >= 0) 
        error->all(FLERR,"Fix property/atom vector name already exists");
      index[nvalue] = atom->add_custom(&arg[iarg][2],0);
      nvalue++;
    } else if (strstr(arg[iarg],"d_") == arg[iarg]) {
      style[nvalue] = DOUBLE;
      int tmp;
      index[nvalue] = atom->find_custom(&arg[iarg][2],tmp);
      if (index[nvalue] >= 0) 
        error->all(FLERR,"Fix property/atom vector name already exists");
      index[nvalue] = atom->add_custom(&arg[iarg][2],1);
      nvalue++;
    } else break;

    iarg++;
  }

  // optional args

  border = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ghost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix property/atom command");
      if (strcmp(arg[iarg+1],"no") == 0) border = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) border = 1;
      else error->all(FLERR,"Illegal fix property/atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix property/atom command");
  }

  if (border) comm_border = nvalue;

  // perform initial allocation of atom-based array
  // register with Atom class

  nmax_old = 0;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);
  if (border) atom->add_callback(2);
}

/* ---------------------------------------------------------------------- */

FixPropertyAtom::~FixPropertyAtom()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);
  if (border) atom->delete_callback(id,2);

  // deallocate per-atom vectors in Atom class
  // set ptrs to NULL, so they no longer exist for Atom class

  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) {
      atom->molecule_flag = 0;
      memory->destroy(atom->molecule);
      atom->molecule = NULL;
    } else if (style[m] == INTEGER) {
      atom->remove_custom(0,index[m]);
    } else if (style[m] == DOUBLE) {
      atom->remove_custom(1,index[m]);
    }
  }

  delete [] style;
  delete [] index;
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtom::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtom::init()
{
  // check again if atom style already defines molecule vector
  //   based on molecular setting
  // this could happen if script changed atom_style
  //   after this fix was already in place

  if (molecule_flag && atom->molecular)
    error->all(FLERR,"Fix property/atom mol when atom_style "
               "already has molecule attribute");
}

/* ----------------------------------------------------------------------
   unpack section of data file
------------------------------------------------------------------------- */

void FixPropertyAtom::read_data_section(char *keyword, int n, char *buf)
{
  int j,m,tagdata;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = atom->count_words(buf);
  *next = '\n';

  if (nwords != nvalue+1) {
    char str[128];
    sprintf(str,"Incorrect %s format in data file",keyword);
    error->all(FLERR,str);
  }

  char **values = new char*[nwords];

  // loop over lines of atom velocities
  // tokenize the line into values
  // if I own atom tag, unpack its values

  int map_tag_max = atom->map_tag_max;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    tagdata = atoi(values[0]);
    if (tagdata <= 0 || tagdata > map_tag_max) {
      char str[128];
      sprintf(str,"Invalid atom ID in %s section of data file",keyword);
      error->one(FLERR,str);
    }

    // assign words in line to per-atom vectors

    if ((m = atom->map(tagdata)) >= 0) {
      for (j = 0; j < nvalue; j++) {
        if (style[j] == MOLECULE) atom->molecule[m] = atoi(values[j+1]);
        else if (style[j] == INTEGER)
          atom->ivector[index[j]][m] = atoi(values[j+1]);
        else if (style[j] == DOUBLE) 
          atom->dvector[index[j]][m] = atof(values[j+1]);
      }
    }

    buf = next + 1;
  }

  delete [] values;
}

/* ---------------------------------------------------------------------- */

bigint FixPropertyAtom::read_data_skip_lines(char *keyword)
{
  return atom->natoms;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyAtom::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) bytes = atom->nmax * sizeof(int);
    else if (style[m] == INTEGER) bytes = atom->nmax * sizeof(int);
    else if (style[m] == DOUBLE) bytes = atom->nmax * sizeof(double);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
   initialize new values to 0,
   since AtomVec class won't do it as atoms are added, 
   e.g. in create_atom() or data_atom()
------------------------------------------------------------------------- */

void FixPropertyAtom::grow_arrays(int nmax)
{
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) {
      memory->grow(atom->molecule,nmax,"atom:molecule");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->molecule[nmax_old],0,nbytes);
    } else if (style[m] == INTEGER) {
      memory->grow(atom->ivector[index[m]],nmax,"atom:ivector");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[m]][nmax_old],0,nbytes);
    } else if (style[m] == DOUBLE) {
      memory->grow(atom->dvector[index[m]],nmax,"atom:dvector");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->dvector[index[m]][nmax_old],0,nbytes);
    }
  }

  nmax_old = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtom::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE)
      atom->molecule[j] = atom->molecule[i];
    else if (style[m] == INTEGER)
      atom->ivector[index[m]][j] = atom->ivector[index[m]][i];
    else if (style[m] == DOUBLE)
      atom->dvector[index[m]][j] = atom->dvector[index[m]][i];
  }
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_border(int n, int *list, double *buf)
{
  int i,j,k;

  int m = 0;
  for (k = 0; k < nvalue; k++) {
    if (style[k] == MOLECULE) {
      int *molecule = atom->molecule;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = molecule[j];
      }
    } else if (style[j] == INTEGER) {
      int *ivector = atom->ivector[index[k]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ivector[j];
      }
    } else if (style[j] == DOUBLE) {
      double *dvector = atom->dvector[index[k]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = dvector[j];
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixPropertyAtom::unpack_border(int n, int first, double *buf)
{
  int i,k,last;

  int m = 0;
  for (k = 0; k < nvalue; k++) {
    if (style[k] == MOLECULE) {
      int *molecule = atom->molecule;
      last = first + n;
      for (i = first; i < last; i++)
        molecule[i] = static_cast<int> (buf[m++]);
    } else if (style[k] == INTEGER) {
      int *ivector = atom->ivector[index[k]];
      last = first + n;
      for (i = first; i < last; i++)
        ivector[i] = static_cast<int> (buf[m++]);
    } else if (style[k] == DOUBLE) {
      double *dvector = atom->dvector[index[k]];
      last = first + n;
      for (i = first; i < last; i++)
        dvector[i] = buf[m++];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) buf[m] = atom->molecule[i];
    else if (style[m] == INTEGER) buf[m] = atom->ivector[index[m]][i];
    else if (style[m] == DOUBLE) buf[m] = atom->dvector[index[m]][i];
  }
  return nvalue;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) 
      atom->molecule[nlocal] = static_cast<int> (buf[m]);
    else if (style[m] == INTEGER) 
      atom->ivector[index[m]][nlocal] = static_cast<int> (buf[m]);
    else if (style[m] == DOUBLE)
      atom->dvector[index[m]][nlocal] = buf[m];
  }
  return nvalue;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_restart(int i, double *buf)
{
  buf[0] = nvalue+1;
  for (int m = 1; m <= nvalue; m++) {
    if (style[m] == MOLECULE) buf[m] = atom->molecule[i];
    else if (style[m] == INTEGER) buf[m] = atom->ivector[index[m]][i];
    else if (style[m] == DOUBLE) buf[m] = atom->dvector[index[m]][i];
  }
  return nvalue+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyAtom::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int i = 0; i < nvalue; i++) {
    if (style[i] == MOLECULE) 
      atom->molecule[nlocal] = static_cast<int> (extra[nlocal][m++]);
    else if (style[m] == INTEGER) 
      atom->ivector[index[m]][nlocal] = static_cast<int> (extra[nlocal][m++]);
    else if (style[m] == DOUBLE)
      atom->dvector[index[m]][nlocal] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::maxsize_restart()
{
  return nvalue+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::size_restart(int nlocal)
{
  return nvalue+1;
}
