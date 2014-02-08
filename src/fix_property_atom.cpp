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

enum{MOLECULE,CHARGE,INTEGER,DOUBLE};

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix property/atom command");

  restart_peratom = 1;
  wd_section = 1;

  int iarg = 3;
  nvalue = narg-iarg;
  style = new int[nvalue];
  index = new int[nvalue];

  molecule_flag = 0;
  q_flag = 0;

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
    } else if (strcmp(arg[iarg],"q") == 0) {
      if (atom->q_flag)
        error->all(FLERR,"Fix property/atom q when atom_style "
                   "already has charge attribute");
      if (q_flag) 
        error->all(FLERR,"Fix property/atom cannot specify q twice");
      style[nvalue] = CHARGE;
      atom->q_flag = q_flag = 1;
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

  // store current atom style

  int n = strlen(atom->atom_style) + 1;
  astyle = new char[n];
  strcpy(astyle,atom->atom_style);

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
    } else if (style[m] == CHARGE) {
      atom->q_flag = 0;
      memory->destroy(atom->q);
      atom->q = NULL;
    } else if (style[m] == INTEGER) {
      atom->remove_custom(0,index[m]);
    } else if (style[m] == DOUBLE) {
      atom->remove_custom(1,index[m]);
    }
  }

  delete [] style;
  delete [] index;
  delete [] astyle;
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
  // error if atom style has changed since fix was defined
  // don't allow this b/c user could change to style that defines molecule,q

  if (strcmp(astyle,atom->atom_style) != 0)
    error->all(FLERR,"Atom style was redefined after using fix property/atom");
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
------------------------------------------------------------------------- */

void FixPropertyAtom::read_data_section(char *keyword, int n, char *buf)
{
  int j,m;
  tagint itag;
  char *next;

  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

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

  // loop over lines of atom info
  // tokenize the line into values
  // if I own atom tag, unpack its values

  tagint map_tag_max = atom->map_tag_max;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    for (j = 1; j < nwords; j++)
      values[j] = strtok(NULL," \t\n\r\f");

    itag = ATOTAGINT(values[0]);
    if (itag <= 0 || itag > map_tag_max) {
      char str[128];
      sprintf(str,"Invalid atom ID in %s section of data file",keyword);
      error->one(FLERR,str);
    }

    // assign words in line to per-atom vectors

    if ((m = atom->map(itag)) >= 0) {
      for (j = 0; j < nvalue; j++) {
        if (style[j] == MOLECULE) atom->molecule[m] = ATOTAGINT(values[j+1]);
        else if (style[j] == CHARGE) atom->q[m] = atof(values[j+1]);
        else if (style[j] == INTEGER)
          atom->ivector[index[j]][m] = atoi(values[j+1]);
        else if (style[j] == DOUBLE) 
          atom->dvector[index[j]][m] = atof(values[j+1]);
      }
    }

    buf = next + 1;
  }

  delete [] values;

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }
}

/* ----------------------------------------------------------------------
   return # of lines in section of data file labeled by keyword
------------------------------------------------------------------------- */

bigint FixPropertyAtom::read_data_skip_lines(char *keyword)
{
  return atom->natoms;
}

/* ----------------------------------------------------------------------
   return size I own for Mth data section
   # of data sections = 1 for this fix
   nx = # of local atoms
   ny = columns = tag + nvalues
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_size(int mth, int &nx, int &ny)
{
  nx = atom->nlocal;
  ny = nvalue + 1;
}

/* ----------------------------------------------------------------------
   pack values for Mth data section into buf
   buf allocated by caller as Nlocal by Nvalues+1
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_pack(int mth, double **buf)
{
  int i;

  // 1st column = atom tag
  // rest of columns = per-atom values

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) buf[i][0] = ubuf(tag[i]).d;

  for (int m = 0; m < nvalue; m++) {
    int mp1 = m+1;
    if (style[m] == MOLECULE) {
      tagint *molecule = atom->molecule;
      for (i = 0; i < nlocal; i++) buf[i][mp1] = ubuf(molecule[i]).d;
    } else if (style[m] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < nlocal; i++) buf[i][mp1] = q[i];
    } else if (style[m] == INTEGER) {
      int *ivec = atom->ivector[index[m]];
      for (i = 0; i < nlocal; i++) buf[i][mp1] = ubuf(ivec[i]).d;
    } else if (style[m] == DOUBLE) {
      double *dvec = atom->dvector[index[m]];
      for (i = 0; i < nlocal; i++) buf[i][mp1] = dvec[i];
    }
  }
}

/* ----------------------------------------------------------------------
   write section keyword for Mth data section to file
   use Molecules or Charges if that is only field, else use fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_keyword(int mth, FILE *fp)
{
  if (nvalue == 1 && style[0] == MOLECULE) fprintf(fp,"\nMolecules\n\n");
  else if (nvalue == 1 && style[0] == CHARGE) fprintf(fp,"\nCharges\n\n");
  else fprintf(fp,"\n%s\n\n",id);
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section(int mth, FILE *fp, 
                                         int n, double **buf, int index)
{
  int m;

  for (int i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT,(tagint) ubuf(buf[i][0]).i);
    for (m = 0; m < nvalue; m++) {
      if (style[m] == MOLECULE)
        fprintf(fp," " TAGINT_FORMAT,(tagint) ubuf(buf[i][m+1]).i);
      else if (style[m] == INTEGER)
        fprintf(fp," %d",(int) ubuf(buf[i][m+1]).i);
      else fprintf(fp," %g",buf[i][m+1]);
    }
    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyAtom::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) bytes = atom->nmax * sizeof(tagint);
    else if (style[m] == CHARGE) bytes = atom->nmax * sizeof(double);
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
      size_t nbytes = (nmax-nmax_old) * sizeof(tagint);
      memset(&atom->molecule[nmax_old],0,nbytes);
    } else if (style[m] == CHARGE) {
      memory->grow(atom->q,nmax,"atom:q");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->q[nmax_old],0,nbytes);
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
    else if (style[m] == CHARGE)
      atom->q[j] = atom->q[i];
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
      tagint *molecule = atom->molecule;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(molecule[j]).d;
      }
    } else if (style[k] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = q[j];
      }
    } else if (style[k] == INTEGER) {
      int *ivector = atom->ivector[index[k]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(ivector[j]).i;
      }
    } else if (style[k] == DOUBLE) {
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
      tagint *molecule = atom->molecule;
      last = first + n;
      for (i = first; i < last; i++)
        molecule[i] = (tagint) ubuf(buf[m++]).i;
    } else if (style[k] == CHARGE) {
      double *q = atom->q;
      last = first + n;
      for (i = first; i < last; i++)
        q[i] = buf[m++];
    } else if (style[k] == INTEGER) {
      int *ivector = atom->ivector[index[k]];
      last = first + n;
      for (i = first; i < last; i++)
        ivector[i] = (int) ubuf(buf[m++]).i;
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
    if (style[m] == MOLECULE) buf[m] = ubuf(atom->molecule[i]).d;
    else if (style[m] == CHARGE) buf[m] = atom->q[i];
    else if (style[m] == INTEGER) buf[m] = ubuf(atom->ivector[index[m]][i]).d;
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
      atom->molecule[nlocal] = (tagint) ubuf(buf[m]).i;
    else if (style[m] == CHARGE)
      atom->q[nlocal] = buf[m];
    else if (style[m] == INTEGER) 
      atom->ivector[index[m]][nlocal] = (int) ubuf(buf[m]).i;
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

  int m = 1;
  for (int j = 0; j < nvalue; j++) {
    if (style[j] == MOLECULE) buf[m++] = ubuf(atom->molecule[i]).d;
    else if (style[j] == CHARGE) buf[m++] = atom->q[i];
    else if (style[j] == INTEGER) buf[m++] = ubuf(atom->ivector[index[j]][i]).d;
    else if (style[j] == DOUBLE) buf[m++] = atom->dvector[index[j]][i];
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
      atom->molecule[nlocal] = (tagint) ubuf(extra[nlocal][m++]).i;
    else if (style[i] == CHARGE)
      atom->q[nlocal] = extra[nlocal][m++];
    else if (style[i] == INTEGER) 
      atom->ivector[index[i]][nlocal] = (int) ubuf(extra[nlocal][m++]).i;
    else if (style[i] == DOUBLE)
      atom->dvector[index[i]][nlocal] = extra[nlocal][m++];
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
