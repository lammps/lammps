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

#include "fix_property_atom.h"
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "fmt/format.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{MOLECULE,CHARGE,RMASS,IVEC,DVEC,IARRAY,DARRAY};

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvalue(0), style(NULL), index(NULL), astyle(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix property/atom command");

  restart_peratom = 1;
  wd_section = 1;

  int iarg = 3;
  nvalue = narg-iarg;
  style = new int[nvalue];
  cols = new int[nvalue];
  index = new int[nvalue];

  molecule_flag = 0;
  q_flag = 0;
  rmass_flag = 0;

  nvalue = 0;
  values_peratom = 0;
  
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mol") == 0) {
      if (atom->molecule_flag)
        error->all(FLERR,"Fix property/atom mol when atom_style "
                   "already has molecule attribute");
      if (molecule_flag)
        error->all(FLERR,"Fix property/atom cannot specify mol twice");
      style[nvalue] = MOLECULE;
      cols[nvalue] = 0;
      atom->molecule_flag = molecule_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg],"q") == 0) {
      if (atom->q_flag)
        error->all(FLERR,"Fix property/atom q when atom_style "
                   "already has charge attribute");
      if (q_flag)
        error->all(FLERR,"Fix property/atom cannot specify q twice");
      style[nvalue] = CHARGE;
      cols[nvalue] = 0;
      atom->q_flag = q_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg],"rmass") == 0) {
      if (atom->rmass_flag)
        error->all(FLERR,"Fix property/atom rmass when atom_style "
                   "already has rmass attribute");
      if (rmass_flag)
        error->all(FLERR,"Fix property/atom cannot specify rmass twice");
      style[nvalue] = RMASS;
      cols[nvalue] = 0;
      atom->rmass_flag = rmass_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strstr(arg[iarg],"i_") == arg[iarg] ||
	       strstr(arg[iarg],"d_") == arg[iarg]) {
      int which = 0;
      if (arg[iarg][0] == 'd') which = 1;
      if (which == 0) style[nvalue] = IVEC;
      else style[nvalue] = DVEC;
      int tmp1,tmp2;
      index[nvalue] = atom->find_custom(&arg[iarg][2],tmp1,tmp2);
      if (index[nvalue] >= 0)
        error->all(FLERR,"Fix property/atom vector name already exists");
      cols[nvalue] = 0;
      index[nvalue] = atom->add_custom(&arg[iarg][2],which,cols[nvalue]);
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strstr(arg[iarg],"i2_") == arg[iarg] || 
	       strstr(arg[iarg],"d2_") == arg[iarg]) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix property/atom command");
      int which = 0;
      if (arg[iarg][0] == 'd') which = 1;
      if (which == 0) style[nvalue] = IARRAY;
      else style[nvalue] = DARRAY;
      int tmp1,tmp2;
      index[nvalue] = atom->find_custom(&arg[iarg][3],tmp1,tmp2);
      if (index[nvalue] >= 0)
        error->all(FLERR,"Fix property/atom array name already exists");
      cols[nvalue] = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      if (cols[nvalue] < 1)
	error->all(FLERR,"Invalid array columns in fix property/atom");
      index[nvalue] = atom->add_custom(&arg[iarg][3],which,cols[nvalue]);
      values_peratom += cols[nvalue];
      nvalue++;
      iarg += 2;
    } else break;
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

  if (border) comm_border = values_peratom;

  // warn if mol or charge keyword used without ghost yes

  if (border == 0) {
    int flag = 0;
    for (int i = 0; i < nvalue; i++)
      if (style[i] == MOLECULE
          || style[i] == CHARGE
          || style[i] == RMASS) flag = 1;
    if (flag && comm->me == 0)
      error->warning(FLERR,"Fix property/atom mol or charge or rmass "
                     "w/out ghost communication");
  }

  // store current atom style

  int n = strlen(atom->atom_style) + 1;
  astyle = new char[n];
  strcpy(astyle,atom->atom_style);

  // perform initial allocation of atom-based array
  // register with Atom class

  nmax_old = 0;
  if (!lmp->kokkos) grow_arrays(atom->nmax);
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

  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) {
      atom->molecule_flag = 0;
      memory->destroy(atom->molecule);
      atom->molecule = NULL;
    } else if (style[nv] == CHARGE) {
      atom->q_flag = 0;
      memory->destroy(atom->q);
      atom->q = NULL;
    } else if (style[nv] == RMASS) {
      atom->rmass_flag = 0;
      memory->destroy(atom->rmass);
      atom->rmass = NULL;
    } else if (style[nv] == IVEC) {
      atom->remove_custom(index[nv],0,cols[nv]);
    } else if (style[nv] == DVEC) {
      atom->remove_custom(index[nv],1,cols[nv]);
    } else if (style[nv] == IARRAY) {
      atom->remove_custom(index[nv],0,cols[nv]);
    } else if (style[nv] == DARRAY) {
      atom->remove_custom(index[nv],1,cols[nv]);
    }
  }

  delete [] style;
  delete [] cols;
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
  // don't allow this because user could change to style that defines molecule,q

  if (strcmp(astyle,atom->atom_style) != 0)
    error->all(FLERR,"Atom style was redefined after using fix property/atom");
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to first atomID field if multiple data files are read
------------------------------------------------------------------------- */

void FixPropertyAtom::read_data_section(char *keyword, int n, char *buf,
                                        tagint id_offset)
{
  int j,k,m,iword,ncol,nv;
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
  int nwords = utils::trim_and_count_words(buf);
  *next = '\n';

  if (nwords != values_peratom+1)
    error->all(FLERR,fmt::format("Incorrect {} format in data file",keyword));

  char **values = new char*[nwords];

  // loop over lines of atom info
  // tokenize the line into values
  // if I own atom tag, unpack its values

  tagint map_tag_max = atom->map_tag_max;

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');

    values[0] = strtok(buf," \t\n\r\f");
    if (values[0] == NULL)
      error->all(FLERR,fmt::format("Too few lines in {} section of data file",keyword));

    int format_ok = 1;
    for (j = 1; j < nwords; j++) {
      values[j] = strtok(NULL," \t\n\r\f");
      if (values[j] == NULL) format_ok = 0;
    }
    if (!format_ok)
      error->all(FLERR,fmt::format("Incorrect {} format in data file",keyword));

    itag = ATOTAGINT(values[0]) + id_offset;
    if (itag <= 0 || itag > map_tag_max)
      error->all(FLERR,fmt::format("Invalid atom ID {} in {} section of "
                                   "data file",itag, keyword));

    // assign words in line to per-atom vectors
    // iword = position in vector of words
    
    if ((m = atom->map(itag)) >= 0) {
      iword = 1;
      for (nv = 0; nv < nvalue; nv++) {
        if (style[nv] == MOLECULE) atom->molecule[m] = ATOTAGINT(values[iword++]);
        else if (style[nv] == CHARGE) atom->q[m] = atof(values[iword++]);
        else if (style[nv] == RMASS) atom->rmass[m] = atof(values[iword++]);
        else if (style[nv] == IVEC)
          atom->ivector[index[nv]][m] = atoi(values[iword++]);
        else if (style[nv] == DVEC)
          atom->dvector[index[nv]][m] = atof(values[iword++]);
        else if (style[nv] == IARRAY) {
	  ncol = cols[nv];
	  for (k = 0; k < ncol; k++)
	    atom->iarray[index[nv]][m][k] = atoi(values[iword++]);
        } else if (style[nv] == DARRAY) {
	  ncol = cols[nv];
	  for (k = 0; k < ncol; k++)
	    atom->darray[index[nv]][m][k] = atof(values[iword++]);
	}
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

bigint FixPropertyAtom::read_data_skip_lines(char * /*keyword*/)
{
  return atom->natoms;
}

/* ----------------------------------------------------------------------
   return size I own for Mth data section
   # of data sections = 1 for this fix
   nx = # of local atoms
   ny = columns = tag + values_peratom
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_size(int /*mth*/, int &nx, int &ny)
{
  nx = atom->nlocal;
  ny = values_peratom + 1;
}

/* ----------------------------------------------------------------------
   pack values for Mth data section into 2d buf
   buf allocated by caller as Nlocal by Nvalues+1
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_pack(int /*mth*/, double **buf)
{
  int i,k,ncol;

  // 1st column = atom tag
  // rest of columns = per-atom values

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) buf[i][0] = ubuf(tag[i]).d;

  int icol = 1;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      for (i = 0; i < nlocal; i++) buf[i][icol] = ubuf(molecule[i]).d;
      icol++;
    } else if (style[nv] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < nlocal; i++) buf[i][icol] = q[i];
      icol++;
    } else if (style[nv] == RMASS) {
      double *rmass = atom->rmass;
      for (i = 0; i < nlocal; i++) buf[i][icol] = rmass[i];
      icol++;
    } else if (style[nv] == IVEC) {
      int *ivec = atom->ivector[index[nv]];
      for (i = 0; i < nlocal; i++) buf[i][icol] = ubuf(ivec[i]).d;
      icol++;
    } else if (style[nv] == DVEC) {
      double *dvec = atom->dvector[index[nv]];
      for (i = 0; i < nlocal; i++) buf[i][icol] = dvec[i];
      icol++;
    } else if (style[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < nlocal; i++)
	for (k = 0; k < ncol; k++)
	  buf[i][icol+k] = ubuf(iarray[i][k]).d;
      icol += ncol;
    } else if (style[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < nlocal; i++)
	for (k = 0; k < ncol; k++)
	  buf[i][icol+k] = ubuf(darray[i][k]).d;
      icol += ncol;
    }
  }
}

/* ----------------------------------------------------------------------
   write section keyword for Mth data section to file
   use Molecules or Charges if that is only field, else use fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_keyword(int /*mth*/, FILE *fp)
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

void FixPropertyAtom::write_data_section(int /*mth*/, FILE *fp,
                                         int n, double **buf, int /*index*/)
{
  int k,icol,ncol,nv;

  for (int i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT,(tagint) ubuf(buf[i][0]).i);
    icol = 1;
    for (nv = 0; nv < nvalue; nv++) {
      if (style[nv] == MOLECULE)
        fprintf(fp," " TAGINT_FORMAT,(tagint) ubuf(buf[i][icol++]).i);
      else if (style[nv] == CHARGE)
	fprintf(fp," %g",buf[i][icol++]);
      else if (style[nv] == RMASS)
	fprintf(fp," %g",buf[i][icol++]);
      else if (style[nv] == IVEC)
        fprintf(fp," %d",(int) ubuf(buf[i][icol++]).i);
      else if (style[nv] == DVEC)
	fprintf(fp," %g",buf[i][icol++]);
      else if (style[nv] == IARRAY) {
	ncol = cols[nv];
	for (k = 0; k < ncol; k++)
	  fprintf(fp," %d",(int) ubuf(buf[i][icol+k]).i);
	icol += ncol;
      } else if (style[nv] == DARRAY) {
	ncol = cols[nv];
	for (k = 0; k < ncol; k++)
	  fprintf(fp," %g",buf[i][icol+k]);
	icol += ncol;
      }
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
    else if (style[m] == RMASS) bytes = atom->nmax * sizeof(double);
    else if (style[m] == IVEC) bytes = atom->nmax * sizeof(int);
    else if (style[m] == DVEC) bytes = atom->nmax * sizeof(double);
    else if (style[m] == IARRAY) bytes = atom->nmax * cols[m] * sizeof(int);
    else if (style[m] == DARRAY) bytes = atom->nmax * cols[m] * sizeof(double);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
   also initialize new values to 0
     since AtomVec class won't do it as atoms are added,
     e.g. in create_atom() or data_atom()
------------------------------------------------------------------------- */

void FixPropertyAtom::grow_arrays(int nmax)
{
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) {
      memory->grow(atom->molecule,nmax,"atom:molecule");
      size_t nbytes = (nmax-nmax_old) * sizeof(tagint);
      memset(&atom->molecule[nmax_old],0,nbytes);
    } else if (style[nv] == CHARGE) {
      memory->grow(atom->q,nmax,"atom:q");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->q[nmax_old],0,nbytes);
    } else if (style[nv] == RMASS) {
      memory->grow(atom->rmass,nmax,"atom:rmass");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->rmass[nmax_old],0,nbytes);
    } else if (style[nv] == IVEC) {
      memory->grow(atom->ivector[index[nv]],nmax,"atom:ivector");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[nv]][nmax_old],0,nbytes);
    } else if (style[nv] == DVEC) {
      memory->grow(atom->dvector[index[nv]],nmax,"atom:dvector");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->dvector[index[nv]][nmax_old],0,nbytes);
    } else if (style[nv] == IARRAY) {
      memory->grow(atom->iarray[index[nv]],nmax,cols[nv],"atom:iarray");
      size_t nbytes = (nmax-nmax_old) * cols[nv] * sizeof(int);
      if (nbytes) memset(&atom->iarray[index[nv]][nmax_old][0],0,nbytes);
    } else if (style[nv] == DARRAY) {
      memory->grow(atom->darray[index[nv]],nmax,cols[nv],"atom:darray");
      size_t nbytes = (nmax-nmax_old) * cols[nv] * sizeof(double);
      if (nbytes) memset(&atom->darray[index[nv]][nmax_old][0],0,nbytes);
    }
  }

  nmax_old = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtom::copy_arrays(int i, int j, int /*delflag*/)
{
  int k,ncol;
  
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE)
      atom->molecule[j] = atom->molecule[i];
    else if (style[nv] == CHARGE)
      atom->q[j] = atom->q[i];
    else if (style[nv] == RMASS)
      atom->rmass[j] = atom->rmass[i];
    else if (style[nv] == IVEC)
      atom->ivector[index[nv]][j] = atom->ivector[index[nv]][i];
    else if (style[nv] == DVEC)
      atom->dvector[index[nv]][j] = atom->dvector[index[nv]][i];
    else if (style[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	atom->iarray[index[nv]][j][k] = atom->iarray[index[nv]][i][k];
    } else if (style[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	atom->darray[index[nv]][j][k] = atom->darray[index[nv]][i][k];
    }
  }
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_border(int n, int *list, double *buf)
{
  int i,j,k,ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(molecule[j]).d;
      }
    } else if (style[nv] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = q[j];
      }
    } else if (style[nv] == RMASS) {
      double *rmass = atom->rmass;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = rmass[j];
      }
    } else if (style[nv] == IVEC) {
      int *ivector = atom->ivector[index[nv]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(ivector[j]).d;
      }
    } else if (style[nv] == DVEC) {
      double *dvector = atom->dvector[index[nv]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = dvector[j];
      }
    } else if (style[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < n; i++) {
        j = list[i];
	for (k = 0; k < ncol; k++)
	  buf[m++] = ubuf(iarray[j][k]).d;
      }
    } else if (style[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < n; i++) {
        j = list[i];
	for (k = 0; k < ncol; k++)
	  buf[m++] = darray[j][k];
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
  int i,k,last,ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      last = first + n;
      for (i = first; i < last; i++)
        molecule[i] = (tagint) ubuf(buf[m++]).i;
    } else if (style[nv] == CHARGE) {
      double *q = atom->q;
      last = first + n;
      for (i = first; i < last; i++)
        q[i] = buf[m++];
    } else if (style[nv] == RMASS) {
      double *rmass = atom->rmass;
      last = first + n;
      for (i = first; i < last; i++)
        rmass[i] = buf[m++];
    } else if (style[nv] == IVEC) {
      int *ivector = atom->ivector[index[nv]];
      last = first + n;
      for (i = first; i < last; i++)
        ivector[i] = (int) ubuf(buf[m++]).i;
    } else if (style[nv] == DVEC) {
      double *dvector = atom->dvector[index[nv]];
      last = first + n;
      for (i = first; i < last; i++)
        dvector[i] = buf[m++];
    } else if (style[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      last = first + n;
      for (i = first; i < last; i++)
	for (k = 0; k < ncol; k++)
	  iarray[i][k] = (int) ubuf(buf[m++]).i;
    } else if (style[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      last = first + n;
      for (i = first; i < last; i++)
	for (k = 0; k < ncol; k++)
	  darray[i][k] = buf[m++];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_exchange(int i, double *buf)
{
  int k,ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) buf[m++] = ubuf(atom->molecule[i]).d;
    else if (style[nv] == CHARGE) buf[m++] = atom->q[i];
    else if (style[nv] == RMASS) buf[m++] = atom->rmass[i];
    else if (style[nv] == IVEC) buf[m++] = ubuf(atom->ivector[index[nv]][i]).d;
    else if (style[nv] == DVEC) buf[m++] = atom->dvector[index[nv]][i];
    else if (style[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) 
	buf[m++] = ubuf(atom->iarray[index[nv]][i][k]).d;
    } else if (style[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) 
	buf[m++] = atom->darray[index[nv]][i][k];
    }
  }
  
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::unpack_exchange(int nlocal, double *buf)
{
  int k,ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE)
      atom->molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
    else if (style[nv] == CHARGE)
      atom->q[nlocal] = buf[m++];
    else if (style[nv] == RMASS)
      atom->rmass[nlocal] = buf[m++];
    else if (style[nv] == IVEC)
      atom->ivector[index[nv]][nlocal] = (int) ubuf(buf[m++]).i;
    else if (style[nv] == DVEC)
      atom->dvector[index[nv]][nlocal] = buf[m++];
    else if (style[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) 
	atom->iarray[index[nv]][nlocal][k] = (int) ubuf(buf[m++]).i;
    } else if (style[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) 
	atom->darray[index[nv]][nlocal][k] = buf[m++];
    }
  }
  
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_restart(int i, double *buf)
{
  int k,ncol;
  
  // pack buf[0] this way because other fixes unpack it
  
  buf[0] = values_peratom+1;

  int m = 1;
  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE) buf[m++] = ubuf(atom->molecule[i]).d;
    else if (style[nv] == CHARGE) buf[m++] = atom->q[i];
    else if (style[nv] == RMASS) buf[m++] = atom->rmass[i];
    else if (style[nv] == IVEC) buf[m++] = ubuf(atom->ivector[index[nv]][i]).d;
    else if (style[nv] == DVEC) buf[m++] = atom->dvector[index[nv]][i];
    else if (style[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	buf[m++] = ubuf(atom->iarray[index[nv]][i][k]).d;
    } else if (style[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	buf[m++] = atom->darray[index[nv]][i][k];
    }
  }

  return values_peratom+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyAtom::unpack_restart(int nlocal, int nth)
{
  int k,ncol;
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  for (int nv = 0; nv < nvalue; nv++) {
    if (style[nv] == MOLECULE)
      atom->molecule[nlocal] = (tagint) ubuf(extra[nlocal][m++]).i;
    else if (style[nv] == CHARGE)
      atom->q[nlocal] = extra[nlocal][m++];
    else if (style[nv] == RMASS)
      atom->rmass[nlocal] = extra[nlocal][m++];
    else if (style[nv] == IVEC)
      atom->ivector[index[nv]][nlocal] = (int) ubuf(extra[nlocal][m++]).i;
    else if (style[nv] == DVEC)
      atom->dvector[index[nv]][nlocal] = extra[nlocal][m++];
    else if (style[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	atom->iarray[index[nv]][nlocal][k] = (int) ubuf(extra[nlocal][m++]).i;
    } else if (style[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
	atom->darray[index[nv]][nlocal][k] = extra[nlocal][m++];
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::maxsize_restart()
{
  return values_peratom+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::size_restart(int /*nlocal*/)
{
  return values_peratom+1;
}
