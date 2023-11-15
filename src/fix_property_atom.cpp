/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_property_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "read_data.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtom::FixPropertyAtom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nvalue(0), styles(nullptr), index(nullptr), astyle(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix property/atom command");

  restart_peratom = 1;
  wd_section = 1;    // can be overwitten using optional arguments

  int iarg = 3;
  nvalue = narg - iarg;
  styles = new int[nvalue];
  cols = new int[nvalue];
  index = new int[nvalue];

  molecule_flag = 0;
  q_flag = 0;
  rmass_flag = 0;
  temperature_flag = 0;
  heatflow_flag = 0;
  nmax_old = 0;

  nvalue = 0;
  values_peratom = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "mol") == 0) {
      if (atom->molecule_flag)
        error->all(FLERR, "Fix property/atom mol when atom_style already has molecule attribute");
      if (molecule_flag) error->all(FLERR, "Fix property/atom cannot specify mol twice");
      styles[nvalue] = MOLECULE;
      cols[nvalue] = 0;
      atom->molecule_flag = molecule_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg], "q") == 0) {
      if (atom->q_flag)
        error->all(FLERR, "Fix property/atom q when atom_style already has charge attribute");
      if (q_flag) error->all(FLERR, "Fix property/atom cannot specify q twice");
      styles[nvalue] = CHARGE;
      cols[nvalue] = 0;
      atom->q_flag = q_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg], "rmass") == 0) {
      if (atom->rmass_flag)
        error->all(FLERR, "Fix property/atom rmass when atom_style already has rmass attribute");
      if (rmass_flag) error->all(FLERR, "Fix property/atom cannot specify rmass twice");
      styles[nvalue] = RMASS;
      cols[nvalue] = 0;
      atom->rmass_flag = rmass_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg], "temperature") == 0) {
      if (atom->temperature_flag)
        error->all(FLERR, "Fix property/atom temperature when atom_style already has temperature attribute");
      if (temperature_flag) error->all(FLERR, "Fix property/atom cannot specify temperature twice");
      styles[nvalue] = TEMPERATURE;
      cols[nvalue] = 0;
      atom->temperature_flag = temperature_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;
    } else if (strcmp(arg[iarg], "heatflow") == 0) {
      if (atom->heatflow_flag)
        error->all(FLERR, "Fix property/atom heatflow when atom_style already has heatflow attribute");
      if (heatflow_flag) error->all(FLERR, "Fix property/atom cannot specify heatflow twice");
      styles[nvalue] = HEATFLOW;
      cols[nvalue] = 0;
      atom->heatflow_flag = heatflow_flag = 1;
      values_peratom++;
      nvalue++;
      iarg++;

      // custom atom vector

    } else if (utils::strmatch(arg[iarg], "^i_")) {
      styles[nvalue] = IVEC;
      int flag, ncols;
      index[nvalue] = atom->find_custom(&arg[iarg][2], flag, ncols);
      if (index[nvalue] >= 0) error->all(FLERR, "Fix property/atom vector name already exists");
      if (ReadData::is_data_section(id))
        error->all(FLERR, "Fix property/atom fix ID must not be a data file section name");
      index[nvalue] = atom->add_custom(&arg[iarg][2], 0, 0);
      cols[nvalue] = 0;
      values_peratom++;
      nvalue++;
      iarg++;

    } else if (utils::strmatch(arg[iarg], "^d_")) {
      styles[nvalue] = DVEC;
      int flag, ncols;
      index[nvalue] = atom->find_custom(&arg[iarg][2], flag, ncols);
      if (index[nvalue] >= 0) error->all(FLERR, "Fix property/atom vector name already exists");
      if (ReadData::is_data_section(id))
        error->all(FLERR, "Fix property/atom fix ID must not be a data file section name");
      index[nvalue] = atom->add_custom(&arg[iarg][2], 1, 0);
      cols[nvalue] = 0;
      values_peratom++;
      nvalue++;
      iarg++;

      // custom atom array

    } else if (utils::strmatch(arg[iarg], "^[id]2_")) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix property/atom command");

      int which, flag, ncols;
      which = atom->find_custom(&arg[iarg][3], flag, ncols);
      if (which >= 0)
        error->all(FLERR, "Fix property/atom array name {} already exists", &arg[iarg][3]);
      if (ReadData::is_data_section(id))
        error->all(FLERR, "Fix property/atom fix ID must not be a data file section name");

      ncols = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      if (ncols < 1)
        error->all(FLERR, "Invalid array columns number {} in fix property/atom", ncols);

      if (arg[iarg][0] == 'i') {
        which = 0;
        styles[nvalue] = IARRAY;
      } else {
        which = 1;
        styles[nvalue] = DARRAY;
      }
      index[nvalue] = atom->add_custom(&arg[iarg][3], which, ncols);
      cols[nvalue] = ncols;
      values_peratom += ncols;
      nvalue++;
      iarg += 2;

      // no match

    } else
      break;
  }

  // optional args

  border = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ghost") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix property/atom command");
      border = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "writedata") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix property/atom command");
      wd_section = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix property/atom command");
  }

  if (border) comm_border = values_peratom;

  // warn if mol, charge, rmass, temperature, or heatflow keyword used without ghost yes

  if (border == 0) {
    int flag = 0;
    for (int i = 0; i < nvalue; i++)
      if (styles[i] == MOLECULE || styles[i] == CHARGE || styles[i] == RMASS ||
      styles[i] == TEMPERATURE || styles[i] == HEATFLOW) flag = 1;
    if (flag && comm->me == 0)
      error->warning(FLERR, "Fix property/atom mol, charge, rmass, temperature, or heatflow w/out ghost communication");
  }

  // store current atom style

  astyle = utils::strdup(atom->atom_style);

  // register with Atom class

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);
  if (border) atom->add_callback(Atom::BORDER);
}


/* ---------------------------------------------------------------------- */

void FixPropertyAtom::post_constructor()
{
  // perform initial allocation of atom-based array

  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

FixPropertyAtom::~FixPropertyAtom()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, Atom::GROW);
  atom->delete_callback(id, Atom::RESTART);
  if (border) atom->delete_callback(id, Atom::BORDER);

  // deallocate per-atom vectors in Atom class
  // set ptrs to a null pointer, so they no longer exist for Atom class

  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      atom->molecule_flag = 0;
      memory->destroy(atom->molecule);
      atom->molecule = nullptr;
    } else if (styles[nv] == CHARGE) {
      atom->q_flag = 0;
      memory->destroy(atom->q);
      atom->q = nullptr;
    } else if (styles[nv] == RMASS) {
      atom->rmass_flag = 0;
      memory->destroy(atom->rmass);
      atom->rmass = nullptr;
    } else if (styles[nv] == TEMPERATURE) {
      atom->temperature_flag = 0;
      memory->destroy(atom->temperature);
      atom->temperature = nullptr;
    } else if (styles[nv] == HEATFLOW) {
      atom->heatflow_flag = 0;
      memory->destroy(atom->heatflow);
      atom->heatflow = nullptr;
    } else if (styles[nv] == IVEC) {
      atom->remove_custom(index[nv], 0, cols[nv]);
    } else if (styles[nv] == DVEC) {
      atom->remove_custom(index[nv], 1, cols[nv]);
    } else if (styles[nv] == IARRAY) {
      atom->remove_custom(index[nv], 0, cols[nv]);
    } else if (styles[nv] == DARRAY) {
      atom->remove_custom(index[nv], 1, cols[nv]);
    }
  }

  delete[] styles;
  delete[] cols;
  delete[] index;
  delete[] astyle;
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

  if (strcmp(astyle, atom->atom_style) != 0)
    error->all(FLERR, "Atom style was redefined after using fix property/atom");
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to first atomID field if multiple data files are read
------------------------------------------------------------------------- */

void FixPropertyAtom::read_data_section(char *keyword, int n, char *buf, tagint id_offset)
{
  int j, k, m, ncol;
  tagint itag;
  char *next;

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  // loop over lines of atom info
  // tokenize the line into values
  // if I own atom tag, assign its values

  tagint map_tag_max = atom->map_tag_max;

  for (int i = 0; i < n; i++) {
    next = strchr(buf, '\n');
    if (!next) error->all(FLERR, "Unexpected end of file while reading data section");
    *next = '\0';

    try {
      ValueTokenizer values(buf);
      if ((int) values.count() != values_peratom + 1)
        error->all(FLERR, "Incorrect format in {} section of data file: {} expected {} and got {}",
                   keyword, buf, values_peratom + 1, values.count());

      itag = values.next_tagint() + id_offset;
      if (itag <= 0 || itag > map_tag_max)
        error->all(FLERR, "Invalid atom ID {} in {} section of data file", itag, keyword);

      // assign words in line to per-atom vectors

      if ((m = atom->map(itag)) >= 0) {
        for (j = 0; j < nvalue; j++) {
          if (styles[j] == MOLECULE) {
            atom->molecule[m] = values.next_tagint();
          } else if (styles[j] == CHARGE) {
            atom->q[m] = values.next_double();
          } else if (styles[j] == RMASS) {
            atom->rmass[m] = values.next_double();
          } else if (styles[j] == TEMPERATURE) {
            atom->temperature[m] = values.next_double();
          } else if (styles[j] == HEATFLOW) {
            atom->heatflow[m] = values.next_double();
          } else if (styles[j] == IVEC) {
            atom->ivector[index[j]][m] = values.next_int();
          } else if (styles[j] == DVEC) {
            atom->dvector[index[j]][m] = values.next_double();
          } else if (styles[j] == IARRAY) {
            ncol = cols[j];
            for (k = 0; k < ncol; k++) atom->iarray[index[j]][m][k] = values.next_int();
          } else if (styles[j] == DARRAY) {
            ncol = cols[j];
            for (k = 0; k < ncol; k++) atom->darray[index[j]][m][k] = values.next_double();
          }
        }
      }
    } catch (TokenizerException &e) {
      error->all(FLERR, "Invalid format in {} section of data file '{}': {}", keyword, buf,
                 e.what());
    }
    buf = next + 1;
  }

  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }
}

/* ----------------------------------------------------------------------
   return # of lines in section of data file labeled by keyword. -1 signals use # of added atoms
------------------------------------------------------------------------- */

bigint FixPropertyAtom::read_data_skip_lines(char * /*keyword*/)
{
  return -1;
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
  int i, k, ncol;

  // 1st column = atom tag
  // rest of columns = per-atom values

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) buf[i][0] = ubuf(tag[i]).d;

  int icol = 1;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      for (i = 0; i < nlocal; i++) buf[i][icol] = ubuf(molecule[i]).d;
      icol++;
    } else if (styles[nv] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < nlocal; i++) buf[i][icol] = q[i];
      icol++;
    } else if (styles[nv] == RMASS) {
      double *rmass = atom->rmass;
      for (i = 0; i < nlocal; i++) buf[i][icol] = rmass[i];
      icol++;
    } else if (styles[nv] == TEMPERATURE) {
      double *temperature = atom->temperature;
      for (i = 0; i < nlocal; i++) buf[i][icol] = temperature[i];
      icol++;
    } else if (styles[nv] == HEATFLOW) {
      double *heatflow = atom->heatflow;
      for (i = 0; i < nlocal; i++) buf[i][icol] = heatflow[i];
      icol++;
    } else if (styles[nv] == IVEC) {
      int *ivec = atom->ivector[index[nv]];
      for (i = 0; i < nlocal; i++) buf[i][icol] = ubuf(ivec[i]).d;
      icol++;
    } else if (styles[nv] == DVEC) {
      double *dvec = atom->dvector[index[nv]];
      for (i = 0; i < nlocal; i++) buf[i][icol] = dvec[i];
      icol++;
    } else if (styles[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < nlocal; i++)
        for (k = 0; k < ncol; k++) buf[i][icol + k] = ubuf(iarray[i][k]).d;
      icol += ncol;
    } else if (styles[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < nlocal; i++)
        for (k = 0; k < ncol; k++) buf[i][icol + k] = darray[i][k];
      icol += ncol;
    }
  }
}

/* ----------------------------------------------------------------------
   write a Molecules or Charges section if that is only field
   manages by this fix instance. Otherwise write everything
   to the section labeled by the fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section_keyword(int /*mth*/, FILE *fp)
{
  if (nvalue == 1 && styles[0] == MOLECULE)
    fprintf(fp, "\nMolecules\n\n");
  else if (nvalue == 1 && styles[0] == CHARGE)
    fprintf(fp, "\nCharges\n\n");
  else {
    fprintf(fp, "\n%s #", id);
    // write column hint as comment
    for (int i = 0; i < nvalue; ++i) {
      if (styles[i] == MOLECULE)
        fputs(" mol", fp);
      else if (styles[i] == CHARGE)
        fputs(" q", fp);
      else if (styles[i] == RMASS)
        fputs(" rmass", fp);
      else if (styles[i] == TEMPERATURE)
        fputs(" temperature", fp);
      else if (styles[i] == HEATFLOW)
        fputs(" heatflow", fp);
      else if (styles[i] == IVEC)
        fprintf(fp, " i_%s", atom->ivname[index[i]]);
      else if (styles[i] == DVEC)
        fprintf(fp, " d_%s", atom->dvname[index[i]]);
      else if (styles[i] == IARRAY)
        fprintf(fp, " i_%s", atom->ianame[index[i]]);
      else if (styles[i] == DARRAY)
        fprintf(fp, " d_%s", atom->daname[index[i]]);
    }
    fputs("\n\n", fp);
  }
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixPropertyAtom::write_data_section(int /*mth*/, FILE *fp, int n, double **buf, int /*index*/)
{
  int k, icol, ncol, nv;
  std::string line;

  for (int i = 0; i < n; i++) {
    line = fmt::format("{}", (tagint) ubuf(buf[i][0]).i);
    icol = 1;
    for (nv = 0; nv < nvalue; nv++) {
      if (styles[nv] == MOLECULE)
        line += fmt::format(" {}", (tagint) ubuf(buf[i][icol++]).i);
      else if (styles[nv] == CHARGE)
        line += fmt::format(" {}", buf[i][icol++]);
      else if (styles[nv] == RMASS)
        line += fmt::format(" {}", buf[i][icol++]);
      else if (styles[nv] == TEMPERATURE)
        line += fmt::format(" {}", buf[i][icol++]);
      else if (styles[nv] == HEATFLOW)
        line += fmt::format(" {}", buf[i][icol++]);
      else if (styles[nv] == IVEC)
        line += fmt::format(" {}", (int) ubuf(buf[i][icol++]).i);
      else if (styles[nv] == DVEC)
        line += fmt::format(" {}", buf[i][icol++]);
      else if (styles[nv] == IARRAY) {
        ncol = cols[nv];
        for (k = 0; k < ncol; k++) line += fmt::format(" {}", (int) ubuf(buf[i][icol + k]).i);
        icol += ncol;
      } else if (styles[nv] == DARRAY) {
        ncol = cols[nv];
        for (k = 0; k < ncol; k++) line += fmt::format(" {}", buf[i][icol + k]);
        icol += ncol;
      }
    }
    fmt::print(fp, line + "\n");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyAtom::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvalue; m++) {
    if (styles[m] == MOLECULE)
      bytes = atom->nmax * sizeof(tagint);
    else if (styles[m] == CHARGE)
      bytes = atom->nmax * sizeof(double);
    else if (styles[m] == RMASS)
      bytes = atom->nmax * sizeof(double);
    else if (styles[m] == TEMPERATURE)
      bytes = atom->nmax * sizeof(double);
    else if (styles[m] == HEATFLOW)
      bytes = atom->nmax * sizeof(double);
    else if (styles[m] == IVEC)
      bytes = atom->nmax * sizeof(int);
    else if (styles[m] == DVEC)
      bytes = atom->nmax * sizeof(double);
    else if (styles[m] == IARRAY)
      bytes = (size_t) atom->nmax * cols[m] * sizeof(int);
    else if (styles[m] == DARRAY)
      bytes = (size_t) atom->nmax * cols[m] * sizeof(double);
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
    if (styles[nv] == MOLECULE) {
      memory->grow(atom->molecule, nmax, "atom:molecule");
      size_t nbytes = (nmax - nmax_old) * sizeof(tagint);
      memset(&atom->molecule[nmax_old], 0, nbytes);
    } else if (styles[nv] == CHARGE) {
      memory->grow(atom->q, nmax, "atom:q");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->q[nmax_old], 0, nbytes);
    } else if (styles[nv] == RMASS) {
      memory->grow(atom->rmass, nmax, "atom:rmass");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->rmass[nmax_old], 0, nbytes);
    } else if (styles[nv] == TEMPERATURE) {
      memory->grow(atom->temperature, nmax, "atom:temperature");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->temperature[nmax_old], 0, nbytes);
    } else if (styles[nv] == HEATFLOW) {
      memory->grow(atom->heatflow, nmax, "atom:heatflow");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->heatflow[nmax_old], 0, nbytes);
    } else if (styles[nv] == IVEC) {
      memory->grow(atom->ivector[index[nv]], nmax, "atom:ivector");
      size_t nbytes = (nmax - nmax_old) * sizeof(int);
      memset(&atom->ivector[index[nv]][nmax_old], 0, nbytes);
    } else if (styles[nv] == DVEC) {
      memory->grow(atom->dvector[index[nv]], nmax, "atom:dvector");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->dvector[index[nv]][nmax_old], 0, nbytes);
    } else if (styles[nv] == IARRAY) {
      memory->grow(atom->iarray[index[nv]], nmax, cols[nv], "atom:iarray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(int);
      if (nbytes) memset(&atom->iarray[index[nv]][nmax_old][0], 0, nbytes);
    } else if (styles[nv] == DARRAY) {
      memory->grow(atom->darray[index[nv]], nmax, cols[nv], "atom:darray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(double);
      if (nbytes) memset(&atom->darray[index[nv]][nmax_old][0], 0, nbytes);
    }
  }

  nmax_old = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyAtom::copy_arrays(int i, int j, int /*delflag*/)
{
  int k, ncol;

  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE)
      atom->molecule[j] = atom->molecule[i];
    else if (styles[nv] == CHARGE)
      atom->q[j] = atom->q[i];
    else if (styles[nv] == RMASS)
      atom->rmass[j] = atom->rmass[i];
    else if (styles[nv] == TEMPERATURE)
      atom->temperature[j] = atom->temperature[i];
    else if (styles[nv] == HEATFLOW)
      atom->heatflow[j] = atom->heatflow[i];
    else if (styles[nv] == IVEC)
      atom->ivector[index[nv]][j] = atom->ivector[index[nv]][i];
    else if (styles[nv] == DVEC)
      atom->dvector[index[nv]][j] = atom->dvector[index[nv]][i];
    else if (styles[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) atom->iarray[index[nv]][j][k] = atom->iarray[index[nv]][i][k];
    } else if (styles[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) atom->darray[index[nv]][j][k] = atom->darray[index[nv]][i][k];
    }
  }
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_border(int n, int *list, double *buf)
{
  int i, j, k, ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(molecule[j]).d;
      }
    } else if (styles[nv] == CHARGE) {
      double *q = atom->q;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = q[j];
      }
    } else if (styles[nv] == RMASS) {
      double *rmass = atom->rmass;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = rmass[j];
      }
    } else if (styles[nv] == TEMPERATURE) {
      double *temperature = atom->temperature;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = temperature[j];
      }
    } else if (styles[nv] == HEATFLOW) {
      double *heatflow = atom->heatflow;
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = heatflow[j];
      }
    } else if (styles[nv] == IVEC) {
      int *ivector = atom->ivector[index[nv]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = ubuf(ivector[j]).d;
      }
    } else if (styles[nv] == DVEC) {
      double *dvector = atom->dvector[index[nv]];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = dvector[j];
      }
    } else if (styles[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < n; i++) {
        j = list[i];
        for (k = 0; k < ncol; k++) buf[m++] = ubuf(iarray[j][k]).d;
      }
    } else if (styles[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      for (i = 0; i < n; i++) {
        j = list[i];
        for (k = 0; k < ncol; k++) buf[m++] = darray[j][k];
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
  int i, k, last, ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      tagint *molecule = atom->molecule;
      last = first + n;
      for (i = first; i < last; i++) molecule[i] = (tagint) ubuf(buf[m++]).i;
    } else if (styles[nv] == CHARGE) {
      double *q = atom->q;
      last = first + n;
      for (i = first; i < last; i++) q[i] = buf[m++];
    } else if (styles[nv] == RMASS) {
      double *rmass = atom->rmass;
      last = first + n;
      for (i = first; i < last; i++) rmass[i] = buf[m++];
    } else if (styles[nv] == TEMPERATURE) {
      double *temperature = atom->temperature;
      last = first + n;
      for (i = first; i < last; i++) temperature[i] = buf[m++];
    } else if (styles[nv] == HEATFLOW) {
      double *heatflow = atom->heatflow;
      last = first + n;
      for (i = first; i < last; i++) heatflow[i] = buf[m++];
    } else if (styles[nv] == IVEC) {
      int *ivector = atom->ivector[index[nv]];
      last = first + n;
      for (i = first; i < last; i++) ivector[i] = (int) ubuf(buf[m++]).i;
    } else if (styles[nv] == DVEC) {
      double *dvector = atom->dvector[index[nv]];
      last = first + n;
      for (i = first; i < last; i++) dvector[i] = buf[m++];
    } else if (styles[nv] == IARRAY) {
      int **iarray = atom->iarray[index[nv]];
      ncol = cols[nv];
      last = first + n;
      for (i = first; i < last; i++)
        for (k = 0; k < ncol; k++) iarray[i][k] = (int) ubuf(buf[m++]).i;
    } else if (styles[nv] == DARRAY) {
      double **darray = atom->darray[index[nv]];
      ncol = cols[nv];
      last = first + n;
      for (i = first; i < last; i++)
        for (k = 0; k < ncol; k++) darray[i][k] = buf[m++];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_exchange(int i, double *buf)
{
  int k, ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE)
      buf[m++] = ubuf(atom->molecule[i]).d;
    else if (styles[nv] == CHARGE)
      buf[m++] = atom->q[i];
    else if (styles[nv] == RMASS)
      buf[m++] = atom->rmass[i];
    else if (styles[nv] == TEMPERATURE)
      buf[m++] = atom->temperature[i];
    else if (styles[nv] == HEATFLOW)
      buf[m++] = atom->heatflow[i];
    else if (styles[nv] == IVEC)
      buf[m++] = ubuf(atom->ivector[index[nv]][i]).d;
    else if (styles[nv] == DVEC)
      buf[m++] = atom->dvector[index[nv]][i];
    else if (styles[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) buf[m++] = ubuf(atom->iarray[index[nv]][i][k]).d;
    } else if (styles[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) buf[m++] = atom->darray[index[nv]][i][k];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyAtom::unpack_exchange(int nlocal, double *buf)
{
  int k, ncol;

  int m = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE)
      atom->molecule[nlocal] = (tagint) ubuf(buf[m++]).i;
    else if (styles[nv] == CHARGE)
      atom->q[nlocal] = buf[m++];
    else if (styles[nv] == RMASS)
      atom->rmass[nlocal] = buf[m++];
    else if (styles[nv] == TEMPERATURE)
      atom->temperature[nlocal] = buf[m++];
    else if (styles[nv] == HEATFLOW)
      atom->heatflow[nlocal] = buf[m++];
    else if (styles[nv] == IVEC)
      atom->ivector[index[nv]][nlocal] = (int) ubuf(buf[m++]).i;
    else if (styles[nv] == DVEC)
      atom->dvector[index[nv]][nlocal] = buf[m++];
    else if (styles[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) atom->iarray[index[nv]][nlocal][k] = (int) ubuf(buf[m++]).i;
    } else if (styles[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) atom->darray[index[nv]][nlocal][k] = buf[m++];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyAtom::pack_restart(int i, double *buf)
{
  int k, ncol;

  // pack buf[0] this way because other fixes unpack it

  buf[0] = values_peratom + 1;

  int m = 1;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE)
      buf[m++] = ubuf(atom->molecule[i]).d;
    else if (styles[nv] == CHARGE)
      buf[m++] = atom->q[i];
    else if (styles[nv] == RMASS)
      buf[m++] = atom->rmass[i];
    else if (styles[nv] == TEMPERATURE)
      buf[m++] = atom->temperature[i];
    else if (styles[nv] == HEATFLOW)
      buf[m++] = atom->heatflow[i];
    else if (styles[nv] == IVEC)
      buf[m++] = ubuf(atom->ivector[index[nv]][i]).d;
    else if (styles[nv] == DVEC)
      buf[m++] = atom->dvector[index[nv]][i];
    else if (styles[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) buf[m++] = ubuf(atom->iarray[index[nv]][i][k]).d;
    } else if (styles[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) buf[m++] = atom->darray[index[nv]][i][k];
    }
  }

  return values_peratom + 1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyAtom::unpack_restart(int nlocal, int nth)
{
  int k, ncol;
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE)
      atom->molecule[nlocal] = (tagint) ubuf(extra[nlocal][m++]).i;
    else if (styles[nv] == CHARGE)
      atom->q[nlocal] = extra[nlocal][m++];
    else if (styles[nv] == RMASS)
      atom->rmass[nlocal] = extra[nlocal][m++];
    else if (styles[nv] == TEMPERATURE)
      atom->temperature[nlocal] = extra[nlocal][m++];
    else if (styles[nv] == HEATFLOW)
      atom->heatflow[nlocal] = extra[nlocal][m++];
    else if (styles[nv] == IVEC)
      atom->ivector[index[nv]][nlocal] = (int) ubuf(extra[nlocal][m++]).i;
    else if (styles[nv] == DVEC)
      atom->dvector[index[nv]][nlocal] = extra[nlocal][m++];
    else if (styles[nv] == IARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++)
        atom->iarray[index[nv]][nlocal][k] = (int) ubuf(extra[nlocal][m++]).i;
    } else if (styles[nv] == DARRAY) {
      ncol = cols[nv];
      for (k = 0; k < ncol; k++) atom->darray[index[nv]][nlocal][k] = extra[nlocal][m++];
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::maxsize_restart()
{
  return values_peratom + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyAtom::size_restart(int /*nlocal*/)
{
  return values_peratom + 1;
}
