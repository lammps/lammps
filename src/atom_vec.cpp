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

#include "atom_vec.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

// peratom variables that are auto-included in corresponding child style field lists
// these fields cannot be specified in the fields strings

const std::vector<std::string> AtomVec::default_grow = {"id", "type", "mask", "image",
                                                        "x",  "v",    "f"};
const std::vector<std::string> AtomVec::default_copy = {"id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_comm = {"x"};
const std::vector<std::string> AtomVec::default_comm_vel = {"x", "v"};
const std::vector<std::string> AtomVec::default_reverse = {"f"};
const std::vector<std::string> AtomVec::default_border = {"id", "type", "mask", "x"};
const std::vector<std::string> AtomVec::default_border_vel = {"id", "type", "mask", "x", "v"};
const std::vector<std::string> AtomVec::default_exchange = {"id",    "type", "mask",
                                                            "image", "x",    "v"};
const std::vector<std::string> AtomVec::default_restart = {"id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_create = {"id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_data_atom = {};
const std::vector<std::string> AtomVec::default_data_vel = {};

/* ---------------------------------------------------------------------- */

AtomVec::AtomVec(LAMMPS *lmp) : Pointers(lmp)
{
  nmax = 0;
  ngrow = 0;

  molecular = Atom::ATOMIC;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 0;
  mass_type = dipole_type = PER_ATOM;
  forceclearflag = 0;
  maxexchange = 0;
  bonus_flag = 0;
  size_forward_bonus = size_border_bonus = 0;

  kokkosable = 0;

  nargcopy = 0;
  argcopy = nullptr;

  tag = nullptr;
  type = mask = nullptr;
  image = nullptr;
  x = v = f = nullptr;

  x_hold = nullptr;
  v_hold = omega_hold = angmom_hold = nullptr;

  threads = nullptr;
}

/* ---------------------------------------------------------------------- */

AtomVec::~AtomVec()
{
  int datatype, cols;
  void *pdata;

  for (int i = 0; i < nargcopy; i++) delete[] argcopy[i];
  delete[] argcopy;

  for (int i = 0; i < ngrow; i++) {
    pdata = mgrow.pdata[i];
    datatype = mgrow.datatype[i];
    cols = mgrow.cols[i];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0)
        memory->destroy(*((double **) pdata));
      else if (cols > 0)
        memory->destroy(*((double ***) pdata));
      else {
        memory->destroy(*((double ***) pdata));
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0)
        memory->destroy(*((int **) pdata));
      else if (cols > 0)
        memory->destroy(*((int ***) pdata));
      else {
        memory->destroy(*((int ***) pdata));
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0)
        memory->destroy(*((bigint **) pdata));
      else if (cols > 0)
        memory->destroy(*((bigint ***) pdata));
      else {
        memory->destroy(*((bigint ***) pdata));
      }
    }
  }

  delete[] threads;
}

/* ----------------------------------------------------------------------
   make copy of args for use by restart & replicate
------------------------------------------------------------------------- */

void AtomVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  if (nargcopy)
    argcopy = new char *[nargcopy];
  else
    argcopy = nullptr;
  for (int i = 0; i < nargcopy; i++) argcopy[i] = utils::strdup(arg[i]);
}

/* ----------------------------------------------------------------------
   no additional args by default
------------------------------------------------------------------------- */

void AtomVec::process_args(int narg, char ** /*arg*/)
{
  if (narg) error->all(FLERR, "Invalid atom_style command");
}

/* ----------------------------------------------------------------------
   pull settings from Domain needed for pack_comm_vel and pack_border_vel
   child classes may override this method, but should also invoke it
------------------------------------------------------------------------- */

void AtomVec::init()
{
  deform_vremap = domain->deform_vremap;
  deform_groupbit = domain->deform_groupbit;
  h_rate = domain->h_rate;

  if (lmp->kokkos != nullptr && !kokkosable)
    error->all(FLERR, "KOKKOS package requires a kokkos enabled atom_style");
}

static constexpr bigint DELTA = 16384;

/* ----------------------------------------------------------------------
   roundup N so it is a multiple of DELTA
   error if N exceeds 32-bit int, since will be used as arg to grow()
------------------------------------------------------------------------- */

bigint AtomVec::roundup(bigint n)
{
  if (n % DELTA) n = n / DELTA * DELTA + DELTA;
  if (n > MAXSMALLINT) error->one(FLERR, "Too many atoms created on one or more procs");
  return n;
}

/* ----------------------------------------------------------------------
   grow nmax so it is a multiple of DELTA
------------------------------------------------------------------------- */

void AtomVec::grow_nmax()
{
  nmax = nmax / DELTA * DELTA;
  nmax += DELTA;
}

static constexpr bigint DELTA_BONUS = 8192;

/* ----------------------------------------------------------------------
   grow nmax_bonus so it is a multiple of DELTA_BONUS
------------------------------------------------------------------------- */

int AtomVec::grow_nmax_bonus(int nmax_bonus)
{
  nmax_bonus = nmax_bonus / DELTA_BONUS * DELTA_BONUS;
  nmax_bonus += DELTA_BONUS;
  return nmax_bonus;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVec::grow(int n)
{
  int datatype, cols, maxcols;
  void *pdata;

  if (n == 0)
    grow_nmax();
  else
    nmax = MAX(n, nmax);
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT) error->one(FLERR, "Per-processor system is too big");

  tag = memory->grow(atom->tag, nmax, "atom:tag");
  type = memory->grow(atom->type, nmax, "atom:type");
  mask = memory->grow(atom->mask, nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");
  x = memory->grow(atom->x, nmax, 3, "atom:x");
  v = memory->grow(atom->v, nmax, 3, "atom:v");
  f = memory->grow(atom->f, nmax * comm->nthreads, 3, "atom:f");

  for (int i = 0; i < ngrow; i++) {
    pdata = mgrow.pdata[i];
    datatype = mgrow.datatype[i];
    cols = mgrow.cols[i];
    const int nthreads = threads[i] ? comm->nthreads : 1;
    if (datatype == Atom::DOUBLE) {
      if (cols == 0)
        memory->grow(*((double **) pdata), nmax * nthreads, "atom:dvec");
      else if (cols > 0)
        memory->grow(*((double ***) pdata), nmax * nthreads, cols, "atom:darray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((double ***) pdata), nmax * nthreads, maxcols, "atom:darray");
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0)
        memory->grow(*((int **) pdata), nmax * nthreads, "atom:ivec");
      else if (cols > 0)
        memory->grow(*((int ***) pdata), nmax * nthreads, cols, "atom:iarray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((int ***) pdata), nmax * nthreads, maxcols, "atom:iarray");
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0)
        memory->grow(*((bigint **) pdata), nmax * nthreads, "atom:bvec");
      else if (cols > 0)
        memory->grow(*((bigint ***) pdata), nmax * nthreads, cols, "atom:barray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((bigint ***) pdata), nmax * nthreads, maxcols, "atom:barray");
      }
    }
  }

  for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
    modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVec::copy(int i, int j, int delflag)
{
  int m, n, datatype, cols, collength, ncols;
  void *pdata, *plength;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  if (ncopy) {
    for (n = 0; n < ncopy; n++) {
      pdata = mcopy.pdata[n];
      datatype = mcopy.datatype[n];
      cols = mcopy.cols[n];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++) array[j][m] = array[i][m];
        } else {
          double **array = *((double ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++) array[j][m] = array[i][m];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++) array[j][m] = array[i][m];
        } else {
          int **array = *((int ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++) array[j][m] = array[i][m];
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++) array[j][m] = array[i][m];
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++) array[j][m] = array[i][m];
        }
      }
    }
  }

  if (bonus_flag) copy_bonus(i, j, delflag);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j, delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m, mm, nn, datatype, cols;
  double dx, dy, dz;
  void *pdata;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
    }
  }

  if (ncomm) {
    for (nn = 0; nn < ncomm; nn++) {
      pdata = mcomm.pdata[nn];
      datatype = mcomm.datatype[nn];
      cols = mcomm.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = vec[j];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_comm_bonus(n, list, &buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_comm_vel(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m, mm, nn, datatype, cols;
  double dx, dy, dz, dvx, dvy, dvz;
  void *pdata;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (ncomm_vel) {
    for (nn = 0; nn < ncomm_vel; nn++) {
      pdata = mcomm_vel.pdata[nn];
      datatype = mcomm_vel.datatype[nn];
      cols = mcomm_vel.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = vec[j];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_comm_bonus(n, list, &buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_comm(int n, int first, double *buf)
{
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }

  if (ncomm) {
    for (nn = 0; nn < ncomm; nn++) {
      pdata = mcomm.pdata[nn];
      datatype = mcomm.datatype[nn];
      cols = mcomm.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = buf[m++];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) vec[i] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) unpack_comm_bonus(n, first, &buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_comm_vel(int n, int first, double *buf)
{
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (ncomm_vel) {
    for (nn = 0; nn < ncomm_vel; nn++) {
      pdata = mcomm_vel.pdata[nn];
      datatype = mcomm_vel.datatype[nn];
      cols = mcomm_vel.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = buf[m++];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) vec[i] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) unpack_comm_bonus(n, first, &buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_reverse(int n, int first, double *buf)
{
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }

  if (nreverse) {
    for (nn = 0; nn < nreverse; nn++) {
      pdata = mreverse.pdata[nn];
      datatype = mreverse.datatype[nn];
      cols = mreverse.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) { buf[m++] = vec[i]; }
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++) buf[m++] = array[i][mm];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) { buf[m++] = ubuf(vec[i]).d; }
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) { buf[m++] = ubuf(vec[i]).d; }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
          }
        }
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_reverse(int n, int *list, double *buf)
{
  int i, j, m, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }

  if (nreverse) {
    for (nn = 0; nn < nreverse; nn++) {
      pdata = mreverse.pdata[nn];
      datatype = mreverse.datatype[nn];
      cols = mreverse.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            vec[j] += buf[m++];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) array[j][mm] += buf[m++];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            vec[j] += (int) ubuf(buf[m++]).i;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) array[j][mm] += (int) ubuf(buf[m++]).i;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            vec[j] += (bigint) ubuf(buf[m++]).i;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) array[j][mm] += (bigint) ubuf(buf[m++]).i;
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m, mm, nn, datatype, cols;
  double dx, dy, dz;
  void *pdata;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    }
  }

  if (nborder) {
    for (nn = 0; nn < nborder; nn++) {
      pdata = mborder.pdata[nn];
      datatype = mborder.datatype[nn];
      cols = mborder.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = vec[j];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_border_bonus(n, list, &buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n, list, &buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_border_vel(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, m, mm, nn, datatype, cols;
  double dx, dy, dz, dvx, dvy, dvz;
  void *pdata;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
      }
    } else {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
      }
    }
  }

  if (nborder_vel) {
    for (nn = 0; nn < nborder_vel; nn++) {
      pdata = mborder_vel.pdata[nn];
      datatype = mborder_vel.datatype[nn];
      cols = mborder_vel.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = vec[j];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            buf[m++] = ubuf(vec[j]).d;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_border_bonus(n, list, &buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n, list, &buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_border(int n, int first, double *buf)
{
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  last = first + n;
  while (last > nmax) grow(0);

  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
  }

  if (nborder) {
    for (nn = 0; nn < nborder; nn++) {
      pdata = mborder.pdata[nn];
      datatype = mborder.datatype[nn];
      cols = mborder.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = buf[m++];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) vec[i] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_border_bonus(n, first, &buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->unpack_border(n, first, &buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_border_vel(int n, int first, double *buf)
{
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;

  m = 0;
  last = first + n;
  while (last > nmax) grow(0);

  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }

  if (nborder_vel) {
    for (nn = 0; nn < nborder_vel; nn++) {
      pdata = mborder_vel.pdata[nn];
      datatype = mborder_vel.datatype[nn];
      cols = mborder_vel.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = buf[m++];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) vec[i] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++) array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_border_bonus(n, first, &buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->unpack_border(n, first, &buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVec::pack_exchange(int i, double *buf)
{
  int mm, nn, datatype, cols, collength, ncols;
  void *pdata, *plength;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  if (nexchange) {
    for (nn = 0; nn < nexchange; nn++) {
      pdata = mexchange.pdata[nn];
      datatype = mexchange.datatype[nn];
      cols = mexchange.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[m++] = vec[i];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (mm = 0; mm < cols; mm++) buf[m++] = array[i][mm];
        } else {
          double **array = *((double ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++) buf[m++] = array[i][mm];
        }
      }
      if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[m++] = ubuf(vec[i]).d;
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
        } else {
          int **array = *((int ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++) buf[m++] = ubuf(array[i][mm]).d;
        }
      }
      if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[m++] = ubuf(vec[i]).d;
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata);
          for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++) buf[m++] = ubuf(array[i][mm]).d;
        }
      }
    }
  }

  if (bonus_flag) m += pack_exchange_bonus(i, &buf[m]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::unpack_exchange(double *buf)
{
  int mm, nn, datatype, cols, collength, ncols;
  void *pdata, *plength;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  if (nexchange) {
    for (nn = 0; nn < nexchange; nn++) {
      pdata = mexchange.pdata[nn];
      datatype = mexchange.datatype[nn];
      cols = mexchange.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[nlocal] = buf[m++];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (mm = 0; mm < cols; mm++) array[nlocal][mm] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[nlocal][collength - 1];
          else
            ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = buf[m++];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[nlocal] = (int) ubuf(buf[m++]).i;
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (mm = 0; mm < cols; mm++) array[nlocal][mm] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[nlocal][collength - 1];
          else
            ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[nlocal] = (bigint) ubuf(buf[m++]).i;
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata);
          for (mm = 0; mm < cols; mm++) array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength)
            ncols = (*((int ***) plength))[nlocal][collength - 1];
          else
            ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_exchange_bonus(nlocal, &buf[m]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->unpack_exchange(nlocal, &buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVec::size_restart()
{
  int i, nn, cols, collength, ncols;
  void *plength;

  // NOTE: need to worry about overflow of returned int N

  int nlocal = atom->nlocal;

  // 11 = length storage + id,type,mask,image,x,v

  int n = 11 * nlocal;

  if (nrestart) {
    for (nn = 0; nn < nrestart; nn++) {
      cols = mrestart.cols[nn];
      if (cols == 0)
        n += nlocal;
      else if (cols > 0)
        n += cols * nlocal;
      else {
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        for (i = 0; i < nlocal; i++) {
          if (collength)
            ncols = (*((int ***) plength))[i][collength - 1];
          else
            ncols = (*((int **) plength))[i];
          n += ncols;
        }
      }
    }
  }

  if (bonus_flag) n += size_restart_bonus();

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++) n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVec::pack_restart(int i, double *buf)
{
  int mm, nn, datatype, cols, collength, ncols;
  void *pdata, *plength;

  // if needed, change values before packing

  pack_restart_pre(i);

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  for (nn = 0; nn < nrestart; nn++) {
    pdata = mrestart.pdata[nn];
    datatype = mrestart.datatype[nn];
    cols = mrestart.cols[nn];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        buf[m++] = vec[i];
      } else if (cols > 0) {
        double **array = *((double ***) pdata);
        for (mm = 0; mm < cols; mm++) buf[m++] = array[i][mm];
      } else {
        double **array = *((double ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[i][collength - 1];
        else
          ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++) buf[m++] = array[i][mm];
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        buf[m++] = ubuf(vec[i]).d;
      } else if (cols > 0) {
        int **array = *((int ***) pdata);
        for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
      } else {
        int **array = *((int ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[i][collength - 1];
        else
          ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++) buf[m++] = ubuf(array[i][mm]).d;
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        buf[m++] = ubuf(vec[i]).d;
      } else if (cols > 0) {
        bigint **array = *((bigint ***) pdata);
        for (mm = 0; mm < cols; mm++) buf[m++] = ubuf(array[i][mm]).d;
      } else {
        bigint **array = *((bigint ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[i][collength - 1];
        else
          ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++) buf[m++] = ubuf(array[i][mm]).d;
      }
    }
  }

  if (bonus_flag) m += pack_restart_bonus(i, &buf[m]);

  // if needed, restore values after packing

  pack_restart_post(i);

  // invoke fixes which store peratom restart info

  for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
    m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVec::unpack_restart(double *buf)
{
  int mm, nn, datatype, cols, collength, ncols;
  void *pdata, *plength;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store) memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  for (nn = 0; nn < nrestart; nn++) {
    pdata = mrestart.pdata[nn];
    datatype = mrestart.datatype[nn];
    cols = mrestart.cols[nn];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = buf[m++];
      } else if (cols > 0) {
        double **array = *((double ***) pdata);
        for (mm = 0; mm < cols; mm++) array[nlocal][mm] = buf[m++];
      } else {
        double **array = *((double ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[nlocal][collength - 1];
        else
          ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = buf[m++];
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        vec[nlocal] = (int) ubuf(buf[m++]).i;
      } else if (cols > 0) {
        int **array = *((int ***) pdata);
        for (mm = 0; mm < cols; mm++) array[nlocal][mm] = (int) ubuf(buf[m++]).i;
      } else {
        int **array = *((int ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[nlocal][collength - 1];
        else
          ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = (int) ubuf(buf[m++]).i;
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = (bigint) ubuf(buf[m++]).i;
      } else if (cols > 0) {
        bigint **array = *((bigint ***) pdata);
        for (mm = 0; mm < cols; mm++) array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
      } else {
        bigint **array = *((bigint ***) pdata);
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        if (collength)
          ncols = (*((int ***) plength))[nlocal][collength - 1];
        else
          ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++) array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
      }
    }
  }

  if (bonus_flag) m += unpack_restart_bonus(nlocal, &buf[m]);

  // if needed, initialize other peratom values

  unpack_restart_init(nlocal);

  // store extra restart info which fixes can unpack when instantiated

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int>(buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVec::create_atom(int itype, double *coord)
{
  int m, n, datatype, cols;
  void *pdata;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  // initialization additional fields

  for (n = 0; n < ncreate; n++) {
    pdata = mcreate.pdata[n];
    datatype = mcreate.datatype[n];
    cols = mcreate.cols[n];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = 0.0;
      } else {
        double **array = *((double ***) pdata);
        for (m = 0; m < cols; m++) array[nlocal][m] = 0.0;
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        vec[nlocal] = 0;
      } else {
        int **array = *((int ***) pdata);
        for (m = 0; m < cols; m++) array[nlocal][m] = 0;
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = 0;
      } else {
        bigint **array = *((bigint ***) pdata);
        for (m = 0; m < cols; m++) array[nlocal][m] = 0;
      }
    }
  }

  // if needed, initialize non-zero peratom values

  create_atom_post(nlocal);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other peratom quantities
------------------------------------------------------------------------- */

void AtomVec::data_atom(double *coord, imageint imagetmp, const std::vector<std::string> &values,
                        std::string &extract)
{
  int m, n, datatype, cols;
  void *pdata;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = imagetmp;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  int ivalue = 0;
  for (n = 0; n < ndata_atom; n++) {
    pdata = mdata_atom.pdata[n];
    datatype = mdata_atom.datatype[n];
    cols = mdata_atom.cols[n];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = utils::numeric(FLERR, values[ivalue++], true, lmp);
      } else {
        double **array = *((double ***) pdata);
        if (array == atom->x) {    // x was already set by coord arg
          ivalue += cols;
          continue;
        }
        for (m = 0; m < cols; m++)                                      \
          array[nlocal][m] = utils::numeric(FLERR, values[ivalue++], true, lmp);
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        if (vec == atom->type) {    // custom treatment of atom types
          extract = values[ivalue++];
          continue;
        }
        vec[nlocal] = utils::inumeric(FLERR, values[ivalue++], true, lmp);
      } else {
        int **array = *((int ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = utils::inumeric(FLERR, values[ivalue++], true, lmp);
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = utils::bnumeric(FLERR, values[ivalue++], true, lmp);
      } else {
        bigint **array = *((bigint ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = utils::bnumeric(FLERR, values[ivalue++], true, lmp);
      }
    }
  }

  // error checks applicable to all styles

  if ((atom->tag_enable && (tag[nlocal] <= 0)) || (!atom->tag_enable && (tag[nlocal] != 0)))
    error->one(FLERR, "Invalid atom ID {} in line {} of Atoms section of data file", tag[nlocal],
               nlocal + 1);

  // if needed, modify unpacked values or initialize other peratom values

  data_atom_post(nlocal);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVec::pack_data(double **buf)
{
  int i, j, m, n, datatype, cols;
  void *pdata;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {

    // if needed, change values before packing

    pack_data_pre(i);

    j = 0;
    for (n = 0; n < ndata_atom; n++) {
      pdata = mdata_atom.pdata[n];
      datatype = mdata_atom.datatype[n];
      cols = mdata_atom.cols[n];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[i][j++] = vec[i];
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = array[i][m];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = ubuf(array[i][m]).d;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = ubuf(array[i][m]).d;
        }
      }
    }

    buf[i][j++] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][j++] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][j++] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;

    // if needed, restore values after packing

    pack_data_post(i);
  }
}

/* ----------------------------------------------------------------------
   write atom info to data file
   id is first field, 3 image flags are final fields
------------------------------------------------------------------------- */

void AtomVec::write_data(FILE *fp, int n, double **buf)
{
  int i, j, m, nn, datatype, cols;

  for (i = 0; i < n; i++) {
    fmt::print(fp, "{}", ubuf(buf[i][0]).i);

    j = 1;
    for (nn = 1; nn < ndata_atom; nn++) {
      datatype = mdata_atom.datatype[nn];
      cols = mdata_atom.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          fmt::print(fp, " {:.16}", buf[i][j++]);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", buf[i][j++]);
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          if (atom->types_style == Atom::LABELS &&
              atom->peratom[mdata_atom.index[nn]].name == "type") {
            fmt::print(fp, " {}", atom->lmap->typelabel[ubuf(buf[i][j++]).i - 1]);
          } else
            fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        }
      }
    }

    fmt::print(fp, " {} {} {}\n", ubuf(buf[i][j]).i, ubuf(buf[i][j + 1]).i, ubuf(buf[i][j + 2]).i);
  }
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVec::data_vel(int ilocal, const std::vector<std::string> &values)
{
  int m, n, datatype, cols;
  void *pdata;

  double **v = atom->v;
  int ivalue = 1;
  v[ilocal][0] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  v[ilocal][1] = utils::numeric(FLERR, values[ivalue++], true, lmp);
  v[ilocal][2] = utils::numeric(FLERR, values[ivalue++], true, lmp);

  if (ndata_vel > 2) {
    for (n = 2; n < ndata_vel; n++) {
      pdata = mdata_vel.pdata[n];
      datatype = mdata_vel.datatype[n];
      cols = mdata_vel.cols[n];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[ilocal] = utils::numeric(FLERR, values[ivalue++], true, lmp);
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::numeric(FLERR, values[ivalue++], true, lmp);
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[ilocal] = utils::inumeric(FLERR, values[ivalue++], true, lmp);
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::inumeric(FLERR, values[ivalue++], true, lmp);
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[ilocal] = utils::bnumeric(FLERR, values[ivalue++], true, lmp);
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::bnumeric(FLERR, values[ivalue++], true, lmp);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVec::pack_vel(double **buf)
{
  int i, j, m, n, datatype, cols;
  void *pdata;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    j = 0;
    for (n = 0; n < ndata_vel; n++) {
      pdata = mdata_vel.pdata[n];
      datatype = mdata_vel.datatype[n];
      cols = mdata_vel.cols[n];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[i][j++] = vec[i];
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = array[i][m];
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = ubuf(array[i][m]).d;
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++) buf[i][j++] = ubuf(array[i][m]).d;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   write velocity info to data file
   id and velocity vector are first 4 fields
------------------------------------------------------------------------- */

void AtomVec::write_vel(FILE *fp, int n, double **buf)
{
  int i, j, m, nn, datatype, cols;

  for (i = 0; i < n; i++) {
    fmt::print(fp, "{}", ubuf(buf[i][0]).i);

    j = 1;
    for (nn = 1; nn < ndata_vel; nn++) {
      datatype = mdata_vel.datatype[nn];
      cols = mdata_vel.cols[nn];
      if (datatype == Atom::DOUBLE) {
        if (cols == 0) {
          fmt::print(fp, " {}", buf[i][j++]);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", buf[i][j++]);
        }
      } else if (datatype == Atom::INT) {
        if (cols == 0) {
          fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        }
      } else if (datatype == Atom::BIGINT) {
        if (cols == 0) {
          fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        } else {
          for (m = 0; m < cols; m++) fmt::print(fp, " {}", ubuf(buf[i][j++]).i);
        }
      }
    }
    fputs("\n", fp);
  }
}

/* ----------------------------------------------------------------------
   pack bond info for data file into buf if non-nullptr
   return count of bonds from this proc
   do not count/pack bonds with bondtype = 0
   if bondtype is negative, flip back to positive
------------------------------------------------------------------------- */

int AtomVec::pack_bond(tagint **buf)
{
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i, j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++) {
        if (bond_type[i][j] == 0) continue;
        if (buf) {
          buf[m][0] = MAX(bond_type[i][j], -bond_type[i][j]);
          buf[m][1] = tag[i];
          buf[m][2] = bond_atom[i][j];
        }
        m++;
      }
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++)
        if (tag[i] < bond_atom[i][j]) {
          if (bond_type[i][j] == 0) continue;
          if (buf) {
            buf[m][0] = MAX(bond_type[i][j], -bond_type[i][j]);
            buf[m][1] = tag[i];
            buf[m][2] = bond_atom[i][j];
          }
          m++;
        }
  }

  return m;
}

/* ----------------------------------------------------------------------
   write bond info to data file
------------------------------------------------------------------------- */

void AtomVec::write_bond(FILE *fp, int n, tagint **buf, int index)
{
  std::string typestr;
  for (int i = 0; i < n; i++) {
    typestr = std::to_string(buf[i][0]);
    if (atom->types_style == Atom::LABELS) typestr = atom->lmap->btypelabel[buf[i][0] - 1];
    fmt::print(fp, "{} {} {} {}\n", index, typestr, buf[i][1], buf[i][2]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   pack angle info for data file into buf if non-nullptr
   return count of angles from this proc
   do not count/pack angles with angletype = 0
   if angletype is negative, flip back to positive
------------------------------------------------------------------------- */

int AtomVec::pack_angle(tagint **buf)
{
  tagint *tag = atom->tag;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i, j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_angle[i]; j++) {
        if (angle_type[i][j] == 0) continue;
        if (buf) {
          buf[m][0] = MAX(angle_type[i][j], -angle_type[i][j]);
          buf[m][1] = angle_atom1[i][j];
          buf[m][2] = angle_atom2[i][j];
          buf[m][3] = angle_atom3[i][j];
        }
        m++;
      }
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_angle[i]; j++)
        if (tag[i] == angle_atom2[i][j]) {
          if (angle_type[i][j] == 0) continue;
          if (buf) {
            buf[m][0] = MAX(angle_type[i][j], -angle_type[i][j]);
            buf[m][1] = angle_atom1[i][j];
            buf[m][2] = angle_atom2[i][j];
            buf[m][3] = angle_atom3[i][j];
          }
          m++;
        }
  }

  return m;
}

/* ----------------------------------------------------------------------
   write angle info to data file
------------------------------------------------------------------------- */

void AtomVec::write_angle(FILE *fp, int n, tagint **buf, int index)
{
  std::string typestr;
  for (int i = 0; i < n; i++) {
    typestr = std::to_string(buf[i][0]);
    if (atom->types_style == Atom::LABELS) typestr = atom->lmap->atypelabel[buf[i][0] - 1];
    fmt::print(fp, "{} {} {} {} {}\n", index, typestr, buf[i][1], buf[i][2], buf[i][3]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   pack dihedral info for data file
------------------------------------------------------------------------- */

int AtomVec::pack_dihedral(tagint **buf)
{
  tagint *tag = atom->tag;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_type = atom->dihedral_type;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i, j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++) {
        if (buf) {
          buf[m][0] = MAX(dihedral_type[i][j], -dihedral_type[i][j]);
          buf[m][1] = dihedral_atom1[i][j];
          buf[m][2] = dihedral_atom2[i][j];
          buf[m][3] = dihedral_atom3[i][j];
          buf[m][4] = dihedral_atom4[i][j];
        }
        m++;
      }
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++)
        if (tag[i] == dihedral_atom2[i][j]) {
          if (buf) {
            buf[m][0] = MAX(dihedral_type[i][j], -dihedral_type[i][j]);
            buf[m][1] = dihedral_atom1[i][j];
            buf[m][2] = dihedral_atom2[i][j];
            buf[m][3] = dihedral_atom3[i][j];
            buf[m][4] = dihedral_atom4[i][j];
          }
          m++;
        }
  }

  return m;
}

/* ----------------------------------------------------------------------
   write dihedral info to data file
------------------------------------------------------------------------- */

void AtomVec::write_dihedral(FILE *fp, int n, tagint **buf, int index)
{
  std::string typestr;
  for (int i = 0; i < n; i++) {
    typestr = std::to_string(buf[i][0]);
    if (atom->types_style == Atom::LABELS) typestr = atom->lmap->dtypelabel[buf[i][0] - 1];
    fmt::print(fp, "{} {} {} {} {} {}\n", index, typestr, buf[i][1], buf[i][2], buf[i][3],
               buf[i][4]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   pack improper info for data file
------------------------------------------------------------------------- */

int AtomVec::pack_improper(tagint **buf)
{
  tagint *tag = atom->tag;
  int *num_improper = atom->num_improper;
  int **improper_type = atom->improper_type;
  tagint **improper_atom1 = atom->improper_atom1;
  tagint **improper_atom2 = atom->improper_atom2;
  tagint **improper_atom3 = atom->improper_atom3;
  tagint **improper_atom4 = atom->improper_atom4;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i, j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_improper[i]; j++) {
        if (buf) {
          buf[m][0] = MAX(improper_type[i][j], -improper_type[i][j]);
          buf[m][1] = improper_atom1[i][j];
          buf[m][2] = improper_atom2[i][j];
          buf[m][3] = improper_atom3[i][j];
          buf[m][4] = improper_atom4[i][j];
        }
        m++;
      }
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_improper[i]; j++)
        if (tag[i] == improper_atom2[i][j]) {
          if (buf) {
            buf[m][0] = MAX(improper_type[i][j], -improper_type[i][j]);
            buf[m][1] = improper_atom1[i][j];
            buf[m][2] = improper_atom2[i][j];
            buf[m][3] = improper_atom3[i][j];
            buf[m][4] = improper_atom4[i][j];
          }
          m++;
        }
  }

  return m;
}

/* ----------------------------------------------------------------------
   write improper info to data file
------------------------------------------------------------------------- */

void AtomVec::write_improper(FILE *fp, int n, tagint **buf, int index)
{
  std::string typestr;
  for (int i = 0; i < n; i++) {
    typestr = std::to_string(buf[i][0]);
    if (atom->types_style == Atom::LABELS) typestr = atom->lmap->itypelabel[buf[i][0] - 1];
    fmt::print(fp, "{} {} {} {} {} {}\n", index, typestr, buf[i][1], buf[i][2], buf[i][3],
               buf[i][4]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   convert info input by read_data from general to restricted triclinic
   atom coords are converted in Atom::data_atoms()
   parent class operates on data from Velocities section of data file
   child classes operate on all other data: Atoms, Ellipsoids, Lines, Triangles, etc
------------------------------------------------------------------------- */

void AtomVec::read_data_general_to_restricted(int nlocal_previous, int nlocal)
{
  int datatype, cols;
  void *pdata;

  for (int n = 1; n < ndata_vel; n++) {
    pdata = mdata_vel.pdata[n];
    datatype = mdata_vel.datatype[n];
    cols = mdata_vel.cols[n];

    // operate on v, omega, angmom
    // no other read_data Velocities fields are Nx3 double arrays

    if (datatype == Atom::DOUBLE) {
      if (cols == 3) {
        double **array = *((double ***) pdata);
        for (int i = nlocal_previous; i < nlocal; i++)
          domain->general_to_restricted_vector(array[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   convert info output by write_data from restricted to general triclinic
   create "hold" copy of original restricted data to restore after data file is written
   parent class only operates on x and data from Velocities section of data file
   child classes operate on all other data: Atoms, Ellipsoids, Lines, Triangles, etc
------------------------------------------------------------------------- */

void AtomVec::write_data_restricted_to_general()
{
  int datatype, cols;
  void *pdata;

  int nlocal = atom->nlocal;

  memory->create(x_hold,nlocal,3,"atomvec:x_hold");
  if (nlocal) memcpy(&x_hold[0][0],&x[0][0],3*nlocal*sizeof(double));
  for (int i = 0; i < nlocal; i++)
    domain->restricted_to_general_coords(x[i]);

  double **omega = atom->omega;
  double **angmom = atom->angmom;

  for (int n = 1; n < ndata_vel; n++) {
    pdata = mdata_vel.pdata[n];
    datatype = mdata_vel.datatype[n];
    cols = mdata_vel.cols[n];

    // operate on v, omega, angmom
    // no other write_data Velocities fields are Nx3 double arrays

    if (datatype == Atom::DOUBLE) {
      if (cols == 3) {
        double **array = *((double ***) pdata);

        if (array == v) {
          memory->create(v_hold,nlocal,3,"atomvec:v_hold");
          if (nlocal) memcpy(&v_hold[0][0],&v[0][0],3*nlocal*sizeof(double));
          for (int i = 0; i < nlocal; i++)
            domain->restricted_to_general_vector(v[i]);
        } else if (array == omega) {
          memory->create(omega_hold,nlocal,3,"atomvec:omega_hold");
          if (nlocal) memcpy(&omega_hold[0][0],&omega[0][0],3*nlocal*sizeof(double));
          for (int i = 0; i < nlocal; i++)
            domain->restricted_to_general_vector(omega[i]);
        } else if (array == angmom) {
          memory->create(angmom_hold,nlocal,3,"atomvec:angmom_hold");
          if (nlocal) memcpy(&angmom_hold[0][0],&angmom[0][0],3*nlocal*sizeof(double));
          for (int i = 0; i < nlocal; i++)
            domain->restricted_to_general_vector(angmom[i]);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   restore info output by write_data to restricted triclinic
   original data is in "hold" arrays
   parent class only operates on x and data from Velocities section of data file
   child classes operate on all other data: Atoms, Ellipsoids, Lines, Triangles, etc
------------------------------------------------------------------------- */

void AtomVec::write_data_restore_restricted()
{
  int nlocal = atom->nlocal;

  if (x_hold) {
    memcpy(&x[0][0],&x_hold[0][0],3*nlocal*sizeof(double));
    memory->destroy(x_hold);
    x_hold = nullptr;
  }

  // operate on v, omega, angmom
  // no other write_data Velocities fields are Nx3 double arrays

  if (v_hold) {
    memcpy(&v[0][0],&v_hold[0][0],3*nlocal*sizeof(double));
    memory->destroy(v_hold);
    v_hold = nullptr;
  }

  if (omega_hold) {
    memcpy(&atom->omega[0][0],&omega_hold[0][0],3*nlocal*sizeof(double));
    memory->destroy(omega_hold);
    omega_hold = nullptr;
  }

  if (angmom_hold) {
    memcpy(&atom->angmom[0][0],&angmom_hold[0][0],3*nlocal*sizeof(double));
    memory->destroy(angmom_hold);
    angmom_hold = nullptr;
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double AtomVec::memory_usage()
{
  int datatype, cols, maxcols;
  void *pdata;

  double bytes = 0;

  bytes += memory->usage(tag, nmax);
  bytes += memory->usage(type, nmax);
  bytes += memory->usage(mask, nmax);
  bytes += memory->usage(image, nmax);
  bytes += memory->usage(x, nmax, 3);
  bytes += memory->usage(v, nmax, 3);
  bytes += memory->usage(f, nmax * comm->nthreads, 3);

  for (int i = 0; i < ngrow; i++) {
    pdata = mgrow.pdata[i];
    datatype = mgrow.datatype[i];
    cols = mgrow.cols[i];
    const int nthreads = threads[i] ? comm->nthreads : 1;
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        bytes += memory->usage(*((double **) pdata), nmax * nthreads);
      } else if (cols > 0) {
        bytes += memory->usage(*((double ***) pdata), nmax * nthreads, cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((double ***) pdata), nmax * nthreads, maxcols);
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        bytes += memory->usage(*((int **) pdata), nmax * nthreads);
      } else if (cols > 0) {
        bytes += memory->usage(*((int ***) pdata), nmax * nthreads, cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((int ***) pdata), nmax * nthreads, maxcols);
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bytes += memory->usage(*((bigint **) pdata), nmax * nthreads);
      } else if (cols > 0) {
        bytes += memory->usage(*((bigint ***) pdata), nmax * nthreads, cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((bigint ***) pdata), nmax * nthreads, maxcols);
      }
    }
  }

  if (bonus_flag) bytes += memory_usage_bonus();

  return bytes;
}

// ----------------------------------------------------------------------
// internal methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVec::setup_fields()
{
  int n, cols;

  if ((fields_data_atom.size() < 1) || (fields_data_atom[0] != "id"))
    error->all(FLERR, "Atom style fields_data_atom must have 'id' as first field");
  if ((fields_data_vel.size() < 2) || (fields_data_vel[0] != "id") || (fields_data_vel[1] != "v"))
    error->all(FLERR, "Atom style fields_data_vel must have 'id' and 'v' as first two fields");

  // process field strings
  // return # of fields and matching index into atom.peratom (in Method struct)

  ngrow = process_fields(fields_grow, default_grow, &mgrow);
  ncopy = process_fields(fields_copy, default_copy, &mcopy);
  ncomm = process_fields(fields_comm, default_comm, &mcomm);
  ncomm_vel = process_fields(fields_comm_vel, default_comm_vel, &mcomm_vel);
  nreverse = process_fields(fields_reverse, default_reverse, &mreverse);
  nborder = process_fields(fields_border, default_border, &mborder);
  nborder_vel = process_fields(fields_border_vel, default_border_vel, &mborder_vel);
  nexchange = process_fields(fields_exchange, default_exchange, &mexchange);
  nrestart = process_fields(fields_restart, default_restart, &mrestart);
  ncreate = process_fields(fields_create, default_create, &mcreate);
  ndata_atom = process_fields(fields_data_atom, default_data_atom, &mdata_atom);
  ndata_vel = process_fields(fields_data_vel, default_data_vel, &mdata_vel);

  // populate field-based data struct for each method to use

  init_method(ngrow, &mgrow);
  init_method(ncopy, &mcopy);
  init_method(ncomm, &mcomm);
  init_method(ncomm_vel, &mcomm_vel);
  init_method(nreverse, &mreverse);
  init_method(nborder, &mborder);
  init_method(nborder_vel, &mborder_vel);
  init_method(nexchange, &mexchange);
  init_method(nrestart, &mrestart);
  init_method(ncreate, &mcreate);
  init_method(ndata_atom, &mdata_atom);
  init_method(ndata_vel, &mdata_vel);

  // create threads data struct for grow and memory_usage to use

  if (ngrow)
    threads = new bool[ngrow];
  else
    threads = nullptr;
  for (int i = 0; i < ngrow; i++) {
    const auto &field = atom->peratom[mgrow.index[i]];
    threads[i] = field.threadflag == 1;
  }

  // set style-specific sizes

  comm_x_only = 1;
  if (ncomm) comm_x_only = 0;
  if (bonus_flag && size_forward_bonus) comm_x_only = 0;

  if (nreverse == 0)
    comm_f_only = 1;
  else
    comm_f_only = 0;

  size_forward = 3;
  for (n = 0; n < ncomm; n++) {
    cols = mcomm.cols[n];
    if (cols == 0)
      size_forward++;
    else
      size_forward += cols;
  }
  if (bonus_flag) size_forward += size_forward_bonus;

  size_reverse = 3;
  for (n = 0; n < nreverse; n++) {
    cols = mreverse.cols[n];
    if (cols == 0)
      size_reverse++;
    else
      size_reverse += cols;
  }

  size_border = 6;
  for (n = 0; n < nborder; n++) {
    cols = mborder.cols[n];
    if (cols == 0)
      size_border++;
    else
      size_border += cols;
  }
  if (bonus_flag) size_border += size_border_bonus;

  size_velocity = 3;
  for (n = 0; n < ncomm_vel; n++) {
    cols = mcomm_vel.cols[n];
    if (cols == 0)
      size_velocity++;
    else
      size_velocity += cols;
  }

  size_data_atom = 0;
  for (n = 0; n < ndata_atom; n++) {
    cols = mdata_atom.cols[n];
    if (atom->peratom[mdata_atom.index[n]].name == "x") xcol_data = size_data_atom + 1;
    if (cols == 0)
      size_data_atom++;
    else
      size_data_atom += cols;
  }

  size_data_vel = 0;
  for (n = 0; n < ndata_vel; n++) {
    cols = mdata_vel.cols[n];
    if (cols == 0)
      size_data_vel++;
    else
      size_data_vel += cols;
  }
}

/* ----------------------------------------------------------------------
   process a single field string
------------------------------------------------------------------------- */

int AtomVec::process_fields(const std::vector<std::string> &words,
                            const std::vector<std::string> &def_words, Method *method)
{
  int nfield = words.size();
  int ndef = def_words.size();

  // process fields one by one, add to index vector

  const auto &peratom = atom->peratom;
  const int nperatom = peratom.size();

  // allocate memory in method
  method->resize(nfield);

  std::vector<int> &index = method->index;
  int match;

  for (int i = 0; i < nfield; i++) {
    const std::string &field = words[i];

    // find field in master Atom::peratom list

    for (match = 0; match < nperatom; match++)
      if (field == peratom[match].name) break;
    if (match == nperatom) error->all(FLERR, "Peratom field {} not recognized", field);
    index[i] = match;

    // error if field appears multiple times

    for (match = 0; match < i; match++)
      if (index[i] == index[match]) error->all(FLERR, "Peratom field {} is repeated", field);

    // error if field is in default str

    for (match = 0; match < ndef; match++)
      if (field == def_words[match]) error->all(FLERR, "Peratom field {} is a default", field);
  }

  return nfield;
}

/* ----------------------------------------------------------------------
   init method data structs for processing fields
------------------------------------------------------------------------- */

void AtomVec::init_method(int nfield, Method *method)
{
  for (int i = 0; i < nfield; i++) {
    const auto &field = atom->peratom[method->index[i]];
    method->pdata[i] = (void *) field.address;
    method->datatype[i] = field.datatype;
    method->cols[i] = field.cols;
    if (method->cols[i] < 0) {
      method->maxcols[i] = field.address_maxcols;
      method->collength[i] = field.collength;
      method->plength[i] = field.address_length;
    }
  }
}

/* ----------------------------------------------------------------------
   Method class members
------------------------------------------------------------------------- */

void AtomVec::Method::resize(int nfield)
{
  pdata.resize(nfield);
  datatype.resize(nfield);
  cols.resize(nfield);
  maxcols.resize(nfield);
  collength.resize(nfield);
  plength.resize(nfield);
  index.resize(nfield);
}
