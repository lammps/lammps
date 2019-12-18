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

#include "atom_vec.h"
#include <cstdio>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 16384
#define DELTA_BONUS 8192

enum{DOUBLE,INT,BIGINT};

/* ---------------------------------------------------------------------- */

AtomVec::AtomVec(LAMMPS *lmp) : Pointers(lmp)
{
  nmax = 0;

  molecular = 0;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 0;
  mass_type = dipole_type = 0;
  forceclearflag = 0;
  maxexchange = 0;
  bonus_flag = 0;
  size_forward_bonus = size_border_bonus = 0;

  kokkosable = 0;

  nargcopy = 0;
  argcopy = NULL;

  nthreads = comm->nthreads;

  // peratom variables auto-included in corresponding child style fields string
  // these fields cannot be specified in the fields string

  default_grow = "id type mask image x v f";
  default_copy = "id type mask image x v";
  default_comm = "x";
  default_comm_vel = "x v";
  default_reverse = "f";
  default_border = "id type mask x";
  default_border_vel = "id type mask x v";
  default_exchange = "id type mask image x v";
  default_restart = "id type mask image x v";
  default_create = "id type mask image x v";
  default_data_atom = "";
  default_data_vel = "";

  // initializations

  init_method(&mgrow);
  init_method(&mcopy);
  init_method(&mcomm);
  init_method(&mcomm_vel);
  init_method(&mreverse);
  init_method(&mborder);
  init_method(&mborder_vel);
  init_method(&mexchange);
  init_method(&mrestart);
  init_method(&mcreate);
  init_method(&mdata_atom);
  init_method(&mdata_vel);
}

/* ---------------------------------------------------------------------- */

AtomVec::~AtomVec()
{
  for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
  delete [] argcopy;

  destroy_method(&mgrow);
  destroy_method(&mcopy);
  destroy_method(&mcomm);
  destroy_method(&mcomm_vel);
  destroy_method(&mreverse);
  destroy_method(&mborder);
  destroy_method(&mborder_vel);
  destroy_method(&mexchange);
  destroy_method(&mrestart);
  destroy_method(&mcreate);
  destroy_method(&mdata_atom);
  destroy_method(&mdata_vel);

  delete [] threads;
}

/* ----------------------------------------------------------------------
   make copy of args for use by restart & replicate
------------------------------------------------------------------------- */

void AtomVec::store_args(int narg, char **arg)
{
  nargcopy = narg;
  argcopy = new char*[nargcopy];
  for (int i = 0; i < nargcopy; i++) {
    int n = strlen(arg[i]) + 1;
    argcopy[i] = new char[n];
    strcpy(argcopy[i],arg[i]);
  }
}

/* ----------------------------------------------------------------------
   no additional args by default
------------------------------------------------------------------------- */

void AtomVec::process_args(int narg, char ** /*arg*/)
{
  if (narg) error->all(FLERR,"Invalid atom_style command");
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

  if (lmp->kokkos != NULL && !kokkosable)
    error->all(FLERR,"KOKKOS package requires a kokkos enabled atom_style");
}

/* ----------------------------------------------------------------------
   grow nmax so it is a multiple of DELTA
------------------------------------------------------------------------- */

void AtomVec::grow_nmax()
{
  nmax = nmax/DELTA * DELTA;
  nmax += DELTA;
}

/* ----------------------------------------------------------------------
   grow nmax_bonus so it is a multiple of DELTA_BONUS
------------------------------------------------------------------------- */

int AtomVec::grow_nmax_bonus(int nmax_bonus)
{
  nmax_bonus = nmax_bonus/DELTA_BONUS * DELTA_BONUS;
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
  int datatype,cols,maxcols;
  void *pdata,*plength;

  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*nthreads,3,"atom:f");

  for (int i = 0; i < ngrow; i++) {
    pdata = mgrow.pdata[i];
    datatype = mgrow.datatype[i];
    cols = mgrow.cols[i];
    if (datatype == DOUBLE) {
      if (cols == 0)
        memory->grow(*((double **) pdata),nmax*threads[i],"atom:dvec");
      else if (cols > 0)
        memory->grow(*((double ***) pdata),nmax*threads[i],cols,"atom:darray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((double ***) pdata),nmax*threads[i],maxcols,"atom:darray");
      }
    } else if (datatype == INT) {
      if (cols == 0)
        memory->grow(*((int **) pdata),nmax*threads[i],"atom:ivec");
      else if (cols > 0)
        memory->grow(*((int ***) pdata),nmax*threads[i],cols,"atom:iarray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((int ***) pdata),nmax*threads[i],maxcols,"atom:iarray");
      }
    } else if (datatype == BIGINT) {
      if (cols == 0)
        memory->grow(*((bigint **) pdata),nmax*threads[i],"atom:bvec");
      else if (cols > 0)
        memory->grow(*((bigint ***) pdata),nmax*threads[i],cols,"atom:barray");
      else {
        maxcols = *(mgrow.maxcols[i]);
        memory->grow(*((int ***) pdata),nmax*threads[i],maxcols,"atom:barray");
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
  int m,n,datatype,cols,collength,ncols;
  void *pdata,*plength;

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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            array[j][m] = array[i][m];
        } else { 
          double **array = *((double ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++)
            array[j][m] = array[i][m];
       }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            array[j][m] = array[i][m];
        } else {
          int **array = *((int ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++)
            array[j][m] = array[i][m];
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[j] = vec[i];
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            array[j][m] = array[i][m];
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mcopy.collength[n];
          plength = mcopy.plength[n];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (m = 0; m < ncols; m++)
            array[j][m] = array[i][m];
        }
      }
    }
  }

  if (bonus_flag) copy_bonus(i,j,delflag);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m,mm,nn,datatype,cols;
  double dx,dy,dz;
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
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
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
      if (datatype == DOUBLE) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == INT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == BIGINT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_comm_bonus(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_comm_vel(int n, int *list, double *buf,
                           int pbc_flag, int *pbc)
{
  int i,j,m,mm,nn,datatype,cols;
  double dx,dy,dz,dvx,dvy,dvz;
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
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
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
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
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
      if (datatype == DOUBLE) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == INT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == BIGINT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_comm_bonus(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_comm(int n, int first, double *buf)
{
  int i,m,last,mm,nn,datatype,cols;
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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++)
            vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = buf[m++];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++)
            vec[i] = ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++)
            vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) unpack_comm_bonus(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last,mm,nn,datatype,cols;
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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++)
            vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = buf[m++];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++)
            vec[i] = ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++)
            vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) unpack_comm_bonus(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_reverse(int n, int first, double *buf)
{
  int i,m,last,mm,nn,datatype,cols;
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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++) {
            buf[m++] = vec[i];
          }
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++)
              buf[m++] = array[i][mm];
          }
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++) {
            buf[m++] = ubuf(vec[i]).d;
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[i][mm]).d;
          }
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++) {
            buf[m++] = ubuf(vec[i]).d;
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++) {
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[i][mm]).d;
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
  int i,j,m,mm,nn,datatype,cols;
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
      if (datatype == DOUBLE) {
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
            for (mm = 0; mm < cols; mm++)
              array[j][mm] += buf[m++];
          }
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            vec[j] += buf[m++];
          }
        } else {
          int **array = *((int ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++)
              array[j][mm] += buf[m++];
          }
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            vec[j] += buf[m++];
          }
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = 0; i < n; i++) {
            j = list[i];
            for (mm = 0; mm < cols; mm++)
              array[j][mm] += buf[m++];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m,mm,nn,datatype,cols;
  double dx,dy,dz;
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
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
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
      if (datatype == DOUBLE) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == INT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == BIGINT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_border_bonus(n,list,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::pack_border_vel(int n, int *list, double *buf, 
                             int pbc_flag, int *pbc)
{
  int i,j,m,mm,nn,datatype,cols;
  double dx,dy,dz,dvx,dvy,dvz;
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
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
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
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
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
      if (datatype == DOUBLE) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = array[j][mm];
          }
        }
      } else if (datatype == INT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      } else if (datatype == BIGINT) {
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
            for (mm = 0; mm < cols; mm++)
              buf[m++] = ubuf(array[j][mm]).d;
          }
        }
      }
    }
  }

  if (bonus_flag) m += pack_border_bonus(n,list,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_border(int n, int first, double *buf)
{
  int i,m,last,mm,nn,datatype,cols;
  void *pdata;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++)
            vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = buf[m++];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++)
            vec[i] = ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++)
            vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_border_bonus(n,first,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVec::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last,mm,nn,datatype,cols;
  void *pdata;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          for (i = first; i < last; i++)
            vec[i] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = buf[m++];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          for (i = first; i < last; i++)
            vec[i] = ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          for (i = first; i < last; i++)
            vec[i] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (i = first; i < last; i++)
            for (mm = 0; mm < cols; mm++)
              array[i][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_border_bonus(n,first,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVec::pack_exchange(int i, double *buf)
{
  int mm,nn,datatype,cols,collength,ncols;
  void *pdata,*plength;

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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[m++] = vec[i];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (mm = 0; mm < cols; mm++)
            buf[m++] = array[i][mm];
        } else {
          double **array = *((double ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++)
            buf[m++] = array[i][mm];
        }
      } if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[m++] = ubuf(vec[i]).d;
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (mm = 0; mm < cols; mm++)
            buf[m++] = ubuf(array[i][mm]).d;
        } else {
          int **array = *((int ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++)
            buf[m++] = ubuf(array[i][mm]).d;
        }
      } if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[m++] = ubuf(vec[i]).d;
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata); 
          for (mm = 0; mm < cols; mm++)
            buf[m++] = ubuf(array[i][mm]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          for (mm = 0; mm < ncols; mm++)
            buf[m++] = ubuf(array[i][mm]).d;
        }
      }
    }
  }

  if (bonus_flag) m += pack_exchange_bonus(i,&buf[m]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVec::unpack_exchange(double *buf)
{
  int mm,nn,datatype,cols,collength,ncols;
  void *pdata,*plength;

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
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[nlocal] = buf[m++];
        } else if (cols > 0) {
          double **array = *((double ***) pdata);
          for (mm = 0; mm < cols; mm++)
            array[nlocal][mm] = buf[m++];
        } else {
          double **array = *((double ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
          else ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++)
            array[nlocal][mm] = buf[m++];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[nlocal] = ubuf(buf[m++]).i;
        } else if (cols > 0) {
          int **array = *((int ***) pdata);
          for (mm = 0; mm < cols; mm++)
            array[nlocal][mm] = (int) ubuf(buf[m++]).i;
        } else {
          int **array = *((int ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
          else ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++)
            array[nlocal][mm] = (int) ubuf(buf[m++]).i;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[nlocal] = (bigint) ubuf(buf[m++]).i;
        } else if (cols > 0) {
          bigint **array = *((bigint ***) pdata);
          for (mm = 0; mm < cols; mm++)
            array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
        } else {
          bigint **array = *((bigint ***) pdata);
          collength = mexchange.collength[nn];
          plength = mexchange.plength[nn];
          if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
          else ncols = (*((int **) plength))[nlocal];
          for (mm = 0; mm < ncols; mm++)
            array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
        }
      }
    }
  }

  if (bonus_flag) m += unpack_exchange_bonus(nlocal,&buf[m]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVec::size_restart()
{
  int i,nn,cols,collength,ncols;
  void *plength;

  // NOTE: need to worry about overflow of returned int N

  int nlocal = atom->nlocal;

  // 11 = length storage + id,type,mask,image,x,v

  int n = 11 * nlocal;

  if (nrestart) {
    for (nn = 0; nn < nrestart; nn++) {
      cols = mrestart.cols[nn];
      if (cols == 0) n += nlocal;
      else if (cols > 0) n += cols*nlocal;
      else {
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        for (i = 0; i < nlocal; i++) {
          if (collength) ncols = (*((int ***) plength))[i][collength-1];
          else ncols = (*((int **) plength))[i];
          n += ncols;
        }
      }
    }
  }

  if (bonus_flag) n += size_restart_bonus();

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVec::pack_restart(int i, double *buf)
{
  int mm,nn,datatype,cols,collength,ncols;
  void *pdata,*plength;

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
    if (datatype == DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        buf[m++] = vec[i];
      } else if (cols > 0) {
        double **array = *((double ***) pdata);
        for (mm = 0; mm < cols; mm++)
          buf[m++] = array[i][mm];
      } else {
        double **array = *((double ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[i][collength-1];
        else ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++)
          buf[m++] = array[i][mm];
      }
    } else if (datatype == INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        buf[m++] = ubuf(vec[i]).d;
      } else if (cols > 0) {
        int **array = *((int ***) pdata);
        for (mm = 0; mm < cols; mm++)
          buf[m++] = ubuf(array[i][mm]).d;
      } else {
        int **array = *((int ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[i][collength-1];
        else ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++)
          buf[m++] = ubuf(array[i][mm]).d;
      }
    } else if (datatype == BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        buf[m++] = ubuf(vec[i]).d;
      } else if (cols > 0) {
        bigint **array = *((bigint ***) pdata);
        for (mm = 0; mm < cols; mm++)
          buf[m++] = ubuf(array[i][mm]).d;
      } else {
        bigint **array = *((bigint ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[i][collength-1];
        else ncols = (*((int **) plength))[i];
        for (mm = 0; mm < ncols; mm++)
          buf[m++] = ubuf(array[i][mm]).d;
      }
    }
  }

  if (bonus_flag) m += pack_restart_bonus(i,&buf[m]);

  // if needed, restore values after packing

  pack_restart_post(i);

  // invoke fixes which store peratom restart info

  for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
    m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVec::unpack_restart(double *buf)
{
  int mm,nn,datatype,cols,collength,ncols;
  void *pdata,*plength;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
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
    if (datatype == DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = buf[m++];
      } else if (cols > 0) {
        double **array = *((double ***) pdata);
        for (mm = 0; mm < cols; mm++)
          array[nlocal][mm] = buf[m++];
      } else {
        double **array = *((double ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
        else ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++)
          array[nlocal][mm] = buf[m++];
      }
    } else if (datatype == INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        vec[nlocal] = ubuf(buf[m++]).i;
      } else if (cols > 0) {
        int **array = *((int ***) pdata);
        for (mm = 0; mm < cols; mm++)
            array[nlocal][mm] = (int) ubuf(buf[m++]).i;
      } else {
        int **array = *((int ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
        else ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++)
            array[nlocal][mm] = (int) ubuf(buf[m++]).i;
      }
    } else if (datatype == BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = (bigint) ubuf(buf[m++]).i;
      } else if (cols > 0) {
        bigint **array = *((bigint ***) pdata);
        for (mm = 0; mm < cols; mm++)
          array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
      } else {
        int **array = *((int ***) pdata);
        collength = mexchange.collength[nn];
        plength = mexchange.plength[nn];
        if (collength) ncols = (*((int ***) plength))[nlocal][collength-1];
        else ncols = (*((int **) plength))[nlocal];
        for (mm = 0; mm < ncols; mm++)
          array[nlocal][mm] = (bigint) ubuf(buf[m++]).i;
      }
    }
  }

  if (bonus_flag) m += unpack_restart_bonus(nlocal,&buf[m]);

  // if needed, initialize other peratom values

  unpack_restart_init(nlocal);

  // store extra restart info which fixes can unpack when instantiated

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
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
  int m,n,datatype,cols;
  void *pdata;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  // initialization additional fields

  for (n = 0; n < ncreate; n++) {
    pdata = mcreate.pdata[n];
    datatype = mcreate.datatype[n];
    cols = mcreate.cols[n];
    if (datatype == DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = 0.0;
      } else {
        double **array = *((double ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0.0;
      }
    } else if (datatype == INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        vec[nlocal] = 0;
      } else {
        int **array = *((int ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0;
      }
    } else if (datatype == BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = 0;
      } else {
        bigint **array = *((bigint ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0;
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

void AtomVec::data_atom(double *coord, imageint imagetmp, char **values)
{
  int m,n,datatype,cols;
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
    if (datatype == DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **) pdata);
        vec[nlocal] = utils::numeric(FLERR,values[ivalue++],true,lmp);
      } else {
        double **array = *((double ***) pdata);
        if (array == atom->x) {      // x was already set by coord arg
          ivalue += cols;
          continue;
        }
        for (m = 0; m < cols; m++)
          array[nlocal][m] = utils::numeric(FLERR,values[ivalue++],true,lmp);
      }
    } else if (datatype == INT) {
      if (cols == 0) {
        int *vec = *((int **) pdata);
        vec[nlocal] = utils::inumeric(FLERR,values[ivalue++],true,lmp);
      } else {
        int **array = *((int ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = utils::inumeric(FLERR,values[ivalue++],true,lmp);
      }
    } else if (datatype == BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **) pdata);
        vec[nlocal] = utils::bnumeric(FLERR,values[ivalue++],true,lmp);
      } else {
        bigint **array = *((bigint ***) pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = utils::bnumeric(FLERR,values[ivalue++],true,lmp);
      }
    }
  }

  // error checks applicable to all styles

  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  // if needed, modify unpacked values or initialize other peratom values

  data_atom_post(nlocal);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVec::pack_data(double **buf)
{
  int i,j,m,n,datatype,cols;
  void *pdata;

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    // if needed, change values before packing

    pack_data_pre(i);

    j = 0;
    for (n = 0; n < ndata_atom; n++) {
      pdata = mdata_atom.pdata[n];
      datatype = mdata_atom.datatype[n];
      cols = mdata_atom.cols[n];
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[i][j++] = vec[i];
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = array[i][m];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = ubuf(array[i][m]).d;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = ubuf(array[i][m]).d;
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
  int i,j,m,nn,datatype,cols;
  void *pdata;

  for (i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT,(tagint) ubuf(buf[i][0]).i);

    j = 1;
    for (nn = 1; nn < ndata_atom; nn++) {
      pdata = mdata_atom.pdata[nn];
      datatype = mdata_atom.datatype[nn];
      cols = mdata_atom.cols[nn];
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          fprintf(fp," %-1.16e",buf[i][j++]);
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            fprintf(fp," %-1.16e",buf[i][j++]);
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          fprintf(fp," %d",(int) ubuf(buf[i][j++]).i);
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            fprintf(fp," %d",(int) ubuf(buf[i][j++]).i);
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          fprintf(fp," " BIGINT_FORMAT,(bigint) ubuf(buf[i][j++]).i);
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            fprintf(fp," " BIGINT_FORMAT,(bigint) ubuf(buf[i][j++]).i);
        }
      }
    }

    fprintf(fp," %d %d %d\n",
            (int) ubuf(buf[i][j]).i,
            (int) ubuf(buf[i][j+1]).i,
            (int) ubuf(buf[i][j+2]).i);
  }
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVec::data_vel(int ilocal, char **values)
{
  int m,n,datatype,cols;
  void *pdata;

  double **v = atom->v;
  v[ilocal][0] = utils::numeric(FLERR,values[0],true,lmp);
  v[ilocal][1] = utils::numeric(FLERR,values[1],true,lmp);
  v[ilocal][2] = utils::numeric(FLERR,values[2],true,lmp);

  if (ndata_vel > 2) {
    int ivalue = 3;
    for (n = 2; n < ndata_vel; n++) {
      pdata = mdata_vel.pdata[n];
      datatype = mdata_vel.datatype[n];
      cols = mdata_vel.cols[n];
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          vec[ilocal] = utils::numeric(FLERR,values[ivalue++],true,lmp);
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::numeric(FLERR,values[ivalue++],true,lmp);
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          vec[ilocal] = utils::inumeric(FLERR,values[ivalue++],true,lmp);
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::inumeric(FLERR,values[ivalue++],true,lmp);
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          vec[ilocal] = utils::bnumeric(FLERR,values[ivalue++],true,lmp);
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            array[ilocal][m] = utils::bnumeric(FLERR,values[ivalue++],true,lmp);
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
  int i,j,m,n,datatype,cols;
  void *pdata;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    j = 0;
    for (n = 0; n < ndata_vel; n++) {
      pdata = mdata_vel.pdata[n];
      datatype = mdata_vel.datatype[n];
      cols = mdata_vel.cols[n];
      if (datatype == DOUBLE) {
        if (cols == 0) {
          double *vec = *((double **) pdata);
          buf[i][j++] = vec[i];
        } else {
          double **array = *((double ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = array[i][m];
        }
      } else if (datatype == INT) {
        if (cols == 0) {
          int *vec = *((int **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          int **array = *((int ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = ubuf(array[i][m]).d;
        }
      } else if (datatype == BIGINT) {
        if (cols == 0) {
          bigint *vec = *((bigint **) pdata);
          buf[i][j++] = ubuf(vec[i]).d;
        } else {
          bigint **array = *((bigint ***) pdata);
          for (m = 0; m < cols; m++)
            buf[i][j++] = ubuf(array[i][m]).d;
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
  int i,j,m,nn,datatype,cols;
  void *pdata;

  for (i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3]);

    if (ndata_vel) {
      j = 4;
      for (nn = 0; nn < ndata_vel; nn++) {
        pdata = mdata_vel.pdata[nn];
        datatype = mdata_vel.datatype[nn];
        cols = mdata_vel.cols[nn];
        if (datatype == DOUBLE) {
          if (cols == 0) {
            double *vec = *((double **) pdata);
            fprintf(fp," %-1.16e",buf[i][j++]);
          } else {
            double **array = *((double ***) pdata);
            for (m = 0; m < cols; m++)
              fprintf(fp," %-1.16e",buf[i][j++]);
          }
        } else if (datatype == INT) {
          if (cols == 0) {
            int *vec = *((int **) pdata);
            fprintf(fp," %d",(int) ubuf(buf[i][j++]).i);
          } else {
            int **array = *((int ***) pdata);
            for (m = 0; m < cols; m++)
              fprintf(fp," %d",(int) ubuf(buf[i][j++]).i);
          }
        } else if (datatype == BIGINT) {
          if (cols == 0) {
            bigint *vec = *((bigint **) pdata);
            fprintf(fp," " BIGINT_FORMAT,(bigint) ubuf(buf[i][j++]).i);
          } else {
            bigint **array = *((bigint ***) pdata);
            for (m = 0; m < cols; m++)
              fprintf(fp," " BIGINT_FORMAT,(bigint) ubuf(buf[i][j++]).i);
          }
        }
      }
    }

    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   pack bond info for data file into buf if non-NULL
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

  int i,j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_bond[i]; j++) {
        if (bond_type[i][j] == 0) continue;
        if (buf) {
          buf[m][0] = MAX(bond_type[i][j],-bond_type[i][j]);
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
            buf[m][0] = MAX(bond_type[i][j],-bond_type[i][j]);
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
  for (int i = 0; i < n; i++) {
    fprintf(fp,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index,buf[i][0],buf[i][1],buf[i][2]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   pack angle info for data file into buf if non-NULL
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

  int i,j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_angle[i]; j++) {
        if (angle_type[i][j] == 0) continue;
        if (buf) {
          buf[m][0] = MAX(angle_type[i][j],-angle_type[i][j]);
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
            buf[m][0] = MAX(angle_type[i][j],-angle_type[i][j]);
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
  for (int i = 0; i < n; i++) {
    fprintf(fp,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " "
            TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index,buf[i][0],buf[i][1],buf[i][2],buf[i][3]);
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

  int i,j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++) {
        if (buf) {
          buf[m][0] = MAX(dihedral_type[i][j],-dihedral_type[i][j]);
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
            buf[m][0] = MAX(dihedral_type[i][j],-dihedral_type[i][j]);
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
  for (int i = 0; i < n; i++) {
    fprintf(fp,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " "
            TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index,buf[i][0],buf[i][1],buf[i][2],buf[i][3],buf[i][4]);
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

  int i,j;
  int m = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_improper[i]; j++) {
        if (buf) {
          buf[m][0] = MAX(improper_type[i][j],-improper_type[i][j]);
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
            buf[m][0] = MAX(improper_type[i][j],-improper_type[i][j]);
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
  for (int i = 0; i < n; i++) {
    fprintf(fp,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " "
            TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index,buf[i][0],buf[i][1],buf[i][2],buf[i][3],buf[i][4]);
    index++;
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVec::memory_usage()
{
  int datatype,cols,index,maxcols;
  void *pdata;

  bigint bytes = 0;

  bytes += memory->usage(tag,nmax);
  bytes += memory->usage(type,nmax);
  bytes += memory->usage(mask,nmax);
  bytes += memory->usage(image,nmax);
  bytes += memory->usage(x,nmax,3);
  bytes += memory->usage(v,nmax,3);
  bytes += memory->usage(f,nmax*nthreads,3);

  for (int i = 0; i < ngrow; i++) {
    pdata = mgrow.pdata[i];
    datatype = mgrow.datatype[i];
    cols = mgrow.cols[i];
    index = mgrow.index[i];
    if (datatype == DOUBLE) {
      if (cols == 0) {
        bytes += memory->usage(*((double **) pdata),nmax*threads[i]);
      } else if (cols > 0) {
        bytes += memory->usage(*((double ***) pdata),nmax*threads[i],cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((double ***) pdata),nmax*threads[i],maxcols);
      }
    } else if (datatype == INT) {
      if (cols == 0) {
        bytes += memory->usage(*((int **) pdata),nmax*threads[i]);
      } else if (cols > 0) {
        bytes += memory->usage(*((int ***) pdata),nmax*threads[i],cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((int ***) pdata),nmax*threads[i],maxcols);
      }
    } else if (datatype == BIGINT) {
      if (cols == 0) {
        bytes += memory->usage(*((bigint **) pdata),nmax*threads[i]);
      } else if (cols > 0) {
        bytes += memory->usage(*((bigint ***) pdata),nmax*threads[i],cols);
      } else {
        maxcols = *(mgrow.maxcols[i]);
        bytes += memory->usage(*((bigint ***) pdata),nmax*threads[i],maxcols);
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
  int n,cols;

  if (strstr(fields_data_atom,"id ") != fields_data_atom)
    error->all(FLERR,"Atom style fields_data_atom must have id as first field");
  if (strstr(fields_data_vel,"id v") != fields_data_vel)
    error->all(FLERR,"Atom style fields_data_vel must have "
               "'id v' as first fields");

  // process field strings
  // return # of fields and matching index into atom->peratom (in Method struct)

  ngrow = process_fields(fields_grow,default_grow,&mgrow);
  ncopy = process_fields(fields_copy,default_copy,&mcopy);
  ncomm = process_fields(fields_comm,default_comm,&mcomm);
  ncomm_vel = process_fields(fields_comm_vel,default_comm_vel,&mcomm_vel);
  nreverse = process_fields(fields_reverse,default_reverse,&mreverse);
  nborder = process_fields(fields_border,default_border,&mborder);
  nborder_vel = process_fields(fields_border_vel,default_border_vel,&mborder_vel);
  nexchange = process_fields(fields_exchange,default_exchange,&mexchange);
  nrestart = process_fields(fields_restart,default_restart,&mrestart);
  ncreate = process_fields(fields_create,default_create,&mcreate);
  ndata_atom = process_fields(fields_data_atom,default_data_atom,&mdata_atom);
  ndata_vel = process_fields(fields_data_vel,default_data_vel,&mdata_vel);

  // populate field-based data struct for each method to use

  create_method(ngrow,&mgrow);
  create_method(ncopy,&mcopy);
  create_method(ncomm,&mcomm);
  create_method(ncomm_vel,&mcomm_vel);
  create_method(nreverse,&mreverse);
  create_method(nborder,&mborder);
  create_method(nborder_vel,&mborder_vel);
  create_method(nexchange,&mexchange);
  create_method(nrestart,&mrestart);
  create_method(ncreate,&mcreate);
  create_method(ndata_atom,&mdata_atom);
  create_method(ndata_vel,&mdata_vel);

  // create threads data struct for grow and memory_usage to use

  threads = new int[ngrow];
  for (int i = 0; i < ngrow; i++) {
    Atom::PerAtom *field = &atom->peratom[mgrow.index[i]];
    if (field->threadflag) threads[i] = nthreads;
    else threads[i] = 1;
  }

  // set style-specific sizes

  comm_x_only = 1;
  if (ncomm) comm_x_only = 0;
  if (bonus_flag && size_forward_bonus) comm_x_only = 0;

  if (nreverse == 0) comm_f_only = 1;
  else comm_f_only = 0;

  size_forward = 3;
  for (n = 0; n < ncomm; n++) {
    cols = mcomm.cols[n];
    if (cols == 0) size_forward++;
    else size_forward += cols;
  }
  if (bonus_flag) size_forward += size_forward_bonus;

  size_reverse = 3;
  for (n = 0; n < nreverse; n++) {
    cols = mreverse.cols[n];
    if (cols == 0) size_reverse++;
    else size_reverse += cols;
  }

  size_border = 6;
  for (n = 0; n < nborder; n++) {
    cols = mborder.cols[n];
    if (cols == 0) size_border++;
    else size_border += cols;
  }
  if (bonus_flag) size_border += size_border_bonus;

  size_velocity = 3;
  for (n = 0; n < ncomm_vel; n++) {
    cols = mcomm_vel.cols[n];
    if (cols == 0) size_velocity++;
    else size_velocity += cols;
  }

  size_data_atom = 0;
  for (n = 0; n < ndata_atom; n++) {
    cols = mdata_atom.cols[n];
    if (strcmp(atom->peratom[mdata_atom.index[n]].name,"x") == 0) 
      xcol_data = size_data_atom + 1;
    if (cols == 0) size_data_atom++;
    else size_data_atom += cols;
  }

  size_data_vel = 0;
  for (n = 0; n < ndata_vel; n++) {
    cols = mdata_vel.cols[n];
    if (cols == 0) size_data_vel++;
    else size_data_vel += cols;
  }
}

/* ----------------------------------------------------------------------
   process a single field string
------------------------------------------------------------------------- */

int AtomVec::process_fields(char *str, const char *default_str, Method *method)
{
  if (str == NULL) {
    method->index = NULL;
    return 0;
  }

  // tokenize words in both strings

  char *copy1,*copy2;
  char **words,**defwords;
  int nfield = tokenize(str,words,copy1);
  int ndef = tokenize((char *) default_str,defwords,copy2);

  // process fields one by one, add to index vector
  
  Atom::PerAtom *peratom = atom->peratom;
  int nperatom = atom->nperatom;

  int *index = new int[nfield];
  int match;

  for (int i = 0; i < nfield; i++) {
    
    // find field in master Atom::peratom list

    for (match = 0; match < nperatom; match++)
      if (strcmp(words[i],peratom[match].name) == 0) break;
    if (match == nperatom) {
      char str[128];
      sprintf(str,"Peratom field %s not recognized",words[i]);
      error->all(FLERR,str);
    }
    index[i] = match;

    // error if field appears multiple times

    for (match = 0; match < i; match++)
      if (index[i] == index[match]) {
	char str[128];
	sprintf(str,"Peratom field %s is repeated",words[i]);
	error->all(FLERR,str);
      }

    // error if field is in default str

    for (match = 0; match < ndef; match++)
      if (strcmp(words[i],defwords[match]) == 0) {
	char str[128];
	sprintf(str,"Peratom field %s is a default",words[i]);
	error->all(FLERR,str);
      }

  }

  delete [] copy1;
  delete [] copy2;
  delete [] words;
  delete [] defwords;

  method->index = index;
  return nfield;
}

/* ----------------------------------------------------------------------
   tokenize str into white-space separated words
   return nwords = number of words
   return words = vector of ptrs to each word
   also return copystr since words points into it, caller will delete copystr
------------------------------------------------------------------------- */

int AtomVec::tokenize(char *str, char **&words, char *&copystr)
{
  int n = strlen(str) + 1;
  copystr = new char[n];
  strcpy(copystr,str);

  int nword = atom->count_words(copystr);
  words = new char*[nword];

  nword = 0;
  char *word = strtok(copystr," ");
  while (word) {
    words[nword++] = word;
    word = strtok(NULL," ");
  }

  return nword;
}

/* ----------------------------------------------------------------------
   create a method data structs for processing fields
------------------------------------------------------------------------- */

void AtomVec::create_method(int nfield, Method *method)
{
  method->pdata = new void*[nfield];
  method->datatype = new int[nfield];
  method->cols = new int[nfield];
  method->maxcols = new int*[nfield];
  method->collength = new int[nfield];
  method->plength = new void*[nfield];

  for (int i = 0; i < nfield; i++) {
    Atom::PerAtom *field = &atom->peratom[method->index[i]];
    method->pdata[i] = (void *) field->address;
    method->datatype[i] = field->datatype;
    method->cols[i] = field->cols;
    if (method->cols[i] < 0) {
      method->maxcols[i] = field->address_maxcols;
      method->collength[i] = field->collength;
      method->plength[i] = field->address_length;
    }
  }
}

/* ----------------------------------------------------------------------
   free memory in a method data structs
------------------------------------------------------------------------- */

void AtomVec::init_method(Method *method)
{
  method->pdata = NULL;
  method->datatype = NULL;
  method->cols = NULL;
  method->maxcols = NULL;
  method->collength = NULL;
  method->plength = NULL;
  method->index = NULL;
}

/* ----------------------------------------------------------------------
   free memory in a method data structs
------------------------------------------------------------------------- */

void AtomVec::destroy_method(Method *method)
{
  delete [] method->pdata;
  delete [] method->datatype;
  delete [] method->cols;
  delete [] method->maxcols;
  delete [] method->collength;
  delete [] method->plength;
  delete [] method->index;
}
