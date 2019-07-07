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

#include <cstdlib>
#include <cstring>
#include "atom_vec_hybrid_kokkos.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::AtomVecHybridKokkos(LAMMPS *lmp) : AtomVecKokkos(lmp) {
  no_comm_vel_flag = 1;
  no_border_vel_flag = 1;
}

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::~AtomVecHybridKokkos()
{
  for (int k = 0; k < nstyles; k++) delete styles[k];
  delete [] styles;
  for (int k = 0; k < nstyles; k++) delete [] keywords[k];
  delete [] keywords;
}

/* ----------------------------------------------------------------------
   process sub-style args
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::process_args(int narg, char **arg)
{
  // build list of all known atom styles

  build_styles();

  // allocate list of sub-styles as big as possibly needed if no extra args

  styles = new AtomVec*[narg];
  keywords = new char*[narg];

  // allocate each sub-style
  // call process_args() with set of args that are not atom style names
  // use known_style() to determine which args these are

  int i,jarg,dummy;

  int iarg = 0;
  nstyles = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"hybrid") == 0)
      error->all(FLERR,"Atom style hybrid cannot have hybrid as an argument");
    for (i = 0; i < nstyles; i++)
      if (strcmp(arg[iarg],keywords[i]) == 0)
        error->all(FLERR,"Atom style hybrid cannot use same atom style twice");
    styles[nstyles] = atom->new_avec(arg[iarg],1,dummy);
    keywords[nstyles] = new char[strlen(arg[iarg])+1];
    strcpy(keywords[nstyles],arg[iarg]);
    jarg = iarg + 1;
    while (jarg < narg && !known_style(arg[jarg])) jarg++;
    styles[nstyles]->process_args(jarg-iarg-1,&arg[iarg+1]);
    iarg = jarg;
    nstyles++;
  }

  // free allstyles created by build_styles()

  for (int i = 0; i < nallstyles; i++) delete [] allstyles[i];
  delete [] allstyles;

  // hybrid settings are MAX or MIN of sub-style settings
  // hybrid sizes are minimal values plus extra values for each sub-style

  molecular = 0;
  comm_x_only = comm_f_only = 1;

  size_forward = 3;
  size_reverse = 3;
  size_border = 6;
  size_data_atom = 5;
  size_data_vel = 4;
  xcol_data = 3;

  for (int k = 0; k < nstyles; k++) {
    if ((styles[k]->molecular == 1 && molecular == 2) ||
        (styles[k]->molecular == 2 && molecular == 1))
      error->all(FLERR,"Cannot mix molecular and molecule template "
                 "atom styles");
    molecular = MAX(molecular,styles[k]->molecular);

    bonds_allow = MAX(bonds_allow,styles[k]->bonds_allow);
    angles_allow = MAX(angles_allow,styles[k]->angles_allow);
    dihedrals_allow = MAX(dihedrals_allow,styles[k]->dihedrals_allow);
    impropers_allow = MAX(impropers_allow,styles[k]->impropers_allow);
    mass_type = MAX(mass_type,styles[k]->mass_type);
    dipole_type = MAX(dipole_type,styles[k]->dipole_type);
    forceclearflag = MAX(forceclearflag,styles[k]->forceclearflag);

    if (styles[k]->molecular == 2) onemols = styles[k]->onemols;

    comm_x_only = MIN(comm_x_only,styles[k]->comm_x_only);
    comm_f_only = MIN(comm_f_only,styles[k]->comm_f_only);
    size_forward += styles[k]->size_forward - 3;
    size_reverse += styles[k]->size_reverse - 3;
    size_border += styles[k]->size_border - 6;
    size_data_atom += styles[k]->size_data_atom - 5;
    size_data_vel += styles[k]->size_data_vel - 4;
  }

  size_velocity = 3;
  if (atom->omega_flag) size_velocity += 3;
  if (atom->angmom_flag) size_velocity += 3;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::init()
{
  AtomVec::init();
  for (int k = 0; k < nstyles; k++) styles[k]->init();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  // sub-styles perform all reallocation
  // turn off nextra_grow so hybrid can do that once below

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;
  for (int k = 0; k < nstyles; k++) styles[k]->grow(nmax);
  atom->nextra_grow = tmp;

  // insure hybrid local ptrs and sub-style ptrs are up to date
  // for sub-styles, do this in case
  //   multiple sub-style reallocs of same array occurred

  grow_reset();

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow_reset()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.d_view;
  h_tag = atomKK->k_tag.h_view;

  type = atomKK->type;
  d_type = atomKK->k_type.d_view;
  h_type = atomKK->k_type.h_view;

  mask = atomKK->mask;
  d_mask = atomKK->k_mask.d_view;
  h_mask = atomKK->k_mask.h_view;

  image = atomKK->image;
  d_image = atomKK->k_image.d_view;
  h_image = atomKK->k_image.h_view;

  x = atomKK->x;
  d_x = atomKK->k_x.d_view;
  h_x = atomKK->k_x.h_view;

  v = atomKK->v;
  d_v = atomKK->k_v.d_view;
  h_v = atomKK->k_v.h_view;

  f = atomKK->f;
  d_f = atomKK->k_f.d_view;
  h_f = atomKK->k_f.h_view;

  v = atomKK->v;
  d_v = atomKK->k_v.d_view;
  h_v = atomKK->k_v.h_view;

  omega = atomKK->omega;
  d_omega = atomKK->k_omega.d_view;
  h_omega = atomKK->k_omega.h_view;

  angmom = atomKK->angmom;
  d_angmom = atomKK->k_angmom.d_view;
  h_angmom = atomKK->k_angmom.h_view;

  for (int k = 0; k < nstyles; k++) styles[k]->grow_reset();
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J for all sub-styles
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::copy(int i, int j, int delflag)
{
  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;
  for (int k = 0; k < nstyles; k++) styles[k]->copy(i,j,delflag);
  atom->nextra_grow = tmp;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::clear_bonus()
{
  for (int k = 0; k < nstyles; k++) styles[k]->clear_bonus();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::force_clear(int n, size_t nbytes)
{
  for (int k = 0; k < nstyles; k++)
    if (styles[k]->forceclearflag) styles[k]->force_clear(n,nbytes);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &k_sendlist,
                     const int & iswap,
                     const DAT::tdual_xfloat_2d &buf,
                     const int &pbc_flag, const int pbc[])
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}
void AtomVecHybridKokkos::unpack_comm_kokkos(const int &n, const int &nfirst,
                        const DAT::tdual_xfloat_2d &buf)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
}
int AtomVecHybridKokkos::pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                   const int & iswap, const int nfirst,
                   const int &pbc_flag, const int pbc[])
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}
int AtomVecHybridKokkos::pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                       DAT::tdual_xfloat_2d buf,int iswap,
                       int pbc_flag, int *pbc, ExecutionSpace space)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}
void AtomVecHybridKokkos::unpack_border_kokkos(const int &n, const int &nfirst,
                          const DAT::tdual_xfloat_2d &buf,
                          ExecutionSpace space)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
}
int AtomVecHybridKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space, int dim,
                         X_FLOAT lo, X_FLOAT hi)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}
int AtomVecHybridKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                           int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                           ExecutionSpace space)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  atomKK->sync(Host,X_MASK);

  int i,j,k,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
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
      buf[m++] = h_x(j,0) + dx;
      buf[m++] = h_x(j,1) + dy;
      buf[m++] = h_x(j,2) + dz;
    }
  }

  // pack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_comm_hybrid(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  atomKK->sync(Host,X_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);

  int i,j,k,m;
  double dx,dy,dz,dvx,dvy,dvz;
  int omega_flag = atom->omega_flag;
  int angmom_flag = atom->angmom_flag;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = h_v(j,0);
      buf[m++] = h_v(j,1);
      buf[m++] = h_v(j,2);
      if (omega_flag) {
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
      if (angmom_flag) {
        buf[m++] = h_angmom(j,0);
        buf[m++] = h_angmom(j,1);
        buf[m++] = h_angmom(j,2);
      }
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
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        if (omega_flag) {
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
        if (angmom_flag) {
          buf[m++] = h_angmom(j,0);
          buf[m++] = h_angmom(j,1);
          buf[m++] = h_angmom(j,2);
        }
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        if (h_mask[i] & deform_groupbit) {
          buf[m++] = h_v(j,0) + dvx;
          buf[m++] = h_v(j,1) + dvy;
          buf[m++] = h_v(j,2) + dvz;
        } else {
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
        }
        if (omega_flag) {
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
        if (angmom_flag) {
          buf[m++] = h_angmom(j,0);
          buf[m++] = h_angmom(j,1);
          buf[m++] = h_angmom(j,2);
        }
      }
    }
  }

  // pack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_comm_hybrid(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
  }

  atomKK->modified(Host,X_MASK);

  // unpack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->unpack_comm_hybrid(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_comm_vel(int n, int first, double *buf)
{
  int i,k,m,last;
  int omega_flag = atom->omega_flag;
  int angmom_flag = atom->angmom_flag;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
    if (omega_flag) {
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    if (angmom_flag) {
      h_angmom(i,0) = buf[m++];
      h_angmom(i,1) = buf[m++];
      h_angmom(i,2) = buf[m++];
    }
  }

  atomKK->modified(Host,X_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);

  // unpack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->unpack_comm_hybrid(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_reverse(int n, int first, double *buf)
{
  atomKK->sync(Host,F_MASK);

  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = h_f(i,0);
    buf[m++] = h_f(i,1);
    buf[m++] = h_f(i,2);
  }

  // pack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_reverse_hybrid(n,first,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_f(j,0) += buf[m++];
    h_f(j,1) += buf[m++];
    h_f(j,2) += buf[m++];
  }

  atomKK->modified(Host,F_MASK);

  // unpack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->unpack_reverse_hybrid(n,list,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  atomKK->sync(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK);

  int i,j,k,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
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
      buf[m++] = h_x(j,0) + dx;
      buf[m++] = h_x(j,1) + dy;
      buf[m++] = h_x(j,2) + dz;
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
    }
  }

  // pack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_border_hybrid(n,list,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  atomKK->sync(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);
  int i,j,k,m;
  double dx,dy,dz,dvx,dvy,dvz;
  int omega_flag = atom->omega_flag;
  int angmom_flag = atom->angmom_flag;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
      buf[m++] = h_v(j,0);
      buf[m++] = h_v(j,1);
      buf[m++] = h_v(j,2);
      if (omega_flag) {
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
      if (angmom_flag) {
        buf[m++] = h_angmom(j,0);
        buf[m++] = h_angmom(j,1);
        buf[m++] = h_angmom(j,2);
      }
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
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag[j]).d;
        buf[m++] = ubuf(h_type[j]).d;
        buf[m++] = ubuf(h_mask[j]).d;
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        if (omega_flag) {
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
        if (angmom_flag) {
          buf[m++] = h_angmom(j,0);
          buf[m++] = h_angmom(j,1);
          buf[m++] = h_angmom(j,2);
        }
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag[j]).d;
        buf[m++] = ubuf(h_type[j]).d;
        buf[m++] = ubuf(h_mask[j]).d;
        if (h_mask[i] & deform_groupbit) {
          buf[m++] = h_v(j,0) + dvx;
          buf[m++] = h_v(j,1) + dvy;
          buf[m++] = h_v(j,2) + dvz;
        } else {
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
        }
        if (omega_flag) {
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
        if (angmom_flag) {
          buf[m++] = h_angmom(j,0);
          buf[m++] = h_angmom(j,1);
          buf[m++] = h_angmom(j,2);
        }
      }
    }
  }

  // pack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_border_hybrid(n,list,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_border(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag[i] = (tagint) ubuf(buf[m++]).i;
    h_type[i] = (int) ubuf(buf[m++]).i;
    h_mask[i] = (int) ubuf(buf[m++]).i;
  }

  atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK);

  // unpack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->unpack_border_hybrid(n,first,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_border_vel(int n, int first, double *buf)
{
  int i,k,m,last;
  int omega_flag = atom->omega_flag;
  int angmom_flag = atom->angmom_flag;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag[i] = (tagint) ubuf(buf[m++]).i;
    h_type[i] = (int) ubuf(buf[m++]).i;
    h_mask[i] = (int) ubuf(buf[m++]).i;
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
    if (omega_flag) {
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    if (angmom_flag) {
      h_angmom(i,0) = buf[m++];
      h_angmom(i,1) = buf[m++];
      h_angmom(i,2) = buf[m++];
    }
  }

  atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);

  // unpack sub-style contributions as contiguous chunks

  for (k = 0; k < nstyles; k++)
    m += styles[k]->unpack_border_hybrid(n,first,&buf[m]);

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   pack each sub-style one after the other
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_exchange(int i, double *buf)
{
  int k,m;

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;

  m = 0;
  for (k = 0; k < nstyles; k++)
    m += styles[k]->pack_exchange(i,&buf[m]);

  atom->nextra_grow = tmp;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for single atom received from another proc
   unpack each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::unpack_exchange(double *buf)
{
  int k,m;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int tmp = atom->nextra_grow;
  atom->nextra_grow = 0;

  m = 0;
  for (k = 0; k < nstyles; k++) {
    m += styles[k]->unpack_exchange(&buf[m]);
    atom->nlocal--;
  }

  atom->nextra_grow = tmp;

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

int AtomVecHybridKokkos::size_restart()
{
  int tmp = atom->nextra_restart;
  atom->nextra_restart = 0;

  int n = 0;
  for (int k = 0; k < nstyles; k++)
    n += styles[k]->size_restart();

  atom->nextra_restart = tmp;

  int nlocal = atom->nlocal;
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (int i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   pack each sub-style one after the other
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_restart(int i, double *buf)
{
  int tmp = atom->nextra_restart;
  atom->nextra_restart = 0;

  int m = 0;
  for (int k = 0; k < nstyles; k++)
    m += styles[k]->pack_restart(i,&buf[m]);

  atom->nextra_restart = tmp;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
   unpack each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int tmp = atom->nextra_store;
  atom->nextra_store = 0;

  int m = 0;
  for (int k = 0; k < nstyles; k++) {
    m += styles[k]->unpack_restart(&buf[m]);
    atom->nlocal--;
  }
  atom->nextra_store = tmp;

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
   create each sub-style one after the other
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  for (int k = 0; k < nstyles; k++) {
    styles[k]->create_atom(itype,coord);
    atom->nlocal--;
  }
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   grow() occurs here so arrays for all sub-styles are grown
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::data_atom(double *coord, imageint imagetmp, char **values)
{
  atomKK->sync(Host,X_MASK|TAG_MASK|TYPE_MASK|IMAGE_MASK|MASK_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  h_tag[nlocal] = ATOTAGINT(values[0]);
  h_type[nlocal] = force->inumeric(FLERR,values[1]);
  if (h_type[nlocal] <= 0 || h_type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom h_type in Atoms section of data file");

  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];

  h_image[nlocal] = imagetmp;
  h_mask[nlocal] = 1;

  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;
  if (atom->omega_flag) {
    h_omega(nlocal,0) = 0.0;
    h_omega(nlocal,1) = 0.0;
    h_omega(nlocal,2) = 0.0;
  }
  if (atom->angmom_flag) {
    h_angmom(nlocal,0) = 0.0;
    h_angmom(nlocal,1) = 0.0;
    h_angmom(nlocal,2) = 0.0;
  }

  atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|IMAGE_MASK|MASK_MASK|V_MASK|OMEGA_MASK/*|ANGMOM_MASK*/);

  // each sub-style parses sub-style specific values

  int m = 5;
  for (int k = 0; k < nstyles; k++)
    m += styles[k]->data_atom_hybrid(nlocal,&values[m]);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::data_vel(int m, char **values)
{
  atomKK->sync(Host,V_MASK);

  h_v(m,0) = force->numeric(FLERR,values[0]);
  h_v(m,1) = force->numeric(FLERR,values[1]);
  h_v(m,2) = force->numeric(FLERR,values[2]);

  atomKK->modified(Host,V_MASK);

  // each sub-style parses sub-style specific values

  int n = 3;
  for (int k = 0; k < nstyles; k++)
    n += styles[k]->data_vel_hybrid(m,&values[n]);
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_data(double **buf)
{
  atomKK->sync(Host,TAG_MASK|TYPE_MASK|X_MASK);

  int k,m;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = ubuf(h_type[i]).d;
    buf[i][2] = h_x(i,0);
    buf[i][3] = h_x(i,1);
    buf[i][4] = h_x(i,2);

    m = 5;
    for (k = 0; k < nstyles; k++)
      m += styles[k]->pack_data_hybrid(i,&buf[i][m]);

    buf[i][m] = ubuf((h_image[i] & IMGMASK) - IMGMAX).d;
    buf[i][m+1] = ubuf((h_image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][m+2] = ubuf((h_image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 h_image flags
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::write_data(FILE *fp, int n, double **buf)
{
  int k,m;

  for (int i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4]);

    m = 5;
    for (k = 0; k < nstyles; k++)
      m += styles[k]->write_data_hybrid(fp,&buf[i][m]);

    fprintf(fp," %d %d %d\n",
            (int) ubuf(buf[i][m]).i,(int) ubuf(buf[i][m+1]).i,
            (int) ubuf(buf[i][m+2]).i);
  }
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_vel(double **buf)
{
  atomKK->sync(Host,V_MASK);

  int k,m;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = h_v(i,0);
    buf[i][2] = h_v(i,1);
    buf[i][3] = h_v(i,2);

    m = 4;
    for (k = 0; k < nstyles; k++)
      m += styles[k]->pack_vel_hybrid(i,&buf[i][m]);
  }
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::write_vel(FILE *fp, int n, double **buf)
{
  int k,m;

  for (int i = 0; i < n; i++) {
    fprintf(fp,TAGINT_FORMAT " %g %g %g",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3]);

    m = 4;
    for (k = 0; k < nstyles; k++)
      m += styles[k]->write_vel_hybrid(fp,&buf[i][m]);

    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   returned value encodes which sub-style and index returned by sub-style
   return -1 if name is unknown to any sub-styles
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::property_atom(char *name)
{
  for (int k = 0; k < nstyles; k++) {
    int index = styles[k]->property_atom(name);
    if (index >= 0) return index*nstyles + k;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_property_atom(int multiindex, double *buf,
                                       int nvalues, int groupbit)
{
  int k = multiindex % nstyles;
  int index = multiindex/nstyles;
  styles[k]->pack_property_atom(index,buf,nvalues,groupbit);
}

/* ----------------------------------------------------------------------
   allstyles = list of all atom styles in this LAMMPS executable
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::build_styles()
{
  nallstyles = 0;
#define ATOM_CLASS
#define AtomStyle(key,Class) nallstyles++;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS

  allstyles = new char*[nallstyles];

  int n;
  nallstyles = 0;
#define ATOM_CLASS
#define AtomStyle(key,Class)                \
  n = strlen(#key) + 1;                     \
  allstyles[nallstyles] = new char[n];      \
  strcpy(allstyles[nallstyles],#key);       \
  nallstyles++;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
}

/* ----------------------------------------------------------------------
   allstyles = list of all known atom styles
------------------------------------------------------------------------- */

int AtomVecHybridKokkos::known_style(char *str)
{
  for (int i = 0; i < nallstyles; i++)
    if (strcmp(str,allstyles[i]) == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecHybridKokkos::memory_usage()
{
  bigint bytes = 0;
  for (int k = 0; k < nstyles; k++) bytes += styles[k]->memory_usage();
  return bytes;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) ((AtomVecKokkos*) styles[k])->sync(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) ((AtomVecKokkos*) styles[k])->sync_overlapping_device(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::modified(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) ((AtomVecKokkos*) styles[k])->modified(space,h_mask);
}
