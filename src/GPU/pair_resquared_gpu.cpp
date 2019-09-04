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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "pair_resquared_gpu.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"
#include "domain.h"
#include "update.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int re_gpu_init(const int ntypes, double **shape, double **well,
                double **cutsq, double **sigma, double **epsilon,
                int **form, double **host_lj1,
                double **host_lj2, double **host_lj3, double **host_lj4,
                double **offset, double *special_lj, const int nlocal,
                const int nall,        const int max_nbors, const int maxspecial,
                const double cell_size,        int &gpu_mode, FILE *screen);
void re_gpu_clear();
int ** re_gpu_compute_n(const int ago, const int inum, const int nall,
                        double **host_x, int *host_type, double *sublo,
                        double *subhi, tagint *tag, int **nspecial, tagint **special,
                        const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        int **ilist, int **jnum, const double cpu_time,
                        bool &success, double **host_quat);
int * re_gpu_compute(const int ago, const int inum, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double **host_quat);
double re_gpu_bytes();

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

/* ---------------------------------------------------------------------- */

PairRESquaredGPU::PairRESquaredGPU(LAMMPS *lmp) : PairRESquared(lmp),
                                                gpu_mode(GPU_FORCE)
{
  reinitflag = 0;
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Pair resquared/gpu requires atom style ellipsoid");
  quat_nmax = 0;
  quat = NULL;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairRESquaredGPU::~PairRESquaredGPU()
{
  re_gpu_clear();
  cpu_time = 0.0;
  memory->destroy(quat);
}

/* ---------------------------------------------------------------------- */

void PairRESquaredGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;

  if (nall > quat_nmax) {
    quat_nmax = static_cast<int>(1.1 * nall);
    memory->grow(quat, quat_nmax, 4, "pair:quat");
  }
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  for (int i=0; i<nall; i++) {
    int qi = ellipsoid[i];
    if (qi > -1) {
      quat[i][0] = bonus[qi].quat[0];
      quat[i][1] = bonus[qi].quat[1];
      quat[i][2] = bonus[qi].quat[2];
      quat[i][3] = bonus[qi].quat[3];
    }
  }

  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = re_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                  atom->type, domain->sublo, domain->subhi,
                                  atom->tag, atom->nspecial, atom->special,
                                  eflag, vflag, eflag_atom, vflag_atom,
                                  host_start, &ilist, &numneigh, cpu_time,
                                  success, quat);
  } else {
    inum = list->inum;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    ilist = re_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                           list->ilist, numneigh, firstneigh, eflag, vflag,
                           eflag_atom, vflag_atom, host_start,
                           cpu_time, success, quat);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  if (host_start < inum) {
    cpu_time = MPI_Wtime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairRESquaredGPU::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with resquared/gpu pair style");
  if (!atom->ellipsoid_flag)
    error->all(FLERR,"Pair resquared/gpu requires atom style ellipsoid");

  // per-type shape precalculations
  // require that atom shapes are identical within each type
  // if shape = 0 for point particle, set shape = 1 as required by Gay-Berne

  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->shape_consistency(i,shape1[i][0],shape1[i][1],shape1[i][2]))
      error->all(FLERR,"Pair resquared/gpu requires atoms with same type have same shape");
    if (setwell[i]) {
      shape2[i][0] = shape1[i][0]*shape1[i][0];
      shape2[i][1] = shape1[i][1]*shape1[i][1];
      shape2[i][2] = shape1[i][2]*shape1[i][2];
      lshape[i] = shape1[i][0]*shape1[i][1]*shape1[i][2];
    }
  }

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }

  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = re_gpu_init(atom->ntypes+1, shape1, well, cutsq, sigma,
                            epsilon, form, lj1, lj2, lj3, lj4, offset,
                            force->special_lj, atom->nlocal,
                            atom->nlocal+atom->nghost, 300, maxspecial,
                            cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
  quat_nmax = static_cast<int>(1.1 * (atom->nlocal + atom->nghost));
  memory->grow(quat, quat_nmax, 4, "pair:quat");
}

/* ---------------------------------------------------------------------- */

double PairRESquaredGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + memory->usage(quat,quat_nmax)+re_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairRESquaredGPU::cpu_compute(int start, int inum, int eflag,
                                   int /* vflag */, int *ilist,
                                   int *numneigh, int **firstneigh)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  int *jlist;
  RE2Vars wi,wj;

  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // not a LJ sphere

    if (lshape[itype] != 0.0) precompute_i(i,wi);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {
        switch (form[itype][jtype]) {

         case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (eflag) one_eng =
              r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
              offset[itype][jtype];
          fforce[0] = r12[0]*forcelj;
          fforce[1] = r12[1]*forcelj;
          fforce[2] = r12[2]*forcelj;
          break;

         case SPHERE_ELLIPSE:
          precompute_i(j,wj);
          one_eng = resquared_lj(j,i,wj,r12,rsq,fforce,rtor,false);
          break;

         case ELLIPSE_SPHERE:
          one_eng = resquared_lj(i,j,wi,r12,rsq,fforce,ttor,true);
          tor[i][0] += ttor[0]*factor_lj;
          tor[i][1] += ttor[1]*factor_lj;
          tor[i][2] += ttor[2]*factor_lj;
          break;

         default:
          precompute_i(j,wj);
          one_eng = resquared_analytic(i,j,wi,wj,r12,rsq,fforce,ttor,rtor);
          tor[i][0] += ttor[0]*factor_lj;
          tor[i][1] += ttor[1]*factor_lj;
          tor[i][2] += ttor[2]*factor_lj;

         break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];

        if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz_full(i,evdwl,0.0,fforce[0],fforce[1],
                                      fforce[2],-r12[0],-r12[1],-r12[2]);
      }
    }
  }
}
