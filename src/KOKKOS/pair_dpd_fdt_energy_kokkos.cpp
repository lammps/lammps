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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_dpd_fdt_energy_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "comm.h"
#include "update.h"
#include "fix.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDfdtEnergyKokkos<DeviceType>::PairDPDfdtEnergyKokkos(LAMMPS *lmp) : PairDPDfdtEnergy(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | TAG_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDPDfdtEnergyKokkos<DeviceType>::~PairDPDfdtEnergyKokkos()
{
  if (allocated) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->destroy_kokkos(k_vatom,vatom);
    k_cutsq = DAT::tdual_ffloat_2d();
    memory->sfree(cutsq);
    eatom = NULL;
    vatom = NULL;
    cutsq = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memory->destroy_kokkos(k_vatom,vatom);
    memory->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairDPDfdtEnergyKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairDPDfdtEnergyKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r = sqrt(rsq);
 if (r < EPSILON) return 0;     // r can be 0.0 in DPD systems
  const F_FLOAT rinv = 1.0/r;
  const F_FLOAT wr = 1.0 - r/cut[itype][jtype];
  const F_FLOAT wd = wr*wr;

 // conservative force = a0 * wr
  return  a0[itype][jtype]*wr*rinv;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairDPDfdtEnergyKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT r = sqrt(rsq);
  if (r < EPSILON) return 0;     // r can be 0.0 in DPD systems
  const F_FLOAT rinv = 1.0/r;
  const F_FLOAT wr = 1.0 - r/cut[itype][jtype];
  const F_FLOAT wd = wr*wr;
     // unshifted eng of conservative term:
     // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
     // eng shifted to 0.0 at cutoff
  return 0.5*a0[itype][jtype]*cut[itype][jtype] * wd;
}


/*
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rinv,wd,wr,factor_dpd;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        rinv = 1.0/r;
        wr = 1.0 - r/cut[itype][jtype];
        wd = wr*wr;

        // conservative force = a0 * wr
        fpair = a0[itype][jtype]*wr;
        fpair *= factor_dpd*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}
*/

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::allocate()
{
  PairDPDfdtEnergy::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memory->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  PairDPDfdtEnergy::settings(2,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDPDfdtEnergyKokkos<DeviceType>::init_style()
{
  PairDPDfdtEnergy::init_style();

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with dpd/fdt/energy/kk");
  }

/*
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd/fdt/energy requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair dpd/fdt/energy requires newton pair on");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->ssa = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"shardlow") == 0)
      neighbor->requests[irequest]->ssa = 1;

  bool eos_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"eos",3) == 0) eos_flag = true;
  if(!eos_flag) error->all(FLERR,"pair_style dpd/fdt/energy requires an EOS to be specified");
*/
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDPDfdtEnergyKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairDPDfdtEnergy::init_one(i,j);

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }
  k_cutsq.h_view(i,j) = cutone*cutone;
  k_cutsq.h_view(j,i) = k_cutsq.h_view(i,j);
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}


namespace LAMMPS_NS {
template class PairDPDfdtEnergyKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairDPDfdtEnergyKokkos<LMPHostType>;
#endif
}

