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

/* ----------------------------------------------------------------------
   Contributing author: Eddy BARRAUD (IFPEN/Sorbonne)
------------------------------------------------------------------------- */

#include "pair_dpd_coul_slater_long_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include "ewald_const.h"
#include "kspace.h"

#include <cmath>
#include "lammps_gpu.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;
using namespace EwaldConst;


/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLongGPU::PairDPDCoulSlaterLongGPU(LAMMPS *lmp) : PairDPDCoulSlaterLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairDPDCoulSlaterLongGPU::~PairDPDCoulSlaterLongGPU()
{
  dpd_coul_slater_long_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairDPDCoulSlaterLongGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  double dtinvsqrt = 1.0 / sqrt(update->dt);

  bool success = true;
  int *ilist, *numneigh, **firstneigh;

  double *q = atom->q;
  dpd_coul_slater_long_gpu_get_extra_data(q);

  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = dpd_coul_slater_long_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
        cpu_time, success, atom->v, dtinvsqrt, seed, update->ntimestep, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    dpd_coul_slater_long_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                    eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time, success, atom->tag,
                    atom->v, dtinvsqrt, seed, update->ntimestep, atom->nlocal, domain->boxlo,
                    domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
  if (host_start < inum) {
    cpu_time = platform::walltime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = platform::walltime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLongGPU::init_style()
{

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd requires ghost atoms store velocity");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/slater/long requires atom attribute q");
  // ensure use of KSpace long-range solver, set g_ewald
  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double mcut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mcut = init_one(i, j);
        mcut *= mcut;
        if (mcut > maxcut) maxcut = mcut;
        cutsq[i][j] = cutsq[j][i] = mcut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      dpd_coul_slater_long_gpu_init(atom->ntypes + 1, cutsq, a0, gamma, sigma,
                   cut_dpd, cut_dpdsq, cut_slatersq, force->special_lj, atom->nlocal,
                   atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen,
                   force->special_coul, force->qqrd2e, g_ewald, lamda);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);

}

/* ---------------------------------------------------------------------- */

double PairDPDCoulSlaterLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + dpd_coul_slater_long_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairDPDCoulSlaterLongGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */, int *ilist,
                             int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
  double r2inv,forcedpd,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  double rsq,r,rinv,dot,wd,randnum,factor_dpd,factor_sqrt;
  int *jlist;
  double slater_term;

  tagint itag, jtag;

  double *q = atom->q;
  double *special_coul = force->special_coul;
  double qqrd2e = force->qqrd2e;

  evdwl = 0.0;
  ecoul = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *special_lj = force->special_lj;
  double dtinvsqrt = 1.0 / sqrt(update->dt);
  int timestep = (int) update->ntimestep;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      factor_sqrt = special_sqrt[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jtag = tag[j];

      // forces if below maximum cutoff
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (evflag) evdwl = ecoul = 0.0;
        // apply DPD force if distance below DPD cutoff
        if (rsq < cut_dpdsq[itype][jtype] && r > EPSILON ) {
          rinv = 1.0 / r;
          delvx = vxtmp - v[j][0];
          delvy = vytmp - v[j][1];
          delvz = vztmp - v[j][2];
          dot = delx * delvx + dely * delvy + delz * delvz;
          wd = 1.0 - r / cut[itype][jtype];

          unsigned int tag1 = itag, tag2 = jtag;
          if (tag1 > tag2) {
            tag1 = jtag;
            tag2 = itag;
          }

          randnum = 0.0;
          saru(tag1, tag2, seed, timestep, randnum);

          // conservative force = a0 * wd
          // drag force = -gamma * wd^2 * (delx dot delv) / r
          // random force = sigma * wd * rnd * dtinvsqrt;
          forcedpd = a0[itype][jtype]*wd;
          forcedpd -= gamma[itype][jtype]*wd*wd*dot*rinv;
          forcedpd *= factor_dpd;
          forcedpd += factor_sqrt*sigma[itype][jtype]*wd*randnum*dtinvsqrt;
          forcedpd *= rinv;

          if (eflag) {
            // eng shifted to 0.0 at cutoff
            evdwl = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype] * wd*wd;
            evdwl *= factor_dpd;
          }

        } else forcedpd = 0.0;

        // apply Slater electrostatic force if distance below Slater cutoff
        // and the two species are charged
        if (rsq < cut_slatersq[itype][jtype]){
          r2inv = 1.0/rsq;
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          slater_term = exp(-2*r/lamda)*(1 + (2*r/lamda*(1+r/lamda)));
          prefactor = qqrd2e * qtmp*q[j]/r;
          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - slater_term);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor*(1-slater_term);
          forcecoul *= r2inv;

          if (eflag) {
            ecoul = prefactor*(erfc - (1 + r/lamda)*exp(-2*r/lamda));
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor*(1.0-(1 + r/lamda)*exp(-2*r/lamda));
          }

        } else forcecoul = 0.0;

        fpair = forcedpd + forcecoul;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (evflag) ev_tally_full(i, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }
}
