/* -------------------------------------------------------------------------
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
   Contributing author: Axel Kohlmeyer (Temple U)
   OpenMP based threading support for LAMMPS
------------------------------------------------------------------------- */

#include "thr_omp.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "math_const.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"

#include <cstring>

using namespace LAMMPS_NS;
using MathConst::THIRD;

/* ---------------------------------------------------------------------- */

ThrOMP::ThrOMP(LAMMPS *ptr, int style) : lmp(ptr), fix(nullptr), thr_style(style), thr_error(0)
{
  // register fix omp with this class
  fix = static_cast<FixOMP *>(lmp->modify->get_fix_by_id("package_omp"));
  if (!fix) lmp->error->all(FLERR, "The 'package omp' command is required for /omp styles");
}

// clang-format off
/* ----------------------------------------------------------------------
   Hook up per thread per atom arrays into the tally infrastructure
   ---------------------------------------------------------------------- */

void ThrOMP::ev_setup_thr(int eflag, int vflag, int nall, double *eatom,
                          double **vatom, double **cvatom, ThrData *thr)
{
  const int tid = thr->get_tid();
  if (tid == 0) thr_error = 0;

  if (thr_style & THR_PAIR) {
    if (eflag & ENERGY_ATOM) {
      thr->eatom_pair = eatom + tid*nall;
      if (nall > 0)
        memset(&(thr->eatom_pair[0]),0,nall*sizeof(double));
    }
    // per-atom virial and per-atom centroid virial are the same for two-body
    // many-body pair styles not yet implemented
    if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
      thr->vatom_pair = vatom + tid*nall;
      if (nall > 0)
        memset(&(thr->vatom_pair[0][0]),0,nall*6*sizeof(double));
    }
    // check cvatom_pair, because can't access centroidstressflag
    if ((vflag & VIRIAL_CENTROID) && cvatom) {
      thr->cvatom_pair = cvatom + tid*nall;
      if (nall > 0)
        memset(&(thr->cvatom_pair[0][0]),0,nall*9*sizeof(double));
    } else {
      thr->cvatom_pair = nullptr;
    }

  }

  if (thr_style & THR_BOND) {
    if (eflag & ENERGY_ATOM) {
      thr->eatom_bond = eatom + tid*nall;
      if (nall > 0)
        memset(&(thr->eatom_bond[0]),0,nall*sizeof(double));
    }
    // per-atom virial and per-atom centroid virial are the same for bonds
    if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
      thr->vatom_bond = vatom + tid*nall;
      if (nall > 0)
        memset(&(thr->vatom_bond[0][0]),0,nall*6*sizeof(double));
    }
  }

  if (thr_style & THR_ANGLE) {
    if (eflag & ENERGY_ATOM) {
      thr->eatom_angle = eatom + tid*nall;
      if (nall > 0)
        memset(&(thr->eatom_angle[0]),0,nall*sizeof(double));
    }
    if (vflag & VIRIAL_ATOM) {
      thr->vatom_angle = vatom + tid*nall;
      if (nall > 0)
        memset(&(thr->vatom_angle[0][0]),0,nall*6*sizeof(double));
    }
    if (vflag & VIRIAL_CENTROID) {
      thr->cvatom_angle = cvatom + tid*nall;
      if (nall > 0)
        memset(&(thr->cvatom_angle[0][0]),0,nall*9*sizeof(double));
    }
  }

  if (thr_style & THR_DIHEDRAL) {
    if (eflag & ENERGY_ATOM) {
      thr->eatom_dihed = eatom + tid*nall;
      if (nall > 0)
        memset(&(thr->eatom_dihed[0]),0,nall*sizeof(double));
    }
    if (vflag & VIRIAL_ATOM) {
      thr->vatom_dihed = vatom + tid*nall;
      if (nall > 0)
        memset(&(thr->vatom_dihed[0][0]),0,nall*6*sizeof(double));
    }
    if (vflag & VIRIAL_CENTROID) {
      thr->cvatom_dihed = cvatom + tid*nall;
      if (nall > 0)
        memset(&(thr->cvatom_dihed[0][0]),0,nall*9*sizeof(double));
    }
  }

  if (thr_style & THR_IMPROPER) {
    if (eflag & ENERGY_ATOM) {
      thr->eatom_imprp = eatom + tid*nall;
      if (nall > 0)
        memset(&(thr->eatom_imprp[0]),0,nall*sizeof(double));
    }
    if (vflag & VIRIAL_ATOM) {
      thr->vatom_imprp = vatom + tid*nall;
      if (nall > 0)
        memset(&(thr->vatom_imprp[0][0]),0,nall*6*sizeof(double));
    }
    if (vflag & VIRIAL_CENTROID) {
      thr->cvatom_imprp = cvatom + tid*nall;
      if (nall > 0)
        memset(&(thr->cvatom_imprp[0][0]),0,nall*9*sizeof(double));
    }
  }
  // nothing to do for THR_KSPACE
}

/* ----------------------------------------------------------------------
   Reduce per thread data into the regular structures
   Reduction of global properties is serialized with a "critical"
   directive, so that only one thread at a time will access the
   global variables. Since we are not synchronized, this should
   come with little overhead. The reduction of per-atom properties
   in contrast is parallelized over threads in the same way as forces.
   ---------------------------------------------------------------------- */

void ThrOMP::reduce_thr(void *style, const int eflag, const int vflag,
                        ThrData *const thr)
{
  const int nlocal = lmp->atom->nlocal;
  const int nghost = lmp->atom->nghost;
  const int nall = nlocal + nghost;
  const int nfirst = lmp->atom->nfirst;
  const int nthreads = lmp->comm->nthreads;
  const int evflag = eflag | vflag;

  const int tid = thr->get_tid();
  double **f = lmp->atom->f;
  double **x = lmp->atom->x;

  int need_force_reduce = 1;

  if (evflag)
    sync_threads();

  switch (thr_style) {

  case THR_PAIR: {

    if (lmp->force->pair->vflag_fdotr) {

      // this is a non-hybrid pair style. compute per thread fdotr
      if (fix->last_pair_hybrid == nullptr) {
        if (lmp->neighbor->includegroup == 0)
          thr->virial_fdotr_compute(x, nlocal, nghost, -1);
        else
          thr->virial_fdotr_compute(x, nlocal, nghost, nfirst);
      } else {
        if (style == fix->last_pair_hybrid) {
          // pair_style hybrid will compute fdotr for us
          // but we first need to reduce the forces
          data_reduce_thr(&(f[0][0]), nall, nthreads, 3, tid);
          fix->did_reduce();
          need_force_reduce = 0;
        }
      }
    }

    if (evflag) {
      auto  const pair = (Pair *)style;

#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          pair->eng_vdwl += thr->eng_vdwl;
          pair->eng_coul += thr->eng_coul;
          thr->eng_vdwl = 0.0;
          thr->eng_coul = 0.0;
        }
        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR))
          for (int i=0; i < 6; ++i) {
            pair->virial[i] += thr->virial_pair[i];
            thr->virial_pair[i] = 0.0;
          }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(pair->eatom[0]), nall, nthreads, 1, tid);
      }
      // per-atom virial and per-atom centroid virial are the same for two-body
      // many-body pair styles not yet implemented
      if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
        data_reduce_thr(&(pair->vatom[0][0]), nall, nthreads, 6, tid);
      }
      // check cvatom_pair, because can't access centroidstressflag
      if ((vflag & VIRIAL_CENTROID) && thr->cvatom_pair) {
        data_reduce_thr(&(pair->cvatom[0][0]), nall, nthreads, 9, tid);
      }
    }
  }
    break;

  case THR_BOND:

    if (evflag) {
      Bond * const bond = lmp->force->bond;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          bond->energy += thr->eng_bond;
          thr->eng_bond = 0.0;
        }

        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR)) {
          for (int i=0; i < 6; ++i) {
            bond->virial[i] += thr->virial_bond[i];
            thr->virial_bond[i] = 0.0;
          }
        }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(bond->eatom[0]), nall, nthreads, 1, tid);
      }
      // per-atom virial and per-atom centroid virial are the same for bonds
      if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
        data_reduce_thr(&(bond->vatom[0][0]), nall, nthreads, 6, tid);
      }

    }
    break;

  case THR_ANGLE:

    if (evflag) {
      Angle * const angle = lmp->force->angle;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          angle->energy += thr->eng_angle;
          thr->eng_angle = 0.0;
        }

        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR)) {
          for (int i=0; i < 6; ++i) {
            angle->virial[i] += thr->virial_angle[i];
            thr->virial_angle[i] = 0.0;
          }
        }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(angle->eatom[0]), nall, nthreads, 1, tid);
      }
      if (vflag & VIRIAL_ATOM) {
        data_reduce_thr(&(angle->vatom[0][0]), nall, nthreads, 6, tid);
      }
      if (vflag & VIRIAL_CENTROID) {
        data_reduce_thr(&(angle->cvatom[0][0]), nall, nthreads, 9, tid);
      }

    }
    break;

  case THR_DIHEDRAL:

    if (evflag) {
      Dihedral * const dihedral = lmp->force->dihedral;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          dihedral->energy += thr->eng_dihed;
          thr->eng_dihed = 0.0;
        }

        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR)) {
          for (int i=0; i < 6; ++i) {
            dihedral->virial[i] += thr->virial_dihed[i];
            thr->virial_dihed[i] = 0.0;
          }
        }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(dihedral->eatom[0]), nall, nthreads, 1, tid);
      }
      if (vflag & VIRIAL_ATOM) {
        data_reduce_thr(&(dihedral->vatom[0][0]), nall, nthreads, 6, tid);
      }
      if (vflag & VIRIAL_CENTROID) {
        data_reduce_thr(&(dihedral->cvatom[0][0]), nall, nthreads, 9, tid);
      }

    }
    break;

  case THR_DIHEDRAL|THR_CHARMM: // special case for CHARMM dihedrals

    if (evflag) {
      Dihedral * const dihedral = lmp->force->dihedral;
      Pair * const pair = lmp->force->pair;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          dihedral->energy += thr->eng_dihed;
          pair->eng_vdwl += thr->eng_vdwl;
          pair->eng_coul += thr->eng_coul;
          thr->eng_dihed = 0.0;
          thr->eng_vdwl = 0.0;
          thr->eng_coul = 0.0;
        }

        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR)) {
          for (int i=0; i < 6; ++i) {
            dihedral->virial[i] += thr->virial_dihed[i];
            pair->virial[i] += thr->virial_pair[i];
            thr->virial_dihed[i] = 0.0;
            thr->virial_pair[i] = 0.0;
          }
        }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(dihedral->eatom[0]), nall, nthreads, 1, tid);
        data_reduce_thr(&(pair->eatom[0]), nall, nthreads, 1, tid);
      }
      if (vflag & VIRIAL_ATOM) {
        data_reduce_thr(&(dihedral->vatom[0][0]), nall, nthreads, 6, tid);
      }
      if (vflag & VIRIAL_CENTROID) {
        data_reduce_thr(&(dihedral->cvatom[0][0]), nall, nthreads, 9, tid);
      }
      // per-atom virial and per-atom centroid virial are the same for two-body
      // many-body pair styles not yet implemented
      if (vflag & (VIRIAL_ATOM | VIRIAL_CENTROID)) {
        data_reduce_thr(&(pair->vatom[0][0]), nall, nthreads, 6, tid);
      }
      // check cvatom_pair, because can't access centroidstressflag
      if ((vflag & VIRIAL_CENTROID) && thr->cvatom_pair) {
        data_reduce_thr(&(pair->cvatom[0][0]), nall, nthreads, 9, tid);
      }
    }
    break;

  case THR_IMPROPER:

    if (evflag) {
      Improper *improper = lmp->force->improper;
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        if (eflag & ENERGY_GLOBAL) {
          improper->energy += thr->eng_imprp;
          thr->eng_imprp = 0.0;
        }

        if (vflag & (VIRIAL_PAIR | VIRIAL_FDOTR)) {
          for (int i=0; i < 6; ++i) {
            improper->virial[i] += thr->virial_imprp[i];
            thr->virial_imprp[i] = 0.0;
          }
        }
      }

      if (eflag & ENERGY_ATOM) {
        data_reduce_thr(&(improper->eatom[0]), nall, nthreads, 1, tid);
      }
      if (vflag & VIRIAL_ATOM) {
        data_reduce_thr(&(improper->vatom[0][0]), nall, nthreads, 6, tid);
      }
      if (vflag & VIRIAL_CENTROID) {
        data_reduce_thr(&(improper->cvatom[0][0]), nall, nthreads, 9, tid);
      }

    }
    break;

  case THR_KSPACE:
    // nothing to do. XXX may need to add support for per-atom info
    break;

  case THR_INTGR:
    // nothing to do
    break;

  default:
    printf("tid:%d unhandled thr_style case %d\n", tid, thr_style);
    break;
  }

  if (style == fix->last_omp_style) {
    if (need_force_reduce) {
      data_reduce_thr(&(f[0][0]), nall, nthreads, 3, tid);
      fix->did_reduce();
    }

    if (lmp->atom->torque)
      data_reduce_thr(&(lmp->atom->torque[0][0]), nall, nthreads, 3, tid);
  }
  thr->timer(Timer::COMM);
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and eng_coul into per thread global and per-atom accumulators
------------------------------------------------------------------------- */

void ThrOMP::e_tally_thr(Pair * const pair, const int i, const int j,
                         const int nlocal, const int newton_pair,
                         const double evdwl, const double ecoul, ThrData * const thr)
{
  if (pair->eflag_global) {
    if (newton_pair) {
      thr->eng_vdwl += evdwl;
      thr->eng_coul += ecoul;
    } else {
      const double evdwlhalf = 0.5*evdwl;
      const double ecoulhalf = 0.5*ecoul;
      if (i < nlocal) {
        thr->eng_vdwl += evdwlhalf;
        thr->eng_coul += ecoulhalf;
      }
      if (j < nlocal) {
        thr->eng_vdwl += evdwlhalf;
        thr->eng_coul += ecoulhalf;
      }
    }
  }
  if (pair->eflag_atom && thr->eatom_pair) {
    const double epairhalf = 0.5 * (evdwl + ecoul);
    if (newton_pair || i < nlocal) thr->eatom_pair[i] += epairhalf;
    if (newton_pair || j < nlocal) thr->eatom_pair[j] += epairhalf;
  }
}

/* helper functions */
static void v_tally(double * const vout, const double * const vin)
{
  vout[0] += vin[0];
  vout[1] += vin[1];
  vout[2] += vin[2];
  vout[3] += vin[3];
  vout[4] += vin[4];
  vout[5] += vin[5];
}

static void v_tally9(double * const vout, const double * const vin)
{
  vout[0] += vin[0];
  vout[1] += vin[1];
  vout[2] += vin[2];
  vout[3] += vin[3];
  vout[4] += vin[4];
  vout[5] += vin[5];
  vout[6] += vin[6];
  vout[7] += vin[7];
  vout[8] += vin[8];
}

static void v_tally(double * const vout, const double scale, const double * const vin)
{
  vout[0] += scale*vin[0];
  vout[1] += scale*vin[1];
  vout[2] += scale*vin[2];
  vout[3] += scale*vin[3];
  vout[4] += scale*vin[4];
  vout[5] += scale*vin[5];
}

/* ----------------------------------------------------------------------
   tally virial into per thread global and per-atom accumulators
------------------------------------------------------------------------- */
void ThrOMP::v_tally_thr(Pair * const pair, const int i, const int j,
                         const int nlocal, const int newton_pair,
                         const double * const v, ThrData * const thr)
{
  if (pair->vflag_global) {
    double * const va = thr->virial_pair;
    if (newton_pair) {
      v_tally(va,v);
    } else {
      if (i < nlocal) v_tally(va,0.5,v);
      if (j < nlocal) v_tally(va,0.5,v);
    }
  }

  if (pair->vflag_atom) {
    if (newton_pair || i < nlocal) {
      double * const va = thr->vatom_pair[i];
      v_tally(va,0.5,v);
    }
    if (newton_pair || j < nlocal) {
      double * const va = thr->vatom_pair[j];
      v_tally(va,0.5,v);
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per thread global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Pair * const pair, const int i, const int j, const int nlocal,
                          const int newton_pair, const double evdwl, const double ecoul,
                          const double fpair, const double delx, const double dely,
                          const double delz, ThrData * const thr)
{

  if (pair->eflag_either)
    e_tally_thr(pair, i, j, nlocal, newton_pair, evdwl, ecoul, thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    v_tally_thr(pair, i, j, nlocal, newton_pair, v, thr);
  }

  if (pair->num_tally_compute > 0) {
    // ev_tally callbacks are not thread safe and thus have to be protected
#if defined(_OPENMP)
#pragma omp critical
#endif
    for (int k=0; k < pair->num_tally_compute; ++k) {
      Compute *c = pair->list_tally_compute[k];
      c->pair_tally_callback(i, j, nlocal, newton_pair,
                             evdwl, ecoul, fpair, delx, dely, delz);
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per thread global and per-atom accumulators
   for full neighbor list. Only tally atom i and also set newton to off.
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_full_thr(Pair * const pair, const int i, const double evdwl,
                               const double ecoul, const double fpair, const double delx,
                               const double dely, const double delz, ThrData * const thr)
{
  if (pair->eflag_either)
    e_tally_thr(pair, i, /*j*/ i+1, /*nlocal*/ i+1, /*newton_pair*/ 0, evdwl, ecoul, thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    v_tally_thr(pair, i, /*j*/ i+1, /*nlocal*/ i+1, /*newton_pair*/ 0, v, thr);
  }

  if (pair->num_tally_compute > 0) {
    // ev_tally callbacks are not thread safe and thus have to be protected
#if defined(_OPENMP)
#pragma omp critical
#endif
    for (int k=0; k < pair->num_tally_compute; ++k) {
      Compute *c = pair->list_tally_compute[k];
      c->pair_tally_callback(i, i+1, i, 0, evdwl, ecoul, fpair, delx, dely, delz);
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_xyz_thr(Pair * const pair, const int i, const int j,
                              const int nlocal, const int newton_pair,
                              const double evdwl, const double ecoul,
                              const double fx, const double fy, const double fz,
                              const double delx, const double dely, const double delz,
                              ThrData * const thr)
{

  if (pair->eflag_either)
    e_tally_thr(pair, i, j, nlocal, newton_pair, evdwl, ecoul, thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = delx*fx;
    v[1] = dely*fy;
    v[2] = delz*fz;
    v[3] = delx*fy;
    v[4] = delx*fz;
    v[5] = dely*fz;

    v_tally_thr(pair, i, j, nlocal, newton_pair, v, thr);
  }
}


/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   for virial, have delx,dely,delz and fx,fy,fz
   called when using full neighbor lists
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_xyz_full_thr(Pair * const pair, const int i,
                                   const double evdwl, const double ecoul,
                                   const double fx, const double fy,
                                   const double fz, const double delx,
                                   const double dely, const double delz,
                                   ThrData * const thr)
{

  if (pair->eflag_either)
    e_tally_thr(pair,i,i,i+1,0,0.5*evdwl,ecoul,thr);

  if (pair->vflag_either) {
    double v[6];
    v[0] = 0.5*delx*fx;
    v[1] = 0.5*dely*fy;
    v[2] = 0.5*delz*fz;
    v[3] = 0.5*delx*fy;
    v[4] = 0.5*delx*fz;
    v[5] = 0.5*dely*fz;

    v_tally_thr(pair,i,i,i+1,0,v,thr);
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally3_thr(Pair * const pair, const int i, const int j, const int k,
                           const double evdwl, const double ecoul,
                           const double * const fj, const double * const fk,
                           const double * const drji, const double * const drki,
                           ThrData * const thr)
{
  if (pair->eflag_either) {
    if (pair->eflag_global) {
      thr->eng_vdwl += evdwl;
      thr->eng_coul += ecoul;
    }
    if (pair->eflag_atom) {
      const double epairthird = THIRD * (evdwl + ecoul);
      thr->eatom_pair[i] += epairthird;
      thr->eatom_pair[j] += epairthird;
      thr->eatom_pair[k] += epairthird;
    }
  }

  if (pair->vflag_either) {
    double v[6];

    v[0] = drji[0]*fj[0] + drki[0]*fk[0];
    v[1] = drji[1]*fj[1] + drki[1]*fk[1];
    v[2] = drji[2]*fj[2] + drki[2]*fk[2];
    v[3] = drji[0]*fj[1] + drki[0]*fk[1];
    v[4] = drji[0]*fj[2] + drki[0]*fk[2];
    v[5] = drji[1]*fj[2] + drki[1]*fk[2];

    if (pair->vflag_global) v_tally(thr->virial_pair,v);

    if (pair->vflag_atom) {
      v_tally(thr->vatom_pair[i],THIRD,v);
      v_tally(thr->vatom_pair[j],THIRD,v);
      v_tally(thr->vatom_pair[k],THIRD,v);
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by AIREBO potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally4_thr(Pair * const pair, const int i, const int j,
                           const int k, const int m, const double evdwl,
                           const double * const fi, const double * const fj,
                           const double * const fk, const double * const drim,
                           const double * const drjm, const double * const drkm,
                           ThrData * const thr)
{
  double v[6];

  if (pair->eflag_either) {
    if (pair->eflag_global) thr->eng_vdwl += evdwl;
    if (pair->eflag_atom) {
      const double epairfourth = 0.25 * evdwl;
      thr->eatom_pair[i] += epairfourth;
      thr->eatom_pair[j] += epairfourth;
      thr->eatom_pair[k] += epairfourth;
      thr->eatom_pair[m] += epairfourth;
    }
  }

  if (pair->vflag_either) {
    v[0] = (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
    v[1] = (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
    v[2] = (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
    v[3] = (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
    v[4] = (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
    v[5] = (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);
    if (pair->vflag_global) v_tally(thr->virial_pair,v);

    if (pair->vflag_atom) {
      v[0] *= 0.25;
      v[1] *= 0.25;
      v[2] *= 0.25;
      v[3] *= 0.25;
      v[4] *= 0.25;
      v[5] *= 0.25;
      v_tally(thr->vatom_pair[i],v);
      v_tally(thr->vatom_pair[j],v);
      v_tally(thr->vatom_pair[k],v);
      v_tally(thr->vatom_pair[m],v);
    }
  }
}

/* ----------------------------------------------------------------------
   tally ecoul and virial into each of n atoms in list
   called by TIP4P potential, newton_pair is always on
   changes v values by dividing by n
 ------------------------------------------------------------------------- */

void ThrOMP::ev_tally_list_thr(Pair * const pair, const int key,
                               const int * const list, const double * const v,
                               const double ecoul, const double alpha,
                               ThrData * const thr)
{
  int i;
  if (pair->eflag_either) {
    if (pair->eflag_global) thr->eng_coul += ecoul;
    if (pair->eflag_atom) {
      if (key == 0) {
        thr->eatom_pair[list[0]] += 0.5*ecoul;
        thr->eatom_pair[list[1]] += 0.5*ecoul;
      } else if (key == 1) {
        thr->eatom_pair[list[0]] += 0.5*ecoul*(1-alpha);
        thr->eatom_pair[list[1]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[2]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[3]] += 0.5*ecoul;
      } else if (key == 2) {
        thr->eatom_pair[list[0]] += 0.5*ecoul;
        thr->eatom_pair[list[1]] += 0.5*ecoul*(1-alpha);
        thr->eatom_pair[list[2]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[3]] += 0.25*ecoul*alpha;
      } else {
        thr->eatom_pair[list[0]] += 0.5*ecoul*(1-alpha);
        thr->eatom_pair[list[1]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[2]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[3]] += 0.5*ecoul*(1-alpha);
        thr->eatom_pair[list[4]] += 0.25*ecoul*alpha;
        thr->eatom_pair[list[5]] += 0.25*ecoul*alpha;
      }
    }
  }

  if (pair->vflag_either) {
    if (pair->vflag_global)
      v_tally(thr->virial_pair,v);

    if (pair->vflag_atom) {
      if (key == 0) {
        for (i = 0; i <= 5; i++) {
          thr->vatom_pair[list[0]][i] += 0.5*v[i];
          thr->vatom_pair[list[1]][i] += 0.5*v[i];
        }
      } else if (key == 1) {
        for (i = 0; i <= 5; i++) {
          thr->vatom_pair[list[0]][i] += 0.5*v[i]*(1-alpha);
          thr->vatom_pair[list[1]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[2]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[3]][i] += 0.5*v[i];
        }
      } else if (key == 2) {
        for (i = 0; i <= 5; i++) {
          thr->vatom_pair[list[0]][i] += 0.5*v[i];
          thr->vatom_pair[list[1]][i] += 0.5*v[i]*(1-alpha);
          thr->vatom_pair[list[2]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[3]][i] += 0.25*v[i]*alpha;
        }
      } else {
        for (i = 0; i <= 5; i++) {
          thr->vatom_pair[list[0]][i] += 0.5*v[i]*(1-alpha);
          thr->vatom_pair[list[1]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[2]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[3]][i] += 0.5*v[i]*(1-alpha);
          thr->vatom_pair[list[4]][i] += 0.25*v[i]*alpha;
          thr->vatom_pair[list[5]][i] += 0.25*v[i]*alpha;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Bond * const bond, const int i, const int j, const int nlocal,
                          const int newton_bond, const double ebond, const double fbond,
                          const double delx, const double dely, const double delz,
                          ThrData * const thr)
{
  if (bond->eflag_either) {
    const double ebondhalf = 0.5*ebond;
    if (newton_bond) {
      if (bond->eflag_global)
        thr->eng_bond += ebond;
      if (bond->eflag_atom) {
        thr->eatom_bond[i] += ebondhalf;
        thr->eatom_bond[j] += ebondhalf;
      }
    } else {
      if (bond->eflag_global) {
        if (i < nlocal) thr->eng_bond += ebondhalf;
        if (j < nlocal) thr->eng_bond += ebondhalf;
      }
      if (bond->eflag_atom) {
        if (i < nlocal) thr->eatom_bond[i] += ebondhalf;
        if (j < nlocal) thr->eatom_bond[j] += ebondhalf;
      }
    }
  }

  if (bond->vflag_either) {
    double v[6];

    v[0] = delx*delx*fbond;
    v[1] = dely*dely*fbond;
    v[2] = delz*delz*fbond;
    v[3] = delx*dely*fbond;
    v[4] = delx*delz*fbond;
    v[5] = dely*delz*fbond;

    if (bond->vflag_global) {
      if (newton_bond)
        v_tally(thr->virial_bond,v);
      else {
        if (i < nlocal)
          v_tally(thr->virial_bond,0.5,v);
        if (j < nlocal)
          v_tally(thr->virial_bond,0.5,v);
      }
    }

    if (bond->vflag_atom) {
      v[0] *= 0.5;
      v[1] *= 0.5;
      v[2] *= 0.5;
      v[3] *= 0.5;
      v[4] *= 0.5;
      v[5] *= 0.5;

      if (newton_bond) {
        v_tally(thr->vatom_bond[i],v);
        v_tally(thr->vatom_bond[j],v);
      } else {
        if (i < nlocal)
          v_tally(thr->vatom_bond[i],v);
        if (j < nlocal)
          v_tally(thr->vatom_bond[j],v);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Angle * const angle, const int i, const int j, const int k,
                          const int nlocal, const int newton_bond, const double eangle,
                          const double * const f1, const double * const f3,
                          const double delx1, const double dely1, const double delz1,
                          const double delx2, const double dely2, const double delz2,
                          ThrData * const thr)
{
  if (angle->eflag_either) {
    const double eanglethird = THIRD*eangle;
    if (newton_bond) {
      if (angle->eflag_global)
        thr->eng_angle += eangle;
      if (angle->eflag_atom) {
        thr->eatom_angle[i] += eanglethird;
        thr->eatom_angle[j] += eanglethird;
        thr->eatom_angle[k] += eanglethird;
      }
    } else {
      if (angle->eflag_global) {
        if (i < nlocal) thr->eng_angle += eanglethird;
        if (j < nlocal) thr->eng_angle += eanglethird;
        if (k < nlocal) thr->eng_angle += eanglethird;
      }
      if (angle->eflag_atom) {
        if (i < nlocal) thr->eatom_angle[i] += eanglethird;
        if (j < nlocal) thr->eatom_angle[j] += eanglethird;
        if (k < nlocal) thr->eatom_angle[k] += eanglethird;
      }
    }
  }

  if (angle->vflag_either) {
    double v[6];

    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (angle->vflag_global) {
      if (newton_bond) {
        v_tally(thr->virial_angle,v);
      } else {
        int cnt = 0;
        if (i < nlocal) ++cnt;
        if (j < nlocal) ++cnt;
        if (k < nlocal) ++cnt;
        v_tally(thr->virial_angle,cnt*THIRD,v);
      }
    }

    if (angle->vflag_atom) {
      v[0] *= THIRD;
      v[1] *= THIRD;
      v[2] *= THIRD;
      v[3] *= THIRD;
      v[4] *= THIRD;
      v[5] *= THIRD;

      if (newton_bond) {
        v_tally(thr->vatom_angle[i],v);
        v_tally(thr->vatom_angle[j],v);
        v_tally(thr->vatom_angle[k],v);
      } else {
        if (i < nlocal) v_tally(thr->vatom_angle[i],v);
        if (j < nlocal) v_tally(thr->vatom_angle[j],v);
        if (k < nlocal) v_tally(thr->vatom_angle[k],v);
      }
    }
  }

  // per-atom centroid virial
  if (angle->cvflag_atom) {
    double f2[3], v1[9], v2[9], v3[9];
    double a1[3], a2[3], a3[3];

    // r0 = (r1+r2+r3)/3
    // rij = ri-rj
    // total virial = r10*f1 + r20*f2 + r30*f3
    // del1: r12
    // del2: r32

    // a1 = r10 = (2*r12 -   r32)/3
    a1[0] = THIRD*(2*delx1-delx2);
    a1[1] = THIRD*(2*dely1-dely2);
    a1[2] = THIRD*(2*delz1-delz2);

    // a2 = r20 = ( -r12 -   r32)/3
    a2[0] = THIRD*(-delx1-delx2);
    a2[1] = THIRD*(-dely1-dely2);
    a2[2] = THIRD*(-delz1-delz2);

    // a3 = r30 = ( -r12 + 2*r32)/3
    a3[0] = THIRD*(-delx1+2*delx2);
    a3[1] = THIRD*(-dely1+2*dely2);
    a3[2] = THIRD*(-delz1+2*delz2);

    f2[0] = - f1[0] - f3[0];
    f2[1] = - f1[1] - f3[1];
    f2[2] = - f1[2] - f3[2];

    v1[0] = a1[0]*f1[0];
    v1[1] = a1[1]*f1[1];
    v1[2] = a1[2]*f1[2];
    v1[3] = a1[0]*f1[1];
    v1[4] = a1[0]*f1[2];
    v1[5] = a1[1]*f1[2];
    v1[6] = a1[1]*f1[0];
    v1[7] = a1[2]*f1[0];
    v1[8] = a1[2]*f1[1];

    v2[0] = a2[0]*f2[0];
    v2[1] = a2[1]*f2[1];
    v2[2] = a2[2]*f2[2];
    v2[3] = a2[0]*f2[1];
    v2[4] = a2[0]*f2[2];
    v2[5] = a2[1]*f2[2];
    v2[6] = a2[1]*f2[0];
    v2[7] = a2[2]*f2[0];
    v2[8] = a2[2]*f2[1];

    v3[0] = a3[0]*f3[0];
    v3[1] = a3[1]*f3[1];
    v3[2] = a3[2]*f3[2];
    v3[3] = a3[0]*f3[1];
    v3[4] = a3[0]*f3[2];
    v3[5] = a3[1]*f3[2];
    v3[6] = a3[1]*f3[0];
    v3[7] = a3[2]*f3[0];
    v3[8] = a3[2]*f3[1];

    if (newton_bond) {
      v_tally9(thr->cvatom_angle[i],v1);
      v_tally9(thr->cvatom_angle[j],v2);
      v_tally9(thr->cvatom_angle[k],v3);
    } else {
      if (i < nlocal) v_tally9(thr->cvatom_angle[i],v1);
      if (j < nlocal) v_tally9(thr->cvatom_angle[j],v2);
      if (k < nlocal) v_tally9(thr->cvatom_angle[k],v3);
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial from 1-3 repulsion of SPICA angle into accumulators
------------------------------------------------------------------------- */

void ThrOMP::ev_tally13_thr(Angle * const angle, const int i1, const int i3,
                            const int nlocal, const int newton_bond,
                            const double epair, const double fpair,
                            const double delx, const double dely,
                            const double delz, ThrData * const thr)
{

  if (angle->eflag_either) {
    const double epairhalf = 0.5 * epair;

    if (angle->eflag_global) {
      if (newton_bond || i1 < nlocal)
        thr->eng_angle += epairhalf;
      if (newton_bond || i3 < nlocal)
        thr->eng_angle += epairhalf;
    }

    if (angle->eflag_atom) {
      if (newton_bond || i1 < nlocal) thr->eatom_angle[i1] += epairhalf;
      if (newton_bond || i3 < nlocal) thr->eatom_angle[i3] += epairhalf;
    }
  }

  if (angle->vflag_either) {
    double v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (angle->vflag_global) {
      double * const va = thr->virial_angle;
      if (newton_bond || i1 < nlocal) v_tally(va,0.5,v);
      if (newton_bond || i3 < nlocal) v_tally(va,0.5,v);
    }

    if (angle->vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        double * const va = thr->vatom_angle[i1];
        v_tally(va,0.5,v);
      }
      if (newton_bond || i3 < nlocal) {
        double * const va = thr->vatom_angle[i3];
        v_tally(va,0.5,v);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Dihedral * const dihed, const int i1, const int i2,
                          const int i3, const int i4, const int nlocal,
                          const int newton_bond, const double edihedral,
                          const double * const f1, const double * const f3,
                          const double * const f4, const double vb1x,
                          const double vb1y, const double vb1z, const double vb2x,
                          const double vb2y, const double vb2z, const double vb3x,
                          const double vb3y, const double vb3z, ThrData * const thr)
{

  if (dihed->eflag_either) {
    if (dihed->eflag_global) {
      if (newton_bond) {
        thr->eng_dihed += edihedral;
      } else {
        const double edihedralquarter = 0.25*edihedral;
        int cnt = 0;
        if (i1 < nlocal) ++cnt;
        if (i2 < nlocal) ++cnt;
        if (i3 < nlocal) ++cnt;
        if (i4 < nlocal) ++cnt;
        thr->eng_dihed += static_cast<double>(cnt)*edihedralquarter;
      }
    }
    if (dihed->eflag_atom) {
      const double edihedralquarter = 0.25*edihedral;
      if (newton_bond) {
        thr->eatom_dihed[i1] += edihedralquarter;
        thr->eatom_dihed[i2] += edihedralquarter;
        thr->eatom_dihed[i3] += edihedralquarter;
        thr->eatom_dihed[i4] += edihedralquarter;
      } else {
        if (i1 < nlocal) thr->eatom_dihed[i1] +=  edihedralquarter;
        if (i2 < nlocal) thr->eatom_dihed[i2] +=  edihedralquarter;
        if (i3 < nlocal) thr->eatom_dihed[i3] +=  edihedralquarter;
        if (i4 < nlocal) thr->eatom_dihed[i4] +=  edihedralquarter;
      }
    }
  }

  if (dihed->vflag_either) {
    double v[6];
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (dihed->vflag_global) {
      if (newton_bond) {
        v_tally(thr->virial_dihed,v);
      } else {
        int cnt = 0;
        if (i1 < nlocal) ++cnt;
        if (i2 < nlocal) ++cnt;
        if (i3 < nlocal) ++cnt;
        if (i4 < nlocal) ++cnt;
        v_tally(thr->virial_dihed,0.25*static_cast<double>(cnt),v);
      }
    }

    v[0] *= 0.25;
    v[1] *= 0.25;
    v[2] *= 0.25;
    v[3] *= 0.25;
    v[4] *= 0.25;
    v[5] *= 0.25;

    if (dihed->vflag_atom) {
      if (newton_bond) {
        v_tally(thr->vatom_dihed[i1],v);
        v_tally(thr->vatom_dihed[i2],v);
        v_tally(thr->vatom_dihed[i3],v);
        v_tally(thr->vatom_dihed[i4],v);
      } else {
        if (i1 < nlocal) v_tally(thr->vatom_dihed[i1],v);
        if (i2 < nlocal) v_tally(thr->vatom_dihed[i2],v);
        if (i3 < nlocal) v_tally(thr->vatom_dihed[i3],v);
        if (i4 < nlocal) v_tally(thr->vatom_dihed[i4],v);
      }
    }
  }

  // per-atom centroid virial
  if (dihed->cvflag_atom) {
    double f2[3], v1[9], v2[9], v3[9], v4[9];
    double a1[3], a2[3], a3[3], a4[3];

    // r0 = (r1+r2+r3+r4)/4
    // rij = ri-rj
    // total virial = r10*f1 + r20*f2 + r30*f3 + r40*f4
    // vb1: r12
    // vb2: r32
    // vb3: r43

    // a1 = r10 = (3*r12 - 2*r32 -   r43)/4
    a1[0] = 0.25*(3*vb1x - 2*vb2x - vb3x);
    a1[1] = 0.25*(3*vb1y - 2*vb2y - vb3y);
    a1[2] = 0.25*(3*vb1z - 2*vb2z - vb3z);

    // a2 = r20 = ( -r12 - 2*r32 -   r43)/4
    a2[0] = 0.25*(-vb1x - 2*vb2x - vb3x);
    a2[1] = 0.25*(-vb1y - 2*vb2y - vb3y);
    a2[2] = 0.25*(-vb1z - 2*vb2z - vb3z);

    // a3 = r30 = ( -r12 + 2*r32 -   r43)/4
    a3[0] = 0.25*(-vb1x + 2*vb2x - vb3x);
    a3[1] = 0.25*(-vb1y + 2*vb2y - vb3y);
    a3[2] = 0.25*(-vb1z + 2*vb2z - vb3z);

    // a4 = r40 = ( -r12 + 2*r32 + 3*r43)/4
    a4[0] = 0.25*(-vb1x + 2*vb2x + 3*vb3x);
    a4[1] = 0.25*(-vb1y + 2*vb2y + 3*vb3y);
    a4[2] = 0.25*(-vb1z + 2*vb2z + 3*vb3z);

    f2[0] = - f1[0] - f3[0] - f4[0];
    f2[1] = - f1[1] - f3[1] - f4[1];
    f2[2] = - f1[2] - f3[2] - f4[2];

    v1[0] = a1[0]*f1[0];
    v1[1] = a1[1]*f1[1];
    v1[2] = a1[2]*f1[2];
    v1[3] = a1[0]*f1[1];
    v1[4] = a1[0]*f1[2];
    v1[5] = a1[1]*f1[2];
    v1[6] = a1[1]*f1[0];
    v1[7] = a1[2]*f1[0];
    v1[8] = a1[2]*f1[1];

    v2[0] = a2[0]*f2[0];
    v2[1] = a2[1]*f2[1];
    v2[2] = a2[2]*f2[2];
    v2[3] = a2[0]*f2[1];
    v2[4] = a2[0]*f2[2];
    v2[5] = a2[1]*f2[2];
    v2[6] = a2[1]*f2[0];
    v2[7] = a2[2]*f2[0];
    v2[8] = a2[2]*f2[1];

    v3[0] = a3[0]*f3[0];
    v3[1] = a3[1]*f3[1];
    v3[2] = a3[2]*f3[2];
    v3[3] = a3[0]*f3[1];
    v3[4] = a3[0]*f3[2];
    v3[5] = a3[1]*f3[2];
    v3[6] = a3[1]*f3[0];
    v3[7] = a3[2]*f3[0];
    v3[8] = a3[2]*f3[1];

    v4[0] = a4[0]*f4[0];
    v4[1] = a4[1]*f4[1];
    v4[2] = a4[2]*f4[2];
    v4[3] = a4[0]*f4[1];
    v4[4] = a4[0]*f4[2];
    v4[5] = a4[1]*f4[2];
    v4[6] = a4[1]*f4[0];
    v4[7] = a4[2]*f4[0];
    v4[8] = a4[2]*f4[1];

    if (newton_bond) {
      v_tally9(thr->cvatom_dihed[i1],v1);
      v_tally9(thr->cvatom_dihed[i2],v2);
      v_tally9(thr->cvatom_dihed[i3],v3);
      v_tally9(thr->cvatom_dihed[i4],v4);
    } else {
      if (i1 < nlocal) v_tally9(thr->cvatom_dihed[i1],v1);
      if (i2 < nlocal) v_tally9(thr->cvatom_dihed[i2],v2);
      if (i3 < nlocal) v_tally9(thr->cvatom_dihed[i3],v3);
      if (i4 < nlocal) v_tally9(thr->cvatom_dihed[i4],v4);
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

void ThrOMP::ev_tally_thr(Improper * const imprp, const int i1, const int i2,
                          const int i3, const int i4, const int nlocal,
                          const int newton_bond, const double eimproper,
                          const double * const f1, const double * const f3,
                          const double * const f4, const double vb1x,
                          const double vb1y, const double vb1z, const double vb2x,
                          const double vb2y, const double vb2z, const double vb3x,
                          const double vb3y, const double vb3z, ThrData * const thr)
{

  if (imprp->eflag_either) {
    if (imprp->eflag_global) {
      if (newton_bond) {
        thr->eng_imprp += eimproper;
      } else {
        const double eimproperquarter = 0.25*eimproper;
        int cnt = 0;
        if (i1 < nlocal) ++cnt;
        if (i2 < nlocal) ++cnt;
        if (i3 < nlocal) ++cnt;
        if (i4 < nlocal) ++cnt;
        thr->eng_imprp += static_cast<double>(cnt)*eimproperquarter;
      }
    }
    if (imprp->eflag_atom) {
      const double eimproperquarter = 0.25*eimproper;
      if (newton_bond) {
        thr->eatom_imprp[i1] += eimproperquarter;
        thr->eatom_imprp[i2] += eimproperquarter;
        thr->eatom_imprp[i3] += eimproperquarter;
        thr->eatom_imprp[i4] += eimproperquarter;
      } else {
        if (i1 < nlocal) thr->eatom_imprp[i1] +=  eimproperquarter;
        if (i2 < nlocal) thr->eatom_imprp[i2] +=  eimproperquarter;
        if (i3 < nlocal) thr->eatom_imprp[i3] +=  eimproperquarter;
        if (i4 < nlocal) thr->eatom_imprp[i4] +=  eimproperquarter;
      }
    }
  }

  if (imprp->vflag_either) {
    double v[6];
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (imprp->vflag_global) {
      if (newton_bond) {
        v_tally(thr->virial_imprp,v);
      } else {
        int cnt = 0;
        if (i1 < nlocal) ++cnt;
        if (i2 < nlocal) ++cnt;
        if (i3 < nlocal) ++cnt;
        if (i4 < nlocal) ++cnt;
        v_tally(thr->virial_imprp,0.25*static_cast<double>(cnt),v);
      }
    }

    v[0] *= 0.25;
    v[1] *= 0.25;
    v[2] *= 0.25;
    v[3] *= 0.25;
    v[4] *= 0.25;
    v[5] *= 0.25;

    if (imprp->vflag_atom) {
      if (newton_bond) {
        v_tally(thr->vatom_imprp[i1],v);
        v_tally(thr->vatom_imprp[i2],v);
        v_tally(thr->vatom_imprp[i3],v);
        v_tally(thr->vatom_imprp[i4],v);
      } else {
        if (i1 < nlocal) v_tally(thr->vatom_imprp[i1],v);
        if (i2 < nlocal) v_tally(thr->vatom_imprp[i2],v);
        if (i3 < nlocal) v_tally(thr->vatom_imprp[i3],v);
        if (i4 < nlocal) v_tally(thr->vatom_imprp[i4],v);
      }
    }
  }

  // per-atom centroid virial
  if (imprp->cvflag_atom) {
    double f2[3], v1[9], v2[9], v3[9], v4[9];
    double a1[3], a2[3], a3[3], a4[3];

    // r0 = (r1+r2+r3+r4)/4
    // rij = ri-rj
    // total virial = r10*f1 + r20*f2 + r30*f3 + r40*f4
    // vb1: r12
    // vb2: r32
    // vb3: r43

    // a1 = r10 = (3*r12 - 2*r32 -   r43)/4
    a1[0] = 0.25*(3*vb1x - 2*vb2x - vb3x);
    a1[1] = 0.25*(3*vb1y - 2*vb2y - vb3y);
    a1[2] = 0.25*(3*vb1z - 2*vb2z - vb3z);

    // a2 = r20 = ( -r12 - 2*r32 -   r43)/4
    a2[0] = 0.25*(-vb1x - 2*vb2x - vb3x);
    a2[1] = 0.25*(-vb1y - 2*vb2y - vb3y);
    a2[2] = 0.25*(-vb1z - 2*vb2z - vb3z);

    // a3 = r30 = ( -r12 + 2*r32 -   r43)/4
    a3[0] = 0.25*(-vb1x + 2*vb2x - vb3x);
    a3[1] = 0.25*(-vb1y + 2*vb2y - vb3y);
    a3[2] = 0.25*(-vb1z + 2*vb2z - vb3z);

    // a4 = r40 = ( -r12 + 2*r32 + 3*r43)/4
    a4[0] = 0.25*(-vb1x + 2*vb2x + 3*vb3x);
    a4[1] = 0.25*(-vb1y + 2*vb2y + 3*vb3y);
    a4[2] = 0.25*(-vb1z + 2*vb2z + 3*vb3z);

    f2[0] = - f1[0] - f3[0] - f4[0];
    f2[1] = - f1[1] - f3[1] - f4[1];
    f2[2] = - f1[2] - f3[2] - f4[2];

    v1[0] = a1[0]*f1[0];
    v1[1] = a1[1]*f1[1];
    v1[2] = a1[2]*f1[2];
    v1[3] = a1[0]*f1[1];
    v1[4] = a1[0]*f1[2];
    v1[5] = a1[1]*f1[2];
    v1[6] = a1[1]*f1[0];
    v1[7] = a1[2]*f1[0];
    v1[8] = a1[2]*f1[1];

    v2[0] = a2[0]*f2[0];
    v2[1] = a2[1]*f2[1];
    v2[2] = a2[2]*f2[2];
    v2[3] = a2[0]*f2[1];
    v2[4] = a2[0]*f2[2];
    v2[5] = a2[1]*f2[2];
    v2[6] = a2[1]*f2[0];
    v2[7] = a2[2]*f2[0];
    v2[8] = a2[2]*f2[1];

    v3[0] = a3[0]*f3[0];
    v3[1] = a3[1]*f3[1];
    v3[2] = a3[2]*f3[2];
    v3[3] = a3[0]*f3[1];
    v3[4] = a3[0]*f3[2];
    v3[5] = a3[1]*f3[2];
    v3[6] = a3[1]*f3[0];
    v3[7] = a3[2]*f3[0];
    v3[8] = a3[2]*f3[1];

    v4[0] = a4[0]*f4[0];
    v4[1] = a4[1]*f4[1];
    v4[2] = a4[2]*f4[2];
    v4[3] = a4[0]*f4[1];
    v4[4] = a4[0]*f4[2];
    v4[5] = a4[1]*f4[2];
    v4[6] = a4[1]*f4[0];
    v4[7] = a4[2]*f4[0];
    v4[8] = a4[2]*f4[1];

    if (newton_bond) {
      v_tally9(thr->cvatom_imprp[i1],v1);
      v_tally9(thr->cvatom_imprp[i2],v2);
      v_tally9(thr->cvatom_imprp[i3],v3);
      v_tally9(thr->cvatom_imprp[i4],v4);
    } else {
      if (i1 < nlocal) v_tally9(thr->cvatom_imprp[i1],v1);
      if (i2 < nlocal) v_tally9(thr->cvatom_imprp[i2],v2);
      if (i3 < nlocal) v_tally9(thr->cvatom_imprp[i3],v3);
      if (i4 < nlocal) v_tally9(thr->cvatom_imprp[i4],v4);
    }
  }

}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void ThrOMP::v_tally2_thr(Pair *const pair, const int i, const int j, const double fpair,
                          const double * const drij, ThrData * const thr)
{
  double v[6];

  v[0] = drij[0]*drij[0]*fpair;
  v[1] = drij[1]*drij[1]*fpair;
  v[2] = drij[2]*drij[2]*fpair;
  v[3] = drij[0]*drij[1]*fpair;
  v[4] = drij[0]*drij[2]*fpair;
  v[5] = drij[1]*drij[2]*fpair;
  if (pair->vflag_global) v_tally(thr->virial_pair,v);

  if (pair->vflag_atom) {
    v[0] *= 0.5;
    v[1] *= 0.5;
    v[2] *= 0.5;
    v[3] *= 0.5;
    v[4] *= 0.5;
    v[5] *= 0.5;
    v_tally(thr->vatom_pair[i],v);
    v_tally(thr->vatom_pair[j],v);
  }
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by RexaFF potential, newton_pair is always on
   fi is magnitude of force on atom i, deli is the direction
   note that the other atom (j) is not updated, due to newton on
------------------------------------------------------------------------- */

void ThrOMP::v_tally2_newton_thr(Pair *const pair, const int i, const double * const fi,
                                 const double * const deli, ThrData * const thr)
{
  double v[6];

  v[0] = deli[0]*fi[0];
  v[1] = deli[1]*fi[1];
  v[2] = deli[2]*fi[2];
  v[3] = deli[0]*fi[1];
  v[4] = deli[0]*fi[2];
  v[5] = deli[1]*fi[2];
  if (pair->vflag_global) v_tally(thr->virial_pair,v);
  if (pair->vflag_atom) v_tally(thr->vatom_pair[i],v);
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void ThrOMP::v_tally3_thr(Pair *const pair, const int i, const int j, const int k,
                          const double * const fi, const double * const fj,
                          const double * const drik, const double * const drjk,
                          ThrData * const thr)
{
  double v[6];

  v[0] = (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = (drik[1]*fi[2] + drjk[1]*fj[2]);
  if (pair->vflag_global) v_tally(thr->virial_pair,v);

  if (pair->vflag_atom) {
    v[0] *= THIRD;
    v[1] *= THIRD;
    v[2] *= THIRD;
    v[3] *= THIRD;
    v[4] *= THIRD;
    v[5] *= THIRD;
    v_tally(thr->vatom_pair[i],v);
    v_tally(thr->vatom_pair[j],v);
    v_tally(thr->vatom_pair[k],v);
  }
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
------------------------------------------------------------------------- */

void ThrOMP::v_tally4_thr(Pair *const pair, const int i, const int j, const int k, const int m,
                          const double * const fi, const double * const fj,
                          const double * const fk, const double * const drim,
                          const double * const drjm, const double * const drkm,
                          ThrData * const thr)
{
  double v[6];

  v[0] = (drim[0]*fi[0] + drjm[0]*fj[0] + drkm[0]*fk[0]);
  v[1] = (drim[1]*fi[1] + drjm[1]*fj[1] + drkm[1]*fk[1]);
  v[2] = (drim[2]*fi[2] + drjm[2]*fj[2] + drkm[2]*fk[2]);
  v[3] = (drim[0]*fi[1] + drjm[0]*fj[1] + drkm[0]*fk[1]);
  v[4] = (drim[0]*fi[2] + drjm[0]*fj[2] + drkm[0]*fk[2]);
  v[5] = (drim[1]*fi[2] + drjm[1]*fj[2] + drkm[1]*fk[2]);
  if (pair->vflag_global) v_tally(thr->virial_pair,v);

  if (pair->vflag_atom) {
    v[0] *= 0.25;
    v[1] *= 0.25;
    v[2] *= 0.25;
    v[3] *= 0.25;
    v[4] *= 0.25;
    v[5] *= 0.25;
    v_tally(thr->vatom_pair[i],v);
    v_tally(thr->vatom_pair[j],v);
    v_tally(thr->vatom_pair[k],v);
    v_tally(thr->vatom_pair[m],v);
  }
}

/* ---------------------------------------------------------------------- */

double ThrOMP::memory_usage_thr()
{
  double bytes=0.0;

  return bytes;
}
