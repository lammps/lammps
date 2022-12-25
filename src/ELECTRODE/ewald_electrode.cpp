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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "ewald_electrode.h"

#include "atom.h"
#include "boundary_correction.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "pair.h"
#include "slab_2d.h"
#include "slab_dipole.h"
#include "update.h"
#include "wire_dipole.h"

#include <algorithm>
#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

EwaldElectrode::EwaldElectrode(LAMMPS *lmp) : Ewald(lmp), boundcorr(nullptr)
{
  eikr_step = -1;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldElectrode::~EwaldElectrode()
{
  delete boundcorr;
}

/* ---------------------------------------------------------------------- */

void EwaldElectrode::init()
{
  if (comm->me == 0) utils::logmesg(lmp, "Ewald/electrode initialization ...\n");

  // error check
  if (domain->triclinic) error->all(FLERR, "Cannot (yet) use ewald/electrode with triclinic box ");
  // triclinic_check();
  if (domain->dimension == 2) error->all(FLERR, "Cannot use ewald/electrode with 2d simulation");

  if (!atom->q_flag) error->all(FLERR, "KSpace style ewald/electrode requires atom attribute q");

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with ewald/electrode");
  if (slabflag) {
    if (wireflag) error->all(FLERR, "Cannot use slab and wire corrections together");
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || domain->boundary[2][0] != 1 ||
        domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab ewald/electrod");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
        domain->boundary[1][0] != 1 || domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire ewald/electrode");
  }

  int *per = domain->periodicity;
  auto equal_periodicity = [per](int a[3]) {
    for (int i = 0; i < 3; i++)
      if (a[i] != per[i]) return false;
    return true;
  };
  int periodicity_2d[] = {1, 1, 0};
  int periodicity_1d[] = {0, 0, 1};
  if (boundcorr != nullptr) delete boundcorr;
  if (slabflag == 1) {
    // EW3Dc dipole correction
    if (!equal_periodicity(periodicity_2d))
      error->all(FLERR, "Two dimensional system must use p p f");
    boundcorr = new SlabDipole(lmp);
  } else if (slabflag == 3) {
    // EW2D
    if (!equal_periodicity(periodicity_2d))
      error->all(FLERR, "Two dimensional system must use p p f");
    boundcorr = new Slab2d(lmp);
  } else if (wireflag == 1) {
    // EW3Dc wire correction
    if (!equal_periodicity(periodicity_1d))
      error->all(FLERR, "One dimensional system must use f f p");
    boundcorr = new WireDipole(lmp);
  } else {
    // dummy BoundaryCorrection for ffield
    boundcorr = new BoundaryCorrection(lmp);
  }

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  triclinic = domain->triclinic;
  pair_check();

  int itmp;
  double *p_cutoff = (double *) force->pair->extract("cut_coul", itmp);
  if (p_cutoff == nullptr) error->all(FLERR, "KSpace style is incompatible with Pair style");
  double cutoff = *p_cutoff;

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0)
    accuracy = accuracy_absolute;
  else
    accuracy = accuracy_relative * two_charge_force;

  // setup K-space resolution

  bigint natoms = atom->natoms;

  // use xprd,yprd,zprd even if triclinic so grid size is the same
  // adjust z dimension for 2d slab Ewald
  // 3d Ewald just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  if (!gewaldflag) {
    if (accuracy <= 0.0) error->all(FLERR, "KSpace accuracy must be > 0");
    if (q2 == 0.0) error->all(FLERR, "Must use 'kspace_modify gewald' for uncharged system");
    g_ewald = accuracy * sqrt(natoms * cutoff * xprd * yprd * zprd) / (2.0 * q2);
    if (g_ewald >= 1.0)
      g_ewald = (1.35 - 0.15 * log(accuracy)) / cutoff;
    else
      g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // setup Ewald coefficients so can print stats

  setup();

  // final RMS accuracy

  double lprx = rms(kxmax_orig, xprd_wire, natoms, q2);
  double lpry = rms(kymax_orig, yprd_wire, natoms, q2);
  double lprz = rms(kzmax_orig, zprd_slab, natoms, q2);
  double lpr = sqrt(lprx * lprx + lpry * lpry + lprz * lprz) / sqrt(3.0);
  double q2_over_sqrt = q2 / sqrt(natoms * cutoff * xprd_wire * yprd_wire * zprd_slab);
  double spr = 2.0 * q2_over_sqrt * exp(-g_ewald * g_ewald * cutoff * cutoff);
  double tpr = estimate_table_accuracy(q2_over_sqrt, spr);
  double estimated_accuracy = sqrt(lpr * lpr + spr * spr + tpr * tpr);

  // stats

  if (comm->me == 0) {
    std::string mesg = fmt::format("  G vector (1/distance) = {:.8g}\n", g_ewald);
    mesg += fmt::format("  estimated absolute RMS force accuracy = {:.8g}\n", estimated_accuracy);
    mesg += fmt::format("  estimated relative force accuracy = {:.8g}\n",
                        estimated_accuracy / two_charge_force);
    mesg += fmt::format("  KSpace vectors: actual max1d max3d = {} {} {}\n", kcount, kmax, kmax3d);
    mesg += fmt::format("                  kxmax kymax kzmax  = {} {} {}\n", kxmax, kymax, kzmax);
    utils::logmesg(lmp, mesg);
  }
}

/* ----------------------------------------------------------------------
   adjust Ewald coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void EwaldElectrode::setup()
{
  // volume-dependent factors

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // adjustment of z dimension for 2d slab Ewald
  // 3d Ewald just uses zprd since slab_volfactor = 1.0

  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;
  volume = xprd_wire * yprd_wire * zprd_slab;

  area = xprd_wire * yprd_wire;

  unitk[0] = 2.0 * MY_PI / xprd_wire;
  unitk[1] = 2.0 * MY_PI / yprd_wire;
  unitk[2] = 2.0 * MY_PI / zprd_slab;

  int kmax_old = kmax;

  if (kewaldflag == 0) {
    // determine kmax
    // function of current box size, accuracy, G_ewald (short-range cutoff)

    bigint natoms = atom->natoms;
    double err;
    kxmax = 1;
    kymax = 1;
    kzmax = 1;

    err = rms(kxmax, xprd_wire, natoms, q2);
    while (err > accuracy) {
      kxmax++;
      err = rms(kxmax, xprd_wire, natoms, q2);
    }

    err = rms(kymax, yprd_wire, natoms, q2);
    while (err > accuracy) {
      kymax++;
      err = rms(kymax, yprd_wire, natoms, q2);
    }

    err = rms(kzmax, zprd_slab, natoms, q2);
    while (err > accuracy) {
      kzmax++;
      err = rms(kzmax, zprd_slab, natoms, q2);
    }

    kmax = MAX(kxmax, kymax);
    kmax = MAX(kmax, kzmax);
    kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;

    double gsqxmx = unitk[0] * unitk[0] * kxmax * kxmax;
    double gsqymx = unitk[1] * unitk[1] * kymax * kymax;
    double gsqzmx = unitk[2] * unitk[2] * kzmax * kzmax;
    gsqmx = MAX(gsqxmx, gsqymx);
    gsqmx = MAX(gsqmx, gsqzmx);

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

  } else {
    kxmax = kx_ewald;
    kymax = ky_ewald;
    kzmax = kz_ewald;

    kxmax_orig = kxmax;
    kymax_orig = kymax;
    kzmax_orig = kzmax;

    kmax = MAX(kxmax, kymax);
    kmax = MAX(kmax, kzmax);
    kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;

    double gsqxmx = unitk[0] * unitk[0] * kxmax * kxmax;
    double gsqymx = unitk[1] * unitk[1] * kymax * kymax;
    double gsqzmx = unitk[2] * unitk[2] * kzmax * kzmax;
    gsqmx = MAX(gsqxmx, gsqymx);
    gsqmx = MAX(gsqmx, gsqzmx);
  }

  gsqmx *= 1.00001;

  // if size has grown, reallocate k-dependent and nlocal-dependent arrays

  if (kmax > kmax_old) {
    deallocate();
    allocate();
    group_allocate_flag = 0;

    memory->destroy(ek);
    memory->destroy3d_offset(cs, -kmax_created);
    memory->destroy3d_offset(sn, -kmax_created);
    nmax = atom->nmax;
    memory->create(ek, nmax, 3, "ewald/electrode:ek");
    memory->create3d_offset(cs, -kmax, kmax, 3, nmax, "ewald/electrode:cs");
    memory->create3d_offset(sn, -kmax, kmax, 3, nmax, "ewald/electrode:sn");
    kmax_created = kmax;
  }

  // pre-compute Ewald coefficients

  coeffs();
}

/* ----------------------------------------------------------------------
   compute the Ewald long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldElectrode::compute(int eflag, int vflag)
{
  // set energy/virial flags
  ev_init(eflag, vflag);

  qsum_qsq(0);

  // return if there are no charges
  if (qsqsum == 0.0) return;

  update_eikr(true);

  MPI_Allreduce(sfacrl, sfacrl_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim, sfacim_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  // K-space portion of electric field
  // double loop over K-vectors and local atoms
  // perform per-atom calculations if needed

  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  int kx, ky, kz;
  double cypz, sypz, exprl, expim, partial, partial_peratom;

  for (int i = 0; i < nlocal; i++) {
    ek[i][0] = 0.0;
    ek[i][1] = 0.0;
    ek[i][2] = 0.0;
  }

  for (int k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];

    for (int i = 0; i < nlocal; i++) {
      cypz = cs[ky][1][i] * cs[kz][2][i] - sn[ky][1][i] * sn[kz][2][i];
      sypz = sn[ky][1][i] * cs[kz][2][i] + cs[ky][1][i] * sn[kz][2][i];
      exprl = cs[kx][0][i] * cypz - sn[kx][0][i] * sypz;
      expim = sn[kx][0][i] * cypz + cs[kx][0][i] * sypz;
      partial = expim * sfacrl_all[k] - exprl * sfacim_all[k];
      ek[i][0] += partial * eg[k][0];
      ek[i][1] += partial * eg[k][1];
      ek[i][2] += partial * eg[k][2];

      if (evflag_atom) {
        partial_peratom = exprl * sfacrl_all[k] + expim * sfacim_all[k];
        if (eflag_atom) eatom[i] += q[i] * ug[k] * partial_peratom;
        if (vflag_atom)
          for (int j = 0; j < 6; j++) vatom[i][j] += ug[k] * vg[k][j] * partial_peratom;
      }
    }
  }

  // convert E-field to force

  const double qscale = qqrd2e * scale;

  for (int i = 0; i < nlocal; i++) {
    if (wireflag != 2) {
      f[i][0] += qscale * q[i] * ek[i][0];
      f[i][1] += qscale * q[i] * ek[i][1];
    }
    if (slabflag != 2) { f[i][2] += qscale * q[i] * ek[i][2]; }
  }

  // sum global energy across KSpace vevs and add in volume-dependent term

  if (eflag_global) {
    for (int k = 0; k < kcount; k++)
      energy += ug[k] * (sfacrl_all[k] * sfacrl_all[k] + sfacim_all[k] * sfacim_all[k]);

    energy -= g_ewald * qsqsum / MY_PIS + MY_PI2 * qsum * qsum / (g_ewald * g_ewald * volume);
    energy *= qscale;
  }

  // global virial

  if (vflag_global) {
    double uk;
    for (int k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k] * sfacrl_all[k] + sfacim_all[k] * sfacim_all[k]);
      for (int j = 0; j < 6; j++) virial[j] += uk * vg[k][j];
    }
    for (int j = 0; j < 6; j++) virial[j] *= qscale;
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    if (eflag_atom) {
      for (int i = 0; i < nlocal; i++) {
        eatom[i] -=
            g_ewald * q[i] * q[i] / MY_PIS + MY_PI2 * q[i] * qsum / (g_ewald * g_ewald * volume);
        eatom[i] *= qscale;
      }
    }

    if (vflag_atom)
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < 6; j++) vatom[i][j] *= q[i] * qscale;
  }

  boundcorr->compute_corr(qsum, eflag_atom, eflag_global, energy, eatom);
}

/* ---------------------------------------------------------------------- */

void EwaldElectrode::eik_dot_r()
{
  int i, k, l, m, n, ic;
  double cstr1, sstr1, cstr2, sstr2, cstr3, sstr3, cstr4, sstr4;
  double sqk, clpm, slpm;

  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  n = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic] * unitk[ic];
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk[ic] * x[i][ic]);
        sn[1][ic][i] = sin(unitk[ic] * x[i][ic]);
        cs[-1][ic][i] = cs[1][ic][i];
        sn[-1][ic][i] = -sn[1][ic][i];
        cstr1 += q[i] * cs[1][ic][i];
        sstr1 += q[i] * sn[1][ic][i];
      }
      if (slabflag != 3 || ic < 2) {    // skip (0, 0, m) for ew2d
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
      }
    }
  }

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m * unitk[ic] * m * unitk[ic];
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cs[m][ic][i] = cs[m - 1][ic][i] * cs[1][ic][i] - sn[m - 1][ic][i] * sn[1][ic][i];
          sn[m][ic][i] = sn[m - 1][ic][i] * cs[1][ic][i] + cs[m - 1][ic][i] * sn[1][ic][i];
          cs[-m][ic][i] = cs[m][ic][i];
          sn[-m][ic][i] = -sn[m][ic][i];
          cstr1 += q[i] * cs[m][ic][i];
          sstr1 += q[i] * sn[m][ic][i];
        }
        if (slabflag != 3 || ic < 2) {    // skip (0, 0, m) for ew2d
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
        }
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (k * unitk[0] * k * unitk[0]) + (l * unitk[1] * l * unitk[1]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 += q[i] * (cs[k][0][i] * cs[l][1][i] - sn[k][0][i] * sn[l][1][i]);
          sstr1 += q[i] * (sn[k][0][i] * cs[l][1][i] + cs[k][0][i] * sn[l][1][i]);
          cstr2 += q[i] * (cs[k][0][i] * cs[l][1][i] + sn[k][0][i] * sn[l][1][i]);
          sstr2 += q[i] * (sn[k][0][i] * cs[l][1][i] - cs[k][0][i] * sn[l][1][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (l * unitk[1] * l * unitk[1]) + (m * unitk[2] * m * unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 += q[i] * (cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i]);
          sstr1 += q[i] * (sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i]);
          cstr2 += q[i] * (cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i]);
          sstr2 += q[i] * (sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (k * unitk[0] * k * unitk[0]) + (m * unitk[2] * m * unitk[2]);
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        cstr2 = 0.0;
        sstr2 = 0.0;
        for (i = 0; i < nlocal; i++) {
          cstr1 += q[i] * (cs[k][0][i] * cs[m][2][i] - sn[k][0][i] * sn[m][2][i]);
          sstr1 += q[i] * (sn[k][0][i] * cs[m][2][i] + cs[k][0][i] * sn[m][2][i]);
          cstr2 += q[i] * (cs[k][0][i] * cs[m][2][i] + sn[k][0][i] * sn[m][2][i]);
          sstr2 += q[i] * (sn[k][0][i] * cs[m][2][i] - cs[k][0][i] * sn[m][2][i]);
        }
        sfacrl[n] = cstr1;
        sfacim[n++] = sstr1;
        sfacrl[n] = cstr2;
        sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k * unitk[0] * k * unitk[0]) + (l * unitk[1] * l * unitk[1]) +
            (m * unitk[2] * m * unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          cstr3 = 0.0;
          sstr3 = 0.0;
          cstr4 = 0.0;
          sstr4 = 0.0;
          for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i];
            slpm = sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i];
            cstr1 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr1 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i];
            slpm = -sn[l][1][i] * cs[m][2][i] + cs[l][1][i] * sn[m][2][i];
            cstr2 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr2 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] + sn[l][1][i] * sn[m][2][i];
            slpm = sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i];
            cstr3 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr3 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);

            clpm = cs[l][1][i] * cs[m][2][i] - sn[l][1][i] * sn[m][2][i];
            slpm = -sn[l][1][i] * cs[m][2][i] - cs[l][1][i] * sn[m][2][i];
            cstr4 += q[i] * (cs[k][0][i] * clpm - sn[k][0][i] * slpm);
            sstr4 += q[i] * (sn[k][0][i] * clpm + cs[k][0][i] * slpm);
          }
          sfacrl[n] = cstr1;
          sfacim[n++] = sstr1;
          sfacrl[n] = cstr2;
          sfacim[n++] = sstr2;
          sfacrl[n] = cstr3;
          sfacim[n++] = sstr3;
          sfacrl[n] = cstr4;
          sfacim[n++] = sstr4;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pre-compute coefficients for each Ewald K-vector
------------------------------------------------------------------------- */

void EwaldElectrode::coeffs()
{
  int k, l, m;
  double sqk, vterm;

  double g_ewald_sq_inv = 1.0 / (g_ewald * g_ewald);
  double preu = 4.0 * MY_PI / volume;

  kcount = 0;

  // (k,0,0), (0,l,0), (0,0,m), skip (0,0) in case of EW2D (slabflag == 3)
  for (m = 1; m <= kmax; m++) {
    sqk = (m * unitk[0]) * (m * unitk[0]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = m;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = 0;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 2.0 * unitk[0] * m * ug[kcount];
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 0.0;
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0 + vterm * (unitk[0] * m) * (unitk[0] * m);
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m * unitk[1]) * (m * unitk[1]);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = m;
      kzvecs[kcount] = 0;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 2.0 * unitk[1] * m * ug[kcount];
      eg[kcount][2] = 0.0;
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0 + vterm * (unitk[1] * m) * (unitk[1] * m);
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m * unitk[2]) * (m * unitk[2]);
    if (sqk <= gsqmx && slabflag != 3) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
      vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[1] * l) * (unitk[1] * l);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 0.0;
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0;
        vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = -l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 0.0;
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0;
        vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = 0.0;
        kcount++;
        ;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[1] * l) * (unitk[1] * l) + (unitk[2] * m) * (unitk[2] * m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 0.0;
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0;
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
        kcount++;

        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = -m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 0.0;
        eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
        eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
        vg[kcount][0] = 1.0;
        vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = 0.0;
        vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
        kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[2] * m) * (unitk[2] * m);
      if (sqk <= gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 0.0;
        eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
        vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0;
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
        vg[kcount][5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = -m;
        ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
        eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
        eg[kcount][1] = 0.0;
        eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
        vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
        vg[kcount][1] = 1.0;
        vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
        vg[kcount][3] = 0.0;
        vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
        vg[kcount][5] = 0.0;
        kcount++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (unitk[0] * k) * (unitk[0] * k) + (unitk[1] * l) * (unitk[1] * l) +
            (unitk[2] * m) * (unitk[2] * m);
        if (sqk <= gsqmx) {
          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
          vterm = -2.0 * (1.0 / sqk + 0.25 * g_ewald_sq_inv);
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = 2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = 2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = -vterm * unitk[1] * l * unitk[2] * m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu * exp(-0.25 * sqk * g_ewald_sq_inv) / sqk;
          eg[kcount][0] = 2.0 * unitk[0] * k * ug[kcount];
          eg[kcount][1] = -2.0 * unitk[1] * l * ug[kcount];
          eg[kcount][2] = -2.0 * unitk[2] * m * ug[kcount];
          vg[kcount][0] = 1.0 + vterm * (unitk[0] * k) * (unitk[0] * k);
          vg[kcount][1] = 1.0 + vterm * (unitk[1] * l) * (unitk[1] * l);
          vg[kcount][2] = 1.0 + vterm * (unitk[2] * m) * (unitk[2] * m);
          vg[kcount][3] = -vterm * unitk[0] * k * unitk[1] * l;
          vg[kcount][4] = -vterm * unitk[0] * k * unitk[2] * m;
          vg[kcount][5] = vterm * unitk[1] * l * unitk[2] * m;
          kcount++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   group-group interactions
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   compute the Ewald total long-range force and energy for groups A and B
 ------------------------------------------------------------------------- */

void EwaldElectrode::compute_group_group(int /*groupbit_A*/, int /*groupbit_B*/, int /*AA_flag*/)
{
  error->all(FLERR, "Cannot (yet) use ewald/electrode style with compute group/group");
}

/* ----------------------------------------------------------------------
   compute b-vector of constant potential approach
 ------------------------------------------------------------------------- */

void EwaldElectrode::compute_vector(double *vec, int sensor_grpbit, int source_grpbit,
                                    bool invert_source)
{
  update_eikr(false);

  int const nlocal = atom->nlocal;
  double *q = atom->q;
  int *mask = atom->mask;
  std::vector<double> q_cos(kcount);
  std::vector<double> q_sin(kcount);

  for (int k = 0; k < kcount; k++) {
    int const kx = kxvecs[k];
    int const ky = kyvecs[k];
    int const kz = kzvecs[k];
    double q_cos_k = 0;
    double q_sin_k = 0;
    for (int i = 0; i < nlocal; i++) {
      bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
      if (!i_in_source) continue;    // only electrode atoms
      double const cos_kxky = cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
      double const sin_kxky = sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];
      double const cos_kr = cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
      double const sin_kr = sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];

      q_cos_k += q[i] * cos_kr;
      q_sin_k += q[i] * sin_kr;
    }
    q_cos[k] = q_cos_k;
    q_sin[k] = q_sin_k;
  }

  MPI_Allreduce(MPI_IN_PLACE, &q_cos.front(), kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &q_sin.front(), kcount, MPI_DOUBLE, MPI_SUM, world);

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & sensor_grpbit)) continue;
    double bi = 0;
    for (int k = 0; k < kcount; k++) {
      int const kx = kxvecs[k];
      int const ky = kyvecs[k];
      int const kz = kzvecs[k];
      double const cos_kxky = cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
      double const sin_kxky = sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];
      double const cos_kr = cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
      double const sin_kr = sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];
      bi += 2 * ug[k] * (cos_kr * q_cos[k] + sin_kr * q_sin[k]);
    }
    vec[i] += bi;
  }
}

/* ----------------------------------------------------------------------
   compute b-vector EW3DC correction of constant potential approach
 ------------------------------------------------------------------------- */

void EwaldElectrode::compute_vector_corr(double *vec, int sensor_grpbit, int source_grpbit,
                                         bool invert_source)
{
  update_eikr(false);
  boundcorr->vector_corr(vec, sensor_grpbit, source_grpbit, invert_source);
}

/* ----------------------------------------------------------------------
   compute individual interactions between all pairs of atoms in group A
   and B. see lammps_gather_atoms_concat() on how all sn and cs have been
   obtained.
 ------------------------------------------------------------------------- */

void EwaldElectrode::compute_matrix(bigint *imat, double **matrix, bool /* timer_flag */)
{
  update_eikr(false);
  int nlocal = atom->nlocal;
  int nprocs = comm->nprocs;

  double *csx, *csy, *csz, *snx, *sny, *snz;
  double *csx_all, *csy_all, *csz_all;
  double *snx_all, *sny_all, *snz_all;
  bigint *jmat, *jmat_local;
  // how many local group atoms owns each proc and how many in total
  bigint ngroup = 0;
  int ngrouplocal = std::count_if(&imat[0], &imat[nlocal], [](int i) {
    return i >= 0;
  });
  MPI_Allreduce(&ngrouplocal, &ngroup, 1, MPI_INT, MPI_SUM, world);

  // gather only subset of local sn and cs on each proc

  memory->create(csx, ngrouplocal * (kxmax + 1), "ewald/electrode:csx");
  memory->create(snx, ngrouplocal * (kxmax + 1), "ewald/electrode:snx");
  memory->create(csy, ngrouplocal * (kymax + 1), "ewald/electrode:csy");
  memory->create(sny, ngrouplocal * (kymax + 1), "ewald/electrode:sny");
  memory->create(snz, ngrouplocal * (kzmax + 1), "ewald/electrode:snz");
  memory->create(csz, ngrouplocal * (kzmax + 1), "ewald/electrode:csz");

  memory->create(jmat_local, ngrouplocal, "ewald/electrode:jmat_local");

  // copy subsets of local sn and cn to new local group arrays
  // beeing as memory efficient as one can possibly be ...

  for (int i = 0, n = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;

    for (int k = 0; k <= kxmax; k++) {
      csx[k + n * (kxmax + 1)] = cs[k][0][i];
      snx[k + n * (kxmax + 1)] = sn[k][0][i];
    }
    for (int k = 0; k <= kymax; k++) {
      csy[k + n * (kymax + 1)] = cs[k][1][i];
      sny[k + n * (kymax + 1)] = sn[k][1][i];
    }
    for (int k = 0; k <= kzmax; k++) {
      csz[k + n * (kzmax + 1)] = cs[k][2][i];
      snz[k + n * (kzmax + 1)] = sn[k][2][i];
    }
    jmat_local[n] = imat[i];
    n++;
  }

  if (((bigint) kxmax + 1) * ngroup > INT_MAX)
    error->all(FLERR, "kmax is too large, integer overflows might occur.");

  memory->create(csx_all, ((bigint) kxmax + 1) * ngroup, "ewald/electrode:csx_all");
  memory->create(snx_all, ((bigint) kxmax + 1) * ngroup, "ewald/electrode:snx_all");
  memory->create(csy_all, ((bigint) kymax + 1) * ngroup, "ewald/electrode:csy_all");
  memory->create(sny_all, ((bigint) kymax + 1) * ngroup, "ewald/electrode:sny_all");
  memory->create(csz_all, ((bigint) kzmax + 1) * ngroup, "ewald/electrode:csz_all");
  memory->create(snz_all, ((bigint) kzmax + 1) * ngroup, "ewald/electrode:snz_all");

  memory->create(jmat, ngroup, "ewald/electrode:jmat");

  int *recvcounts, *displs;
  memory->create(recvcounts, nprocs, "ewald/electrode:recvcounts");
  memory->create(displs, nprocs, "ewald/electrode:displs");

  // gather subsets global cs and sn
  int n = (kxmax + 1) * ngrouplocal;

  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(csx, n, MPI_DOUBLE, csx_all, recvcounts, displs, MPI_DOUBLE, world);
  MPI_Allgatherv(&snx[0], n, MPI_DOUBLE, snx_all, recvcounts, displs, MPI_DOUBLE, world);
  n = (kymax + 1) * ngrouplocal;
  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&csy[0], n, MPI_DOUBLE, csy_all, recvcounts, displs, MPI_DOUBLE, world);
  MPI_Allgatherv(&sny[0], n, MPI_DOUBLE, sny_all, recvcounts, displs, MPI_DOUBLE, world);

  n = (kzmax + 1) * ngrouplocal;
  MPI_Allgather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&csz[0], n, MPI_DOUBLE, csz_all, recvcounts, displs, MPI_DOUBLE, world);
  MPI_Allgatherv(&snz[0], n, MPI_DOUBLE, snz_all, recvcounts, displs, MPI_DOUBLE, world);

  // gather subsets global matrix indexing

  MPI_Allgather(&ngrouplocal, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
  MPI_Allgatherv(&jmat_local[0], ngrouplocal, MPI_LMP_BIGINT, jmat, recvcounts, displs,
                 MPI_LMP_BIGINT, world);

  memory->destroy(displs);
  memory->destroy(recvcounts);

  memory->destroy(jmat_local);

  // aij for each atom pair in groups; first loop over i,j then over k to
  // reduce memory access
  for (int i = 0; i < nlocal; i++) {
    if (imat[i] < 0) continue;

    for (bigint j = 0; j < ngroup; j++) {
      // matrix is symmetric, skip upper triangular matrix
      if (jmat[j] > imat[i]) continue;

      double aij = 0.0;

      for (int k = 0; k < kcount; k++) {
        // local  indexing  cs[k_idim][idim][i]       <>
        // csx_all[i+k*ngrouplocal+displs[comm->me]]]

        // anyway, use local sn and cs for simplicity

        int const kx = kxvecs[k];
        int const ky = kyvecs[k];
        int const kz = kzvecs[k];
        int const sign_ky = (ky > 0) - (ky < 0);
        int const sign_kz = (kz > 0) - (kz < 0);

        double cos_kxky = cs[kx][0][i] * cs[ky][1][i] - sn[kx][0][i] * sn[ky][1][i];
        double sin_kxky = sn[kx][0][i] * cs[ky][1][i] + cs[kx][0][i] * sn[ky][1][i];

        double const cos_kxkykz_i = cos_kxky * cs[kz][2][i] - sin_kxky * sn[kz][2][i];
        double const sin_kxkykz_i = sin_kxky * cs[kz][2][i] + cos_kxky * sn[kz][2][i];

        // global indexing  csx_all[kx+j*(kxmax+1)]  <>  csx_all[kx][j]

        int const kxj = kx + j * (kxmax + 1);
        int const kyj = abs(ky) + j * (kymax + 1);
        int const kzj = abs(kz) + j * (kzmax + 1);

        cos_kxky = csx_all[kxj] * csy_all[kyj] - snx_all[kxj] * sny_all[kyj] * sign_ky;
        sin_kxky = snx_all[kxj] * csy_all[kyj] + csx_all[kxj] * sny_all[kyj] * sign_ky;

        double const cos_kxkykz_j = cos_kxky * csz_all[kzj] - sin_kxky * snz_all[kzj] * sign_kz;
        double const sin_kxkykz_j = sin_kxky * csz_all[kzj] + cos_kxky * snz_all[kzj] * sign_kz;

        aij += 2.0 * ug[k] * (cos_kxkykz_i * cos_kxkykz_j + sin_kxkykz_i * sin_kxkykz_j);
      }
      matrix[imat[i]][jmat[j]] += aij;
      if (imat[i] != jmat[j]) matrix[jmat[j]][imat[i]] += aij;
    }
  }

  memory->destroy(jmat);
  memory->destroy(csx_all);
  memory->destroy(snx_all);
  memory->destroy(csy_all);
  memory->destroy(sny_all);
  memory->destroy(csz_all);
  memory->destroy(snz_all);
  memory->destroy(csx);
  memory->destroy(snx);
  memory->destroy(csy);
  memory->destroy(sny);
  memory->destroy(csz);
  memory->destroy(snz);
}

/* ----------------------------------------------------------------------
   compute individual corrections between all pairs of atoms in group A
   and B. see lammps_gather_atoms_concat() on how all sn and cs have been
   obtained.
 ------------------------------------------------------------------------- */

void EwaldElectrode::compute_matrix_corr(bigint *imat, double **matrix)
{
  update_eikr(false);
  boundcorr->matrix_corr(imat, matrix);
}

/* ---------------------------------------------------------------------- */

void EwaldElectrode::update_eikr(bool enforce_update)
{
  if (eikr_step < update->ntimestep || enforce_update) {
    // extend size of per-atom arrays if necessary
    if (atom->nmax > nmax) {
      memory->destroy(ek);
      memory->destroy3d_offset(cs, -kmax_created);
      memory->destroy3d_offset(sn, -kmax_created);
      nmax = atom->nmax;
      memory->create(ek, nmax, 3, "ewald/electrode:ek");
      memory->create3d_offset(cs, -kmax, kmax, 3, nmax, "ewald/electrode:cs");
      memory->create3d_offset(sn, -kmax, kmax, 3, nmax, "ewald/electrode:sn");
      kmax_created = kmax;
    }
    eikr_step = update->ntimestep;
    eik_dot_r();
  }
}
