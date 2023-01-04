// clang-format off
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

#include "kspace.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

KSpace::KSpace(LAMMPS *lmp) : Pointers(lmp)
{
  order_allocated = 0;
  energy = 0.0;
  virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;

  triclinic_support = 1;
  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag =
    dipoleflag = spinflag = 0;
  compute_flag = 1;
  group_group_enable = 0;
  stagger_flag = 0;

  order = 5;
  gridflag = 0;
  gewaldflag = 0;
  minorder = 2;
  overlap_allowed = 1;
  fftbench = 0;

  // default to using MPI collectives for FFT/remap only on IBM BlueGene

#ifdef __bg__
  collective_flag = 1;
#else
  collective_flag = 0;
#endif

  kewaldflag = 0;

  order_6 = 5;
  gridflag_6 = 0;
  gewaldflag_6 = 0;
  auto_disp_flag = 0;

  conp_one_step = true;
  slabflag = wireflag = 0;
  differentiation_flag = 0;
  slab_volfactor = 1;
  wire_volfactor = 1;
  suffix_flag = Suffix::NONE;
  adjust_cutoff_flag = 1;
  scalar_pressure_flag = 0;
  warn_nonneutral = 1;
  warn_nocharge = 1;

  accuracy_absolute = -1.0;
  accuracy_real_6 = -1.0;
  accuracy_kspace_6 = -1.0;

  neighrequest_flag = 1;
  mixflag = 0;

  splittol = 1.0e-6;

  maxeatom = maxvatom = 0;
  eatom = nullptr;
  vatom = nullptr;
  centroidstressflag = CENTROID_NOTAVAIL;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  copymode = 0;

  memory->create(gcons,7,7,"kspace:gcons");
  gcons[2][0] = 15.0 / 8.0;
  gcons[2][1] = -5.0 / 4.0;
  gcons[2][2] = 3.0 / 8.0;
  gcons[3][0] = 35.0 / 16.0;
  gcons[3][1] = -35.0 / 16.0;
  gcons[3][2] = 21.0 / 16.0;
  gcons[3][3] = -5.0 / 16.0;
  gcons[4][0] = 315.0 / 128.0;
  gcons[4][1] = -105.0 / 32.0;
  gcons[4][2] = 189.0 / 64.0;
  gcons[4][3] = -45.0 / 32.0;
  gcons[4][4] = 35.0 / 128.0;
  gcons[5][0] = 693.0 / 256.0;
  gcons[5][1] = -1155.0 / 256.0;
  gcons[5][2] = 693.0 / 128.0;
  gcons[5][3] = -495.0 / 128.0;
  gcons[5][4] = 385.0 / 256.0;
  gcons[5][5] = -63.0 / 256.0;
  gcons[6][0] = 3003.0 / 1024.0;
  gcons[6][1] = -3003.0 / 512.0;
  gcons[6][2] = 9009.0 / 1024.0;
  gcons[6][3] = -2145.0 / 256.0;
  gcons[6][4] = 5005.0 / 1024.0;
  gcons[6][5] = -819.0 / 512.0;
  gcons[6][6] = 231.0 / 1024.0;

  memory->create(dgcons,7,6,"kspace:dgcons");
  dgcons[2][0] = -5.0 / 2.0;
  dgcons[2][1] = 3.0 / 2.0;
  dgcons[3][0] = -35.0 / 8.0;
  dgcons[3][1] = 21.0 / 4.0;
  dgcons[3][2] = -15.0 / 8.0;
  dgcons[4][0] = -105.0 / 16.0;
  dgcons[4][1] = 189.0 / 16.0;
  dgcons[4][2] = -135.0 / 16.0;
  dgcons[4][3] = 35.0 / 16.0;
  dgcons[5][0] = -1155.0 / 128.0;
  dgcons[5][1] = 693.0 / 32.0;
  dgcons[5][2] = -1485.0 / 64.0;
  dgcons[5][3] = 385.0 / 32.0;
  dgcons[5][4] = -315.0 / 128.0;
  dgcons[6][0] = -3003.0 / 256.0;
  dgcons[6][1] = 9009.0 / 256.0;
  dgcons[6][2] = -6435.0 / 128.0;
  dgcons[6][3] = 5005.0 / 128.0;
  dgcons[6][4] = -4095.0 / 256.0;
  dgcons[6][5] = 693.0 / 256.0;
}

/* ---------------------------------------------------------------------- */

KSpace::~KSpace()
{
  if (copymode) return;

  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(gcons);
  memory->destroy(dgcons);
}

/* ----------------------------------------------------------------------
   calculate this in init() so that units are finalized
------------------------------------------------------------------------- */

void KSpace::two_charge()
{
  two_charge_force = force->qqr2e *
    (force->qelectron * force->qelectron) /
    (force->angstrom * force->angstrom);
}

/* ---------------------------------------------------------------------- */

void KSpace::triclinic_check()
{
  if (domain->triclinic && triclinic_support != 1)
    error->all(FLERR,"KSpace style does not yet support triclinic geometries");
}

/* ---------------------------------------------------------------------- */

void KSpace::compute_dummy(int eflag, int vflag)
{
  ev_init(eflag,vflag);
}

/* ----------------------------------------------------------------------
   check that pair style is compatible with long-range solver
------------------------------------------------------------------------- */

void KSpace::pair_check()
{
  if (force->pair == nullptr)
    error->all(FLERR,"KSpace solver requires a pair style");

  if (ewaldflag && !force->pair->ewaldflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (pppmflag && !force->pair->pppmflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (msmflag && !force->pair->msmflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (dispersionflag && !force->pair->dispersionflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (dipoleflag && !force->pair->dipoleflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (spinflag && !force->pair->spinflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (tip4pflag && !force->pair->tip4pflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");

  if (force->pair->dispersionflag && !dispersionflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (force->pair->tip4pflag && !tip4pflag)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for bitwise settings of eflag/vflag
   set the following flags, values are otherwise set to 0:
     evflag       != 0 if any bits of eflag or vflag are set
     eflag_global != 0 if ENERGY_GLOBAL bit of eflag set
     eflag_atom   != 0 if ENERGY_ATOM bit of eflag set
     eflag_either != 0 if eflag_global or eflag_atom is set
     vflag_global != 0 if VIRIAL_PAIR or VIRIAL_FDOTR bit of vflag set
     vflag_atom   != 0 if VIRIAL_ATOM bit of vflag set
                       no current support for centroid stress
     vflag_either != 0 if vflag_global or vflag_atom is set
     evflag_atom  != 0 if eflag_atom or vflag_atom is set
------------------------------------------------------------------------- */

void KSpace::ev_setup(int eflag, int vflag, int alloc)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_either = vflag;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  vflag_atom = vflag & VIRIAL_ATOM;

  if (eflag_atom || vflag_atom) evflag_atom = 1;
  else evflag_atom = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom,maxeatom,"kspace:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom,maxvatom,6,"kspace:vatom");
    }
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (tip4pflag) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (tip4pflag) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   compute qsum,qsqsum,q2 and give error/warning if not charge neutral
   called initially, when particle count changes, when charges are changed
------------------------------------------------------------------------- */

void KSpace::qsum_qsq(int warning_flag)
{
  const double * const q = atom->q;
  const int nlocal = atom->nlocal;
  double qsum_local(0.0), qsqsum_local(0.0);

#if defined(_OPENMP)
#pragma omp parallel for default(shared) reduction(+:qsum_local,qsqsum_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    qsum_local += q[i];
    qsqsum_local += q[i]*q[i];
  }

  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsum_local,&qsqsum,1,MPI_DOUBLE,MPI_SUM,world);

  if ((qsqsum == 0.0) && (comm->me == 0) && warn_nocharge && warning_flag) {
    error->warning(FLERR,"Using kspace solver on system with no charge");
    warn_nocharge = 0;
  }

  q2 = qsqsum * force->qqrd2e;

  // not yet sure of the correction needed for non-neutral systems
  // so issue warning or error

  if (fabs(qsum) > SMALL) {
    std::string message = fmt::format("System is not charge neutral, net "
                                      "charge = {:.8}",qsum);
    if (!warn_nonneutral) error->all(FLERR,message);
    if (warn_nonneutral == 1 && comm->me == 0) error->warning(FLERR,message);
    warn_nonneutral = 2;
  }
}

/* ----------------------------------------------------------------------
   estimate the accuracy of the short-range coulomb tables
------------------------------------------------------------------------- */

double KSpace::estimate_table_accuracy(double q2_over_sqrt, double spr)
{
  double table_accuracy = 0.0;
  int nctb = force->pair->ncoultablebits;
  if (comm->me == 0) {
    if (nctb)
      error->message(FLERR,"  using {}-bit tables for long-range coulomb",nctb);
    else
      error->message(FLERR,"  using polynomial approximation for long-range coulomb");
  }

  if (nctb) {
    double empirical_precision[17];
    empirical_precision[6] =  6.99E-03;
    empirical_precision[7] =  1.78E-03;
    empirical_precision[8] =  4.72E-04;
    empirical_precision[9] =  1.17E-04;
    empirical_precision[10] = 2.95E-05;
    empirical_precision[11] = 7.41E-06;
    empirical_precision[12] = 1.76E-06;
    empirical_precision[13] = 9.28E-07;
    empirical_precision[14] = 7.46E-07;
    empirical_precision[15] = 7.32E-07;
    empirical_precision[16] = 7.30E-07;
    if (nctb <= 6) table_accuracy = empirical_precision[6];
    else if (nctb <= 16) table_accuracy = empirical_precision[nctb];
    else table_accuracy = empirical_precision[16];
    table_accuracy *= q2_over_sqrt;
    if ((table_accuracy > spr) && (comm->me == 0))
      error->warning(FLERR,"For better accuracy use 'pair_modify table 0'");
  }

  return table_accuracy;
}

/* ----------------------------------------------------------------------
   convert box coords vector to transposed triclinic lamda (0-1) coords
   vector, lamda = [(H^-1)^T] * v, does not preserve vector magnitude
   v and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void KSpace::x2lamdaT(double *v, double *lamda)
{
  double *h_inv = domain->h_inv;
  double lamda_tmp[3];

  lamda_tmp[0] = h_inv[0]*v[0];
  lamda_tmp[1] = h_inv[5]*v[0] + h_inv[1]*v[1];
  lamda_tmp[2] = h_inv[4]*v[0] + h_inv[3]*v[1] + h_inv[2]*v[2];

  lamda[0] = lamda_tmp[0];
  lamda[1] = lamda_tmp[1];
  lamda[2] = lamda_tmp[2];
}

/* ----------------------------------------------------------------------
   convert lamda (0-1) coords vector to transposed box coords vector
   lamda = (H^T) * v, does not preserve vector magnitude
   v and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void KSpace::lamda2xT(double *lamda, double *v)
{
  double h[6];
  h[0] = domain->h[0];
  h[1] = domain->h[1];
  h[2] = domain->h[2];
  h[3] = fabs(domain->h[3]);
  h[4] = fabs(domain->h[4]);
  h[5] = fabs(domain->h[5]);
  double v_tmp[3];

  v_tmp[0] = h[0]*lamda[0];
  v_tmp[1] = h[5]*lamda[0] + h[1]*lamda[1];
  v_tmp[2] = h[4]*lamda[0] + h[3]*lamda[1] + h[2]*lamda[2];

  v[0] = v_tmp[0];
  v[1] = v_tmp[1];
  v[2] = v_tmp[2];
}

/* ----------------------------------------------------------------------
   convert triclinic lamda (0-1) coords vector to box coords vector
   v = H * lamda, does not preserve vector magnitude
   lamda and v can point to same 3-vector
------------------------------------------------------------------------- */

void KSpace::lamda2xvector(double *lamda, double *v)
{
  double *h = domain->h;

  v[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2];
  v[1] = h[1]*lamda[1] + h[3]*lamda[2];
  v[2] = h[2]*lamda[2];
}

/* ----------------------------------------------------------------------
   modify parameters of the KSpace style
------------------------------------------------------------------------- */

void KSpace::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mesh") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      nx_pppm = nx_msm_max = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      ny_pppm = ny_msm_max = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      nz_pppm = nz_msm_max = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      if (nx_pppm == 0 && ny_pppm == 0 && nz_pppm == 0)
        gridflag = 0;
      else if (nx_pppm <= 0 || ny_pppm <= 0 || nz_pppm <= 0)
        error->all(FLERR,"Kspace_modify mesh parameters must be all "
                   "zero or all positive");
      else gridflag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"mesh/disp") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      nx_pppm_6 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      ny_pppm_6 = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      nz_pppm_6 = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      if (nx_pppm_6 == 0 && ny_pppm_6 == 0 && nz_pppm_6 == 0)
        gridflag_6 = 0;
      else if (nx_pppm_6 <= 0 || ny_pppm_6 <= 0 || nz_pppm_6 == 0)
        error->all(FLERR,"Kspace_modify mesh/disp parameters must be all zero or all positive");
      else gridflag_6 = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"order") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      order = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"order/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      order_6 = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"minorder") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      minorder = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (minorder < 2) error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"overlap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      overlap_allowed = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      accuracy_absolute = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      g_ewald = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (g_ewald == 0.0) gewaldflag = 0;
      else gewaldflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      g_ewald_6 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (g_ewald_6 == 0.0) gewaldflag_6 = 0;
      else gewaldflag_6 = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"slab") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"nozforce") == 0) {
        slabflag = 2;
      } else if (strcmp(arg[iarg+1],"ew2d") == 0) {
        slabflag = 3;
      } else {
        slabflag = 1;
        slab_volfactor = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (slab_volfactor <= 1.0)
          error->all(FLERR,"Bad kspace_modify slab parameter");
        if (slab_volfactor < 2.0 && comm->me == 0)
          error->warning(FLERR,"Kspace_modify slab param < 2.0 may "
                         "cause unphysical behavior");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"wire") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"noxyforce") == 0) {
        wireflag = 2;
      } else {
        wireflag = 1;
        wire_volfactor = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (wire_volfactor <= 1.0)
          error->all(FLERR,"Bad kspace_modify slab parameter");
        if (wire_volfactor < 2.0 && comm->me == 0)
          error->warning(FLERR,"Kspace_modify slab param < 2.0 may "
                         "cause unphysical behavior");
      }
      warn_nonneutral = 0; // can't use wire correction with non-neutral system
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "amat") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal kspace_modify command");
      if (!pppmflag) error->all(FLERR, "Illegal kspace_modify command 'amat'"
                                      "available for pppm/conp, only");
      if (strcmp(arg[iarg + 1], "twostep") == 0) {
        conp_one_step = false;
      } else if (strcmp(arg[iarg + 1], "onestep") == 0) {
        conp_one_step = true;
      } else {
        error->all(FLERR, "Illegal kspace_modify command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"compute") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      compute_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"fftbench") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      fftbench = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"collective") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      collective_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"diff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"ad") == 0) differentiation_flag = 1;
      else if (strcmp(arg[iarg+1],"ik") == 0) differentiation_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff/adjust") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      adjust_cutoff_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"kmax/ewald") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      kx_ewald = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      ky_ewald = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      kz_ewald = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      if (kx_ewald < 0 || ky_ewald < 0 || kz_ewald < 0)
        error->all(FLERR,"Bad kspace_modify kmax/ewald parameter");
      if (kx_ewald > 0 && ky_ewald > 0 && kz_ewald > 0)
        kewaldflag = 1;
      else
        kewaldflag = 0;
      iarg += 4;
    } else if (strcmp(arg[iarg],"mix/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"pair") == 0) mixflag = 0;
      else if (strcmp(arg[iarg+1],"geom") == 0) mixflag = 1;
      else if (strcmp(arg[iarg+1],"none") == 0) mixflag = 2;
      else error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"force/disp/real") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      accuracy_real_6 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"force/disp/kspace") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      accuracy_kspace_6 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"eigtol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      splittol = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (splittol >= 1.0)
        error->all(FLERR,"Kspace_modify eigtol must be smaller than one");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pressure/scalar") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      scalar_pressure_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"disp/auto") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      auto_disp_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal kspace_modify command");
      iarg += n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void *KSpace::extract(const char *str)
{
  if (strcmp(str,"scale") == 0) return (void *) &scale;
  return nullptr;
}
