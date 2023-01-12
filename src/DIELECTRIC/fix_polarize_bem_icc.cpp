/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
     Trung Nguyen and Monica Olvera de la Cruz (Northwestern)
   References:
     1) Tyagi, Suzen, Sega, Barbosa, Kantorovich, Holm, J. Chem. Phys. 2010,
        132, 154112
     2) Barros, Sinkovits, Luijten, J. Chem. Phys 2014, 140, 064903

   Implement a boundary-element solver for image charge computation (ICC)
   using successive overrelaxation
------------------------------------------------------------------------- */

#include "fix_polarize_bem_icc.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "msm_dielectric.h"
#include "pair_coul_cut_dielectric.h"
#include "pair_coul_long_dielectric.h"
#include "pair_lj_cut_coul_cut_dielectric.h"
#include "pair_lj_cut_coul_debye_dielectric.h"
#include "pair_lj_cut_coul_long_dielectric.h"
#include "pair_lj_cut_coul_msm_dielectric.h"
#include "pppm_dielectric.h"
#include "random_park.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

FixPolarizeBEMICC::FixPolarizeBEMICC(LAMMPS *_lmp, int narg, char **arg) : Fix(_lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR, "Illegal fix polarize/bem/icc command");

  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR, "Fix polarize requires atom style dielectric");

  // parse required arguments

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) error->all(FLERR, "Illegal fix polarize/bem/icc command");
  double tol = utils::numeric(FLERR, arg[4], false, lmp);
  tol_abs = tol_rel = tol;

  itr_max = 20;
  omega = 0.7;
  randomized = 0;
  ave_charge = 0;

  efield_pair = nullptr;
  efield_kspace = nullptr;

  comm_forward = 1;
  kspaceflag = 0;

  global_freq = 1;
  vector_flag = 1;
  size_vector = 2;
  extvector = 0;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMICC::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::init()
{
  int ncount = group->count(igroup);
  if (comm->me == 0) {
    utils::logmesg(lmp, "BEM/ICC solver for {} induced charges\n", ncount);
    utils::logmesg(lmp, " using pair style {}\n", force->pair_style);
    if (force->kspace) utils::logmesg(lmp, " using kspace style {}\n", force->kspace_style);
  }

  // initialize random induced charges with zero sum

  if (randomized) {

    int i;
    double *q = atom->q;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    auto random = new RanPark(lmp, seed_charge + comm->me);
    for (i = 0; i < 100; i++) random->uniform();
    double sum, tmp = 0;
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      q[i] = ave_charge * (random->uniform() - 0.5);
      tmp += q[i];
    }
    MPI_Allreduce(&tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    sum /= (double) ncount;

    tmp = 0;
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      q[i] -= sum;
      tmp += q[i];
    }
    MPI_Allreduce(&tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, world);

    delete random;
  }
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::setup(int /*vflag*/)
{
  // check if the pair styles in use are compatible

  if (strcmp(force->pair_style, "lj/cut/coul/long/dielectric") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulLongDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/long/dielectric/omp") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulLongDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/msm/dielectric") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulMSMDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/cut/dielectric") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulCutDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/cut/dielectric/omp") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulCutDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/debye/dielectric") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulDebyeDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/debye/dielectric/omp") == 0)
    efield_pair = (dynamic_cast<PairLJCutCoulDebyeDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "coul/long/dielectric") == 0)
    efield_pair = (dynamic_cast<PairCoulLongDielectric *>(force->pair))->efield;
  else if (strcmp(force->pair_style, "coul/cut/dielectric") == 0)
    efield_pair = (dynamic_cast<PairCoulCutDielectric *>(force->pair))->efield;
  else
    error->all(FLERR, "Pair style not compatible with fix polarize/bem/icc");

  // check if kspace is used for force computation

  if (force->kspace) {
    kspaceflag = 1;
    if (strcmp(force->kspace_style, "pppm/dielectric") == 0)
      efield_kspace = (dynamic_cast<PPPMDielectric *>(force->kspace))->efield;
    else if (strcmp(force->kspace_style, "msm/dielectric") == 0)
      efield_kspace = (dynamic_cast<MSMDielectric *>(force->kspace))->efield;
    else
      error->all(FLERR, "Kspace style not compatible with fix polarize/bem/icc");
  } else {
    if (kspaceflag == 1) {    // users specified kspace yes but there is no kspace pair style
      error->warning(FLERR, "No Kspace pair style available for fix polarize/bem/icc");
      kspaceflag = 0;
    }
  }

  // NOTE: epsilon0e2q converts (epsilon0 * efield) to the unit of (charge unit / squared distance unit)
  // efield is computed by pair and kspace styles in the unit of energy unit / charge unit / distance unit
  // for units real efield is in the unit of kcal/mol/e/A
  // converting from (F/m) (kcal/mol/e/A)  to e/A^2 (1 e = 1.6e-19 C, 1 m = 1e+10 A)
  // epsilon0e2q = 8.854187812813e-12 (C^2/N/m^2) * (4184 Nm/6.023e+23) /e/A
  //             = 8.854187812813e-12 * (4184/6.023e+23) * (1/1.6e-19)^2 e^2 / (1e+10 A) /e/A
  //             = 0.000240263377163643 e/A^2

  // for units metal efield is in the unit of eV/e/A
  // converting from (F/m) (eV/e/A)  to e/A^2 (1 V = 1 Nm/C)
  // epsilon0e2q = 8.854187812813e-12 (C^2/N/m^2) * (1 e Nm/C) /e/A
  //             = 8.854187812813e-12 * 1/1.6e-19 e^2 / (1e+10 A) /e/A
  //             = 0.00553386738300813 e/A^2

  // for units si efield is in the unit of J/C/m
  // converting from (F/m) (J/C/m) to C/m^2
  // epsilon0e2q = 8.854187812813e-12 (C^2/N/m^2) * (1 Nm/C/m)
  //             = 8.854187812813e-12 C/m^2

  // for units nano efield is in the unit of attogram nm^2/ns^2/e/nm
  // converting from (F/m) (attogram nm^2/ns^2/e/nm) to e/nm^2
  // epsilon0e2q = 8.854187812813e-12 (C^2/N/m^2) * (1e-21 kg nm^2 / (1e-18s^2) / e / nm)
  //             = 8.854187812813e-12 (C^2/N/m^2) * (1e-21 kg 1e-9 m / (1e-18s^2) / e)
  //             = 8.854187812813e-12 (1/1.6e-19)^2 (1e-21 * 1e-9 / (1e-18)) e / (1e+18 nm^2)
  //             = 0.000345866711328125 e/nm^2

  epsilon0e2q = 1.0;
  if (strcmp(update->unit_style, "real") == 0)
    epsilon0e2q = 0.000240263377163643;
  else if (strcmp(update->unit_style, "metal") == 0)
    epsilon0e2q = 0.00553386738300813;
  else if (strcmp(update->unit_style, "si") == 0)
    epsilon0e2q = 8.854187812813e-12;
  else if (strcmp(update->unit_style, "nano") == 0)
    epsilon0e2q = 0.000345866711328125;
  else if (strcmp(update->unit_style, "lj") != 0)
    error->all(FLERR, "Only unit styles 'lj', 'real', 'metal', 'si' and 'nano' are supported");

  compute_induced_charges();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::pre_force(int)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  compute_induced_charges();

  // make sure forces are reset to zero before actual forces are computed

  force_clear();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::compute_induced_charges()
{
  double *q_scaled = atom->q_scaled;
  double *q = atom->q;
  double **norm = atom->mu;
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *epsilon = atom->epsilon;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int eflag = 1;
  int vflag = 0;
  int itr;

  // use Eq. (64) in Barros et al. to initialize the induced charges
  // Note: area[i] is included here to ensure correct charge unit
  //   for direct use in force/efield compute
  // for induced charges q_real stores the free surface charge
  //   q_real are read from the data file
  // Note that the electrical fields here are due to the rescaled real charges,
  //   and also multiplied by epsilon[i]
  // for the interface particles assume that epsilon[i] = em[i]
  // NOTE: 13Dec2022 pair and kspace with dielectric suffix operate on q_scaled

  force_clear();
  force->pair->compute(eflag, vflag);
  if (kspaceflag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double Ex = efield_pair[i][0];
    double Ey = efield_pair[i][1];
    double Ez = efield_pair[i][2];
    if (kspaceflag) {
      Ex += efield_kspace[i][0];
      Ey += efield_kspace[i][1];
      Ez += efield_kspace[i][2];
    }

    // divide (Ex,Ey,Ez) by epsilon[i] here
    double ndotE = epsilon0e2q * (Ex * norm[i][0] + Ey * norm[i][1] + Ez * norm[i][2]) /
        epsilon[i] / (2 * MY_PI);
    double q_free = q[i];
    double q_bound = 0;
    q_bound = (1.0 / em[i] - 1) * q_free - (ed[i] / (2 * em[i])) * ndotE * area[i];
    q_scaled[i] = q_free + q_bound;
  }

  // communicate q_scaled between neighboring procs

  comm->forward_comm(this);

  // iterate until convergence

  for (itr = 0; itr < itr_max; itr++) {

    force_clear();
    force->pair->compute(eflag, vflag);
    if (kspaceflag) force->kspace->compute(eflag, vflag);
    if (force->newton) comm->reverse_comm();

    double tol = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;

      double q_free = q[i];
      double qtmp = q_scaled[i] - q_free;
      double Ex = efield_pair[i][0];
      double Ey = efield_pair[i][1];
      double Ez = efield_pair[i][2];
      if (kspaceflag) {
        Ex += efield_kspace[i][0];
        Ey += efield_kspace[i][1];
        Ez += efield_kspace[i][2];
      }

      // Eq. (69) in Barros et al.,  sigma_f[i] = q_real[i] / area[i]
      // note the area[i] is included here to ensure correct charge unit
      // for direct use in force/efield compute

      double ndotE = epsilon0e2q * (Ex * norm[i][0] + Ey * norm[i][1] + Ez * norm[i][2]) /
          (4 * MY_PI) / epsilon[i];
      double q_bound = q_scaled[i] - q_free;
      q_bound = (1 - omega) * q_bound +
          omega * ((1.0 / em[i] - 1) * q_free - (ed[i] / em[i]) * ndotE * area[i]);
      q_scaled[i] = q_free + q_bound;

      // Eq. (11) in Tyagi et al., with f from Eq. (6)
      // NOTE: Tyagi et al. defined the normal vector n_i pointing
      // from the medium containg the ions toward the other medium,
      // which makes the normal vector direction depend on the ion position
      // Also, since Tyagi et al. chose epsilon_1 for the uniform dielectric constant
      // of the equivalent system, there is epsilon_1 in f in Eq. (6).
      // Here we are using (q/epsilon_i) for the real charges and 1.0 for the equivalent system
      // hence there's no epsilon_1 in the factor f

      //double dot = (Ex*norm[i][0] + Ey*norm[i][1] + Ez*norm[i][2]);
      //double f = (ed[i] / (2 * em[i])) / (2*MY_PI);
      //q[i] = (1 - omega) * q[i] - omega * epsilon0 * f * dot * area[i];

      double delta = fabs(qtmp - q_bound);
      double r = (fabs(qtmp) > 0) ? delta / fabs(qtmp) : 0;
      if (tol < r) tol = r;
    }

    // communicate q_scaled for efield compute in the next iteration

    comm->forward_comm(this);

    MPI_Allreduce(&tol, &rho, 1, MPI_DOUBLE, MPI_MAX, world);
    if (itr > 0 && rho < tol_rel) break;
  }

  iterations = itr;

  // compute the total induced charges of the interface particles
  // for interface particles: set the charge to be the sum of unscaled (free) charges and induced charges

  double tmp = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double q_bound = q_scaled[i] - q[i];
    tmp += q_bound;
    q[i] = q_scaled[i];
  }

  // ensure sum of all induced charges being zero

  int ncount = group->count(igroup);
  double sum = 0;
  MPI_Allreduce(&tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
  double qboundave = sum/(double)ncount;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    q[i] -=  qboundave;
  }
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::force_clear()
{
  int nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&atom->f[0][0], 0, 3 * nbytes);
    if (torqueflag) memset(&atom->torque[0][0], 0, 3 * nbytes);
    if (extraflag) atom->avec->force_clear(0, nbytes);
  }
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMICC::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "itr_max") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      itr_max = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "omega") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      omega = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kspace") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      kspaceflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dielectrics") == 0) {
      if (iarg + 6 > narg) error->all(FLERR, "Illegal fix_modify command");
      double epsiloni = -1, areai = -1;
      double qunscaledi = 0;
      int set_charge = 0;
      double ediff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      double emean = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strcmp(arg[iarg + 3], "NULL") != 0)
        epsiloni = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strcmp(arg[iarg + 4], "NULL") != 0)
        areai = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strcmp(arg[iarg + 5], "NULL") != 0) {
        qunscaledi = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
        set_charge = 1;
      }
      set_dielectric_params(ediff, emean, epsiloni, areai, set_charge, qunscaledi);

      iarg += 6;
    } else if (strcmp(arg[iarg], "rand") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix_modify command");
      ave_charge = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      seed_charge = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      randomized = 1;
      iarg += 3;
    } else
      error->all(FLERR, "Illegal fix_modify command");
  }

  return iarg;
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMICC::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                         int * /*pbc*/)
{
  int m;
  for (m = 0; m < n; m++) buf[m] = atom->q_scaled[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMICC::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) atom->q_scaled[i] = buf[m];
}


/* ----------------------------------------------------------------------
   set dielectric params for the atoms in the group
------------------------------------------------------------------------- */

void FixPolarizeBEMICC::set_dielectric_params(double ediff, double emean, double epsiloni,
                                              double areai, int set_charge, double qvalue)
{
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *q = atom->q;
  double *epsilon = atom->epsilon;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ed[i] = ediff;
      em[i] = emean;
      if (areai > 0) area[i] = areai;
      if (epsiloni > 0) epsilon[i] = epsiloni;
      if (set_charge) q[i] = qvalue;
    }
  }
}

/* ----------------------------------------------------------------------
   return the actual number of iterations
      and current relative error
------------------------------------------------------------------------- */

double FixPolarizeBEMICC::compute_vector(int n)
{
  if (n == 0)
    return iterations;
  else if (n == 1)
    return rho;
  else
    return 0;
}
