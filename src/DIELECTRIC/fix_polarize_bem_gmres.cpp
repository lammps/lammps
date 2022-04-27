/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
     Trung Nguyen and Monica Olvera de la Cruz (Northwestern)

   Implement a boundary element solver for surface induced charges
   using the GMRES algorithm

   Reference: Barros, Sinkovits, Luijten, J. Chem. Phys 2014, 140, 064903

   GMRES solver: Original C version by Lili Ju, C++ version by John Burkardt.
   The version adapted here is implemented where A is an operator acting on x,
   and the residual (b-Ax) is computed directly.

   References:
     1) R. Barrett, M. Berry, T. Chan, J. Demmel, J. Donato, J. Dongarra,
     V. Eijkhout, R. Pozo, C. Romine, H. van der Vorst,
     Templates for the Solution of Linear Systems: Building Blocks for
     Iterative Methods, SIAM, 1994, ISBN: 0898714710, LC: QA297.8.T45.

     2) T. Kelley, Iterative Methods for Linear and Nonlinear Equations,
     SIAM, 2004, ISBN: 0898713528, LC: QA297.8.K45.
------------------------------------------------------------------------- */

#include "fix_polarize_bem_gmres.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
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

FixPolarizeBEMGMRES::FixPolarizeBEMGMRES(LAMMPS *_lmp, int narg, char **arg) :
    Fix(_lmp, narg, arg), q_backup(nullptr), c(nullptr), g(nullptr), h(nullptr), r(nullptr),
    s(nullptr), v(nullptr), y(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal fix polarize/bem/gmres command");

  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR, "Fix polarize requires atom style dielectric");

  // parse required arguments

  nevery = utils::numeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) error->all(FLERR, "Illegal fix polarize/bem/gmres command");
  double tol = utils::numeric(FLERR, arg[4], false, lmp);
  tol_abs = tol_rel = tol;

  itr_max = 20;
  mr = 0;
  randomized = 0;
  ave_charge = 0;

  efield_pair = nullptr;
  efield_kspace = nullptr;

  comm_forward = 1;
  nmax = 0;
  allocated = 0;
  kspaceflag = 0;

  induced_charge_idx = nullptr;
  induced_charges = nullptr;
  rhs = nullptr;
  buffer = nullptr;
  tag2mat = nullptr;
  mat2tag = nullptr;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;

  FixPolarizeBEMGMRES::grow_arrays(atom->nmax);
  atom->add_callback(0);    // to ensure to work with atom->sort()

  // output the residual and actual number of iterations

  global_freq = 1;
  vector_flag = 1;
  size_vector = 2;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

FixPolarizeBEMGMRES::~FixPolarizeBEMGMRES()
{
  memory->destroy(q_backup);
  memory->destroy(induced_charge_idx);
  memory->destroy(induced_charges);
  memory->destroy(rhs);
  memory->destroy(buffer);
  memory->destroy(mat2tag);
  memory->destroy(tag2mat);

  if (allocated) FixPolarizeBEMGMRES::deallocate();
  atom->delete_callback(id, 0);
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMGMRES::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::init()
{
  // mapping induced charge matrix/vector to atom tags and vice versa

  int i, maxtag;
  double *q = atom->q;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint max_tag = -1;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) max_tag = MAX(max_tag, tag[i]);

  tagint itmp;
  MPI_Allreduce(&max_tag, &itmp, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  maxtag = (int) itmp;

  int *ncount;
  memory->create(ncount, maxtag + 1, "polarize:ncount");
  for (i = 0; i <= maxtag; i++) ncount[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) ncount[tag[i]]++;

  memory->create(tag2mat, maxtag + 1, "polarize:tag2mat");
  MPI_Allreduce(ncount, tag2mat, maxtag + 1, MPI_INT, MPI_SUM, world);

  num_induced_charges = 0;
  for (i = 0; i <= maxtag; i++)
    if (tag2mat[i])
      tag2mat[i] = num_induced_charges++;
    else
      tag2mat[i] = -1;

  memory->create(mat2tag, num_induced_charges, "polarize:mat2tag");

  num_induced_charges = 0;
  for (i = 0; i <= maxtag; i++)
    if (tag2mat[i] >= 0) mat2tag[num_induced_charges++] = i;

  for (i = 0; i < nlocal; i++) {
    induced_charge_idx[i] = -1;
    if (mask[i] & groupbit) induced_charge_idx[i] = tag2mat[tag[i]];
  }

  memory->destroy(ncount);

  // allocate memory for the solver

  memory->create(induced_charges, num_induced_charges, "polarize:induced_charges");
  memory->create(rhs, num_induced_charges, "polarize:rhs");
  memory->create(buffer, num_induced_charges, "polarize:buffer");

  mat_dim = num_induced_charges;
  if (mr > mat_dim - 1 || mr <= 0) mr = mat_dim - 1;

  if (allocated == 0) {
    allocate();
    allocated = 1;
  }

  // initialize random induced charges with zero sum

  if (randomized) {

    auto random = new RanPark(lmp, seed_charge + comm->me);
    for (i = 0; i < 100; i++) random->uniform();
    double sum, tmp = 0;
    for (i = 0; i < nlocal; i++) {
      if (induced_charge_idx[i] < 0) continue;
      q[i] = ave_charge * (random->uniform() - 0.5);
      tmp += q[i];
    }
    MPI_Allreduce(&tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    sum /= (double) num_induced_charges;

    tmp = 0;
    for (i = 0; i < nlocal; i++) {
      if (induced_charge_idx[i] < 0) continue;
      q[i] -= sum;
      tmp += q[i];
    }
    MPI_Allreduce(&tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, world);

    if (comm->me == 0) utils::logmesg(lmp, "ave induced charge q = {:.8}\n", sum);
    delete random;
  }

  if (comm->me == 0)
    utils::logmesg(lmp, "GMRES solver for {} induced charges using maximum {} q-vectors\n",
                   num_induced_charges, mr);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::setup(int /*vflag*/)
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
    error->all(FLERR, "Pair style not compatible with fix polarize/bem/gmres");

  if (kspaceflag) {
    if (force->kspace) {
      if (strcmp(force->kspace_style, "pppm/dielectric") == 0)
        efield_kspace = (dynamic_cast<PPPMDielectric *>(force->kspace))->efield;
      else if (strcmp(force->kspace_style, "msm/dielectric") == 0)
        efield_kspace = (dynamic_cast<MSMDielectric *>(force->kspace))->efield;
      else
        error->all(FLERR, "Kspace style not compatible with fix polarize/bem/gmres");
    } else
      error->all(FLERR, "No Kspace style available for fix polarize/bem/gmres");
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

  first = 1;
  compute_induced_charges();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::pre_force(int)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  compute_induced_charges();

  // make sure forces are reset to zero before actual forces are computed

  force_clear();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::compute_induced_charges()
{
  double *q = atom->q;
  double *q_real = atom->q_unscaled;
  double **norm = atom->mu;
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *epsilon = atom->epsilon;
  int nlocal = atom->nlocal;
  int eflag = 0;
  int vflag = 0;

  // compute the right hand side (vector b) of Eq. (40) according to Eq. (42)
  // keep the scaled real charges intact here to compute efield for the right hand side (b)
  //   and backup all the charges
  // for induced charges q_real stores the free surface charge
  // set the induced charges to be zero to compute the right hand side (b)
  // the current value can be accessed via induced_charges[induced_charge_idx[i]]

  for (int i = 0; i < nlocal; i++) {
    q_backup[i] = q[i];
    if (induced_charge_idx[i] >= 0) q[i] = 0;
  }

  comm->forward_comm(this);

  // note here q[i] are the bound charges including area
  // so that kspace solver can be used directly with the charge values
  // for the moment, require that newton off and full neighbor list for pair
  // Note that in the definition of the electrical fields in Equations (41) and (53)
  //   in Ref. Barros et al there is a factor 1/(4pi), and
  //   that these are the electrical field is due to the rescaled real charges
  // Note: the right-hand side (b) is in the unit of charge density

  force_clear();
  force->pair->compute(eflag, vflag);
  if (kspaceflag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();

  for (int i = 0; i < num_induced_charges; i++) buffer[i] = 0;

  for (int i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) continue;
    int idx = induced_charge_idx[i];
    if (ed[i] == 0) {
      buffer[idx] = 0;
      continue;
    }
    double Ex = efield_pair[i][0];
    double Ey = efield_pair[i][1];
    double Ez = efield_pair[i][2];
    if (kspaceflag) {
      Ex += efield_kspace[i][0];
      Ey += efield_kspace[i][1];
      Ez += efield_kspace[i][2];
    }
    double ndotE = epsilon0e2q * (Ex * norm[i][0] + Ey * norm[i][1] + Ez * norm[i][2]) / epsilon[i];
    double sigma_f = q_real[i] / area[i];
    buffer[idx] = (1 - em[i]) * sigma_f - ed[i] * ndotE / (4 * MY_PI);
  }

  MPI_Allreduce(buffer, rhs, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);

  // compute the initial residual r before iteration
  // while it seems that assigning induced charges to the last values
  //   could improve convergence, it's not necessarily the cases
  //   where the system is evolving or diffusive
  // to be defensive here, reset induced charges to zeros
  //   and initial residual r equal to rhs

  for (int i = 0; i < num_induced_charges; i++) {
    induced_charges[i] = 0;
    r[i] = rhs[i];
  }

  // get the norm of the right hand side vector

  normb = sqrt(vec_dot(rhs, rhs, num_induced_charges));
  if (normb < tol_abs) return;

  // use the GMRES solver to solve for the induced charges

  gmres_solve(induced_charges, r);

  // set the particle charges in the group to be the induced charges
  // restore the charges of the real particles (that are not in the group)

  for (int i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] >= 0) {
      int idx = induced_charge_idx[i];
      q[i] = induced_charges[idx] * area[i] + q_real[i];
    } else {
      q[i] = q_backup[i];
    }
  }

  comm->forward_comm(this);

  if (first) first = 0;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::gmres_solve(double *x, double *r)
{
  int i, j, k, k_copy, n, itr;
  double av, htmp, mu, rho_tol;
  double delta = 1.0e-03;

  n = mat_dim;

  // compute the relative tolerance
  // rho = norm(r)
  rho = sqrt(vec_dot(r, r, n));
  rho_tol = rho * tol_rel;

  // the outer loop to itr_max
  // let the compiler optimize the 1d loops

  for (itr = 1; itr <= itr_max; itr++) {

    // the first vector v (i.e. v[0]) is the updated residual normalized

    for (i = 0; i < n; i++) v[i + 0 * n] = r[i] / rho;

    g[0] = rho;
    for (i = 1; i <= mr; i++) g[i] = 0.0;

    // fill up h with zero

    memset(h, 0, (size_t) (mr + 1) * mr * sizeof(double));

    // the inner loop k = 1..(n-1)
    // build up the k-th Krylov space,
    // actually build the Q_k matrix of size n by k,
    //   whose the columns are k vectors v(1)...v(k)
    // remember that v[0] is computed from the updated residual as above

    for (k = 1; k <= mr; k++) {
      k_copy = k;

      // compute v(k) <- a * v(k-1)
      // here is the tricky part: v(k-1) plays a role as "charges"
      // matvec(a, v+(k-1)*n, v+k*n, n);

      apply_operator(v + (k - 1) * n, v + k * n, n);

      // compute the norm of the vector v(k)

      av = sqrt(vec_dot(v + k * n, v + k * n, n));

      // Arnoldi iteration to find v's
      // orthogonalize the k vectors v(1) . . . v(k)

      for (j = 1; j <= k; j++) {
        h[(j - 1) + (k - 1) * (mr + 1)] = vec_dot(v + k * n, v + (j - 1) * n, n);
        for (i = 0; i < n; i++)
          v[i + k * n] = v[i + k * n] - h[(j - 1) + (k - 1) * (mr + 1)] * v[i + (j - 1) * n];
      }

      // compute the norm of the newly created vector v(k)

      h[k + (k - 1) * (mr + 1)] = sqrt(vec_dot(v + k * n, v + k * n, n));

      // if the norm is close to zero, repeat the above orthogonalization

      if ((av + delta * h[k + (k - 1) * (mr + 1)]) == av) {
        for (j = 1; j <= k; j++) {
          htmp = vec_dot(v + k * n, v + (j - 1) * n, n);
          h[(j - 1) + (k - 1) * (mr + 1)] = h[(j - 1) + (k - 1) * (mr + 1)] + htmp;
          for (i = 0; i < n; i++) v[i + k * n] = v[i + k * n] - htmp * v[i + (j - 1) * n];
        }
        h[k + (k - 1) * (mr + 1)] = sqrt(vec_dot(v + k * n, v + k * n, n));
      }

      // if the norm of v(k) is nonzero, normalize v(k)

      if (h[k + (k - 1) * (mr + 1)] != 0.0) {
        for (i = 0; i < n; i++) { v[i + k * n] = v[i + k * n] / h[k + (k - 1) * (mr + 1)]; }
      }

      // if k is not the first iteration,
      // find the vector y that minimizes the norm of the residual
      //   using the least square method

      if (k > 1) {

        // update y(i-1) <- h(k-1, i-1) for i = 1...(k+1)

        for (i = 1; i <= k + 1; i++) y[i - 1] = h[(i - 1) + (k - 1) * (mr + 1)];

        // apply the Given rotation to y[j-1] and y[j] for j = 1..(k-1)

        for (j = 1; j <= k - 1; j++) mult_givens(c[j - 1], s[j - 1], j - 1, y);

        // update h(k-1, i-1) <- y(i-1) for i = 1..(k_1)

        for (i = 1; i <= k + 1; i++) h[i - 1 + (k - 1) * (mr + 1)] = y[i - 1];
      }

      // compute cosine and sine terms of the Given rotations

      mu = sqrt(h[(k - 1) + (k - 1) * (mr + 1)] * h[(k - 1) + (k - 1) * (mr + 1)] +
                h[k + (k - 1) * (mr + 1)] * h[k + (k - 1) * (mr + 1)]);
      c[k - 1] = h[(k - 1) + (k - 1) * (mr + 1)] / mu;
      s[k - 1] = -h[k + (k - 1) * (mr + 1)] / mu;

      // update h(k-1,k-1) and set h(k-1,k) to zero

      h[(k - 1) + (k - 1) * (mr + 1)] =
          c[k - 1] * h[(k - 1) + (k - 1) * (mr + 1)] - s[k - 1] * h[k + (k - 1) * (mr + 1)];
      h[k + (k - 1) * (mr + 1)] = 0;

      // apply the Givens rotation to g[k-1] and g[k]

      mult_givens(c[k - 1], s[k - 1], k - 1, g);

      // compute the norm of the residual

      rho = fabs(g[k]);

      if (rho <= rho_tol && rho <= tol_abs) break;
    }

    k = k_copy - 1;

    // compute the estimate y from h

    y[k] = g[k] / h[k + k * (mr + 1)];
    for (i = k; i >= 1; i--) {
      y[i - 1] = g[i - 1];
      for (j = i + 1; j <= k + 1; j++)
        y[i - 1] = y[i - 1] - h[(i - 1) + (j - 1) * (mr + 1)] * y[j - 1];
      y[i - 1] = y[i - 1] / h[(i - 1) + (i - 1) * (mr + 1)];
    }

    // update x at the current iteration: x <- Q(n by k) * y (k by 1)

    for (i = 1; i <= n; i++) {
      for (j = 1; j <= k + 1; j++) x[i - 1] = x[i - 1] + v[(i - 1) + (j - 1) * n] * y[j - 1];
    }

    // update the residual with the updated induced charges (x)

    update_residual(x, r, n);

    // rho = norm(r)

    rho = sqrt(vec_dot(r, r, n));

    // Barros et al. suggested the condition: norm(r) < EPSILON norm(b)

    if (rho < tol_rel * normb) break;

    // general GMRES convergence criteria

    if (rho <= rho_tol && rho <= tol_abs) break;
  }

  iterations = itr;
}

/* ----------------------------------------------------------------------
  compute the result of operator A on a given vector w
  matvec(A, v(k-1), v(k), n);
------------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::apply_operator(double *w, double *Aw, int /*n*/)
{
  int i;
  double *q = atom->q;
  double **norm = atom->mu;
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *epsilon = atom->epsilon;
  int nlocal = atom->nlocal;
  int eflag = 0;
  int vflag = 0;

  // set the induced charges to be w
  // the real charges are set to zero: Aw only involves sigma_b (not sigma_f)
  // need not to revert the induced charges
  //   because update_residual() will set the induced charges anyway

  for (i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) {
      q[i] = 0;
    } else {
      int idx = induced_charge_idx[i];
      q[i] = w[idx] * area[i];
    }
  }

  comm->forward_comm(this);

  // compute the electrical field due to w*area: y = A (w*area)

  force_clear();
  force->pair->compute(eflag, vflag);
  if (kspaceflag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();

  // now efield is the electrical field due to induced charges only
  // Note that in the definition of the electrical fields in Equations (41) and (53)
  // in Ref. Barros et al there is a factor 1/(4pi).

  for (i = 0; i < num_induced_charges; i++) buffer[i] = 0;

  for (i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) continue;

    int idx = induced_charge_idx[i];
    double Ex = efield_pair[i][0];
    double Ey = efield_pair[i][1];
    double Ez = efield_pair[i][2];
    if (kspaceflag) {
      Ex += efield_kspace[i][0];
      Ey += efield_kspace[i][1];
      Ez += efield_kspace[i][2];
    }
    double ndotE = epsilon0e2q * (Ex * norm[i][0] + Ey * norm[i][1] + Ez * norm[i][2]) / epsilon[i];
    buffer[idx] = em[i] * w[idx] + ed[i] * ndotE / (4 * MY_PI);
  }

  MPI_Allreduce(buffer, Aw, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
  need to turn the real charges back on
  set the induced charges to be w
  compute the new residual in r = b - Ax without directly computing A x
  using Eq. (60) in Barros et al.
------------------------------------------------------------------------ */

void FixPolarizeBEMGMRES::update_residual(double *w, double *r, int /*n*/)
{
  int i;
  double *q = atom->q;
  double *q_real = atom->q_unscaled;
  double **norm = atom->mu;
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *epsilon = atom->epsilon;
  int nlocal = atom->nlocal;
  int eflag = 0;
  int vflag = 0;

  // compute the Coulombic forces and electrical field E
  //   due to both ions and induced charges
  // note here q[i] = the bound charges including area + free surface charges
  //   so that kspace solver can be used directly

  for (i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) {
      q[i] = q_backup[i];
    } else {
      int idx = induced_charge_idx[i];
      q[i] = w[idx] * area[i] + q_real[i];
    }
  }

  comm->forward_comm(this);

  force_clear();
  force->pair->compute(eflag, vflag);
  if (kspaceflag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();

  // compute the residual according to Eq. (60) in Barros et al.
  // Note: in the definition of the electrical fields in Equations (41) and (53)
  // in Ref. Barros et al there is a factor 1/(4pi).

  for (i = 0; i < num_induced_charges; i++) buffer[i] = 0;

  for (i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) continue;

    int idx = induced_charge_idx[i];
    if (ed[i] == 0) {
      buffer[idx] = 0;
      continue;
    }
    double Ex = efield_pair[i][0];
    double Ey = efield_pair[i][1];
    double Ez = efield_pair[i][2];
    if (kspaceflag) {
      Ex += efield_kspace[i][0];
      Ey += efield_kspace[i][1];
      Ez += efield_kspace[i][2];
    }
    double ndotE = epsilon0e2q * (Ex * norm[i][0] + Ey * norm[i][1] + Ez * norm[i][2]) /
        epsilon[i] / (4 * MY_PI);
    double sigma_f = q_real[i] / area[i];
    buffer[idx] = (1 - em[i]) * sigma_f - em[i] * w[idx] - ed[i] * ndotE;
  }

  MPI_Allreduce(buffer, r, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::force_clear()
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

double FixPolarizeBEMGMRES::vec_dot(const double *a1, const double *a2, int n)
{
  double value = 0.0;
  for (int i = 0; i < n; i++) value += (a1[i] * a2[i]);
  return value;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixPolarizeBEMGMRES::memory_usage()
{
  double bytes = 0;
  bytes += mat_dim * sizeof(double);                   // induced_charges
  bytes += mat_dim * sizeof(double);                   // buffer
  bytes += mat_dim * sizeof(double);                   // rhs
  bytes += atom->nmax * sizeof(double);                // induced_charge_idx
  bytes += atom->nmax * sizeof(double);                // q_backup
  bytes += mr * sizeof(double);                        // c
  bytes += (mr + 1) * sizeof(double);                  // g
  bytes += (double) (mr + 1) * mr * sizeof(double);    // h
  bytes += mat_dim * sizeof(double);                   // r
  bytes += (double) mr * (mr + 1) * sizeof(double);    // s
  bytes += mat_dim * sizeof(double);                   // v
  bytes += (double) (mr + 1) * mr * sizeof(double);    // y
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::allocate()
{
  memory->create(c, mr, "polarize:c");
  memory->create(g, mr + 1, "polarize:g");
  memory->create(h, (mr + 1) * mr, "polarize:h");
  memory->create(r, mat_dim, "polarize:r");
  memory->create(s, mr, "polarize:s");
  memory->create(v, mat_dim * (mr + 1), "polarize:v");
  memory->create(y, mr + 1, "polarize:y");
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::deallocate()
{
  memory->destroy(c);
  memory->destroy(g);
  memory->destroy(h);
  memory->destroy(r);
  memory->destroy(s);
  memory->destroy(v);
  memory->destroy(y);
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMGMRES::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "itr_max") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      itr_max = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "mr") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      mr = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kspace") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      kspaceflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dielectrics") == 0) {
      if (iarg + 6 > narg) error->all(FLERR, "Illegal fix_modify command");
      double epsiloni = -1, areai = -1;
      double qreali = 0;
      int set_charge = 0;
      double ediff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      double emean = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strcmp(arg[iarg + 3], "NULL") != 0)
        epsiloni = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strcmp(arg[iarg + 4], "NULL") != 0)
        areai = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strcmp(arg[iarg + 5], "NULL") != 0) {
        qreali = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
        set_charge = 1;
      }
      set_dielectric_params(ediff, emean, epsiloni, areai, set_charge, qreali);

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

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::grow_arrays(int n)
{
  if (n > nmax) nmax = n;
  memory->grow(induced_charge_idx, nmax, "polarize:induced_charge_idx");
  memory->grow(q_backup, nmax, "polarize:q_backup");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::copy_arrays(int i, int j, int /*delflag*/)
{
  induced_charge_idx[j] = induced_charge_idx[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::set_arrays(int i)
{
  induced_charge_idx[i] = -1;
}

/* ---------------------------------------------------------------------- */

int FixPolarizeBEMGMRES::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                           int * /*pbc*/)
{
  int m;
  for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPolarizeBEMGMRES::pack_exchange(int i, double *buf)
{
  buf[0] = ubuf(induced_charge_idx[i]).d;
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixPolarizeBEMGMRES::unpack_exchange(int nlocal, double *buf)
{
  induced_charge_idx[nlocal] = (int) ubuf(buf[0]).i;
  return 1;
}

/* ----------------------------------------------------------------------
   return the number of iterations for convergence
      and current residual
------------------------------------------------------------------------- */

double FixPolarizeBEMGMRES::compute_vector(int n)
{
  if (n == 0)
    return iterations;
  else if (n == 1)
    return rho;
  else
    return 0;
}

/* ----------------------------------------------------------------------
   set dielectric params for the atoms in the group
------------------------------------------------------------------------- */

void FixPolarizeBEMGMRES::set_dielectric_params(double ediff, double emean, double epsiloni,
                                                double areai, int set_charge, double qvalue)
{
  double *area = atom->area;
  double *ed = atom->ed;
  double *em = atom->em;
  double *q_unscaled = atom->q_unscaled;
  double *epsilon = atom->epsilon;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ed[i] = ediff;
      em[i] = emean;
      if (areai > 0) area[i] = areai;
      if (epsiloni > 0) epsilon[i] = epsiloni;
      if (set_charge) q_unscaled[i] = qvalue;
    }
  }
}
