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
     based on the original implementation by Honghao Li (Northwestern)
   Reference:
     Jadhao, Solis, Olvera de la Cruz, J. Chem. Phys. 138, 054119, 2013

   Solve the following eq. for induced charges on fixed sharp interfaces:
     (Rww + Rww^T) w = q Rwq
   at every time step, the vector (q Rwq) is computed, and so
       w = [Rww + Rww^T)^(-1)] (q Rwq)
   NOTE: Oct 7, 2019: switch to using a conjugate gradient solver
------------------------------------------------------------------------- */

#include "fix_polarize_functional.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "msm_dielectric.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_coul_cut_dielectric.h"
#include "pair_coul_long_dielectric.h"
#include "pair_lj_cut_coul_cut_dielectric.h"
#include "pair_lj_cut_coul_long_dielectric.h"
#include "pair_lj_cut_coul_msm_dielectric.h"
#include "pppm_dielectric.h"
#include "random_park.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;
using namespace MathConst;
using namespace MathSpecial;

enum { REAL2SCALED = 0, SCALED2REAL = 1 };

#define EPSILON 1e-6

//#define _POLARIZE_DEBUG

/* ---------------------------------------------------------------------- */

FixPolarizeFunctional::FixPolarizeFunctional(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR, "Illegal fix polarize/functional command");

  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR, "Fix polarize/functional requires atom style dielectric");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery < 0) error->all(FLERR, "Illegal fix polarize/functional command");

  tolerance = EPSILON;
  if (narg == 5) tolerance = utils::numeric(FLERR, arg[4], false, lmp);

  comm_forward = 1;
  nmax = 0;
  allocated = 0;

  induced_charge_idx = nullptr;
  induced_charges = nullptr;
  tag2mat = nullptr;
  mat2tag = nullptr;
  tag2mat_ions = nullptr;
  mat2tag_ions = nullptr;
  ion_idx = nullptr;

  rhs1 = nullptr;
  rhs2 = nullptr;
  buffer1 = nullptr;
  buffer2 = nullptr;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;

  Rww = nullptr;
  inverse_matrix = nullptr;
  G1ww = nullptr;
  G2ww = nullptr;
  G3ww = nullptr;
  ndotGww = nullptr;

  qiRqwVector = nullptr;
  G1qw_real = nullptr;
  sum2G2wq = nullptr;

  sum1G2qw = nullptr;
  sum1G1qw_epsilon = nullptr;
  sum2ndotGwq_epsilon = nullptr;

  efield_pair = nullptr;
  efield_kspace = nullptr;

  includingG3ww = 1;

  cg_r = cg_p = cg_Ap = nullptr;
  cg_A = nullptr;

  grow_arrays(atom->nmax);
  atom->add_callback(0);    // to ensure to work with atom->sort()
}

/* ---------------------------------------------------------------------- */

FixPolarizeFunctional::~FixPolarizeFunctional()
{
  memory->destroy(mat2tag);
  memory->destroy(tag2mat);
  memory->destroy(tag2mat_ions);
  memory->destroy(mat2tag_ions);
  memory->destroy(ion_idx);
  memory->destroy(induced_charge_idx);
  memory->destroy(induced_charges);
  memory->destroy(rhs1);
  memory->destroy(rhs2);
  memory->destroy(buffer1);
  memory->destroy(buffer2);

  if (allocated) deallocate();
  atom->delete_callback(id, 0);
}

/* ---------------------------------------------------------------------- */

int FixPolarizeFunctional::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::init()
{
  // mapping induced charge matrix/vector to atom tags and vice versa

  int i, maxtag;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  tagint max_tag;
  tagint itmp;
  int *ncount = nullptr;

  // induced charge arrays setup

  max_tag = -1;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) max_tag = MAX(max_tag, tag[i]);

  MPI_Allreduce(&max_tag, &itmp, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  maxtag = (int) itmp;

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

  // ion arrays setup

  max_tag = -1;
  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) max_tag = MAX(max_tag, tag[i]);

  MPI_Allreduce(&max_tag, &itmp, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  maxtag = (int) itmp;

  memory->create(ncount, maxtag + 1, "polarize:ncount");
  for (i = 0; i <= maxtag; i++) ncount[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) ncount[tag[i]]++;

  memory->create(tag2mat_ions, maxtag + 1, "polarize:tag2mat_ions");
  MPI_Allreduce(ncount, tag2mat_ions, maxtag + 1, MPI_INT, MPI_SUM, world);

  num_ions = 0;
  for (i = 0; i <= maxtag; i++)
    if (tag2mat_ions[i])
      tag2mat_ions[i] = num_ions++;
    else
      tag2mat_ions[i] = -1;

  memory->create(mat2tag_ions, num_ions, "polarize:mat2tag_ions");
  memory->create(rhs1, num_induced_charges, "polarize:rhs1");
  memory->create(rhs2, num_induced_charges, "polarize:rhs2");
  int buffer_size = (num_induced_charges > num_ions) ? num_induced_charges : num_ions;
  memory->create(buffer1, buffer_size, num_induced_charges, "polarize:buffer1");
  memory->create(buffer2, num_induced_charges, num_induced_charges, "polarize:buffer2");
  memory->create(induced_charges, num_induced_charges, "polarize:induced_charges");

  num_ions = 0;
  for (i = 0; i <= maxtag; i++)
    if (tag2mat_ions[i] >= 0) mat2tag_ions[num_ions++] = i;

  for (i = 0; i < nlocal; i++) {
    ion_idx[i] = -1;
    if (!(mask[i] & groupbit)) ion_idx[i] = tag2mat_ions[tag[i]];
  }

  memory->destroy(ncount);

  if (allocated == 0) {
    nmax = atom->nmax;
    allocate();
    allocated = 1;
  }

  // need a full neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;

  if (force->kspace)
    g_ewald = force->kspace->g_ewald;
  else
    g_ewald = 0.01;

  if (comm->me == 0)
    utils::logmesg(lmp, "Direct solver using a variational approach for {} induced charges\n",
                   num_induced_charges);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::setup(int /*vflag*/)
{
  // check if the pair styles in use are compatible

  if (strcmp(force->pair_style, "lj/cut/coul/long/dielectric") == 0)
    efield_pair = ((PairLJCutCoulLongDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/long/dielectric/omp") == 0)
    efield_pair = ((PairLJCutCoulLongDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/msm/dielectric") == 0)
    efield_pair = ((PairLJCutCoulMSMDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/cut/dielectric") == 0)
    efield_pair = ((PairLJCutCoulCutDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "lj/cut/coul/cut/dielectric/omp") == 0)
    efield_pair = ((PairLJCutCoulCutDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "coul/long/dielectric") == 0)
    efield_pair = ((PairCoulLongDielectric *) force->pair)->efield;
  else if (strcmp(force->pair_style, "coul/cut/dielectric") == 0)
    efield_pair = ((PairCoulCutDielectric *) force->pair)->efield;
  else
    error->all(FLERR, "Pair style not compatible with fix polarize/functional");

  if (force->kspace) {

    kspaceflag = 1;
    if (strcmp(force->kspace_style, "pppm/dielectric") == 0)
      efield_kspace = ((PPPMDielectric *) force->kspace)->efield;
    else if (strcmp(force->kspace_style, "msm/dielectric") == 0)
      efield_kspace = ((MSMDielectric *) force->kspace)->efield;
    else
      error->all(FLERR, "Kspace style not compatible with fix polarize/functional");

  } else {

    if (kspaceflag == 1) {    // users specified kspace yes
      error->warning(FLERR, "No Kspace style available for fix polarize/functional");
      kspaceflag = 0;
    }
  }

  update_induced_charges();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::setup_pre_force(int /*vflag*/)
{
  // calculate Rww before the run (assuming that the interface is fixed for now)
  // otherwise this should be done every time step in pre_force()

  calculate_Rww_cutoff();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::pre_force(int)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  // solve for the induced charges

  update_induced_charges();
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::update_induced_charges()
{
  // convert all ions from scaled charges (q) to real q by multiplying with epsilon

  charge_rescaled(SCALED2REAL);

  // compute the right hand side vector qiRwVector

  calculate_qiRqw_cutoff();

  // conjugate gradient solver for w from Rww * w = -qRqw

  for (int i = 0; i < num_induced_charges; i++)
    for (int j = 0; j < num_induced_charges; j++) cg_A[i][j] = Rww[i][j] + Rww[j][i];

  for (int i = 0; i < num_induced_charges; i++) induced_charges[i] = 0;

  cg_solver(cg_A, qiRqwVector, induced_charges, num_induced_charges);

  // assign charges to the particles in the group

  double *q = atom->q;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (induced_charge_idx[i] < 0) continue;
    int idx = induced_charge_idx[i];
    q[i] = -induced_charges[idx] / (4 * MY_PI);
  }

  // revert to scaled charges to calculate forces

  charge_rescaled(REAL2SCALED);
}

/* ----------------------------------------------------------------------
  scaled2real = 1: convert ion charges from scaled values (divided by epsilon) to real values
              = 0:                          real values to scaled values
------------------------------------------------------------------------- */

void FixPolarizeFunctional::charge_rescaled(int scaled2real)
{
  double *q = atom->q;
  double *q_real = atom->q_unscaled;
  double *epsilon = atom->epsilon;
  int nlocal = atom->nlocal;

  if (scaled2real) {
    for (int i = 0; i < nlocal; i++)
      if (induced_charge_idx[i] < 0) q[i] = q_real[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (induced_charge_idx[i] < 0) q[i] = q_real[i] / epsilon[i];
  }

  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::allocate()
{
  // initialize all data
  // interface terms, all matrix of M*M
  memory->create(inverse_matrix, num_induced_charges, num_induced_charges, "fix:inverse_matrix");
  memory->create(Rww, num_induced_charges, num_induced_charges, "fix:Rww");
  memory->create(G1ww, num_induced_charges, num_induced_charges, "fix:G1ww");
  memory->create(ndotGww, num_induced_charges, num_induced_charges, "fix:ndotGww");
  memory->create(G2ww, num_induced_charges, num_induced_charges, "fix:G2ww");
  memory->create(G3ww, num_induced_charges, num_induced_charges, "fix:G3ww");

  // each step, qw, qq terms, temp data

  memory->create(qiRqwVector, num_induced_charges, "fix:qiRqwVector");
  memory->create(sum2G2wq, num_induced_charges, "fix:sum2G2wq");
  memory->create(G1qw_real, num_ions, num_induced_charges, "fix:G1qw_real");

  // arrays of M
  memory->create(sum1G2qw, num_induced_charges, "fix:sum1G2qw");
  memory->create(sum1G1qw_epsilon, num_induced_charges, "fix:sum1G1qw_epsilon");
  memory->create(sum2ndotGwq_epsilon, num_induced_charges, "fix:sum2ndotGwq_epsilon");

  memory->create(cg_r, num_induced_charges, "polarize:cg_r");
  memory->create(cg_p, num_induced_charges, "polarize:cg_p");
  memory->create(cg_Ap, num_induced_charges, "polarize:cg_Ap");
  memory->create(cg_A, num_induced_charges, num_induced_charges, "polarize:cg_A");
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::deallocate()
{
  memory->destroy(inverse_matrix);
  memory->destroy(Rww);
  memory->destroy(G1ww);
  memory->destroy(G2ww);
  memory->destroy(G3ww);
  memory->destroy(ndotGww);

  memory->destroy(qiRqwVector);
  memory->destroy(sum2G2wq);
  memory->destroy(G1qw_real);

  memory->destroy(sum1G2qw);
  memory->destroy(sum1G1qw_epsilon);
  memory->destroy(sum2ndotGwq_epsilon);

  memory->destroy(cg_r);
  memory->destroy(cg_p);
  memory->destroy(cg_Ap);
  memory->destroy(cg_A);
}

/* ---------------------------------------------------------------------- */

int FixPolarizeFunctional::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kspace") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix_modify command");
      if (strcmp(arg[iarg + 1], "yes") == 0)
        kspaceflag = 1;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        kspaceflag = 0;
      else
        error->all(FLERR, "Illegal fix_modify command for fix polarize/functional");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dielectrics") == 0) {
      if (iarg + 6 > narg) error->all(FLERR, "Illegal fix_modify command");
      double epsiloni = -1, areai = -1;
      double q_unscaled = 0;
      int set_charge = 0;
      double ediff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      double emean = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strcmp(arg[iarg + 3], "nullptr") != 0)
        epsiloni = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strcmp(arg[iarg + 4], "nullptr") != 0)
        areai = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strcmp(arg[iarg + 5], "nullptr") != 0) {
        q_unscaled = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
        set_charge = 1;
      }
      set_dielectric_params(ediff, emean, epsiloni, areai, set_charge, q_unscaled);

      iarg += 6;
    } else
      error->all(FLERR, "Illegal fix_modify command");
  }

  return iarg;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPolarizeFunctional::pack_exchange(int i, double *buf)
{
  buf[0] = ubuf(induced_charge_idx[i]).d;
  buf[1] = ubuf(ion_idx[i]).d;
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixPolarizeFunctional::unpack_exchange(int nlocal, double *buf)
{
  induced_charge_idx[nlocal] = (int) ubuf(buf[0]).i;
  ion_idx[nlocal] = (int) ubuf(buf[1]).i;
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixPolarizeFunctional::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m;
  for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPolarizeFunctional::grow_arrays(int n)
{
  if (n > nmax) nmax = n;
  memory->grow(induced_charge_idx, nmax, "fix:induced_charge_idx");
  memory->grow(ion_idx, nmax, "fix:ion_idx");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPolarizeFunctional::copy_arrays(int i, int j, int /*delflag*/)
{
  induced_charge_idx[j] = induced_charge_idx[i];
  ion_idx[j] = ion_idx[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPolarizeFunctional::set_arrays(int i)
{
  induced_charge_idx[i] = -1;
  ion_idx[i] = 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixPolarizeFunctional::memory_usage()
{
  double bytes = 0;
  bytes += square(num_induced_charges) * sizeof(double);            // inverse_matrix
  bytes += square(num_induced_charges) * sizeof(double);            // Rww
  bytes += square(num_induced_charges) * sizeof(double);            // G1ww
  bytes += square(num_induced_charges) * sizeof(double);            // ndotGww
  bytes += square(num_induced_charges) * sizeof(double);            // G2ww
  bytes += square(num_induced_charges) * sizeof(double);            // G3ww
  bytes += num_induced_charges * sizeof(double);                    // qiRqwVector
  bytes += num_induced_charges * sizeof(double);                    // sum2G2wq
  bytes += num_induced_charges * sizeof(double);                    // sum1G2qw
  bytes += num_induced_charges * sizeof(double);                    // sum1G1qw_epsilon
  bytes += num_induced_charges * sizeof(double);                    // sum2ndotGwq_epsilon
  bytes += (double)num_ions * num_induced_charges * sizeof(double); // G1qw_real
  bytes += nmax * sizeof(int);                                      // induced_charge_idx
  bytes += nmax * sizeof(int);                                      // ion_idx
  bytes += num_induced_charges * sizeof(double);                    // induced_charges
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::calculate_Rww_cutoff()
{
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *area = atom->area;
  double *curvature = atom->curvature;
  double **norm = atom->mu;
  double *ed = atom->ed;
  double *em = atom->em;

  // invoke full neighbor list

  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate G1ww, gradG1ww, ndotG1ww
  // fill up buffer1 with local G1ww and buffer2 with local ndotGww
  // seperate into two loops to let the compiler optimize/or later vectorization

  for (int i = 0; i < num_induced_charges; i++)
    for (int j = 0; j < num_induced_charges; j++) buffer1[i][j] = 0;

  for (int i = 0; i < num_induced_charges; i++)
    for (int j = 0; j < num_induced_charges; j++) buffer2[i][j] = 0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = induced_charge_idx[i];
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {

          // interface particles: k can be ghost atoms
          double delx = xtmp - x[k][0];
          double dely = ytmp - x[k][1];
          double delz = ztmp - x[k][2];
          domain->minimum_image(delx, dely, delz);
          int mk = tag2mat[tag[k]];

          // G1ww[mi][mk] = calculate_greens_ewald(delx, dely, delz);
          buffer1[mi][mk] = calculate_greens_ewald(delx, dely, delz);

          // gradG1ww is vector, directly change it in the function
          double gradG1ww[3];
          calculate_grad_greens_ewald(gradG1ww, delx, dely, delz);

          // use mu to store the normal vector of interface vertex
          buffer2[mi][mk] = MathExtra::dot3(norm[i], gradG1ww) / (4 * MY_PI);
        }
      }

      // special treatment for the diagonal terms,
      // even though in the above loop there is mk == mi

      buffer1[mi][mi] = calculate_greens_ewald_self_vertex(area[i]);
      buffer2[mi][mi] = calculate_ndotgreens_ewald_self_vertex(area[i], curvature[i]) / (4 * MY_PI);
    }
  }

  MPI_Allreduce(buffer1[0], G1ww[0], num_induced_charges * num_induced_charges, MPI_DOUBLE, MPI_SUM,
                world);
  MPI_Allreduce(buffer2[0], ndotGww[0], num_induced_charges * num_induced_charges, MPI_DOUBLE,
                MPI_SUM, world);

  // calculate G2ww
  // fill up buffer1 with local G2ww

  for (int i = 0; i < num_induced_charges; i++)
    for (int j = 0; j < num_induced_charges; j++) buffer1[i][j] = 0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = induced_charge_idx[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;

        if (mask[k] & groupbit) {
          // interface particles: k can be ghost atoms
          int mk = tag2mat[tag[k]];
          double temp = 0;
          for (int ll = 0; ll < jnum; ll++) {
            int l = jlist[ll] & NEIGHMASK;
            if (mask[l] & groupbit) {
              // interface particles: l can be ghost atoms
              int ml = tag2mat[tag[l]];
              temp += G1ww[mi][ml] * ndotGww[ml][mk] * area[l] * ed[l];
            }
          }
          //G2ww[mi][mk] = temp;
          buffer1[mi][mk] = temp;
        }
      }

      double temp = 0;
      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {
          // interface particles: k can be ghost atoms
          int mk = tag2mat[tag[k]];
          temp += G1ww[mi][mk] * ndotGww[mk][mi] * area[k] * ed[k];
        }
      }
      //G2ww[mi][mi] = temp;
      buffer1[mi][mi] = temp;
    }
  }

  MPI_Allreduce(buffer1[0], G2ww[0], num_induced_charges * num_induced_charges, MPI_DOUBLE, MPI_SUM,
                world);

  // calculate G3ww and Rww
  // G3ww is implemented as in _exact(), but can be optionally excluded
  // due to its minor contribution
  // fill up buffer1 with local G3ww

  for (int i = 0; i < num_induced_charges; i++)
    for (int j = 0; j < num_induced_charges; j++) buffer1[i][j] = 0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = induced_charge_idx[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;

        if (mask[k] & groupbit) {

          // interface particles: k can be ghost atoms

          int mk = tag2mat[tag[k]];

          double a1 = em[i] * (em[k] - 1.0);
          double a2 = 1.0 - em[i] - em[k];

          // The first term (w/ G1ww) contributes the most to Rww
          // the second term (w/ G2ww) includes certain correction

          //Rww[mi][mk] = a1 * G1ww[mi][mk] + a2 * G2ww[mi][mk];
          buffer1[mi][mk] = a1 * G1ww[mi][mk] + a2 * G2ww[mi][mk];

          if (includingG3ww) {
            double temp = 0;
            for (int ll = 0; ll < jnum; ll++) {
              int l = jlist[ll] & NEIGHMASK;
              if (mask[l] & groupbit) {
                // interface particles: l can be ghost atoms
                int ml = tag2mat[tag[l]];
                temp += (ndotGww[ml][mi]) * G2ww[ml][mk] * area[l] * ed[l];
              }
            }
            G3ww[mi][mk] = temp;
            //Rww[mi][mk] += G3ww[mi][mk];
            buffer1[mi][mk] += G3ww[mi][mk];
          }
        }
      }

      if (includingG3ww) {
        double temp = 0;
        for (int ll = 0; ll < jnum; ll++) {
          int l = jlist[ll] & NEIGHMASK;
          if (mask[l] & groupbit) {
            // interface particles: l can be ghost atoms
            int ml = tag2mat[tag[l]];
            temp += (ndotGww[ml][mi]) * G2ww[ml][mi] * area[l] * ed[l];
          }
        }
        G3ww[mi][mi] = temp;
        //Rww[mi][mi] += G3ww[mi][mi];
        buffer1[mi][mi] += G3ww[mi][mi];
      }

      // including the diagonal term
      double a1 = em[i] * (em[i] - 1.0);
      double a2 = 1.0 - em[i] - em[i];

      // The first term (w/ G1ww) contributes the most to Rww
      // the second term (w/ G2ww) includes certain correction
      //Rww[mi][mi] = a1 * G1ww[mi][mi] + a2 * G2ww[mi][mi];

      buffer1[mi][mi] = a1 * G1ww[mi][mi] + a2 * G2ww[mi][mi];
    }
  }

  MPI_Allreduce(buffer1[0], Rww[0], num_induced_charges * num_induced_charges, MPI_DOUBLE, MPI_SUM,
                world);

#ifdef _POLARIZE_DEBUG
  if (comm->me == 0) {
    FILE *fp = fopen("Rww-functional.txt", "w");
    for (int i = 0; i < num_induced_charges; i++)
      fprintf(fp, "%d %g %g %g\n", i, Rww[i][i], Rww[i][num_induced_charges / 2],
              Rww[num_induced_charges / 2][i]);
    fclose(fp);
  }
#endif
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::calculate_qiRqw_cutoff()
{
  int ii, i, k, kk, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, r;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  double **x = atom->x;
  double *q = atom->q_unscaled;
  double *epsilon = atom->epsilon;
  double *area = atom->area;
  double **norm = atom->mu;
  double *ed = atom->ed;
  double *em = atom->em;

  // invoke full neighbor list

  int inum, *ilist, *jlist, *numneigh, **firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate G1qw_real
  // fill up buffer1 with local G1qw_real

  for (int i = 0; i < num_ions; i++)
    for (int j = 0; j < num_induced_charges; j++) buffer1[i][j] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) {
      // ion particles
      int mi = ion_idx[i];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (kk = 0; kk < jnum; kk++) {
        k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {
          // interface particles: k can be ghost atoms
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          domain->minimum_image(delx, dely, delz);
          r = sqrt(delx * delx + dely * dely + delz * delz);

          int mk = tag2mat[tag[k]];
          //G1qw_real[mi][mk] = greens_real(r);
          buffer1[mi][mk] = greens_real(r);
        }
      }
    }
  }

  MPI_Allreduce(&buffer1[0][0], &G1qw_real[0][0], num_ions * num_induced_charges, MPI_DOUBLE,
                MPI_SUM, world);

  // the following loop need the above results: gradG1wq_real
  // calculate sum1G1qw_epsilon and sum2ndotGwq_epsilon
  // fill up rhs1 with local sum1G1qw_epsilon and rhs2 with local sum2ndotGwq_epsilon

  memset(rhs1, 0, num_induced_charges * sizeof(double));
  memset(rhs2, 0, num_induced_charges * sizeof(double));

  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];    // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      int mk = induced_charge_idx[k];
      xtmp = x[k][0];
      ytmp = x[k][1];
      ztmp = x[k][2];
      jlist = firstneigh[k];
      jnum = numneigh[k];

      double tempndotG[3] = {0.0, 0.0, 0.0};
      double temp_sum1 = 0;
      for (ii = 0; ii < jnum; ii++) {
        i = jlist[ii] & NEIGHMASK;
        if (!(mask[i] & groupbit)) {
          // ions particles: i can be ghost atoms
          delx = x[i][0] - xtmp;
          dely = x[i][1] - ytmp;
          delz = x[i][2] - ztmp;
          domain->minimum_image(delx, dely, delz);

          int mi = tag2mat_ions[tag[i]];    //ion_idx[i];

          //calculate_grad_greens_real(gradG1wq_real[mk][mi], delx, dely, delz);
          double gradG1wq[3];
          calculate_grad_greens_real(gradG1wq, delx, dely, delz);
          MathExtra::scale3(-1.0, gradG1wq);

          tempndotG[0] += gradG1wq[0] * (q[i] / epsilon[i]);
          tempndotG[1] += gradG1wq[1] * (q[i] / epsilon[i]);
          tempndotG[2] += gradG1wq[2] * (q[i] / epsilon[i]);
          temp_sum1 += G1qw_real[mi][mk] * q[i] / epsilon[i];
        }
      }

      //sum1G1qw_epsilon[mk] = temp_sum1;// + ewaldDielectric->sum1G1qw_k_epsilon[mk];
      rhs1[mk] = temp_sum1;

      //sum2ndotGwq_epsilon[mk] = MathExtra::dot3(norm[k], tempndotG);
      rhs2[mk] = MathExtra::dot3(norm[k], tempndotG);
    }
  }

  MPI_Allreduce(rhs1, sum1G1qw_epsilon, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(rhs2, sum2ndotGwq_epsilon, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);

  // calculate G2, gradient G2
  // sum2G2wq and sum1G2qw

  //  for (int i = 0; i < num_induced_charges; i++) rhs1[i] = rhs2[i] = 0;
  memset(rhs1, 0, num_induced_charges * sizeof(double));
  memset(rhs2, 0, num_induced_charges * sizeof(double));

  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];    // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      int mk = induced_charge_idx[k];
      jlist = firstneigh[k];
      jnum = numneigh[k];

      double tempwq = 0;
      double temp = 0;
      for (ii = 0; ii < jnum; ii++) {
        i = jlist[ii] & NEIGHMASK;
        if (mask[i] & groupbit) {
          // interface particles: i can be ghost atoms
          int mi = tag2mat[tag[i]];
          tempwq += G1ww[mk][mi] * (sum2ndotGwq_epsilon[mi]) * area[i] * ed[i];
          temp += sum1G1qw_epsilon[mi] * (ndotGww[mi][mk]) * area[i] * ed[i];
        }
      }

      // add the corresponding self terms
      tempwq += G1ww[mk][mk] * (sum2ndotGwq_epsilon[mk]) * area[k] * ed[k];
      temp += sum1G1qw_epsilon[mk] * (ndotGww[mk][mk]) * area[k] * ed[k];

      //sum2G2wq[mk] = tempwq;
      rhs1[mk] = tempwq;
      //sum1G2qw[mk] = temp;
      rhs2[mk] = temp;
    }
  }

  MPI_Allreduce(rhs1, sum2G2wq, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(rhs2, sum1G2qw, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);

  // calculate G3, gradient G3
  // fill up rhs1 with local qiRqwVector

  memset(rhs1, 0, num_induced_charges * sizeof(double));

  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];    // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      int mk = induced_charge_idx[k];
      jlist = firstneigh[k];
      jnum = numneigh[k];

      double sum1G3qw = 0;
      double qiRwwVectorTemp1 = 0;
      for (ii = 0; ii < jnum; ii++) {
        i = jlist[ii] & NEIGHMASK;
        if (mask[i] & groupbit) {
          // interface particles: i can be ghost atoms
          int mi = tag2mat[tag[i]];
          sum1G3qw += sum2ndotGwq_epsilon[mi] * G2ww[mi][mk] * area[i] * ed[i];
        } else {
          // ions particles: i can be ghost atoms
          int mi = tag2mat_ions[tag[i]];    //ion_idx[i];
          qiRwwVectorTemp1 += q[i] * (1.0 - em[k] / epsilon[i]) * G1qw_real[mi][mk];
        }
      }

      // add the diagonal term

      sum1G3qw += sum2ndotGwq_epsilon[mk] * G2ww[mk][mk] * area[k] * ed[k];

      // qiRwwVectorTemp2 is a significant contribution, of which sum2G2wq is significant
      double qiRwwVectorTemp2 = (1.0 - 2.0 * em[k]) * sum2G2wq[mk] + sum1G2qw[mk] + 2.0 * sum1G3qw;

      // qiRqwVector[mk] = qiRwwVectorTemp1 + qiRwwVectorTemp2;
      rhs1[mk] = qiRwwVectorTemp1 + qiRwwVectorTemp2;
    }
  }

  MPI_Allreduce(rhs1, qiRqwVector, num_induced_charges, MPI_DOUBLE, MPI_SUM, world);

#ifdef _POLARIZE_DEBUG
  if (comm->me == 0) {
    FILE *fp = fopen("qRqw-functional.txt", "w");
    for (int i = 0; i < num_induced_charges; i++) fprintf(fp, "%d %g\n", i, qiRqwVector[i]);
    fclose(fp);
  }
#endif
}

/* ----------------------------------------------------------------------
   set dielectric params for the atom in the group
------------------------------------------------------------------------- */

void FixPolarizeFunctional::set_dielectric_params(double ediff, double emean, double epsiloni,
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

/* ----------------------------------------------------------------------
   real Green's function
------------------------------------------------------------------------ */

double FixPolarizeFunctional::greens_real(double r)
{
  return erfc(g_ewald * r) / r;
}

/* ---------------------------------------------------------------------- */

double FixPolarizeFunctional::grad_greens_real_factor(double r)
{
  double alpharij = g_ewald * r;
  double factor = erfc(alpharij) + 2.0 * alpharij / MY_PIS * exp(-(alpharij * alpharij));
  double r3 = r * r * r;
  return (factor * (-1.0 / r3));
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::calculate_grad_greens_real(double *vec, double dx, double dy, double dz)
{
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  double real = grad_greens_real_factor(r);
  vec[0] = real * dx;
  vec[1] = real * dy;
  vec[2] = real * dz;
}

/* ---------------------------------------------------------------------- */

double FixPolarizeFunctional::calculate_greens_ewald(double dx, double dy, double dz)
{
  // excluding the reciprocal term
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  return greens_real(r);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::calculate_grad_greens_ewald(double *vec, double dx, double dy,
                                                        double dz)
{
  // real part of grad greens, excluding the reciprocal term
  calculate_grad_greens_real(vec, dx, dy, dz);
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::calculate_matrix_multiply_vector(double **matrix, double *in_vec,
                                                             double *out_vec, int M)
{
  for (int k = 0; k < M; ++k) {
    double temp = 0.0;
    for (int l = 0; l < M; ++l) { temp += matrix[k][l] * in_vec[l]; }
    out_vec[k] = temp;
  }
}

/* ---------------------------------------------------------------------- */

double FixPolarizeFunctional::calculate_greens_ewald_self_vertex(double area)
{
  // excluding the reciprocal term
  double corr = 2.0 * MY_PIS / sqrt(area);
  double self_energy = -2.0 * g_ewald / MY_PIS;
  return corr + self_energy;
}

/* ---------------------------------------------------------------------- */

double FixPolarizeFunctional::calculate_ndotgreens_ewald_self_vertex(double area, double curvature)
{
  // this term is important, cannot be set to zero
  // curvature = 1 / R, minus if norm is inverse of R to center.

  return curvature * MY_PIS / sqrt(area);
}

/* ----------------------------------------------------------------------
  compute the inner product between two vectors x and y: x^t * y
    where ^t is the transpose operator
-- ---------------------------------------------------------------------- */

double FixPolarizeFunctional::inner_product(double *x, double *y, int N)
{
  double t = 0;
  for (int i = 0; i < N; i++) t += x[i] * y[i];
  return t;
}

/* ---------------------------------------------------------------------- */

void FixPolarizeFunctional::cg_solver(double **A, double *b, double *x, int N)
{
  calculate_matrix_multiply_vector(A, x, cg_p, N);
  for (int i = 0; i < N; i++) {
    cg_r[i] = b[i] - cg_p[i];
    cg_p[i] = cg_r[i];
  }
  double rsq = inner_product(cg_r, cg_r, N);

  // maximum number of iterations do not exceed N
  for (int k = 0; k < N; k++) {

    // Ap = A * p
    calculate_matrix_multiply_vector(A, cg_p, cg_Ap, N);

    // pAp = p^t * Ap
    double pAp = inner_product(cg_p, cg_Ap, N);

    // alpha = r^t * r / pAp
    double alpha = rsq / pAp;

    // x = x + alpha * p
    // r = r - alpha * Ap
    for (int i = 0; i < N; i++) {
      x[i] = x[i] + alpha * cg_p[i];
      cg_r[i] = cg_r[i] - alpha * cg_Ap[i];
    }

    double rsq_new = inner_product(cg_r, cg_r, N);
    if (rsq_new < tolerance) break;

    // beta = rsq_new / rsq
    double beta = rsq_new / rsq;
    for (int i = 0; i < N; i++) cg_p[i] = cg_r[i] + beta * cg_p[i];
    rsq = rsq_new;
  }
}
