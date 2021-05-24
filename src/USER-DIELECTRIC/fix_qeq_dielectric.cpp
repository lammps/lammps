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
   Contributing authors:
     full-matrix: Honghao Li (Northwestern University)
     neighbor list: Trung Nguyen  (Northwestern University)
   Ref.: Jadhao, Solis, Olvera de la Cruz, J. Chem. Phys. 138, 054119, 2013
------------------------------------------------------------------------- */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include "fix_qeq_dielectric.h"
#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "kspace.h"
#include "group.h"
#include "pair.h"
#include "respa.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "citeme.h"
#include "error.h"
#include <gsl/gsl_linalg.h>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

#define _POLARIZE_DEBUG

/* ---------------------------------------------------------------------- */

FixQEqDielectric::FixQEqDielectric(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix qeq/dielectric command");

  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR,"fix qeq/dielectric requires atom style dielectric");

  comm_forward = 1;

  // register with Atom class
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // kspace info
//  ewaldDielectric = NULL;

  tags_interface = NULL;
  tags_ions = NULL;

  // interface terms
  Rww = NULL;
  inverse_matrix = NULL;
  G1ww = NULL;
  G2ww = NULL;
  G3ww = NULL;
  ndotGww = NULL;
  gradG1ww = NULL;

  qiRqwVector = NULL;
  G1qq_real = NULL;
  G1qw_real = NULL;
  gradG1wq_real = NULL;
  sum2G2wq = NULL;

  sum1G2qw = NULL;
  sum1G3qw = NULL;
  sum1G1qw_epsilon = NULL;
  sum2ndotGwq_epsilon = NULL;

  q_induced_charges = NULL;

  full = 1;
  includingG3ww = 1;
}

/* ---------------------------------------------------------------------- */

FixQEqDielectric::~FixQEqDielectric()
{
  // unregister callbacks to this fix from Atom class

  if (copymode) return;

  atom->delete_callback(id,0);

//  if (ewaldDielectric) delete ewaldDielectric;
//  ewaldDielectric = NULL;

  memory->destroy(tags_interface);
  memory->destroy(tags_ions);

  memory->destroy(inverse_matrix);
  memory->destroy(Rww);
  memory->destroy(G1ww);
  memory->destroy(G2ww);
  memory->destroy(G3ww);
  memory->destroy(ndotGww);
  memory->destroy3d_offset(gradG1ww,0);

  memory->destroy(qiRqwVector);
  memory->destroy(sum2G2wq);
  memory->destroy(G1qq_real);
  memory->destroy(G1qw_real);
  memory->destroy3d_offset(gradG1wq_real,0);

  memory->destroy(sum1G2qw);
  memory->destroy(sum1G3qw);
  memory->destroy(sum1G1qw_epsilon);
  memory->destroy(sum2ndotGwq_epsilon);

  memory->destroy(q_induced_charges);

}

/* ---------------------------------------------------------------------- */

int FixQEqDielectric::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::init()
{
  if (atom->map_style == 0)
    error->all(FLERR,"Fix qeq/dielectric requires an atom map, see atom_modify");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/dielectric group has no atoms");

  if (comm->nprocs > 1) error->all(FLERR,"Fix qeq/dielectric works with 1 MPI for now");

  // need a full neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // TODO: some data needed later, some just temp should be defined and release in setup_pre_force.
  // initial interfaces terms
  n_induced_charges = ngroup;
  n_ions = atom->nlocal - n_induced_charges;
  if (force->kspace) g_ewald = force->kspace->g_ewald;
  else g_ewald = 0.01;

  // kspace
//  int narg = 1;
//  char *arg[10] = {"0.0000"}; // this does not matter, the accuracy_relative is manually changed to kspace accuracy_relative in constructor.
//  ewaldDielectric = new EwaldDielectric(n_induced_charges, n_ions, lmp, narg, arg);
  // tags, after ewaldDielectric created, tags will be transfer to that class too.
  memory->create(tags_interface, n_induced_charges, "fix:tags_interface");
  memory->create(tags_ions, n_ions, "fix:tags_ions");

  // initialize all data
  // interface terms, all matrix of M*M
  memory->create(inverse_matrix, n_induced_charges, n_induced_charges, "fix:inverse_matrix");
  memory->create(Rww, n_induced_charges, n_induced_charges, "fix:Rww");
  memory->create(G1ww, n_induced_charges, n_induced_charges, "fix:G1ww");
  memory->create(ndotGww, n_induced_charges, n_induced_charges, "fix:ndotGww");
  memory->create(G2ww, n_induced_charges, n_induced_charges, "fix:G2ww");
  memory->create(G3ww, n_induced_charges, n_induced_charges, "fix:G3ww");
  memory->create3d_offset(gradG1ww,0,n_induced_charges,n_induced_charges,3,"fix:gradG1ww");

  // each step, qw, qq terms, temp data
  memory->create(qiRqwVector, n_induced_charges, "fix:qiRqwVector");
  memory->create(sum2G2wq, n_induced_charges, "fix:sum2G2wq");
  memory->create(G1qq_real, n_ions, n_ions, "fix:G1qq_real");
  memory->create(G1qw_real, n_ions, n_induced_charges, "fix:G1qw_real");
  memory->create3d_offset(gradG1wq_real,0,n_induced_charges,n_ions,3,"fix:gradG1wq_real");
  // array of M
  memory->create(sum1G2qw, n_induced_charges, "fix:sum1G2qw");
  memory->create(sum1G3qw, n_induced_charges, "fix:sum1G3qw");
  memory->create(sum1G1qw_epsilon, n_induced_charges, "fix:sum1G1qw_epsilon");
  memory->create(sum2ndotGwq_epsilon, n_induced_charges, "fix:sum2ndotGwq_epsilon");

  memory->create(q_induced_charges, n_induced_charges, "fix:q_induced_charges");

  setup_tags();
//  print_all_properties();
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::setup_pre_force(int vflag)
{
//  printf("========= calculate structure factors and kspace terms ==========\n");
  // calculate structure factor and kspace terms first
  // TODO this step seems not necessary, check it later
  // since the 3-loop calculate explicit kspace term G(sk, sm), so the structure factor stuff is not necessary
//  ewaldDielectric->calculate_structure_factors();
//  ewaldDielectric->calculate_kspace_terms();

  printf("==== Start setup_pre_force in fix_qeq_dielectric ======\n");
  clock_t t = clock();
  #ifdef _OPENMP
  double startTime = omp_get_wtime();
  #endif

  if (full) calculate_Rww_full();
  else calculate_Rww_cutoff();

  // calculate inverse matrix (Rww + Rww^T)^(-1),
  printf("==== Calculate the inverse of matrix Rww ======\n");
  calculate_inverse_matrix(Rww, inverse_matrix, n_induced_charges);

  pre_force(vflag);

  printf("==== End setup_pre_force in fix_qeq_dielectric ======\n");

}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_Rww_full()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *area = avec->area;
  double *curvature = avec->curvature;
  double **norm = avec->mu;
  double *ed = avec->ed;
  double *em = avec->em;

  // calculate interface term Rww and the inverse matrix (Rww + Rww^T)^(-1);

  #pragma omp parallel
  {

  // calculate G1ww, gradG1ww, ndotG1ww
  // =========== Loop over all nlocal ===============, we need all pair G^Kspace(i, j)
  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    // nk: all nlocal index in atom[] starts with n---, k is the index of all the matrix
    int nk = get_index_interface(k);
    double xtmp = x[nk][0];
    double ytmp = x[nk][1];
    double ztmp = x[nk][2];

    for (int l = 0; l < n_induced_charges; ++l) {
      int nl = get_index_interface(l);

      // ik = il, dr = 0, different equation
      if (nk == nl) {
        G1ww[k][l] = calculate_greens_ewald_self_vertex(area[nk]);
        ndotGww[k][l] = calculate_ndotgreens_ewald_self_vertex(area[nk], curvature[nk]) / (4*MY_PI);
        continue;
      } else {
          // k != l
        double delx = xtmp - x[nl][0];
        double dely = ytmp - x[nl][1];
        double delz = ztmp - x[nl][2];

        domain->minimum_image(delx,dely,delz);

        G1ww[k][l] = calculate_greens_ewald(delx, dely, delz);
        // gradG1ww is vector, directly change it in the function
        calculate_grad_greens_ewald(gradG1ww[k][l], delx, dely, delz);
        // use mu to store the normal vector of interface vertex
        ndotGww[k][l] = MathExtra::dot3(norm[nk], gradG1ww[k][l]) / (4*MY_PI);
      }
    }
  }

  // calculate G2ww
  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    for (int l = 0; l < n_induced_charges; ++l) {
      double temp = 0;
      for (int m = 0; m < n_induced_charges; ++m) {
        int nm = get_index_interface(m);
        temp += G1ww[k][m] * ndotGww[m][l] * area[nm] * ed[nm];
      }
      G2ww[k][l] = temp;
    }
  }

  // calculate G3ww and Rww
  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);
    for (int l = 0; l < n_induced_charges; ++l) {
      int nl = get_index_interface(l);
      double temp = 0;
      for (int m = 0; m < n_induced_charges; ++m) {
        int nm = get_index_interface(m);
        temp += (ndotGww[m][k]) * G2ww[m][l] * area[nm] * ed[nm];
      }
      G3ww[k][l] = temp;
      double a1 = em[nk] * (em[nl] - 1.0);
      double a2 = 1.0 - em[nk] - em[nl];
      // The first term (w/ G1ww) contributes the most to Rww
      // the second term (w/ G2ww) includes certain correction
      // the third term (w/ G3ww) corresponds to a minor contribution
      Rww[k][l] = a1 * G1ww[k][l] + a2 * G2ww[k][l] + G3ww[k][l];
    }
  }
  } // end of the whole parallel region

  #ifdef _POLARIZE_DEBUG
  if (comm->me == 0) {
    FILE* fp = fopen("Rww-qeq-full.txt", "w");
    for (int i = 0; i < n_induced_charges; i++)
      fprintf(fp, "%d %g\n", i, Rww[i][i]);
    fclose(fp);
  }
  #endif
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_Rww_cutoff()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *area = avec->area;
  double *curvature = avec->curvature;
  double **norm = avec->mu;
  double *ed = avec->ed;
  double *em = avec->em;

  // invoke full neighbor list (will copy or build if necessary)

  int inum,jnum,*ilist,*jlist,*numneigh,**firstneigh;

  inum = list->inum;               // number of entries in the neighbor list
  ilist = list->ilist;             // ilist[ii] gives the atom index of the entry ii in the list
  numneigh = list->numneigh;       // numneigh[i] gives the number of neighbors of local atom i
  firstneigh = list->firstneigh;   // firstneigh[i] gives the pointer to the neighbors of i

  // calculate G1ww, gradG1ww, ndotG1ww

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = get_matrix_index_from_local_index(i);
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {
          // interface particles
          double delx = xtmp - x[k][0];
          double dely = ytmp - x[k][1];
          double delz = ztmp - x[k][2];
          domain->minimum_image(delx,dely,delz);
          int mk = get_matrix_index_from_local_index(k);

          G1ww[mi][mk] = calculate_greens_ewald(delx, dely, delz);
          // gradG1ww is vector, directly change it in the function
          calculate_grad_greens_ewald(gradG1ww[mi][mk], delx, dely, delz);
          // use mu to store the normal vector of interface vertex
          ndotGww[mi][mk] = MathExtra::dot3(norm[i], gradG1ww[mi][mk]) / (4*MY_PI);
        }
      }

      // special treatment for the diagonal terms because in the above loop there is no mk == mi

      G1ww[mi][mi] = calculate_greens_ewald_self_vertex(area[i]);
      ndotGww[mi][mi] = calculate_ndotgreens_ewald_self_vertex(area[i], curvature[i]) / (4*MY_PI);
    }
  }

  // calculate G2ww

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = get_matrix_index_from_local_index(i);
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;

        if (mask[k] & groupbit) {
          // interface particles
          int mk = get_matrix_index_from_local_index(k);
          double temp = 0;
          for (int ll = 0; ll < jnum; ll++) {
            int l = jlist[ll] & NEIGHMASK;
            if (mask[l] & groupbit) {
              // interface particles
              int ml = get_matrix_index_from_local_index(l);
              temp += G1ww[mi][ml] * ndotGww[ml][mk] * area[l] * ed[l];
            }
          }
          G2ww[mi][mk] = temp;
        }
      }

      // including the diagonal term
      double temp = 0;
      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {
          // interface particles
          int mk = get_matrix_index_from_local_index(k);
          temp += G1ww[mi][mk] * ndotGww[mk][mi] * area[k] * ed[k];
        }
      }
      G2ww[mi][mi] = temp;
    }
  }

  // calculate G3ww and Rww
/*
  // G3ww is implemented as in _exact(), but can be optionally excluded
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);
    for (int l = 0; l < n_induced_charges; ++l) {
      int nl = get_index_interface(l);
      double a1 = em[nk] * (em[nl] - 1.0);
      double a2 = 1.0 - em[nk] - em[nl];
      // The first term (w/ G1ww) contributes the most to Rww
      // the second term (w/ G2ww) includes certain correction
      Rww[k][l] = a1 * G1ww[k][l] + a2 * G2ww[k][l];

      // the third term (w/ G3ww) corresponds to a minor contribution
      if (includingG3ww) {
        double temp = 0;
        for (int m = 0; m < n_induced_charges; ++m) {
          int nm = get_index_interface(m);
          temp += (ndotGww[m][k]) * G2ww[m][l] * area[nm] * ed[nm];
        }
        G3ww[k][l] = temp;
        Rww[k][l] += G3ww[k][l];
      }
    }
  }
*/

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      // interface particles
      int mi = get_matrix_index_from_local_index(i);
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (int kk = 0; kk < jnum; kk++) {
        int k = jlist[kk] & NEIGHMASK;

        if (mask[k] & groupbit) {
          // interface particles
          int mk = get_matrix_index_from_local_index(k);

          double a1 = em[mi] * (em[mk] - 1.0);
          double a2 = 1.0 - em[mi] - em[mk];
          // The first term (w/ G1ww) contributes the most to Rww
          // the second term (w/ G2ww) includes certain correction
          Rww[mi][mk] = a1 * G1ww[mi][mk] + a2 * G2ww[mi][mk];

          if (includingG3ww) {
            double temp = 0;
            for (int ll = 0; ll < jnum; ll++) {
              int l = jlist[ll] & NEIGHMASK;
              if (mask[l] & groupbit) {
                // interface particles
                int ml = get_matrix_index_from_local_index(l);
                temp += (ndotGww[ml][mi]) * G2ww[ml][mk] * area[l] * ed[l];
              }
            }
            G3ww[mi][mk] = temp;
            Rww[mi][mk] += G3ww[mi][mk];
          }
        }
      }

      if (includingG3ww) {
        double temp = 0;
        for (int ll = 0; ll < jnum; ll++) {
          int l = jlist[ll] & NEIGHMASK;
          if (mask[l] & groupbit) {
            // interface particles
            int ml = get_matrix_index_from_local_index(l);
            temp += (ndotGww[ml][mi]) * G2ww[ml][mi] * area[l] * ed[l];
          }
        }
        G3ww[mi][mi] = temp;
        Rww[mi][mi] += G3ww[mi][mi];
      }

      // including the diagonal term
      double a1 = em[i] * (em[i] - 1.0);
      double a2 = 1.0 - em[i] - em[i];
      // The first term (w/ G1ww) contributes the most to Rww
      // the second term (w/ G2ww) includes certain correction
      Rww[mi][mi] = a1 * G1ww[mi][mi] + a2 * G2ww[mi][mi];
      if (includingG3ww) Rww[mi][mi] += G3ww[mi][mi];
    }
  }

  #ifdef _POLARIZE_DEBUG
  if (comm->me == 0) {
    FILE* fp = fopen("Rww-qeq-cutoff.txt", "w");
    for (int i = 0; i < n_induced_charges; i++)
      fprintf(fp, "%d %g\n", i, Rww[i][i]);
    fclose(fp);
  }
  #endif
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::pre_force(int vflag)
{
//  setup_tags_local();

  // very important! change all ions q to real q to calculate the structure factor,
  // later will change to scaled q to calculate force.
  change_all_ions_q_to_real_q();

  // calculate structure factor and kspace terms first
//  ewaldDielectric->calculate_structure_factors();
//  ewaldDielectric->calculate_kspace_terms();

  // main calculation
  if (full) calculate_qiRqw_full();
  else calculate_qiRqw_cutoff();

  // compute induced charges
  compute_induced_charges();
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_qiRqw_full()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *q = avec->q_real;
  double *epsilon = avec->epsilon;
  double *area = avec->area;
  double **norm = avec->mu;
  double *ed = avec->ed;
  double *em = avec->em;

  // compute Green's functions
  // compute qRqw

  #pragma omp parallel
  {

  // calculate G1, gradient G1
  #pragma omp for
  for (int i = 0; i < n_ions; i++) {
    int ni = get_index_ions(i);
    double xtmp = x[ni][0];
    double ytmp = x[ni][1];
    double ztmp = x[ni][2];

    for (int k = 0; k < n_induced_charges; ++k) {
      int nk = get_index_interface(k);
      double delx = xtmp - x[nk][0];
      double dely = ytmp - x[nk][1];
      double delz = ztmp - x[nk][2];
      domain->minimum_image(delx,dely,delz);
      double r = sqrt(delx * delx + dely * dely + delz * delz);
      G1qw_real[i][k] = greens_real(r);
    }
  }

  // the following loop need the above results,
  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);
    double xtmp = x[nk][0];
    double ytmp = x[nk][1];
    double ztmp = x[nk][2];
    double tempndotG[3] = {0.0, 0.0, 0.0};
    double temp_sum1 = 0;
    for (int i = 0; i < n_ions; i++) {
      int ni = get_index_ions(i);
      // posvecqw[i][k]
      double delx = x[ni][0] - xtmp;
      double dely = x[ni][1] - ytmp;
      double delz = x[ni][2] - ztmp;
      domain->minimum_image(delx,dely,delz);
      calculate_grad_greens_real(gradG1wq_real[k][i], delx, dely, delz);
      MathExtra::scale3(-1.0, gradG1wq_real[k][i]);

      tempndotG[0] += gradG1wq_real[k][i][0] * (q[ni] / epsilon[ni]);
      tempndotG[1] += gradG1wq_real[k][i][1] * (q[ni] / epsilon[ni]);
      tempndotG[2] += gradG1wq_real[k][i][2] * (q[ni] / epsilon[ni]);
      temp_sum1 += G1qw_real[i][k] * q[ni] / epsilon[ni];
    }
    sum1G1qw_epsilon[k] = temp_sum1; // + ewaldDielectric->sum1G1qw_k_epsilon[k];
//      double temp_sum2 = 0;
//      for (int ic = 0; ic < 3; ++ic) {
//        temp_sum2 += norm[nk][ic] * (tempndotG[ic] + ewaldDielectric->sum2gradG1wq_k_epsilon[k][ic]);
//        temp_sum2 += norm[nk][ic] * (tempndotG[ic]);
//      }
    double temp_sum2 = MathExtra::dot3(norm[nk], tempndotG);
    sum2ndotGwq_epsilon[k] = temp_sum2;
  }

  // calculate G2, gradient G2
  #pragma omp for nowait
  for (int k = 0; k < n_induced_charges; ++k) {
    double tempwq = 0;
    for (int m = 0; m < n_induced_charges; ++m) {
      int nm = get_index_interface(m);
      tempwq += G1ww[k][m] * (sum2ndotGwq_epsilon[m]) * area[nm] * ed[nm];
    }
    sum2G2wq[k] = tempwq;
  }

  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    double temp = 0;
    for (int m = 0; m < n_induced_charges; ++m) {
      int nm = get_index_interface(m);
      temp += sum1G1qw_epsilon[m] * (ndotGww[m][k]) * area[nm] * ed[nm];
    }
    sum1G2qw[k] = temp;
  }

  // calculate G3, gradient G3
  #pragma omp for
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);

    double qiRwwVectorTemp1 = 0;
    for (int i = 0; i < n_ions; ++i) {
      int ni = get_index_ions(i);
      qiRwwVectorTemp1 += q[ni] * (1.0 - em[nk] / epsilon[ni]) * G1qw_real[i][k];
    }
//  qiRwwVectorTemp1 += ewaldDielectric->sum1G1qw_k[k] - em[nk] * ewaldDielectric->sum1G1qw_k_epsilon[k];

    double temp = 0;
    for (int m = 0; m < n_induced_charges; ++m) {
      int nm = get_index_interface(m);
      temp += sum2ndotGwq_epsilon[m] * G2ww[m][k] * area[nm] * ed[nm];
    }
    sum1G3qw[k] = temp;

    // qiRwwVectorTemp2 is a significant contribution, of which sum2G2wq is significant
    double qiRwwVectorTemp2 = (1.0 - 2.0 * em[nk]) * sum2G2wq[k] + sum1G2qw[k] + 2.0 * sum1G3qw[k];
    qiRqwVector[k] = qiRwwVectorTemp1 + qiRwwVectorTemp2;
  }

  } // end of parallel region
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_qiRqw_cutoff()
{
  int ii, i, mi, k, kk, mk, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, r;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double *q = avec->q_real;
  double *epsilon = avec->epsilon;
  double *area = avec->area;
  double **norm = avec->mu;
  double *ed = avec->ed;
  double *em = avec->em;

  // invoke full neighbor list (will copy or build if necessary)

  int inum,*ilist,*jlist,*numneigh,**firstneigh;

  inum = list->inum;               // number of entries in the neighbor list
  ilist = list->ilist;             // ilist[ii] gives the atom index of the entry ii in the list
  numneigh = list->numneigh;       // numneigh[i] gives the number of neighbors of local atom i
  firstneigh = list->firstneigh;   // firstneigh[i] gives the pointer to the neighbors of i

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) {
      // ion particles
      mi = get_matrix_index_from_local_index(i);
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (kk = 0; kk < jnum; kk++) {
        k = jlist[kk] & NEIGHMASK;
        if (mask[k] & groupbit) {
          // interface particles
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          domain->minimum_image(delx,dely,delz);
          r = sqrt(delx * delx + dely * dely + delz * delz);
          mk = get_matrix_index_from_local_index(k);
          G1qw_real[mi][mk] = greens_real(r);
        }
      }
    }
  }

  // the following loop need the above results,
  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];  // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      mk = get_matrix_index_from_local_index(k);
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
          // ions particles
          delx = x[i][0] - xtmp;
          dely = x[i][1] - ytmp;
          delz = x[i][2] - ztmp;
          domain->minimum_image(delx,dely,delz);

          mi = get_matrix_index_from_local_index(i);
          calculate_grad_greens_real(gradG1wq_real[mk][mi], delx, dely, delz);
          MathExtra::scale3(-1.0, gradG1wq_real[mk][mi]);

          tempndotG[0] += gradG1wq_real[mk][mi][0] * (q[i] / epsilon[i]);
          tempndotG[1] += gradG1wq_real[mk][mi][1] * (q[i] / epsilon[i]);
          tempndotG[2] += gradG1wq_real[mk][mi][2] * (q[i] / epsilon[i]);
          temp_sum1 += G1qw_real[mi][mk] * q[i] / epsilon[i];
        }
      }
      sum1G1qw_epsilon[mk] = temp_sum1;// + ewaldDielectric->sum1G1qw_k_epsilon[mk];
//      double temp_sum2 = 0;
//      for (int ic = 0; ic < 3; ++ic) {
//        temp_sum2 += norm[k][ic] * (tempndotG[ic] + ewaldDielectric->sum2gradG1wq_k_epsilon[mk][ic]);
//        temp_sum2 += norm[k][ic] * (tempndotG[ic]);
//      }
      double temp_sum2 = MathExtra::dot3(norm[k], tempndotG);
      sum2ndotGwq_epsilon[mk] = temp_sum2;
    }
  }

  // calculate G2, gradient G2
  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];  // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      mk = get_matrix_index_from_local_index(k);
      jlist = firstneigh[k];
      jnum = numneigh[k];

      double tempwq = 0;
      double temp = 0;
      for (ii = 0; ii < jnum; ii++) {
        i = jlist[ii] & NEIGHMASK;
        if (mask[i] & groupbit) {
          // interface particles
          mi = get_matrix_index_from_local_index(i);
          tempwq += G1ww[mk][mi] * (sum2ndotGwq_epsilon[mi]) * area[i] * ed[i];
          temp += sum1G1qw_epsilon[mi] * (ndotGww[mi][mk]) * area[i] * ed[i];
        }
      }

      // add the corresponding self terms
      tempwq += G1ww[mk][mk] * (sum2ndotGwq_epsilon[mk]) * area[k] * ed[k];
      temp += sum1G1qw_epsilon[mk] * (ndotGww[mk][mk]) * area[k] * ed[k];

      sum2G2wq[mk] = tempwq;
      sum1G2qw[mk] = temp;
    }
  }

  // calculate G3, gradient G3
  for (kk = 0; kk < inum; kk++) {
    k = ilist[kk];  // k is local index
    if (mask[k] & groupbit) {
      // interface particles
      mk = get_matrix_index_from_local_index(k);

      // qiRwwVectorTemp2 is a significant contribution, of which sum2G2wq is significant
      double qiRwwVectorTemp1 = 0;
      double qiRwwVectorTemp2 = (1.0 - 2.0 * em[k]) * sum2G2wq[mk];

      if (includingG3ww) {
        jlist = firstneigh[k];
        jnum = numneigh[k];
        double temp = 0;
        for (ii = 0; ii < jnum; ii++) {
          i = jlist[ii] & NEIGHMASK;
          if (mask[i] & groupbit) {
            // interface particles
            mi = get_matrix_index_from_local_index(i);
            temp += sum2ndotGwq_epsilon[mi] * G2ww[mi][mk] * area[i] * ed[i];
          } else {
            // ions particles
            mi = get_matrix_index_from_local_index(i);
            qiRwwVectorTemp1 += q[i] * (1.0 - em[k] / epsilon[i]) * G1qw_real[mi][mk];
          }
        }
        // qiRwwVectorTemp1 += ewaldDielectric->sum1G1qw_k[mk] - em[k] * ewaldDielectric->sum1G1qw_k_epsilon[mk];

        // add the corresponding self term
        temp += sum2ndotGwq_epsilon[mk] * G2ww[mk][mk] * area[k] * ed[k];
        sum1G3qw[mk] = temp;

        qiRwwVectorTemp2 += sum1G2qw[mk] + 2.0 * sum1G3qw[mk];
      }
      
//      double qiRwwVectorTemp2 = (1.0 - 2.0 * em[k]) * sum2G2wq[mk] + + sum1G2qw[mk] + 2.0 * sum1G3qw[mk];
      qiRqwVector[mk] = qiRwwVectorTemp1 + qiRwwVectorTemp2;
    }
  }
  #ifdef _POLARIZE_DEBUG
  if (comm->me == 0) {
    FILE* fp = fopen("qRqw-qeq.txt", "w");
    for (int i = 0; i < n_induced_charges; i++)
      fprintf(fp, "%d %g\n", i, qiRqwVector[i]);
    fclose(fp);
  }
  #endif
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::pre_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::compute_induced_charges()
{
  // inverse_matrix * qRqw, notice the w = - inverse_matrix * qRqw, do not forget the minus sign.
  calculate_matrix_multiply_vector(inverse_matrix, qiRqwVector,
    q_induced_charges, n_induced_charges);

  // map q induced charges back into atom, from global id tag = k to nlocal index of nk.
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);
    atom->q[nk] = -q_induced_charges[k] / (4*MY_PI);
  }

  // important! change to scaled q to calculate force.
  change_all_ions_q_to_scaled_q();
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::change_all_ions_q_to_real_q()
{
  for (int i = 0; i < n_ions; ++i) {
    int ni = get_index_ions(i);
    atom->q[ni] = avec->q_real[ni];
  }
  // communicate the induced charges
  // for single MPI, this updates the ghost atoms' charges
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::change_all_ions_q_to_scaled_q()
{
  for (int i = 0; i < n_ions; ++i) {
    int ni = get_index_ions(i);
    atom->q[ni] = avec->q_real[ni] / avec->epsilon[ni];
  }
  // communicate the induced charges
  // for single MPI, this updates the ghost atoms' charges
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

int FixQEqDielectric::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"full") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) full = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) full = 0;
      else error->all(FLERR,"Illegal fix_modify command for fix qeq/dielectric");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dielectrics") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix_modify command");
      double epsiloni=-1, areai=-1;
      double ediff = force->numeric(FLERR,arg[iarg+1]);
      double emean = force->numeric(FLERR,arg[iarg+2]);
      if (strcmp(arg[iarg+3],"NULL") != 0) epsiloni = force->numeric(FLERR,arg[iarg+3]);
      if (strcmp(arg[iarg+4],"NULL") != 0) areai = force->numeric(FLERR,arg[iarg+4]);

      set_dielectric_params(ediff, emean, epsiloni, areai);

      iarg += 5;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return iarg;
}

/* ---------------------------------------------------------------------- */

int FixQEqDielectric::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;
  for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEqDielectric::memory_usage()
{
  double bytes = 0;

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqDielectric::grow_arrays(int nmax)
{

}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqDielectric::copy_arrays(int i, int j, int delflag)
{

}

/* ----------------------------------------------------------------------
   set dielectric params for the atom in the group
------------------------------------------------------------------------- */

void FixQEqDielectric::set_dielectric_params(double ediff, double emean,
   double epsiloni, double areai)
{
  int i;
  double *area = avec->area;
  double *ed = avec->ed;
  double *em = avec->em;
  double *epsilon = avec->epsilon;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      ed[i] = ediff;
      em[i] = emean;
      if (areai > 0) area[i] = areai;
      if (epsiloni > 0) epsilon[i] = epsiloni;
    }
  }
}

/* ----------------------------------------------------------------------
  setup mapping from tags to indices in matrices/vectors
------------------------------------------------------------------------- */
void FixQEqDielectric::setup_tags()
{
  // This function should not be put into init(), because the indices would change, so put it in pre_force_setup().
  // TODO: this function should be called each time the local index changes.
//  if (atom->nlocal != atom->natoms)
//    error->all(FLERR, "Single MPI only");

  //TODO: suppose there are only 2 types of atoms, one for interface and one for ions.

  tagint* tag = atom->tag;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;

  int i_interface = 0;
  int i_ions = 0;
  int t;
  for (int i = 0; i < nlocal; ++i) {
    // interface atoms
    if (mask[i] & groupbit) {
      t = tag[i];
      tags_interface[i_interface] = t;
      if (t >= matrix_index_from_global_id.size())
        matrix_index_from_global_id.resize(t + 1);
      matrix_index_from_global_id[t] = i_interface;
      i_interface++;
    } else {
      // ions atoms
      t = tag[i];
      tags_ions[i_ions] = t;
      if (t >= matrix_index_from_global_id.size())
        matrix_index_from_global_id.resize(t + 1);
      matrix_index_from_global_id[t] = i_ions;
      i_ions++;
    }
  }

//  assert(i_interface == n_induced_charges);
//  assert(i_ions == n_ions);

  // ewald calculations need the information of all tags
//  ewaldDielectric->tags_interface = tags_interface;
//  ewaldDielectric->tags_ions = tags_ions;
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::print_all_properties()
{
  // set epsilon_ions, em, ed, and q of interfaces
  double *q_real = avec->q_real;
  double *q = atom->q;
  double *epsilon = avec->epsilon;
  double *area = avec->area;
  double *curvature = avec->curvature;
  double **norm = avec->mu;
  double *ed = avec->ed;
  double *em = avec->em;
  // debug file
  FILE *debugFile = fopen("_debug.log", "w");
  fprintf(debugFile, "============  Ions Info =============\n");
  for (int i = 0; i < n_ions; ++i) {
    int ni = get_index_ions(i);
    fprintf(debugFile, "id: %d, local index: %d, real q: %f, rescaled q: %f, epsilon: %f\n",
      atom->tag[ni], ni, q_real[ni], atom->q[ni], epsilon[ni]);
  }
  fprintf(debugFile, "============  Interface Info =============\n");
  for (int k = 0; k < n_induced_charges; ++k) {
    int nk = get_index_interface(k);
    fprintf(debugFile, "id: %d, local index: %d, xyz: %f %f %f, q: %f, scaled q: %f, area: %f, "
      "curvature: %f, em: %f, ed: %f, norm: %f %f %f\n",
      atom->tag[nk], nk, atom->x[nk][0], atom->x[nk][1], atom->x[nk][2],
      q_real[nk], atom->q[nk], area[nk], curvature[nk], em[nk], ed[nk],
      norm[nk][0], norm[nk][1], norm[nk][2]);
  }
  fclose(debugFile);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::setup_tags_local()
{
  // TODO: this function should be called each time the local index changes.
  //TODO: suppose there are only 2 types of atoms, one for interface and one for ions.
  tags_interface_local.clear();
  tags_ions_local.clear();
  int i_interface = 0;
  int i_ions = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    // interface atoms
    if ((atom->mask[i] & groupbit)) {
      tags_interface_local.push_back(atom->tag[i]);
      i_interface++;
    } else {
      // ions atoms
      tags_ions_local.push_back(atom->tag[i]);
      i_ions++;
    }
  }
  n_induced_charges_local = i_interface;
  n_ions_local = i_ions;

  // ewald calculations need the information of all tags
//  ewaldDielectric->tags_interface_local = &tags_interface_local;
//  ewaldDielectric->tags_ions_local = &tags_ions_local;
}

/* ---------------------------------------------------------------------- 
   real Green's function
------------------------------------------------------------------------ */

double FixQEqDielectric::greens_real(double r)
{
  return erfc(g_ewald * r) / r;
}

/* ---------------------------------------------------------------------- */

double FixQEqDielectric::grad_greens_real_factor(double r)
{
  double alpharij = g_ewald * r;
  double factor = erfc(alpharij) + 2.0 * alpharij / MY_PIS * exp(-(alpharij * alpharij));
  double r3 = r*r*r;
  return factor * (-1.0 / r3);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_grad_greens_real(double *vec, double dx, double dy, double dz)
{
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  double real = grad_greens_real_factor(r);
  vec[0] = real * dx;
  vec[1] = real * dy;
  vec[2] = real * dz;
}

/* ---------------------------------------------------------------------- */

double FixQEqDielectric::calculate_greens_ewald(double dx, double dy, double dz)
{
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  return greens_real(r);// + calculate_greens_ewald_reciprocal(dx, dy, dz);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_grad_greens_ewald(double *vec, double dx, double dy, double dz)
 {
  // real part of grad greens, must be first, because the vec is reset
  calculate_grad_greens_real(vec, dx, dy, dz);
  // kspace part, must be 2nd, because the vec is self additive.
//  calculate_grad_greens_ewald_reciprocal(vec, dx, dy, dz);
}

/* ---------------------------------------------------------------------- 
   given the matrix Rww, calculate inverse matrix (Rww + Rww^T)^(-1),
------------------------------------------------------------------------ */

void FixQEqDielectric::calculate_inverse_matrix(double **Rww, double **inverse_matrix, int M)
{
  // convel rt into gsmatrix
  gsl_matrix *matrixRww = gsl_matrix_alloc(M, M);
  for (int k = 0; k < M; ++k) {
    for (int l = 0; l < M; ++l) {
      gsl_matrix_set(matrixRww, k, l, Rww[k][l]);
    }
  }

  // QT: = Rww + Rww ^ T   Matrix
  gsl_matrix *QT = gsl_matrix_alloc(M, M);
  gsl_matrix_transpose_memcpy(QT, matrixRww);
  gsl_matrix_add(QT, matrixRww);


  int sign;
  gsl_permutation *p = gsl_permutation_alloc(M);
  gsl_linalg_LU_decomp(QT, p, &sign);
  gsl_linalg_LU_invert(QT, p, matrixRww); // here matrixRww to store inverse_matrix

  // convert back
  for (int k = 0; k < M; ++k) {
    for (int l = 0; l < M; ++l) {
      inverse_matrix[k][l] = gsl_matrix_get(matrixRww, k, l);
    }
  }

  gsl_matrix_free(QT);
  gsl_permutation_free(p);
  gsl_matrix_free(matrixRww);
}

/* ---------------------------------------------------------------------- */

void FixQEqDielectric::calculate_matrix_multiply_vector(double **matrix,
  double *in_vec, double *out_vec, int M)
{
  #pragma parallel omp for
  for (int k = 0; k < M; ++k) {
    double temp = 0.0;
    for (int l = 0; l < M; ++l) {
      temp += matrix[k][l] * in_vec[l];
    }
    out_vec[k] = temp;
  }
}

/* ---------------------------------------------------------------------- */

double FixQEqDielectric::calculate_greens_ewald_self_vertex(double area)
{
  // this term is very important, cannot be set to zero. see <Bugs in lammps implementation>.pptx
  double corr = 2.0 * MY_PIS / sqrt(area);
  double self_energy = -2.0 * g_ewald / MY_PIS;
  return corr + self_energy;// + calculate_greens_ewald_reciprocal(0, 0, 0);
}

/* ---------------------------------------------------------------------- */

double FixQEqDielectric::calculate_ndotgreens_ewald_self_vertex(double area, double curvature)
{
  // this term is important, cannot be set to zero. see <Bugs in lammps implementation>.pptx
  // curvature = 1 / R, minus if norm is inverse of R to center.
  // Honghao Li's result, the same with Erik's paper, J Chem Phys 140 064903 (2014)
  return curvature * MY_PIS / sqrt(area);
}

