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
   Contributing author:
   Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

   Hybrid & sub-group capabilities added by Ray Shan (Materials Design)

   OpenMP based threading support for fix qeq/reax/omp added
   by Hasan Metin Aktulga (MSU), Chris Knight (ALCF), Paul Coffman (ALCF),
   Kurt O'Hearn (MSU), Ray Shan (Materials Design), Wei Jiang (ALCF)

   Integration of the pair_style reax/c/omp into the User-OMP package
   by Axel Kohlmeyer (Temple U.)

   Please cite the related publication:
   H. M. Aktulga, C. Knight, P. Coffman, K. A. O'Hearn, T. R. Shan,
   W. Jiang, "Optimizing the performance of reactive molecular dynamics
   simulations for multi-core architectures", International Journal of
   High Performance Computing Applications, to appear.
 ------------------------------------------------------------------------- */

#include "fix_qeq_reax_omp.h"
#include <mpi.h>
#include <cmath>
#include "pair_reaxc.h"
#include "atom.h"
#include "comm.h"
#include "neigh_list.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "reaxc_defs.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

#define EV_TO_KCAL_PER_MOL 14.4
//#define DANGER_ZONE     0.95
//#define LOOSE_ZONE      0.7
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN_NBRS 100

/* ---------------------------------------------------------------------- */

FixQEqReaxOMP::FixQEqReaxOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEqReax(lmp, narg, arg)
{
  if (narg<8 || narg>9) error->all(FLERR,"Illegal fix qeq/reax/omp command");

  b_temp = NULL;

  // ASPC: Kolafa, J. Comp. Chem., 25(3), 335 (2003)
  do_aspc = 0;
  aspc_order = 1;
  // Must be consistent with nprev to store history: nprev = aspc_order + 2
  aspc_order_max = nprev - 2;
  aspc_omega = 0.0;
  aspc_b = NULL;
}

FixQEqReaxOMP::~FixQEqReaxOMP()
{
  memory->destroy(b_temp);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::post_constructor()
{
  pertype_parameters(pertype_option);
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::allocate_storage()
{
  FixQEqReax::allocate_storage();

  // dual CG support
  int size = nmax;
  if (dual_enabled) size*= 2;
  memory->create(b_temp, comm->nthreads, size, "qeq/reax/omp:b_temp");
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::deallocate_storage()
{
  memory->destroy(b_temp);

  FixQEqReax::deallocate_storage();
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::init()
{
  FixQEqReax::init();

  // APSC setup
  if (do_aspc) {
    memory->create(aspc_b, aspc_order_max+2, "qeq/reax/aspc_b");

    // Calculate damping factor
    double o = double(aspc_order);
    aspc_omega = (o+2.0) / (2*o+3.0);

    // Calculate B coefficients
    double c = (4.0 * o + 6.0) / (o + 3.0);
    aspc_b[0] = c;

    double n =  1.0;
    double d =  4.0;
    double s = -1.0;
    double f =  2.0;

    for (int i=1; i<aspc_order_max+2; i++) {
      c*= (o + n) / (o + d);
      aspc_b[i] = s * f * c;

      s *= -1.0;
      f += 1.0;
      n -= 1.0;
      d += 1.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::compute_H()
{
  int inum, *ilist, *numneigh, **firstneigh;
  double SMALL = 0.0001;

  int *type = atom->type;
  tagint * tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;

  if (reaxc) {
    inum = reaxc->list->inum;
    ilist = reaxc->list->ilist;
    numneigh = reaxc->list->numneigh;
    firstneigh = reaxc->list->firstneigh;
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }
  int ai, num_nbrs;

  // sumscan of the number of neighbors per atom to determine the offsets
  // most likely, we are overallocating. desirable to work on this part
  // to reduce the memory footprint of the far_nbrs list.

  num_nbrs = 0;

  for (int itr_i = 0; itr_i < inum; ++itr_i) {
    ai = ilist[itr_i];
    H.firstnbr[ai] = num_nbrs;
    num_nbrs += numneigh[ai];
  }
  m_fill = num_nbrs;

  // fill in the H matrix

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i, j, ii, jj, mfill, jnum, flag;
    int *jlist;
    double dx, dy, dz, r_sqr;

    mfill = 0;

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        jlist = firstneigh[i];
        jnum = numneigh[i];
        mfill = H.firstnbr[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];

          dx = x[j][0] - x[i][0];
          dy = x[j][1] - x[i][1];
          dz = x[j][2] - x[i][2];
          r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

          flag = 0;
          if (r_sqr <= SQR(swb)) {
            if (j < n) flag = 1;
            else if (tag[i] < tag[j]) flag = 1;
            else if (tag[i] == tag[j]) {
              if (dz > SMALL) flag = 1;
              else if (fabs(dz) < SMALL) {
                if (dy > SMALL) flag = 1;
                else if (fabs(dy) < SMALL && dx > SMALL) flag = 1;
              }
            }
          }

          if (flag) {
            H.jlist[mfill] = j;
            H.val[mfill] = calculate_H( sqrt(r_sqr), shld[type[i]][type[j]] );
            mfill++;
          }
        }

        H.numnbrs[i] = mfill - H.firstnbr[i];
      }
    }

    if (mfill >= H.m) {
      char str[128];
      sprintf(str,"H matrix size has been exceeded: mfill=%d H.m=%d\n",
              mfill, H.m);
      error->warning(FLERR,str);
      error->all(FLERR,"Fix qeq/reax/omp has insufficient QEq matrix size");
    }
  } // omp

}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::init_storage()
{
  int NN;

  if (reaxc) NN = reaxc->list->inum + reaxc->list->gnum;
  else NN = list->inum + list->gnum;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < NN; i++) {
    Hdia_inv[i] = 1. / eta[atom->type[i]];
    b_s[i] = -chi[atom->type[i]];
    b_t[i] = -1.0;
    b_prc[i] = 0;
    b_prm[i] = 0;
    s[i] = t[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::pre_force(int /* vflag */)
{

#ifdef OMP_TIMING
  double endTimeBase, startTimeBase, funcstartTimeBase;
  funcstartTimeBase = MPI_Wtime();
#endif

  double t_start, t_end;

  if (update->ntimestep % nevery) return;
  if (comm->me == 0) t_start = MPI_Wtime();

  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;

  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) reallocate_storage();
  if (n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

#ifdef OMP_TIMING
  startTimeBase = MPI_Wtime();
#endif

  init_matvec();

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEINITMVINDEX] += (endTimeBase-startTimeBase);
  startTimeBase = endTimeBase;
#endif

  if (dual_enabled) {
    matvecs = dual_CG(b_s, b_t, s, t); // OMP_TIMING inside dual_CG
  } else {
    matvecs_s = CG(b_s, s);     // CG on s - parallel

#ifdef OMP_TIMING
    endTimeBase = MPI_Wtime();
    ompTimingData[COMPUTECG1INDEX] += (endTimeBase-startTimeBase);
    ompTimingCount[COMPUTECG1INDEX]++;
    ompTimingCGCount[COMPUTECG1INDEX]+= matvecs_s;
    startTimeBase = endTimeBase;
#endif

    matvecs_t = CG(b_t, t);     // CG on t - parallel

#ifdef OMP_TIMING
    endTimeBase = MPI_Wtime();
    ompTimingData[COMPUTECG2INDEX] += (endTimeBase-startTimeBase);
    ompTimingCount[COMPUTECG2INDEX]++;
    ompTimingCGCount[COMPUTECG2INDEX]+= matvecs_t;
    startTimeBase = endTimeBase;
#endif

  } // if (dual_enabled)

#ifdef OMP_TIMING
  startTimeBase = MPI_Wtime();
#endif

  calculate_Q();

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTECALCQINDEX] += (endTimeBase-startTimeBase);
#endif

  if (comm->me == 0) {
    t_end = MPI_Wtime();
    qeq_time = t_end - t_start;
  }

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEQEQINDEX] += (endTimeBase-funcstartTimeBase);
#endif
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::init_matvec()
{
#ifdef OMP_TIMING
  long endTimeBase, startTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  /* fill-in H matrix */
  compute_H();

  int nn,i;
  int *ilist;

  if (reaxc) {
    nn = reaxc->list->inum;
    ilist = reaxc->list->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }

  // Should really be more careful with initialization and first (aspc_order+2) MD steps
  if (do_aspc) {

    double m_aspc_omega = 1.0 - aspc_omega;
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) private(i)
#endif
    for (int ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {

        /* init pre-conditioner for H and init solution vectors */
        Hdia_inv[i] = 1. / eta[ atom->type[i] ];
        b_s[i]      = -chi[ atom->type[i] ];
        b_t[i]      = -1.0;

        // Predictor Step
        double tp = 0.0;
        double sp = 0.0;
        for (int j=0; j<aspc_order+2; j++) {
          tp+= aspc_b[j] * t_hist[i][j];
          sp+= aspc_b[j] * s_hist[i][j];
        }

        // Corrector Step
        t[i] = aspc_omega * t_hist[i][0] + m_aspc_omega * tp;
        s[i] = aspc_omega * s_hist[i][0] + m_aspc_omega * sp;
      }
    }

  } else {

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) private(i)
#endif
    for (int ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {

        /* init pre-conditioner for H and init solution vectors */
        Hdia_inv[i] = 1. / eta[ atom->type[i] ];
        b_s[i]      = -chi[ atom->type[i] ];
        b_t[i]      = -1.0;

        /* linear extrapolation for s & t from previous solutions */
        //s[i] = 2 * s_hist[i][0] - s_hist[i][1];
        //t[i] = 2 * t_hist[i][0] - t_hist[i][1];

        /* quadratic extrapolation for s & t from previous solutions */
        //s[i] = s_hist[i][2] + 3 * ( s_hist[i][0] - s_hist[i][1] );
        t[i] = t_hist[i][2] + 3 * ( t_hist[i][0] - t_hist[i][1] );

        /* cubic extrapolation for s & t from previous solutions */
        s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
        //t[i] = 4*(t_hist[i][0]+t_hist[i][2])-(6*t_hist[i][1]+t_hist[i][3]);
      }
    }
  }

  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector( t );

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEMVCOMPINDEX] += (long) (endTimeBase-startTimeBase);
#endif
}

/* ---------------------------------------------------------------------- */

int FixQEqReaxOMP::CG( double *b, double *x)
{
  int  i, ii, imax;
  double alpha, beta, b_norm;
  double sig_old, sig_new;

  double my_buf[2], buf[2];

  int nn, jj;
  int *ilist;
  if (reaxc) {
    nn = reaxc->list->inum;
    ilist = reaxc->list->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }

  imax = 200;

  pack_flag = 1;
  sparse_matvec( &H, x, q );
  comm->reverse_comm_fix( this); //Coll_Vector( q );

  double tmp1, tmp2;
  tmp1 = tmp2 = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) private(i) reduction(+:tmp1,tmp2)
#endif
  for (jj = 0; jj < nn; ++jj) {
    i = ilist[jj];
    if (atom->mask[i] & groupbit) {
      r[i] = b[i] - q[i];
      d[i] = r[i] * Hdia_inv[i]; //pre-condition

      tmp1 += b[i] * b[i];
      tmp2 += r[i] * d[i];
    }
  }

  my_buf[0] = tmp1;
  my_buf[1] = tmp2;

  MPI_Allreduce(&my_buf, &buf, 2, MPI_DOUBLE, MPI_SUM, world);

  b_norm = sqrt(buf[0]);
  sig_new = buf[1];

  for (i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i) {
    comm->forward_comm_fix(this); //Dist_vector( d );
    sparse_matvec( &H, d, q );
    comm->reverse_comm_fix(this); //Coll_vector( q );

    tmp1 = 0.0;
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) private(ii) reduction(+:tmp1)
#endif
      for (jj = 0; jj < nn; jj++) {
        ii = ilist[jj];
        if (atom->mask[ii] & groupbit) tmp1 += d[ii] * q[ii];
      }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp master
#endif
      {
        MPI_Allreduce(&tmp1, &tmp2, 1, MPI_DOUBLE, MPI_SUM, world);

        alpha = sig_new / tmp2;
        tmp1 = 0.0;
      }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50) private(ii) reduction(+:tmp1)
#endif
      for (jj = 0; jj < nn; jj++) {
        ii = ilist[jj];
        if (atom->mask[ii] & groupbit) {
          x[ii] += alpha * d[ii];
          r[ii] -= alpha * q[ii];

          // pre-conditioning
          p[ii] = r[ii] * Hdia_inv[ii];
          tmp1 += r[ii] * p[ii];
        }
      }
    } // omp parallel

    sig_old = sig_new;

    MPI_Allreduce(&tmp1, &tmp2, 1, MPI_DOUBLE, MPI_SUM, world);

    sig_new = tmp2;
    beta = sig_new / sig_old;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) private(ii)
#endif
    for (jj = 0; jj < nn; jj++) {
      ii = ilist[jj];
      if (atom->mask[ii] & groupbit) d[ii] = p[ii] + beta * d[ii];
    }
  }

  if (i >= imax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq/reax CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }

  return i;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::sparse_matvec( sparse_matrix *A, double *x, double *b)
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i, j, itr_j;
    int nn, NN, ii;
    int *ilist;
    int nthreads = comm->nthreads;
#if defined(_OPENMP)
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    if (reaxc) {
      nn = reaxc->list->inum;
      NN = reaxc->list->inum + reaxc->list->gnum;
      ilist = reaxc->list->ilist;
    } else {
      nn = list->inum;
      NN = list->inum + list->gnum;
      ilist = list->ilist;
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) b[i] = eta[ atom->type[i] ] * x[i];
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = nn; ii < NN; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) b[i] = 0;
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i)
      for (int t=0; t<nthreads; t++) b_temp[t][i] = 0.0;

    // Wait for b accumulated and b_temp zeroed.
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
          j = A->jlist[itr_j];
          b[i] += A->val[itr_j] * x[j];

          b_temp[tid][j] += A->val[itr_j] * x[i];
        }
      }
    }

    // Wait till b_temp accumulated
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i)
      for (int t = 0; t < nthreads; ++t) b[i] += b_temp[t][i];

  } //end omp parallel
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::calculate_Q()
{
  int i;
  double *q = atom->q;

  int nn;
  int *ilist;

  if (reaxc) {
    nn = reaxc->list->inum;
    ilist = reaxc->list->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }

  double tmp1, tmp2;
  tmp1 = tmp2 = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) private(i) reduction(+:tmp1,tmp2)
#endif
  for (int ii = 0; ii < nn; ii++) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      tmp1 += s[i];
      tmp2 += t[i];
    }
  }

  double my_buf[2], buf[2];
  buf[0] = 0.0;
  buf[1] = 0.0;

  my_buf[0] = tmp1;
  my_buf[1] = tmp2;

  MPI_Allreduce(&my_buf,&buf,2,MPI_DOUBLE,MPI_SUM,world);

  double u = buf[0] / buf[1];

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) private(i)
#endif
  for (int ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      // backup s & t
      for (int k = nprev-1; k > 0; --k) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm_fix( this); //Dist_vector( atom->q );
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::vector_sum( double* dest, double c, double* v,
                                double d, double* y, int k)
{
  int i;
  int *ilist;

  if (reaxc) ilist = reaxc->list->ilist;
  else ilist = list->ilist;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) private(i)
#endif
  for (int ii=0; ii<k; ii++) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) dest[i] = c * v[i] + d * y[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::vector_add( double* dest, double c, double* v, int k)
{
  int i;
  int *ilist;

  if (reaxc) ilist = reaxc->list->ilist;
  else ilist = list->ilist;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) private(i)
#endif
  for (int ii=0; ii<k; ii++) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) dest[i] += c * v[i];
  }
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
/* dual CG support                                                        */
/* ---------------------------------------------------------------------- */

int FixQEqReaxOMP::dual_CG( double *b1, double *b2, double *x1, double *x2)
{

#ifdef OMP_TIMING
  double endTimeBase, startTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  int  i, imax;
  double alpha_s, alpha_t, beta_s, beta_t, b_norm_s, b_norm_t;
  double sig_old_s, sig_old_t, sig_new_s, sig_new_t;

  double my_buf[4], buf[4];

  int nn, ii, jj;
  int *ilist;
  if (reaxc) {
    nn = reaxc->list->inum;
    ilist = reaxc->list->ilist;
  } else {
    nn = list->inum;
    ilist = list->ilist;
  }

  imax = 200;

  pack_flag = 5; // forward 2x d and reverse 2x q
  dual_sparse_matvec( &H, x1, x2, q );
  comm->reverse_comm_fix( this); //Coll_Vector( q );

  double tmp1, tmp2, tmp3, tmp4;
  tmp1 = tmp2 = tmp3 = tmp4 = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) private(i) reduction(+:tmp1,tmp2,tmp3,tmp4)
#endif
  for (jj = 0; jj < nn; ++jj) {
    i = ilist[jj];
    if (atom->mask[i] & groupbit) {
      int indxI = 2 * i;
      r[indxI  ] = b1[i] - q[indxI  ];
      r[indxI+1] = b2[i] - q[indxI+1];

      d[indxI  ] = r[indxI  ] * Hdia_inv[i]; //pre-condition
      d[indxI+1] = r[indxI+1] * Hdia_inv[i];

      tmp1 += b1[i] * b1[i];
      tmp2 += b2[i] * b2[i];

      tmp3 += r[indxI  ] * d[indxI  ];
      tmp4 += r[indxI+1] * d[indxI+1];
    }
  }

  my_buf[0] = tmp1;
  my_buf[1] = tmp2;
  my_buf[2] = tmp3;
  my_buf[3] = tmp4;

  MPI_Allreduce(&my_buf, &buf, 4, MPI_DOUBLE, MPI_SUM, world);

  b_norm_s = sqrt(buf[0]);
  b_norm_t = sqrt(buf[1]);

  sig_new_s = buf[2];
  sig_new_t = buf[3];

  for (i = 1; i < imax; ++i) {
    comm->forward_comm_fix(this); //Dist_vector( d );
    dual_sparse_matvec( &H, d, q );
    comm->reverse_comm_fix(this); //Coll_vector( q );

    tmp1 = tmp2 = 0.0;
#if defined(_OPENMP)
#pragma omp parallel
#endif
    {

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) private(ii) reduction(+:tmp1,tmp2)
#endif
      for (jj = 0; jj < nn; jj++) {
        ii = ilist[jj];
        if (atom->mask[ii] & groupbit) {
          int indxI = 2 * ii;
          tmp1 += d[indxI  ] * q[indxI  ];
          tmp2 += d[indxI+1] * q[indxI+1];
        }
      }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp master
#endif
      {
        my_buf[0] = tmp1;
        my_buf[1] = tmp2;

        MPI_Allreduce(&my_buf, &buf, 2, MPI_DOUBLE, MPI_SUM, world);

        alpha_s = sig_new_s / buf[0];
        alpha_t = sig_new_t / buf[1];

        tmp1 = tmp2 = 0.0;
      }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50) private(ii) reduction(+:tmp1,tmp2)
#endif
      for (jj = 0; jj < nn; jj++) {
        ii = ilist[jj];
        if (atom->mask[ii] & groupbit) {
          int indxI = 2 * ii;
          x1[ii] += alpha_s * d[indxI  ];
          x2[ii] += alpha_t * d[indxI+1];

          r[indxI  ] -= alpha_s * q[indxI  ];
          r[indxI+1] -= alpha_t * q[indxI+1];

          // pre-conditioning
          p[indxI  ] = r[indxI  ] * Hdia_inv[ii];
          p[indxI+1] = r[indxI+1] * Hdia_inv[ii];

          tmp1 += r[indxI  ] * p[indxI  ];
          tmp2 += r[indxI+1] * p[indxI+1];
        }
      }
    } // omp parallel

    my_buf[0] = tmp1;
    my_buf[1] = tmp2;

    sig_old_s = sig_new_s;
    sig_old_t = sig_new_t;

    MPI_Allreduce(&my_buf, &buf, 2, MPI_DOUBLE, MPI_SUM, world);

    sig_new_s = buf[0];
    sig_new_t = buf[1];

    if (sqrt(sig_new_s)/b_norm_s <= tolerance
        || sqrt(sig_new_t)/b_norm_t <= tolerance) break;

    beta_s = sig_new_s / sig_old_s;
    beta_t = sig_new_t / sig_old_t;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) private(ii)
#endif
    for (jj = 0; jj < nn; jj++) {
      ii = ilist[jj];
      if (atom->mask[ii] & groupbit) {
        int indxI = 2 * ii;

        d[indxI  ] = p[indxI  ] + beta_s * d[indxI  ];
        d[indxI+1] = p[indxI+1] + beta_t * d[indxI+1];
      }
    }
  }

  i++;
  matvecs_s = matvecs_t = i; // The plus one makes consistent with count from CG()
  matvecs = i;

  // Timing info for iterating s&t together
#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTECG1INDEX] += (endTimeBase-startTimeBase);
  ompTimingCount[COMPUTECG1INDEX]++;
  ompTimingCGCount[COMPUTECG1INDEX]+= i;
  startTimeBase = endTimeBase;
#endif

  // If necessary, converge other system
  if (sqrt(sig_new_s)/b_norm_s > tolerance) {
    pack_flag = 2;
    comm->forward_comm_fix(this); // x1 => s

    i+= CG(b1, x1);
    matvecs_s = i;
  }
  else if (sqrt(sig_new_t)/b_norm_t > tolerance) {
    pack_flag = 3;
    comm->forward_comm_fix(this); // x2 => t

    i+= CG(b2, x2);
    matvecs_t = i;
  }

  // Timing info for remainder of s or t
#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTECG2INDEX] += (endTimeBase-startTimeBase);
  ompTimingCount[COMPUTECG2INDEX]++;
  ompTimingCGCount[COMPUTECG2INDEX]+= i - matvecs;
  startTimeBase = endTimeBase;
#endif

  if ( i >= imax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq/reax CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }

  return i;
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::dual_sparse_matvec( sparse_matrix *A, double *x1, double *x2, double *b)
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i, j, itr_j;
    int nn, NN, ii;
    int *ilist;
    int indxI, indxJ;

    int nthreads = comm->nthreads;
#if defined(_OPENMP)
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    if (reaxc) {
      nn = reaxc->list->inum;
      NN = reaxc->list->inum + reaxc->list->gnum;
      ilist = reaxc->list->ilist;
    } else {
      nn = list->inum;
      NN = list->inum + list->gnum;
      ilist = list->ilist;
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        b[indxI  ] = eta[ atom->type[i] ] * x1[i];
        b[indxI+1] = eta[ atom->type[i] ] * x2[i];
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = nn; ii < NN; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        b[indxI]   = 0;
        b[indxI+1] = 0;
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i) {
      indxI = 2 * i;
      for (int t=0; t<nthreads; t++) {
        b_temp[t][indxI  ] = 0.0;
        b_temp[t][indxI+1] = 0.0;
      }
    }

    // Wait for b accumulated and b_temp zeroed
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
          j = A->jlist[itr_j];
          indxJ = 2 * j;
          b[indxI  ] += A->val[itr_j] * x1[j];
          b[indxI+1] += A->val[itr_j] * x2[j];

          b_temp[tid][indxJ  ] += A->val[itr_j] * x1[i];
          b_temp[tid][indxJ+1] += A->val[itr_j] * x2[i];
        }
      }
    }

    // Wait till b_temp accumulated
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i) {
      indxI = 2 * i;
      for (int t = 0; t < nthreads; ++t) {
        b[indxI  ] += b_temp[t][indxI  ];
        b[indxI+1] += b_temp[t][indxI+1];
      }
    }

  } // omp parallel
}

/* ---------------------------------------------------------------------- */

void FixQEqReaxOMP::dual_sparse_matvec( sparse_matrix *A, double *x, double *b )
{
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i, j, itr_j;
    int nn, NN, ii;
    int *ilist;
    int indxI, indxJ;

    int nthreads = comm->nthreads;
#if defined(_OPENMP)
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    if (reaxc) {
      nn = reaxc->list->inum;
      NN = reaxc->list->inum + reaxc->list->gnum;
      ilist = reaxc->list->ilist;
    } else {
      nn = list->inum;
      NN = list->inum + list->gnum;
      ilist = list->ilist;
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        b[indxI  ] = eta[ atom->type[i] ] * x[indxI  ];
        b[indxI+1] = eta[ atom->type[i] ] * x[indxI+1];
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = nn; ii < NN; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        b[indxI]   = 0;
        b[indxI+1] = 0;
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i) {
      indxI = 2 * i;
      for (int t=0; t<nthreads; t++) {
        b_temp[t][indxI  ] = 0.0;
        b_temp[t][indxI+1] = 0.0;
      }
    }

    // Wait for b accumulated and b_temp zeroed
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (ii = 0; ii < nn; ++ii) {
      i = ilist[ii];
      if (atom->mask[i] & groupbit) {
        indxI = 2 * i;
        for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
          j = A->jlist[itr_j];
          indxJ = 2 * j;
          b[indxI  ] += A->val[itr_j] * x[indxJ  ];
          b[indxI+1] += A->val[itr_j] * x[indxJ+1];

          b_temp[tid][indxJ  ] += A->val[itr_j] * x[indxI  ];
          b_temp[tid][indxJ+1] += A->val[itr_j] * x[indxI+1];
        }
      }
    }

    // Wait till b_temp accumulated
#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (i = 0; i < NN; ++i) {
      indxI = 2 * i;
      for (int t = 0; t < nthreads; ++t) {
        b[indxI  ] += b_temp[t][indxI  ];
        b[indxI+1] += b_temp[t][indxI+1];
      }
    }
  } // omp parallel
}
