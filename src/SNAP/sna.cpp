/* ----------------------------------------------------------------------
   Authors: Aidan Thompson and Christian Trott, Sandia National Labs, 2012
   Property of Sandia National Labs: Not for External Distribution
------------------------------------------------------------------------- */

#include "sna.h"
#include "math.h"
#include "math_const.h"
#include "math_extra.h"
#include "string.h"
#include "stdlib.h"
#include "openmp_snap.h"

#include "memory.h"
#include "error.h"
#include "comm.h"
#include "atom.h"

using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;

/* ----------------------------------------------------------------------

   this implementation is based on the method outlined
   in Bartok[1], using formulae from VMK[2].

   for the Clebsch-Gordan coefficients, we
   convert the VMK half-integral labels
   a, b, c, alpha, beta, gamma
   to array offsets j1, j2, j, m1, m2, m
   using the following relations:

   j1 = 2*a
   j2 = 2*b
   j =  2*c

   m1 = alpha+a      2*alpha = 2*m1 - j1
   m2 = beta+b    or 2*beta = 2*m2 - j2
   m =  gamma+c      2*gamma = 2*m - j

   in this way: 

   -a <= alpha <= a 
   -b <= beta <= b 
   -c <= gamma <= c 

   becomes:

   0 <= m1 <= j1
   0 <= m2 <= j2
   0 <= m <= j

   and the requirement that
   a+b+c be integral implies that
   j1+j2+j must be even. 
   The requirement that:
   
   gamma = alpha+beta 

   becomes:

   2*m - j = 2*m1 - j1 + 2*m2 - j2

   Similarly, for the Wigner U-functions U(J,m,m') we
   convert the half-integral labels J,m,m' to
   array offsets j,ma,mb:

   j = 2*J
   ma = J+m
   mb = J+m'

   so that:

   0 <= j <= 2*Jmax
   0 <= ma, mb <= j.

   For the bispectrum components B(J1,J2,J) we convert to:
   
   j1 = 2*J1
   j2 = 2*J2
   j = 2*J

   and the requirement:

   |J1-J2| <= J <= J1+J2, for j1+j2+j integral

   becomes:

   |j1-j2| <= j <= j1+j2, for j1+j2+j even integer

   or 
   
   j = |j1-j2|, |j1-j2|+2,...,j1+j2-2,j1+j2

   [1] Albert Bartok-Partay, "Gaussian Approximation..."
   Doctoral Thesis, Cambrindge University, (2009)

   [2] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
   "Quantum Theory of Angular Momentum," World Scientific (1988)

------------------------------------------------------------------------- */

SNA::SNA(LAMMPS* lmp, double rfac0_in,
         int twojmax_in, int diagonalstyle_in, int use_shared_arrays_in,
         double rmin0_in, int switch_flag_in) : Pointers(lmp)
{
  wself = 1.0;

  use_shared_arrays = use_shared_arrays_in;
  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;

  twojmax = twojmax_in;
  diagonalstyle = diagonalstyle_in;

  ncoeff = compute_ncoeff();

  create_twojmax_arrays();

  bvec = NULL;
  dbvec = NULL;
  memory->create(bvec, ncoeff, "pair:bvec");
  memory->create(dbvec, ncoeff, 3, "pair:dbvec");
  rij = NULL;
  inside = NULL;
  wj = NULL;
  rcutij = NULL;
  nmax = 0;
  idxj = NULL;

  timers = new double[20];

  for(int i = 0; i < 20; i++) timers[i] = 0;

  print = 0;
  counter = 0;

  build_indexlist();

}

/* ---------------------------------------------------------------------- */

SNA::~SNA()
{
  if(!use_shared_arrays) {
    destroy_twojmax_arrays();
    memory->destroy(rij);
    memory->destroy(inside);
    memory->destroy(wj);
    memory->destroy(rcutij);
    memory->destroy(bvec);
    memory->destroy(dbvec);
  }
  delete[] idxj;
}

void SNA::build_indexlist()
{
  if(diagonalstyle == 0) {
    int idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
          idxj_count++;

    // indexList can be changed here

    idxj = new SNA_LOOPINDICES[idxj_count];
    idxj_max = idxj_count;

    idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
          idxj[idxj_count].j1 = j1;
          idxj[idxj_count].j2 = j2;
          idxj[idxj_count].j = j;
          idxj_count++;
        }
  }

  if(diagonalstyle == 1) {
    int idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j = 0; j <= MIN(twojmax, 2 * j1); j += 2) {
        idxj_count++;
      }

    // indexList can be changed here

    idxj = new SNA_LOOPINDICES[idxj_count];
    idxj_max = idxj_count;

    idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j = 0; j <= MIN(twojmax, 2 * j1); j += 2) {
        idxj[idxj_count].j1 = j1;
        idxj[idxj_count].j2 = j1;
        idxj[idxj_count].j = j;
        idxj_count++;
      }
  }

  if(diagonalstyle == 2) {
    int idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++) {
      idxj_count++;
    }

    // indexList can be changed here

    idxj = new SNA_LOOPINDICES[idxj_count];
    idxj_max = idxj_count;

    idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++) {
      idxj[idxj_count].j1 = j1;
      idxj[idxj_count].j2 = j1;
      idxj[idxj_count].j = j1;
      idxj_count++;
    }
  }

  if(diagonalstyle == 3) {
    int idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) idxj_count++;

    // indexList can be changed here

    idxj = new SNA_LOOPINDICES[idxj_count];
    idxj_max = idxj_count;

    idxj_count = 0;

    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
	  if (j >= j1) {
	    idxj[idxj_count].j1 = j1;
	    idxj[idxj_count].j2 = j2;
	    idxj[idxj_count].j = j;
	    idxj_count++;
	  }
  }

}
/* ---------------------------------------------------------------------- */

void SNA::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
  // if(comm->me == 0) {
  //   if(screen) print_clebsch_gordan(screen);
  //   if(logfile) print_clebsch_gordan(logfile);
  // }
}


void SNA::grow_rij(int newnmax)
{
  if(newnmax <= nmax) return;

  nmax = newnmax;

  if(!use_shared_arrays) {
    memory->destroy(rij);
    memory->destroy(inside);
    memory->destroy(wj);
    memory->destroy(rcutij);
    memory->create(rij, nmax, 3, "pair:rij");
    memory->create(inside, nmax, "pair:inside");
    memory->create(wj, nmax, "pair:wj"); 
    memory->create(rcutij, nmax, "pair:rcutij");
 }
}
/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

void SNA::compute_ui(int jnum)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma 
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  zero_uarraytot();
  addself_uarraytot(wself);

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  for(int j = 0; j < jnum; j++) {
    x = rij[j][0];
    y = rij[j][1];
    z = rij[j][2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);

    compute_uarray(x, y, z, z0, r);
    add_uarraytot(r, wj[j], rcutij[j]);
  }

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[0] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

void SNA::compute_ui_omp(int jnum, int sub_threads)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma 
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  zero_uarraytot();
  addself_uarraytot(wself);

  for(int j = 0; j < jnum; j++) {
    x = rij[j][0];
    y = rij[j][1];
    z = rij[j][2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);
    theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);
    omp_set_num_threads(sub_threads);

#if defined(_OPENMP)
#pragma omp parallel shared(x,y,z,z0,r,sub_threads) default(none)
#endif
    {
      compute_uarray_omp(x, y, z, z0, r, sub_threads);
    }
    add_uarraytot(r, wj[j], rcutij[j]);
  }


}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

void SNA::compute_zi()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        for ma = 0,...,j
  //          for mb = 0,...,jmid
  //            z(j1,j2,j,ma,mb) = 0
  //            for ma1 = Max(0,ma+(j1-j2-j)/2),Min(j1,ma+(j1+j2-j)/2)
  //              sumb1 = 0
  //              ma2 = ma-ma1+(j1+j2-j)/2;
  //              for mb1 = Max(0,mb+(j1-j2-j)/2),Min(j1,mb+(j1+j2-j)/2)
  //                mb2 = mb-mb1+(j1+j2-j)/2;
  //                sumb1 += cg(j1,mb1,j2,mb2,j) *
  //                  u(j1,ma1,mb1) * u(j2,ma2,mb2)
  //              z(j1,j2,j,ma,mb) += sumb1*cg(j1,ma1,j2,ma2,j)

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  // compute_dbidrj() requires full j1/j2/j chunk of z elements
  // use zarray j1/j2 symmetry

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++) {
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
	double sumb1_r, sumb1_i;
	int ma2, mb2;
	for(int mb = 0; 2*mb <= j; mb++)
	  for(int ma = 0; ma <= j; ma++) {
	    zarray_r[j1][j2][j][ma][mb] = 0.0;
	    zarray_i[j1][j2][j][ma][mb] = 0.0;
	    
	    for(int ma1 = MAX(0, (2 * ma - j - j2 + j1) / 2);
		ma1 <= MIN(j1, (2 * ma - j + j2 + j1) / 2); ma1++) {
	      sumb1_r = 0.0;
	      sumb1_i = 0.0;

	      ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;

	      for(int mb1 = MAX(0, (2 * mb - j - j2 + j1) / 2);
              mb1 <= MIN(j1, (2 * mb - j + j2 + j1) / 2); mb1++) {

		mb2 = (2 * mb - j - (2 * mb1 - j1) + j2) / 2;
		sumb1_r += cgarray[j1][j2][j][mb1][mb2] * 
		  (uarraytot_r[j1][ma1][mb1] * uarraytot_r[j2][ma2][mb2] -
		   uarraytot_i[j1][ma1][mb1] * uarraytot_i[j2][ma2][mb2]);
		sumb1_i += cgarray[j1][j2][j][mb1][mb2] * 
		  (uarraytot_r[j1][ma1][mb1] * uarraytot_i[j2][ma2][mb2] +
		   uarraytot_i[j1][ma1][mb1] * uarraytot_r[j2][ma2][mb2]);
	      } // end loop over mb1

	      zarray_r[j1][j2][j][ma][mb] +=
		sumb1_r * cgarray[j1][j2][j][ma1][ma2];
	      zarray_i[j1][j2][j][ma][mb] +=
		sumb1_i * cgarray[j1][j2][j][ma1][ma2];
	    } // end loop over ma1
	  } // end loop over ma, mb
      } // end loop over j
    } // end loop over j1, j2

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[1] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif
}

void SNA::compute_zi_omp(int sub_threads)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        for ma = 0,...,j
  //          for mb = 0,...,j
  //            z(j1,j2,j,ma,mb) = 0
  //            for ma1 = Max(0,ma+(j1-j2-j)/2),Min(j1,ma+(j1+j2-j)/2)
  //              sumb1 = 0
  //              ma2 = ma-ma1+(j1+j2-j)/2;
  //              for mb1 = Max(0,mb+(j1-j2-j)/2),Min(j1,mb+(j1+j2-j)/2)
  //                mb2 = mb-mb1+(j1+j2-j)/2;
  //                sumb1 += cg(j1,mb1,j2,mb2,j) *
  //                  u(j1,ma1,mb1) * u(j2,ma2,mb2)
  //              z(j1,j2,j,ma,mb) += sumb1*cg(j1,ma1,j2,ma2,j)

  if(omp_in_parallel())
    omp_set_num_threads(sub_threads);

  // compute_dbidrj() requires full j1/j2/j chunk of z elements
  // use zarray j1/j2 symmetry

#if defined(_OPENMP)
#pragma omp parallel for schedule(auto) default(none)
#endif
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {

    double sumb1_r, sumb1_i;
    int ma2, mb2;

    for(int ma = 0; ma <= j; ma++)
      for(int mb = 0; mb <= j; mb++) {
        zarray_r[j1][j2][j][ma][mb] = 0.0;
        zarray_i[j1][j2][j][ma][mb] = 0.0;

        for(int ma1 = MAX(0, (2 * ma - j - j2 + j1) / 2);
            ma1 <= MIN(j1, (2 * ma - j + j2 + j1) / 2); ma1++) {
          sumb1_r = 0.0;
          sumb1_i = 0.0;

          ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;

          for(int mb1 = MAX(0, (2 * mb - j - j2 + j1) / 2);
              mb1 <= MIN(j1, (2 * mb - j + j2 + j1) / 2); mb1++) {

            mb2 = (2 * mb - j - (2 * mb1 - j1) + j2) / 2;
            sumb1_r += cgarray[j1][j2][j][mb1][mb2] *
	      (uarraytot_r[j1][ma1][mb1] * uarraytot_r[j2][ma2][mb2] -
	       uarraytot_i[j1][ma1][mb1] * uarraytot_i[j2][ma2][mb2]);
            sumb1_i += cgarray[j1][j2][j][mb1][mb2] *
	      (uarraytot_r[j1][ma1][mb1] * uarraytot_i[j2][ma2][mb2] +
	       uarraytot_i[j1][ma1][mb1] * uarraytot_r[j2][ma2][mb2]);
          }

          zarray_r[j1][j2][j][ma][mb] +=
            sumb1_r * cgarray[j1][j2][j][ma1][ma2];
          zarray_i[j1][j2][j][ma][mb] +=
            sumb1_i * cgarray[j1][j2][j][ma1][ma2];
        }
      }
  }
}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
------------------------------------------------------------------------- */

void SNA::compute_bi()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++) {
      for(int j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        barray[j1][j2][j] = 0.0;

	for(int mb = 0; 2*mb < j; mb++) {
	  for(int ma = 0; ma <= j; ma++) {
            barray[j1][j2][j] +=
              uarraytot_r[j][ma][mb] * zarray_r[j1][j2][j][ma][mb] +
	      uarraytot_i[j][ma][mb] * zarray_i[j1][j2][j][ma][mb];
	  }
	}

	// For j even, special treatment for middle column

	if (j%2 == 0) {
	  int mb = j/2;
	  for(int ma = 0; ma < mb; ma++)
	    barray[j1][j2][j] +=
	      uarraytot_r[j][ma][mb] * zarray_r[j1][j2][j][ma][mb] +
	      uarraytot_i[j][ma][mb] * zarray_i[j1][j2][j][ma][mb];
	  int ma = mb;
	  barray[j1][j2][j] +=
	    (uarraytot_r[j][ma][mb] * zarray_r[j1][j2][j][ma][mb] +
	     uarraytot_i[j][ma][mb] * zarray_i[j1][j2][j][ma][mb])*0.5;
	}

        barray[j1][j2][j] *= 2.0;
      }
    }

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[2] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   copy Bi array to a vector
------------------------------------------------------------------------- */

void SNA::copy_bi2bvec()
{
  int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++)
    if(diagonalstyle == 0) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2) {
          bvec[ncount] = barray[j1][j2][j];
          ncount++;
        }
    } else if(diagonalstyle == 1) {
      j2 = j1;
      for(j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        bvec[ncount] = barray[j1][j2][j];
        ncount++;
      }
    } else if(diagonalstyle == 2) {
      j = j2 = j1;
      bvec[ncount] = barray[j1][j2][j];
      ncount++;
    } else if(diagonalstyle == 3) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
	  if (j >= j1) {
	    bvec[ncount] = barray[j1][j2][j];
	    ncount++;
	  }
    }
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

void SNA::compute_duidrj(double* rij, double wj, double rcut)
{
  double rsq, r, x, y, z, z0, theta0, cs, sn;
  double dz0dr;

  x = rij[0];
  y = rij[1];
  z = rij[2];
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  double rscale0 = rfac0 * MY_PI / (rcut - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  compute_duarray(x, y, z, z0, r, dz0dr, wj, rcut);

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[3] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   calculate derivative of Bi w.r.t. atom j
   variant using indexlist for j1,j2,j
   variant not using symmetry relation
------------------------------------------------------------------------- */

void SNA::compute_dbidrj_nonsymm()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        dbdr(j1,j2,j) = 0
  //        for ma = 0,...,j
  //          for mb = 0,...,j
  //            dzdr = 0
  //            for ma1 = Max(0,ma+(j1-j2-j)/2),Min(j1,ma+(j1+j2-j)/2)
  //              sumb1 = 0
  //              ma2 = ma-ma1+(j1+j2-j)/2;
  //              for mb1 = Max(0,mb+(j1-j2-j)/2),Min(j1,mb+(j1+j2-j)/2)
  //                mb2 = mb-mb1+(j1+j2-j)/2;
  //                sumb1 += cg(j1,mb1,j2,mb2,j) *
  //                  (dudr(j1,ma1,mb1) * u(j2,ma2,mb2) +
  //                  u(j1,ma1,mb1) * dudr(j2,ma2,mb2))
  //              dzdr += sumb1*cg(j1,ma1,j2,ma2,j)
  //            dbdr(j1,j2,j) +=
  //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb) +
  //              Conj(u(j,ma,mb))*dzdr

  double* dbdr;
  double* dudr_r, *dudr_i;
  double sumb1_r[3], sumb1_i[3], dzdr_r[3], dzdr_i[3];
  int ma2;

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  for(int JJ = 0; JJ < idxj_max; JJ++) {
    const int j1 = idxj[JJ].j1;
    const int j2 = idxj[JJ].j2;
    const int j = idxj[JJ].j;

    dbdr = dbarray[j1][j2][j];
    dbdr[0] = 0.0;
    dbdr[1] = 0.0;
    dbdr[2] = 0.0;

    double** *j1duarray_r = duarray_r[j1];
    double** *j2duarray_r = duarray_r[j2];
    double** *j1duarray_i = duarray_i[j1];
    double** *j2duarray_i = duarray_i[j2];
    double** j1uarraytot_r = uarraytot_r[j1];
    double** j2uarraytot_r = uarraytot_r[j2];
    double** j1uarraytot_i = uarraytot_i[j1];
    double** j2uarraytot_i = uarraytot_i[j2];
    double** j1j2jcgarray = cgarray[j1][j2][j];

    for(int ma = 0; ma <= j; ma++)
      for(int mb = 0; mb <= j; mb++) {
        dzdr_r[0] = 0.0;
        dzdr_r[1] = 0.0;
        dzdr_r[2] = 0.0;
        dzdr_i[0] = 0.0;
        dzdr_i[1] = 0.0;
        dzdr_i[2] = 0.0;

        const int max_mb1 = MIN(j1, (2 * mb - j + j2 + j1) / 2) + 1;
        const int max_ma1 = MIN(j1, (2 * ma - j + j2 + j1) / 2) + 1;

        for(int ma1 = MAX(0, (2 * ma - j - j2 + j1) / 2);
            ma1 < max_ma1; ma1++) {

          ma2 = (2 * ma - j - (2 * ma1 - j1) + j2) / 2;
          sumb1_r[0] = 0.0;
          sumb1_r[1] = 0.0;
          sumb1_r[2] = 0.0;
          sumb1_i[0] = 0.0;
          sumb1_i[1] = 0.0;
          sumb1_i[2] = 0.0;

          //inside loop 54 operations (mul and add)
          for(int mb1 = MAX(0, (2 * mb - j - j2 + j1) / 2),
              mb2 = mb + (j1 + j2 - j) / 2 - mb1;
              mb1 < max_mb1; mb1++, mb2--) {

            double* dudr1_r, *dudr1_i, *dudr2_r, *dudr2_i;

            dudr1_r = j1duarray_r[ma1][mb1];
            dudr2_r = j2duarray_r[ma2][mb2];
            dudr1_i = j1duarray_i[ma1][mb1];
            dudr2_i = j2duarray_i[ma2][mb2];

            const double cga_mb1mb2 = j1j2jcgarray[mb1][mb2];
            const double uat_r_ma2mb2 = cga_mb1mb2 * j2uarraytot_r[ma2][mb2];
            const double uat_r_ma1mb1 = cga_mb1mb2 * j1uarraytot_r[ma1][mb1];
            const double uat_i_ma2mb2 = cga_mb1mb2 * j2uarraytot_i[ma2][mb2];
            const double uat_i_ma1mb1 = cga_mb1mb2 * j1uarraytot_i[ma1][mb1];

            for(int k = 0; k < 3; k++) {
              sumb1_r[k] += dudr1_r[k] * uat_r_ma2mb2;
              sumb1_r[k] -= dudr1_i[k] * uat_i_ma2mb2;
              sumb1_i[k] += dudr1_r[k] * uat_i_ma2mb2;
              sumb1_i[k] += dudr1_i[k] * uat_r_ma2mb2;

              sumb1_r[k] += dudr2_r[k] * uat_r_ma1mb1;
              sumb1_r[k] -= dudr2_i[k] * uat_i_ma1mb1;
              sumb1_i[k] += dudr2_r[k] * uat_i_ma1mb1;
              sumb1_i[k] += dudr2_i[k] * uat_r_ma1mb1;
            }
          } // end loop over mb1,mb2

          // dzdr += sumb1*cg(j1,ma1,j2,ma2,j)

          dzdr_r[0] += sumb1_r[0] * j1j2jcgarray[ma1][ma2];
          dzdr_r[1] += sumb1_r[1] * j1j2jcgarray[ma1][ma2];
          dzdr_r[2] += sumb1_r[2] * j1j2jcgarray[ma1][ma2];
          dzdr_i[0] += sumb1_i[0] * j1j2jcgarray[ma1][ma2];
          dzdr_i[1] += sumb1_i[1] * j1j2jcgarray[ma1][ma2];
          dzdr_i[2] += sumb1_i[2] * j1j2jcgarray[ma1][ma2];
        } // end loop over ma1,ma2

        // dbdr(j1,j2,j) +=
        //   Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb) +
        //   Conj(u(j,ma,mb))*dzdr

        dudr_r = duarray_r[j][ma][mb];
        dudr_i = duarray_i[j][ma][mb];

        for(int k = 0; k < 3; k++)
          dbdr[k] +=
            (dudr_r[k] * zarray_r[j1][j2][j][ma][mb] +
             dudr_i[k] * zarray_i[j1][j2][j][ma][mb]) +
            (uarraytot_r[j][ma][mb] * dzdr_r[k] +
             uarraytot_i[j][ma][mb] * dzdr_i[k]);
      } //end loop over ma mb

  } //end loop over j1 j2 j

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[4] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   calculate derivative of Bi w.r.t. atom j
   variant using indexlist for j1,j2,j
   variant using symmetry relation
------------------------------------------------------------------------- */

void SNA::compute_dbidrj()
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        zdb = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            zdb +=
  //              Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
  //        dbdr(j1,j2,j) += 2*zdb
  //        zdb = 0
  //        for mb1 = 0,...,j1mid
  //          for ma1 = 0,...,j1
  //            zdb +=
  //              Conj(dudr(j1,ma1,mb1))*z(j1,j2,j,ma1,mb1)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j1+1)
  //        zdb = 0
  //        for mb2 = 0,...,j2mid
  //          for ma2 = 0,...,j2
  //            zdb +=
  //              Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)
  //        dbdr(j1,j2,j) += 2*zdb*(j+1)/(j2+1)

  double* dbdr;
  double* dudr_r, *dudr_i;
  double sumzdu_r[3];
  double** jjjzarray_r;
  double** jjjzarray_i;
  double jjjmambzarray_r; 
  double jjjmambzarray_i;

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &starttime);
#endif

  for(int JJ = 0; JJ < idxj_max; JJ++) {
    const int j1 = idxj[JJ].j1;
    const int j2 = idxj[JJ].j2;
    const int j = idxj[JJ].j;

    dbdr = dbarray[j1][j2][j];
    dbdr[0] = 0.0;
    dbdr[1] = 0.0;
    dbdr[2] = 0.0;

    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    for(int k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j1 >= j2) {
      jjjzarray_r = zarray_r[j1][j2][j];
      jjjzarray_i = zarray_i[j1][j2][j];
    } else {
      jjjzarray_r = zarray_r[j2][j1][j];
      jjjzarray_i = zarray_i[j2][j1][j];
    }

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {

        dudr_r = duarray_r[j][ma][mb];
        dudr_i = duarray_i[j][ma][mb];
	jjjmambzarray_r = jjjzarray_r[ma][mb];
	jjjmambzarray_i = jjjzarray_i[ma][mb];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {
      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        dudr_r = duarray_r[j][ma][mb];
	dudr_i = duarray_i[j][ma][mb];
	jjjmambzarray_r = jjjzarray_r[ma][mb];
	jjjmambzarray_i = jjjzarray_i[ma][mb];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;
      }
      int ma = mb;
      dudr_r = duarray_r[j][ma][mb];
      dudr_i = duarray_i[j][ma][mb];
      jjjmambzarray_r = jjjzarray_r[ma][mb];
      jjjmambzarray_i = jjjzarray_i[ma][mb];
      for(int k = 0; k < 3; k++)
	sumzdu_r[k] +=
	  (dudr_r[k] * jjjmambzarray_r +
	   dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if jeven

    for(int k = 0; k < 3; k++)
      dbdr[k] += 2.0*sumzdu_r[k];

    // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

    double j1fac = (j+1)/(j1+1.0);

    for(int k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j >= j2) {
      jjjzarray_r = zarray_r[j][j2][j1];
      jjjzarray_i = zarray_i[j][j2][j1];
    } else {
      jjjzarray_r = zarray_r[j2][j][j1];
      jjjzarray_i = zarray_i[j2][j][j1];
    }

    for(int mb1 = 0; 2*mb1 < j1; mb1++)
      for(int ma1 = 0; ma1 <= j1; ma1++) {

        dudr_r = duarray_r[j1][ma1][mb1];
        dudr_i = duarray_i[j1][ma1][mb1];
	jjjmambzarray_r = jjjzarray_r[ma1][mb1];
	jjjmambzarray_i = jjjzarray_i[ma1][mb1];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma1 mb1

    // For j1 even, handle middle column

    if (j1%2 == 0) {
      int mb1 = j1/2;
      for(int ma1 = 0; ma1 < mb1; ma1++) {
        dudr_r = duarray_r[j1][ma1][mb1];
	dudr_i = duarray_i[j1][ma1][mb1];
	jjjmambzarray_r = jjjzarray_r[ma1][mb1];
	jjjmambzarray_i = jjjzarray_i[ma1][mb1];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;
      }
      int ma1 = mb1;
      dudr_r = duarray_r[j1][ma1][mb1];
      dudr_i = duarray_i[j1][ma1][mb1];
      jjjmambzarray_r = jjjzarray_r[ma1][mb1];
      jjjmambzarray_i = jjjzarray_i[ma1][mb1];
      for(int k = 0; k < 3; k++)
	sumzdu_r[k] +=
	  (dudr_r[k] * jjjmambzarray_r +
	   dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if j1even

    for(int k = 0; k < 3; k++)
      dbdr[k] += 2.0*sumzdu_r[k]*j1fac;

    // Sum over Conj(dudr(j2,ma2,mb2))*z(j1,j,j2,ma2,mb2)

    double j2fac = (j+1)/(j2+1.0);

    for(int k = 0; k < 3; k++)
      sumzdu_r[k] = 0.0;

    // use zarray j1/j2 symmetry (optional)

    if (j1 >= j) {
      jjjzarray_r = zarray_r[j1][j][j2];
      jjjzarray_i = zarray_i[j1][j][j2];
    } else {
      jjjzarray_r = zarray_r[j][j1][j2];
      jjjzarray_i = zarray_i[j][j1][j2];
    }

    for(int mb2 = 0; 2*mb2 < j2; mb2++)
      for(int ma2 = 0; ma2 <= j2; ma2++) {

        dudr_r = duarray_r[j2][ma2][mb2];
        dudr_i = duarray_i[j2][ma2][mb2];
	jjjmambzarray_r = jjjzarray_r[ma2][mb2];
	jjjmambzarray_i = jjjzarray_i[ma2][mb2];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;

      } //end loop over ma2 mb2

    // For j2 even, handle middle column

    if (j2%2 == 0) {
      int mb2 = j2/2;
      for(int ma2 = 0; ma2 < mb2; ma2++) {
        dudr_r = duarray_r[j2][ma2][mb2];
	dudr_i = duarray_i[j2][ma2][mb2];
	jjjmambzarray_r = jjjzarray_r[ma2][mb2];
	jjjmambzarray_i = jjjzarray_i[ma2][mb2];
        for(int k = 0; k < 3; k++)
          sumzdu_r[k] +=
            dudr_r[k] * jjjmambzarray_r +
	    dudr_i[k] * jjjmambzarray_i;
      }
      int ma2 = mb2;
      dudr_r = duarray_r[j2][ma2][mb2];
      dudr_i = duarray_i[j2][ma2][mb2];
      jjjmambzarray_r = jjjzarray_r[ma2][mb2];
      jjjmambzarray_i = jjjzarray_i[ma2][mb2];
      for(int k = 0; k < 3; k++)
	sumzdu_r[k] +=
	  (dudr_r[k] * jjjmambzarray_r +
	   dudr_i[k] * jjjmambzarray_i)*0.5;
    } // end if j2even

    for(int k = 0; k < 3; k++)
      dbdr[k] += 2.0*sumzdu_r[k]*j2fac;

  } //end loop over j1 j2 j

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME, &endtime);
  timers[4] += (endtime.tv_sec - starttime.tv_sec + 1.0 *
                (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
#endif

}

/* ----------------------------------------------------------------------
   copy Bi derivatives into a vector
------------------------------------------------------------------------- */

void SNA::copy_dbi2dbvec()
{
  int ncount, j1, j2, j;

  ncount = 0;

  for(j1 = 0; j1 <= twojmax; j1++) {
    if(diagonalstyle == 0) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2) {
          dbvec[ncount][0] = dbarray[j1][j2][j][0];
          dbvec[ncount][1] = dbarray[j1][j2][j][1];
          dbvec[ncount][2] = dbarray[j1][j2][j][2];
          ncount++;
        }
    } else if(diagonalstyle == 1) {
      j2 = j1;
      for(j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2) {
        dbvec[ncount][0] = dbarray[j1][j2][j][0];
        dbvec[ncount][1] = dbarray[j1][j2][j][1];
        dbvec[ncount][2] = dbarray[j1][j2][j][2];
        ncount++;
      }
    } else if(diagonalstyle == 2) {
      j = j2 = j1;
      dbvec[ncount][0] = dbarray[j1][j2][j][0];
      dbvec[ncount][1] = dbarray[j1][j2][j][1];
      dbvec[ncount][2] = dbarray[j1][j2][j][2];
      ncount++;
    } else if(diagonalstyle == 3) {
      for(j2 = 0; j2 <= j1; j2++)
        for(j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
	  if (j >= j1) {
	    dbvec[ncount][0] = dbarray[j1][j2][j][0];
	    dbvec[ncount][1] = dbarray[j1][j2][j][1];
	    dbvec[ncount][2] = dbarray[j1][j2][j][2];
	    ncount++;
	  }
    }
  }
}

/* ---------------------------------------------------------------------- */

void SNA::zero_uarraytot()
{
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        uarraytot_r[j][ma][mb] = 0.0;
        uarraytot_i[j][ma][mb] = 0.0;
      }
}

/* ---------------------------------------------------------------------- */

void SNA::addself_uarraytot(double wself_in)
{
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++) {
      uarraytot_r[j][ma][ma] = wself_in;
      uarraytot_i[j][ma][ma] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

void SNA::add_uarraytot(double r, double wj, double rcut)
{
  double sfac;

  sfac = compute_sfac(r, rcut);

  sfac *= wj;

  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        uarraytot_r[j][ma][mb] +=
          sfac * uarray_r[j][ma][mb];
        uarraytot_i[j][ma][mb] +=
          sfac * uarray_i[j][ma][mb];
      }
}

void SNA::add_uarraytot_omp(double r, double wj, double rcut)
{
  double sfac;

  sfac = compute_sfac(r, rcut);

  sfac *= wj;

#if defined(_OPENMP)
#pragma omp for
#endif
  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        uarraytot_r[j][ma][mb] +=
          sfac * uarray_r[j][ma][mb];
        uarraytot_i[j][ma][mb] +=
          sfac * uarray_i[j][ma][mb];
      }
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

void SNA::compute_uarray(double x, double y, double z,
                         double z0, double r)
{
  double r0inv;
  double a_r, b_r, a_i, b_i;
  double rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2

  uarray_r[0][0][0] = 1.0;
  uarray_i[0][0][0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {

    // fill in left side of matrix layer from previous layer

    for (int mb = 0; 2*mb <= j; mb++) {
      uarray_r[j][0][mb] = 0.0;
      uarray_i[j][0][mb] = 0.0;

      for (int ma = 0; ma < j; ma++) {
	rootpq = rootpqarray[j - ma][j - mb];
        uarray_r[j][ma][mb] +=
          rootpq *
          (a_r * uarray_r[j - 1][ma][mb] + 
	   a_i * uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma][mb] +=
          rootpq *
          (a_r * uarray_i[j - 1][ma][mb] - 
	   a_i * uarray_r[j - 1][ma][mb]);

	rootpq = rootpqarray[ma + 1][j - mb];
        uarray_r[j][ma + 1][mb] =
          -rootpq *
          (b_r * uarray_r[j - 1][ma][mb] + 
	   b_i * uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma + 1][mb] =
          -rootpq *
          (b_r * uarray_i[j - 1][ma][mb] - 
	   b_i * uarray_r[j - 1][ma][mb]);
      }
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    int mbpar = -1;
    for (int mb = 0; 2*mb <= j; mb++) {
      mbpar = -mbpar;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
    	mapar = -mapar;
    	if (mapar == 1) {
    	  uarray_r[j][j-ma][j-mb] = uarray_r[j][ma][mb];
    	  uarray_i[j][j-ma][j-mb] = -uarray_i[j][ma][mb];
    	} else {
    	  uarray_r[j][j-ma][j-mb] = -uarray_r[j][ma][mb];
    	  uarray_i[j][j-ma][j-mb] = uarray_i[j][ma][mb];
    	}	  
      }
    }
  }
}

void SNA::compute_uarray_omp(double x, double y, double z,
                             double z0, double r, int sub_threads)
{
  double r0inv;
  double a_r, b_r, a_i, b_i;
  double rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2

  uarray_r[0][0][0] = 1.0;
  uarray_i[0][0][0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
#if defined(_OPENMP)
#pragma omp for
#endif
    for (int mb = 0; mb < j; mb++) {
      uarray_r[j][0][mb] = 0.0;
      uarray_i[j][0][mb] = 0.0;

      for (int ma = 0; ma < j; ma++) {
	rootpq = rootpqarray[j - ma][j - mb];
        uarray_r[j][ma][mb] +=
	  rootpq *
          (a_r * uarray_r[j - 1][ma][mb] + 
	   a_i * uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma][mb] +=
	  rootpq *
          (a_r * uarray_i[j - 1][ma][mb] - 
	   a_i * uarray_r[j - 1][ma][mb]);

	rootpq = rootpqarray[ma + 1][j - mb];
        uarray_r[j][ma + 1][mb] =
	  -rootpq *
          (b_r * uarray_r[j - 1][ma][mb] + 
	   b_i * uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma + 1][mb] =
	  -rootpq *
          (b_r * uarray_i[j - 1][ma][mb] - 
	   b_i * uarray_r[j - 1][ma][mb]);
      }
    }

    int mb = j;
    uarray_r[j][0][mb] = 0.0;
    uarray_i[j][0][mb] = 0.0;

#if defined(_OPENMP)
#pragma omp for
#endif
    for (int ma = 0; ma < j; ma++) {
      rootpq = rootpqarray[j - ma][mb];
      uarray_r[j][ma][mb] +=
	rootpq *
        (b_r * uarray_r[j - 1][ma][mb - 1] - 
	 b_i * uarray_i[j - 1][ma][mb - 1]);
      uarray_i[j][ma][mb] +=
	rootpq *
        (b_r * uarray_i[j - 1][ma][mb - 1] + 
	 b_i * uarray_r[j - 1][ma][mb - 1]);

      rootpq = rootpqarray[ma + 1][mb];
      uarray_r[j][ma + 1][mb] =
	rootpq *
        (a_r * uarray_r[j - 1][ma][mb - 1] - 
	 a_i * uarray_i[j - 1][ma][mb - 1]);
      uarray_i[j][ma + 1][mb] =
	rootpq *
        (a_r * uarray_i[j - 1][ma][mb - 1] + 
	 a_i * uarray_r[j - 1][ma][mb - 1]);
    }
  }
}

/* ----------------------------------------------------------------------
   compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

void SNA::compute_duarray(double x, double y, double z,
                          double z0, double r, double dz0dr,
			  double wj, double rcut)
{
  double r0inv;
  double a_r, a_i, b_r, b_i;
  double da_r[3], da_i[3], db_r[3], db_i[3];
  double dz0[3], dr0inv[3], dr0invdr;
  double rootpq;

  double rinv = 1.0 / r;
  double ux = x * rinv;
  double uy = y * rinv;
  double uz = z * rinv;

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

  dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  for (int k = 0; k < 3; k++) {
    da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }

  da_i[2] += -r0inv;

  for (int k = 0; k < 3; k++) {
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }

  db_i[0] += -r0inv;
  db_r[1] += r0inv;

  uarray_r[0][0][0] = 1.0;
  duarray_r[0][0][0][0] = 0.0;
  duarray_r[0][0][0][1] = 0.0;
  duarray_r[0][0][0][2] = 0.0;
  uarray_i[0][0][0] = 0.0;
  duarray_i[0][0][0][0] = 0.0;
  duarray_i[0][0][0][1] = 0.0;
  duarray_i[0][0][0][2] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    for (int mb = 0; 2*mb <= j; mb++) {
      uarray_r[j][0][mb] = 0.0;
      duarray_r[j][0][mb][0] = 0.0;
      duarray_r[j][0][mb][1] = 0.0;
      duarray_r[j][0][mb][2] = 0.0;
      uarray_i[j][0][mb] = 0.0;
      duarray_i[j][0][mb][0] = 0.0;
      duarray_i[j][0][mb][1] = 0.0;
      duarray_i[j][0][mb][2] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[j - ma][j - mb];
        uarray_r[j][ma][mb] += rootpq *
                               (a_r *  uarray_r[j - 1][ma][mb] +
                                a_i *  uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma][mb] += rootpq *
                               (a_r *  uarray_i[j - 1][ma][mb] -
                                a_i *  uarray_r[j - 1][ma][mb]);

        for (int k = 0; k < 3; k++) {
          duarray_r[j][ma][mb][k] +=
            rootpq * (da_r[k] * uarray_r[j - 1][ma][mb] +
                      da_i[k] * uarray_i[j - 1][ma][mb] +
                      a_r * duarray_r[j - 1][ma][mb][k] +
                      a_i * duarray_i[j - 1][ma][mb][k]);
          duarray_i[j][ma][mb][k] +=
            rootpq * (da_r[k] * uarray_i[j - 1][ma][mb] -
                      da_i[k] * uarray_r[j - 1][ma][mb] +
                      a_r * duarray_i[j - 1][ma][mb][k] -
                      a_i * duarray_r[j - 1][ma][mb][k]);
        }

	rootpq = rootpqarray[ma + 1][j - mb];
        uarray_r[j][ma + 1][mb] =
          -rootpq * (b_r *  uarray_r[j - 1][ma][mb] +
                     b_i *  uarray_i[j - 1][ma][mb]);
        uarray_i[j][ma + 1][mb] =
          -rootpq * (b_r *  uarray_i[j - 1][ma][mb] -
                     b_i *  uarray_r[j - 1][ma][mb]);

        for (int k = 0; k < 3; k++) {
          duarray_r[j][ma + 1][mb][k] =
            -rootpq * (db_r[k] * uarray_r[j - 1][ma][mb] +
                       db_i[k] * uarray_i[j - 1][ma][mb] +
                       b_r * duarray_r[j - 1][ma][mb][k] +
                       b_i * duarray_i[j - 1][ma][mb][k]);
          duarray_i[j][ma + 1][mb][k] =
            -rootpq * (db_r[k] * uarray_i[j - 1][ma][mb] -
                       db_i[k] * uarray_r[j - 1][ma][mb] +
                       b_r * duarray_i[j - 1][ma][mb][k] -
                       b_i * duarray_r[j - 1][ma][mb][k]);
        }
      }
    }

    int mbpar = -1;
    for (int mb = 0; 2*mb <= j; mb++) {
      mbpar = -mbpar;
      int mapar = -mbpar;
      for (int ma = 0; ma <= j; ma++) {
    	mapar = -mapar;
    	if (mapar == 1) {
    	  uarray_r[j][j-ma][j-mb] = uarray_r[j][ma][mb];
    	  uarray_i[j][j-ma][j-mb] = -uarray_i[j][ma][mb];
    	  for (int k = 0; k < 3; k++) {
    	    duarray_r[j][j-ma][j-mb][k] = duarray_r[j][ma][mb][k];
    	    duarray_i[j][j-ma][j-mb][k] = -duarray_i[j][ma][mb][k];
    	  }
    	} else {
    	  uarray_r[j][j-ma][j-mb] = -uarray_r[j][ma][mb];
    	  uarray_i[j][j-ma][j-mb] = uarray_i[j][ma][mb];
    	  for (int k = 0; k < 3; k++) {
    	    duarray_r[j][j-ma][j-mb][k] = -duarray_r[j][ma][mb][k];
    	    duarray_i[j][j-ma][j-mb][k] = duarray_i[j][ma][mb][k];
    	  }
    	}
      }
    }
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;

  for (int j = 0; j <= twojmax; j++)
    for (int ma = 0; ma <= j; ma++)
      for (int mb = 0; mb <= j; mb++) {
        duarray_r[j][ma][mb][0] = dsfac * uarray_r[j][ma][mb] * ux +
                                  sfac * duarray_r[j][ma][mb][0];
        duarray_i[j][ma][mb][0] = dsfac * uarray_i[j][ma][mb] * ux +
                                  sfac * duarray_i[j][ma][mb][0];
        duarray_r[j][ma][mb][1] = dsfac * uarray_r[j][ma][mb] * uy +
                                  sfac * duarray_r[j][ma][mb][1];
        duarray_i[j][ma][mb][1] = dsfac * uarray_i[j][ma][mb] * uy +
                                  sfac * duarray_i[j][ma][mb][1];
        duarray_r[j][ma][mb][2] = dsfac * uarray_r[j][ma][mb] * uz +
                                  sfac * duarray_r[j][ma][mb][2];
        duarray_i[j][ma][mb][2] = dsfac * uarray_i[j][ma][mb] * uz +
                                  sfac * duarray_i[j][ma][mb][2];
      }
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

double SNA::memory_usage()
{
  int jdim = twojmax + 1;
  double bytes;
  bytes = jdim * jdim * jdim * jdim * jdim * sizeof(double);
  bytes += 2 * jdim * jdim * jdim * sizeof(complex<double>);
  bytes += 2 * jdim * jdim * jdim * sizeof(double);
  bytes += jdim * jdim * jdim * 3 * sizeof(complex<double>);
  bytes += jdim * jdim * jdim * 3 * sizeof(double);
  bytes += ncoeff * sizeof(double);
  bytes += jdim * jdim * jdim * jdim * jdim * sizeof(complex<double>);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void SNA::create_twojmax_arrays()
{
  int jdim = twojmax + 1;

  memory->create(cgarray, jdim, jdim, jdim, jdim, jdim,
                 "sna:cgarray");
  memory->create(rootpqarray, jdim+1, jdim+1,
                 "sna:rootpqarray");
  memory->create(barray, jdim, jdim, jdim,
                 "sna:barray");
  memory->create(dbarray, jdim, jdim, jdim, 3,
                 "sna:dbarray");

  memory->create(duarray_r, jdim, jdim, jdim, 3,
                 "sna:duarray");
  memory->create(duarray_i, jdim, jdim, jdim, 3,
                 "sna:duarray");

  memory->create(uarray_r, jdim, jdim, jdim,
                 "sna:uarray");
  memory->create(uarray_i, jdim, jdim, jdim,
                 "sna:uarray");

  if(!use_shared_arrays) {
    memory->create(uarraytot_r, jdim, jdim, jdim,
                   "sna:uarraytot");
    memory->create(zarray_r, jdim, jdim, jdim, jdim, jdim,
                   "sna:zarray");
    memory->create(uarraytot_i, jdim, jdim, jdim,
                   "sna:uarraytot");
    memory->create(zarray_i, jdim, jdim, jdim, jdim, jdim,
                   "sna:zarray");
  }

}

/* ---------------------------------------------------------------------- */

void SNA::destroy_twojmax_arrays()
{
  memory->destroy(cgarray);
  memory->destroy(rootpqarray);
  memory->destroy(barray);

  memory->destroy(dbarray);

  memory->destroy(duarray_r);
  memory->destroy(duarray_i);

  memory->destroy(uarray_r);
  memory->destroy(uarray_i);

  if(!use_shared_arrays) {
    memory->destroy(uarraytot_r);
    memory->destroy(zarray_r);
    memory->destroy(uarraytot_i);
    memory->destroy(zarray_i);
  }
}


/* ----------------------------------------------------------------------
   store n! and return it as needed
------------------------------------------------------------------------- */

double SNA::factorial(int n)
{
  const int nmax = 167; //  Largest n supported
  double nfac_table[nmax+1] = {
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
    355687428096000,
    6.402373705728e+15,
    1.21645100408832e+17,
    2.43290200817664e+18,
    5.10909421717094e+19,
    1.12400072777761e+21,
    2.5852016738885e+22,
    6.20448401733239e+23,
    1.5511210043331e+25,
    4.03291461126606e+26,
    1.08888694504184e+28,
    3.04888344611714e+29,
    8.8417619937397e+30,
    2.65252859812191e+32,
    8.22283865417792e+33,
    2.63130836933694e+35,
    8.68331761881189e+36,
    2.95232799039604e+38,
    1.03331479663861e+40,
    3.71993326789901e+41,
    1.37637530912263e+43,
    5.23022617466601e+44,
    2.03978820811974e+46,
    8.15915283247898e+47,
    3.34525266131638e+49,
    1.40500611775288e+51,
    6.04152630633738e+52,
    2.65827157478845e+54,
    1.1962222086548e+56,
    5.50262215981209e+57,
    2.58623241511168e+59,
    1.24139155925361e+61,
    6.08281864034268e+62,
    3.04140932017134e+64,
    1.55111875328738e+66,
    8.06581751709439e+67,
    4.27488328406003e+69,
    2.30843697339241e+71,
    1.26964033536583e+73,
    7.10998587804863e+74,
    4.05269195048772e+76,
    2.35056133128288e+78,
    1.3868311854569e+80,
    8.32098711274139e+81,
    5.07580213877225e+83,
    3.14699732603879e+85,
    1.98260831540444e+87,
    1.26886932185884e+89,
    8.24765059208247e+90,
    5.44344939077443e+92,
    3.64711109181887e+94,
    2.48003554243683e+96,
    1.71122452428141e+98,
    1.19785716699699e+100,
    8.50478588567862e+101,
    6.12344583768861e+103,
    4.47011546151268e+105,
    3.30788544151939e+107,
    2.48091408113954e+109,
    1.88549470166605e+111,
    1.45183092028286e+113,
    1.13242811782063e+115,
    8.94618213078297e+116,
    7.15694570462638e+118,
    5.79712602074737e+120,
    4.75364333701284e+122,
    3.94552396972066e+124,
    3.31424013456535e+126,
    2.81710411438055e+128,
    2.42270953836727e+130,
    2.10775729837953e+132,
    1.85482642257398e+134,
    1.65079551609085e+136,
    1.48571596448176e+138,
    1.3520015276784e+140,
    1.24384140546413e+142,
    1.15677250708164e+144,
    1.08736615665674e+146,
    1.03299784882391e+148,
    9.91677934870949e+149,
    9.61927596824821e+151,
    9.42689044888324e+153,
    9.33262154439441e+155,
    9.33262154439441e+157,
    9.42594775983835e+159,
    9.61446671503512e+161,
    9.90290071648618e+163,
    1.02990167451456e+166,
    1.08139675824029e+168,
    1.14628056373471e+170,
    1.22652020319614e+172,
    1.32464181945183e+174,
    1.44385958320249e+176,
    1.58824554152274e+178,
    1.76295255109024e+180,
    1.97450685722107e+182,
    2.23119274865981e+184,
    2.54355973347219e+186,
    2.92509369349301e+188,
    3.3931086844519e+190,
    3.96993716080872e+192,
    4.68452584975429e+194,
    5.5745857612076e+196,
    6.68950291344912e+198,
    8.09429852527344e+200,
    9.8750442008336e+202,
    1.21463043670253e+205,
    1.50614174151114e+207,
    1.88267717688893e+209,
    2.37217324288005e+211,
    3.01266001845766e+213,
    3.8562048236258e+215,
    4.97450422247729e+217,
    6.46685548922047e+219,
    8.47158069087882e+221,
    1.118248651196e+224,
    1.48727070609069e+226,
    1.99294274616152e+228,
    2.69047270731805e+230,
    3.65904288195255e+232,
    5.01288874827499e+234,
    6.91778647261949e+236,
    9.61572319694109e+238,
    1.34620124757175e+241,
    1.89814375907617e+243,
    2.69536413788816e+245,
    3.85437071718007e+247,
    5.5502938327393e+249,
    8.04792605747199e+251,
    1.17499720439091e+254,
    1.72724589045464e+256,
    2.55632391787286e+258,
    3.80892263763057e+260,
    5.71338395644585e+262,
    8.62720977423323e+264,
    1.31133588568345e+267,
    2.00634390509568e+269,
    3.08976961384735e+271,
    4.78914290146339e+273,
    7.47106292628289e+275,
    1.17295687942641e+278,
    1.85327186949373e+280,
    2.94670227249504e+282,
    4.71472363599206e+284,
    7.59070505394721e+286,
    1.22969421873945e+289,
    2.0044015765453e+291,
    3.28721858553429e+293,
    5.42391066613159e+295,
    9.00369170577843e+297,
    1.503616514865e+300,
  };

  if(n < 0 || n > nmax) {
    char str[128];
    sprintf(str, "Invalid argument to factorial %d", n);
    error->all(FLERR, str);
  }

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

double SNA::deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

void SNA::init_clebsch_gordan()
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= twojmax; j2++)
      for (int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2)
        for (int m1 = 0; m1 <= j1; m1 += 1) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2 += 1) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) continue;

	    sum = 0.0;

	    for (int z = MAX(0, MAX(-(j - j2 + aa2) 
				   / 2, -(j - j1 - bb2) / 2));
		z <= MIN((j1 + j2 - j) / 2,
			 MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
		z++) {
	      ifac = z % 2 ? -1 : 1;
	      sum += ifac /
		(factorial(z) *
		 factorial((j1 + j2 - j) / 2 - z) *
		 factorial((j1 - aa2) / 2 - z) *
		 factorial((j2 + bb2) / 2 - z) *
		 factorial((j - j2 + aa2) / 2 + z) *
		 factorial((j - j1 - bb2) / 2 + z));
	    }

	    cc2 = 2 * m - j;
	    dcg = deltacg(j1, j2, j);
	    sfaccg = sqrt(factorial((j1 + aa2) / 2) *
			factorial((j1 - aa2) / 2) *
			factorial((j2 + bb2) / 2) *
			factorial((j2 - bb2) / 2) *
			factorial((j  + cc2) / 2) *
			factorial((j  - cc2) / 2) *
			(j + 1));

	    cgarray[j1][j2][j][m1][m2] = sum * dcg * sfaccg;
	  }
	}
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax 
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

void SNA::init_rootpqarray()
{
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p][q] = sqrt(static_cast<double>(p)/q);
}

/* ----------------------------------------------------------------------
   a = j/2
------------------------------------------------------------------------- */

void SNA::jtostr(char* str, int j)
{
  if(j % 2 == 0)
    sprintf(str, "%d", j / 2);
  else
    sprintf(str, "%d/2", j);
}

/* ----------------------------------------------------------------------
   aa = m - j/2
------------------------------------------------------------------------- */

void SNA::mtostr(char* str, int j, int m)
{
  if(j % 2 == 0)
    sprintf(str, "%d", m - j / 2);
  else
    sprintf(str, "%d/2", 2 * m - j);
}

/* ----------------------------------------------------------------------
   list values of Clebsch-Gordan coefficients
   using notation of VMK Table 8.11
------------------------------------------------------------------------- */

void SNA::print_clebsch_gordan(FILE* file)
{
  char stra[20], strb[20], strc[20], straa[20], strbb[20], strcc[20];
  int m, aa2, bb2;

  fprintf(file, "a, aa, b, bb, c, cc, c(a,aa,b,bb,c,cc) \n");

  for (int j1 = 0; j1 <= twojmax; j1++) {
    jtostr(stra, j1);

    for (int j2 = 0; j2 <= twojmax; j2++) {
      jtostr(strb, j2);

      for (int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        jtostr(strc, j);

        for (int m1 = 0; m1 <= j1; m1 += 1) {
          mtostr(straa, j1, m1);
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2 += 1) {
            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) continue;

            mtostr(strbb, j2, m2);
            mtostr(strcc, j, m);

            fprintf(file, "%s\t%s\t%s\t%s\t%s\t%s\t%g\n",
                    stra, straa, strb, strbb, strc, strcc,
                    cgarray[j1][j2][j][m1][m2]);
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int SNA::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    if(diagonalstyle == 0) {
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          ncount++;
    } else if(diagonalstyle == 1) {
      int j2 = j1;

      for (int j = abs(j1 - j2);
          j <= MIN(twojmax, j1 + j2); j += 2)
        ncount++;
    } else if(diagonalstyle == 2) {
      ncount++;
    } else if(diagonalstyle == 3) {
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = abs(j1 - j2);
            j <= MIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) ncount++;
    }

  return ncount;
}

/* ---------------------------------------------------------------------- */

double SNA::compute_sfac(double r, double rcut)
{
  if (switch_flag == 0) return 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double SNA::compute_dsfac(double r, double rcut)
{
  if (switch_flag == 0) return 0.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      double rcutfac = MY_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}

