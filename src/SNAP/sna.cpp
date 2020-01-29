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
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#include "sna.h"
#include <cmath>
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

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

SNA::SNA(LAMMPS* lmp, double rfac0_in, int twojmax_in,
         double rmin0_in, int switch_flag_in, int bzero_flag_in,
         int alloy_flag_in, int wselfall_flag_in, int nelements_in) : Pointers(lmp)
{
  wself = 1.0;

  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;
  bzero_flag = bzero_flag_in;
  bnorm_flag = alloy_flag_in;
  alloy_flag = alloy_flag_in;
  wselfall_flag = wselfall_flag_in;
  nelements = nelements_in;

  twojmax = twojmax_in;

  compute_ncoeff();

  rij = NULL;
  inside = NULL;
  wj = NULL;
  rcutij = NULL;
  element = NULL;
  nmax = 0;
  idxz = NULL;
  idxb = NULL;
  ulist_r_ij = NULL;
  ulist_i_ij = NULL;

  build_indexlist();
  create_twojmax_arrays();

  if (bzero_flag) {
    double www = wself*wself*wself;
    for(int j = 0; j <= twojmax; j++)
      if (bnorm_flag)
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
  }

}

/* ---------------------------------------------------------------------- */

SNA::~SNA()
{
  memory->destroy(rij);
  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  delete[] idxz;
  delete[] idxb;
  destroy_twojmax_arrays();
}

void SNA::build_indexlist()
{

  // index list for cglist

  int jdim = twojmax + 1;
  memory->create(idxcg_block, jdim, jdim, jdim,
                 "sna:idxcg_block");

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j1][j2][j] = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;

  // index list for uarray
  // need to include both halves

  memory->create(idxu_block, jdim,
                 "sna:idxu_block");

  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;

  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
  idxb = new SNA_BINDICES[idxb_max];

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count].j1 = j1;
          idxb[idxb_count].j2 = j2;
          idxb[idxb_count].j = j;
          idxb_count++;
        }

  // reverse index list for beta and b

  memory->create(idxb_block, jdim, jdim, jdim,
                 "sna:idxb_block");
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j1][j2][j] = idxb_count;
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  idxz = new SNA_ZINDICES[idxz_max];

  memory->create(idxz_block, jdim, jdim, jdim,
                 "sna:idxz_block");

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j1][j2][j] = idxz_count;

        // find right beta[jjb] entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            idxz[idxz_count].j1 = j1;
            idxz[idxz_count].j2 = j2;
            idxz[idxz_count].j = j;
            idxz[idxz_count].ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count].ma2max = (2 * ma - j - (2 * idxz[idxz_count].ma1min - j1) + j2) / 2;
            idxz[idxz_count].na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count].ma1min + 1;
            idxz[idxz_count].mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count].mb2max = (2 * mb - j - (2 * idxz[idxz_count].mb1min - j1) + j2) / 2;
            idxz[idxz_count].nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count].mb1min + 1;
            idxz[idxz_count].ma = ma;
            idxz[idxz_count].mb = mb;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count].jju = jju;

            idxz_count++;
          }
      }
}

/* ---------------------------------------------------------------------- */

void SNA::init()
{
  init_clebsch_gordan();
  //   print_clebsch_gordan();
  init_rootpqarray();
}


void SNA::grow_rij(int newnmax)
{
  if(newnmax <= nmax) return;

  nmax = newnmax;

  memory->destroy(rij);
  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  memory->create(rij, nmax, 3, "pair:rij");
  memory->create(inside, nmax, "pair:inside");
  memory->create(wj, nmax, "pair:wj");
  memory->create(rcutij, nmax, "pair:rcutij");
  memory->create(element, nmax, "sna:element");
  memory->create(ulist_r_ij, nmax, idxu_max, "sna:ulist_ij");
  memory->create(ulist_i_ij, nmax, idxu_max, "sna:ulist_ij");
}

/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

void SNA::compute_ui(int jnum, int ielem)
{
  double rsq, r, x, y, z, z0, theta0;

  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  zero_uarraytot(ielem);
  addself_uarraytot(wself, ielem);

  for(int j = 0; j < jnum; j++) {
    x = rij[j][0];
    y = rij[j][1];
    z = rij[j][2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);

    compute_uarray(x, y, z, z0, r, j);
    add_uarraytot(r, wj[j], rcutij[j], j, element[j]);
  }

}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

void SNA::compute_zi()
{

  int idouble = 0;
  double * zptr_r;
  double * zptr_i;
  for(int elem1 = 0; elem1 < nelements; elem1++)
    for(int elem2 = 0; elem2 < nelements; elem2++) {

      zptr_r = &zlist_r[idouble*idxz_max];
      zptr_i = &zlist_i[idouble*idxz_max];

      for (int jjz = 0; jjz < idxz_max; jjz++) {
        const int j1 = idxz[jjz].j1;
        const int j2 = idxz[jjz].j2;
        const int j = idxz[jjz].j;
        const int ma1min = idxz[jjz].ma1min;
        const int ma2max = idxz[jjz].ma2max;
        const int na = idxz[jjz].na;
        const int mb1min = idxz[jjz].mb1min;
        const int mb2max = idxz[jjz].mb2max;
        const int nb = idxz[jjz].nb;

        const double *cgblock = cglist + idxcg_block[j1][j2][j];

        zptr_r[jjz] = 0.0;
        zptr_i[jjz] = 0.0;

        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;
        for (int ib = 0; ib < nb; ib++) {

          double suma1_r = 0.0;
          double suma1_i = 0.0;

          const double *u1_r = &ulisttot_r[elem1*idxu_max+jju1];
          const double *u1_i = &ulisttot_i[elem1*idxu_max+jju1];
          const double *u2_r = &ulisttot_r[elem2*idxu_max+jju2];
          const double *u2_i = &ulisttot_i[elem2*idxu_max+jju2];

          int ma1 = ma1min;
          int ma2 = ma2max;
          int icga = ma1min * (j2 + 1) + ma2max;

          for (int ia = 0; ia < na; ia++) {
            suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
            suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
            ma1++;
            ma2--;
            icga += j2;
          } // end loop over ia

          if (bnorm_flag){
            zptr_r[jjz] += cgblock[icgb] * suma1_r/(j+1);
            zptr_i[jjz] += cgblock[icgb] * suma1_i/(j+1);
          }
          else {
            zptr_r[jjz] += cgblock[icgb] * suma1_r;
            zptr_i[jjz] += cgblock[icgb] * suma1_i;
          }
          jju1 += j1 + 1;
          jju2 -= j2 + 1;
          icgb += j2;
        } // end loop over ib
      } // end loop over jjz
      idouble++;
    }
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices
------------------------------------------------------------------------- */

void SNA::compute_yi(const double* beta)
{
  int jju;
  double betaj;
  int itriple;
  int jelem;
  double temp;

  for(int ielem1 = 0; ielem1 < nelements; ielem1++)
    for(int j = 0; j <= twojmax; j++) {
      jju = idxu_block[j];
      for(int mb = 0; 2*mb <= j; mb++)
        for(int ma = 0; ma <= j; ma++) {
          ylist_r[ielem1*idxu_max+jju] = 0.0;
          ylist_i[ielem1*idxu_max+jju] = 0.0;
          jju++;
        } // end loop over ma, mb
    } // end loop over j

  for(int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {
        for (int jjz = 0; jjz < idxz_max; jjz++) {
          const int j1 = idxz[jjz].j1;
          const int j2 = idxz[jjz].j2;
          const int j = idxz[jjz].j;
          const int ma1min = idxz[jjz].ma1min;
          const int ma2max = idxz[jjz].ma2max;
          const int na = idxz[jjz].na;
          const int mb1min = idxz[jjz].mb1min;
          const int mb2max = idxz[jjz].mb2max;
          const int nb = idxz[jjz].nb;
          const int ma = idxz[jjz].ma;
          const int mb = idxz[jjz].mb;

          const double *cgblock = cglist + idxcg_block[j1][j2][j];

          double ztmp_r = 0.0;
          double ztmp_i = 0.0;

          int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
          int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
          int icgb = mb1min * (j2 + 1) + mb2max;
          for (int ib = 0; ib < nb; ib++) {

            double suma1_r = 0.0;
            double suma1_i = 0.0;

            const double *u1_r = &ulisttot_r[elem1*idxu_max+jju1];
            const double *u1_i = &ulisttot_i[elem1*idxu_max+jju1];
            const double *u2_r = &ulisttot_r[elem2*idxu_max+jju2];
            const double *u2_i = &ulisttot_i[elem2*idxu_max+jju2];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
              suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
              suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
              ma1++;
              ma2--;
              icga += j2;
            } // end loop over ia

            if (bnorm_flag) {
              ztmp_r += cgblock[icgb] * suma1_r / (j+1);
              ztmp_i += cgblock[icgb] * suma1_i / (j+1);
            } else {
              ztmp_r += cgblock[icgb] * suma1_r;
              ztmp_i += cgblock[icgb] * suma1_i;
            }

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib

          // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
          // find right y_list[jju] and beta[jjb] entries
          // multiply and divide by j+1 factors
          // account for multiplicity of 1, 2, or 3

          jju = idxz[jjz].jju;
          for(int elem3 = 0; elem3 < nelements; elem3++) {
          // pick out right beta value
          if (alloy_flag) {
            if (j >= j1) {
              const int jjb = idxb_block[j1][j2][j];
              itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
              if (j1 ==j && j2 == j && elem1 == elem2 && elem1 == elem3) betaj = 3 * beta[itriple];
              else if (j1 == j && j2 == j && elem2 == elem3 && elem1!=elem3)
                  betaj =
                      2 * (beta[itriple] + beta[((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb]);
              else if (j1 == j && j2 == j && elem1 == elem2 && elem1 != elem3)
                betaj = beta[itriple] + beta[((elem2 * nelements + elem1) * nelements + elem3) * idxb_max + jjb] +
                        beta[((elem1 * nelements + elem3) * nelements + elem2) * idxb_max + jjb];
              else if (j1 ==j && (elem1 == elem3 || elem2 == elem3)) // this line covers quite a few cases
                betaj = 2 * beta[itriple];
              else if (j1 == j && j2 != j && elem2 == elem1 && elem1 != elem3)
                betaj = beta[itriple] + beta[((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb];
              else
                betaj = beta[((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb];
            } else if (j >= j2) {
              const int jjb = idxb_block[j][j2][j1];
              itriple = ((elem3 * nelements + elem2) * nelements + elem1)*idxb_max + jjb;
              if (j2 == j) {
                if (elem3 == elem2)
                  betaj = 2 * beta[itriple];
                else if (elem1 == elem2)
                  betaj = beta[itriple] + beta[((elem1 * nelements + elem3) * nelements + elem2) * idxb_max + jjb];
                else if (elem1 == elem3)
                  betaj = beta[itriple] + beta[((elem2 * nelements + elem1) * nelements + elem3) * idxb_max + jjb];
              } else
                betaj = beta[itriple];
            } else {
              const int jjb = idxb_block[j2][j][j1];
              itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
              betaj = beta[itriple];
            }
          } else {
            if (j >= j1) {
              const int jjb = idxb_block[j1][j2][j];
              if (j1 == j) {
                if (j2 == j) betaj = 3*beta[jjb];
                else betaj = 2*beta[jjb];
              } else betaj = beta[jjb];
            } else if (j >= j2) {
              const int jjb = idxb_block[j][j2][j1];
              if (j2 == j) betaj = 2*beta[jjb];
              else betaj = beta[jjb];
            } else {
              const int jjb = idxb_block[j2][j][j1];
              betaj = beta[jjb];
            }
          }

          if (!bnorm_flag && j1 > j)
            betaj *= (j1 + 1) / (j + 1.0);

          ylist_r[elem3 * idxu_max + jju] += betaj * ztmp_r;
          ylist_i[elem3 * idxu_max + jju] += betaj * ztmp_i;
          //if (elem3==0 && j ==4 &&  ma==2 &&  mb==0)
          //  fprintf(screen, "%i %i %i %i %i %i %f %f %f\n", elem1, elem2, elem3, j, j1, j2, betaj, ztmp_r, ylist_r[elem3 * idxu_max + jju]);

        }
      } // end loop over jjz
    }

}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

void SNA::compute_deidrj(double* dedr)
{

  for(int k = 0; k < 3; k++)
    dedr[k] = 0.0;

  int ielem = elem_duarray;
  for(int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {

        double* dudr_r = dulist_r[jju];
        double* dudr_i = dulist_i[jju];
        double jjjmambyarray_r = ylist_r[ielem*idxu_max+jju];
        double jjjmambyarray_i = ylist_i[ielem*idxu_max+jju];

        for(int k = 0; k < 3; k++)
          dedr[k] +=
            dudr_r[k] * jjjmambyarray_r +
            dudr_i[k] * jjjmambyarray_i;
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        double* dudr_r = dulist_r[jju];
        double* dudr_i = dulist_i[jju];
        double jjjmambyarray_r = ylist_r[ielem*idxu_max+jju];
        double jjjmambyarray_i = ylist_i[ielem*idxu_max+jju];

        for(int k = 0; k < 3; k++)
          dedr[k] +=
            dudr_r[k] * jjjmambyarray_r +
            dudr_i[k] * jjjmambyarray_i;
        jju++;
      }

      double* dudr_r = dulist_r[jju];
      double* dudr_i = dulist_i[jju];
      double jjjmambyarray_r = ylist_r[ielem*idxu_max+jju];
      double jjjmambyarray_i = ylist_i[ielem*idxu_max+jju];

      for(int k = 0; k < 3; k++)
        dedr[k] +=
          (dudr_r[k] * jjjmambyarray_r +
           dudr_i[k] * jjjmambyarray_i)*0.5;
      // jju++;

    } // end if jeven

  } // end loop over j

  for(int k = 0; k < 3; k++)
    dedr[k] *= 2.0;

}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
------------------------------------------------------------------------- */

void SNA::compute_bi(int ielem) {
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  int itriple = 0;
  int idouble = 0;
  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {

      double *zptr_r = &zlist_r[idouble*idxz_max];
      double *zptr_i = &zlist_i[idouble*idxz_max];

      for (int elem3 = 0; elem3 < nelements; elem3++) {
        for (int jjb = 0; jjb < idxb_max; jjb++) {
          const int j1 = idxb[jjb].j1;
          const int j2 = idxb[jjb].j2;
          const int j = idxb[jjb].j;

          int jjz = idxz_block[j1][j2][j];
          int jju = idxu_block[j];
          double sumzu = 0.0;
          for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju] * zptr_r[jjz] +
                       ulisttot_i[elem3*idxu_max+jju] * zptr_i[jjz];
              jjz++;
              jju++;
            } // end loop over ma, mb

          // For j even, handle middle column

          if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju] * zptr_r[jjz] +
                       ulisttot_i[elem3*idxu_max+jju] * zptr_i[jjz];
              jjz++;
              jju++;
            }

            sumzu += 0.5 * (ulisttot_r[elem3*idxu_max+jju] * zptr_r[jjz] +
                            ulisttot_i[elem3*idxu_max+jju] * zptr_i[jjz]);
          } // end if jeven

          blist[itriple*idxb_max+jjb] = 2.0 * sumzu;

        }
        itriple++;
      }
      idouble++;
    }

  // apply bzero shift

  if (bzero_flag) {
    if (!wselfall_flag) {
      itriple = (ielem*nelements+ielem)*nelements+ielem;
      for (int jjb = 0; jjb < idxb_max; jjb++) {
        const int j = idxb[jjb].j;
        blist[itriple*idxb_max+jjb] -= bzero[j];
      } // end loop over JJ
    } else {
      int itriple = 0;
      for(int elem1 = 0; elem1 < nelements; elem1++)
        for(int elem2 = 0; elem2 < nelements; elem2++) {
          for(int elem3 = 0; elem3 < nelements; elem3++) {
            for (int jjb = 0; jjb < idxb_max; jjb++) {
              const int j = idxb[jjb].j;
              blist[itriple*idxb_max+jjb] -= bzero[j];
            } // end loop over JJ
            itriple++;
          } // end loop over elem3
        } // end loop over elem1,elem2
    }
  }
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
  //              Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
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
  double* zptr_r;
  double* zptr_i;
  int jjz, jju;

  int idouble;
  int itriple;

  // set all the derivatives to zero once

  for(int jjb = 0; jjb < idxb_max; jjb++) {

    for(int elem1 = 0; elem1 < nelements; elem1++)
      for(int elem2 = 0; elem2 < nelements; elem2++)
        for(int elem3 = 0; elem3 < nelements; elem3++) {

          itriple = (elem1 * nelements + elem2) * nelements + elem3;

          dbdr = dblist[itriple*idxb_max+jjb];
          dbdr[0] = 0.0;
          dbdr[1] = 0.0;
          dbdr[2] = 0.0;
        }

  }

  int elem3 = elem_duarray;

  for(int jjb = 0; jjb < idxb_max; jjb++) {
    const int j1 = idxb[jjb].j1;
    const int j2 = idxb[jjb].j2;
    const int j = idxb[jjb].j;


    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    for(int elem1 = 0; elem1 < nelements; elem1++)
      for(int elem2 = 0; elem2 < nelements; elem2++) {

        jjz = idxz_block[j1][j2][j];
        jju = idxu_block[j];
        idouble = elem1*nelements+elem2;
        itriple = (elem1*nelements+elem2)*nelements+elem3;
        dbdr = dblist[itriple*idxb_max+jjb];
        zptr_r = &zlist_r[idouble*idxz_max];
        zptr_i = &zlist_i[idouble*idxz_max];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j even, handle middle column

        if (j % 2 == 0) {
          int mb = j / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = dulist_r[jju];
          dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if jeven

        for (int k = 0; k < 3; k++)
          dbdr[k] += 2.0 * sumzdu_r[k];
        // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

        double j1fac = (j + 1) / (j1 + 1.0);

        idouble = elem1*nelements+elem2;
        itriple = (elem3*nelements+elem2)*nelements+elem1;
        dbdr = dblist[itriple*idxb_max+jjb];
        jjz = idxz_block[j][j2][j1];
        jju = idxu_block[j1];
        zptr_r = &zlist_r[idouble*idxz_max];
        zptr_i = &zlist_i[idouble*idxz_max];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j1; mb++)
          for (int ma = 0; ma <= j1; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j1 even, handle middle column

        if (j1 % 2 == 0) {
          int mb = j1 / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = dulist_r[jju];
          dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j1even

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
            dbdr[k] += 2.0 * sumzdu_r[k];
          else
            dbdr[k] += 2.0 * sumzdu_r[k] * j1fac;

        // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)

        double j2fac = (j + 1) / (j2 + 1.0);

        idouble = elem2*nelements+elem1;
        itriple = (elem1*nelements+elem3)*nelements+elem2;
        dbdr = dblist[itriple*idxb_max+jjb];
        jjz = idxz_block[j][j1][j2];
        jju = idxu_block[j2];
        zptr_r = &zlist_r[idouble*idxz_max];
        zptr_i = &zlist_i[idouble*idxz_max];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j2; mb++)
          for (int ma = 0; ma <= j2; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j2 even, handle middle column

        if (j2 % 2 == 0) {
          int mb = j2 / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = dulist_r[jju];
            dudr_i = dulist_i[jju];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = dulist_r[jju];
          dudr_i = dulist_i[jju];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j2even

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
            dbdr[k] += 2.0 * sumzdu_r[k];
          else
            dbdr[k] += 2.0 * sumzdu_r[k] * j2fac;

      }
  } //end loop over j1 j2 j

}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

void SNA::compute_duidrj(double* rij, double wj, double rcut, int jj, int jelem)
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

  elem_duarray = jelem;
  compute_duarray(x, y, z, z0, r, dz0dr, wj, rcut, jj);
}

/* ---------------------------------------------------------------------- */

void SNA::zero_uarraytot(int ielem)
{
  for(int jelem = 0; jelem < nelements; jelem++)
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        ulisttot_r[jelem*idxu_max+jju] = 0.0;
        ulisttot_i[jelem*idxu_max+jju] = 0.0;

        // utot(j,ma,ma) = wself, sometimes
        if (jelem == ielem || wselfall_flag)
          if (ma==mb)
          ulisttot_r[jelem*idxu_max+jju] = wself; ///// double check this
        jju++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void SNA::addself_uarraytot(double wself_in, int ielem)
{
  // for (int j = 0; j <= twojmax; j++) {
  //   int jju = idxu_block[j];
  //   for (int ma = 0; ma <= j; ma++) {
  //     ulisttot_r[ielem*idxu_max+jju] = wself_in;
  //     ulisttot_i[ielem*idxu_max+jju] = 0.0;
  //     jju += j+2;
  //   }
  // }
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

void SNA::add_uarraytot(double r, double wj, double rcut, int jj, int ielem)
{
  double sfac;

  sfac = compute_sfac(r, rcut);

  sfac *= wj;

  double* ulist_r = ulist_r_ij[jj];
  double* ulist_i = ulist_i_ij[jj];

  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        ulisttot_r[ielem*idxu_max+jju] +=
          sfac * ulist_r[jju];
        ulisttot_i[ielem*idxu_max+jju] +=
          sfac * ulist_i[jju];
        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

void SNA::compute_uarray(double x, double y, double z,
                         double z0, double r, int jj)
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


  double* ulist_r = ulist_r_ij[jj];
  double* ulist_i = ulist_i_ij[jj];

  ulist_r[0] = 1.0;
  ulist_i[0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    for (int mb = 0; 2*mb <= j; mb++) {
      ulist_r[jju] = 0.0;
      ulist_i[jju] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[j - ma][j - mb];
        ulist_r[jju] +=
          rootpq *
          (a_r * ulist_r[jjup] +
           a_i * ulist_i[jjup]);
        ulist_i[jju] +=
          rootpq *
          (a_r * ulist_i[jjup] -
           a_i * ulist_r[jjup]);

        rootpq = rootpqarray[ma + 1][j - mb];
        ulist_r[jju+1] =
          -rootpq *
          (b_r * ulist_r[jjup] +
           b_i * ulist_i[jjup]);
        ulist_i[jju+1] =
          -rootpq *
          (b_r * ulist_i[jjup] -
           b_i * ulist_r[jjup]);
        jju++;
        jjup++;
      }
      jju++;
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        if (mapar == 1) {
          ulist_r[jjup] = ulist_r[jju];
          ulist_i[jjup] = -ulist_i[jju];
        } else {
          ulist_r[jjup] = -ulist_r[jju];
          ulist_i[jjup] = ulist_i[jju];
        }
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  }
}

/* ----------------------------------------------------------------------
   Compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

void SNA::compute_duarray(double x, double y, double z,
                          double z0, double r, double dz0dr,
                          double wj, double rcut, int jj)
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

  double* ulist_r = ulist_r_ij[jj];
  double* ulist_i = ulist_i_ij[jj];

  dulist_r[0][0] = 0.0;
  dulist_r[0][1] = 0.0;
  dulist_r[0][2] = 0.0;
  dulist_i[0][0] = 0.0;
  dulist_i[0][1] = 0.0;
  dulist_i[0][2] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    for (int mb = 0; 2*mb <= j; mb++) {
      dulist_r[jju][0] = 0.0;
      dulist_r[jju][1] = 0.0;
      dulist_r[jju][2] = 0.0;
      dulist_i[jju][0] = 0.0;
      dulist_i[jju][1] = 0.0;
      dulist_i[jju][2] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[j - ma][j - mb];
        for (int k = 0; k < 3; k++) {
          dulist_r[jju][k] +=
            rootpq * (da_r[k] * ulist_r[jjup] +
                      da_i[k] * ulist_i[jjup] +
                      a_r * dulist_r[jjup][k] +
                      a_i * dulist_i[jjup][k]);
          dulist_i[jju][k] +=
            rootpq * (da_r[k] * ulist_i[jjup] -
                      da_i[k] * ulist_r[jjup] +
                      a_r * dulist_i[jjup][k] -
                      a_i * dulist_r[jjup][k]);
        }

        rootpq = rootpqarray[ma + 1][j - mb];
        for (int k = 0; k < 3; k++) {
          dulist_r[jju+1][k] =
            -rootpq * (db_r[k] * ulist_r[jjup] +
                       db_i[k] * ulist_i[jjup] +
                       b_r * dulist_r[jjup][k] +
                       b_i * dulist_i[jjup][k]);
          dulist_i[jju+1][k] =
            -rootpq * (db_r[k] * ulist_i[jjup] -
                       db_i[k] * ulist_r[jjup] +
                       b_r * dulist_i[jjup][k] -
                       b_i * dulist_r[jjup][k]);
        }
        jju++;
        jjup++;
      }
      jju++;
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju+(j+1)*(j+1)-1;
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        if (mapar == 1) {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup][k] = dulist_r[jju][k];
            dulist_i[jjup][k] = -dulist_i[jju][k];
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup][k] = -dulist_r[jju][k];
            dulist_i[jjup][k] = dulist_i[jju][k];
          }
        }
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  }

  double sfac = compute_sfac(r, rcut);
  double dsfac = compute_dsfac(r, rcut);

  sfac *= wj;
  dsfac *= wj;
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist_r[jju][0] = dsfac * ulist_r[jju] * ux +
                                  sfac * dulist_r[jju][0];
        dulist_i[jju][0] = dsfac * ulist_i[jju] * ux +
                                  sfac * dulist_i[jju][0];
        dulist_r[jju][1] = dsfac * ulist_r[jju] * uy +
                                  sfac * dulist_r[jju][1];
        dulist_i[jju][1] = dsfac * ulist_i[jju] * uy +
                                  sfac * dulist_i[jju][1];
        dulist_r[jju][2] = dsfac * ulist_r[jju] * uz +
                                  sfac * dulist_r[jju][2];
        dulist_i[jju][2] = dsfac * ulist_i[jju] * uz +
                                  sfac * dulist_i[jju][2];
        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

double SNA::memory_usage()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  double bytes;

  bytes = 0;

  bytes += jdimpq*jdimpq * sizeof(double);               // pqarray
  bytes += idxcg_max * sizeof(double);                   // cglist

  bytes += nmax * idxu_max * sizeof(double) * 2;         // ulist_ij
  bytes += idxu_max * nelements * sizeof(double) * 2;    // ulisttot
  bytes += idxu_max * 3 * sizeof(double) * 2;            // dulist

  bytes += idxz_max * ndoubles * sizeof(double) * 2;     // zlist
  bytes += idxb_max * ntriples * sizeof(double);         // blist
  bytes += idxb_max * ntriples * 3 * sizeof(double);     // dblist
  bytes += idxu_max * nelements * sizeof(double) * 2;    // ylist

  bytes += jdim * jdim * jdim * sizeof(int);             // idxcg_block
  bytes += jdim * sizeof(int);                           // idxu_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxz_block
  bytes += jdim * jdim * jdim * sizeof(int);             // idxb_block

  bytes += idxz_max * sizeof(SNA_ZINDICES);              // idxz
  bytes += idxb_max * sizeof(SNA_BINDICES);              // idxb

  if (bzero_flag)
  bytes += jdim * sizeof(double);                        // bzero

  bytes += nmax * 3 * sizeof(double);                    // rij
  bytes += nmax * sizeof(int);                           // inside
  bytes += nmax * sizeof(double);                        // wj
  bytes += nmax * sizeof(double);                        // rcutij

  return bytes;
}
/* ---------------------------------------------------------------------- */

void SNA::create_twojmax_arrays()
{
  int jdimpq = twojmax + 2;
  memory->create(rootpqarray, jdimpq, jdimpq,
                 "sna:rootpqarray");
  memory->create(cglist, idxcg_max, "sna:cglist");
  memory->create(ulisttot_r, idxu_max*nelements, "sna:ulisttot");
  memory->create(ulisttot_i, idxu_max*nelements, "sna:ulisttot");
  memory->create(dulist_r, idxu_max, 3, "sna:dulist");
  memory->create(dulist_i, idxu_max, 3, "sna:dulist");
  memory->create(zlist_r, idxz_max*ndoubles, "sna:zlist");
  memory->create(zlist_i, idxz_max*ndoubles, "sna:zlist");
  memory->create(blist, idxb_max*ntriples, "sna:blist");
  memory->create(dblist, idxb_max*ntriples, 3, "sna:dblist");
  memory->create(ylist_r, idxu_max*nelements, "sna:ylist");
  memory->create(ylist_i, idxu_max*nelements, "sna:ylist");

  if (bzero_flag)
    memory->create(bzero, twojmax+1,"sna:bzero");
  else
    bzero = NULL;

}

/* ---------------------------------------------------------------------- */

void SNA::destroy_twojmax_arrays()
{
  memory->destroy(rootpqarray);
  memory->destroy(cglist);
  memory->destroy(ulisttot_r);
  memory->destroy(ulisttot_i);
  memory->destroy(dulist_r);
  memory->destroy(dulist_i);
  memory->destroy(zlist_r);
  memory->destroy(zlist_i);
  memory->destroy(blist);
  memory->destroy(dblist);
  memory->destroy(ylist_r);
  memory->destroy(ylist_i);

  memory->destroy(idxcg_block);
  memory->destroy(idxu_block);
  memory->destroy(idxz_block);
  memory->destroy(idxb_block);

  if (bzero_flag)
    memory->destroy(bzero);

}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

double SNA::factorial(int n)
{
  if (n < 0 || n > nmaxfactorial) {
    char str[128];
    sprintf(str, "Invalid argument to factorial %d", n);
    error->all(FLERR, str);
  }

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

const double SNA::nfac_table[] = {
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
  1.503616514865e+300, // nmaxfactorial = 167
};

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

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

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

            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}

/* ----------------------------------------------------------------------
   print out values of Clebsch-Gordan coefficients
   format and notation follows VMK Table 8.11
------------------------------------------------------------------------- */

void SNA::print_clebsch_gordan()
{
  if (comm->me) return;

  int aa2, bb2, cc2;
  for (int j = 0; j <= twojmax; j += 1) {
    printf("c = %g\n",j/2.0);
    printf("a alpha b beta C_{a alpha b beta}^{c alpha+beta}\n");
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        if (j1-j2 <= j && j1+j2 >= j && (j1+j2+j)%2 == 0) {
          int idxcg_count = idxcg_block[j1][j2][j];
          for (int m1 = 0; m1 <= j1; m1++) {
            aa2 = 2*m1-j1;
            for (int m2 = 0; m2 <= j2; m2++) {
              bb2 = 2*m2-j2;
              double cgtmp = cglist[idxcg_count];
              cc2 = aa2+bb2;
              if (cc2 >= -j && cc2 <= j)
                if (j1 != j2 || (aa2 > bb2 && aa2 >= -bb2) || (aa2 == bb2 && aa2 >= 0))
                  printf("%4g %4g %4g %4g %10.6g\n",
                         j1/2.0,aa2/2.0,j2/2.0,bb2/2.0,cgtmp);
              idxcg_count++;
            }
          }
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

/* ---------------------------------------------------------------------- */

void SNA::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  ndoubles = nelements*nelements;
  ntriples = nelements*nelements*nelements;
  if (alloy_flag)
    ncoeff = ncount*ntriples;
  else
    ncoeff = ncount;
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

