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

/* ----------------------------------------------------------------------
   Contributing authors: W. Michael Brown, Intel
------------------------------------------------------------------------- */

#if defined(__AVX512F__)
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)

#include "sna_intel.h"

#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"

#include <cmath>

using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;
using namespace ip_simd;

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
   Doctoral Thesis, Cambridge University, (2009)

   [2] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
   "Quantum Theory of Angular Momentum," World Scientific (1988)

------------------------------------------------------------------------- */

SNAIntel::SNAIntel(LAMMPS* lmp, double rfac0_in, int twojmax_in,
                   double rmin0_in, int switch_flag_in, int bzero_flag_in,
                   int chem_flag_in, int bnorm_flag_in, int wselfall_flag_in,
                   int nelements_in, int switch_inner_flag_in) : Pointers(lmp)
{
  wself = 1.0;

  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;
  switch_inner_flag = switch_inner_flag_in;
  bzero_flag = bzero_flag_in;
  chem_flag = chem_flag_in;
  bnorm_flag = bnorm_flag_in;
  wselfall_flag = wselfall_flag_in;

  if (bnorm_flag != chem_flag)
    lmp->error->warning(FLERR, "bnormflag and chemflag are not equal."
                        "This is probably not what you intended");

  if (chem_flag)
    nelements = nelements_in;
  else
    nelements = 1;

  twojmax = twojmax_in;

  compute_ncoeff();

  rij = nullptr;
  inside = nullptr;
  wj = nullptr;
  rcutij = nullptr;
  sinnerij = nullptr;
  dinnerij = nullptr;
  element = nullptr;
  nmax = 0;
  idxz = nullptr;
  idxb = nullptr;
  ulist_r_ij = nullptr;
  ulist_i_ij = nullptr;

  build_indexlist();
  create_twojmax_arrays();

  if (bzero_flag) {
    double www = wself*wself*wself;
    for (int j = 0; j <= twojmax; j++)
      if (bnorm_flag)
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
  }

}

/* ---------------------------------------------------------------------- */

SNAIntel::~SNAIntel()
{
  memory->destroy(rij);
  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(sinnerij);
  memory->destroy(dinnerij);
  if (chem_flag) memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  delete[] idxz;
  delete[] idxb;
  destroy_twojmax_arrays();
}

void SNAIntel::build_indexlist()
{

  // index list for cglist

  int jdim = twojmax + 1;
  memory->create(idxcg_block, jdim, jdim, jdim,
                 "sna:idxcg_block");

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
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

  for (int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;

  // index list for beta and B

  int idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
  idxb = new SNA_BINDICES[idxb_max];

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
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
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j1][j2][j] = idxb_count;
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  idxz = new SNA_ZINDICES[idxz_max];

  memory->create(idxz_block, jdim, jdim, jdim,
                 "sna:idxz_block");

  idxz_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
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
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count].jju = jju;

            idxz_count++;
          }
      }
}

/* ---------------------------------------------------------------------- */

void SNAIntel::init()
{
  init_clebsch_gordan();
  //   print_clebsch_gordan();
  init_rootpqarray();
}

void SNAIntel::grow_rij(int newnmax)
{
  if (newnmax <= nmax) return;

  nmax = newnmax;

  memory->destroy(rij);
  memory->destroy(inside);
  memory->destroy(wj);
  memory->destroy(rcutij);
  memory->destroy(sinnerij);
  memory->destroy(dinnerij);
  if (chem_flag) memory->destroy(element);
  memory->destroy(ulist_r_ij);
  memory->destroy(ulist_i_ij);
  memory->create(rij, nmax, 3, "pair:rij");
  memory->create(inside, nmax, "pair:inside");
  memory->create(wj, nmax, "pair:wj");
  memory->create(rcutij, nmax, "pair:rcutij");
  memory->create(sinnerij, nmax, "pair:sinnerij");
  memory->create(dinnerij, nmax, "pair:dinnerij");
  if (chem_flag) memory->create(element, nmax, "sna:element");
  memory->create(ulist_r_ij, nmax, idxu_max, "sna:ulist_ij");
  memory->create(ulist_i_ij, nmax, idxu_max, "sna:ulist_ij");
}

/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

void SNAIntel::compute_ui(const SNA_IVEC &jnum, const SNA_IVEC &ielem,
                          const int max_jnum)
{
  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  zero_uarraytot(ielem);

  for (int j = 0; j < max_jnum; j++) {
    const SNA_DVEC x = rij[j][0];
    const SNA_DVEC y = rij[j][1];
    const SNA_DVEC z = rij[j][2];
    const SNA_DVEC rcut = rcutij[j];
    const SNA_DVEC rsq = x * x + y * y + z * z;
    const SNA_DVEC r = SIMD_sqrt(rsq);
    const SNA_DVEC rscale0 = SIMD_rcp(rcut - rmin0) * rfac0 * MY_PI;
    const SNA_DVEC theta0 = (r - rmin0) * rscale0;
    const SNA_DVEC z0 = r * SIMD_rcp(SIMD_tan(theta0));

    compute_uarray(x, y, z, z0, r, j, jnum);
    add_uarraytot(r, j, jnum);
  }

}

/* ----------------------------------------------------------------------
   pick out right beta value
------------------------------------------------------------------------- */

double SNAIntel::choose_beta(const int j, const int j1, const int j2,
                             const int elem1, const int elem2, const int elem3,
                             int &itriple)
{
  double bfactor;
  if (j >= j1) {
    const int jjb = idxb_block[j1][j2][j];
    itriple = ((elem1 * nelements + elem2) * nelements + elem3) *
      idxb_max + jjb;
    if (j1 == j) {
      if (j2 == j)
        bfactor = 3.0;
      else
        bfactor = 2.0;
    } else
      bfactor = 1.0;
  } else if (j >= j2) {
    const int jjb = idxb_block[j][j2][j1];
    itriple = ((elem3 * nelements + elem2) * nelements + elem1) *
      idxb_max + jjb;
    if (j2 == j)
      bfactor = 2.0;
    else
      bfactor = 1.0;
  } else {
    const int jjb = idxb_block[j2][j][j1];
    itriple = ((elem2 * nelements + elem3) * nelements + elem1) *
      idxb_max + jjb;
    bfactor = 1.0;
  }

  if (!bnorm_flag && j1 > j)
    bfactor *= (1.0 + j1) / (1.0 + j);

  return bfactor;
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist indices
------------------------------------------------------------------------- */

template <int COMPUTE_YI>
void SNAIntel::compute_zi_or_yi(const SNA_DVEC* beta)
{
  if (COMPUTE_YI) {
    memset(ylist_r,0,idxu_max*nelements*sizeof(SNA_DVEC));
    memset(ylist_i,0,idxu_max*nelements*sizeof(SNA_DVEC));
  }

  double *zlist_rp = (double *)zlist_r;
  double *zlist_ip = (double *)zlist_i;

  int zlist_i = 0;

  for (int elem1 = 0; elem1 < nelements; elem1++)
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

        const double *cgblock = cglist + idxcg_block[j1][j2][j];

        SNA_DVEC ztmp_r = 0.0;
        SNA_DVEC ztmp_i = 0.0;

        const double *u_r = (double *)ulisttot_r;
        const double *u_i = (double *)ulisttot_i;

        int jju1 = elem1 * idxu_max + idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = elem2 * idxu_max + idxu_block[j2] + (j2 + 1) * mb2max;
        jju1 *= vector_width();
        jju2 *= vector_width();
        int icgb = mb1min * (j2 + 1) + mb2max;
        for (int ib = 0; ib < nb; ib++) {

          SNA_DVEC suma1_r = 0.0;
          SNA_DVEC suma1_i = 0.0;

          int ma1 = ma1min * vector_width();
          int ma2 = ma2max * vector_width();
          int icga = ma1min * (j2 + 1) + ma2max;

          for (int ia = 0; ia < na; ia++) {
            const SNA_DVEC u1_r = SIMD_load(u_r + jju1 + ma1);
            const SNA_DVEC u2_r = SIMD_load(u_r + jju2 + ma2);
            const SNA_DVEC u1_i = SIMD_load(u_i + jju1 + ma1);
            const SNA_DVEC u2_i = SIMD_load(u_i + jju2 + ma2);
            suma1_r += (u1_r*u2_r - u1_i*u2_i) * cgblock[icga];
            suma1_i += (u1_r*u2_i + u1_i*u2_r) * cgblock[icga];
            ma1+= vector_width();
            ma2-= vector_width();
            icga += j2;
          } // end loop over ia

          ztmp_r += suma1_r * cgblock[icgb];
          ztmp_i += suma1_i * cgblock[icgb];

          jju1 += (j1 + 1) * vector_width();
          jju2 -= (j2 + 1) * vector_width();
          icgb += j2;
        } // end loop over ib

        // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
        // find right y_list[jju] and beta[jjb] entries
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        if (bnorm_flag) {
          ztmp_i *= SIMD_rcp(SIMD_set(static_cast<double>(j+1)));
          ztmp_r *= SIMD_rcp(SIMD_set(static_cast<double>(j+1)));
        }

        if (COMPUTE_YI) {
          int jju = idxz[jjz].jju;
          for (int elem3 = 0; elem3 < nelements; elem3++) {
            int itriple;
            double bfactor = choose_beta(j, j1, j2, elem1, elem2, elem3,
                                         itriple);
            const SNA_DVEC betaj = beta[itriple] * bfactor;
            const int i = elem3 * idxu_max + jju;
            SIMD_store(&(ylist_r[i]), SIMD_load(ylist_r + i) + betaj * ztmp_r);
            SIMD_store(&(ylist_i[i]), SIMD_load(ylist_i + i) + betaj * ztmp_i);
          }
        } else {
          SIMD_store(zlist_rp + zlist_i, ztmp_r);
          SIMD_store(zlist_ip + zlist_i, ztmp_i);
          zlist_i += vector_width();
        }
      }// end loop over jjz
    }
}

/* ----------------------------------------------------------------------
   compute Yi from Zi
------------------------------------------------------------------------- */

void SNAIntel::compute_yi_from_zi(const SNA_DVEC* beta)
{
  memset(ylist_r,0,idxu_max*nelements*sizeof(SNA_DVEC));
  memset(ylist_i,0,idxu_max*nelements*sizeof(SNA_DVEC));

  double *zlist_rp = (double *)zlist_r;
  double *zlist_ip = (double *)zlist_i;

  int zlist_i = 0;

  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {
      for (int jjz = 0; jjz < idxz_max; jjz++) {
        const int j1 = idxz[jjz].j1;
        const int j2 = idxz[jjz].j2;
        const int j = idxz[jjz].j;

        const SNA_DVEC ztmp_r = SIMD_load(zlist_rp + zlist_i);
        const SNA_DVEC ztmp_i = SIMD_load(zlist_ip + zlist_i);
        zlist_i += vector_width();

        int jju = idxz[jjz].jju;
        for (int elem3 = 0; elem3 < nelements; elem3++) {
          int itriple;
          double bfactor = choose_beta(j, j1, j2, elem1, elem2, elem3,
                                       itriple);
          const SNA_DVEC betaj = beta[itriple] * bfactor;
          const int i = elem3 * idxu_max + jju;
          SIMD_store(&(ylist_r[i]), SIMD_load(ylist_r + i) + betaj * ztmp_r);
          SIMD_store(&(ylist_i[i]), SIMD_load(ylist_i + i) + betaj * ztmp_i);
        }
      } // end loop over jjz
    }
}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

void SNAIntel::compute_deidrj_e(const int jj, const SNA_IVEC &jnum,
                                SNA_DVEC* dedr)
{
  double *ylist_rp = (double *)ylist_r;
  double *ylist_ip = (double *)ylist_i;
  double *dulist_rp = (double *)(dulist_r[0]);
  double *dulist_ip = (double *)(dulist_i[0]);

  for (int k = 0; k < 3; k++)
    dedr[k] = SIMD_set(0.0);

  SNA_IVEC jelem;
  if (chem_flag) jelem = SIMD_load(element + jj);
  else jelem = SIMD256_set(0);

  SIMD_mask m(jj < jnum);

  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j] * vector_width();
    int jju3 = jju * 3;
    SNA_IVEC i = jelem*idxu_max*vector_width() + jju + SIMD256_count();

    for (int mb = 0; 2*mb < j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        SNA_DVEC jjjmambyarray_r = SIMD_gather(m, ylist_rp, i);
        SNA_DVEC jjjmambyarray_i = SIMD_gather(m, ylist_ip, i);
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
          SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
          SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
          dedr[k] = SIMD_add(m, dedr[k], du);
          jju3 += vector_width();
        }
        i = i + vector_width();
      }

    if (j%2 == 0) {
      int mb = j / 2;
      for (int ma = 0; ma < mb; ma++) {
        SNA_DVEC jjjmambyarray_r = SIMD_gather(m, ylist_rp, i);
        SNA_DVEC jjjmambyarray_i = SIMD_gather(m, ylist_ip, i);
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
          SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
          SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
          dedr[k] = SIMD_add(m, dedr[k], du);
          jju3 += vector_width();
        }
        i = i + vector_width();
      }

      SNA_DVEC jjjmambyarray_r = SIMD_gather(m, ylist_rp, i);
      SNA_DVEC jjjmambyarray_i = SIMD_gather(m, ylist_ip, i);
      for (int k = 0; k < 3; k++) {
        SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
        SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
        SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
        dedr[k] = SIMD_fma(m, SIMD_set(0.5), du, dedr[k]);
        jju3 += vector_width();
      }
    } // if j%2
  } // for j

  for (int k = 0; k < 3; k++)
    dedr[k] = dedr[k] * 2.0;
}

/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

void SNAIntel::compute_deidrj(const int jj, const SNA_IVEC &jnum,
                              SNA_DVEC* dedr)
{
  double *ylist_rp = (double *)ylist_r;
  double *ylist_ip = (double *)ylist_i;
  double *dulist_rp = (double *)(dulist_r[0]);
  double *dulist_ip = (double *)(dulist_i[0]);

  for (int k = 0; k < 3; k++)
    dedr[k] = SIMD_set(0.0);

  SIMD_mask m(jj < jnum);

  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j] * vector_width();
    int jju3 = jju * 3;

    for (int mb = 0; 2*mb < j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        SNA_DVEC jjjmambyarray_r = SIMD_load(ylist_rp + jju);
        SNA_DVEC jjjmambyarray_i = SIMD_load(ylist_ip + jju);
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
          SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
          SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
          dedr[k] = SIMD_add(m, dedr[k], du);
          jju3 += vector_width();
        }
        jju += vector_width();
      }

    if (j%2 == 0) {
      int mb = j / 2;
      for (int ma = 0; ma < mb; ma++) {
        SNA_DVEC jjjmambyarray_r = SIMD_load(ylist_rp + jju);
        SNA_DVEC jjjmambyarray_i = SIMD_load(ylist_ip + jju);
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
          SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
          SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
          dedr[k] = SIMD_add(m, dedr[k], du);
          jju3 += vector_width();
        }
        jju += vector_width();
      }

      SNA_DVEC jjjmambyarray_r = SIMD_load(ylist_rp + jju);
      SNA_DVEC jjjmambyarray_i = SIMD_load(ylist_ip + jju);
      for (int k = 0; k < 3; k++) {
        SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
        SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
        SNA_DVEC du = du_r * jjjmambyarray_r + du_i * jjjmambyarray_i;
        dedr[k] = SIMD_fma(m, SIMD_set(0.5), du, dedr[k]);
        jju3 += vector_width();
      }
    } // if j%2
  } // for j

  for (int k = 0; k < 3; k++)
    dedr[k] = dedr[k] * 2.0;
}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
------------------------------------------------------------------------- */

void SNAIntel::compute_bi(const SNA_IVEC &ielem) {
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  double *ulisttot_rp = (double *)ulisttot_r;
  double *ulisttot_ip = (double *)ulisttot_i;
  double *blistp = (double *)blist;

  int itriple = 0;
  int idouble = 0;
  for (int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {

      double *zlist_rp = (double *)(zlist_r + idouble*idxz_max);
      double *zlist_ip = (double *)(zlist_i + idouble*idxz_max);

      for (int elem3 = 0; elem3 < nelements; elem3++) {
        for (int jjb = 0; jjb < idxb_max; jjb++) {
          const int j1 = idxb[jjb].j1;
          const int j2 = idxb[jjb].j2;
          const int j = idxb[jjb].j;

          int jjz = idxz_block[j1][j2][j] * vector_width();
          int jju = (elem3 * idxu_max + idxu_block[j]) * vector_width();
          SNA_DVEC sumzu(0.0);
          for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
              const SNA_DVEC utot_r = SIMD_load(ulisttot_rp + jju);
              const SNA_DVEC utot_i = SIMD_load(ulisttot_ip + jju);
              const SNA_DVEC z_r = SIMD_load(zlist_rp + jjz);
              const SNA_DVEC z_i = SIMD_load(zlist_ip + jjz);
              sumzu = sumzu + utot_r * z_r + utot_i * z_i;
              jjz += vector_width();
              jju += vector_width();
            } // end loop over ma, mb

          // For j even, handle middle column

          if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
              const SNA_DVEC utot_r = SIMD_load(ulisttot_rp + jju);
              const SNA_DVEC utot_i = SIMD_load(ulisttot_ip + jju);
              const SNA_DVEC z_r = SIMD_load(zlist_rp + jjz);
              const SNA_DVEC z_i = SIMD_load(zlist_ip + jjz);
              sumzu = sumzu + utot_r * z_r + utot_i * z_i;
              jjz += vector_width();
              jju += vector_width();
            }

            const SNA_DVEC utot_r = SIMD_load(ulisttot_rp + jju);
            const SNA_DVEC utot_i = SIMD_load(ulisttot_ip + jju);
            const SNA_DVEC z_r = SIMD_load(zlist_rp + jjz);
            const SNA_DVEC z_i = SIMD_load(zlist_ip + jjz);
            sumzu = sumzu + (utot_r * z_r + utot_i * z_i) * 0.5;
          } // end if jeven

          SIMD_store(blistp + (itriple*idxb_max+jjb) * vector_width(),
                     sumzu * 2.0);
        }
        itriple++;
      }
      idouble++;
    }

  // apply bzero shift

  if (bzero_flag) {
    if (!wselfall_flag) {
      SNA_IVEC itriplev = (ielem*nelements+ielem)*nelements+ielem;
      for (int jjb = 0; jjb < idxb_max; jjb++) {
        const int j = idxb[jjb].j;
        SNA_IVEC i = (itriplev*idxb_max+jjb) * vector_width() + SIMD256_count();
        SIMD_scatter(blistp, i, SIMD_gather(blistp, i) - bzero[j]);
      } // end loop over JJ
    } else {
      int itriple = 0;
      for (int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {
          for (int elem3 = 0; elem3 < nelements; elem3++) {
            for (int jjb = 0; jjb < idxb_max; jjb++) {
              const int j = idxb[jjb].j;
              int i = (itriple*idxb_max+jjb) * vector_width();
              SIMD_store(blistp + i, SIMD_load(blistp + i) - bzero[j]);
            } // end loop over JJ
            itriple++;
          } // end loop over elem3
        } // end loop over elem1,elem2
    }
  }
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

void SNAIntel::compute_duidrj(const int jj, const SNA_IVEC &jnum)
{
  const SNA_DVEC x = rij[jj][0];
  const SNA_DVEC y = rij[jj][1];
  const SNA_DVEC z = rij[jj][2];
  const SNA_DVEC rcut = rcutij[jj];
  const SNA_DVEC rsq = x * x + y * y + z * z;
  const SNA_DVEC r = SIMD_sqrt(rsq);
  const SNA_DVEC rscale0 = SIMD_rcp(rcut - rmin0) * rfac0 * MY_PI;
  const SNA_DVEC theta0 = (r - rmin0) * rscale0;
  const SNA_DVEC z0 = r * SIMD_rcp(SIMD_tan(theta0));
  const SNA_DVEC dz0dr = z0 * SIMD_rcp(r) - (r*rscale0) * (rsq + z0 * z0) *
    SIMD_rcp(rsq);
  compute_duarray(x, y, z, z0, r, dz0dr, wj[jj], rcut, jj, jnum);
}

/* ---------------------------------------------------------------------- */

void SNAIntel::zero_uarraytot(const SNA_IVEC &ielem)
{
  double *ulisttot_rp = (double *)ulisttot_r;
  double *ulisttot_ip = (double *)ulisttot_i;
  for (int jelem = 0; jelem < nelements; jelem++)
    for (int j = 0; j <= twojmax; j++) {
      int jju = (jelem * idxu_max + idxu_block[j]) * vector_width();
      for (int mb = 0; mb <= j; mb++) {
        for (int ma = 0; ma <= j; ma++) {
          SIMD_store(ulisttot_rp + jju, SIMD_set(0.0));
          SIMD_store(ulisttot_ip + jju, SIMD_set(0.0));

          // utot(j,ma,ma) = wself, sometimes
          if (ma == mb) {
            if (wselfall_flag || nelements == 1)
              SIMD_store(ulisttot_rp + jju, SIMD_set(wself));
            else {
              SIMD_mask m(ielem == jelem);
              SIMD_store(ulisttot_rp + jju,
                         SIMD_zero_masked(~m, SIMD_set(wself)));
            }
          }
          jju += vector_width();
        }
      }
    }
}



/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

void SNAIntel::add_uarraytot(const SNA_DVEC &r, const int jj,
                             const SNA_IVEC &jnum)
{
  SNA_DVEC sfac = compute_sfac(r, rcutij[jj], sinnerij[jj], dinnerij[jj]);
  sfac *= wj[jj];

  double *ulisttot_rp = (double *)ulisttot_r;
  double *ulisttot_ip = (double *)ulisttot_i;
  const double* ulist_r = (double *)(ulist_r_ij[jj]);
  const double* ulist_i = (double *)(ulist_i_ij[jj]);

  SIMD_mask m(jj < jnum);

  if (chem_flag && nelements > 1) {
    SNA_IVEC jelem = SIMD_load(element+jj);
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j] * vector_width();
      SNA_IVEC i = jelem*idxu_max*vector_width() + jju + SIMD256_count();
      for (int mb = 0; mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          SNA_DVEC utot_r = SIMD_gather(m, ulisttot_rp, i);
          SNA_DVEC utot_i = SIMD_gather(m, ulisttot_ip, i);
          utot_r = SIMD_fma(m, sfac, SIMD_load(ulist_r + jju), utot_r);
          utot_i = SIMD_fma(m, sfac, SIMD_load(ulist_i + jju), utot_i);
          SIMD_scatter(m, ulisttot_rp, i, utot_r);
          SIMD_scatter(m, ulisttot_ip, i, utot_i);
          jju += vector_width();
          i = i + vector_width();
        }
    }
  } else {
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j] * vector_width();
      for (int mb = 0; mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          SNA_DVEC utot_r = SIMD_load(ulisttot_rp + jju);
          SNA_DVEC utot_i = SIMD_load(ulisttot_ip + jju);
          utot_r = SIMD_fma(m, sfac, SIMD_load(ulist_r + jju), utot_r);
          utot_i = SIMD_fma(m, sfac, SIMD_load(ulist_i + jju), utot_i);
          SIMD_store(ulisttot_rp + jju, utot_r);
          SIMD_store(ulisttot_ip + jju, utot_i);
          jju += vector_width();
        }
    }
  }
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

void SNAIntel::compute_uarray(const SNA_DVEC &x, const SNA_DVEC &y,
                              const SNA_DVEC &z, const SNA_DVEC &z0,
                              const SNA_DVEC &r, const int jj,
                              const SNA_IVEC &jnum)
{
  // compute Cayley-Klein parameters for unit quaternion

  const SNA_DVEC r0inv = SIMD_invsqrt(r * r + z0 * z0);
  const SNA_DVEC a_r = z0 * r0inv;
  const SNA_DVEC a_i = -z * r0inv;
  const SNA_DVEC b_r = y * r0inv;
  const SNA_DVEC b_i = -x * r0inv;

  // VMK Section 4.8.2

  double *ulist_rp = (double *)(ulist_r_ij[jj]);
  double *ulist_ip = (double *)(ulist_i_ij[jj]);

  SIMD_store(ulist_rp, SIMD_set(1.0));
  SIMD_store(ulist_ip, SIMD_set(0.0));

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j] * vector_width();
    int jjup = idxu_block[j-1] * vector_width();

    // fill in left side of matrix layer from previous layer

    for (int mb = 0; 2*mb <= j; mb++) {
      SIMD_store(ulist_rp + jju, SIMD_set(0.0));
      SIMD_store(ulist_ip + jju, SIMD_set(0.0));

      for (int ma = 0; ma < j; ma++) {
        double rootpq = rootpqarray[j - ma][j - mb];
        SNA_DVEC u_r = SIMD_load(ulist_rp + jju);
        SNA_DVEC u_i = SIMD_load(ulist_ip + jju);
        const SNA_DVEC up_r = SIMD_load(ulist_rp + jjup);
        const SNA_DVEC up_i = SIMD_load(ulist_ip + jjup);

        SNA_DVEC u_ro, u_io;

        u_ro = a_r * up_r + a_i * up_i;
        u_r = SIMD_fma(SIMD_set(rootpq), u_ro, u_r);
        SIMD_store(ulist_rp + jju, u_r);
        u_io = a_r * up_i - a_i * up_r;
        u_i = SIMD_fma(SIMD_set(rootpq), u_io, u_i);
        SIMD_store(ulist_ip + jju, u_i);

        jju += vector_width();

        rootpq = -rootpqarray[ma + 1][j - mb];
        u_r = (b_r * up_r + b_i * up_i) * rootpq;
        SIMD_store(ulist_rp + jju, u_r);
        u_i = (b_r * up_i - b_i * up_r) * rootpq;
        SIMD_store(ulist_ip + jju, u_i);

        jjup += vector_width();
      }
      jju += vector_width();
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = (jju+(j+1)*(j+1)-1) * vector_width();
    jju *=  vector_width();
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        if (mapar == 1) {
          SIMD_store(ulist_rp + jjup, SIMD_load(ulist_rp + jju));
          SIMD_store(ulist_ip + jjup, -SIMD_load(ulist_ip + jju));
        } else {
          SIMD_store(ulist_rp + jjup, -SIMD_load(ulist_rp + jju));
          SIMD_store(ulist_ip + jjup, SIMD_load(ulist_ip + jju));
        }
        mapar = -mapar;
        jju += vector_width();
        jjup -= vector_width();
      }
      mbpar = -mbpar;
    }
  }
}

/* ----------------------------------------------------------------------
   Compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

void SNAIntel::compute_duarray(const SNA_DVEC &x, const SNA_DVEC &y,
                               const SNA_DVEC &z, const SNA_DVEC &z0,
                               const SNA_DVEC &r, const SNA_DVEC &dz0dr,
                               const SNA_DVEC &wj, const SNA_DVEC &rcut,
                               const int jj, const SNA_IVEC &jnum)
{
  const SNA_DVEC rinv = SIMD_rcp(r);
  const SNA_DVEC r0inv = SIMD_invsqrt(r * r + z0 * z0);
  SNA_DVEC up[3];
  up[0] = x * rinv;
  up[1] = y * rinv;
  up[2] = z * rinv;
  const SNA_DVEC a_r = z0 * r0inv;
  const SNA_DVEC a_i = -z * r0inv;
  const SNA_DVEC b_r = y * r0inv;
  const SNA_DVEC b_i = -x * r0inv;
  const SNA_DVEC dr0invdr = -SIMD_pow(r0inv, 3.0) * (r + z0 * dz0dr);

  SNA_DVEC dr0inv[3], da_r[3], da_i[3];
  for (int k = 0; k < 3; k++) {
    dr0inv[k] = dr0invdr * up[k];
    da_r[k] = dz0dr * up[k] * r0inv + z0 * dr0inv[k];
    da_i[k] = -z * dr0inv[k];
  }
  da_i[2] += -r0inv;

  double *ulist_rp = (double *)(ulist_r_ij[jj]);
  double *ulist_ip = (double *)(ulist_i_ij[jj]);
  double *dulist_rp = (double *)(dulist_r[0]);
  double *dulist_ip = (double *)(dulist_i[0]);

  SNA_DVEC db_r[3], db_i[3];
  for (int k = 0; k < 3; k++) {
    SIMD_store(dulist_rp + k * vector_width(), SIMD_set(0.0));
    SIMD_store(dulist_ip + k * vector_width(), SIMD_set(0.0));
    db_r[k] = y * dr0inv[k];
    db_i[k] = -x * dr0inv[k];
  }
  db_i[0] -= r0inv;
  db_r[1] += r0inv;

  for (int j = 1; j <= twojmax; j++) {
    int jju3 = idxu_block[j] * 3 * vector_width();
    int jjup = idxu_block[j-1] * vector_width();
    int jjup3 = jjup * 3;
    for (int mb = 0; 2*mb <= j; mb++) {
      for (int k = 0; k < 3; k++) {
        SIMD_store(dulist_rp + jju3 + k * vector_width(), SIMD_set(0.0));
        SIMD_store(dulist_ip + jju3 + k * vector_width(), SIMD_set(0.0));
      }

      for (int ma = 0; ma < j; ma++) {
        const double rootpq = rootpqarray[j - ma][j - mb];
        const double mrootpq = -rootpqarray[ma + 1][j - mb];
        const SNA_DVEC up_r = SIMD_load(ulist_rp + jjup);
        const SNA_DVEC up_i = SIMD_load(ulist_ip + jjup);
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = SIMD_load(dulist_rp + jju3);
          SNA_DVEC du_i = SIMD_load(dulist_ip + jju3);
          const SNA_DVEC dup_r = SIMD_load(dulist_rp + jjup3);
          const SNA_DVEC dup_i = SIMD_load(dulist_ip + jjup3);

          SNA_DVEC du_ro, du_io;

          du_ro = (da_r[k]*up_r + da_i[k]*up_i + a_r*dup_r + a_i*dup_i);
          du_r = SIMD_fma(SIMD_set(rootpq), du_ro, du_r);
          SIMD_store(dulist_rp + jju3, du_r);

          du_io = (da_r[k]*up_i - da_i[k]*up_r + a_r*dup_i - a_i*dup_r);
          du_i = SIMD_fma(SIMD_set(rootpq), du_io, du_i);
          SIMD_store(dulist_ip + jju3, du_i);

          du_r = (db_r[k]*up_r + db_i[k]*up_i + b_r*dup_r + b_i*dup_i);
          SIMD_store(dulist_rp + jju3 + 3 * vector_width(), du_r * mrootpq);

          du_i = (db_r[k]*up_i - db_i[k]*up_r + b_r*dup_i - b_i*dup_r);
          SIMD_store(dulist_ip + jju3 + 3 * vector_width(), du_i * mrootpq);

          jju3 += vector_width();
          jjup3 += vector_width();
        }
        jjup += vector_width();
      } // for ma
      jju3 += 3 * vector_width();
    } // for mb

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    SNA_DVEC *du_r_p = dulist_r[0];
    SNA_DVEC *du_i_p = dulist_i[0];

    int jju = idxu_block[j];
    jjup = (jju+(j+1)*(j+1)-1) * 3 * vector_width();
    jju *=  3 * vector_width();
    int mbpar = 1;
    for (int mb = 0; 2*mb <= j; mb++) {
      int mapar = mbpar;
      for (int ma = 0; ma <= j; ma++) {
        if (mapar == 1) {
          for (int k = 0; k < 3; k++) {
            SIMD_store(dulist_rp + jjup, SIMD_load(dulist_rp + jju));
            SIMD_store(dulist_ip + jjup, -SIMD_load(dulist_ip + jju));
            jju += vector_width();
            jjup += vector_width();
          }
        } else {
          for (int k = 0; k < 3; k++) {
            SIMD_store(dulist_rp + jjup, -SIMD_load(dulist_rp + jju));
            SIMD_store(dulist_ip + jjup, SIMD_load(dulist_ip + jju));
            jju += vector_width();
            jjup += vector_width();
          }
        }
        mapar = -mapar;
        jjup -= 6 * vector_width();
      } // for ma
      mbpar = -mbpar;
    } // for mb
  } // for j

  SNA_DVEC dsfac;
  SNA_DVEC sfac = compute_sfac_dsfac(r, rcut, sinnerij[jj], dinnerij[jj],
                                      dsfac);
  sfac = sfac * wj;
  dsfac = dsfac * wj;

  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j] * vector_width();
    int jju3 = jju * 3;
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        const SNA_DVEC ur_dsfac = dsfac * SIMD_load(ulist_rp + jju);
        const SNA_DVEC ui_dsfac = dsfac * SIMD_load(ulist_ip + jju);
        jju += vector_width();
        for (int k = 0; k < 3; k++) {
          SNA_DVEC du_r = ur_dsfac * up[k] + sfac * SIMD_load(dulist_rp+jju3);
          SIMD_store(dulist_rp + jju3, du_r);
          SNA_DVEC du_i = ui_dsfac * up[k] + sfac * SIMD_load(dulist_ip+jju3);
          SIMD_store(dulist_ip + jju3, du_i);
          jju3 += vector_width();
        }
      }
  }
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

double SNAIntel::memory_usage()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  double bytes;

  bytes = 0;

  bytes += (double)jdimpq*jdimpq * sizeof(double);               // pqarray
  bytes += (double)idxcg_max * sizeof(double);                   // cglist

  bytes += (double)nmax * idxu_max * sizeof(SNA_DVEC) * 2;       // ulist_ij
  bytes += (double)idxu_max * nelements * sizeof(SNA_DVEC) * 2;  // ulisttot
  bytes += (double)idxu_max * 3 * sizeof(SNA_DVEC) * 2;          // dulist

  bytes += (double)idxz_max * ndoubles * sizeof(SNA_DVEC) * 2;   // zlist
  bytes += (double)idxb_max * ntriples * sizeof(SNA_DVEC);       // blist
  bytes += (double)idxb_max * ntriples * 3 * sizeof(double);     // dblist
  bytes += (double)idxu_max * nelements * sizeof(SNA_DVEC) * 2;  // ylist

  bytes += (double)jdim * jdim * jdim * sizeof(int);             // idxcg_block
  bytes += (double)jdim * sizeof(int);                           // idxu_block
  bytes += (double)jdim * jdim * jdim * sizeof(int);             // idxz_block
  bytes += (double)jdim * jdim * jdim * sizeof(int);             // idxb_block

  bytes += (double)idxz_max * sizeof(SNA_ZINDICES);              // idxz
  bytes += (double)idxb_max * sizeof(SNA_BINDICES);              // idxb

  if (bzero_flag)
  bytes += (double)jdim * sizeof(double);                        // bzero

  bytes += (double)nmax * 3 * sizeof(SNA_DVEC);                  // rij
  bytes += (double)nmax * sizeof(SNA_IVEC);                      // inside
  bytes += (double)nmax * sizeof(SNA_DVEC);                      // wj
  bytes += (double)nmax * sizeof(SNA_DVEC);                      // rcutij
  bytes += (double)nmax * sizeof(SNA_DVEC);                      // sinnerij
  bytes += (double)nmax * sizeof(SNA_DVEC);                      // dinnerij
  if (chem_flag) bytes += (double)nmax * sizeof(SNA_IVEC);       // element

  return bytes;
}

/* ---------------------------------------------------------------------- */

void SNAIntel::create_twojmax_arrays()
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
    bzero = nullptr;

}

/* ---------------------------------------------------------------------- */

void SNAIntel::destroy_twojmax_arrays()
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
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

double SNAIntel::deltacg(int j1, int j2, int j)
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

void SNAIntel::init_clebsch_gordan()
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j) {
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

void SNAIntel::print_clebsch_gordan()
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

void SNAIntel::init_rootpqarray()
{
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p][q] = sqrt(static_cast<double>(p)/q);
}

/* ---------------------------------------------------------------------- */

void SNAIntel::compute_ncoeff()
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
  if (chem_flag)
    ncoeff = ncount*ntriples;
  else
    ncoeff = ncount;
}

/* ---------------------------------------------------------------------- */

double SNAIntel::compute_sfac(double r, double rcut, double sinner, double dinner)
{
  double sfac;

  // calculate sfac = sfac_outer

  if (switch_flag == 0) sfac = 1.0;
  else if (r <= rmin0) sfac = 1.0;
  else if (r > rcut) sfac = 0.0;
  else {
    double rcutfac = MY_PI / (rcut - rmin0);
    sfac = 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
  }

  // calculate sfac *= sfac_inner, rarely visited

  if (switch_inner_flag == 1 && r < sinner + dinner) {
    if (r > sinner - dinner) {
      double rcutfac = MY_PI2 / dinner;
      sfac *= 0.5 * (1.0 - cos(MY_PI2 + (r - sinner) * rcutfac));
    } else sfac = 0.0;
  }

  return sfac;
}

/* ---------------------------------------------------------------------- */

SNA_DVEC SNAIntel::compute_sfac(const SNA_DVEC &r, const SNA_DVEC &rcut,
                                const SNA_DVEC &sinner, const SNA_DVEC &dinner)
{
  // calculate sfac = sfac_outer

  // if (switch_flag == 0 || r <= rmin0)
  SNA_DVEC sfac = SIMD_set(1.0);
  if (switch_flag != 0) {
    // r <= rcut && r > rmin0
    const SIMD_mask i(r > rmin0);
    const SIMD_mask m(r <= rcut);
    const SNA_DVEC rcutfac = SIMD_rcp(rcut - rmin0) * MY_PI;
    const SNA_DVEC sfac_m = (SIMD_cos((r - rmin0) * rcutfac) + 1.0) * 0.5;
    sfac = SIMD_set(sfac, m & i, sfac_m);
    // (r > rcut) && (r> rmin0)
    sfac = SIMD_zero_masked(m | i, sfac);
  }

  // calculate sfac *= sfac_inner, rarely visited

  if (switch_inner_flag == 1) {
    const SIMD_mask m(r < sinner + dinner);
    // if any(m)
    const SIMD_mask i(r > sinner - dinner);
    const SNA_DVEC rcutfac = SIMD_rcp(dinner) * MY_PI2;
    const SNA_DVEC sfac_m = (SIMD_set(1.0) - SIMD_cos((r-sinner) * rcutfac +
                                                      MY_PI2)) * 0.5;
    sfac = SIMD_set(sfac, m & i, sfac_m);
    sfac = SIMD_zero_masked((~m) | i, sfac);
  }

  return sfac;
}

/* ---------------------------------------------------------------------- */

SNA_DVEC SNAIntel::compute_sfac_dsfac(const SNA_DVEC & r,
                                      const SNA_DVEC & rcut,
                                      const SNA_DVEC & sinner,
                                      const SNA_DVEC & dinner,
                                      SNA_DVEC &dsfac)
{
  // calculate sfac = sfac_outer

  // if (switch_flag == 0 || r <= rmin0)
  SNA_DVEC sfac = SIMD_set(1.0);
  dsfac = SIMD_set(0.0);
  if (switch_flag != 0) {
    // r <= rcut && r > rmin0
    const SIMD_mask i(r > rmin0);
    const SIMD_mask m(r <= rcut);
    const SNA_DVEC rcutfac = SIMD_rcp(rcut - rmin0) * MY_PI;
    const SNA_DVEC trig_arg = (r - rmin0) * rcutfac;
    const SNA_DVEC sfac_m = (SIMD_cos(trig_arg) + 1.0) * 0.5;
    const SNA_DVEC dsfac_m = SIMD_sin(trig_arg) * rcutfac * -0.5;
    sfac = SIMD_set(sfac, m & i, sfac_m);
    dsfac = SIMD_set(dsfac, m & i, dsfac_m);
    // (r > rcut) && (r> rmin0)
    sfac = SIMD_zero_masked(m | i, sfac);
  }

  // calculate sfac *= sfac_inner, rarely visited

  if (switch_inner_flag == 1) {
    const SIMD_mask m(r < sinner + dinner);
    const SIMD_mask i(r > sinner - dinner);
    if (any(m & i)) {
      const SNA_DVEC rcutfac = SIMD_rcp(dinner) * MY_PI2;
      const SNA_DVEC trig_arg = (r - sinner) * rcutfac + MY_PI2;
      const SNA_DVEC sfac_inner = (SIMD_set(1.0) - SIMD_cos(trig_arg)) * 0.5;
      const SNA_DVEC dsfac_inner = rcutfac * 0.5 * SIMD_sin(trig_arg);
      dsfac = SIMD_set(dsfac, m & i, dsfac * sfac_inner +
                       sfac * dsfac_inner);
      sfac = SIMD_set(sfac, m & i, sfac_inner);
    }
    sfac = SIMD_zero_masked((~m) | i, sfac);
    dsfac = SIMD_zero_masked((~m) | i, dsfac);
  }

  return sfac;
}

template void SNAIntel::compute_zi_or_yi<1>(const SNA_DVEC *);
template void SNAIntel::compute_zi_or_yi<0>(const SNA_DVEC *);

#endif
#endif
