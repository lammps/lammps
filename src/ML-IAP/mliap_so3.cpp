/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Contributing authors: Byungkyun Kang (University of Nevada, Las Vegas)
 ------------------------------------------------------------------------- */

#include "mliap_so3.h"

#include "error.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"

#include <cmath>

#include "mliap_so3_math.h"

using namespace SO3Math;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 1.0e-8

/* ---------------------------------------------------------------------- */

MLIAP_SO3::MLIAP_SO3(LAMMPS *lmp, double vrcut, int vlmax, int vnmax, double valpha) : Pointers(lmp)
{
  m_rcut = vrcut;
  m_alpha = valpha;
  m_lmax = vlmax;
  m_nmax = vnmax;
  compute_ncoeff();

  m_Nmax = (m_nmax + m_lmax + 1) * 10;
  m_numYlms = (m_lmax + 1) * (m_lmax + 1);

  m_ellpl1 = nullptr;
  m_ellm1 = nullptr;
  m_pfac = nullptr;
  m_dfac0 = nullptr;
  m_dfac1 = nullptr;
  m_dfac2 = nullptr;
  m_dfac3 = nullptr;
  m_dfac4 = nullptr;
  m_dfac5 = nullptr;

  m_g_array = nullptr;
  m_w = nullptr;
  m_rootpq = nullptr;
  m_idxu_block = nullptr;
  m_idxylm = nullptr;

  m_rip_array = nullptr;
  m_rip_darray = nullptr;

  m_sbes_array = nullptr;
  m_sbes_darray = nullptr;

  m_plist_r = nullptr;
  m_plist_i = nullptr;

  m_clist_r = nullptr;
  m_clist_i = nullptr;

  m_ulist_r = nullptr;
  m_ulist_i = nullptr;

  m_Ylms_r = nullptr;
  m_Ylms_i = nullptr;

  m_dYlm_r = nullptr;
  m_dYlm_i = nullptr;

  m_dplist_r = nullptr;
  m_dplist_i = nullptr;

  m_dclist_r = nullptr;
  m_dclist_i = nullptr;

  m_tempdp_r = nullptr;

  m_clisttot_r = nullptr;
  m_clisttot_i = nullptr;

  m_init_arrays = 0;
  m_dfac_l1 = m_dfac_l2 = 0;
  m_pfac_l1 = m_pfac_l2 = 0;
  m_idxu_count = m_idxy_count = 0;
  alloc_init = alloc_arrays = 0.0;
}

/* ---------------------------------------------------------------------- */

MLIAP_SO3::~MLIAP_SO3()
{

  memory->destroy(m_ellpl1);
  memory->destroy(m_ellm1);
  memory->destroy(m_pfac);
  memory->destroy(m_dfac0);
  memory->destroy(m_pfac);
  memory->destroy(m_dfac0);
  memory->destroy(m_dfac1);
  memory->destroy(m_dfac2);
  memory->destroy(m_dfac3);
  memory->destroy(m_dfac4);
  memory->destroy(m_dfac5);
  memory->destroy(m_w);
  memory->destroy(m_g_array);

  memory->destroy(m_rootpq);
  memory->destroy(m_idxu_block);
  memory->destroy(m_idxylm);

  memory->destroy(m_rip_array);
  memory->destroy(m_rip_darray);

  memory->destroy(m_sbes_array);
  memory->destroy(m_sbes_darray);

  memory->destroy(m_plist_r);
  memory->destroy(m_plist_i);

  memory->destroy(m_clist_r);
  memory->destroy(m_clist_i);

  memory->destroy(m_ulist_r);
  memory->destroy(m_ulist_i);

  memory->destroy(m_Ylms_r);
  memory->destroy(m_Ylms_i);
  memory->destroy(m_dYlm_r);
  memory->destroy(m_dYlm_i);

  memory->destroy(m_dplist_r);
  memory->destroy(m_dplist_i);

  memory->destroy(m_dclist_r);
  memory->destroy(m_dclist_i);

  memory->destroy(m_tempdp_r);

  memory->destroy(m_clisttot_r);
  memory->destroy(m_clisttot_i);
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::compute_ncoeff()
{

  ncoeff = m_nmax * (m_nmax + 1) * (m_lmax + 1) / 2;
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::init()
{
  int i, totali;

  totali = m_lmax + 1;
  memory->destroy(m_ellpl1);
  memory->create(m_ellpl1, totali, "MLIAP_SO3:m_ellpl1");
  memory->destroy(m_ellm1);
  memory->create(m_ellm1, totali, "MLIAP_SO3:m_ellm1");
  alloc_init = 2.0 * totali * sizeof(double);

  for (int l = 1; l < m_lmax + 1; l++) {
    m_ellpl1[l] = get_sum(0, l + 2, 1, 2);
    m_ellm1[l] = get_sum(0, l, 1, 2);
  }

  double pfac1 = 1.0 / 4.0 / MY_PI;
  m_pfac_l1 = m_lmax + 2;
  m_pfac_l2 = (m_lmax + 2) * (m_lmax + 2) + 1;
  totali = m_pfac_l1 * m_pfac_l2;
  memory->destroy(m_pfac);
  memory->create(m_pfac, totali, "MLIAP_SO3:m_pfac");
  alloc_init += totali * sizeof(double);
  for (int l = 0; l < m_lmax + 2; l++)
    for (int m = -l; m < l + 1; m++)
      m_pfac[l * m_pfac_l2 + m] = sqrt((2.0 * l + 1.0) * pfac1) * powsign(m);

  m_dfac_l1 = m_lmax + 1;
  m_dfac_l2 = m_numYlms + 1;
  totali = m_dfac_l1 * m_dfac_l2;
  memory->destroy(m_dfac0);
  memory->create(m_dfac0, totali, "MLIAP_SO3:m_dfac0");
  memory->destroy(m_dfac1);
  memory->create(m_dfac1, totali, "MLIAP_SO3:m_dfac1");
  memory->destroy(m_dfac2);
  memory->create(m_dfac2, totali, "MLIAP_SO3:m_dfac2");
  memory->destroy(m_dfac3);
  memory->create(m_dfac3, totali, "MLIAP_SO3:m_dfac3");
  memory->destroy(m_dfac4);
  memory->create(m_dfac4, totali, "MLIAP_SO3:m_dfac4");
  memory->destroy(m_dfac5);
  memory->create(m_dfac5, totali, "MLIAP_SO3:m_dfac5");
  alloc_init += 6.0 * totali * sizeof(double);

  for (int l = 1; l < m_lmax + 1; l++)
    for (int m = -l; m < l + 1; m++) {
      m_dfac0[l * m_dfac_l2 + m] =
          -sqrt(((l + 1.0) * (l + 1.0) - m * m) / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      m_dfac1[l * m_dfac_l2 + m] =
          sqrt((l * l - m * m) / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
      m_dfac2[l * m_dfac_l2 + m] =
          -sqrt((l + m + 1.0) * (l + m + 2.0) / 2.0 / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      m_dfac3[l * m_dfac_l2 + m] =
          sqrt((l - m - 1.0) * (l - m) / 2.0 / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
      m_dfac4[l * m_dfac_l2 + m] =
          -sqrt((l - m + 1.0) * (l - m + 2.0) / 2.0 / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      m_dfac5[l * m_dfac_l2 + m] =
          sqrt((l + m - 1.0) * (l + m) / 2.0 / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
    }

  totali = m_nmax * m_nmax;
  memory->destroy(m_w);
  memory->create(m_w, totali, "MLIAP_SO3:w");
  alloc_init += totali * sizeof(double);

  for (i = 0; i < totali; i++) m_w[i] = 0.0;

  compute_W(m_nmax, m_w);

  totali = m_nmax * m_Nmax;
  memory->destroy(m_g_array);
  memory->create(m_g_array, totali, "MLIAP_SO3:g_array");
  alloc_init += totali * sizeof(double);

  for (i = 0; i < totali; i++) m_g_array[i] = 0.0;

  init_garray(m_nmax, m_lmax, m_rcut, m_alpha, m_w, m_nmax, m_g_array, m_Nmax);

  int twolmax;
  twolmax = 2 * (m_lmax + 1);
  int m_ldim = twolmax + 1;
  totali = m_ldim * m_ldim;
  memory->destroy(m_rootpq);
  memory->create(m_rootpq, totali, "MLIAP_SO3:rootpq");
  alloc_init += totali * sizeof(double);

  for (i = 0; i < totali; i++) m_rootpq[i] = 0.0;

  for (int p = 1; p < m_ldim; p++)
    for (int q = 1; q < m_ldim; q++) m_rootpq[p * m_ldim + q] = sqrt(static_cast<double>(p) / q);

  memory->destroy(m_idxu_block);
  memory->create(m_idxu_block, m_ldim, "MLIAP_SO3:idxu_bloc");
  alloc_init += totali * sizeof(double);

  for (i = 0; i < m_ldim; i++) m_idxu_block[i] = 0;

  totali = square(m_lmax + 2);
  memory->destroy(m_idxylm);
  memory->create(m_idxylm, totali, "MLIAP_SO3:idxylm");
  alloc_init += totali * sizeof(double);

  for (i = 0; i < totali; i++) m_idxylm[i] = 0;

  m_idxu_count = m_idxy_count = 0;

  for (int l = 0; l < m_ldim; l++) {
    m_idxu_block[l] = m_idxu_count;
    for (int mb = 0; mb < l + 1; mb++)
      for (int ma = 0; ma < l + 1; ma++) {
        if (l % 2 == 0 && ma == l / 2) {
          m_idxylm[m_idxy_count] = m_idxu_count;
          m_idxy_count += 1;
        }
        m_idxu_count += 1;
      }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::init_arrays(int nlocal, int ncoefs)
{
  int totali = nlocal * ncoefs;
  memory->destroy(m_plist_r);
  memory->create(m_plist_r, totali, "MLIAP_SO3:m_plist_r");
  memory->destroy(m_plist_i);
  memory->create(m_plist_i, totali, "MLIAP_SO3:m_plist_i");
  alloc_arrays = 2.0 * totali * sizeof(double);

  totali = m_nmax * m_numYlms;
  memory->destroy(m_clist_r);
  memory->create(m_clist_r, totali, "MLIAP_SO3:m_clist_r");
  memory->destroy(m_clist_i);
  memory->create(m_clist_i, totali, "MLIAP_SO3:m_clist_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = m_idxu_count;
  memory->destroy(m_ulist_r);
  memory->create(m_ulist_r, totali, "MLIAP_SO3:m_ulist_r");
  memory->destroy(m_ulist_i);
  memory->create(m_ulist_i, totali, "MLIAP_SO3:m_ulist_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = (m_lmax + 2) * (m_lmax + 2);
  memory->destroy(m_Ylms_r);
  memory->create(m_Ylms_r, totali, "MLIAP_SO3:m_Ylms_r");
  memory->destroy(m_Ylms_i);
  memory->create(m_Ylms_i, totali, "MLIAP_SO3:m_Ylms_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = m_numYlms * 3;
  memory->destroy(m_dYlm_r);
  memory->create(m_dYlm_r, totali, "MLIAP_SO3:m_dYlm_r");
  memory->destroy(m_dYlm_i);
  memory->create(m_dYlm_i, totali, "MLIAP_SO3:m_dYlm_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = m_nmax * m_numYlms * 3;
  memory->destroy(m_dclist_r);
  memory->create(m_dclist_r, totali, "MLIAP_SO3:m_dclist_r");
  memory->destroy(m_dclist_i);
  memory->create(m_dclist_i, totali, "MLIAP_SO3:m_dclist_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = 3 * m_nmax * (m_nmax + 1) * (m_lmax + 1) / 2;
  memory->destroy(m_tempdp_r);
  memory->create(m_tempdp_r, totali, "MLIAP_SO3:m_tempdp_r");
  alloc_arrays += totali * sizeof(double);

  totali = m_nmax * m_numYlms;
  memory->destroy(m_clisttot_r);
  memory->create(m_clisttot_r, totali, "MLIAP_SO3:m_clisttot_r");
  memory->destroy(m_clisttot_i);
  memory->create(m_clisttot_i, totali, "MLIAP_SO3:m_clisttot_i");
  alloc_arrays += 2.0 * totali * sizeof(double);
  m_init_arrays = 1;
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::memory_usage()
{
  return alloc_init + alloc_arrays;
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::compute_W(int nmax, double *arr)
{
  int alpha, beta, temp1, temp2;

  for (alpha = 1; alpha < nmax + 1; alpha++) {
    temp1 = (2 * alpha + 5) * (2 * alpha + 6) * (2 * alpha + 7);
    for (beta = 1; beta < alpha + 1; beta++) {
      temp2 = (2 * beta + 5) * (2 * beta + 6) * (2 * beta + 7);
      arr[(alpha - 1) * nmax + beta - 1] =
          sqrt(temp1 * temp2) / (5 + alpha + beta) / (6 + alpha + beta) / (7 + alpha + beta);
      arr[(beta - 1) * nmax + alpha - 1] = arr[(alpha - 1) * nmax + beta - 1];
    }
  }

  int i, j, k, n = nmax;
  auto outeig = new double[n];
  auto outeigvec = new double[n * n];
  auto arrinv = new double[n * n];

  auto sqrtD = new double[n * n];
  auto tempM = new double[n * n];

  auto temparr = new double *[n];
  auto tempvl = new double *[n];
  auto tempout = new double[n];

  int info;

  info = invert_matrix(n, arr, arrinv);
  if (info != 0) error->all(FLERR, "Invert matrix Error in W calculation!");

  for (i = 0; i < n; i++) {
    temparr[i] = new double[n];
    tempvl[i] = new double[n];
    for (j = 0; j < n; j++) temparr[i][j] = arrinv[i * n + j];
  }

  jacobin(n, temparr, tempout, tempvl);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) outeigvec[i * n + j] = tempvl[i][j];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (i == j)
        sqrtD[i * n + j] = sqrt(tempout[i]);
      else
        sqrtD[i * n + j] = 0.0;

  double dtemp;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      dtemp = 0;
      for (k = 0; k < n; k++) dtemp += outeigvec[i * n + k] * sqrtD[k * n + j];

      tempM[i * n + j] = dtemp;
    }

  info = invert_matrix(n, outeigvec, arrinv);
  if (info != 0) error->all(FLERR, "Invert matrix Error in W calculation!");

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      dtemp = 0;
      for (k = 0; k < n; k++) dtemp += tempM[i * n + k] * arrinv[k * n + j];

      arr[i * n + j] = dtemp;
    }

  delete[] outeig;
  delete[] outeigvec;
  delete[] arrinv;

  delete[] sqrtD;
  delete[] tempM;

  for (i = 0; i < n; i++) {
    delete[] temparr[i];
    delete[] tempvl[i];
  }
  delete[] temparr;
  delete[] tempvl;
  delete[] tempout;
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::compute_pi(int nmax, int lmax, double *clisttot_r, double *clisttot_i, int lcl2,
                           double *plist_r, double *plist_i, int lpl2, int indpl)
{
  int n1, n2, j, l, m, i = 0;
  double norm;
  for (n1 = 0; n1 < nmax; n1++)
    for (n2 = 0; n2 < n1 + 1; n2++) {
      j = 0;
      for (l = 0; l < lmax + 1; l++) {
        norm = 2.0 * sqrt(2.0) * MY_PI / sqrt(2.0 * l + 1.0);

        for (m = -l; m < l + 1; m++) {

          plist_r[lpl2 * indpl + i] += (clisttot_r[lcl2 * n1 + j] * clisttot_r[lcl2 * n2 + j] +
                                        clisttot_i[lcl2 * n1 + j] * clisttot_i[lcl2 * n2 + j]) *
              norm;
          plist_i[lpl2 * indpl + i] += (-clisttot_r[lcl2 * n1 + j] * clisttot_i[lcl2 * n2 + j] +
                                        clisttot_i[lcl2 * n1 + j] * clisttot_r[lcl2 * n2 + j]) *
              norm;

          j += 1;
        }
        i += 1;
      }
    }
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::phi(double r, int alpha, double rcut)
{
  return powint((rcut - r), (alpha + 2)) /
      sqrt(2 * powint(rcut, (2 * alpha + 7)) / (2 * alpha + 5) / (2 * alpha + 6) / (2 * alpha + 7));
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::compute_g(double r, int n, int nmax, double rcut, double *w, int lw1)
{
  double Sum;
  Sum = 0.0;
  int alpha;

  for (alpha = 1; alpha < nmax + 1; alpha++)
    Sum += w[(n - 1) * lw1 + alpha - 1] * phi(r, alpha, rcut);

  return Sum;
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::Cosine(double Rij, double Rc)
{

  return 0.5 * (cos(MY_PI * Rij / Rc) + 1.0);
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::CosinePrime(double Rij, double Rc)
{

  return -0.5 * MY_PI / Rc * sin(MY_PI * Rij / Rc);
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::compute_sfac(double r, double rcut)
{

  if (r > rcut)
    return 0.0;
  else
    return Cosine(r, rcut);
}

/* ---------------------------------------------------------------------- */

double MLIAP_SO3::compute_dsfac(double r, double rcut)
{

  if (r > rcut)
    return 0.0;
  else
    return CosinePrime(r, rcut);
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::compute_dpidrj(int nmax, int lmax, double *clisttot_r, double *clisttot_i,
                               int lctot2, double *dclist_r, double *dclist_i, int ldcli2,
                               int ldcli3, double *dplist_r, int dpli2)
{
  double temp_r;
  double norm;
  int i, n1, n2, j, l, m, ii;
  i = 0;
  for (n1 = 0; n1 < nmax; n1++)
    for (n2 = 0; n2 < n1 + 1; n2++) {
      j = 0;
      for (l = 0; l < lmax + 1; l++) {
        norm = 2.0 * sqrt(2.0) * MY_PI / sqrt(2.0 * l + 1.0);
        for (m = -l; m < l + 1; m++) {
          for (ii = 0; ii < 3; ii++) {

            temp_r = dclist_r[(n1 * ldcli2 + j) * ldcli3 + ii] * clisttot_r[n2 * lctot2 + j] +
                dclist_i[(n1 * ldcli2 + j) * ldcli3 + ii] * clisttot_i[n2 * lctot2 + j];

            temp_r += clisttot_r[n1 * lctot2 + j] * dclist_r[(n2 * ldcli2 + j) * ldcli3 + ii] +
                clisttot_i[n1 * lctot2 + j] * dclist_i[(n2 * ldcli2 + j) * ldcli3 + ii];

            temp_r *= norm;

            dplist_r[i * dpli2 + ii] += temp_r;
          }
          j += 1;
        }
        i += 1;
      }
    }
}

/* ---------------------------------------------------------------------- */

int MLIAP_SO3::get_sum(int istart, int iend, int id, int imult)
{
  int ires = 0;
  int i;

  for (i = istart; i < iend; i = i + id) ires += i * imult;

  return ires;
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::compute_uarray_recursive(double x, double y, double z, double r, int twol,
                                         double *ulist_r, double *ulist_i, int *idxu_block,
                                         double *rootpqarray)
{
  int l, llu, llup, mb, ma, mbpar, mapar;
  double rootpq;
  int ldim = twol + 1;

  double theta, phi, atheta, btheta;

  double aphi_r, aphi_i, a_r, a_i, b_r, b_i;

  theta = acos(z / r);
  phi = atan2(y, x);

  atheta = cos(theta / 2);
  btheta = sin(theta / 2);

  aphi_r = cos(phi / 2);
  aphi_i = sin(phi / 2);

  a_r = atheta * aphi_r;
  a_i = atheta * aphi_i;
  b_r = btheta * aphi_r;
  b_i = btheta * aphi_i;

  ulist_r[0] = 1.0;
  ulist_i[0] = 0.0;

  for (l = 1; l < ldim; l++) {

    llu = idxu_block[l];
    llup = idxu_block[l - 1];
    mb = 0;

    while (2 * mb <= l) {

      ulist_r[llu] = 0.0;
      ulist_i[llu] = 0.0;
      for (ma = 0; ma < l; ma++) {

        rootpq = rootpqarray[(l - ma) * ldim + l - mb];

        ulist_r[llu] += rootpq * (a_r * ulist_r[llup] + a_i * ulist_i[llup]);
        ulist_i[llu] += rootpq * (a_r * ulist_i[llup] - a_i * ulist_r[llup]);

        rootpq = rootpqarray[(ma + 1) * ldim + l - mb];

        ulist_r[llu + 1] += -rootpq * (b_r * ulist_r[llup] + b_i * ulist_i[llup]);
        ulist_i[llu + 1] += -rootpq * (b_r * ulist_i[llup] - b_i * ulist_r[llup]);

        llu += 1;
        llup += 1;
      }

      llu += 1;
      mb += 1;
    }

    llu = idxu_block[l];
    llup = llu + (l + 1) * (l + 1) - 1;
    mbpar = 1;
    mb = 0;

    while (2 * mb <= l) {
      mapar = mbpar;
      for (ma = 0; ma < l + 1; ma++) {
        if (mapar == 1) {

          ulist_r[llup] = ulist_r[llu];
          ulist_i[llup] = -ulist_i[llu];

        } else {

          ulist_r[llup] = -ulist_r[llu];
          ulist_i[llup] = ulist_i[llu];
        }
        mapar = -mapar;
        llu += 1;
        llup -= 1;
      }
      mbpar = -mbpar;
      mb += 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::init_garray(int nmax, int lmax, double rcut, double alpha, double *w, int lw1,
                            double *g_array, int lg2)
{
  int i, n, Nmax = (nmax + lmax + 1) * 10;
  double x, xi;

  for (i = 1; i < Nmax + 1; i++) {
    // roots of Chebyshev polynomial of degree N
    x = cos((2 * i - 1) * MY_PI / 2 / Nmax);
    // transform the interval [-1,1] to [0, rcut]
    xi = rcut / 2 * (x + 1);
    for (n = 1; n < nmax + 1; n++)
      // r**2*g(n)(r)*e^(-alpha*r**2)
      g_array[(n - 1) * lg2 + i - 1] = rcut / 2 * MY_PI / Nmax * sqrt(1 - x * x) * xi * xi *
          compute_g(xi, n, nmax, rcut, w, lw1) * exp(-alpha * xi * xi);
  }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::get_sbes_array(int nlocal, int *numneighs, double **rij, int lmax, double rcut,
                               double alpha)
{
  int i, j;
  bigint neighbor;

  double x, y, z, ri, xi, rb;
  double sa, sb;
  double pfac1, pfac2, pfac3, pfac4;
  double exts;

  pfac1 = alpha * rcut;
  pfac4 = rcut / 2.0;
  pfac3 = MY_PI / 2.0 / m_Nmax;

  bigint ipair = 0;
  bigint gindex;
  int findex = m_Nmax * (m_lmax + 1);
  const bigint mindex = m_lmax + 1;

  for (int ii = 0; ii < nlocal; ii++) {

    for (neighbor = 0; neighbor < numneighs[ii]; neighbor++) {

      x = rij[ipair][0];
      y = rij[ipair][1];
      z = rij[ipair][2];
      ipair++;

      ri = sqrt(x * x + y * y + z * z);

      if (ri < SMALL) continue;

      pfac2 = pfac1 * ri;
      gindex = (ipair - 1) * findex;

      for (i = 1; i < m_Nmax + 1; i++) {
        const bigint i1mindex = (bigint) (i - 1) * mindex;

        x = cos((2 * i - 1) * pfac3);
        xi = pfac4 * (x + 1);
        rb = pfac2 * (x + 1);

        sa = sinh(rb) / rb;
        sb = (cosh(rb) - sa) / rb;

        m_sbes_array[gindex + i1mindex + 0] = sa;
        m_sbes_array[gindex + i1mindex + 1] = sb;

        for (j = 2; j < lmax + 1; j++)
          m_sbes_array[gindex + i1mindex + j] = m_sbes_array[gindex + i1mindex + j - 2] -
              (2 * j - 1) / rb * m_sbes_array[gindex + i1mindex + j - 1];

        exts = m_sbes_array[gindex + i1mindex + j - 2] -
            (2 * j - 1) / rb * m_sbes_array[gindex + i1mindex + j - 1];

        m_sbes_darray[gindex + i1mindex + 0] = sb;

        for (j = 1; j < lmax; j++)
          m_sbes_darray[gindex + i1mindex + j] = xi *
              (j * m_sbes_array[gindex + i1mindex + j - 1] +
               (j + 1) * m_sbes_array[gindex + i1mindex + j + 1]) /
              (2 * j + 1);

        m_sbes_darray[gindex + i1mindex + j] =
            xi * (j * m_sbes_array[gindex + i1mindex + j - 1] + (j + 1) * exts) / (2 * j + 1);
        m_sbes_darray[gindex + i1mindex + 0] = xi * sb;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::get_rip_array(int nlocal, int *numneighs, double **rij, int nmax, int lmax,
                              double alpha)
{
  int i, l, n;
  double integrald, integral = 0.0;
  bigint neighbor;

  double x, y, z, ri, expfac;

  bigint ipair = 0;

  for (int ii = 0; ii < nlocal; ii++)

    for (neighbor = 0; neighbor < numneighs[ii]; neighbor++) {

      x = rij[ipair][0];
      y = rij[ipair][1];
      z = rij[ipair][2];
      ipair++;

      ri = sqrt(x * x + y * y + z * z);

      if (ri < SMALL) continue;

      expfac = 4 * MY_PI * exp(-alpha * ri * ri);

      for (n = 1; n < nmax + 1; n++)
        for (l = 0; l < lmax + 1; l++) {

          integral = 0.0;
          integrald = 0.0;
          for (i = 0; i < m_Nmax; i++) {
            integral += m_g_array[(n - 1) * m_Nmax + i] *
                m_sbes_array[(ipair - 1) * m_Nmax * (m_lmax + 1) + (bigint) i * (m_lmax + 1) + l];
            integrald += m_g_array[(n - 1) * m_Nmax + i] *
                m_sbes_darray[(ipair - 1) * m_Nmax * (m_lmax + 1) + (bigint) i * (m_lmax + 1) + l];
          }

          m_rip_array[(ipair - 1) * m_nmax * (m_lmax + 1) + (bigint) (n - 1) * (m_lmax + 1) + l] =
              integral * expfac;
          m_rip_darray[(ipair - 1) * m_nmax * (m_lmax + 1) + (bigint) (n - 1) * (m_lmax + 1) + l] =
              integrald * expfac;
        }
    }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::spectrum(int nlocal, int *numneighs, int *jelems, double *wjelem, double **rij,
                         int nmax, int lmax, double rcut, double alpha, int ncoefs)
{
  init_arrays(nlocal, ncoefs);

  bigint totaln = 0;
  bigint totali;
  double Ylm_r, Ylm_i;
  int weight, neighbor;
  double x, y, z, r;
  double r_int;
  int twolmax = 2 * (lmax + 1);
  int findex;
  bigint gindex;
  bigint ipair = 0;
  double sfac;

  findex = m_nmax * (m_lmax + 1);

  for (int i = 0; i < nlocal; i++) totaln += numneighs[i];

  totali = totaln * m_Nmax * (m_lmax + 1);
  memory->destroy(m_sbes_array);
  memory->create(m_sbes_array, totali, "MLIAP_SO3:m_sbes_array");
  memory->destroy(m_sbes_darray);
  memory->create(m_sbes_darray, totali, "MLIAP_SO3:m_sbes_darray");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = totaln * m_nmax * (m_lmax + 1);
  memory->destroy(m_rip_array);
  memory->create(m_rip_array, totali, "MLIAP_SO3:m_rip_array");
  memory->destroy(m_rip_darray);
  memory->create(m_rip_darray, totali, "MLIAP_SO3:m_rip_darray");
  alloc_arrays += 2.0 * totali * sizeof(double);

  totali = totaln * ncoefs * 3;
  memory->destroy(m_dplist_r);
  memory->create(m_dplist_r, totali, "MLIAP_SO3:m_dplist_r");
  memory->destroy(m_dplist_i);
  memory->create(m_dplist_i, totali, "MLIAP_SO3:m_dplist_i");
  alloc_arrays += 2.0 * totali * sizeof(double);

  get_sbes_array(nlocal, numneighs, rij, lmax, rcut, alpha);

  get_rip_array(nlocal, numneighs, rij, nmax, lmax, alpha);

  totali = (bigint) nlocal * ncoefs;
  for (int i = 0; i < totali; i++) {
    m_plist_r[i] = 0.0;
    m_plist_i[i] = 0.0;
  }

  for (int ii = 0; ii < nlocal; ii++) {
    totali = (bigint) nmax * m_numYlms;

    for (bigint ti = 0; ti < totali; ti++) {
      m_clisttot_r[ti] = 0.0;
      m_clisttot_i[ti] = 0.0;
    }

    for (neighbor = 0; neighbor < numneighs[ii]; neighbor++) {
      const int jelem = jelems[ipair];
      weight = wjelem[jelem];

      x = rij[ipair][0];
      y = rij[ipair][1];
      z = rij[ipair][2];
      ipair++;

      r = sqrt(x * x + y * y + z * z);

      if (r < SMALL) continue;
      totali = (bigint) nmax * m_numYlms;
      for (bigint ti = 0; ti < totali; ti++) {
        m_clist_r[ti] = 0.0;
        m_clist_i[ti] = 0.0;
      }
      for (bigint ti = 0; ti < m_idxu_count; ti++) {
        m_ulist_r[ti] = 0.0;
        m_ulist_i[ti] = 0.0;
      }

      compute_uarray_recursive(x, y, z, r, twolmax, m_ulist_r, m_ulist_i, m_idxu_block, m_rootpq);

      sfac = compute_sfac(r, rcut);

      gindex = (ipair - 1) * findex;
      for (int n = 1; n < nmax + 1; n++) {
        int i = 0;
        for (int l = 0; l < lmax + 1; l++) {
          r_int = m_rip_array[gindex + (bigint) (n - 1) * (m_lmax + 1) + l];

          for (int m = -l; m < l + 1; m++) {

            Ylm_r = (m_ulist_r[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            m_clist_r[(n - 1) * m_numYlms + i] += r_int * Ylm_r * sfac;
            Ylm_i = (m_ulist_i[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            m_clist_i[(n - 1) * m_numYlms + i] += r_int * Ylm_i * sfac;
            i += 1;
          }
        }
      }

      totali = (bigint) nmax * m_numYlms;
      for (bigint tn = 0; tn < totali; tn++) {
        m_clist_r[tn] = m_clist_r[tn] * double(weight);
        m_clist_i[tn] = m_clist_i[tn] * double(weight);
      }

      for (bigint tn = 0; tn < totali; tn++) {
        m_clisttot_r[tn] += m_clist_r[tn];
        m_clisttot_i[tn] += m_clist_i[tn];
      }
    }
    compute_pi(nmax, lmax, m_clisttot_r, m_clisttot_i, m_numYlms, m_plist_r, m_plist_i, ncoefs, ii);
  }
}

/* ---------------------------------------------------------------------- */

void MLIAP_SO3::spectrum_dxdr(int nlocal, int *numneighs, int *jelems, double *wjelem, double **rij,
                              int nmax, int lmax, double rcut, double alpha, bigint npairs,
                              int ncoefs)
{
  bigint totaln = 0;
  bigint totali;
  double dr_int[3];
  double Ylm_r, Ylm_i;

  double rvec[3];
  double dexpfac[3];
  double dfact[6];

  double xcov0_r, xcov0_i, xcovpl1_r, xcovpl1_i, xcovm1_r, xcovm1_i;
  double comj_i;
  double r_int;
  double r_int_temp;

  double oneofr;
  int findex;
  bigint gindex;

  int numps, weight, neighbor;

  double x, y, z, r;

  bigint ipair = 0;
  bigint idpair = 0;
  double sfac, dsfac, dsfac_arr[3];

  findex = m_nmax * (m_lmax + 1);

  for (int i = 0; i < nlocal; i++) totaln += numneighs[i];

  totali = totaln * m_Nmax * (m_lmax + 1);
  memory->destroy(m_sbes_array);
  memory->create(m_sbes_array, totali, "MLIAP_SO3:m_sbes_array");
  memory->destroy(m_sbes_darray);
  memory->create(m_sbes_darray, totali, "MLIAP_SO3:m_sbes_darray");

  totali = totaln * m_nmax * (m_lmax + 1);
  memory->destroy(m_rip_array);
  memory->create(m_rip_array, totali, "MLIAP_SO3:m_rip_array");
  memory->destroy(m_rip_darray);
  memory->create(m_rip_darray, totali, "MLIAP_SO3:m_rip_darray");

  totali = totaln * ncoefs * 3;
  memory->destroy(m_dplist_r);
  memory->create(m_dplist_r, totali, "MLIAP_SO3:m_dplist_r");
  memory->destroy(m_dplist_i);
  memory->create(m_dplist_i, totali, "MLIAP_SO3:m_dplist_i");

  totali = npairs * ncoefs * 3;

  for (int i = 0; i < totali; i++) {
    m_dplist_r[i] = 0.0;
    m_dplist_i[i] = 0.0;
  }

  numps = nmax * (nmax + 1) * (lmax + 1) / 2;

  get_sbes_array(nlocal, numneighs, rij, lmax, rcut, alpha);

  get_rip_array(nlocal, numneighs, rij, nmax, lmax, alpha);

  int twolmax = 2 * (lmax + 1);

  for (int ii = 0; ii < nlocal; ii++) {

    totali = (bigint) nmax * m_numYlms;
    for (bigint ti = 0; ti < totali; ti++) {
      m_clisttot_r[ti] = 0.0;
      m_clisttot_i[ti] = 0.0;
    }

    for (neighbor = 0; neighbor < numneighs[ii]; neighbor++) {

      const int jelem = jelems[ipair];
      weight = wjelem[jelem];

      x = rij[ipair][0];
      y = rij[ipair][1];
      z = rij[ipair][2];
      ipair++;

      r = sqrt(x * x + y * y + z * z);

      if (r < SMALL) continue;
      totali = (bigint) nmax * m_numYlms;

      for (bigint ti = 0; ti < totali; ti++) {
        m_clist_r[ti] = 0.0;
        m_clist_i[ti] = 0.0;
      }

      for (bigint ti = 0; ti < m_idxu_count; ti++) {
        m_ulist_r[ti] = 0.0;
        m_ulist_i[ti] = 0.0;
      }

      compute_uarray_recursive(x, y, z, r, twolmax, m_ulist_r, m_ulist_i, m_idxu_block, m_rootpq);

      sfac = compute_sfac(r, rcut);

      gindex = (ipair - 1) * findex;
      for (int n = 1; n < nmax + 1; n++) {
        int i = 0;
        for (int l = 0; l < lmax + 1; l++) {
          r_int = m_rip_array[gindex + (bigint) (n - 1) * (m_lmax + 1) + l];

          for (int m = -l; m < l + 1; m++) {

            Ylm_r = (m_ulist_r[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            m_clist_r[(n - 1) * m_numYlms + i] += r_int * Ylm_r * sfac;
            Ylm_i = (m_ulist_i[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            m_clist_i[(n - 1) * m_numYlms + i] += r_int * Ylm_i * sfac;
            i += 1;
          }
        }
      }

      totali = (bigint) nmax * m_numYlms;
      for (bigint tn = 0; tn < totali; tn++) {
        m_clist_r[tn] = m_clist_r[tn] * double(weight);
        m_clist_i[tn] = m_clist_i[tn] * double(weight);
      }

      for (bigint tn = 0; tn < totali; tn++) {
        m_clisttot_r[tn] += m_clist_r[tn];
        m_clisttot_i[tn] += m_clist_i[tn];
      }
    }

    for (neighbor = 0; neighbor < numneighs[ii]; neighbor++) {

      const int jelem = jelems[idpair];
      weight = wjelem[jelem];

      x = rij[idpair][0];
      y = rij[idpair][1];
      z = rij[idpair][2];
      idpair++;

      r = sqrt(x * x + y * y + z * z);
      if (r < SMALL) continue;

      totali = (bigint) nmax * m_numYlms * 3;
      for (bigint tn = 0; tn < totali; tn++) {
        m_dclist_r[tn] = 0.0;
        m_dclist_i[tn] = 0.0;
      }

      for (int ti = 0; ti < m_idxu_count; ti++) {
        m_ulist_r[ti] = 0.0;
        m_ulist_i[ti] = 0.0;
      }

      compute_uarray_recursive(x, y, z, r, twolmax, m_ulist_r, m_ulist_i, m_idxu_block, m_rootpq);

      /////////  compute_carray_wD  ////////
      {
        rvec[0] = x;
        rvec[1] = y;
        rvec[2] = z;
        totali = ((bigint) lmax + 2) * (lmax + 2);
        for (bigint tn = 0; tn < totali; tn++) {
          m_Ylms_r[tn] = 0.0;
          m_Ylms_i[tn] = 0.0;
        }

        int i = 0;
        for (int l = 0; l < lmax + 2; l++) {
          for (int m = -l; m < l + 1; m++) {
            m_Ylms_r[i] = (m_ulist_r[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            m_Ylms_i[i] = (m_ulist_i[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
            i += 1;
          }
        }

        totali = ((bigint) lmax + 1) * (lmax + 1) * 3;
        for (bigint tn = 0; tn < totali; tn++) {
          m_dYlm_r[tn] = 0.0;
          m_dYlm_i[tn] = 0.0;
        }

        comj_i = 1.0 / sqrt(2.0);
        oneofr = 1.0 / r;

        i = 1;
        for (int l = 1; l < lmax + 1; l++) {
          for (int m = -l; m < l + 1; m++) {

            dfact[0] = m_dfac0[l * m_dfac_l2 + m] * oneofr;
            dfact[1] = m_dfac1[l * m_dfac_l2 + m] * oneofr;
            dfact[2] = m_dfac2[l * m_dfac_l2 + m] * oneofr;
            dfact[3] = m_dfac3[l * m_dfac_l2 + m] * oneofr;
            dfact[4] = m_dfac4[l * m_dfac_l2 + m] * oneofr;
            dfact[5] = m_dfac5[l * m_dfac_l2 + m] * oneofr;

            xcov0_r = dfact[0] * m_Ylms_r[m_ellpl1[l] + m];
            xcov0_i = dfact[0] * m_Ylms_i[m_ellpl1[l] + m];
            if (abs(m) <= l - 1.0) {
              xcov0_r += dfact[1] * m_Ylms_r[m_ellm1[l] + m];
              xcov0_i += dfact[1] * m_Ylms_i[m_ellm1[l] + m];
            }
            xcovpl1_r = dfact[2] * m_Ylms_r[m_ellpl1[l] + m + 1];
            xcovpl1_i = dfact[2] * m_Ylms_i[m_ellpl1[l] + m + 1];
            if (abs(m + 1) <= l - 1.0) {
              xcovpl1_r -= dfact[3] * m_Ylms_r[m_ellm1[l] + m + 1];
              xcovpl1_i -= dfact[3] * m_Ylms_i[m_ellm1[l] + m + 1];
            }
            xcovm1_r = dfact[4] * m_Ylms_r[m_ellpl1[l] + m - 1];
            xcovm1_i = dfact[4] * m_Ylms_i[m_ellpl1[l] + m - 1];
            if (fabs(m - 1.0) <= l - 1.0) {
              xcovm1_r -= dfact[5] * m_Ylms_r[m_ellm1[l] + m - 1];
              xcovm1_i -= dfact[5] * m_Ylms_i[m_ellm1[l] + m - 1];
            }
            m_dYlm_r[i * 3 + 0] = 1.0 / sqrt(2.0) * (xcovm1_r - xcovpl1_r);
            m_dYlm_r[i * 3 + 1] = -comj_i * (xcovm1_i + xcovpl1_i);
            m_dYlm_r[i * 3 + 2] = xcov0_r;

            m_dYlm_i[i * 3 + 0] = 1.0 / sqrt(2.0) * (xcovm1_i - xcovpl1_i);
            m_dYlm_i[i * 3 + 1] = comj_i * (xcovm1_r + xcovpl1_r);
            m_dYlm_i[i * 3 + 2] = xcov0_i;

            i += 1;
          }
        }

        for (int ii = 0; ii < 3; ii++) dexpfac[ii] = -2.0 * alpha * rvec[ii];

        sfac = compute_sfac(r, rcut);
        dsfac = compute_dsfac(r, rcut);
        for (int ii = 0; ii < 3; ii++) { dsfac_arr[ii] = dsfac * rvec[ii] / r; }

        for (int n = 1; n < nmax + 1; n++) {
          int i = 0;
          for (int l = 0; l < lmax + 1; l++) {
            r_int = m_rip_array[(idpair - 1) * m_nmax * (m_lmax + 1) +
                                (bigint) (n - 1) * (m_lmax + 1) + l];
            r_int_temp = m_rip_darray[(idpair - 1) * m_nmax * (m_lmax + 1) +
                                      (bigint) (n - 1) * (m_lmax + 1) + l];

            for (int ii = 0; ii < 3; ii++) dr_int[ii] = r_int_temp * 2.0 * alpha * rvec[ii] / r;

            for (int m = -l; m < l + 1; m++) {

              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 0] +=
                  (r_int * m_Ylms_r[i] * dexpfac[0] + dr_int[0] * m_Ylms_r[i] +
                   r_int * m_dYlm_r[i * 3 + 0]) *
                  sfac;
              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 1] +=
                  (r_int * m_Ylms_r[i] * dexpfac[1] + dr_int[1] * m_Ylms_r[i] +
                   r_int * m_dYlm_r[i * 3 + 1]) *
                  sfac;
              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 2] +=
                  (r_int * m_Ylms_r[i] * dexpfac[2] + dr_int[2] * m_Ylms_r[i] +
                   r_int * m_dYlm_r[i * 3 + 2]) *
                  sfac;

              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 0] +=
                  (r_int * m_Ylms_i[i] * dexpfac[0] + dr_int[0] * m_Ylms_i[i] +
                   r_int * m_dYlm_i[i * 3 + 0]) *
                  sfac;
              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 1] +=
                  (r_int * m_Ylms_i[i] * dexpfac[1] + dr_int[1] * m_Ylms_i[i] +
                   r_int * m_dYlm_i[i * 3 + 1]) *
                  sfac;
              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 2] +=
                  (r_int * m_Ylms_i[i] * dexpfac[2] + dr_int[2] * m_Ylms_i[i] +
                   r_int * m_dYlm_i[i * 3 + 2]) *
                  sfac;

              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 0] += (r_int * m_Ylms_r[i]) * dsfac_arr[0];
              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 1] += (r_int * m_Ylms_r[i]) * dsfac_arr[1];
              m_dclist_r[((n - 1) * m_numYlms + i) * 3 + 2] += (r_int * m_Ylms_r[i]) * dsfac_arr[2];

              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 0] += (r_int * m_Ylms_i[i]) * dsfac_arr[0];
              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 1] += (r_int * m_Ylms_i[i]) * dsfac_arr[1];
              m_dclist_i[((n - 1) * m_numYlms + i) * 3 + 2] += (r_int * m_Ylms_i[i]) * dsfac_arr[2];

              i += 1;
            }
          }
        }
      }
      /////// end compute_carray_wD //////////////////

      totali = (bigint) nmax * m_numYlms * 3;
      for (bigint tn = 0; tn < totali; tn++) {
        m_dclist_r[tn] = m_dclist_r[tn] * double(weight);
        m_dclist_i[tn] = m_dclist_i[tn] * double(weight);
      }

      totali = (bigint) numps * 3;
      for (bigint ti = 0; ti < totali; ti++) m_tempdp_r[ti] = 0.0;

      compute_dpidrj(nmax, lmax, m_clisttot_r, m_clisttot_i, m_numYlms, m_dclist_r, m_dclist_i,
                     m_numYlms, 3, m_tempdp_r, 3);

      for (bigint tn = 0; tn < totali; tn++)
        m_dplist_r[((idpair - 1) * (numps * 3)) + tn] += m_tempdp_r[tn];

    }    //for(neighbor=0;neighbor<numneighs[ii];neighbor++){

  }    //for (int ii = 0; ii < nlocal; ii++) {
}
