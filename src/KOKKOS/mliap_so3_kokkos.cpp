// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
 ------------------------------------------------------------------------- */

#include "mliap_so3_kokkos.h"

#include "error.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "mliap_so3_math.h"

#include <cmath>

using namespace SO3Math;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecialKokkos;

#define SMALL 1.0e-8

/* ---------------------------------------------------------------------- */

template <class DeviceType>
MLIAP_SO3Kokkos<DeviceType>::MLIAP_SO3Kokkos(LAMMPS *lmp, double vrcut, int vlmax, int vnmax, double valpha) : Pointers(lmp)
{
  m_rcut = vrcut;
  m_alpha = valpha;
  m_lmax = vlmax;
  m_nmax = vnmax;
  compute_ncoeff();

  m_Nmax = (m_nmax + m_lmax + 1) * 10;
  m_numYlms = (m_lmax + 1) * (m_lmax + 1);

  m_init_arrays = 0;
  m_dfac_l1 = m_dfac_l2 = 0;
  m_pfac_l1 = m_pfac_l2 = 0;
  m_idxu_count = m_idxy_count = 0;
  alloc_init = alloc_arrays = 0.0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
MLIAP_SO3Kokkos<DeviceType>::~MLIAP_SO3Kokkos()
{
  memoryKK->destroy_kokkos(m_ellpl1);
  memoryKK->destroy_kokkos(m_ellm1);
  memoryKK->destroy_kokkos(m_pfac);
  memoryKK->destroy_kokkos(m_Ylms);
  memoryKK->destroy_kokkos(m_dfac0);
  memoryKK->destroy_kokkos(m_dfac1);
  memoryKK->destroy_kokkos(m_dfac2);
  memoryKK->destroy_kokkos(m_dfac3);
  memoryKK->destroy_kokkos(m_dfac4);
  memoryKK->destroy_kokkos(m_dfac5);
  memoryKK->destroy_kokkos(m_w);
  memoryKK->destroy_kokkos(m_g_array);

  memoryKK->destroy_kokkos(m_rootpq);
  memoryKK->destroy_kokkos(m_idxu_block);
  memoryKK->destroy_kokkos(m_idxylm);

  memoryKK->destroy_kokkos(m_rip_array);
  memoryKK->destroy_kokkos(m_rip_darray);

  memoryKK->destroy_kokkos(m_sbes_array);
  memoryKK->destroy_kokkos(m_sbes_darray);

  memoryKK->destroy_kokkos(m_plist_r);

  memoryKK->destroy_kokkos(m_ulist_r);
  memoryKK->destroy_kokkos(m_ulist_i);

  memoryKK->destroy_kokkos(m_dYlm_r);
  memoryKK->destroy_kokkos(m_dYlm_i);

  memoryKK->destroy_kokkos(k_dplist_r);

  memoryKK->destroy_kokkos(m_dclist);

  memoryKK->destroy_kokkos(m_clisttot_r);
  memoryKK->destroy_kokkos(m_clisttot_i);

  t_numneighs = int_1d();
  t_jelems = int_1d();
  t_wjelem = float_1d();
  t_rij = float_2d();
  t_ij = int_1d();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::compute_ncoeff()
{
  ncoeff = m_nmax * (m_nmax + 1) * (m_lmax + 1) / 2;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::init()
{
  int totali;

  totali = m_lmax + 1;
  memoryKK->destroy_kokkos(m_ellpl1);
  memoryKK->create_kokkos(m_ellpl1, totali, "MLIAP_SO3Kokkos:m_ellpl1");
  memoryKK->destroy_kokkos(m_ellm1);
  memoryKK->create_kokkos(m_ellm1, totali, "MLIAP_SO3Kokkos:m_ellm1");
  alloc_init = 2.0 * totali * sizeof(double);
  using range=Kokkos::RangePolicy<DeviceType>;
  auto ellpl1 = m_ellpl1, ellm1 = m_ellm1;

  Kokkos::parallel_for(range(0,m_lmax), KOKKOS_LAMBDA (int ll) {
    int l=ll+1;
    ellpl1[l] = get_sum(0, l + 2, 1, 2);
    ellm1[l]  = get_sum(0, l, 1, 2);
  });

  double pfac1 = 1.0 / 4.0 / MY_PI;
  m_pfac_l1 = m_lmax + 2;
  m_pfac_l2 = (m_lmax + 2) * (m_lmax + 2) + 1;
  totali = m_pfac_l1 * m_pfac_l2;
  memoryKK->destroy_kokkos(m_pfac);
  memoryKK->create_kokkos(m_pfac, totali, "MLIAP_SO3Kokkos:m_pfac");
  memoryKK->destroy_kokkos(m_Ylms);
  memoryKK->create_kokkos(m_Ylms, totali, "MLIAP_SO3Kokkos:m_Ylms");
  alloc_init += 2 * totali * sizeof(double);

  auto pfac = m_pfac, Ylms = m_Ylms;
  auto pfac_l2 = m_pfac_l2, lmax = m_lmax;
  // Serial but just to make sure run with device memory
  Kokkos::parallel_for(range(0,1), KOKKOS_LAMBDA (int ) {
    int i=0;
    for (int l = 0; l < lmax + 2; l++)
      for (int m = -l; m < l + 1; m++) {
        pfac[l * pfac_l2 + m] = sqrt((2.0 * l + 1.0) * pfac1) * powsign(m);
        Ylms[i] = pfac[l * pfac_l2 + m];
        i += 1;
      }
  });

  m_dfac_l1 = m_lmax + 1;
  m_dfac_l2 = m_numYlms + 1;
  totali = m_dfac_l1 * m_dfac_l2;
  memoryKK->destroy_kokkos(m_dfac0);
  memoryKK->create_kokkos(m_dfac0, totali, "MLIAP_SO3Kokkos:m_dfac0");
  memoryKK->destroy_kokkos(m_dfac1);
  memoryKK->create_kokkos(m_dfac1, totali, "MLIAP_SO3Kokkos:m_dfac1");
  memoryKK->destroy_kokkos(m_dfac2);
  memoryKK->create_kokkos(m_dfac2, totali, "MLIAP_SO3Kokkos:m_dfac2");
  memoryKK->destroy_kokkos(m_dfac3);
  memoryKK->create_kokkos(m_dfac3, totali, "MLIAP_SO3Kokkos:m_dfac3");
  memoryKK->destroy_kokkos(m_dfac4);
  memoryKK->create_kokkos(m_dfac4, totali, "MLIAP_SO3Kokkos:m_dfac4");
  memoryKK->destroy_kokkos(m_dfac5);
  memoryKK->create_kokkos(m_dfac5, totali, "MLIAP_SO3Kokkos:m_dfac5");
  alloc_init += 6.0 * totali * sizeof(double);

  auto dfac0 = m_dfac0,dfac1 = m_dfac1,dfac2 = m_dfac2,dfac3 = m_dfac3,dfac4 = m_dfac4,dfac5 = m_dfac5;
  auto dfac_l2 = m_dfac_l2;

  Kokkos::parallel_for(range(0,m_lmax), KOKKOS_LAMBDA (int ll) {
    int l = ll+1;
    for (int m = -l; m < l + 1; m++) {
      dfac0[l * dfac_l2 + m] =
          -sqrt(((l + 1.0) * (l + 1.0) - m * m) / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      dfac1[l * dfac_l2 + m] =
          sqrt((l * l - m * m) / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
      dfac2[l * dfac_l2 + m] =
          -sqrt((l + m + 1.0) * (l + m + 2.0) / 2.0 / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      dfac3[l * dfac_l2 + m] =
          sqrt((l - m - 1.0) * (l - m) / 2.0 / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
      dfac4[l * dfac_l2 + m] =
          -sqrt((l - m + 1.0) * (l - m + 2.0) / 2.0 / (2.0 * l + 1.0) / (2.0 * l + 3.0)) * l;
      dfac5[l * dfac_l2 + m] =
          sqrt((l + m - 1.0) * (l + m) / 2.0 / (2.0 * l - 1.0) / (2.0 * l + 1.0)) * (l + 1.0);
    }
  });

  totali = m_nmax * m_nmax;
  memoryKK->destroy_kokkos(m_w);
  memoryKK->create_kokkos(m_w, totali, "MLIAP_SO3Kokkos:w");
  alloc_init += totali * sizeof(double);

  totali = m_nmax * m_Nmax;
  memoryKK->destroy_kokkos(m_g_array);
  memoryKK->create_kokkos(m_g_array, totali, "MLIAP_SO3Kokkos:g_array");
  alloc_init += totali * sizeof(double);

  {
    Kokkos::View<double*, LMPHostType> w("w", m_nmax * m_nmax), g_array("g_array", m_nmax * m_Nmax);
    compute_W(m_nmax, w.data());
    init_garray(m_nmax, m_lmax, m_rcut, m_alpha, w.data(), m_nmax, g_array.data(), m_Nmax);
    Kokkos::deep_copy(m_w,w);
    Kokkos::deep_copy(m_g_array,g_array);
  }

  int twolmax;
  twolmax = 2 * (m_lmax + 1);
  int m_ldim = twolmax + 1;
  totali = m_ldim * m_ldim;
  memoryKK->destroy_kokkos(m_rootpq);
  memoryKK->create_kokkos(m_rootpq, totali, "MLIAP_SO3Kokkos:rootpq");
  alloc_init += totali * sizeof(double);

  auto rootpq=m_rootpq;
  auto ldim=m_ldim;
  Kokkos::parallel_for(range(1,m_ldim), KOKKOS_LAMBDA (int p) {
    for (int q = 1; q < ldim; q++)
      rootpq[p * ldim + q] = sqrt(static_cast<double>(p) / q);
  });

  memoryKK->destroy_kokkos(m_idxu_block);
  memoryKK->create_kokkos(m_idxu_block, m_ldim, "MLIAP_SO3Kokkos:idxu_bloc");
  alloc_init += totali * sizeof(double);

  totali = square(m_lmax + 2);
  memoryKK->destroy_kokkos(m_idxylm);
  memoryKK->create_kokkos(m_idxylm, totali, "MLIAP_SO3Kokkos:idxylm");
  alloc_init += totali * sizeof(double);

  Kokkos::View<int*, LMPHostType> idxu_block("idxu_block",m_ldim);
  Kokkos::View<int*, LMPHostType> idxylm("idxylm",totali);
  m_idxu_count = m_idxy_count = 0;
  for (int l = 0; l < m_ldim; l++) {
    idxu_block[l] = m_idxu_count;
    for (int mb = 0; mb < l + 1; mb++)
      for (int ma = 0; ma < l + 1; ma++) {
        if (l % 2 == 0 && ma == l / 2) {
          idxylm[m_idxy_count] = m_idxu_count;
          m_idxy_count += 1;
        }
        m_idxu_count += 1;
      }
  }
  Kokkos::deep_copy(m_idxu_block,idxu_block);
  Kokkos::deep_copy(m_idxylm,idxylm);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::init_arrays(int nlocal, int ncoefs)
{

  {
    // sum up the temp memory space per particle
    int sum_of_temps=0;
    //ulist_r and i
    sum_of_temps += m_idxu_count*2;
    //m_dYlm_r
    sum_of_temps += m_numYlms * 3 * 2;
    //m_clisttot_r
    sum_of_temps += m_nmax * m_numYlms * 2;
    // k_dclist
    sum_of_temps += m_nmax * m_numYlms * 3 * 2;
    m_chunk_size = m_temp_memory_size/(sum_of_temps * sizeof(double));
  }

  int totali = nlocal * ncoefs;
  if ( nlocal > (int)m_plist_r.extent(0)) {
    memoryKK->destroy_kokkos(m_plist_r);
    memoryKK->create_kokkos(m_plist_r, nlocal, ncoefs, "MLIAP_SO3Kokkos:m_plist_r");
    alloc_arrays = totali * sizeof(double);
  }

  int num_of_temp = std::min(nlocal, m_chunk_size);
  if ((int)m_ulist_r.extent(0) < num_of_temp ) {
    totali = m_idxu_count;
    memoryKK->destroy_kokkos(m_ulist_r);
    memoryKK->create_kokkos(m_ulist_r, num_of_temp, totali, "MLIAP_SO3Kokkos:m_ulist_r");
    memoryKK->destroy_kokkos(m_ulist_i);
    memoryKK->create_kokkos(m_ulist_i, num_of_temp, totali, "MLIAP_SO3Kokkos:m_ulist_i");
    alloc_arrays += 2.0 * totali * num_of_temp * sizeof(double);

    totali = m_numYlms * 3;
    memoryKK->destroy_kokkos(m_dYlm_r);
    memoryKK->create_kokkos(m_dYlm_r, num_of_temp, m_numYlms, 3, "MLIAP_SO3Kokkos:m_dYlm_r");
    memoryKK->destroy_kokkos(m_dYlm_i);
    memoryKK->create_kokkos(m_dYlm_i, num_of_temp, m_numYlms, 3, "MLIAP_SO3Kokkos:m_dYlm_i");
    alloc_arrays += 2.0 * m_numYlms * 3 * num_of_temp * sizeof(double);

    memoryKK->destroy_kokkos(m_dclist);
    memoryKK->create_kokkos(m_dclist, num_of_temp, m_nmax, m_numYlms, 3, "MLIAP_SO3Kokkos:k_dclist_r");
    alloc_arrays += m_nmax * m_numYlms * 3 * num_of_temp* sizeof(double);

    memoryKK->destroy_kokkos(m_clisttot_r);
    memoryKK->create_kokkos(m_clisttot_r, num_of_temp, m_nmax, m_numYlms, "MLIAP_SO3Kokkos:m_clisttot_r");
    memoryKK->destroy_kokkos(m_clisttot_i);
    memoryKK->create_kokkos(m_clisttot_i, num_of_temp, m_nmax, m_numYlms, "MLIAP_SO3Kokkos:m_clisttot_i");
    alloc_arrays += 2.0 * m_nmax * m_numYlms * num_of_temp * sizeof(double);
    m_init_arrays = 1;
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double MLIAP_SO3Kokkos<DeviceType>::memory_usage()
{
  return alloc_init + alloc_arrays;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::compute_W(int nmax, double *arr)
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
template <class DeviceType>
template <typename ViewType>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::compute_pi(int nmax, int lmax, ViewType clisttot_r, ViewType clisttot_i, int /*lcl2*/,
                           float_2d plist_r, int indpl) const
{
  int n1, n2, j, l, m, i = 0;
  double norm;
  for (n1 = 0; n1 < nmax; n1++)
    for (n2 = 0; n2 < n1 + 1; n2++) {
      j = 0;
      for (l = 0; l < lmax + 1; l++) {
        norm = 2.0 * sqrt(2.0) * MY_PI / sqrt(2.0 * l + 1.0);

        for (m = -l; m < l + 1; m++) {

          plist_r(indpl, i) += (clisttot_r(n1, j) * clisttot_r(n2, j) +
                                clisttot_i(n1, j) * clisttot_i(n2, j)) *
              norm;
          j += 1;
        }
        i += 1;
      }
    }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double MLIAP_SO3Kokkos<DeviceType>::phi(double r, int alpha, double rcut)
{
  return powint((rcut - r), (alpha + 2)) /
      sqrt(2 * powint(rcut, (2 * alpha + 7)) / (2 * alpha + 5) / (2 * alpha + 6) / (2 * alpha + 7));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double MLIAP_SO3Kokkos<DeviceType>::compute_g(double r, int n, int nmax, double rcut, double *w, int lw1)
{
  double Sum;
  Sum = 0.0;
  int alpha;

  for (alpha = 1; alpha < nmax + 1; alpha++)
    Sum += w[(n - 1) * lw1 + alpha - 1] * phi(r, alpha, rcut);

  return Sum;
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double MLIAP_SO3Kokkos<DeviceType>::Cosine(double Rij, double Rc) const
{

  return 0.5 * (cos(MY_PI * Rij / Rc) + 1.0);
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double MLIAP_SO3Kokkos<DeviceType>::CosinePrime(double Rij, double Rc) const
{

  return -0.5 * MY_PI / Rc * sin(MY_PI * Rij / Rc);
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double MLIAP_SO3Kokkos<DeviceType>::compute_sfac(double r, double rcut) const
{
  if (r > rcut)
    return 0.0;
  else
    return Cosine(r, rcut);
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION
double MLIAP_SO3Kokkos<DeviceType>::compute_dsfac(double r, double rcut) const
{
  if (r > rcut)
    return 0.0;
  else
    return CosinePrime(r, rcut);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
int MLIAP_SO3Kokkos<DeviceType>::get_sum(int istart, int iend, int id, int imult)
{
  int ires = 0;
  int i;

  for (i = istart; i < iend; i = i + id) ires += i * imult;

  return ires;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <typename UlistView>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::compute_uarray_recursive(double x, double y, double z, double r, int twol,
                                               UlistView ulist_r, UlistView ulist_i, int_1d idxu_block,
                                               float_1d rootpqarray) const
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

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::init_garray(int nmax, int lmax, double rcut, double alpha, double *w, int lw1,
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

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::operator() (const MLIAPSO3GetSBESArrayTag&, int ii) const{
   int ipair = t_ij(ii);
   for (int neighbor = 0; neighbor < t_numneighs[ii]; neighbor++) {
     double x = t_rij(ipair, 0);
     double y = t_rij(ipair, 1);
     double z = t_rij(ipair, 2);
     ipair++;

     double ri = sqrt(x * x + y * y + z * z);

     if (ri < SMALL) continue;

     const double pfac1 = t_alpha * t_rcut;
     const double pfac4 = t_rcut / 2.0;
     const double pfac3 = MY_PI / 2.0 / m_Nmax;
     double pfac2 = pfac1 * ri;
     const int findex = m_Nmax * (m_lmax + 1);
     const int gindex = (ipair - 1) * findex;
     const int mindex = m_lmax + 1;

     for (int i = 1; i < m_Nmax + 1; i++) {
       const bigint i1mindex = (bigint) (i - 1) * mindex;

       x = cos((2 * i - 1) * pfac3);
       double xi = pfac4 * (x + 1);
       double rb = pfac2 * (x + 1);

       double sa = sinh(rb) / rb;
       double sb = (cosh(rb) - sa) / rb;

       m_sbes_array[gindex + i1mindex + 0] = sa;
       m_sbes_array[gindex + i1mindex + 1] = sb;

       int j;
       for (j = 2; j < t_lmax + 1; j++)
         m_sbes_array[gindex + i1mindex + j] = m_sbes_array[gindex + i1mindex + j - 2] -
             (2 * j - 1) / rb * m_sbes_array[gindex + i1mindex + j - 1];

       double exts = m_sbes_array[gindex + i1mindex + j - 2] -
           (2 * j - 1) / rb * m_sbes_array[gindex + i1mindex + j - 1];

       m_sbes_darray[gindex + i1mindex + 0] = sb;

       for (j = 1; j < t_lmax; j++)
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

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::operator() (const MLIAPSO3GetRipArrayTag&, int ii) const{
   int ipair = t_ij(ii);
   for (int neighbor = 0; neighbor < t_numneighs[ii]; neighbor++) {

     double x = t_rij(ipair, 0);
     double y = t_rij(ipair, 1);
     double z = t_rij(ipair, 2);
     ipair++;

     double ri = sqrt(x * x + y * y + z * z);
     if (ri < SMALL) continue;

     double expfac = 4 * MY_PI * exp(-t_alpha * ri * ri);
     for (int n = 1; n < t_nmax + 1; n++)
       for (int l = 0; l < t_lmax + 1; l++) {
         double integral = 0.0, integrald = 0.0;
         for (int i = 0; i < m_Nmax; i++) {
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

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::spectrum(int nlocal, DAT::tdual_int_1d numneighs, DAT::tdual_int_1d jelems, DAT::tdual_float_1d wjelem,
                               DAT::tdual_float_2d rij,  DAT::tdual_int_1d k_ij,
                               int nmax, int lmax, double rcut, double alpha, int totaln, int ncoefs)
{
  init_arrays(nlocal, ncoefs);

  bigint totali;

  totali = totaln * m_Nmax * (m_lmax + 1);
  if ( totali > (int)m_sbes_array.extent(0)) {
    memoryKK->realloc_kokkos(m_sbes_array, "MLIAP_SO3Kokkos:m_sbes_array", totali);
    memoryKK->realloc_kokkos(m_sbes_darray, "MLIAP_SO3Kokkos:m_sbes_darray", totali);
    alloc_arrays += 2.0 * totali * sizeof(double);
  }

  totali = totaln * m_nmax * (m_lmax + 1);
  if ( totali > (int)m_rip_array.extent(0)) {
    memoryKK->realloc_kokkos(m_rip_array, "MLIAP_SO3Kokkos:m_rip_array", totali);
    memoryKK->realloc_kokkos(m_rip_darray, "MLIAP_SO3Kokkos:m_rip_darray", totali);
    alloc_arrays += 2.0 * totali * sizeof(double);
  }

  totali = totaln * ncoefs * 3;
  if ( totali > (int)k_dplist_r.extent(0)) {
    memoryKK->realloc_kokkos(k_dplist_r, "MLIAP_SO3Kokkos:m_dplist_r", (int)totaln, ncoefs, 3);
    alloc_arrays += 2.0 * totali * sizeof(double);
  }

  t_numneighs = numneighs.template view<DeviceType>();
  t_jelems = jelems.template view<DeviceType>();
  t_wjelem = wjelem.template view<DeviceType>();
  t_rij = rij.template view<DeviceType>();
  t_ij = k_ij.template view<DeviceType>();
  t_nmax = nmax;
  t_lmax = lmax;
  t_rcut = rcut;
  t_alpha = alpha;

  {
    Kokkos::RangePolicy<DeviceType,MLIAPSO3GetSBESArrayTag> range(0,nlocal);
    Kokkos::parallel_for(range, *this);
  }

  {
    Kokkos::RangePolicy<DeviceType,MLIAPSO3GetRipArrayTag> range(0,nlocal);
    Kokkos::parallel_for(range, *this);
  }

  Kokkos::deep_copy(m_plist_r,0.);
  for (int start=0; start < nlocal; start += m_chunk_size) {
    int end=std::min(nlocal, start+m_chunk_size);
    Kokkos::deep_copy(m_clisttot_r, 0);
    Kokkos::deep_copy(m_clisttot_i, 0);
    {
      Kokkos::RangePolicy<DeviceType,MLIAPSO3SpectrumTag> range(start,end);
      Kokkos::parallel_for(range, *this);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::operator() (const MLIAP_SO3Kokkos<DeviceType>::MLIAPSO3SpectrumTag&, int ii) const {
  int ii_chunk = ii%m_chunk_size;
  auto ulist_r = Kokkos::subview(m_ulist_r,ii_chunk,Kokkos::ALL);
  auto ulist_i = Kokkos::subview(m_ulist_i,ii_chunk,Kokkos::ALL);
  auto clisttot_r=Kokkos::subview(m_clisttot_r,ii_chunk,Kokkos::ALL, Kokkos::ALL);
  auto clisttot_i=Kokkos::subview(m_clisttot_i,ii_chunk,Kokkos::ALL, Kokkos::ALL);

  const int twolmax = 2 * (t_lmax + 1);
  const int findex = m_nmax * (m_lmax + 1);
  int ipair = t_ij(ii);
  for (int neighbor = 0; neighbor < t_numneighs[ii]; neighbor++) {
    const int jelem = t_jelems[ipair];
    int weight = t_wjelem[jelem];

    double x = t_rij(ipair, 0);
    double y = t_rij(ipair, 1);
    double z = t_rij(ipair, 2);
    ipair++;

    double r = sqrt(x * x + y * y + z * z);

    if (r < SMALL) continue;
    for (int ti = 0; ti < m_idxu_count; ti++) {
      ulist_r[ti] = 0.0;
      ulist_i[ti] = 0.0;
    }

    compute_uarray_recursive(x, y, z, r, twolmax, ulist_r, ulist_i, m_idxu_block, m_rootpq);

    double sfac_weight = compute_sfac(r, t_rcut)*double(weight);

    int gindex = (ipair - 1) * findex;
    for (int n = 1; n < t_nmax + 1; n++) {
      int i = 0;
      for (int l = 0; l < t_lmax + 1; l++) {
        double r_int = m_rip_array[gindex + (bigint) (n - 1) * (m_lmax + 1) + l];

        for (int m = -l; m < l + 1; m++) {

          double Ylm_r = (ulist_r[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
          clisttot_r((n - 1), i) += r_int * Ylm_r * sfac_weight;
          double Ylm_i = (ulist_i[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
          clisttot_i((n - 1), i) += r_int * Ylm_i * sfac_weight;
          i += 1;
        }
      }
    }
  }
  compute_pi(t_nmax, t_lmax, clisttot_r, clisttot_i, m_numYlms, m_plist_r, ii);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAP_SO3Kokkos<DeviceType>::spectrum_dxdr(int nlocal, DAT::tdual_int_1d numneighs, DAT::tdual_int_1d jelems, DAT::tdual_float_1d wjelem,
                                    DAT::tdual_float_2d rij, DAT::tdual_int_1d k_ij,
                                    int nmax, int lmax, double rcut, double alpha, bigint totaln,
                                    int ncoefs)
{
  bigint totali;

  if ( nlocal > (int)m_clisttot_r.extent(0)){
    memoryKK->destroy_kokkos(m_clisttot_r);
    memoryKK->create_kokkos(m_clisttot_r, nlocal, m_nmax, m_numYlms, "MLIAP_SO3Kokkos:m_clisttot_r");
    memoryKK->destroy_kokkos(m_clisttot_i);
    memoryKK->create_kokkos(m_clisttot_i, nlocal, m_nmax, m_numYlms, "MLIAP_SO3Kokkos:m_clisttot_i");
    int num_of_temp = std::min(nlocal, m_chunk_size);
    int delta=num_of_temp-m_ulist_r.extent(0);
    if (delta > 0 ){
      memoryKK->destroy_kokkos(m_ulist_r);
      memoryKK->create_kokkos(m_ulist_r, num_of_temp, m_idxu_count, "MLIAP_SO3Kokkos:m_ulist_r");
      memoryKK->destroy_kokkos(m_ulist_i);
      memoryKK->create_kokkos(m_ulist_i, num_of_temp, m_idxu_count, "MLIAP_SO3Kokkos:m_ulist_i");
      alloc_arrays += 2.0 * m_idxu_count * delta * sizeof(double);
      memoryKK->destroy_kokkos(m_dYlm_r);
      memoryKK->create_kokkos(m_dYlm_r, num_of_temp, m_numYlms, 3, "MLIAP_SO3Kokkos:m_dYlm_r");
      memoryKK->destroy_kokkos(m_dYlm_i);
      memoryKK->create_kokkos(m_dYlm_i, num_of_temp, m_numYlms, 3, "MLIAP_SO3Kokkos:m_dYlm_i");
      alloc_arrays += 2.0 * m_numYlms * 3 * delta * sizeof(double);
    }
  }

  totali = totaln * m_Nmax * (m_lmax + 1);
  if ( totali > (int)m_sbes_array.extent(0)) {
    memoryKK->destroy_kokkos(m_sbes_array);
    memoryKK->create_kokkos(m_sbes_array, totali, "MLIAP_SO3Kokkos:m_sbes_array");
    memoryKK->destroy_kokkos(m_sbes_darray);
    memoryKK->create_kokkos(m_sbes_darray, totali, "MLIAP_SO3Kokkos:m_sbes_darray");

    totali = totaln * m_nmax * (m_lmax + 1);
    memoryKK->destroy_kokkos(m_rip_array);
    memoryKK->create_kokkos(m_rip_array, totali, "MLIAP_SO3Kokkos:m_rip_array");
    memoryKK->destroy_kokkos(m_rip_darray);
    memoryKK->create_kokkos(m_rip_darray, totali, "MLIAP_SO3Kokkos:m_rip_darray");

    memoryKK->destroy_kokkos(k_dplist_r);
    memoryKK->create_kokkos(k_dplist_r, (int)totaln, ncoefs, 3, "MLIAP_SO3Kokkos:m_dplist_r");
  }

  t_numneighs=numneighs.template view<DeviceType>();
  t_jelems=jelems.template view<DeviceType>();
  t_wjelem=wjelem.template view<DeviceType>();
  t_rij=rij.template view<DeviceType>();
  t_ij = k_ij.template view<DeviceType>();
  t_nmax=nmax;
  t_lmax=lmax;
  t_rcut=rcut;
  t_alpha=alpha;

  {
    Kokkos::RangePolicy<DeviceType,MLIAPSO3GetSBESArrayTag> range(0,nlocal);
    Kokkos::parallel_for(range, *this);
  }

  {
    Kokkos::RangePolicy<DeviceType,MLIAPSO3GetRipArrayTag> range(0,nlocal);
    Kokkos::parallel_for(range, *this);
  }

  Kokkos::deep_copy(k_dplist_r, 0.0);
  for (int start=0; start < nlocal; start += m_chunk_size) {
    int end=std::min(nlocal, start+m_chunk_size);
    Kokkos::deep_copy(m_clisttot_r, 0);
    Kokkos::deep_copy(m_clisttot_i, 0);
    {
      Kokkos::RangePolicy<DeviceType,MLIAPSO3SpectrumDXDRTag> range(start,end);
      Kokkos::parallel_for(range, *this);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void MLIAP_SO3Kokkos<DeviceType>::operator() (const MLIAP_SO3Kokkos<DeviceType>::MLIAPSO3SpectrumDXDRTag&, int ii) const {
  //TO-DO Need to move m_ulist_r, m_ulist_i into local shared memory
  int ii_chunk = ii%m_chunk_size;
  auto ulist_r = Kokkos::subview(m_ulist_r,ii_chunk,Kokkos::ALL);
  auto ulist_i = Kokkos::subview(m_ulist_i,ii_chunk,Kokkos::ALL);
  auto clisttot_r=Kokkos::subview(m_clisttot_r,ii_chunk,Kokkos::ALL, Kokkos::ALL);
  auto clisttot_i=Kokkos::subview(m_clisttot_i,ii_chunk,Kokkos::ALL, Kokkos::ALL);
  auto dYlm_r = Kokkos::subview(m_dYlm_r,ii_chunk,Kokkos::ALL,Kokkos::ALL);
  auto dYlm_i = Kokkos::subview(m_dYlm_i,ii_chunk,Kokkos::ALL,Kokkos::ALL);
  auto dclist = Kokkos::subview(m_dclist,ii_chunk,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
  int twolmax = 2 * (t_lmax + 1);
  int findex = m_nmax * (m_lmax + 1);
  int ipair = t_ij(ii);
  for (int neighbor = 0; neighbor < t_numneighs[ii]; neighbor++) {

    const int jelem = t_jelems[ipair];
    int weight = t_wjelem[jelem];

    double x = t_rij(ipair, 0);
    double y = t_rij(ipair, 1);
    double z = t_rij(ipair, 2);
    ipair++;

    double r = sqrt(x * x + y * y + z * z);

    if (r < SMALL) continue;

    for (bigint ti = 0; ti < m_idxu_count; ti++) {
      ulist_r[ti] = 0.0;
      ulist_i[ti] = 0.0;
    }

    compute_uarray_recursive(x, y, z, r, twolmax, ulist_r, ulist_i, m_idxu_block, m_rootpq);

    double sfac_weight = compute_sfac(r, t_rcut)*double(weight);
    bigint gindex = (ipair - 1) * findex;
    for (int n = 0; n < t_nmax; n++) {
      int i = 0;
      for (int l = 0; l < t_lmax + 1; l++) {
        double r_int = m_rip_array[gindex + n * (m_lmax + 1) + l];

        for (int m = -l; m < l + 1; m++) {
          double Ylm_r = (ulist_r[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
          clisttot_r(n, i) += r_int * Ylm_r * sfac_weight;
          double Ylm_i = (ulist_i[m_idxylm[i]]) * m_pfac[l * m_pfac_l2 + m];
          clisttot_i(n, i) += r_int * Ylm_i * sfac_weight;
          i += 1;
        }
      }
    }

  }

  ipair = t_ij(ii);
  for (int neighbor = 0; neighbor < t_numneighs[ii]; neighbor++) {
    const int jelem = t_jelems[ipair];
    int weight = t_wjelem[jelem];

    double x = t_rij(ipair, 0);
    double y = t_rij(ipair, 1);
    double z = t_rij(ipair, 2);
    ipair++;
    auto dplist_r=Kokkos::subview(k_dplist_r,ipair-1,Kokkos::ALL, Kokkos::ALL);

    double r = sqrt(x * x + y * y + z * z);
    if (r < SMALL) continue;

    for (int ti = 0; ti < m_idxu_count; ti++) {
      ulist_r[ti] = 0.0;
      ulist_i[ti] = 0.0;
    }

    compute_uarray_recursive(x, y, z, r, twolmax, ulist_r, ulist_i, m_idxu_block, m_rootpq);

    /////////  compute_carray_wD  ////////
    {
      double rvec[3];
      rvec[0] = x;
      rvec[1] = y;
      rvec[2] = z;

      for (int i=0;i<m_numYlms;++i)
        for (int j=0;j<3;++j)
          dYlm_r(i,j) = dYlm_i(i,j) = 0.0;

      double comj_i = 1.0 / sqrt(2.0);
      double oneofr = 1.0 / r;

      double dexpfac[3];
      for (int ii = 0; ii < 3; ii++) dexpfac[ii] = -2.0 * t_alpha * rvec[ii];

      double sfac = compute_sfac(r, t_rcut);
      double dsfac = compute_dsfac(r, t_rcut);
      double dsfac_arr[3];
      for (int ii = 0; ii < 3; ii++) { dsfac_arr[ii] = dsfac * rvec[ii] / r; }

      int i = 1;
      for (int l = 1; l < t_lmax + 1; l++) {
        for (int m = -l; m < l + 1; m++) {
          double dfact[6];
          double xcov0_r, xcov0_i, xcovpl1_r, xcovpl1_i, xcovm1_r, xcovm1_i;
          dfact[0] = m_dfac0[l * m_dfac_l2 + m] * oneofr;
          dfact[1] = m_dfac1[l * m_dfac_l2 + m] * oneofr;
          dfact[2] = m_dfac2[l * m_dfac_l2 + m] * oneofr;
          dfact[3] = m_dfac3[l * m_dfac_l2 + m] * oneofr;
          dfact[4] = m_dfac4[l * m_dfac_l2 + m] * oneofr;
          dfact[5] = m_dfac5[l * m_dfac_l2 + m] * oneofr;

          int idx=m_ellpl1[l] + m;
          xcov0_r = dfact[0] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
          xcov0_i = dfact[0] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          if (abs(m) <= l - 1.0) {
            idx=m_ellm1[l] + m;
            xcov0_r += dfact[1] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
            xcov0_i += dfact[1] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          }
          idx=m_ellpl1[l] + m + 1;
          xcovpl1_r = dfact[2] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
          xcovpl1_i = dfact[2] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          if (abs(m + 1) <= l - 1.0) {
            idx=m_ellm1[l] + m + 1;
            xcovpl1_r -= dfact[3] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
            xcovpl1_i -= dfact[3] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          }
          idx=m_ellpl1[l] + m - 1;
          xcovm1_r = dfact[4] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
          xcovm1_i = dfact[4] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          if (fabs(m - 1.0) <= l - 1.0) {
            idx=m_ellm1[l] + m - 1;
            xcovm1_r -= dfact[5] * m_Ylms[idx]*(ulist_r[m_idxylm[idx]]);
            xcovm1_i -= dfact[5] * m_Ylms[idx]*(ulist_i[m_idxylm[idx]]);
          }
          dYlm_r(i, 0) = 1.0 / sqrt(2.0) * (xcovm1_r - xcovpl1_r);
          dYlm_r(i, 1) = -comj_i * (xcovm1_i + xcovpl1_i);
          dYlm_r(i, 2) = xcov0_r;

          dYlm_i(i, 0) = 1.0 / sqrt(2.0) * (xcovm1_i - xcovpl1_i);
          dYlm_i(i, 1) = comj_i * (xcovm1_r + xcovpl1_r);
          dYlm_i(i, 2) = xcov0_i;

          i += 1;
        }
      }

      // Break up the loops into real and imaginary,
      //Real loop
      int dp_indx = 0;
      for (int n = 0; n < t_nmax; n++) {
        int i = 0;
        for (int l = 0; l < t_lmax + 1; l++) {
          double r_int = m_rip_array[(ipair - 1) * m_nmax * (m_lmax + 1) +
                             (bigint) n * (m_lmax + 1) + l];
          double r_int_temp = m_rip_darray[(ipair - 1) * m_nmax * (m_lmax + 1) +
                                    (bigint) n * (m_lmax + 1) + l];
          double dr_int[3];
          for (int ii = 0; ii < 3; ii++) dr_int[ii] = r_int_temp * 2.0 * t_alpha * rvec[ii] / r;

          for (int m = -l; m < l + 1; m++) {
            double Ylms_r = m_Ylms[i]*(ulist_r[m_idxylm[i]]);
            dclist(n, i,  0) =
                (r_int * Ylms_r * dexpfac[0] + dr_int[0] * Ylms_r +
                 r_int * dYlm_r(i, 0)) *
                sfac;
            dclist(n, i,  1) =
                (r_int * Ylms_r * dexpfac[1] + dr_int[1] * Ylms_r +
                 r_int * dYlm_r(i, 1)) *
                sfac;
            dclist(n, i,  2) =
                (r_int * Ylms_r * dexpfac[2] + dr_int[2] * Ylms_r +
                 r_int * dYlm_r(i, 2)) *
                sfac;

            dclist(n, i,  0) += (r_int * Ylms_r) * dsfac_arr[0];
            dclist(n, i,  1) += (r_int * Ylms_r) * dsfac_arr[1];
            dclist(n, i,  2) += (r_int * Ylms_r) * dsfac_arr[2];

            for (int k=0;k<3;++k){
              dclist(n,i,k) *= double(weight);
            }
            i += 1;
          }
        }

        for (int n2 = 0; n2 < n + 1; n2++) {
          int j = 0;
          for (int l = 0; l < t_lmax + 1; l++) {
            double norm = 2.0 * sqrt(2.0) * MY_PI / sqrt(2.0 * l + 1.0);
            for (int m = -l; m < l + 1; m++) {
              for (int idim = 0; idim < 3; idim++) {
                double temp_r;
                temp_r = dclist(n,j,idim) * clisttot_r(n2, j);
                temp_r += clisttot_r(n, j) * dclist(n2,j,idim);
                temp_r *= norm;
                dplist_r(dp_indx, idim) += temp_r;
              }
              j += 1;
            }
            dp_indx += 1;
          }
        }

      }

      // Imaginary loop
      dp_indx = 0;
      for (int n = 0; n < t_nmax; n++) {
        int i = 0;
        for (int l = 0; l < t_lmax + 1; l++) {
          double r_int = m_rip_array[(ipair - 1) * m_nmax * (m_lmax + 1) +
                             (bigint) n * (m_lmax + 1) + l];
          double r_int_temp = m_rip_darray[(ipair - 1) * m_nmax * (m_lmax + 1) +
                                    (bigint) n * (m_lmax + 1) + l];
          double dr_int[3];
          for (int ii = 0; ii < 3; ii++) dr_int[ii] = r_int_temp * 2.0 * t_alpha * rvec[ii] / r;

          for (int m = -l; m < l + 1; m++) {
            double Ylms_i = m_Ylms[i]*(ulist_i[m_idxylm[i]]);
            dclist(n, i,  0) =
                (r_int * Ylms_i * dexpfac[0] + dr_int[0] * Ylms_i +
                 r_int * dYlm_i(i, 0)) *
                sfac;
            dclist(n, i,  1) =
                (r_int * Ylms_i * dexpfac[1] + dr_int[1] * Ylms_i +
                 r_int * dYlm_i(i, 1)) *
                sfac;
            dclist(n, i,  2) =
                (r_int * Ylms_i * dexpfac[2] + dr_int[2] * Ylms_i +
                 r_int * dYlm_i(i, 2)) *
                sfac;

            dclist(n, i,  0) += (r_int * Ylms_i) * dsfac_arr[0];
            dclist(n, i,  1) += (r_int * Ylms_i) * dsfac_arr[1];
            dclist(n, i,  2) += (r_int * Ylms_i) * dsfac_arr[2];

            for (int k=0;k<3;++k){
              dclist(n,i,k) *= double(weight);
            }
            i += 1;
          }
        }
        for (int n2 = 0; n2 < n + 1; n2++) {
          int j = 0;
          for (int l = 0; l < t_lmax + 1; l++) {
            double norm = 2.0 * sqrt(2.0) * MY_PI / sqrt(2.0 * l + 1.0);
            for (int m = -l; m < l + 1; m++) {
              for (int idim = 0; idim < 3; idim++) {
                double temp_r;
                temp_r = dclist(n,j,idim) * clisttot_i(n2, j);
                temp_r += clisttot_i(n, j) * dclist(n2,j,idim);
                temp_r *= norm;
                dplist_r(dp_indx, idim) += temp_r;
              }
              j += 1;
            }
            dp_indx += 1;
          }
        }

      }
    } /////// end compute_carray_wD //////////////////
  }    //for(neighbor=0;neighbor<numneighs[ii];neighbor++){
}

namespace LAMMPS_NS {
template class MLIAP_SO3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MLIAP_SO3Kokkos<LMPHostType>;
#endif
}
