// clang-format off
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
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "fix_eos_table_rx_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "comm.h"
#include <cmath>
#include "atom_masks.h"

#define MAXLINE 1024

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEOStableRXKokkos<DeviceType>::FixEOStableRXKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEOStableRX(lmp, narg, arg)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  update_table = 1;
  h_table = new TableHost();
  d_table = new TableDevice();

  k_error_flag = DAT::tdual_int_scalar("fix:error_flag");
  k_warning_flag = DAT::tdual_int_scalar("fix:warning_flag");

  k_dHf = DAT::tdual_float_1d("fix:dHf",nspecies);
  k_energyCorr = DAT::tdual_float_1d("fix:energyCorr",nspecies);
  k_tempCorrCoeff = DAT::tdual_float_1d("fix:tempCorrCoeff",nspecies);
  k_moleculeCorrCoeff = DAT::tdual_float_1d("fix:moleculeCorrCoeff",nspecies);
  for (int n = 0; n < nspecies; n++) {
    k_dHf.h_view(n) = dHf[n];
    k_energyCorr.h_view(n) = energyCorr[n];
    k_tempCorrCoeff.h_view(n) = tempCorrCoeff[n];
    k_moleculeCorrCoeff.h_view(n) = moleculeCorrCoeff[n];
  }

  k_dHf.modify<LMPHostType>();
  k_dHf.sync<DeviceType>();
  d_dHf = k_dHf.view<DeviceType>();

  k_energyCorr.modify<LMPHostType>();
  k_energyCorr.sync<DeviceType>();
  d_energyCorr = k_energyCorr.view<DeviceType>();

  k_tempCorrCoeff.modify<LMPHostType>();
  k_tempCorrCoeff.sync<DeviceType>();
  d_tempCorrCoeff = k_tempCorrCoeff.view<DeviceType>();

  k_moleculeCorrCoeff.modify<LMPHostType>();
  k_moleculeCorrCoeff.sync<DeviceType>();
  d_moleculeCorrCoeff = k_moleculeCorrCoeff.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEOStableRXKokkos<DeviceType>::~FixEOStableRXKokkos()
{
  if (copymode) return;

  delete h_table;
  delete d_table;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::setup(int /*vflag*/)
{
  if (update_table)
    create_kokkos_tables();

  copymode = 1;

  int nlocal = atom->nlocal;
  boltz = force->boltz;
  mask = atomKK->k_mask.view<DeviceType>();
  uCond = atomKK->k_uCond.view<DeviceType>();
  uMech = atomKK->k_uMech.view<DeviceType>();
  uChem = atomKK->k_uChem.view<DeviceType>();
  dpdTheta= atomKK->k_dpdTheta.view<DeviceType>();
  uCG = atomKK->k_uCG.view<DeviceType>();
  uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();

  if (!this->restart_reset) {
    atomKK->sync(execution_space,MASK_MASK | UCHEM_MASK | UCG_MASK | UCGNEW_MASK);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXSetup>(0,nlocal),*this);
    atomKK->modified(execution_space,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);
  }

  // Communicate the updated momenta and velocities to all nodes
  atomKK->sync(Host,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);
  comm->forward_comm_fix(this);
  atomKK->modified(Host,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);

  atomKK->sync(execution_space,MASK_MASK | UCOND_MASK | UMECH_MASK | UCHEM_MASK | DPDTHETA_MASK | DVECTOR_MASK);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXTemperatureLookup>(0,nlocal),*this);
  atomKK->modified(execution_space,DPDTHETA_MASK);

  error_check();

  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::operator()(TagFixEOStableRXSetup, const int &i) const {
  if (mask[i] & groupbit) {
    const double duChem = uCG[i] - uCGnew[i];
    uChem[i] += duChem;
    uCG[i] = 0.0;
    uCGnew[i] = 0.0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::operator()(TagFixEOStableRXTemperatureLookup, const int &i) const {
  if (mask[i] & groupbit)
    temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::init()
{
  if (update_table)
    create_kokkos_tables();

  copymode = 1;

  int nlocal = atom->nlocal;
  boltz = force->boltz;
  mask = atomKK->k_mask.view<DeviceType>();
  uCond = atomKK->k_uCond.view<DeviceType>();
  uMech = atomKK->k_uMech.view<DeviceType>();
  uChem = atomKK->k_uChem.view<DeviceType>();
  dpdTheta= atomKK->k_dpdTheta.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();

  if (this->restart_reset) {
    atomKK->sync(execution_space,MASK_MASK | UCOND_MASK | UMECH_MASK | UCHEM_MASK | DPDTHETA_MASK | DVECTOR_MASK);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXTemperatureLookup>(0,nlocal),*this);
    atomKK->modified(execution_space,DPDTHETA_MASK);
  } else {
    atomKK->sync(execution_space,MASK_MASK | UCOND_MASK | UMECH_MASK | UCHEM_MASK | DPDTHETA_MASK | DVECTOR_MASK);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXInit>(0,nlocal),*this);
    atomKK->modified(execution_space,UCOND_MASK | UMECH_MASK | UCHEM_MASK);
  }

  error_check();

  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::operator()(TagFixEOStableRXInit, const int &i) const {
  double tmp;
  if (mask[i] & groupbit) {
    if (dpdTheta[i] <= 0.0)
      k_error_flag.template view<DeviceType>()() = 1;
    energy_lookup(i,dpdTheta[i],tmp);
    uCond[i] = 0.0;
    uMech[i] = tmp;
    uChem[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::post_integrate()
{
  if (update_table)
    create_kokkos_tables();

  copymode = 1;

  int nlocal = atom->nlocal;
  boltz = force->boltz;
  mask = atomKK->k_mask.view<DeviceType>();
  uCond = atomKK->k_uCond.view<DeviceType>();
  uMech = atomKK->k_uMech.view<DeviceType>();
  uChem = atomKK->k_uChem.view<DeviceType>();
  dpdTheta= atomKK->k_dpdTheta.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();

  atomKK->sync(execution_space,MASK_MASK | UCOND_MASK | UMECH_MASK | UCHEM_MASK | DPDTHETA_MASK | DVECTOR_MASK);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXTemperatureLookup2>(0,nlocal),*this);
  atomKK->modified(execution_space,DPDTHETA_MASK);

  error_check();

  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::operator()(TagFixEOStableRXTemperatureLookup2, const int &i) const {
  if (mask[i] & groupbit) {
    temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
    if (dpdTheta[i] <= 0.0)
      k_error_flag.template view<DeviceType>()() = 1;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::end_of_step()
{
  if (update_table)
    create_kokkos_tables();

  copymode = 1;

  int nlocal = atom->nlocal;
  boltz = force->boltz;
  mask = atomKK->k_mask.view<DeviceType>();
  uCond = atomKK->k_uCond.view<DeviceType>();
  uMech = atomKK->k_uMech.view<DeviceType>();
  uChem = atomKK->k_uChem.view<DeviceType>();
  dpdTheta= atomKK->k_dpdTheta.view<DeviceType>();
  uCG = atomKK->k_uCG.view<DeviceType>();
  uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();


  // Communicate the ghost uCGnew
  atomKK->sync(Host,UCG_MASK | UCGNEW_MASK);
  comm->reverse_comm_fix(this);
  atomKK->modified(Host,UCG_MASK | UCGNEW_MASK);

  atomKK->sync(execution_space,MASK_MASK | UCHEM_MASK | UCG_MASK | UCGNEW_MASK);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXSetup>(0,nlocal),*this);
  atomKK->modified(execution_space,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);

  // Communicate the updated momenta and velocities to all nodes
  atomKK->sync(Host,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);
  comm->forward_comm_fix(this);
  atomKK->modified(Host,UCHEM_MASK | UCG_MASK | UCGNEW_MASK);

  atomKK->sync(execution_space,MASK_MASK | UCOND_MASK | UMECH_MASK | UCHEM_MASK | DPDTHETA_MASK | DVECTOR_MASK);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEOStableRXTemperatureLookup2>(0,nlocal),*this);
  atomKK->modified(execution_space,DPDTHETA_MASK);

  error_check();

  copymode = 0;
}

/* ----------------------------------------------------------------------
   calculate potential ui at temperature thetai
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::energy_lookup(int id, double thetai, double &ui) const
{
  int itable, nPG;
  double fraction, uTmp, nMolecules, nTotal, nTotalPG;
  double tolerance = 1.0e-10;

  ui = 0.0;
  nTotal = 0.0;
  nTotalPG = 0.0;
  nPG = 0;

  if (rx_flag) {
    for (int ispecies = 0; ispecies < nspecies; ispecies++) {
      nTotal += dvector(ispecies,id);
      if (fabs(d_moleculeCorrCoeff[ispecies]) > tolerance) {
        nPG++;
        nTotalPG += dvector(ispecies,id);
      }
    }
  } else {
    nTotal = 1.0;
  }

  for (int ispecies=0;ispecies<nspecies;ispecies++) {
    //Table *tb = &tables[ispecies];
    //thetai = MAX(thetai,tb->lo);
    thetai = MAX(thetai,d_table_const.lo(ispecies));
    //thetai = MIN(thetai,tb->hi);
    thetai = MIN(thetai,d_table_const.hi(ispecies));

    if (tabstyle == LINEAR) {
      //itable = static_cast<int> ((thetai - tb->lo) * tb->invdelta);
      itable = static_cast<int> ((thetai - d_table_const.lo(ispecies)) * d_table_const.invdelta(ispecies));
      //fraction = (thetai - tb->r[itable]) * tb->invdelta;
      fraction = (thetai - d_table_const.r(ispecies,itable)) * d_table_const.invdelta(ispecies);
      //uTmp = tb->e[itable] + fraction*tb->de[itable];
      uTmp = d_table_const.e(ispecies,itable) + fraction*d_table_const.de(ispecies,itable);

      uTmp += d_dHf[ispecies];
      uTmp += d_tempCorrCoeff[ispecies]*thetai; // temperature correction
      uTmp += d_energyCorr[ispecies]; // energy correction
      if (nPG > 0) ui += d_moleculeCorrCoeff[ispecies]*nTotalPG/double(nPG); // molecule correction

      if (rx_flag) nMolecules = dvector(ispecies,id);
      else nMolecules = 1.0;
      ui += nMolecules*uTmp;
    }
  }
  ui = ui - double(nTotal+1.5)*boltz*thetai;
}

/* ----------------------------------------------------------------------
   calculate temperature thetai at energy ui
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEOStableRXKokkos<DeviceType>::temperature_lookup(int id, double ui, double &thetai) const
{
  //Table *tb = &tables[0];

  int it;
  double t1,t2,u1,u2,f1,f2;
  double maxit = 100;
  double temp;
  double delta = 0.001;
  double tolerance = 1.0e-10;
  int lo = d_table_const.lo(0);
  int hi = d_table_const.hi(0);

  // Store the current thetai in t1
  t1 = MAX(thetai,lo);
  t1 = MIN(t1,hi);
  if (t1==hi) delta = -delta;

  // Compute u1 at thetai
  energy_lookup(id,t1,u1);

  // Compute f1
  f1 = u1 - ui;

  // Compute guess of t2
  t2 = (1.0 + delta)*t1;

  // Compute u2 at t2
  energy_lookup(id,t2,u2);

  // Compute f1
  f2 = u2 - ui;

  // Apply the Secant Method
  for (it=0; it<maxit; it++) {
    if (fabs(f2-f1) < MY_EPSILON) {
      if (std::isnan(f1) || std::isnan(f2)) k_error_flag.template view<DeviceType>()() = 2;
      temp = t1;
      temp = MAX(temp,lo);
      temp = MIN(temp,hi);
      k_warning_flag.template view<DeviceType>()() = 1;
      break;
    }
    temp = t2 - f2*(t2-t1)/(f2-f1);
    if (fabs(temp-t2) < tolerance) break;
    f1 = f2;
    t1 = t2;
    t2 = temp;
    energy_lookup(id,t2,u2);
    f2 = u2 - ui;
  }
  if (it==maxit) {
    if (std::isnan(f1) || std::isnan(f2) || std::isnan(ui) || std::isnan(thetai) || std::isnan(t1) || std::isnan(t2))
      k_error_flag.template view<DeviceType>()() = 2;
    else
      k_error_flag.template view<DeviceType>()() = 3;
  }
  thetai = temp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixEOStableRXKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int ii,jj,m;
  HAT::t_efloat_1d h_uChem = atomKK->k_uChem.h_view;
  HAT::t_efloat_1d h_uCG = atomKK->k_uCG.h_view;
  HAT::t_efloat_1d h_uCGnew = atomKK->k_uCGnew.h_view;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    buf[m++] = h_uChem[jj];
    buf[m++] = h_uCG[jj];
    buf[m++] = h_uCGnew[jj];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;
  HAT::t_efloat_1d h_uChem = atomKK->k_uChem.h_view;
  HAT::t_efloat_1d h_uCG = atomKK->k_uCG.h_view;
  HAT::t_efloat_1d h_uCGnew = atomKK->k_uCGnew.h_view;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++) {
    h_uChem[ii]  = buf[m++];
    h_uCG[ii]    = buf[m++];
    h_uCGnew[ii] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixEOStableRXKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  HAT::t_efloat_1d h_uCG = atomKK->k_uCG.h_view;
  HAT::t_efloat_1d h_uCGnew = atomKK->k_uCGnew.h_view;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = h_uCG[i];
    buf[m++] = h_uCGnew[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  HAT::t_efloat_1d h_uCG = atomKK->k_uCG.h_view;
  HAT::t_efloat_1d h_uCGnew = atomKK->k_uCGnew.h_view;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    h_uCG[j] += buf[m++];
    h_uCGnew[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::error_check()
{
  k_error_flag.template modify<DeviceType>();
  k_error_flag.template sync<LMPHostType>();
  if (k_error_flag.h_view() == 1)
    error->one(FLERR,"Internal temperature <= zero");
  else if (k_error_flag.h_view() == 2)
    error->one(FLERR,"NaN detected in secant solver.");
  else if (k_error_flag.h_view() == 3)
    error->one(FLERR,"Maxit exceeded in secant solver.");

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (k_warning_flag.h_view()) {
    error->warning(FLERR,"Secant solver did not converge because table bounds were exceeded.");
    k_warning_flag.h_view() = 0;
    k_warning_flag.template modify<LMPHostType>();
    k_warning_flag.template sync<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEOStableRXKokkos<DeviceType>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(d_table->lo,h_table->lo,ntables,"Table::lo");
  memoryKK->create_kokkos(d_table->hi,h_table->hi,ntables,"Table::hi");
  memoryKK->create_kokkos(d_table->invdelta,h_table->invdelta,ntables,"Table::invdelta");

  if (tabstyle == LINEAR) {
    memoryKK->create_kokkos(d_table->r,h_table->r,ntables,tablength,"Table::r");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,tlm1,"Table::de");
  }

  for (int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->lo[i] = tb->lo;
    h_table->hi[i] = tb->hi;
    h_table->invdelta[i] = tb->invdelta;

    for (int j = 0; j < (int)h_table->r.extent(1); j++)
      h_table->r(i,j) = tb->r[j];
    for (int j = 0; j < (int)h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for (int j = 0; j < (int)h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
  }

  Kokkos::deep_copy(d_table->lo,h_table->lo);
  Kokkos::deep_copy(d_table->hi,h_table->hi);
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  Kokkos::deep_copy(d_table->r,h_table->r);
  Kokkos::deep_copy(d_table->e,h_table->e);
  Kokkos::deep_copy(d_table->de,h_table->de);

  d_table_const.lo = d_table->lo;
  d_table_const.hi = d_table->hi;
  d_table_const.invdelta = d_table->invdelta;
  d_table_const.r = d_table->r;
  d_table_const.e = d_table->e;
  d_table_const.de = d_table->de;

  update_table = 0;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixEOStableRXKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixEOStableRXKokkos<LMPHostType>;
#endif
}
