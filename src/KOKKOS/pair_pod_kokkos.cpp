// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   aE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ngoc Cuong Nguyen (MIT)
------------------------------------------------------------------------- */

#include "pair_pod_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"

#include <cstring>
#include <chrono>

#include "eapod.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using MathSpecial::powint;

enum{FS,FS_SHIFTEDSCALED};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairPODKokkos<DeviceType>::PairPODKokkos(LAMMPS *lmp) : PairPOD(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  ni = 0;
  nimax = 0;
  nij = 0;
  nijmax = 0;
  atomBlockSize = 2048;
  nAtomBlocks = 0;
  timing = 0;
  for (int i=0; i<100; i++) comptime[i] = 0;

  host_flag = (execution_space == Host);
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairPODKokkos<DeviceType>::~PairPODKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPODKokkos<DeviceType>::init_style()
{
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR,"Pair style pod/kk can currently only run on a single "
                         "CPU thread");

    PairPOD::init_style();
    return;
  }

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style POD requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style POD requires newton pair on");

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair pace/kk");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPODKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairPOD::init_one(i,j);

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPODKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "pair_coeff", error);

  PairPOD::coeff(narg,arg); // create a PairPOD object

  copy_from_pod_class(PairPOD::fastpodptr); // copy parameters and arrays from pod class

  int n = atom->ntypes + 1;
  MemKK::realloc_kokkos(d_map, "pair_pod:map", n);

  MemKK::realloc_kokkos(k_cutsq, "pair_pod:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();

  MemKK::realloc_kokkos(k_scale, "pair_pod:scale", n, n);
  d_scale = k_scale.template view<DeviceType>();

  // Set up element lists

  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int i = 1; i <= atom->ntypes; i++)
    h_map(i) = map[i];

  Kokkos::deep_copy(d_map,h_map);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPODKokkos<DeviceType>::allocate()
{
  PairPOD::allocate();
}

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& max_neighs) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPODKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  int newton_pair = force->newton_pair;
  if (newton_pair == false)
    error->all(FLERR,"PairPODKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  maxneigh = 0;
  if (host_flag) {
    inum = list->inum;
    d_numneigh = typename ArrayTypes<DeviceType>::t_int_1d("pair_pod:numneigh",inum);
    for (int i=0; i<inum; i++) d_numneigh(i) = list->numneigh[i];
    d_ilist = typename ArrayTypes<DeviceType>::t_int_1d("pair_pod:ilist",inum);
    for (int i=0; i<inum; i++) d_ilist(i) = list->ilist[i];

    int maxn = 0;
    for (int i=0; i<inum; i++)
      if (maxn < list->numneigh[i]) maxn = list->numneigh[i];
    MemoryKokkos::realloc_kokkos(d_neighbors,"neighlist:neighbors",inum,maxn);
    for (int i=0; i<inum; i++) {
      int gi = list->ilist[i];
      int m = list->numneigh[gi];
      if (maxneigh<m) maxneigh = m;
      for (int l = 0; l < m; l++) {           // loop over each atom around atom i
        d_neighbors(gi, l) = list->firstneigh[gi][l];
      }
    }
  }
  else {
    NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
    d_numneigh = k_list->d_numneigh;
    d_neighbors = k_list->d_neighbors;
    d_ilist = k_list->d_ilist;
    inum = list->inum;
    int maxneighs;
    Kokkos::parallel_reduce("PairPODKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneighs));
    maxneigh = maxneighs;
  }

  auto begin = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();

  // determine the number of atom blocks and divide atoms into blocks
  nAtomBlocks = calculateNumberOfIntervals(inum, atomBlockSize);
  if (nAtomBlocks > 100) nAtomBlocks = 100;
  divideInterval(atomBlocks, inum, nAtomBlocks);

  int nmax = 0;
  for (int block=0; block<nAtomBlocks; block++) {
    int n = atomBlocks[block+1] - atomBlocks[block];
    if (nmax < n) nmax = n;
  }
  grow_atoms(nmax);
  grow_pairs(nmax*maxneigh);

  rcutsq = rcut*rcut;
  for (int block=0; block<nAtomBlocks; block++) {
    int gi1 = atomBlocks[block]-1;
    int gi2 = atomBlocks[block+1]-1;
    ni = gi2 - gi1; // total number of atoms in the current atom block

    begin = std::chrono::high_resolution_clock::now();
    // calculate the total number of pairs (i,j) in the current atom block
    nij = NeighborCount(numij, rcutsq, gi1, ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[0] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    begin = std::chrono::high_resolution_clock::now();
    // obtain the neighbors within rcut
    NeighborList(rij, numij, typeai, idxi, ai, aj, ti, tj, rcutsq, gi1, ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[1] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    // compute atomic energy and force for the current atom block
    begin = std::chrono::high_resolution_clock::now();
    blockatom_energyforce(ei, fij, ni, nij);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[2] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    begin = std::chrono::high_resolution_clock::now();
    // tally atomic energy to global energy
    tallyenergy(ei, gi1, ni);

    // tally atomic force to global force
    tallyforce(fij, ai, aj, nij);

    // tally atomic stress
    if (vflag) {
      tallystress(fij, rij, ai, aj, nij);
    }
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[3] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    //savedatafordebugging();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::copy_from_pod_class(EAPOD *podptr)
{
  nelements = podptr->nelements; // number of elements
  onebody = podptr->onebody;   // one-body descriptors
  besseldegree = podptr->besseldegree; // degree of Bessel functions
  inversedegree = podptr->inversedegree; // degree of inverse functions
  nbesselpars = podptr->nbesselpars;  // number of Bessel parameters
  nCoeffPerElement = podptr->nCoeffPerElement; // number of coefficients per element = (nl1 + Mdesc*nClusters)
  ns = podptr->ns;      // number of snapshots for radial basis functions
  nl1 = podptr->nl1;  // number of one-body descriptors
  nl2 = podptr->nl2;  // number of two-body descriptors
  nl3 = podptr->nl3;  // number of three-body descriptors
  nl4 = podptr->nl4;  // number of four-body descriptors
  nl23 = podptr->nl23; // number of two-body x three-body descriptors
  nl33 = podptr->nl33; // number of three-body x three-body descriptors
  nl34 = podptr->nl34; // number of three-body x four-body descriptors
  nl44 = podptr->nl44; // number of four-body x four-body descriptors
  nl = podptr->nl;   // number of local descriptors
  nrbf2 = podptr->nrbf2;
  nrbf3 = podptr->nrbf3;
  nrbf4 = podptr->nrbf4;
  nrbfmax = podptr->nrbfmax; // number of radial basis functions
  nabf3 = podptr->nabf3;     // number of three-body angular basis functions
  nabf4 = podptr->nabf4;     // number of four-body angular basis functions
  K3 = podptr->K3;           // number of three-body monomials
  K4 = podptr->K4;           // number of four-body monomials
  Q4 = podptr->Q4;           // number of four-body monomial coefficients
  nClusters = podptr->nClusters; // number of environment clusters
  nComponents = podptr->nComponents; // number of principal components
  Mdesc = podptr->Mdesc; // number of base descriptors

  rin = podptr->rin;
  rcut = podptr->rcut;
  rmax = rcut - rin;

  MemKK::realloc_kokkos(besselparams, "pair_pod:besselparams", 3);
  auto h_besselparams = Kokkos::create_mirror_view(besselparams);
  h_besselparams[0] = podptr->besselparams[0];
  h_besselparams[1] = podptr->besselparams[1];
  h_besselparams[2] = podptr->besselparams[2];
  Kokkos::deep_copy(besselparams, h_besselparams);

  MemKK::realloc_kokkos(elemindex, "pair_pod:elemindex", nelements*nelements);
  auto h_elemindex = Kokkos::create_mirror_view(elemindex);
  for (int i=0; i<nelements*nelements; i++) h_elemindex[i] = podptr->elemindex[i];
  Kokkos::deep_copy(elemindex, h_elemindex);

  MemKK::realloc_kokkos(Phi, "pair_pod:Phi", ns*ns);
  auto h_Phi = Kokkos::create_mirror_view(Phi);
  for (int i=0; i<ns*ns; i++) h_Phi[i] = podptr->Phi[i];
  Kokkos::deep_copy(Phi, h_Phi);

  MemKK::realloc_kokkos(coefficients, "pair_pod:coefficients", nCoeffPerElement * nelements);
  auto h_coefficients = Kokkos::create_mirror_view(coefficients);
  for (int i=0; i<nCoeffPerElement * nelements; i++) h_coefficients[i] = podptr->coeff[i];
  Kokkos::deep_copy(coefficients, h_coefficients);

  if (nClusters > 1) {
    MemKK::realloc_kokkos(Proj, "pair_pod:Proj",  Mdesc * nComponents * nelements);
    auto h_Proj = Kokkos::create_mirror_view(Proj);
    for (int i=0; i<Mdesc * nComponents * nelements; i++) h_Proj[i] = podptr->Proj[i];
    Kokkos::deep_copy(Proj, h_Proj);

    MemKK::realloc_kokkos(Centroids, "pair_pod:Centroids",  nClusters * nComponents * nelements);
    auto h_Centroids = Kokkos::create_mirror_view(Centroids);
    for (int i=0; i<nClusters * nComponents * nelements; i++) h_Centroids[i] = podptr->Centroids[i];
    Kokkos::deep_copy(Centroids, h_Centroids);
  }

  MemKK::realloc_kokkos(pn3, "pair_pod:pn3", nabf3+1); // array stores the number of monomials for each degree
  MemKK::realloc_kokkos(pq3, "pair_pod:pq3", K3*2); // array needed for the recursive computation of the angular basis functions
  MemKK::realloc_kokkos(pc3, "pair_pod:pc3", K3);   // array needed for the computation of the three-body descriptors
  MemKK::realloc_kokkos(pa4, "pair_pod:pa4", nabf4+1); // this array is a subset of the array {0, 1, 4, 10, 19, 29, 47, 74, 89, 119, 155, 209, 230, 275, 335, 425, 533, 561, 624, 714, 849, 949, 1129, 1345}
  MemKK::realloc_kokkos(pb4, "pair_pod:pb4", Q4*3); // array stores the indices of the monomials needed for the computation of the angular basis functions
  MemKK::realloc_kokkos(pc4, "pair_pod:pc4", Q4);   // array of monomial coefficients needed for the computation of the four-body descriptors

  auto h_pn3 = Kokkos::create_mirror_view(pn3);
  for (int i=0; i<nabf3+1; i++) h_pn3[i] = podptr->pn3[i];
  Kokkos::deep_copy(pn3, h_pn3);

  auto h_pq3 = Kokkos::create_mirror_view(pq3);
  for (int i = 0; i < K3*2; i++) h_pq3[i] = podptr->pq3[i];
  Kokkos::deep_copy(pq3, h_pq3);

  auto h_pc3 = Kokkos::create_mirror_view(pc3);
  for (int i = 0; i < K3; i++) h_pc3[i] = podptr->pc3[i];
  Kokkos::deep_copy(pc3, h_pc3);

  auto h_pa4 = Kokkos::create_mirror_view(pa4);
  for (int i = 0; i < nabf4+1; i++) h_pa4[i] = podptr->pa4[i];
  Kokkos::deep_copy(pa4, h_pa4);

  auto h_pb4 = Kokkos::create_mirror_view(pb4);
  for (int i = 0; i < Q4*3; i++) h_pb4[i] = podptr->pb4[i];
  Kokkos::deep_copy(pb4, h_pb4);

  auto h_pc4 = Kokkos::create_mirror_view(pc4);
  for (int i = 0; i < Q4; i++) h_pc4[i] = podptr->pc4[i];
  Kokkos::deep_copy(pc4, h_pc4);

  MemKK::realloc_kokkos(ind33l, "pair_pod:ind33l", nl33);
  MemKK::realloc_kokkos(ind33r, "pair_pod:ind33r", nl33);
  MemKK::realloc_kokkos(ind34l, "pair_pod:ind34l", nl34);
  MemKK::realloc_kokkos(ind34r, "pair_pod:ind34r", nl34);
  MemKK::realloc_kokkos(ind44l, "pair_pod:ind44l", nl44);
  MemKK::realloc_kokkos(ind44r, "pair_pod:ind44r", nl44);

  auto h_ind33l = Kokkos::create_mirror_view(ind33l);
  for (int i = 0; i < nl33; i++) h_ind33l[i] = podptr->ind33l[i];
  Kokkos::deep_copy(ind33l, h_ind33l);

  auto h_ind33r = Kokkos::create_mirror_view(ind33r);
  for (int i = 0; i < nl33; i++) h_ind33r[i] = podptr->ind33r[i];
  Kokkos::deep_copy(ind33r, h_ind33r);

  auto h_ind34l = Kokkos::create_mirror_view(ind34l);
  for (int i = 0; i < nl34; i++) h_ind34l[i] = podptr->ind34l[i];
  Kokkos::deep_copy(ind34l, h_ind34l);

  auto h_ind34r = Kokkos::create_mirror_view(ind34r);
  for (int i = 0; i < nl34; i++) h_ind34r[i] = podptr->ind34r[i];
  Kokkos::deep_copy(ind34r, h_ind34r);

  auto h_ind44l = Kokkos::create_mirror_view(ind44l);
  for (int i = 0; i < nl44; i++) h_ind44l[i] = podptr->ind44l[i];
  Kokkos::deep_copy(ind44l, h_ind44l);

  auto h_ind44r = Kokkos::create_mirror_view(ind44r);
  for (int i = 0; i < nl44; i++) h_ind44r[i] = podptr->ind44r[i];
  Kokkos::deep_copy(ind44r, h_ind44r);
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::divideInterval(int *intervals, int N, int M)
{
  int intervalSize = N / M; // Basic size of each interval
  int remainder = N % M;    // Remainder to distribute
  intervals[0] = 1;         // Start of the first interval
  for (int i = 1; i <= M; i++) {
    intervals[i] = intervals[i - 1] + intervalSize + (remainder > 0 ? 1 : 0);
    if (remainder > 0) {
      remainder--;
    }
  }
}

template<class DeviceType>
int PairPODKokkos<DeviceType>::calculateNumberOfIntervals(int N, int intervalSize)
{
  int M = N / intervalSize;
  if (N % intervalSize != 0) {
    M++; // Add an additional interval to cover the remainder
  }

  return M;
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::grow_atoms(int Ni)
{
  if (Ni > nimax) {
    nimax = Ni;
    MemKK::realloc_kokkos(numij, "pair_pod:numij", nimax+1);
    MemKK::realloc_kokkos(ei, "pair_pod:ei", nimax);
    MemKK::realloc_kokkos(typeai, "pair_pod:typeai", nimax);
    int n = nimax * nelements * K3 * nrbfmax;
    MemKK::realloc_kokkos(sumU, "pair_pod:sumU", n);
    MemKK::realloc_kokkos(forcecoeff, "pair_pod:forcecoeff", n);
    MemKK::realloc_kokkos(bd, "pair_pod:bd", nimax * Mdesc);
    MemKK::realloc_kokkos(cb, "pair_pod:cb", nimax * Mdesc);
    if (nClusters > 1)MemKK::realloc_kokkos(pd, "pair_pod:pd", nimax * (1 + nComponents + 3*nClusters));
    Kokkos::deep_copy(numij, 0);
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::grow_pairs(int Nij)
{
  if (Nij > nijmax) {
    nijmax = Nij;
    MemKK::realloc_kokkos(rij, "pair_pod:r_ij", 3 * nijmax);
    MemKK::realloc_kokkos(fij, "pair_pod:f_ij", 3 * nijmax);
    MemKK::realloc_kokkos(idxi, "pair_pod:idxi", nijmax);
    MemKK::realloc_kokkos(ai, "pair_pod:ai", nijmax);
    MemKK::realloc_kokkos(aj, "pair_pod:aj", nijmax);
    MemKK::realloc_kokkos(ti, "pair_pod:ti", nijmax);
    MemKK::realloc_kokkos(tj, "pair_pod:tj", nijmax);
    MemKK::realloc_kokkos(rbf, "pair_pod:rbf", nijmax * nrbfmax);
    MemKK::realloc_kokkos(rbfx, "pair_pod:rbfx", nijmax * nrbfmax);
    MemKK::realloc_kokkos(rbfy, "pair_pod:rbfy", nijmax * nrbfmax);
    MemKK::realloc_kokkos(rbfz, "pair_pod:rbfz", nijmax * nrbfmax);
    int kmax = (K3 > ns) ? K3 : ns;
    MemKK::realloc_kokkos(abf, "pair_pod:abf", nijmax * kmax);
    MemKK::realloc_kokkos(abfx, "pair_pod:abfx", nijmax * kmax);
    MemKK::realloc_kokkos(abfy, "pair_pod:abfy", nijmax * kmax);
    MemKK::realloc_kokkos(abfz, "pair_pod:abfz", nijmax * kmax);
  }
}

template<class DeviceType>
int PairPODKokkos<DeviceType>::NeighborCount(t_pod_1i l_numij, double l_rcutsq, int gi1, int Ni)
{
  // create local shadow views for KOKKOS_LAMBDA to pass them into parallel_for
  auto l_ilist = d_ilist;
  auto l_x = x;
  auto l_numneigh = d_numneigh;
  auto l_neighbors = d_neighbors;

  // compute number of pairs for each atom i
  Kokkos::parallel_for("NeighborCount", Kokkos::TeamPolicy<>(Ni, Kokkos::AUTO), KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
    int i = team.league_rank();
    int gi = l_ilist(gi1 + i);
    double xi0 = l_x(gi, 0);
    double xi1 = l_x(gi, 1);
    double xi2 = l_x(gi, 2);
    int jnum = l_numneigh(gi);
    int ncount = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
        [&] (const int jj, int& count) {
      int j = l_neighbors(gi,jj);
      j &= NEIGHMASK;
      double delx = xi0 - l_x(j,0);
      double dely = xi1 - l_x(j,1);
      double delz = xi2 - l_x(j,2);
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < l_rcutsq) count++;
    },ncount);

    l_numij(i+1) = ncount;
  });

  // accumalative sum
  Kokkos::parallel_scan("InclusivePrefixSum", Ni + 1, KOKKOS_LAMBDA(int i, int& update, const bool final) {
    if (i > 0) {
      update += l_numij(i);
      if (final) {
        l_numij(i) = update;
      }
    }
  });

  int total_neighbors = 0;
  Kokkos::deep_copy(Kokkos::View<int,Kokkos::HostSpace>(&total_neighbors), Kokkos::subview(l_numij, Ni));

  return total_neighbors;
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::NeighborList(t_pod_1d l_rij, t_pod_1i l_numij,  t_pod_1i l_typeai,
  t_pod_1i l_idxi, t_pod_1i l_ai, t_pod_1i l_aj, t_pod_1i l_ti, t_pod_1i l_tj, double l_rcutsq, int gi1, int Ni)
{
  // create local shadow views for KOKKOS_LAMBDA to pass them into parallel_for
  auto l_ilist = d_ilist;
  auto l_x = x;
  auto l_numneigh = d_numneigh;
  auto l_neighbors = d_neighbors;
  auto l_map = d_map;
  auto l_type = type;

  Kokkos::parallel_for("NeighborList", Kokkos::TeamPolicy<>(Ni, Kokkos::AUTO), KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
    int i = team.league_rank();
    int gi = l_ilist(gi1 + i);
    double xi0 = l_x(gi, 0);
    double xi1 = l_x(gi, 1);
    double xi2 = l_x(gi, 2);
    int itype = l_map(l_type(gi)) + 1; //map[atomtypes[gi]] + 1;
    l_typeai(i) = itype;
    int jnum = l_numneigh(gi);
    int nij0 = l_numij(i);
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,jnum),
        [&] (const int jj, int& offset, bool final) {
      int gj = l_neighbors(gi,jj);
      gj &= NEIGHMASK;
      double delx = l_x(gj,0) - xi0;
      double dely = l_x(gj,1) - xi1;
      double delz = l_x(gj,2) - xi2;
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= l_rcutsq) return;
      if (final) {
        int nij1 = nij0 + offset;
        l_rij(nij1 * 3 + 0) = delx;
        l_rij(nij1 * 3 + 1) = dely;
        l_rij(nij1 * 3 + 2) = delz;
        l_idxi(nij1) = i;
        l_ai(nij1) = gi;
        l_aj(nij1) = gj;
        l_ti(nij1) = itype;
        l_tj(nij1) = l_map(l_type(gj)) + 1; //map[atomtypes[gj)) + 1;
      }
      offset++;
    });
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::radialbasis(t_pod_1d rbft, t_pod_1d rbftx, t_pod_1d rbfty, t_pod_1d rbftz,
    t_pod_1d l_rij, t_pod_1d l_besselparams, double l_rin, double l_rmax, int l_besseldegree,
    int l_inversedegree, int l_nbesselpars, int Nij)
{
  Kokkos::parallel_for("ComputeRadialBasis", Nij, KOKKOS_LAMBDA(int n) {
    double xij1 = l_rij(0+3*n);
    double xij2 = l_rij(1+3*n);
    double xij3 = l_rij(2+3*n);

    double dij = sqrt(xij1*xij1 + xij2*xij2 + xij3*xij3);
    double dr1 = xij1/dij;
    double dr2 = xij2/dij;
    double dr3 = xij3/dij;

    double r = dij - l_rin;
    double y = r/l_rmax;
    double y2 = y*y;

    double y3 = 1.0 - y2*y;
    double y4 = y3*y3 + 1e-6;
    double y5 = sqrt(y4);
    double y6 = exp(-1.0/y5);
    double y7 = y4*sqrt(y4);

    // Calculate the final cutoff function as y6/exp(-1)
    double fcut = y6/exp(-1.0);

    // Calculate the derivative of the final cutoff function
    double dfcut = ((3.0/(l_rmax*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;

    // Calculate fcut/r, fcut/r^2, and dfcut/r
    double f1 = fcut/r;
    double f2 = f1/r;
    double df1 = dfcut/r;

    double alpha = l_besselparams(0);
    double t1 = (1.0-exp(-alpha));
    double t2 = exp(-alpha*r/l_rmax);
    double x0 =  (1.0 - t2)/t1;
    double dx0 = (alpha/l_rmax)*t2/t1;

    alpha = l_besselparams(1);
    t1 = (1.0-exp(-alpha));
    t2 = exp(-alpha*r/l_rmax);
    double x1 =  (1.0 - t2)/t1;
    double dx1 = (alpha/l_rmax)*t2/t1;

    alpha = l_besselparams(2);
    t1 = (1.0-exp(-alpha));
    t2 = exp(-alpha*r/l_rmax);
    double x2 =  (1.0 - t2)/t1;
    double dx2 = (alpha/l_rmax)*t2/t1;

    for (int i=0; i<l_besseldegree; i++) {
      double a = (i+1)*MY_PI;
      double b = (sqrt(2.0/(l_rmax))/(i+1));
      double af1 = a*f1;

      double sinax = sin(a*x0);
      int idxni = n + Nij*i;
      rbft(idxni) = b*f1*sinax;
      double drbftdr = b*(df1*sinax - f2*sinax + af1*cos(a*x0)*dx0);
      rbftx(idxni) = drbftdr*dr1;
      rbfty(idxni) = drbftdr*dr2;
      rbftz(idxni) = drbftdr*dr3;

      sinax = sin(a*x1);
      idxni = n + Nij*i + Nij*l_besseldegree*1;
      rbft(idxni) = b*f1*sinax;
      drbftdr = b*(df1*sinax - f2*sinax + af1*cos(a*x1)*dx1);
      rbftx(idxni) = drbftdr*dr1;
      rbfty(idxni) = drbftdr*dr2;
      rbftz(idxni) = drbftdr*dr3;

      sinax = sin(a*x2);
      idxni = n + Nij*i + Nij*l_besseldegree*2;
      rbft(idxni) = b*f1*sinax;
      drbftdr = b*(df1*sinax - f2*sinax + af1*cos(a*x2)*dx2);
      rbftx(idxni) = drbftdr*dr1;
      rbfty(idxni) = drbftdr*dr2;
      rbftz(idxni) = drbftdr*dr3;
    }

    // Calculate fcut/dij and dfcut/dij
    f1 = fcut/dij;
    double a = 1.0;
    for (int i=0; i<l_inversedegree; i++) {
      int p = l_besseldegree*l_nbesselpars + i;
      int idxni = n + Nij*p;
      a = a*dij;

      rbft(idxni) = fcut/a;

      double drbftdr = (dfcut - (i+1.0)*f1)/a;
      rbftx(idxni) = drbftdr*dr1;
      rbfty(idxni) = drbftdr*dr2;
      rbftz(idxni) = drbftdr*dr3;
    }
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::matrixMultiply(t_pod_1d a, t_pod_1d b, t_pod_1d c, int r1, int c1, int c2)
{
    Kokkos::parallel_for("MatrixMultiply", r1 * c2, KOKKOS_LAMBDA(int idx) {
        int j = idx / r1;  // Calculate column index
        int i = idx % r1;  // Calculate row index
        double sum = 0.0;
        for (int k = 0; k < c1; ++k) {
            sum += a(i + r1*k) * b(k + c1*j);  // Manually calculate the 1D index
        }
        c(i + r1*j) = sum;  // Manually calculate the 1D index for c
    });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::angularbasis(t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
        t_pod_1d l_rij, t_pod_1i l_pq3, int l_K3, int N)
{
  Kokkos::parallel_for("AngularBasis", N, KOKKOS_LAMBDA(int j) {
    double x = l_rij(j*3 + 0);
    double y = l_rij(j*3 + 1);
    double z = l_rij(j*3 + 2);

    double xx = x*x;
    double yy = y*y;
    double zz = z*z;
    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double dij = sqrt(xx + yy + zz);
    const double u = x / dij;
    const double v = y / dij;
    const double w = z / dij;

    double dij3 = dij * dij * dij;
    const double dudx = (yy + zz) / dij3;
    const double dudy = -xy / dij3;
    const double dudz = -xz / dij3;

    const double dvdx = -xy / dij3;
    const double dvdy = (xx + zz) / dij3;
    const double dvdz = -yz / dij3;

    const double dwdx = -xz / dij3;
    const double dwdy = -yz / dij3;
    const double dwdz = (xx + yy) / dij3;

    int idxa = j;
    l_abf(idxa) = 1.0;
    l_abfx(idxa) = 0.0;
    l_abfy(idxa) = 0.0;
    l_abfz(idxa) = 0.0;

    // Loop over all angular basis functions
    for (int n=1; n<l_K3; n++) {
      // Get indices for angular basis function
      int d = l_pq3(n + l_K3);
      int mj = j + N*(l_pq3(n)-1);
      idxa = j + N*n;
      // Calculate angular basis function and its derivatives using recursion relation
      if (d==1) {
        l_abf(idxa) = l_abf(mj)*u;
        l_abfx(idxa) = l_abfx(mj)*u + l_abf(mj);
        l_abfy(idxa) = l_abfy(mj)*u;
        l_abfz(idxa) = l_abfz(mj)*u;
      }
      else if (d==2) {
        l_abf(idxa) = l_abf(mj)*v;
        l_abfx(idxa) = l_abfx(mj)*v;
        l_abfy(idxa) = l_abfy(mj)*v + l_abf(mj);
        l_abfz(idxa) = l_abfz(mj)*v;
      }
      else if (d==3) {
        l_abf(idxa) = l_abf(mj)*w;
        l_abfx(idxa) = l_abfx(mj)*w;
        l_abfy(idxa) = l_abfy(mj)*w;
        l_abfz(idxa) = l_abfz(mj)*w + l_abf(mj);
      }
    }
    for (int n=1; n<l_K3; n++) {
      idxa = j + N*n;
      x = l_abfx(idxa);
      y = l_abfy(idxa);
      z = l_abfz(idxa);
      l_abfx(idxa) = x*dudx + y*dvdx + z*dwdx;
      l_abfy(idxa) = x*dudy + y*dvdy + z*dwdy;
      l_abfz(idxa) = x*dudz + y*dvdz + z*dwdz;
    }
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::radialangularsum(t_pod_1d l_sumU, t_pod_1d l_rbf, t_pod_1d l_abf, t_pod_1i l_tj,
    t_pod_1i l_numij, const int l_nelements, const int l_nrbf3, const int l_K3, const int Ni, const int Nij)
{
  int totalIterations = l_nrbf3 * l_K3 * Ni;
  if (l_nelements==1) {
    Kokkos::parallel_for("RadialAngularSum", totalIterations, KOKKOS_LAMBDA(int idx) {
      int k = idx % l_K3;
      int temp = idx / l_K3;
      int m = temp % l_nrbf3;
      int i = temp / l_nrbf3;
      int kmi = k + l_K3*m + l_K3*l_nrbf3*i;

      int start = l_numij(i);
      int nj = l_numij(i+1)-start;
      double sum=0.0;
      for (int j=0; j<nj; j++) {
        int n = start + j;
        sum += l_rbf(n + Nij * m) * l_abf(n + Nij * k);
      }
      l_sumU(kmi) = sum;
    });
  }
  else {
    Kokkos::parallel_for("RadialAngularSum", totalIterations, KOKKOS_LAMBDA(int idx) {
      int k = idx % l_K3;
      int temp = idx / l_K3;
      int m = temp % l_nrbf3;
      int i = temp / l_nrbf3;
      int kmi = l_nelements*k + l_nelements*l_K3*m + l_nelements*l_K3*l_nrbf3*i;
      int start = l_numij(i);
      int nj = l_numij(i+1)-start;

      double tm[10];
      for (int j=0; j<l_nelements; j++) tm[j] = 0;
      for (int j=0; j<nj; j++) {
        int n = start + j;
        int ia = n + Nij * k;
        int ib = n + Nij * m;
        int tn = l_tj(n) - 1; // offset the atom type by 1, since atomtype is 1-based
        tm[tn] += l_rbf(ib) * l_abf(ia);
      }
      for (int j=0; j<l_nelements; j++) l_sumU(j + kmi) = tm[j];
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::twobodydesc(t_pod_1d d2,  t_pod_1d l_rbf, t_pod_1i l_idxi, t_pod_1i l_tj,
        int l_nrbf2, const int Ni, const int Nij)
{
  int totalIterations = l_nrbf2 * Nij;
  Kokkos::parallel_for("twobodydesc", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx / l_nrbf2; // pair index
    int m = idx % l_nrbf2; // rbd index
    int i2 = n + Nij * m; // Index of the radial basis function for atom n and RBF m
    Kokkos::atomic_add(&d2(l_idxi(n) + Ni * (m + l_nrbf2 * (l_tj(n) - 1))), l_rbf(i2)); // Add the radial basis function to the corresponding descriptor
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::twobody_forces(t_pod_1d fij, t_pod_1d cb2, t_pod_1d l_rbfx, t_pod_1d l_rbfy,
        t_pod_1d l_rbfz, t_pod_1i l_idxi, t_pod_1i l_tj, int l_nrbf2, const int Ni, const int Nij)
{
  int totalIterations = l_nrbf2 * Nij;
  Kokkos::parallel_for("twobody_forces", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx / l_nrbf2; // pair index
    int m = idx % l_nrbf2; // rbd index
    int i2 = n + Nij * m; // Index of the radial basis function for atom n and RBF m
    int i1 = 3*n;
    double c = cb2(l_idxi(n) + Ni*m + Ni*l_nrbf2*(l_tj(n) - 1));
    Kokkos::atomic_add(&fij(0 + i1), c*l_rbfx(i2)); // Add the derivative with respect to x to the corresponding descriptor derivative
    Kokkos::atomic_add(&fij(1 + i1), c*l_rbfy(i2)); // Add the derivative with respect to y to the corresponding descriptor derivative
    Kokkos::atomic_add(&fij(2 + i1), c*l_rbfz(i2)); // Add the derivative with respect to z to the corresponding descriptor derivative
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::threebodydesc(t_pod_1d d3, t_pod_1d l_sumU, t_pod_1i l_pc3, t_pod_1i l_pn3,
        int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, const int Ni)
{
  int totalIterations = l_nrbf3 * Ni;
  Kokkos::parallel_for("ThreeBodyDesc", totalIterations, KOKKOS_LAMBDA(int idx) {
    int m = idx % l_nrbf3;
    int i = idx / l_nrbf3;
    int nmi = l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3*i;
    for (int p = 0; p < l_nabf3; p++) {
      int n1 = l_pn3(p);
      int n2 = l_pn3(p + 1);
      int nn = n2 - n1;
      int ipm = i + Ni * (p + l_nabf3 * m);
      int k = 0;
      for (int i1 = 0; i1 < l_nelements; i1++) {
        for (int i2 = i1; i2 < l_nelements; i2++) {
          double tmp=0;
          for (int q = 0; q < nn; q++) {
            tmp += l_pc3(n1 + q) * l_sumU(i1 + l_nelements * (n1 + q) + nmi) * l_sumU(i2 + l_nelements * (n1 + q) + nmi);
          }
          d3(ipm + totalIterations * l_nabf3 * k) = tmp;
          k += 1;
        }
      }
    }
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::threebody_forces(t_pod_1d fij, t_pod_1d cb3, t_pod_1d l_rbf, t_pod_1d l_rbfx,
    t_pod_1d l_rbfy, t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
    t_pod_1d l_sumU, t_pod_1i l_idxi, t_pod_1i l_tj, t_pod_1i l_pc3, t_pod_1i l_pn3, t_pod_1i l_elemindex,
    int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, int Ni, int Nij)
{
  int totalIterations = l_nrbf3 * Nij;
  if (l_nelements==1) {
    Kokkos::parallel_for("threebody_forces1", totalIterations, KOKKOS_LAMBDA(int idx) {
      int j = idx / l_nrbf3;       // Calculate j using integer division
      int m = idx % l_nrbf3;       // Calculate m using modulo operation
      int idxR = j + Nij * m;  // Pre-compute the index for rbf
      double rbfBase = l_rbf(idxR);
      double rbfxBase = l_rbfx(idxR);
      double rbfyBase = l_rbfy(idxR);
      double rbfzBase = l_rbfz(idxR);
      double fx = 0;
      double fy = 0;
      double fz = 0;
      for (int p = 0; p < l_nabf3; p++) {
        double c3 = 2.0 * cb3(l_idxi(j) + Ni*p + Ni*l_nabf3*m);
        int n1 = l_pn3(p);
        int nn = l_pn3(p + 1) - n1;
        int idxU = l_K3 * m + l_K3*l_nrbf3*l_idxi(j);
        for (int q = 0; q < nn; q++) {
          int idxNQ = n1 + q;  // Combine n1 and q into a single index for pc3 and sumU
          double f = c3 * l_pc3(idxNQ) * l_sumU(idxNQ + idxU);
          int idxA = j + Nij*idxNQ;  // Pre-compute the index for abf
          double abfA = l_abf(idxA);

          // Use the pre-computed indices to update dd3
          fx += f * (l_abfx(idxA) * rbfBase + rbfxBase * abfA);
          fy += f * (l_abfy(idxA) * rbfBase + rbfyBase * abfA);
          fz += f * (l_abfz(idxA) * rbfBase + rbfzBase * abfA);
        }
      }
      int ii = 3 * j;  // Pre-compute the base index for dd3
      Kokkos::atomic_add(&fij(0 + ii), fx); // Add the derivative with respect to x to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(1 + ii), fy); // Add the derivative with respect to y to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(2 + ii), fz); // Add the derivative with respect to z to the corresponding descriptor derivative
    });
  }
  else {
    int N3 = Ni *  l_nabf3 * l_nrbf3;
    Kokkos::parallel_for("threebody_forces2", totalIterations, KOKKOS_LAMBDA(int idx) {
      int j = idx / l_nrbf3;  // Derive the original j value
      int m = idx % l_nrbf3;  // Derive the original m value
      int i2 = l_tj(j) - 1;
      int idxK = l_nelements * l_K3 * m + l_nelements*l_K3*l_nrbf3*l_idxi(j);
      int idxR = j + Nij * m;  // Pre-compute the index for rbf
      double rbfBase = l_rbf(idxR);
      double rbfxBase = l_rbfx(idxR);
      double rbfyBase = l_rbfy(idxR);
      double rbfzBase = l_rbfz(idxR);
      double fx = 0;
      double fy = 0;
      double fz = 0;
      for (int p = 0; p < l_nabf3; p++) {
        int n1 = l_pn3(p);
        int nn = l_pn3(p + 1) - n1;
        int jmp = l_idxi(j) + Ni*(p + l_nabf3*m);
        for (int i1 = 0; i1 < l_nelements; i1++) {
          double c3 = (i1 == i2) ? 2.0 : 1.0;
          c3 = c3 * cb3(jmp + N3*l_elemindex(i2 + l_nelements * i1));
          for (int q = 0; q < nn; q++) {
            int idxNQ = n1 + q;  // Combine n1 and q into a single index
            int idxA = j + Nij*idxNQ;  // Pre-compute the index for abf
            double abfA = l_abf(idxA);
            double f = c3 * l_pc3(idxNQ) * l_sumU(i1 + l_nelements * idxNQ + idxK);
            fx += f * (l_abfx(idxA) * rbfBase + rbfxBase * abfA);
            fy += f * (l_abfy(idxA) * rbfBase + rbfyBase * abfA);
            fz += f * (l_abfz(idxA) * rbfBase + rbfzBase * abfA);
          }
        }
      }
      int ii = 3 * j;  // Pre-compute the base index for dd3
      Kokkos::atomic_add(&fij(0 + ii), fx); // Add the derivative with respect to x to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(1 + ii), fy); // Add the derivative with respect to y to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(2 + ii), fz); // Add the derivative with respect to z to the corresponding descriptor derivative
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::threebody_forcecoeff(t_pod_1d fb3, t_pod_1d cb3,
    t_pod_1d l_sumU, t_pod_1i l_pc3, t_pod_1i l_pn3, t_pod_1i l_elemindex,
    int l_nelements, int l_nrbf3, int l_nabf3, int l_K3, int Ni)
{
  int totalIterations = l_nrbf3 * Ni;
  if (l_nelements==1) {
    Kokkos::parallel_for("threebody_forcecoeff1", totalIterations, KOKKOS_LAMBDA(int idx) {
      int i = idx / l_nrbf3;       // Calculate j using integer division
      int m = idx % l_nrbf3;       // Calculate m using modulo operation
      for (int p = 0; p < l_nabf3; p++) {
        double c3 = 2.0 * cb3(i + Ni*p + Ni*l_nabf3*m);
        int n1 = l_pn3(p);
        int nn = l_pn3(p + 1) - n1;
        int idxU = l_K3 * m + l_K3*l_nrbf3*i;
        for (int q = 0; q < nn; q++) {
          int idxNQ = n1 + q;  // Combine n1 and q into a single index for pc3 and sumU
          fb3(idxNQ + idxU) += c3 * l_pc3(idxNQ) * l_sumU(idxNQ + idxU);
        }
      }
    });
  }
  else {
    int N3 = Ni *  l_nabf3 * l_nrbf3;
    Kokkos::parallel_for("threebody_forcecoeff2", totalIterations, KOKKOS_LAMBDA(int idx) {
      int i = idx / l_nrbf3;  // Derive the original j value
      int m = idx % l_nrbf3;  // Derive the original m value
      for (int p = 0; p < l_nabf3; p++) {
        int n1 = l_pn3(p);
        int nn = l_pn3(p + 1) - n1;
        int jmp = i + Ni*(p + l_nabf3*m);
        for (int q = 0; q < nn; q++) {
          int k = n1 + q;  // Combine n1 and q into a single index
          int idxU = l_nelements * k + l_nelements * l_K3 * m + l_nelements*l_K3*l_nrbf3*i;
          for (int i1 = 0; i1 < l_nelements; i1++) {
            double tm = l_pc3[k] * l_sumU[i1 + idxU];
            for (int i2 = i1; i2 < l_nelements; i2++) {
              int em = l_elemindex[i2 + l_nelements * i1];
              double t1 = tm * cb3[jmp + N3*em]; // Ni *  nabf3 * nrbf3 * nelements*(nelements+1)/2
              fb3[i2 + idxU] += t1;   // K3*nrbf3*Ni
              fb3[i1 + idxU] += l_pc3[k] * cb3[jmp + N3*em] * l_sumU[i2 + idxU];
            }
          }
        }
      }
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::fourbodydesc(t_pod_1d d4,  t_pod_1d l_sumU, t_pod_1i l_pa4, t_pod_1i l_pb4,
    t_pod_1i l_pc4, int l_nelements, int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni)
{
  int totalIterations = l_nrbf4 * Ni;
  Kokkos::parallel_for("fourbodydesc", totalIterations, KOKKOS_LAMBDA(int idx) {
    int m = idx % l_nrbf4;
    int i = idx / l_nrbf4;
    int idxU = l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * i;
    for (int p = 0; p < l_nabf4; p++) {
      int n1 = l_pa4(p);
      int n2 = l_pa4(p + 1);
      int nn = n2 - n1;
      int k = 0;
      for (int i1 = 0; i1 < l_nelements; i1++) {
        for (int i2 = i1; i2 < l_nelements; i2++) {
          for (int i3 = i2; i3 < l_nelements; i3++) {
            double tmp = 0.0;
            for (int q = 0; q < nn; q++) {
              int c = l_pc4(n1 + q);
              int j1 = l_pb4(n1 + q);
              int j2 = l_pb4(n1 + q + l_Q4);
              int j3 = l_pb4(n1 + q + 2 * l_Q4);
              tmp += c * l_sumU(idxU + i1 + l_nelements * j1) * l_sumU(idxU + i2 + l_nelements * j2) * l_sumU(idxU + i3 + l_nelements * j3);
            }
            int kk = p + l_nabf4 * m + l_nabf4 * l_nrbf4 * k;
            d4(i + Ni * kk) = tmp;
            k += 1;
          }
        }
      }
    }
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::fourbody_forces(t_pod_1d fij, t_pod_1d cb4, t_pod_1d l_rbf, t_pod_1d l_rbfx,
    t_pod_1d l_rbfy, t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
    t_pod_1d l_sumU, t_pod_1i l_idxi, t_pod_1i l_tj, t_pod_1i l_pa4, t_pod_1i l_pb4, t_pod_1i l_pc4,
    int l_nelements, int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni, int Nij)
{
  int totalIterations = l_nrbf4 * Nij;
  if (l_nelements==1) {
    Kokkos::parallel_for("fourbody_forces1", totalIterations, KOKKOS_LAMBDA(int idx) {
      int j = idx / l_nrbf4;  // Derive the original j value
      int m = idx % l_nrbf4;  // Derive the original m value
      int idxU = l_K3 * m + l_K3*l_nrbf3*l_idxi(j);
      int baseIdxJ = j + Nij * m;  // Pre-compute the index for rbf
      double rbfBase = l_rbf(baseIdxJ);
      double rbfxBase = l_rbfx(baseIdxJ);
      double rbfyBase = l_rbfy(baseIdxJ);
      double rbfzBase = l_rbfz(baseIdxJ);
      double fx = 0;
      double fy = 0;
      double fz = 0;
      for (int p = 0; p < l_nabf4; p++) {
        int n1 = l_pa4(p);
        int n2 = l_pa4(p + 1);
        int nn = n2 - n1;
        double c4 = cb4(l_idxi(j) + Ni*p + Ni*l_nabf4*m);
        for (int q = 0; q < nn; q++) {
          int idxNQ = n1 + q;  // Combine n1 and q into a single index
          double c = c4 * l_pc4[idxNQ];
          int j1 = l_pb4(idxNQ);
          int j2 = l_pb4(idxNQ + l_Q4);
          int j3 = l_pb4(idxNQ + 2 * l_Q4);
          double c1 = l_sumU(idxU + j1);
          double c2 = l_sumU(idxU + j2);
          double c3 = l_sumU(idxU + j3);
          double t12 = c * c1 * c2;
          double t13 = c * c1 * c3;
          double t23 = c * c2 * c3;

          // Pre-calculate commonly used indices
          int baseIdxJ3 = j + Nij * j3; // Common index for j3 terms
          int baseIdxJ2 = j + Nij * j2; // Common index for j2 terms
          int baseIdxJ1 = j + Nij * j1; // Common index for j1 terms

          // Temporary variables to store repeated calculations
          double abfBaseJ1 = l_abf(baseIdxJ1);
          double abfBaseJ2 = l_abf(baseIdxJ2);
          double abfBaseJ3 = l_abf(baseIdxJ3);
          // Update dd4 using pre-computed indices
          fx += t12 * (l_abfx(baseIdxJ3) * rbfBase + rbfxBase * abfBaseJ3)
                            + t13 * (l_abfx(baseIdxJ2) * rbfBase + rbfxBase * abfBaseJ2)
                            + t23 * (l_abfx(baseIdxJ1) * rbfBase + rbfxBase * abfBaseJ1);
          fy += t12 * (l_abfy(baseIdxJ3) * rbfBase + rbfyBase * abfBaseJ3)
                            + t13 * (l_abfy(baseIdxJ2) * rbfBase + rbfyBase * abfBaseJ2)
                            + t23 * (l_abfy(baseIdxJ1) * rbfBase + rbfyBase * abfBaseJ1);
          fz += t12 * (l_abfz(baseIdxJ3) * rbfBase + rbfzBase * abfBaseJ3)
                            + t13 * (l_abfz(baseIdxJ2) * rbfBase + rbfzBase * abfBaseJ2)
                            + t23 * (l_abfz(baseIdxJ1) * rbfBase + rbfzBase * abfBaseJ1);
        }
      }
      int ii = 3 * j;  // Pre-compute the base index for dd3
      Kokkos::atomic_add(&fij(0 + ii), fx); // Add the derivative with respect to x to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(1 + ii), fy); // Add the derivative with respect to y to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(2 + ii), fz); // Add the derivative with respect to z to the corresponding descriptor derivative
    });
  }
  else {
    int N3 = Ni * l_nabf4 * l_nrbf4;
    Kokkos::parallel_for("fourbody_forces2", totalIterations, KOKKOS_LAMBDA(int idx) {
      int j = idx / l_nrbf4;  // Derive the original j value
      int m = idx % l_nrbf4;  // Derive the original m value
      int idxM = j + Nij * m;
      double rbfM = l_rbf(idxM);
      double rbfxM = l_rbfx(idxM);
      double rbfyM = l_rbfy(idxM);
      double rbfzM = l_rbfz(idxM);
      int typej = l_tj(j) - 1;
      double fx = 0;
      double fy = 0;
      double fz = 0;
      for (int p = 0; p < l_nabf4; p++)  {
        int n1 = l_pa4(p);
        int n2 = l_pa4(p + 1);
        int nn = n2 - n1;
        int jpm = l_idxi(j) + Ni*p + Ni*l_nabf4*m;
        int k = 0;
        for (int i1 = 0; i1 < l_nelements; i1++) {
          for (int i2 = i1; i2 < l_nelements; i2++) {
            for (int i3 = i2; i3 < l_nelements; i3++) {
              for (int q = 0; q < nn; q++) {
                double c = l_pc4(n1 + q) * cb4(jpm + N3*k);
                int j1 = l_pb4(n1 + q);
                int j2 = l_pb4(n1 + q + l_Q4);
                int j3 = l_pb4(n1 + q + 2 * l_Q4);

                int idx1 = i1 + l_nelements * j1 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * l_idxi(j);
                int idx2 = i2 + l_nelements * j2 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * l_idxi(j);
                int idx3 = i3 + l_nelements * j3 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * l_idxi(j);
                double c1 = l_sumU(idx1);
                double c2 = l_sumU(idx2 );
                double c3 = l_sumU(idx3);
                double t12 = c*(c1 * c2);
                double t13 = c*(c1 * c3);
                double t23 = c*(c2 * c3);

                int idxJ3 = j + Nij * j3;
                int idxJ2 = j + Nij * j2;
                int idxJ1 = j + Nij * j1;
                double abfJ1 = l_abf(idxJ1);
                double abfJ2 = l_abf(idxJ2);
                double abfJ3 = l_abf(idxJ3);
                double abfxJ1 = l_abfx(idxJ1);
                double abfxJ2 = l_abfx(idxJ2);
                double abfxJ3 = l_abfx(idxJ3);
                double abfyJ1 = l_abfy(idxJ1);
                double abfyJ2 = l_abfy(idxJ2);
                double abfyJ3 = l_abfy(idxJ3);
                double abfzJ1 = l_abfz(idxJ1);
                double abfzJ2 = l_abfz(idxJ2);
                double abfzJ3 = l_abfz(idxJ3);

                // Compute contributions for each condition
                if (typej == i3) {
                    fx += t12 * (abfxJ3 * rbfM + rbfxM * abfJ3);
                    fy += t12 * (abfyJ3 * rbfM + rbfyM * abfJ3);
                    fz += t12 * (abfzJ3 * rbfM + rbfzM * abfJ3);
                }
                if (typej == i2) {
                    fx += t13 * (abfxJ2 * rbfM + rbfxM * abfJ2);
                    fy += t13 * (abfyJ2 * rbfM + rbfyM * abfJ2);
                    fz += t13 * (abfzJ2 * rbfM + rbfzM * abfJ2);
                }
                if (typej == i1) {
                    fx += t23 * (abfxJ1 * rbfM + rbfxM * abfJ1);
                    fy += t23 * (abfyJ1 * rbfM + rbfyM * abfJ1);
                    fz += t23 * (abfzJ1 * rbfM + rbfzM * abfJ1);
                }
              }
              k += 1;
            }
          }
        }
      }
      int ii = 3 * j;  // Pre-compute the base index for dd3
      Kokkos::atomic_add(&fij(0 + ii), fx); // Add the derivative with respect to x to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(1 + ii), fy); // Add the derivative with respect to y to the corresponding descriptor derivative
      Kokkos::atomic_add(&fij(2 + ii), fz); // Add the derivative with respect to z to the corresponding descriptor derivative
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::fourbody_forcecoeff(t_pod_1d fb4, t_pod_1d cb4,
    t_pod_1d l_sumU, t_pod_1i l_pa4, t_pod_1i l_pb4, t_pod_1i l_pc4, int l_nelements,
        int l_nrbf3, int l_nrbf4, int l_nabf4, int l_K3, int l_Q4, int Ni)
{
  int totalIterations = l_nrbf4 * Ni;
  if (l_nelements==1) {
    Kokkos::parallel_for("fourbody_forcecoeff1", totalIterations, KOKKOS_LAMBDA(int idx) {
      int i = idx / l_nrbf4;  // Derive the original j value
      int m = idx % l_nrbf4;  // Derive the original m value
      int idxU = l_K3 * m + l_K3*l_nrbf3*i;
      for (int p = 0; p < l_nabf4; p++) {
        int n1 = l_pa4(p);
        int n2 = l_pa4(p + 1);
        int nn = n2 - n1;
        double c4 = cb4(i + Ni*p + Ni*l_nabf4*m);
        for (int q = 0; q < nn; q++) {
          int idxNQ = n1 + q;  // Combine n1 and q into a single index
          double c = c4 * l_pc4[idxNQ];
          int j1 = idxU + l_pb4(idxNQ);
          int j2 = idxU + l_pb4(idxNQ + l_Q4);
          int j3 = idxU + l_pb4(idxNQ + 2 * l_Q4);
          double c1 = l_sumU(j1);
          double c2 = l_sumU(j2);
          double c3 = l_sumU(j3);
          fb4[j3] += c * c1 * c2;
          fb4[j2] += c * c1 * c3;
          fb4[j1] += c * c2 * c3;
        }
      }
    });
  }
  else {
    int N3 = Ni * l_nabf4 * l_nrbf4;
    Kokkos::parallel_for("fourbody_forcecoeff2", totalIterations, KOKKOS_LAMBDA(int idx) {
      int i = idx / l_nrbf4;  // Derive the original j value
      int m = idx % l_nrbf4;  // Derive the original m value
      for (int p = 0; p < l_nabf4; p++)  {
        int n1 = l_pa4(p);
        int n2 = l_pa4(p + 1);
        int nn = n2 - n1;
        int jpm = i + Ni*p + Ni*l_nabf4*m;
        for (int q = 0; q < nn; q++) {
          double c = l_pc4(n1 + q);
          int j1 = l_pb4(n1 + q);
          int j2 = l_pb4(n1 + q + l_Q4);
          int j3 = l_pb4(n1 + q + 2 * l_Q4);
          int idx1 = l_nelements * j1 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * i;
          int idx2 = l_nelements * j2 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * i;
          int idx3 = l_nelements * j3 + l_nelements * l_K3 * m + l_nelements * l_K3 * l_nrbf3 * i;
          int k = 0;
          for (int i1 = 0; i1 < l_nelements; i1++) {
            double c1 = l_sumU[idx1 + i1];
            for (int i2 = i1; i2 < l_nelements; i2++) {
              double c2 = l_sumU[idx2 + i2];
              for (int i3 = i2; i3 < l_nelements; i3++) {
                double c3 = l_sumU[idx3 + i3];
                double c4 = c * cb4[jpm + N3*k];
                fb4[idx3 + i3] += c4*(c1 * c2);
                fb4[idx2 + i2] += c4*(c1 * c3);
                fb4[idx1 + i1] += c4*(c2 * c3);
                k += 1;
              }
            }
          }
        }
      }
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::allbody_forces(t_pod_1d fij, t_pod_1d l_forcecoeff, t_pod_1d l_rbf, t_pod_1d l_rbfx,
    t_pod_1d l_rbfy, t_pod_1d l_rbfz, t_pod_1d l_abf, t_pod_1d l_abfx, t_pod_1d l_abfy, t_pod_1d l_abfz,
    t_pod_1i l_idxi, t_pod_1i l_tj, int l_nelements, int l_nrbf3, int l_K3, int Nij)
{
  int totalIterations = l_nrbf3 * Nij;
  Kokkos::parallel_for("allbody_forces", totalIterations, KOKKOS_LAMBDA(int idx) {
    int j = idx / l_nrbf3;       // Calculate j using integer division
    int m = idx % l_nrbf3;       // Calculate m using modulo operation
    int i2 = l_tj(j) - 1;
    int idxR = j + Nij * m;  // Pre-compute the index for rbf
    double rbfBase = l_rbf(idxR);
    double rbfxBase = l_rbfx(idxR);
    double rbfyBase = l_rbfy(idxR);
    double rbfzBase = l_rbfz(idxR);
    double fx = 0;
    double fy = 0;
    double fz = 0;
    for (int k = 0; k < l_K3; k++) {
      int idxU = l_nelements * k + l_nelements * l_K3 * m + l_nelements*l_K3*l_nrbf3*l_idxi[j];
      double fc = l_forcecoeff[i2 + idxU];
      int idxA = j + Nij*k;  // Pre-compute the index for abf
      double abfA = l_abf[idxA];
      double abfxA = l_abfx[idxA];
      double abfyA = l_abfy[idxA];
      double abfzA = l_abfz[idxA];
      fx += fc * (abfxA * rbfBase + rbfxBase * abfA); // K3*nrbf3*Nij
      fy += fc * (abfyA * rbfBase + rbfyBase * abfA);
      fz += fc * (abfzA * rbfBase + rbfzBase * abfA);
    }
    int ii = 3 * j;  // Pre-compute the base index for dd3
    Kokkos::atomic_add(&fij(0 + ii), fx); // Add the derivative with respect to x to the corresponding descriptor derivative
    Kokkos::atomic_add(&fij(1 + ii), fy); // Add the derivative with respect to y to the corresponding descriptor derivative
    Kokkos::atomic_add(&fij(2 + ii), fz); // Add the derivative with respect to z to the corresponding descriptor derivative
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::crossdesc(t_pod_1d d12, t_pod_1d d1, t_pod_1d d2, t_pod_1i ind1, t_pod_1i ind2, int n12, int Ni)
{
  int totalIterations = n12 * Ni;
  Kokkos::parallel_for("crossdesc", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx % Ni;
    int i = idx / Ni;

    d12(n + Ni * i) = d1(n + Ni * ind1(i)) * d2(n + Ni * ind2(i));
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::crossdesc_reduction(t_pod_1d cb1, t_pod_1d cb2, t_pod_1d c12, t_pod_1d d1,
        t_pod_1d d2, t_pod_1i ind1, t_pod_1i ind2, int n12, int Ni)
{
  int totalIterations = n12 * Ni;
  Kokkos::parallel_for("crossdesc_reduction", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx % Ni; // Ni
    int m = idx / Ni; // n12
    int k1 = ind1(m); // dd1
    int k2 = ind2(m); // dd2
    int m1 = n + Ni * k1; // d1
    int m2 = n + Ni * k2; // d2
    double c = c12(n + Ni * m);
    Kokkos::atomic_add(&cb1(m1), c * d2(m2));
    Kokkos::atomic_add(&cb2(m2), c * d1(m1));
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::set_array_to_zero(t_pod_1d a, int N)
{
  Kokkos::parallel_for("initialize_array", N, KOKKOS_LAMBDA(int i) {
    a(i) = 0.0;
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::blockatom_base_descriptors(t_pod_1d bd, int Ni, int Nij)
{
  auto begin = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();

  auto d2 = Kokkos::subview(bd, std::make_pair(0, Ni * nl2));
  auto d3 = Kokkos::subview(bd, std::make_pair(Ni * nl2, Ni * (nl2 + nl3)));
  auto d4 = Kokkos::subview(bd, std::make_pair(Ni * (nl2 + nl3), Ni * (nl2 + nl3 + nl4)));
  auto d33 = Kokkos::subview(bd, std::make_pair(Ni * (nl2 + nl3 + nl4), Ni * (nl2 + nl3 + nl4 + nl33)));
  auto d34 = Kokkos::subview(bd, std::make_pair(Ni * (nl2 + nl3 + nl4 + nl33), Ni * (nl2 + nl3 + nl4 + nl33 + nl34)));
  auto d44 = Kokkos::subview(bd, std::make_pair(Ni * (nl2 + nl3 + nl4 + nl33 + nl34), Ni * (nl2 + nl3 + nl4 + nl33 + nl34 + nl44)));

  begin = std::chrono::high_resolution_clock::now();
  radialbasis(abf, abfx, abfy, abfz, rij, besselparams, rin, rmax,
        besseldegree, inversedegree, nbesselpars, Nij);
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[10] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  begin = std::chrono::high_resolution_clock::now();
  matrixMultiply(abf,  Phi, rbf, Nij, ns,  nrbfmax);
  matrixMultiply(abfx, Phi, rbfx, Nij, ns,  nrbfmax);
  matrixMultiply(abfy, Phi, rbfy, Nij, ns,  nrbfmax);
  matrixMultiply(abfz, Phi, rbfz, Nij, ns,  nrbfmax);
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[11] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  begin = std::chrono::high_resolution_clock::now();
  set_array_to_zero(d2, Ni*nl2);
  twobodydesc(d2, rbf, idxi, tj, nrbf2, Ni, Nij);
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[12] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  if ((nl3 > 0) && (Nij>1)) {
    begin = std::chrono::high_resolution_clock::now();
    angularbasis(abf, abfx, abfy, abfz, rij, pq3, K3, Nij);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[13] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    begin = std::chrono::high_resolution_clock::now();
    set_array_to_zero(sumU, nelements * nrbf3 * K3 * Ni);
    radialangularsum(sumU, rbf, abf, tj, numij, nelements, nrbf3, K3, Ni, Nij);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[14] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

    begin = std::chrono::high_resolution_clock::now();
    //set_array_to_zero(d3, Ni*nl3);
    threebodydesc(d3, sumU, pc3, pn3, nelements, nrbf3, nabf3, K3, Ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[15] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;
  }

  if ((nl4 > 0) && (Nij>2)) {
    begin = std::chrono::high_resolution_clock::now();
    //set_array_to_zero(d4, Ni*nl4);
    fourbodydesc(d4, sumU, pa4, pb4, pc4, nelements, nrbf3, nrbf4, nabf4, K3, Q4, Ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[16] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;
  }

  if ((nl33>0) && (Nij>3)) {
    begin = std::chrono::high_resolution_clock::now();
    crossdesc(d33, d3, d3, ind33l, ind33r, nl33, Ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[17] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;
  }

  if ((nl34>0) && (Nij>4)) {
    begin = std::chrono::high_resolution_clock::now();
    crossdesc(d34, d3, d4, ind34l, ind34r, nl34, Ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[18] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;
  }

  if ((nl44>0) && (Nij>5)) {
    begin = std::chrono::high_resolution_clock::now();
    crossdesc(d44, d4, d4, ind44l, ind44r, nl44, Ni);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    comptime[19] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::blockatom_base_coefficients(t_pod_1d ei, t_pod_1d cb, t_pod_1d B, int Ni)
{
  auto cefs = coefficients;
  auto tyai = typeai;
  int nDes = Mdesc;
  int nCoeff = nCoeffPerElement;

  Kokkos::parallel_for("atomic_energies", Ni, KOKKOS_LAMBDA(int n) {
    int nc = nCoeff*(tyai[n]-1);
    ei[n] = cefs[0 + nc];
    for (int m=0; m<nDes; m++)
      ei[n] += cefs[1 + m + nc]*B[n + Ni*m];
  });

  int totalIterations = Ni*nDes;
  Kokkos::parallel_for("base_coefficients", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx % Ni;
    int m = idx / Ni;
    int nc = nCoeff*(tyai[n]-1);
    cb[n + Ni*m] = cefs[1 + m + nc];
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::blockatom_environment_descriptors(t_pod_1d ei, t_pod_1d cb, t_pod_1d B, int Ni)
{
  auto P = Kokkos::subview(pd, std::make_pair(0, Ni * nClusters));
  auto cp = Kokkos::subview(pd, std::make_pair(Ni * nClusters, 2 * Ni * nClusters));
  auto D = Kokkos::subview(pd, std::make_pair(2 * Ni * nClusters, 3 * Ni * nClusters));
  auto pca = Kokkos::subview(pd, std::make_pair(3 * Ni * nClusters, 3 * Ni * nClusters + Ni * nComponents));
  auto sumD = Kokkos::subview(pd, std::make_pair(3 * Ni * nClusters + Ni * nComponents, 3 * Ni * nClusters + Ni * nComponents + Ni));

  auto proj = Proj;
  auto cent = Centroids;
  auto cefs = coefficients;
  auto tyai = typeai;

  int nCom = nComponents;
  int nCls = nClusters;
  int nDes = Mdesc;
  int nCoeff = nCoeffPerElement;

  int totalIterations = Ni*nCom;
  Kokkos::parallel_for("pca", totalIterations, KOKKOS_LAMBDA(int idx) {
    int i = idx % Ni;
    int k = idx / Ni;
    double sum = 0.0;
    int typei = tyai[i]-1;
    for (int m = 0; m < nDes; m++) {
      sum += proj[k + nCom*m + nCom*nDes*typei] * B[i + Ni*m];
    }
    pca[i + Ni*k] = sum;
  });

  totalIterations = Ni*nCls;
  Kokkos::parallel_for("inverse_square_distances", totalIterations, KOKKOS_LAMBDA(int idx) {
    int i = idx % Ni;
    int j = idx / Ni;
    int typei = tyai[i]-1;
    double sum = 1e-20;
    for (int k = 0; k < nCom; k++) {
      double c = cent[k + j * nCom + nCls*nCom*typei];
      double p = pca[i + Ni*k];
      sum += (p - c) * (p - c);
    }
    D[i + Ni*j] = 1.0 / sum;
  });

  Kokkos::parallel_for("Probabilities", Ni, KOKKOS_LAMBDA(int i) {
    double sum = 0;
    for (int j = 0; j < nCls; j++) sum += D[i + Ni*j];
    sumD[i] = sum;
    for (int j = 0; j < nCls; j++) P[i + Ni*j] = D[i + Ni*j]/sum;
  });

  Kokkos::parallel_for("atomic_energies", Ni, KOKKOS_LAMBDA(int n) {
    int nc = nCoeff*(tyai[n]-1);
    ei[n] = cefs[0 + nc];
    for (int k = 0; k<nCls; k++)
      for (int m=0; m<nDes; m++)
        ei[n] += cefs[1 + m + nDes*k + nc]*B[n + Ni*m]*P[n + Ni*k];
  });

  Kokkos::parallel_for("env_coefficients", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx % Ni;
    int k = idx / Ni;
    int nc = nCoeff*(tyai[n]-1);
    double sum = 0;
    for (int m = 0; m<nDes; m++)
      sum += cefs[1 + m + k*nDes + nc]*B[n + Ni*m];
    cp[n + Ni*k] = sum;
  });

  totalIterations = Ni*nDes;
  Kokkos::parallel_for("base_coefficients", totalIterations, KOKKOS_LAMBDA(int idx) {
    int n = idx % Ni;
    int m = idx / Ni;
    int nc = nCoeff*(tyai[n]-1);
    double sum = 0.0;
    for (int k = 0; k<nCls; k++)
      sum += cefs[1 + m + k*nDes + nc]*P[n + Ni*k];
    cb[n + Ni*m] = sum;
  });

  Kokkos::parallel_for("base_env_coefficients", totalIterations, KOKKOS_LAMBDA(int idx) {
    int i = idx % Ni;
    int m = idx / Ni;
    int typei = tyai[i]-1;
    double S1 = 1/sumD[i];
    double S2 = sumD[i]*sumD[i];
    double sum = 0.0;
    for (int j=0; j<nCls; j++) {
      double dP_dB = 0.0;
      for (int k = 0; k < nCls; k++) {
        double dP_dD = -D[i + Ni*j] / S2;
        if (k==j) dP_dD += S1;
        double dD_dB = 0.0;
        double D2 = 2 * D[i + Ni*k] * D[i + Ni*k];
        for (int n = 0; n < nCom; n++) {
          double dD_dpca = D2 * (cent[n + k * nCom + nCls*nCom*typei] - pca[i + Ni*n]);
          dD_dB += dD_dpca * proj[n + m * nCom + nCom*nDes*typei];
        }
        dP_dB += dP_dD * dD_dB;
      }
      sum += cp[i + Ni*j]*dP_dB;
    }
    cb[i + Ni*m] += sum;
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::blockatom_energyforce(t_pod_1d l_ei, t_pod_1d l_fij, int Ni, int Nij)
{
  auto begin = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();

  // calculate base descriptors and their derivatives with respect to atom coordinates
  begin = std::chrono::high_resolution_clock::now();
  blockatom_base_descriptors(bd, Ni, Nij);
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[4] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  begin = std::chrono::high_resolution_clock::now();
  if (nClusters > 1) {
    blockatom_environment_descriptors(l_ei, cb, bd, Ni);
  }
  else {
    blockatom_base_coefficients(l_ei, cb, bd, Ni);
  }
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[5] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  begin = std::chrono::high_resolution_clock::now();
  auto d3 = Kokkos::subview(bd, std::make_pair(Ni * nl2, Ni * (nl2 + nl3)));
  auto d4 = Kokkos::subview(bd, std::make_pair(Ni * (nl2 + nl3), Ni * (nl2 + nl3 + nl4)));
  auto cb2 = Kokkos::subview(cb, std::make_pair(0, Ni * nl2));
  auto cb3 = Kokkos::subview(cb, std::make_pair(Ni * nl2, Ni * (nl2 + nl3)));
  auto cb4 = Kokkos::subview(cb, std::make_pair(Ni * (nl2 + nl3), Ni * (nl2 + nl3 + nl4)));
  auto cb33 = Kokkos::subview(cb, std::make_pair(Ni * (nl2 + nl3 + nl4), Ni * (nl2 + nl3 + nl4 + nl33)));
  auto cb34 = Kokkos::subview(cb, std::make_pair(Ni * (nl2 + nl3 + nl4 + nl33), Ni * (nl2 + nl3 + nl4 + nl33 + nl34)));
  auto cb44 = Kokkos::subview(cb, std::make_pair(Ni * (nl2 + nl3 + nl4 + nl33 + nl34), Ni * (nl2 + nl3 + nl4 + nl33 + nl34 + nl44)));

  if ((nl33>0) && (Nij>3)) {
    crossdesc_reduction(cb3, cb3, cb33, d3, d3, ind33l, ind33r, nl33, Ni);
  }
  if ((nl34>0) && (Nij>4)) {
    crossdesc_reduction(cb3, cb4, cb34, d3, d4, ind34l, ind34r, nl34, Ni);
  }
  if ((nl44>0) && (Nij>5)) {
    crossdesc_reduction(cb4, cb4, cb44, d4, d4, ind44l, ind44r, nl44, Ni);
  }
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[6] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  begin = std::chrono::high_resolution_clock::now();
  set_array_to_zero(l_fij, 3*Nij);
  if ((nl2 > 0) && (Nij>0)) twobody_forces(l_fij, cb2, rbfx, rbfy, rbfz, idxi, tj, nrbf2, Ni, Nij);
  Kokkos::fence();
  end = std::chrono::high_resolution_clock::now();
  comptime[7] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;

  set_array_to_zero(forcecoeff, nelements * nrbf3 * K3 * Ni);
  if ((nl3 > 0) && (Nij>1)) threebody_forcecoeff(forcecoeff, cb3, sumU,
          pc3, pn3, elemindex, nelements, nrbf3, nabf3, K3, Ni);
  if ((nl4 > 0) && (Nij>2)) fourbody_forcecoeff(forcecoeff, cb4, sumU,
      pa4, pb4, pc4, nelements, nrbf3, nrbf4, nabf4, K3, Q4, Ni);
  if ((nl3 > 0) && (Nij>1)) allbody_forces(l_fij, forcecoeff, rbf, rbfx, rbfy, rbfz, abf, abfx, abfy, abfz,
          idxi, tj, nelements, nrbf3, K3, Nij);
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::tallyforce(t_pod_1d l_fij, t_pod_1i l_ai, t_pod_1i l_aj, int Nij)
{
  auto l_f = f;
  Kokkos::parallel_for("TallyForce", Nij, KOKKOS_LAMBDA(int n) {
    int im = l_ai(n);
    int jm = l_aj(n);
    int n3 = 3*n;
    double fx = l_fij(n3 + 0);
    double fy = l_fij(n3 + 1);
    double fz = l_fij(n3 + 2);
    Kokkos::atomic_add(&l_f(im, 0), fx);
    Kokkos::atomic_add(&l_f(im, 1), fy);
    Kokkos::atomic_add(&l_f(im, 2), fz);
    Kokkos::atomic_sub(&l_f(jm, 0), fx);
    Kokkos::atomic_sub(&l_f(jm, 1), fy);
    Kokkos::atomic_sub(&l_f(jm, 2), fz);
  });
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::tallyenergy(t_pod_1d l_ei, int istart, int Ni)
{
  auto l_eatom = d_eatom;

  // For global energy tally
  if (eflag_global) {
    double local_eng_vdwl = 0.0;
    Kokkos::parallel_reduce("GlobalEnergyTally", Ni, KOKKOS_LAMBDA(int k, E_FLOAT& update) {
      update += l_ei(k);
    }, local_eng_vdwl);

    // Update global energy on the host after the parallel region
    eng_vdwl += local_eng_vdwl;
  }

  // For per-atom energy tally
  if (eflag_atom) {
    Kokkos::parallel_for("PerAtomEnergyTally", Ni, KOKKOS_LAMBDA(int k) {
      l_eatom(istart + k) += l_ei(k);
    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::tallystress(t_pod_1d l_fij, t_pod_1d l_rij, t_pod_1i l_ai, t_pod_1i l_aj, int Nij)
{
  auto l_vatom = d_vatom;

  if (vflag_global) {
    for (int j=0; j<3; j++) {
      F_FLOAT sum = 0.0;
      Kokkos::parallel_reduce("GlobalStressTally", Nij, KOKKOS_LAMBDA(int k, F_FLOAT& update) {
        int k3 = 3*k;
        update += l_rij(j + k3) * l_fij(j + k3);
      }, sum);
      virial[j] -= sum;
    }

    F_FLOAT sum = 0.0;
    Kokkos::parallel_reduce("GlobalStressTally", Nij, KOKKOS_LAMBDA(int k, F_FLOAT& update) {
      int k3 = 3*k;
      update += l_rij(k3) * l_fij(1 + k3);
    }, sum);
    virial[3] -= sum;

    sum = 0.0;
    Kokkos::parallel_reduce("GlobalStressTally", Nij, KOKKOS_LAMBDA(int k, F_FLOAT& update) {
      int k3 = 3*k;
      update += l_rij(k3) * l_fij(2 + k3);
    }, sum);
    virial[4] -= sum;

    sum = 0.0;
    Kokkos::parallel_reduce("GlobalStressTally", Nij, KOKKOS_LAMBDA(int k, F_FLOAT& update) {
      int k3 = 3*k;
      update += l_rij(1+k3) * l_fij(2+k3);
    }, sum);
    virial[5] -= sum;
  }

  if (vflag_atom) {
    Kokkos::parallel_for("PerAtomStressTally", Nij, KOKKOS_LAMBDA(int k) {
      int i = l_ai(k);
      int j = l_aj(k);
      int k3 = 3*k;
      double v_local[6];
      v_local[0] = -l_rij(k3) * l_fij(k3 + 0);
      v_local[1] = -l_rij(k3 + 1) * l_fij(k3 + 1);
      v_local[2] = -l_rij(k3 + 2) * l_fij(k3 + 2);
      v_local[3] = -l_rij(k3 + 0) * l_fij(k3 + 1);
      v_local[4] = -l_rij(k3 + 0) * l_fij(k3 + 2);
      v_local[5] = -l_rij(k3 + 1) * l_fij(k3 + 2);

      for (int d = 0; d < 6; ++d) {
        Kokkos::atomic_add(&l_vatom(i, d), 0.5 * v_local[d]);
      }

      for (int d = 0; d < 6; ++d) {
        Kokkos::atomic_add(&l_vatom(j, d), 0.5 * v_local[d]);
      }

    });
  }
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::savematrix2binfile(std::string filename, t_pod_1d d_A, int nrows, int ncols)
{
  auto A = Kokkos::create_mirror_view(d_A);
  Kokkos::deep_copy(A, d_A);

  FILE *fp = fopen(filename.c_str(), "wb");
  double sz[2];
  sz[0] = (double) nrows;
  sz[1] = (double) ncols;
  fwrite( reinterpret_cast<char*>( sz ), sizeof(double) * (2), 1, fp);
  fwrite( reinterpret_cast<char*>( A.data() ), sizeof(double) * (nrows*ncols), 1, fp);
  fclose(fp);
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::saveintmatrix2binfile(std::string filename, t_pod_1i d_A, int nrows, int ncols)
{
  auto A = Kokkos::create_mirror_view(d_A);
  Kokkos::deep_copy(A, d_A);

  FILE *fp = fopen(filename.c_str(), "wb");
  int sz[2];
  sz[0] = nrows;
  sz[1] = ncols;
  fwrite( reinterpret_cast<char*>( sz ), sizeof(int) * (2), 1, fp);
  fwrite( reinterpret_cast<char*>( A.data() ), sizeof(int) * (nrows*ncols), 1, fp);
  fclose(fp);
}

template<class DeviceType>
void PairPODKokkos<DeviceType>::savedatafordebugging()
{
  saveintmatrix2binfile("podkktypeai.bin", typeai, ni, 1);
  saveintmatrix2binfile("podkknumij.bin", numij, ni+1, 1);
  saveintmatrix2binfile("podkkai.bin", ai, nij, 1);
  saveintmatrix2binfile("podkkaj.bin", aj, nij, 1);
  saveintmatrix2binfile("podkkti.bin", ti, nij, 1);
  saveintmatrix2binfile("podkktj.bin", tj, nij, 1);
  saveintmatrix2binfile("podkkidxi.bin", idxi, nij, 1);
  savematrix2binfile("podkkrbf.bin", rbf, nrbfmax, nij);
  savematrix2binfile("podkkrbfx.bin", rbfx, nrbfmax, nij);
  savematrix2binfile("podkkrbfy.bin", rbfy, nrbfmax, nij);
  savematrix2binfile("podkkrbfz.bin", rbfz, nrbfmax, nij);
  int kmax = (K3 > ns) ? K3 : ns;
  savematrix2binfile("podkkabf.bin", abf,   kmax, nij);
  savematrix2binfile("podkkabfx.bin", abfx, kmax, nij);
  savematrix2binfile("podkkabfy.bin", abfy, kmax, nij);
  savematrix2binfile("podkkabfz.bin", abfz, kmax, nij);
  savematrix2binfile("podkkbd.bin", bd, ni, Mdesc);
  savematrix2binfile("podkksumU.bin", sumU, nelements * K3 * nrbfmax, ni);
  savematrix2binfile("podkkrij.bin", rij, 3, nij);
  savematrix2binfile("podkkfij.bin", fij, 3, nij);
  savematrix2binfile("podkkei.bin", ei, ni, 1);

  error->all(FLERR, "Save data and stop the run for debugging");
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPODKokkos<DeviceType>::memory_usage()
{
  double bytes = 0;

  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairPODKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairPODKokkos<LMPHostType>;
#endif
}
