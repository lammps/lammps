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
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_snap_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "kokkos.h"
#include "sna.h"

#define MAXLINE 1024
#define MAXWORD 3

using namespace LAMMPS_NS;

// Outstanding issues with quadratic term
// 1. there seems to a problem with compute_optimized energy calc
// it does not match compute_regular, even when quadratic coeffs = 0

//static double t1 = 0.0;
//static double t2 = 0.0;
//static double t3 = 0.0;
//static double t4 = 0.0;
//static double t5 = 0.0;
//static double t6 = 0.0;
//static double t7 = 0.0;
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::PairSNAPKokkos(LAMMPS *lmp) : PairSNAP(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  vector_length = 8;
  k_cutsq = tdual_fparams("PairSNAPKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkos<DeviceType>::~PairSNAPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  if (diagonalstyle != 3)
    error->all(FLERR,"Must use diagonal style = 3 with pair snap/kk");

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->request(this,instance_me);

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) { // still need atomics, even though using a full neigh list
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else {
    error->all(FLERR,"Must use half neighbor list style with pair snap/kk");
  }
}

/* ---------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  int newton_pair = force->newton_pair;
  if (newton_pair == false)
    error->all(FLERR,"PairSNAPKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  int inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  /*
  for (int i = 0; i < nlocal; i++) {
    typename t_neigh_list::t_neighs neighs_i = neigh_list.get_neighs(i);
    const int num_neighs = neighs_i.get_num_neighs();
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }*/
  int max_neighs = 0;
  Kokkos::parallel_reduce("PairSNAPKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Experimental::Max<int>(max_neighs));

  snaKK.nmax = max_neighs;

  T_INT team_scratch_size = snaKK.size_team_scratch_arrays();
  T_INT thread_scratch_size = snaKK.size_thread_scratch_arrays();

  //printf("Sizes: %i %i\n",team_scratch_size/1024,thread_scratch_size/1024);
  int team_size_max = Kokkos::TeamPolicy<DeviceType>::team_size_max(*this);
  int vector_length = 8;
#ifdef KOKKOS_ENABLE_CUDA
  int team_size = 32;//max_neighs;
  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
#else
  int team_size = 1;
#endif

  EV_FLOAT ev;

  if (eflag) {
    if (neighflag == HALF) {
      typename Kokkos::TeamPolicy<DeviceType, TagPairSNAP<HALF,1> > policy(inum,team_size,vector_length);
      Kokkos::parallel_reduce(policy
          .set_scratch_size(1,Kokkos::PerThread(thread_scratch_size))
          .set_scratch_size(1,Kokkos::PerTeam(team_scratch_size))
        ,*this,ev);
    } else if (neighflag == HALFTHREAD) {
      typename Kokkos::TeamPolicy<DeviceType, TagPairSNAP<HALFTHREAD,1> > policy(inum,team_size,vector_length);
      Kokkos::parallel_reduce(policy
          .set_scratch_size(1,Kokkos::PerThread(thread_scratch_size))
          .set_scratch_size(1,Kokkos::PerTeam(team_scratch_size))
        ,*this,ev);
    }
  } else {
    if (neighflag == HALF) {
      typename Kokkos::TeamPolicy<DeviceType, TagPairSNAP<HALF,0> > policy(inum,team_size,vector_length);
      Kokkos::parallel_for(policy
          .set_scratch_size(1,Kokkos::PerThread(thread_scratch_size))
          .set_scratch_size(1,Kokkos::PerTeam(team_scratch_size))
        ,*this);
    } else if (neighflag == HALFTHREAD) {
      typename Kokkos::TeamPolicy<DeviceType, TagPairSNAP<HALFTHREAD,0> > policy(inum,team_size,vector_length);
      Kokkos::parallel_for(policy
          .set_scratch_size(1,Kokkos::PerThread(thread_scratch_size))
          .set_scratch_size(1,Kokkos::PerTeam(team_scratch_size))
        ,*this);
    }
  }

//static int step =0;
//step++;
//if (step%10==0)
//        printf(" %e %e %e %e %e (%e %e): %e\n",t1,t2,t3,t4,t5,t6,t7,t1+t2+t3+t4+t5);

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::allocate()
{
  PairSNAP::allocate();

  int n = atom->ntypes;
  d_map = Kokkos::View<T_INT*, DeviceType>("PairSNAPKokkos::map",n+1);
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairSNAP::init_one(i,j);
  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSNAPKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairSNAP::coeff(narg,arg);

  // Set up element lists

  d_radelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:radelem",nelements);
  d_wjelem = Kokkos::View<F_FLOAT*, DeviceType>("pair:wjelem",nelements);
  d_coeffelem = Kokkos::View<F_FLOAT**, Kokkos::LayoutRight, DeviceType>("pair:coeffelem",nelements,ncoeffall);

  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_wjelem = Kokkos::create_mirror_view(d_wjelem);
  auto h_coeffelem = Kokkos::create_mirror_view(d_coeffelem);
  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int ielem = 0; ielem < nelements; ielem++) {
    h_radelem(ielem) = radelem[ielem];
    h_wjelem(ielem) = wjelem[ielem];
    for (int jcoeff = 0; jcoeff < ncoeffall; jcoeff++) {
      h_coeffelem(ielem,jcoeff) = coeffelem[ielem][jcoeff];
    }
  }

  for (int i = 1; i <= atom->ntypes; i++) {
    h_map(i) = map[i];
  }

  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_wjelem,h_wjelem);
  Kokkos::deep_copy(d_coeffelem,h_coeffelem);
  Kokkos::deep_copy(d_map,h_map);

  // deallocate non-kokkos sna

  if (sna) {
    for (int tid = 0; tid<nthreads; tid++)
      delete sna[tid];
    delete [] sna;
    sna = NULL;
  }

  // allocate memory for per OpenMP thread data which
  // is wrapped into the sna class

  snaKK = SNAKokkos<DeviceType>(rfac0,twojmax,
                  diagonalstyle,use_shared_arrays,
                  rmin0,switchflag,bzeroflag);
    //if (!use_shared_arrays)
  snaKK.grow_rij(nmax);
  snaKK.init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAP<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType, TagPairSNAP<NEIGHFLAG,EVFLAG> >::member_type& team, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int ii = team.league_rank();
  const int i = d_ilist[ii];
  SNAKokkos<DeviceType> my_sna(snaKK,team);
  const double x_i = x(i,0);
  const double y_i = x(i,1);
  const double z_i = x(i,2);
  const int type_i = type[i];
  const int elem_i = d_map[type_i];
  const double radi = d_radelem[elem_i];

  const int num_neighs = d_numneigh[i];

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // wj = weights for neighbors of I within cutoff
  // rcutij = cutoffs for neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  //Kokkos::Timer timer;
  int ninside = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,num_neighs),
      [&] (const int jj, int& count) {
    Kokkos::single(Kokkos::PerThread(team), [&] (){
      T_INT j = d_neighbors(i,jj);
      const F_FLOAT dx = x(j,0) - x_i;
      const F_FLOAT dy = x(j,1) - y_i;
      const F_FLOAT dz = x(j,2) - z_i;

      const int type_j = type(j);
      const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
      const int elem_j = d_map[type_j];

      if ( rsq < rnd_cutsq(type_i,type_j) )
       count++;
    });
  },ninside);

  //t1 += timer.seconds(); timer.reset();

  if (team.team_rank() == 0)
  Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_neighs),
      [&] (const int jj, int& offset, bool final){
  //for (int jj = 0; jj < num_neighs; jj++) {
    T_INT j = d_neighbors(i,jj);
    const F_FLOAT dx = x(j,0) - x_i;
    const F_FLOAT dy = x(j,1) - y_i;
    const F_FLOAT dz = x(j,2) - z_i;

    const int type_j = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
    const int elem_j = d_map[type_j];

    if ( rsq < rnd_cutsq(type_i,type_j) ) {
      if (final) {
        my_sna.rij(offset,0) = dx;
        my_sna.rij(offset,1) = dy;
        my_sna.rij(offset,2) = dz;
        my_sna.inside[offset] = j;
        my_sna.wj[offset] = d_wjelem[elem_j];
        my_sna.rcutij[offset] = (radi + d_radelem[elem_j])*rcutfac;
      }
      offset++;
    }
  });

  //t2 += timer.seconds(); timer.reset();

  team.team_barrier();
  // compute Ui, Zi, and Bi for atom I
  my_sna.compute_ui(team,ninside);
  //t3 += timer.seconds(); timer.reset();
  team.team_barrier();
  my_sna.compute_zi(team);
  //t4 += timer.seconds(); timer.reset();
  team.team_barrier();

  if (quadraticflag) {
    my_sna.compute_bi(team);
    team.team_barrier();
    my_sna.copy_bi2bvec(team);
    team.team_barrier();
  }

  // for neighbors of I within cutoff:
  // compute dUi/drj and dBi/drj
  // Fij = dEi/dRj = -dEi/dRi => add to Fi, subtract from Fj

  Kokkos::View<double*,Kokkos::LayoutRight,DeviceType,Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    d_coeffi(d_coeffelem,elem_i,Kokkos::ALL);

  Kokkos::parallel_for (Kokkos::TeamThreadRange(team,ninside),
      [&] (const int jj) {
  //for (int jj = 0; jj < ninside; jj++) {
    int j = my_sna.inside[jj];

    //Kokkos::Timer timer2;
    my_sna.compute_duidrj(team,&my_sna.rij(jj,0),
                           my_sna.wj[jj],my_sna.rcutij[jj]);
    //t6 += timer2.seconds(); timer2.reset();
    my_sna.compute_dbidrj(team);
    //t7 += timer2.seconds(); timer2.reset();
    my_sna.copy_dbi2dbvec(team);

    Kokkos::single(Kokkos::PerThread(team), [&] (){
    F_FLOAT fij[3];

    fij[0] = 0.0;
    fij[1] = 0.0;
    fij[2] = 0.0;

    // linear contributions

    for (int k = 1; k <= ncoeff; k++) {
      double bgb = d_coeffi[k];
      fij[0] += bgb*my_sna.dbvec(k-1,0);
      fij[1] += bgb*my_sna.dbvec(k-1,1);
      fij[2] += bgb*my_sna.dbvec(k-1,2);
    }

    if (quadraticflag) {

      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        double bveci = my_sna.bvec[icoeff];
        double fack = d_coeffi[k]*bveci;
        double dbvecix = my_sna.dbvec(icoeff,0);
        double dbveciy = my_sna.dbvec(icoeff,1);
        double dbveciz = my_sna.dbvec(icoeff,2);
        fij[0] += fack*dbvecix;
        fij[1] += fack*dbveciy;
        fij[2] += fack*dbveciz;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          double facki = d_coeffi[k]*bveci;
          double fackj = d_coeffi[k]*my_sna.bvec[jcoeff];
          fij[0] += facki*my_sna.dbvec(jcoeff,0)+fackj*dbvecix;
          fij[1] += facki*my_sna.dbvec(jcoeff,1)+fackj*dbveciy;
          fij[2] += facki*my_sna.dbvec(jcoeff,2)+fackj*dbveciz;
          k++;
        }
      }
    }

    // Hard-coded ZBL potential
    //const double dx = my_sna.rij(jj,0);
    //const double dy = my_sna.rij(jj,1);
    //const double dz = my_sna.rij(jj,2);
    //const double fdivr = -1.5e6/pow(dx*dx + dy*dy + dz*dz,7.0);
    //fij[0] += dx*fdivr;
    //fij[1] += dy*fdivr;
    //fij[2] += dz*fdivr;

    //OK
    //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf SNAP-COMPARE: FIJ\n"
    //    ,x(i,0),x(i,1),x(i,2),x(j,0),x(j,1),x(j,2),fij[0],fij[1],fij[2] );
    a_f(i,0) += fij[0];
    a_f(i,1) += fij[1];
    a_f(i,2) += fij[2];
    a_f(j,0) -= fij[0];
    a_f(j,1) -= fij[1];
    a_f(j,2) -= fij[2];

    // tally global and per-atom virial contribution

    if (EVFLAG) {
      if (vflag_either) {
        v_tally_xyz<NEIGHFLAG>(ev,i,j,
          fij[0],fij[1],fij[2],
          -my_sna.rij(jj,0),-my_sna.rij(jj,1),
          -my_sna.rij(jj,2));
      }
    }

    });
  });
  //t5 += timer.seconds(); timer.reset();

  // tally energy contribution

  if (EVFLAG) {
    if (eflag_either) {

      if (!quadraticflag) {
        my_sna.compute_bi(team);
        team.team_barrier();
        my_sna.copy_bi2bvec(team);
        team.team_barrier();
      }

      // E = beta.B + 0.5*B^t.alpha.B
      // coeff[k] = beta[k-1] or
      // coeff[k] = alpha_ii or
      // coeff[k] = alpha_ij = alpha_ji, j != i

      Kokkos::single(Kokkos::PerTeam(team), [&] () {

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k

        double evdwl = d_coeffi[0];

        // linear contributions
        // could use thread vector range on this loop

        for (int k = 1; k <= ncoeff; k++)
          evdwl += d_coeffi[k]*my_sna.bvec[k-1];

        // quadratic contributions

        if (quadraticflag) {
          int k = ncoeff+1;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bveci = my_sna.bvec[icoeff];
            evdwl += 0.5*d_coeffi[k++]*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              evdwl += d_coeffi[k++]*bveci*my_sna.bvec[jcoeff];
            }
          }
        }

        //ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
        if (eflag_global) ev.evdwl += evdwl;
        if (eflag_atom) d_eatom[i] += evdwl;
      });
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::operator() (TagPairSNAP<NEIGHFLAG,EVFLAG>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAP<NEIGHFLAG,EVFLAG> >::member_type& team) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSNAP<NEIGHFLAG,EVFLAG>(), team, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const E_FLOAT v0 = delx*fx;
  const E_FLOAT v1 = dely*fy;
  const E_FLOAT v2 = delz*fz;
  const E_FLOAT v3 = delx*fy;
  const E_FLOAT v4 = delx*fz;
  const E_FLOAT v5 = dely*fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i,0) += 0.5*v0;
    a_vatom(i,1) += 0.5*v1;
    a_vatom(i,2) += 0.5*v2;
    a_vatom(i,3) += 0.5*v3;
    a_vatom(i,4) += 0.5*v4;
    a_vatom(i,5) += 0.5*v5;
    a_vatom(j,0) += 0.5*v0;
    a_vatom(j,1) += 0.5*v1;
    a_vatom(j,2) += 0.5*v2;
    a_vatom(j,3) += 0.5*v3;
    a_vatom(j,4) += 0.5*v4;
    a_vatom(j,5) += 0.5*v5;
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

template<class DeviceType>
double PairSNAPKokkos<DeviceType>::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);
  bytes += n*n*sizeof(double);
  bytes += 3*nmax*sizeof(double);
  bytes += nmax*sizeof(int);
  bytes += (2*ncoeffall)*sizeof(double);
  bytes += (ncoeff*3)*sizeof(double);
  bytes += snaKK.memory_usage();
  return bytes;
}
