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

#include "comm_kokkos.h"
#include "kokkos.h"
#include "atom.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "atom_vec_kokkos.h"
#include "domain.h"
#include "atom_masks.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 10000
#define BUFEXTRA 1000

enum{SINGLE,MULTI};

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space
------------------------------------------------------------------------- */

CommKokkos::CommKokkos(LAMMPS *lmp) : CommBrick(lmp)
{
  sendlist = NULL;  // need to free this since parent allocated?
  k_sendlist = ArrayTypes<LMPDeviceType>::tdual_int_2d();

  // error check for disallow of OpenMP threads?

  // initialize comm buffers & exchange memory

  maxsend = BUFMIN;
  k_buf_send = ArrayTypes<LMPDeviceType>::
    tdual_xfloat_2d("comm:k_buf_send",(maxsend+BUFEXTRA+5)/6,6);
  buf_send = k_buf_send.view<LMPHostType>().ptr_on_device();

  maxrecv = BUFMIN;
  k_buf_recv = ArrayTypes<LMPDeviceType>::
    tdual_xfloat_2d("comm:k_buf_recv",(maxrecv+5)/6,6);
  buf_recv = k_buf_recv.view<LMPHostType>().ptr_on_device();

  k_exchange_sendlist = ArrayTypes<LMPDeviceType>::
    tdual_int_1d("comm:k_exchange_sendlist",100);
  k_exchange_copylist = ArrayTypes<LMPDeviceType>::
    tdual_int_1d("comm:k_exchange_copylist",100);
  k_count = ArrayTypes<LMPDeviceType>::tdual_int_1d("comm:k_count",1);
  k_sendflag = ArrayTypes<LMPDeviceType>::tdual_int_1d("comm:k_sendflag",100);

  // next line is bogus?

  memory->create(maxsendlist,maxswap,"comm:maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
  }
  memory->create_kokkos(k_sendlist,sendlist,maxswap,BUFMIN,"comm:sendlist");
}

/* ---------------------------------------------------------------------- */

CommKokkos::~CommKokkos()
{
  memory->destroy_kokkos(k_sendlist,sendlist);
  memory->destroy_kokkos(k_buf_send,buf_send);
  memory->destroy_kokkos(k_buf_recv,buf_recv);
}

/* ---------------------------------------------------------------------- */

void CommKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  exchange_comm_classic = lmp->kokkos->exchange_comm_classic;
  forward_comm_classic = lmp->kokkos->forward_comm_classic;
  exchange_comm_on_host = lmp->kokkos->exchange_comm_on_host;
  forward_comm_on_host = lmp->kokkos->forward_comm_on_host;

  CommBrick::init();
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(int dummy)
{

 if (!forward_comm_classic) {
    if (forward_comm_on_host) forward_comm_device<LMPHostType>(dummy);
    else forward_comm_device<LMPDeviceType>(dummy);
    return;
  }

  k_sendlist.sync<LMPHostType>();

  if (comm_x_only) {
    atomKK->sync(Host,X_MASK);
    atomKK->modified(Host,X_MASK);
  } else if (ghost_velocity) {
    atomKK->sync(Host,X_MASK | V_MASK);
    atomKK->modified(Host,X_MASK | V_MASK);
  } else {
    atomKK->sync(Host,ALL_MASK);
    atomKK->modified(Host,ALL_MASK);
  }

  CommBrick::forward_comm(dummy);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::forward_comm_device(int dummy)
{
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVecKokkos *avec = (AtomVecKokkos *) atom->avec;
  double **x = atom->x;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  k_sendlist.sync<DeviceType>();

  for (int iswap = 0; iswap < nswap; iswap++) {

    if (sendproc[iswap] != me) {
      if (comm_x_only) {
        atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
        if (size_forward_recv[iswap]) buf = x[firstrecv[iswap]];
        else buf = NULL;

        if (size_forward_recv[iswap]) {
            buf = atomKK->k_x.view<DeviceType>().ptr_on_device() + 
              firstrecv[iswap]*atomKK->k_x.view<DeviceType>().dimension_1();
            MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        }
        n = avec->pack_comm_kokkos(sendnum[iswap],k_sendlist,
                                   iswap,k_buf_send,pbc_flag[iswap],pbc[iswap]);

        if (n) {
          MPI_Send(k_buf_send.view<DeviceType>().ptr_on_device(),
                   n,MPI_DOUBLE,sendproc[iswap],0,world);
        }

        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
        atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::
                         space,X_MASK);
      } else if (ghost_velocity) {
        error->all(FLERR,"Ghost velocity forward comm not yet "
                   "implemented with Kokkos");
        if (size_forward_recv[iswap])
          MPI_Irecv(k_buf_recv.view<LMPHostType>().ptr_on_device(),
                    size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
      } else {
        if (size_forward_recv[iswap])
          MPI_Irecv(k_buf_recv.view<DeviceType>().ptr_on_device(),
                    size_forward_recv[iswap],MPI_DOUBLE,
                    recvproc[iswap],0,world,&request);
        n = avec->pack_comm_kokkos(sendnum[iswap],k_sendlist,iswap,
                                   k_buf_send,pbc_flag[iswap],pbc[iswap]);
        if (n)
          MPI_Send(k_buf_send.view<DeviceType>().ptr_on_device(),n,
                   MPI_DOUBLE,sendproc[iswap],0,world);
        if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
        avec->unpack_comm_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_recv);
      }

    } else {
      if (!ghost_velocity) {
        if (sendnum[iswap])
          n = avec->pack_comm_self(sendnum[iswap],k_sendlist,iswap,
                                   firstrecv[iswap],pbc_flag[iswap],pbc[iswap]);
      } else if (ghost_velocity) {
        error->all(FLERR,"Ghost velocity forward comm not yet "
                   "implemented with Kokkos");
        n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommKokkos::exchange()
{
  if (!exchange_comm_classic) {
    if (exchange_comm_on_host) exchange_device<LMPHostType>();
    else exchange_device<LMPDeviceType>();
    return;
  }

  atomKK->sync(Host,ALL_MASK);
  atomKK->modified(Host,ALL_MASK);

  CommBrick::exchange();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct BuildExchangeListFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  X_FLOAT _lo,_hi;
  typename AT::t_x_array _x;

  int _nlocal,_dim;
  typename AT::t_int_1d _nsend;
  typename AT::t_int_1d _sendlist;
  typename AT::t_int_1d _sendflag;


  BuildExchangeListFunctor(
      const typename AT::tdual_x_array x,
      const typename AT::tdual_int_1d sendlist,
      typename AT::tdual_int_1d nsend,
      typename AT::tdual_int_1d sendflag,int nlocal, int dim,
                X_FLOAT lo, X_FLOAT hi):
                _x(x.template view<DeviceType>()),
                _sendlist(sendlist.template view<DeviceType>()),
                _nsend(nsend.template view<DeviceType>()),
                _sendflag(sendflag.template view<DeviceType>()),
                _nlocal(nlocal),_dim(dim),
                _lo(lo),_hi(hi){
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if (_x(i,_dim) < _lo || _x(i,_dim) >= _hi) {
      const int mysend=Kokkos::atomic_fetch_add(&_nsend(0),1);
      if(mysend<_sendlist.dimension_0()) {
        _sendlist(mysend) = i;
        _sendflag(i) = 1;
      }
    } else
      _sendflag(i) = 0;
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::exchange_device()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *sublo,*subhi,*buf;
  MPI_Request request;
  MPI_Status status;
  AtomVecKokkos *avec = (AtomVecKokkos *) atom->avec;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()
  // clear ghost count and any ghost bonus data internal to AtomVec

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,ALL_MASK);

  // loop over dimensions
  for (int dim = 0; dim < 3; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;

    if (true) {
      if (k_sendflag.h_view.dimension_0()<nlocal) k_sendflag.resize(nlocal);
      k_count.h_view(0) = k_exchange_sendlist.h_view.dimension_0();
      while (k_count.h_view(0)>=k_exchange_sendlist.h_view.dimension_0()) {
        k_count.h_view(0) = 0;
        k_count.modify<LMPHostType>();
        k_count.sync<DeviceType>();

        BuildExchangeListFunctor<DeviceType> 
          f(atomKK->k_x,k_exchange_sendlist,k_count,k_sendflag,
            nlocal,dim,lo,hi);
        Kokkos::parallel_for(nlocal,f);
        DeviceType::fence();
        k_exchange_sendlist.modify<DeviceType>();
        k_sendflag.modify<DeviceType>();
        k_count.modify<DeviceType>();

        k_count.sync<LMPHostType>();
        if (k_count.h_view(0)>=k_exchange_sendlist.h_view.dimension_0()) {
          k_exchange_sendlist.resize(k_count.h_view(0)*1.1);
          k_exchange_copylist.resize(k_count.h_view(0)*1.1);
          k_count.h_view(0)=k_exchange_sendlist.h_view.dimension_0();
        }
      }
      k_exchange_sendlist.sync<LMPHostType>();
      k_sendflag.sync<LMPHostType>();

      int sendpos = nlocal-1;
      nlocal -= k_count.h_view(0);
      for(int i = 0; i < k_count.h_view(0); i++) {
        if (k_exchange_sendlist.h_view(i)<nlocal) {
          while (k_sendflag.h_view(sendpos)) sendpos--;
          k_exchange_copylist.h_view(i) = sendpos;
          sendpos--;
        } else
        k_exchange_copylist.h_view(i) = -1;
      }

      k_exchange_copylist.modify<LMPHostType>();
      k_exchange_copylist.sync<DeviceType>();
      nsend = 
        avec->pack_exchange_kokkos(k_count.h_view(0),k_buf_send,
                                   k_exchange_sendlist,k_exchange_copylist,
                                   ExecutionSpaceFromDevice<DeviceType>::
                                   space,dim,lo,hi);
      DeviceType::fence();

    } else {
      while (i < nlocal) {
        if (x[i][dim] < lo || x[i][dim] >= hi) {
          if (nsend > maxsend) grow_send_kokkos(nsend,1);
          nsend += avec->pack_exchange(i,&buf_send[nsend]);
          avec->copy(nlocal-1,i,1);
          nlocal--;
        } else i++;
      }
    }
    atom->nlocal = nlocal;

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    if (procgrid[dim] == 1) {
      nrecv = nsend;
      buf = buf_send;
      if (nrecv) {
        atom->nlocal=avec->
          unpack_exchange_kokkos(k_buf_send,nrecv,atom->nlocal,dim,lo,hi,
                                 ExecutionSpaceFromDevice<DeviceType>::space);
        DeviceType::fence();
      }
    } else {
      MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                   &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                     &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
        nrecv += nrecv2;
      }
      if (nrecv > maxrecv) grow_recv_kokkos(nrecv);

      MPI_Irecv(k_buf_recv.view<DeviceType>().ptr_on_device(),nrecv1,
                MPI_DOUBLE,procneigh[dim][1],0,
                world,&request);
      MPI_Send(k_buf_send.view<DeviceType>().ptr_on_device(),nsend,
               MPI_DOUBLE,procneigh[dim][0],0,world);
      MPI_Wait(&request,&status);

      if (procgrid[dim] > 2) {
        MPI_Irecv(k_buf_recv.view<DeviceType>().ptr_on_device()+nrecv1,
                  nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                  world,&request);
        MPI_Send(k_buf_send.view<DeviceType>().ptr_on_device(),nsend,
                 MPI_DOUBLE,procneigh[dim][1],0,world);
        MPI_Wait(&request,&status);
      }

      buf = buf_recv;
      if (nrecv) {
        atom->nlocal = avec->
          unpack_exchange_kokkos(k_buf_recv,nrecv,atom->nlocal,dim,lo,hi,
                                 ExecutionSpaceFromDevice<DeviceType>::space);
        DeviceType::fence();
      }
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list

  }

  atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::space,ALL_MASK);

  if (atom->firstgroupname) {
    /* this is not yet implemented with Kokkos */
    atomKK->sync(Host,ALL_MASK);
    atom->first_reorder();
    atomKK->modified(Host,ALL_MASK);
  }
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate, so don't need to explicitly
     call communicate routine on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommKokkos::borders()
{
  if (!exchange_comm_classic) {
    if (exchange_comm_on_host) borders_device<LMPHostType>();
    else borders_device<LMPDeviceType>();
    return;
  }

  atomKK->sync(Host,ALL_MASK);
  k_sendlist.modify<LMPHostType>();
  atomKK->modified(Host,ALL_MASK);

  CommBrick::borders();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct BuildBorderListFunctor {
	typedef DeviceType device_type;
	typedef ArrayTypes<DeviceType> AT;
  X_FLOAT lo,hi;
  typename AT::t_x_array x;
  int iswap,maxsendlist;
  int nfirst,nlast,dim;
  typename AT::t_int_2d sendlist;
  typename AT::t_int_1d nsend;

  BuildBorderListFunctor(typename AT::tdual_x_array _x, 
                         typename AT::tdual_int_2d _sendlist,
                         typename AT::tdual_int_1d _nsend,int _nfirst, 
                         int _nlast, int _dim,
                         X_FLOAT _lo, X_FLOAT _hi, int _iswap, 
                         int _maxsendlist):
    x(_x.template view<DeviceType>()),
    sendlist(_sendlist.template view<DeviceType>()),
    nsend(_nsend.template view<DeviceType>()),
    nfirst(_nfirst),nlast(_nlast),dim(_dim),
    lo(_lo),hi(_hi),iswap(_iswap),maxsendlist(_maxsendlist){}


  KOKKOS_INLINE_FUNCTION
  void operator() (DeviceType dev) const {
    const int chunk = ((nlast - nfirst + dev.league_size() - 1 ) / 
                       dev.league_size());
    const int teamstart = chunk*dev.league_rank() + nfirst;
    const int teamend = (teamstart + chunk) < nlast?(teamstart + chunk):nlast;
    int mysend = 0;
    for (int i=teamstart + dev.team_rank(); i<teamend; i+=dev.team_size()) {
      if (x(i,dim) >= lo && x(i,dim) <= hi) mysend++;
    }
    const int my_store_pos = dev.team_scan(mysend,&nsend(0));

    if (my_store_pos+mysend < maxsendlist) {
    mysend = my_store_pos;
      for(int i=teamstart + dev.team_rank(); i<teamend; i+=dev.team_size()){
        if (x(i,dim) >= lo && x(i,dim) <= hi) {
          sendlist(iswap,mysend++) = i;
        }
      }
    }
  }

  size_t shmem_size() const { return 1000u;}
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::borders_device() {
  int i,n,itype,iswap,dim,ineed,twoneed,smax,rmax;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  MPI_Status status;
  AtomVecKokkos *avec = (AtomVecKokkos *) atom->avec;

  ExecutionSpace exec_space = ExecutionSpaceFromDevice<DeviceType>::space;
  k_sendlist.modify<DeviceType>();
  atomKK->sync(exec_space,ALL_MASK);

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in list for use in future timesteps

      x = atom->x;
      if (style == SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
      else sendflag = 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (style == SINGLE) {
            typename ArrayTypes<DeviceType>::tdual_int_1d total_send("TS",1);
            total_send.h_view(0) = 0;
            if(exec_space == Device) {
              total_send.template modify<DeviceType>();
              total_send.template sync<LMPDeviceType>();
            }
            BuildBorderListFunctor<DeviceType> f(atomKK->k_x,k_sendlist,
                total_send,nfirst,nlast,dim,lo,hi,iswap,maxsendlist[iswap]);
            Kokkos::ParallelWorkRequest config((nlast-nfirst+127)/128,128);
            Kokkos::parallel_for(config,f);
            DeviceType::fence();
            total_send.template modify<DeviceType>();
            total_send.template sync<LMPHostType>();

            if(total_send.h_view(0) >= maxsendlist[iswap]) {
              grow_list(iswap,total_send.h_view(0));
              total_send.h_view(0) = 0;
              if(exec_space == Device) {
                total_send.template modify<LMPHostType>();
                total_send.template sync<LMPDeviceType>();
              }
              BuildBorderListFunctor<DeviceType> f(atomKK->k_x,k_sendlist,
                  total_send,nfirst,nlast,dim,lo,hi,iswap,maxsendlist[iswap]);
              Kokkos::ParallelWorkRequest config((nlast-nfirst+127)/128,128);
              Kokkos::parallel_for(config,f);
              DeviceType::fence();
              total_send.template modify<DeviceType>();
              total_send.template sync<LMPHostType>();
            }
            nsend = total_send.h_view(0);
          } else {
            error->all(FLERR,"Required border comm not yet "
                       "implemented with Kokkos\n");
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }

        } else {
          error->all(FLERR,"Required border comm not yet "
                     "implemented with Kokkos\n");
          if (style == SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            for (i = atom->nlocal; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }

      // pack up list of border atoms

      if (nsend*size_border > maxsend)
        grow_send_kokkos(nsend*size_border,0);
      if (ghost_velocity) {
        error->all(FLERR,"Required border comm not yet "
                   "implemented with Kokkos\n");
        n = avec->pack_border_vel(nsend,sendlist[iswap],buf_send,
                                  pbc_flag[iswap],pbc[iswap]);
      }
      else
        n = avec->
          pack_border_kokkos(nsend,k_sendlist,k_buf_send,iswap,
                             pbc_flag[iswap],pbc[iswap],exec_space);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,&status);
        if (nrecv*size_border > maxrecv) grow_recv_kokkos(nrecv*size_border);
        if (nrecv) MPI_Irecv(k_buf_recv.view<DeviceType>().ptr_on_device(),
                             nrecv*size_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(k_buf_send.view<DeviceType>().ptr_on_device(),n,
                        MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,&status);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

      // unpack buffer

      if (ghost_velocity) {
        error->all(FLERR,"Required border comm not yet "
                   "implemented with Kokkos\n");
        avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);
      }
      else
        if (sendproc[iswap] != me)
          avec->unpack_border_kokkos(nrecv,atom->nlocal+atom->nghost,
                                     k_buf_recv,exec_space);
        else
          avec->unpack_border_kokkos(nrecv,atom->nlocal+atom->nghost,
                                     k_buf_send,exec_space);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send_kokkos(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv_kokkos(max);

  // reset global->local map

  if (map_style) atom->map_set();
  if (exec_space == Host) k_sendlist.sync<LMPDeviceType>();
  atomKK->modified(exec_space,ALL_MASK);
  DeviceType::fence();
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommKokkos::grow_send_kokkos(int n, int flag, ExecutionSpace space)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  int maxsend_border = (maxsend+BUFEXTRA+5)/atom->avec->size_border + 2;
  if (flag) {
    if(space == Device)
      k_buf_send.modify<LMPDeviceType>();
    else
      k_buf_send.modify<LMPHostType>();

    k_buf_send.resize(maxsend_border,atom->avec->size_border);
    buf_send = k_buf_send.view<LMPHostType>().ptr_on_device();
  }
  else {
    k_buf_send = ArrayTypes<LMPDeviceType>::
      tdual_xfloat_2d("comm:k_buf_send",maxsend_border,atom->avec->size_border);
    buf_send = k_buf_send.view<LMPHostType>().ptr_on_device();
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommKokkos::grow_recv_kokkos(int n, ExecutionSpace space)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  int maxrecv_border = (maxrecv+BUFEXTRA+5)/atom->avec->size_border + 2;
  k_buf_recv = ArrayTypes<LMPDeviceType>::
    tdual_xfloat_2d("comm:k_buf_recv",maxrecv_border,atom->avec->size_border);
  buf_recv = k_buf_recv.view<LMPHostType>().ptr_on_device();
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommKokkos::grow_list(int iswap, int n)
{
  int size = static_cast<int> (BUFFACTOR * n);

  memory->grow_kokkos(k_sendlist,sendlist,maxswap,size,"comm:sendlist");

  for(int i=0;i<maxswap;i++) {
    maxsendlist[i]=size; sendlist[i]=&k_sendlist.view<LMPHostType>()(i,0);
  }
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void CommKokkos::grow_swap(int n)
{
  free_swap();
  allocate_swap(n);
  if (style == MULTI) {
    free_multi();
    allocate_multi(n);
  }

  maxswap = n;
  int size = MAX(k_sendlist.d_view.dimension_1(),BUFMIN);

  memory->grow_kokkos(k_sendlist,sendlist,maxswap,size,"comm:sendlist");

  memory->grow(maxsendlist,n,"comm:maxsendlist");
  for (int i=0;i<maxswap;i++) maxsendlist[i]=size;
}
