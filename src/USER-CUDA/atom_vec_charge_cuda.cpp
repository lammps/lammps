/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "atom_vec_charge_cuda.h"
#include "comm_cuda_cu.h"
#include "atom_vec_charge_cuda_cu.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "user_cuda.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFEXTRA 1000
#define NCUDAEXCHANGE 12 //nextra x y z vx vy vz tag type mask image q

#define BUF_CFLOAT double
/* ---------------------------------------------------------------------- */

AtomVecChargeCuda::AtomVecChargeCuda(LAMMPS *lmp) : AtomVecCharge(lmp)
{
   cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

   maxsend=0;
   cudable=true;
   cuda_init_done=false;
   max_nsend=0;
   cu_copylist=NULL;
   copylist=NULL;
   copylist2=NULL;
}

void AtomVecChargeCuda::grow_copylist(int new_max_nsend)
{
  max_nsend=new_max_nsend;
  delete cu_copylist;
  delete [] copylist2;
  if(copylist) CudaWrapper_FreePinnedHostData((void*) copylist);
  copylist = (int*) CudaWrapper_AllocPinnedHostData(max_nsend*sizeof(int),false);
  copylist2 = new int[max_nsend];
  cu_copylist = new cCudaData<int, int, xx > (copylist, max_nsend);
}

void AtomVecChargeCuda::grow_send(int n,double** buf_send,int flag)  //need to be able to grow the comm send_buffer since the array sahll be copied from the gpu in whole
{
  int old_maxsend=*maxsend+BUFEXTRA;
  *maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
  {
    if(cuda->pinned)
    {
      double* tmp = new double[old_maxsend];
      memcpy((void*) tmp,(void*) *buf_send,old_maxsend*sizeof(double));
      if(*buf_send) CudaWrapper_FreePinnedHostData((void*) (*buf_send));
      *buf_send = (double*) CudaWrapper_AllocPinnedHostData((*maxsend+BUFEXTRA)*sizeof(double),false);
      memcpy(*buf_send,tmp,old_maxsend*sizeof(double));
      delete [] tmp;
    }
    else
    {
     *buf_send = (double *)
      memory->srealloc(*buf_send,(*maxsend+BUFEXTRA)*sizeof(double),
                       "comm:buf_send");
    }
  }
  else {
   if(cuda->pinned)
    {
      if(*buf_send) CudaWrapper_FreePinnedHostData((void*) (*buf_send));
      *buf_send = (double*) CudaWrapper_AllocPinnedHostData((*maxsend+BUFEXTRA)*sizeof(double),false);
    }
    else
    {
      memory->sfree(*buf_send);
      *buf_send = (double *) memory->smalloc((*maxsend+BUFEXTRA)*sizeof(double),
                                          "comm:buf_send");
    }
  }
}

void AtomVecChargeCuda::grow_both(int n)
{
  if(cuda->finished_setup)
  cuda->downloadAll();
  AtomVecCharge::grow(n);
  if(cuda->finished_setup)
  {
    cuda->checkResize();
    cuda->uploadAll();
  }
}

int AtomVecChargeCuda::pack_comm(int n, int* iswap, double *buf,
                             int pbc_flag, int *pbc) //usually this should not be called since comm->communicate handles the communication if only positions are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
          return AtomVecCharge::pack_comm(n,iswap,buf,pbc_flag,pbc);

        int m = Cuda_CommCuda_PackComm(&cuda->shared_data,n,*iswap,(void*) buf,pbc,pbc_flag);
        if((sizeof(X_CFLOAT)!=sizeof(double)) && m)
          m=(m+1)*sizeof(X_CFLOAT)/sizeof(double);
        return m;
}

int AtomVecChargeCuda::pack_comm_vel(int n, int* iswap, double *buf,
                             int pbc_flag, int *pbc) //usually this should not be called since comm->communicate handles the communication if only positions are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
          return AtomVecCharge::pack_comm_vel(n,iswap,buf,pbc_flag,pbc);

        int m = Cuda_CommCuda_PackCommVel(&cuda->shared_data,n,*iswap,(void*) buf,pbc,pbc_flag);
        if((sizeof(X_CFLOAT)!=sizeof(double)) && m)
          m=(m+1)*sizeof(X_CFLOAT)/sizeof(double);
        return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeCuda::unpack_comm(int n, int first, double *buf) //usually this should not be called since comm->communicate handles the communication if only positions are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
           {AtomVecCharge::unpack_comm(n,first,buf); return;}

  Cuda_CommCuda_UnpackComm(&cuda->shared_data,n,first,(void*)buf);
}

void AtomVecChargeCuda::unpack_comm_vel(int n, int first, double *buf) //usually this should not be called since comm->communicate handles the communication if only positions are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
           {AtomVecCharge::unpack_comm_vel(n,first,buf); return;}

  Cuda_CommCuda_UnpackCommVel(&cuda->shared_data,n,first,(void*)buf);
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeCuda::pack_reverse(int n, int first, double *buf) //usually this should not be called since comm->communicate handles the communication if only forces are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
          return AtomVecCharge::pack_reverse(n,first,buf);

  int i,m,last;
  cuda->cu_f->download();
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  cuda->cu_f->upload();
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeCuda::unpack_reverse(int n, int *list, double *buf)//usually this should not be called since comm->communicate handles the communication if only forces are exchanged
{
  if(not cuda->finished_setup || cuda->oncpu)
          {AtomVecCharge::unpack_reverse(n,list,buf); return;}

  int i,j,m;

  m = 0;
  cuda->cu_f->download();
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
  cuda->cu_f->upload();
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeCuda::pack_border(int n, int *iswap, double *buf,
                               int pbc_flag, int *pbc)
{
 if(not cuda->finished_setup || cuda->oncpu)
          return AtomVecCharge::pack_border(n,iswap,buf,pbc_flag,pbc);

        int m = Cuda_AtomVecChargeCuda_PackBorder(&cuda->shared_data,n,*iswap,(void*) buf,pbc,pbc_flag);

  return m;
}

int AtomVecChargeCuda::pack_border_vel(int n, int *iswap, double *buf,
                               int pbc_flag, int *pbc)
{
 if(not cuda->finished_setup || cuda->oncpu)
          return AtomVecCharge::pack_border_vel(n,iswap,buf,pbc_flag,pbc);

        int m = Cuda_AtomVecChargeCuda_PackBorderVel(&cuda->shared_data,n,*iswap,(void*) buf,pbc,pbc_flag);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeCuda::unpack_border(int n, int first, double *buf)
{
  if(not cuda->finished_setup || cuda->oncpu)
           {AtomVecCharge::unpack_border(n,first,buf); return;}
  while(atom->nghost+atom->nlocal+n>=cuda->shared_data.atom.nmax) //ensure there is enough space on device to unpack data
  {
          grow_both(0);
  }
  int flag=Cuda_AtomVecChargeCuda_UnpackBorder(&cuda->shared_data,n,first,(void*)buf);
  if(flag) {printf(" # CUDA: Error: Failed to unpack Border atoms (This might be a bug).\n");}
}

void AtomVecChargeCuda::unpack_border_vel(int n, int first, double *buf)
{
  if(not cuda->finished_setup || cuda->oncpu)
           {AtomVecCharge::unpack_border_vel(n,first,buf); return;}
  while(atom->nghost+atom->nlocal+n>=cuda->shared_data.atom.nmax) //ensure there is enough space on device to unpack data
  {
          grow_both(0);
  }
  int flag=Cuda_AtomVecChargeCuda_UnpackBorderVel(&cuda->shared_data,n,first,(void*)buf);
  if(flag) {printf(" # CUDA: Error: Failed to unpack Border atoms (This might be a bug).\n");}
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */


int AtomVecChargeCuda::pack_exchange(int dim, double *buf)
{
  if(cuda->oncpu)
          return AtomVecCharge::pack_exchange(dim,buf);

  if(not cuda_init_done||domain->box_change)
  {
          Cuda_AtomVecChargeCuda_Init(&cuda->shared_data);
          cuda_init_done=true;
  }
  double** buf_pointer=(double**) buf;
  if(*maxsend<atom->nghost || *buf_pointer==NULL)
  {
          grow_send(atom->nghost>*maxsend?atom->nghost:*maxsend,buf_pointer,0);
          *maxsend=atom->nghost>*maxsend?atom->nghost:*maxsend;
  }

  if(max_nsend==0) grow_copylist(200);

  int nsend_atoms = Cuda_AtomVecChargeCuda_PackExchangeList(&cuda->shared_data,*maxsend,dim,*buf_pointer);

  if(nsend_atoms>max_nsend) grow_copylist(nsend_atoms+100);
  if(nsend_atoms*NCUDAEXCHANGE>*maxsend)
  {
          grow_send((int) (nsend_atoms+100)*NCUDAEXCHANGE,buf_pointer,0);
          Cuda_AtomVecChargeCuda_PackExchangeList(&cuda->shared_data,*maxsend,dim,*buf_pointer);
  }

  int nlocal=atom->nlocal-nsend_atoms;

  for(int i=0;i<nsend_atoms;i++) copylist2[i]=1;
  for(int j=1;j<nsend_atoms+1;j++)
  {
          int i = static_cast <int> ((*buf_pointer)[j]);
          if(i>=nlocal) copylist2[i-nlocal]=-1;
  }

  int actpos=0;
  for(int j=1;j<nsend_atoms+1;j++)
  {
          int i = static_cast <int> ((*buf_pointer)[j]);
          if(i<nlocal)
          {
            while(copylist2[actpos]==-1) actpos++;
              copylist[j-1]=nlocal+actpos;
            actpos++;
          }
  }
  cu_copylist->upload();

  cuda->shared_data.atom.nlocal=nlocal;

  int m = Cuda_AtomVecChargeCuda_PackExchange(&cuda->shared_data,nsend_atoms,*buf_pointer,cu_copylist->dev_data());

  if (atom->nextra_grow)
  for(int j=0;j<nsend_atoms;j++)
  {
      int i=static_cast <int> ((*buf_pointer)[j+1]);
      int nextra=0;
      for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {

        int dm = modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&((*buf_pointer)[m]));
        m+=dm;
                  nextra+=dm;
                  if(i<nlocal)modify->fix[atom->extra_grow[iextra]]->copy_arrays(copylist[j],i,1);
        if(m>*maxsend)  grow_send(m,buf_pointer,1);
      }
      (*buf_pointer)[j+1] = nextra;
  }

  (*buf_pointer)[0] = nsend_atoms;
  atom->nlocal-=nsend_atoms;
  cuda->shared_data.atom.update_nlocal=2;

  if(m==1) return 0;//m is at least 1 in cuda since buf[0] contains number of atoms
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeCuda::unpack_exchange(double *buf)
{
  if(cuda->oncpu)
          return AtomVecCharge::unpack_exchange(buf);
  double *sublo,*subhi;

  int dim=cuda->shared_data.exchange_dim;
  if(domain->box_change)
  Cuda_AtomVecChargeCuda_Init(&cuda->shared_data);
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  int mfirst=0;
  for(int pi=0;pi<(comm->procgrid[dim]>2?2:1);pi++)
  {
  int nlocal = atom->nlocal;
  int nsend_atoms=static_cast<int> (buf[0]);
  if(nsend_atoms>max_nsend) grow_copylist(nsend_atoms+100);

  if (nlocal+nsend_atoms+atom->nghost>=atom->nmax) grow_both(nlocal+nsend_atoms*2+atom->nghost);
  int naccept = Cuda_AtomVecChargeCuda_UnpackExchange(&cuda->shared_data,nsend_atoms,buf,cu_copylist->dev_data());
  cu_copylist->download();
  int m = nsend_atoms*NCUDAEXCHANGE + 1;
  nlocal+=naccept;
  if (atom->nextra_grow)
  for(int j=0;j<nsend_atoms;j++)
  {
    if(copylist[j]>-1)
    {
                    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
                                      m += modify->fix[atom->extra_grow[iextra]]->
                                        unpack_exchange(copylist[j],&buf[m]);
    }
    else
    m+=static_cast <int> (buf[j+1]);
  }
  cuda->shared_data.atom.nlocal=nlocal;
  cuda->shared_data.atom.update_nlocal=2;
  atom->nlocal=nlocal;
  mfirst+=m;
  buf=&buf[m];
  }
  return mfirst;
}
