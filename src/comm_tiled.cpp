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

#include "string.h"
#include "lmptype.h"
#include "comm_tiled.h"
#include "comm_brick.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

enum{SINGLE,MULTI};               // same as in Comm
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(LAMMPS *lmp) : Comm(lmp)
{
  error->all(FLERR,"Comm_style tiled is not yet supported");

  style = 1;
  layout = LAYOUT_UNIFORM;
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommTiled::CommTiled(LAMMPS *lmp, Comm *oldcomm) : Comm(*oldcomm)
{
  error->all(FLERR,"Comm_style tiled is not yet supported");

  style = 1;
  layout = oldcomm->layout;
  copy_arrays(oldcomm);
  init_buffers();
}

/* ---------------------------------------------------------------------- */

CommTiled::~CommTiled()
{
  free_swap();

  if (sendlist) for (int i = 0; i < nswap; i++) memory->destroy(sendlist[i]);
  memory->sfree(sendlist);
  memory->destroy(maxsendlist);

  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
   NOTE: if this is identical to CommBrick, put it into Comm ??
------------------------------------------------------------------------- */

void CommTiled::init_buffers()
{
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");

  nswap = 2 * domain->dimension;
  allocate_swap(nswap);

  //sendlist = (int **) memory->smalloc(nswap*sizeof(int *),"comm:sendlist");
  //memory->create(maxsendlist,nswap,"comm:maxsendlist");
  //for (int i = 0; i < nswap; i++) {
  //  maxsendlist[i] = BUFMIN;
  //  memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
  //}
}

/* ----------------------------------------------------------------------
   NOTE: if this is nearly identical to CommBrick, put it into Comm ??
------------------------------------------------------------------------- */

void CommTiled::init()
{
  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // temporary restrictions

  if (triclinic) 
    error->all(FLERR,"Cannot yet use comm_style tiled with triclinic box");
  if (domain->xperiodic || domain->yperiodic || 
      (domain->dimension == 2 && domain->zperiodic))
    error->all(FLERR,"Cannot yet use comm_style tiled with periodic box");
  if (mode == MULTI)
    error->all(FLERR,"Cannot yet use comm_style tiled with multi-mode comm");

  // comm_only = 1 if only x,f are exchanged in forward/reverse comm
  // comm_x_only = 0 if ghost_velocity since velocities are added

  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;
  if (ghost_velocity) comm_x_only = 0;

  // set per-atom sizes for forward/reverse/border comm
  // augment by velocity and fix quantities if needed

  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;

  if (ghost_velocity) size_forward += atom->avec->size_velocity;
  if (ghost_velocity) size_border += atom->avec->size_velocity;

  for (int i = 0; i < modify->nfix; i++)
    size_border += modify->fix[i]->comm_border;

  // maxexchange = max # of datums/atom in exchange communication
  // maxforward = # of datums in largest forward communication
  // maxreverse = # of datums in largest reverse communication
  // query pair,fix,compute,dump for their requirements
  // pair style can force reverse comm even if newton off

  maxexchange = BUFMIN + maxexchange_fix;
  maxforward = MAX(size_forward,size_border);
  maxreverse = size_reverse;

  if (force->pair) maxforward = MAX(maxforward,force->pair->comm_forward);
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse);

  for (int i = 0; i < modify->nfix; i++) {
    maxforward = MAX(maxforward,modify->fix[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->fix[i]->comm_reverse);
  }

  for (int i = 0; i < modify->ncompute; i++) {
    maxforward = MAX(maxforward,modify->compute[i]->comm_forward);
    maxreverse = MAX(maxreverse,modify->compute[i]->comm_reverse);
  }

  for (int i = 0; i < output->ndump; i++) {
    maxforward = MAX(maxforward,output->dump[i]->comm_forward);
    maxreverse = MAX(maxreverse,output->dump[i]->comm_reverse);
  }

  if (force->newton == 0) maxreverse = 0;
  if (force->pair) maxreverse = MAX(maxreverse,force->pair->comm_reverse_off);
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size and tiling
   sets nsendproc, nrecvproc, sendproc, recvproc
   sets sendother, sendself, pbc_flag, pbc, sendbox
------------------------------------------------------------------------- */

void CommTiled::setup()
{
  int i;

  int dimension;
  int *periodicity;
  double *prd,*sublo,*subhi,*boxlo,*boxhi;

  double cut = MAX(neighbor->cutneighmax,cutghostuser);

  dimension = domain->dimension;
  periodicity = domain->periodicity;
  prd = domain->prd;
  sublo = domain->sublo;
  subhi = domain->subhi;
  boxlo = domain->boxlo;
  boxhi = domain->boxhi;
  cutghost[0] = cutghost[1] = cutghost[2] = cut;
  
  if ((periodicity[0] && cut > prd[0]) ||
      (periodicity[1] && cut > prd[1]) ||
      (dimension == 3 && periodicity[2] && cut > prd[2]))
    error->all(FLERR,"Communication cutoff for comm_style tiled "
               "cannot exceed periodic box length");

  // NOTE: allocate overlap (to Nprocs?)
  // NOTE: allocate 2nd dim of sendproc, recvproc, sendbox
  // NOTE: func pointers for box_drop and box_other
  // NOTE: write box_drop and box_other methods
  // NOTE: for tiled, must do one-time gather of RCB cuts and proc boxes

  int *overlap;
  int noverlap,noverlap1,indexme;
  double lo1[3],hi1[3],lo2[3],hi2[3];
  int one,two;

  nswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    for (int iswap = 0; iswap < 2; iswap++) {

      // ghost box in lower direction

      one = 1;
      lo1[0] = sublo[0]; lo1[1] = sublo[1]; lo1[2] = sublo[2];
      hi1[0] = subhi[0]; hi1[1] = subhi[1]; hi1[2] = subhi[2];
      if (iswap == 0) {
        lo1[idim] = sublo[idim] - cut;
        hi1[idim] = sublo[idim];
      } else {
        lo1[idim] = subhi[idim];
        hi1[idim] = subhi[idim] + cut;
      }
      
      two = 0;
      if (iswap == 0 && periodicity[idim] && lo1[idim] < boxlo[idim]) two = 1;
      if (iswap == 1 && periodicity[idim] && hi1[idim] > boxhi[idim]) two = 1;

      if (two) {
        lo2[0] = sublo[0]; lo2[1] = sublo[1]; lo2[2] = sublo[2];
        hi2[0] = subhi[0]; hi2[1] = subhi[1]; hi2[2] = subhi[2];
        if (iswap == 0) {
          lo2[idim] = lo1[idim] + prd[idim];
          hi2[idim] = hi1[idim] + prd[idim];
          if (sublo[idim] == boxlo[idim]) {
            one = 0;
            hi2[idim] = boxhi[idim];
          }
        } else {
          lo2[idim] = lo1[idim] - prd[idim];
          hi2[idim] = hi1[idim] - prd[idim];
          if (subhi[idim] == boxhi[idim]) {
            one = 0;
            lo2[idim] = boxlo[idim];
          }
        }
      }

      indexme = -1;
      noverlap = 0;
      if (one) {
        if (layout == LAYOUT_UNIFORM) 
          box_drop_uniform(idim,lo1,hi1,noverlap,overlap,indexme);
        else if (layout == LAYOUT_NONUNIFORM) 
          box_drop_nonuniform(idim,lo1,hi1,noverlap,overlap,indexme);
        else
          box_drop_tiled(lo1,hi1,0,nprocs-1,noverlap,overlap,indexme);
      }

      noverlap1 = noverlap;
      if (two) {
        if (layout == LAYOUT_UNIFORM) 
          box_drop_uniform(idim,lo2,hi2,noverlap,overlap,indexme);
        else if (layout == LAYOUT_NONUNIFORM) 
          box_drop_nonuniform(idim,lo2,hi2,noverlap,overlap,indexme);
        else 
          box_drop_tiled(lo2,hi2,0,nprocs-1,noverlap,overlap,indexme);
      }

      // if this (self) proc is in overlap list, move it to end of list
      
      if (indexme >= 0) {
        int tmp = overlap[noverlap-1];
        overlap[noverlap-1] = overlap[indexme];
        overlap[indexme] = tmp;
      }

      // overlap how has list of noverlap procs
      // includes PBC effects

      if (overlap[noverlap-1] == me) sendself[nswap] = 1;
      else sendself[nswap] = 0;
      if (noverlap-sendself[nswap]) sendother[nswap] = 1;
      else sendother[nswap] = 0;

      nsendproc[nswap] = noverlap;
      for (i = 0; i < noverlap; i++) sendproc[nswap][i] = overlap[i];
      if (iswap == 0) {
        nrecvproc[nswap+1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[nswap+1][i] = overlap[i];
      } else {
        nrecvproc[nswap-1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[nswap-1][i] = overlap[i];
      }

      // compute sendbox for each of my sends
      // obox = intersection of ghostbox with other proc's sub-domain
      // sbox = what I need to send to other proc
      //      = sublo to MIN(sublo+cut,subhi) in idim, for iswap = 0
      //      = MIN(subhi-cut,sublo) to subhi in idim, for iswap = 1
      //      = obox in other 2 dims
      // if sbox touches sub-box boundaries in lower dims,
      //   extend sbox in those lower dims to include ghost atoms
      
      double oboxlo[3],oboxhi[3],sbox[6];

      for (i = 0; i < noverlap; i++) {
        pbc_flag[nswap][i] = 0;
        pbc[nswap][i][0] = pbc[nswap][i][1] = pbc[nswap][i][2] =
          pbc[nswap][i][3] = pbc[nswap][i][4] = pbc[nswap][i][5] = 0;
        
        if (layout == LAYOUT_UNIFORM) 
          box_other_uniform(overlap[i],oboxlo,oboxhi);
        else if (layout == LAYOUT_NONUNIFORM) 
          box_other_nonuniform(overlap[i],oboxlo,oboxhi);
        else
          box_other_tiled(overlap[i],oboxlo,oboxhi);
        
        if (i < noverlap1) {
          sbox[0] = MAX(oboxlo[0],lo1[0]);
          sbox[1] = MAX(oboxlo[1],lo1[1]);
          sbox[2] = MAX(oboxlo[2],lo1[2]);
          sbox[3] = MIN(oboxhi[0],hi1[0]);
          sbox[4] = MIN(oboxhi[1],hi1[1]);
          sbox[5] = MIN(oboxhi[2],hi1[2]);
        } else {
          pbc_flag[nswap][i] = 1;
          pbc[nswap][i][idim] = 1;
          sbox[0] = MAX(oboxlo[0],lo2[0]);
          sbox[1] = MAX(oboxlo[1],lo2[1]);
          sbox[2] = MAX(oboxlo[2],lo2[2]);
          sbox[3] = MIN(oboxhi[0],hi2[0]);
          sbox[4] = MIN(oboxhi[1],hi2[1]);
          sbox[5] = MIN(oboxhi[2],hi2[2]);
        }

        if (iswap == 0) {
          sbox[idim] = sublo[idim];
          sbox[3+idim] = MIN(sublo[idim]+cut,subhi[idim]);
        } else {
          sbox[idim] = MAX(subhi[idim]-cut,sublo[idim]);
          sbox[3+idim] = subhi[idim];
        }

        if (idim >= 1) {
          if (sbox[0] == sublo[0]) sbox[0] -= cut;
          if (sbox[4] == subhi[0]) sbox[4] += cut;
        }
        if (idim == 2) {
          if (sbox[1] == sublo[1]) sbox[1] -= cut;
          if (sbox[5] == subhi[1]) sbox[5] += cut;
        }
        
        memcpy(sendbox[nswap][i],sbox,6*sizeof(double));
      }

      nswap++;
    }
  }


  // reallocate requests and statuses to max of any swap

}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::forward_comm(int dummy)
{
  int i,irecv,n,nsend,nrecv;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (comm_x_only) {
      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(x[firstrecv[iswap][i]],size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        avec->pack_comm(sendnum[iswap][nsend],sendlist[iswap][nsend],
                        x[firstrecv[iswap][nrecv]],pbc_flag[iswap][nsend],
                        pbc[iswap][nsend]);
      }

      if (sendother[iswap]) MPI_Waitall(nrecv,requests,statuses);

    } else if (ghost_velocity) {
      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm_vel(sendnum[iswap][i],sendlist[iswap][i],
                                  buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        avec->pack_comm_vel(sendnum[iswap][nsend],sendlist[iswap][nsend],
                            buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        avec->unpack_comm_vel(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],
                              buf_send);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,&status);
          avec->unpack_comm_vel(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsendproc[iswap]; i++) {
          n = avec->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                              buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        avec->pack_comm(sendnum[iswap][nsend],sendlist[iswap][nsend],
                        buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        avec->unpack_comm(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],
                          buf_send);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,&status);
          avec->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                            &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiled::reverse_comm()
{
  int i,irecv,n,nsend,nrecv;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **f = atom->f;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (comm_f_only) {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++)
          MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i],MPI_DOUBLE,
                    sendproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nrecv; i++)
          MPI_Send(f[firstrecv[iswap][i]],size_reverse_send[iswap][i],
                   MPI_DOUBLE,recvproc[iswap][i],0,world);
      }

      if (sendself[iswap]) {
        avec->unpack_reverse(sendnum[iswap][nsend],sendlist[iswap][nsend],
                             f[firstrecv[iswap][nrecv]]);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,&status);
          avec->unpack_reverse(sendnum[iswap][irecv],sendlist[iswap][irecv],
                               &buf_recv[reverse_recv_offset[iswap][irecv]]);
        }
      }
      
    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++)
          MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i],MPI_DOUBLE,
                    sendproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nrecv; i++) {
          n = avec->pack_reverse(recvnum[iswap][i],firstrecv[iswap][i],
                                 buf_send);
          MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        avec->pack_reverse(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],
                           buf_send);
        avec->unpack_reverse(sendnum[iswap][nsend],sendlist[iswap][nsend],
                             buf_send);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,&status);
          avec->unpack_reverse(sendnum[iswap][irecv],sendlist[iswap][irecv],
                               &buf_recv[reverse_recv_offset[iswap][irecv]]);
        }
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

void CommTiled::exchange()
{
  // loop over atoms
  //   if not outside my box, continue
  //   find which proc it is in
  //   find which one of my touching procs it is, else lost
  //   make sure all atoms are "lost" that should be (e.g. outside non-PBC)
  //   add to list to send to that proc
  // loop over touching procs
  //   send buffer to them
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created per swap/proc that will be made
   as list is made, actually do communication
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommTiled::borders()
{
  int i,n,irecv,ngroup,nlast,nsend,nrecv,ncount,rmaxswap;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *bbox;
  double **x;
  MPI_Status status;
  AtomVec *avec = atom->avec;

  // smax = max size of single send in a swap/proc
  // rmax = max size of recvs from all procs for a swap

  int smax = 0;
  int rmax = 0;

  // loop over all swaps in all dimensions

  for (int iswap = 0; iswap < nswap; iswap++) {

    // find atoms within rectangles using <= and >=
    // for x-dim swaps, check owned atoms
    // for yz-dim swaps, check owned and ghost atoms
    // store sent atom indices in list for use in future timesteps
    // NOTE: assume SINGLE mode, add back in logic for MULTI mode later

    x = atom->x;

    for (i = 0; i < nsendproc[iswap]; i++) {
      bbox = sendbox[iswap][i];
      xlo = bbox[0]; xhi = bbox[1];
      ylo = bbox[2]; yhi = bbox[3];
      zlo = bbox[4]; zhi = bbox[5];

      ngroup = atom->nfirst;
      if (iswap < 2) nlast = atom->nlocal;
      else nlast = atom->nlocal + atom->nghost;

      ncount = 0;
      for (i = 0; i < ngroup; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (ncount == maxsendlist[iswap][i]) grow_list(iswap,i,ncount);
          sendlist[iswap][i][ncount++] = i;
        }
      for (i = atom->nlocal; i < nlast; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (ncount == maxsendlist[iswap][i]) grow_list(iswap,i,ncount);
          sendlist[iswap][i][ncount++] = i;
        }
      sendnum[iswap][i] = ncount;
      smax = MAX(smax,ncount);
    }

    // send sendnum counts to procs who recv from me except self
    // copy data to self if sendself is set

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (sendother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&recvnum[iswap][i],1,MPI_INT,
                  recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsend; i++)
        MPI_Send(&sendnum[iswap][i],1,MPI_INT,sendproc[iswap][i],0,world);
    }
    if (sendself[iswap]) recvnum[iswap][nrecv] = sendnum[iswap][nsend];
    if (sendother[iswap]) MPI_Waitall(nrecv,requests,statuses);

    // setup other per swap/proc values from sendnum and recvnum

    rmaxswap = 0;
    for (i = 0; i < nrecvproc[iswap]; i++) {
      rmaxswap += recvnum[iswap][i];
      size_forward_recv[iswap][i] = recvnum[iswap][i]*size_forward;
      size_reverse_send[iswap][i] = recvnum[iswap][i]*size_reverse;
      size_reverse_recv[iswap][i] = sendnum[iswap][i]*size_reverse;
      if (i == 0) {
        firstrecv[iswap][0] = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][i] = firstrecv[iswap][i-1] + recvnum[iswap][i-1];
        forward_recv_offset[iswap][i] = 
          forward_recv_offset[iswap][i-1] + recvnum[iswap][i-1];
      }
    }
    rmax = MAX(rmax,rmaxswap);

    // insure send/recv buffers are large enough for border comm

    if (smax*size_border > maxsend) grow_send(smax*size_border,0);
    if (rmax*size_border > maxrecv) grow_recv(rmax*size_border);

    // swap atoms with other procs using pack_border(), unpack_border()

    if (ghost_velocity) {
      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    recvnum[iswap][i]*size_border,
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsend; i++) {
          n = avec->pack_border_vel(sendnum[iswap][i],sendlist[iswap][i],
                                    buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                   sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        n = avec->pack_border_vel(sendnum[iswap][nsend],sendlist[iswap][nsend],
                                  buf_send,pbc_flag[iswap][nsend],
                                  pbc[iswap][nsend]);
        avec->unpack_border_vel(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],
                                buf_send);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,&status);
          avec->unpack_border(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                              &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                    recvnum[iswap][i]*size_border,
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        for (i = 0; i < nsend; i++) {
          n = avec->pack_border(sendnum[iswap][i],sendlist[iswap][i],
                                buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                   sendproc[iswap][i],0,world);
        }
      }

      if (sendself[iswap]) {
        n = avec->pack_border(sendnum[iswap][nsend],sendlist[iswap][nsend],
                              buf_send,pbc_flag[iswap][nsend],
                              pbc[iswap][nsend]);
        avec->unpack_border(recvnum[iswap][nsend],firstrecv[iswap][nsend],
                            buf_send);
      }

      if (sendother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,&status);
          avec->unpack_border(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                              &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }
    }

    // increment ghost atoms

    n = nrecvproc[iswap];
    atom->nghost += forward_recv_offset[iswap][n-1] + recvnum[iswap][n-1];
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_pair(Pair *pair)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = pair->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                            buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        pair->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                          &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = pair->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                          buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      pair->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_pair(Pair *pair)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = pair->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                    buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        pair->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                  &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = pair->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                  buf_send);
      pair->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_fix(Fix *fix)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = fix->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                           buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        fix->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                         &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = fix->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                         buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      fix->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_fix(Fix *fix)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = fix->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                   buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        fix->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                 &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = fix->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                 buf_send);
      fix->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
   NOTE: complicated b/c don't know # to recv a priori
------------------------------------------------------------------------- */

void CommTiled::forward_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   n = total datums for all atoms, allows for variable number/atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_variable_fix(Fix *fix)
{
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_compute(Compute *compute)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = compute->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                               buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        compute->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                             &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = compute->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                             buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      compute->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_compute(Compute *compute)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = compute->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                   buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        compute->
          unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                              &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = compute->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                 buf_send);
      compute->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],
                                   buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::forward_comm_dump(Dump *dump)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nrecvproc[iswap]; i++)
        MPI_Irecv(&buf_recv[forward_recv_offset[iswap][i]],
                  size_forward_recv[iswap][i],
                  MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nsendproc[iswap]; i++) {
        n = dump->pack_comm(sendnum[iswap][i],sendlist[iswap][i],
                            buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
        MPI_Send(buf_send,n*sendnum[iswap][i],MPI_DOUBLE,
                 sendproc[iswap][i],0,world);
      }
      for (i = 0; i < nrecvproc[iswap]; i++) {
        MPI_Waitany(nrecvproc[iswap],requests,&irecv,&status);
        dump->unpack_comm(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                          &buf_recv[forward_recv_offset[iswap][irecv]]);
      }

    } else {
      n = dump->pack_comm(sendnum[iswap][0],sendlist[iswap][0],
                          buf_send,pbc_flag[iswap][0],pbc[iswap][0]);
      dump->unpack_comm(recvnum[iswap][0],firstrecv[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   n = constant number of datums per atom
------------------------------------------------------------------------- */

void CommTiled::reverse_comm_dump(Dump *dump)
{
  int i,irecv,n;
  MPI_Status status;

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap][0] != me) {
      for (i = 0; i < nsendproc[iswap]; i++)
        MPI_Irecv(&buf_recv[reverse_recv_offset[iswap][i]],
                  size_reverse_recv[iswap][i],MPI_DOUBLE,
                  sendproc[iswap][i],0,world,&requests[i]);
      for (i = 0; i < nrecvproc[iswap]; i++) {
        n = dump->pack_reverse_comm(recvnum[iswap][i],firstrecv[iswap][i],
                                    buf_send);
        MPI_Send(buf_send,n*recvnum[iswap][i],MPI_DOUBLE,
                 recvproc[iswap][i],0,world);
      }
      for (i = 0; i < nsendproc[iswap]; i++) {
        MPI_Waitany(nsendproc[iswap],requests,&irecv,&status);
        dump->unpack_reverse_comm(sendnum[iswap][irecv],sendlist[iswap][irecv],
                                  &buf_recv[reverse_recv_offset[iswap][irecv]]);
      }

    } else {
      n = dump->pack_reverse_comm(recvnum[iswap][0],firstrecv[iswap][0],
                                  buf_send);
      dump->unpack_reverse_comm(sendnum[iswap][0],sendlist[iswap][0],buf_send);
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication of N values in array
------------------------------------------------------------------------- */

void CommTiled::forward_comm_array(int n, double **array)
{
}

/* ----------------------------------------------------------------------
   exchange info provided with all 6 stencil neighbors
------------------------------------------------------------------------- */

int CommTiled::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  int nrecv = n;
  return nrecv;
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   box is owned by me and extends in dim
------------------------------------------------------------------------- */

void CommTiled::box_drop_uniform(int dim, double *lo, double *hi, 
                                 int &noverlap, int *overlap, int &indexme)
{

}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
------------------------------------------------------------------------- */

void CommTiled::box_drop_nonuniform(int dim, double *lo, double *hi, 
                                    int &noverlap, int *overlap, int &indexme)
{
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
------------------------------------------------------------------------- */

void CommTiled::box_drop_tiled(double *lo, double *hi, 
                               int proclower, int procupper,
                               int &noverlap, int *overlap, int &indexme)
{
  // end recursion when partition is a single proc
  // add proc to overlap list

  if (proclower == procupper) {
    if (proclower == me) indexme = noverlap;
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use > and < criteria so does not include a box it only touches
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  double cut = tree[procmid].cut;
  int dim = tree[procmid].dim;
  
  if (lo[dim] < cut) 
    box_drop_tiled(lo,hi,proclower,procmid-1,noverlap,overlap,indexme);
  if (hi[dim] > cut)
    box_drop_tiled(lo,hi,procmid,procupper,noverlap,overlap,indexme);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommTiled::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,maxsend+bufextra,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiled::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiled::grow_list(int iswap, int iwhich, int n)
{
  maxsendlist[iswap][iwhich] = static_cast<int> (BUFFACTOR * n);
  memory->grow(sendlist[iswap][iwhich],maxsendlist[iswap][iwhich],
               "comm:sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void CommTiled::allocate_swap(int n)
{
  memory->create(sendnum,n,"comm:sendnum");
  memory->create(recvnum,n,"comm:recvnum");
  memory->create(sendproc,n,"comm:sendproc");
  memory->create(recvproc,n,"comm:recvproc");
  memory->create(size_forward_recv,n,"comm:size");
  memory->create(size_reverse_send,n,"comm:size");
  memory->create(size_reverse_recv,n,"comm:size");
  memory->create(firstrecv,n,"comm:firstrecv");
  memory->create(pbc_flag,n,"comm:pbc_flag");
  memory->create(pbc,n,6,"comm:pbc");
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void CommTiled::free_swap()
{
  memory->destroy(sendnum);
  memory->destroy(recvnum);
  memory->destroy(sendproc);
  memory->destroy(recvproc);
  memory->destroy(size_forward_recv);
  memory->destroy(size_reverse_send);
  memory->destroy(size_reverse_recv);
  memory->destroy(firstrecv);
  memory->destroy(pbc_flag);
  memory->destroy(pbc);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint CommTiled::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
