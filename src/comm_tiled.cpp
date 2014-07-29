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

// NOTE: change this to 16 after debugged

#define DELTA_PROCS 1

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
  memory->destroy(overlap);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
  deallocate_swap(nswap);
}

/* ----------------------------------------------------------------------
   initialize comm buffers and other data structs local to CommTiled
   NOTE: if this is identical to CommBrick, put it into Comm ??
------------------------------------------------------------------------- */

void CommTiled::init_buffers()
{
  memory->create(overlap,nprocs,"comm:overlap");

  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+bufextra,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");

  nswap = 2 * domain->dimension;
  allocate_swap(nswap);
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
------------------------------------------------------------------------- */

void CommTiled::setup()
{
  int i;

  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  double *prd = domain->prd;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  // set function pointers

  if (layout != LAYOUT_TILED) {
    box_drop = &CommTiled::box_drop_brick;
    box_other = &CommTiled::box_other_brick;
  } else {
    box_drop = &CommTiled::box_drop_tiled;
    box_other = &CommTiled::box_other_tiled;
  }

  // if RCB decomp has just changed, gather needed global RCB info

  if (rcbnew) {
    rcbnew = 0;
    memcpy(&rcbinfo[me].mysplit[0][0],&mysplit[0][0],6*sizeof(double));
    rcbinfo[me].cut = rcbcut;
    rcbinfo[me].dim = rcbcutdim;
    MPI_Allgather(&rcbinfo[me],sizeof(RCBinfo),MPI_CHAR,
                  rcbinfo,sizeof(RCBinfo),MPI_CHAR,world);
  }

  // check that cutoff < any periodic box length

  double cut = MAX(neighbor->cutneighmax,cutghostuser);
  cutghost[0] = cutghost[1] = cutghost[2] = cut;
  
  if ((periodicity[0] && cut > prd[0]) ||
      (periodicity[1] && cut > prd[1]) ||
      (dimension == 3 && periodicity[2] && cut > prd[2]))
    error->all(FLERR,"Communication cutoff for comm_style tiled "
               "cannot exceed periodic box length");

  // loop over 6 swap directions
  // determine which procs I will send to and receive from in each swap
  // done by intersecting ghost box with all proc sub-boxes it overlaps
  // sets nsendproc, nrecvproc, sendproc, recvproc
  // sets sendother, sendself, pbc_flag, pbc, sendbox

  int noverlap,noverlap1,indexme;
  double lo1[3],hi1[3],lo2[3],hi2[3];
  int one,two;

  nswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    for (int iswap = 0; iswap < 2; iswap++) {

      // one = first ghost box in same periodic image
      // two = second ghost box wrapped across periodic boundary
      // either may not exist

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

      // noverlap = # of overlaps of box1/2 with procs via box_drop()
      // overlap = list of overlapping procs
      // if overlap with self, indexme = index of me in list

      indexme = -1;
      noverlap = 0;
      if (one) (this->*box_drop)(idim,lo1,hi1,noverlap,indexme);
      noverlap1 = noverlap;
      if (two) (this->*box_drop)(idim,lo2,hi2,noverlap,indexme);

      // if self is in overlap list, move it to end of list

      if (indexme >= 0) {
        int tmp = overlap[noverlap-1];
        overlap[noverlap-1] = overlap[indexme];
        overlap[indexme] = tmp;
      }

      // reallocate 2nd dimensions of all send/recv arrays, based on noverlap
      // # of sends of this swap = # of recvs of nswap +/- 1

      if (noverlap > nprocmax[iswap]) {
        int oldmax = nprocmax[iswap];
        while (nprocmax[iswap] < noverlap) nprocmax[iswap] += DELTA_PROCS;
        grow_swap_send(iswap,nprocmax[iswap],oldmax);
        if (iswap == 0) grow_swap_recv(iswap+1,nprocmax[iswap]);
        else grow_swap_recv(iswap-1,nprocmax[iswap]);
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
        
        (this->*box_other)(idim,iswap,overlap[i],oboxlo,oboxhi);
        
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
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

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
  int i,m,n,irecv,nlocal,nlast,nsend,nrecv,ncount,rmaxswap;
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
    //       and for ngroup when bordergroup is set

    x = atom->x;

    for (m = 0; m < nsendproc[iswap]; m++) {
      bbox = sendbox[iswap][m];
      xlo = bbox[0]; ylo = bbox[1]; zlo = bbox[2];
      xhi = bbox[3]; yhi = bbox[4]; zhi = bbox[5];

      nlocal = atom->nlocal;
      if (iswap < 2) nlast = atom->nlocal;
      else nlast = atom->nlocal + atom->nghost;

      ncount = 0;
      for (i = 0; i < nlocal; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (ncount == maxsendlist[iswap][m]) grow_list(iswap,i,ncount);
          sendlist[iswap][m][ncount++] = i;
        }
      for (i = atom->nlocal; i < nlast; i++)
        if (x[i][0] >= xlo && x[i][0] <= xhi &&
            x[i][1] >= ylo && x[i][1] <= yhi &&
            x[i][2] >= zlo && x[i][2] <= zhi) {
          if (ncount == maxsendlist[iswap][m]) grow_list(iswap,i,ncount);
          sendlist[iswap][m][ncount++] = i;
        }
      sendnum[iswap][m] = ncount;
      smax = MAX(smax,ncount);
    }

    // send sendnum counts to procs who recv from me except self
    // copy data to self if sendself is set

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (sendother[iswap]) {
      for (m = 0; m < nrecv; m++)
        MPI_Irecv(&recvnum[iswap][m],1,MPI_INT,
                  recvproc[iswap][m],0,world,&requests[m]);
      for (m = 0; m < nsend; m++)
        MPI_Send(&sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],0,world);
    }
    if (sendself[iswap]) recvnum[iswap][nrecv] = sendnum[iswap][nsend];
    if (sendother[iswap]) MPI_Waitall(nrecv,requests,statuses);

    // setup other per swap/proc values from sendnum and recvnum

    rmaxswap = 0;
    for (m = 0; m < nrecvproc[iswap]; m++) {
      rmaxswap += recvnum[iswap][m];
      size_forward_recv[iswap][m] = recvnum[iswap][m]*size_forward;
      size_reverse_send[iswap][m] = recvnum[iswap][m]*size_reverse;
      size_reverse_recv[iswap][m] = sendnum[iswap][m]*size_reverse;
      if (m == 0) {
        firstrecv[iswap][0] = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][m] = firstrecv[iswap][m-1] + recvnum[iswap][m-1];
        forward_recv_offset[iswap][m] = 
          forward_recv_offset[iswap][m-1] + recvnum[iswap][m-1];
      }
    }
    rmax = MAX(rmax,rmaxswap);

    // insure send/recv buffers are large enough for border comm

    if (smax*size_border > maxsend) grow_send(smax*size_border,0);
    if (rmax*size_border > maxrecv) grow_recv(rmax*size_border);

    // swap atoms with other procs using pack_border(), unpack_border()

    if (ghost_velocity) {
      if (sendother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][m]],
                    recvnum[iswap][m]*size_border,
                    MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        for (m = 0; m < nsend; m++) {
          n = avec->pack_border_vel(sendnum[iswap][m],sendlist[iswap][m],
                                    buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
          MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                   sendproc[iswap][m],0,world);
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
        for (m = 0; m < nrecv; m++) {
          MPI_Waitany(nrecv,requests,&irecv,&status);
          avec->unpack_border(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                              &buf_recv[forward_recv_offset[iswap][irecv]]);
        }
      }

    } else {
      if (sendother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[forward_recv_offset[iswap][m]],
                    recvnum[iswap][m]*size_border,
                    MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
        for (m = 0; m < nsend; m++) {
          n = avec->pack_border(sendnum[iswap][m],sendlist[iswap][m],
                                buf_send,pbc_flag[iswap][m],pbc[iswap][m]);
          MPI_Send(buf_send,n*size_border,MPI_DOUBLE,
                   sendproc[iswap][m],0,world);
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
        for (m = 0; m < nrecv; m++) {
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

// NOTE: remaining forward/reverse methods still need to be updated

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

void CommTiled::box_drop_brick(int idim, double *lo, double *hi, 
                               int &noverlap, int &indexme)
{
  // NOTE: this is not triclinic compatible

  double *prd = domain->prd;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  int index,dir;
  if (hi[idim] == sublo[idim]) {
    index = myloc[idim] - 1;
    dir = -1;
  } else if (lo[idim] == subhi[idim]) {
    index = myloc[idim] + 1;
    dir = 1;
  } else if (hi[idim] == boxhi[idim]) {
    index = procgrid[idim] - 1;
    dir = -1;
  } else if (lo[idim] == boxlo[idim]) {
    index = 0;
    dir = 1;
  }

  int other1,other2,proc;
  double lower,upper;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  while (1) {
    lower = boxlo[idim] + prd[idim]*split[index];
    if (index < procgrid[idim]-1) 
      upper = boxlo[idim] + prd[idim]*split[index+1];
    else upper = boxhi[idim];
    if (lower >= hi[idim] || upper <= lo[idim]) break;

    if (idim == 0) proc = grid2proc[index][other1][other2];
    else if (idim == 1) proc = grid2proc[other1][index][other2];
    else proc = grid2proc[other1][other2][idim];

    if (proc == me) indexme = noverlap;
    overlap[noverlap++] = proc;
    index += dir;
    if (index < 0 || index >= procgrid[idim]) break;
  }
}

/* ----------------------------------------------------------------------
   determine overlap list of Noverlap procs the lo/hi box overlaps
   overlap = non-zero area in common between box and proc sub-domain
   recursive method for traversing an RCB tree of cuts
   no need to split lo/hi box as recurse b/c OK if box extends outside RCB box
------------------------------------------------------------------------- */

void CommTiled::box_drop_tiled(int idim, double *lo, double *hi, 
                               int &noverlap, int &indexme)
{
  box_drop_tiled_recurse(lo,hi,0,nprocs-1,noverlap,indexme);
}

void CommTiled::box_drop_tiled_recurse(double *lo, double *hi, 
                                       int proclower, int procupper,
                                       int &noverlap, int &indexme)
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
  double cut = rcbinfo[procmid].cut;
  int idim = rcbinfo[procmid].dim;
  
  if (lo[idim] < cut) 
    box_drop_tiled_recurse(lo,hi,proclower,procmid-1,noverlap,indexme);
  if (hi[idim] > cut)
    box_drop_tiled_recurse(lo,hi,procmid,procupper,noverlap,indexme);
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
------------------------------------------------------------------------- */

void CommTiled::box_other_brick(int idim, int iswap,
                                int proc, double *lo, double *hi)
{
  double *prd = domain->prd;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  lo[0] = sublo[0]; lo[1] = sublo[1]; lo[2] = sublo[2]; 
  hi[0] = subhi[0]; hi[1] = subhi[1]; hi[2] = subhi[2]; 

  int other1,other2,oproc;
  double *split;

  if (idim == 0) {
    other1 = myloc[1]; other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0]; other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0]; other2 = myloc[1];
    split = zsplit;
  }

  int dir = -1;
  if (iswap) dir = 1;
  int index = myloc[idim];
  int n = procgrid[idim];

  for (int i = 0; i < n; i++) {
    index += dir;
    if (index < 0) index = n-1;
    else if (index >= n) index = 0;

    if (idim == 0) oproc = grid2proc[index][other1][other2];
    else if (idim == 1) oproc = grid2proc[other1][index][other2];
    else oproc = grid2proc[other1][other2][idim];

    if (proc == oproc) {
      lo[idim] = boxlo[idim] + prd[idim]*split[index];
      if (split[index+1] < 1.0) 
        hi[idim] = boxlo[idim] + prd[idim]*split[index+1];
      else hi[idim] = boxhi[idim];
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   return other box owned by proc as lo/hi corner pts
------------------------------------------------------------------------- */

void CommTiled::box_other_tiled(int idim, int iswap,
                                int proc, double *lo, double *hi)
{
  double *prd = domain->prd;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double (*split)[2] = rcbinfo[proc].mysplit;

  lo[0] = boxlo[0] + prd[0]*split[0][0];
  if (split[0][1] < 1.0) hi[0] = boxlo[0] + prd[0]*split[0][1];
  else hi[0] = boxhi[0];

  lo[1] = boxlo[1] + prd[1]*split[1][0];
  if (split[1][1] < 1.0) hi[1] = boxlo[1] + prd[1]*split[1][1];
  else hi[1] = boxhi[1];

  lo[2] = boxlo[2] + prd[2]*split[2][0];
  if (split[2][1] < 1.0) hi[2] = boxlo[2] + prd[2]*split[2][1];
  else hi[2] = boxhi[2];
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
               "comm:sendlist[i][j]");
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void CommTiled::allocate_swap(int n)
{
  nsendproc = new int[n];
  nrecvproc = new int[n];
  sendother = new int[n];
  sendself = new int[n];
  nprocmax = new int[n];

  sendproc = new int*[n];
  recvproc = new int*[n];
  sendnum = new int*[n];
  recvnum = new int*[n];
  size_forward_recv = new int*[n];
  firstrecv = new int*[n];
  size_reverse_send = new int*[n];
  size_reverse_recv = new int*[n];
  forward_recv_offset = new int*[n];
  reverse_recv_offset = new int*[n];

  pbc_flag = new int*[n];
  pbc = new int**[n];
  sendbox = new double**[n];
  maxsendlist = new int*[n];
  sendlist = new int**[n];

  for (int i = 0; i < n; i++) {
    sendproc[i] = recvproc[i] = NULL;
    sendnum[i] = recvnum[i] = NULL;
    size_forward_recv[i] = firstrecv[i] = NULL;
    size_reverse_send[i] = size_reverse_recv[i] = NULL;
    forward_recv_offset[i] = reverse_recv_offset[i] = NULL;

    pbc_flag[i] = NULL;
    pbc[i] = NULL;
    sendbox[i] = NULL;
    maxsendlist[i] = NULL;
    sendlist[i] = NULL;
  }

  maxreqstat = 0;
  requests = NULL;
  statuses = NULL;

  for (int i = 0; i < n; i++) {
    nprocmax[i] = DELTA_PROCS;
    grow_swap_send(i,DELTA_PROCS,0);
    grow_swap_recv(i,DELTA_PROCS);
  }
}

/* ----------------------------------------------------------------------
   grow info for swap I, to allow for N procs to communicate with
   ditto for complementary recv for swap I+1 or I-1, as invoked by caller
------------------------------------------------------------------------- */

void CommTiled::grow_swap_send(int i, int n, int nold)
{
  delete [] sendproc[i];
  sendproc[i] = new int[n];
  delete [] sendnum[i];
  sendnum[i] = new int[n];

  delete [] size_reverse_recv[i];
  size_reverse_recv[i] = new int[n];
  delete [] reverse_recv_offset[i];
  reverse_recv_offset[i] = new int[n];

  delete [] pbc_flag[i];
  pbc_flag[i] = new int[n];
  memory->destroy(pbc[i]);
  memory->create(pbc[i],n,6,"comm:pbc_flag");
  memory->destroy(sendbox[i]);
  memory->create(sendbox[i],n,6,"comm:sendbox");

  delete [] maxsendlist[i];
  maxsendlist[i] = new int[n];

  for (int j = 0; j < nold; j++) memory->destroy(sendlist[i][j]);
  delete [] sendlist[i];
  sendlist[i] = new int*[n];
  for (int j = 0; j < n; j++) {
    maxsendlist[i][j] = BUFMIN;
    memory->create(sendlist[i][j],BUFMIN,"comm:sendlist[i][j]");
  }
}

void CommTiled::grow_swap_recv(int i, int n)
{
  delete [] recvproc[i];
  recvproc[i] = new int[n];
  delete [] recvnum[i];
  recvnum[i] = new int[n];

  delete [] size_forward_recv[i];
  size_forward_recv[i] = new int[n];
  delete [] firstrecv[i];
  firstrecv[i] = new int[n];
  delete [] forward_recv_offset[i];
  forward_recv_offset[i] = new int[n];

  delete [] size_reverse_send[i];
  size_reverse_send[i] = new int[n];

  if (n > maxreqstat) {
    maxreqstat = n;
    delete [] requests;
    delete [] statuses;
    requests = new MPI_Request[n];
    statuses = new MPI_Status[n];
  }
}

/* ----------------------------------------------------------------------
   deallocate swap info
------------------------------------------------------------------------- */

void CommTiled::deallocate_swap(int n)
{
  delete [] nsendproc;
  delete [] nrecvproc;
  delete [] sendother;
  delete [] sendself;

  for (int i = 0; i < n; i++) {
    delete [] sendproc[i];
    delete [] recvproc[i];
    delete [] sendnum[i];
    delete [] recvnum[i];
    delete [] size_forward_recv[i];
    delete [] firstrecv[i];
    delete [] size_reverse_send[i];
    delete [] size_reverse_recv[i];
    delete [] forward_recv_offset[i];
    delete [] reverse_recv_offset[i];

    delete [] pbc_flag[i];
    memory->destroy(pbc[i]);
    memory->destroy(sendbox[i]);
    delete [] maxsendlist[i];

    for (int j = 0; j < nprocmax[i]; j++) memory->destroy(sendlist[i][j]);
    delete [] sendlist[i];
  }

  delete [] sendproc;
  delete [] recvproc;
  delete [] sendnum;
  delete [] recvnum;
  delete [] size_forward_recv;
  delete [] firstrecv;
  delete [] size_reverse_send;
  delete [] size_reverse_recv;
  delete [] forward_recv_offset;
  delete [] reverse_recv_offset;

  delete [] pbc_flag;
  delete [] pbc;
  delete [] sendbox;
  delete [] maxsendlist;
  delete [] sendlist;

  delete [] requests;
  delete [] statuses;

  delete [] nprocmax;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint CommTiled::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}
