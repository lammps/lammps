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
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "stdlib.h"
#include "string.h"
#include "fix_ave_atom.h"
#include "fix_reaxc_bonds.h"
#include "atom.h"
#include "update.h"
#include "pair_reax_c.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCBonds::FixReaxCBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix reax/c/bonds command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  nmax = nint(atom->nlocal*1.05);
  ntypes = atom->ntypes;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  global_freq = nfreq = atoi(arg[5]);
  
  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix reax/c/bonds command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix reax/c/bonds command");

  if (me == 0) {
    fp = fopen(arg[6],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/c/bonds file %s",arg[6]);
      error->one(FLERR,str);
    }
  }

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for fix reax/c bonds");

  allocate();

  nvalid = nextvalid();
}

/* ---------------------------------------------------------------------- */

FixReaxCBonds::~FixReaxCBonds()
{
  MPI_Comm_rank(world,&me);

  memory->destroy(sbo);
  memory->destroy(nlp);
  memory->destroy(avq);
  memory->destroy(numneigh);
  memory->destroy(neighid);
  memory->destroy(tmpid);
  memory->destroy(abo);
  memory->destroy(tmpabo);

  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::init()
{
  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/bonds without "
		  "pair_style reax/c");

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::end_of_step()
{
  // skip if not step that requires calculating bonds
  bigint ntimestep = update->ntimestep;

  if (ntimestep != nvalid) return;

  Output_ReaxC_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::Output_ReaxC_Bonds(bigint ntimestep, FILE *fp)

{
  int i, j, k, itype, jtype, iatom, itag, jtag;
  int b, nbuf, nbuf_local, inode;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);

  repeat = nrepeat;
  // zero out average BO for next Nfreq
  if (irepeat == 0)
    for (i = 0; i < nmax; i++) {
      sbo[i] = nlp[i] = avq[i] = 0.0;
      for (j = 0; j < MAXBOND; j++) {
	tmpid[i][j] = 0;
        tmpabo[i][j] = 0.0;
      }
    }

  // Accumulate bonding information for each nvalid
  GatherBond( system, lists);

  // done if irepeat < nrepeat, else reset irepeat and nvalid
  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }
  irepeat = numbonds = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;

  // Determine actual bonding based on averaged bond order
  FindBond( system, lists, numbonds);

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);
  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memory->create(buf,nbuf,"reax/c/bonds:buf");
  for (i = 0; i < nbuf; i ++) buf[i] = 0.0;

  // Pass information to buffer
  PassBuffer( system, buf, nbuf_local);

  // Receive information from buffer for output
  RecvBuffer( system, buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memory->destroy(buf);

}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::GatherBond( reax_system *system, reax_list *lists)
{
  int *ilist, i, ii, inum;
  int j, pj, nj, jtag, jtype;
  double bo_tmp,bo_cut;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  bond_data *bo_ij;
  bo_cut = 0.10; //reaxc->control->bg_cut;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;

    for( pj = Start_Index(i, reaxc->lists); pj < End_Index(i, reaxc->lists); ++pj ) {
      bo_ij = &( reaxc->lists->select.bond_list[pj] );
      j = bo_ij->nbr;
      jtag = reaxc->system->my_atoms[j].orig_id;
      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp > bo_cut) {
	here:;
	if (jtag != tmpid[i][nj] && tmpid[i][nj] != 0) {
	  nj ++;
	  if (nj > MAXBOND) error->all(FLERR,"Increase MAXBOND value");
	  goto here;
	}
	tmpid[i][nj] = jtag;
	tmpabo[i][nj] += bo_tmp;
        nj ++;
      }

    }
    sbo[i] += reaxc->workspace->total_bond_order[i];
    nlp[i] += reaxc->workspace->nlp[i];
    avq[i] += atom->q[i];
  }

}
/* ---------------------------------------------------------------------- */

void FixReaxCBonds::FindBond( reax_system *system, reax_list *lists,
		int &numbonds)
{
  int *ilist, i, ii, inum;
  int j, pj, nj, jtag, jtype;
  double bo_tmp,bo_cut;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  bo_cut = reaxc->control->bg_cut;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    sbo[i] /= repeat;
    nlp[i] /= repeat;
    avq[i] /= repeat;
    numneigh[i] = 0;
    for (j = 0; j < MAXBOND; j++) {
      tmpabo[i][j] /= repeat;
      neighid[i][j] = 0;
      abo[i][j] = 0.0;
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;

    for (j = 0; j < MAXBOND; j++){
      if (tmpabo[i][j] > bo_cut) {
	neighid[i][nj] = tmpid[i][j];
	abo[i][nj] = tmpabo[i][j];
	nj ++;
      }
    }
    numneigh[i] = nj;
    if (nj > numbonds) numbonds = nj;
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::PassBuffer( reax_system *system, double *buf, 
		int &nbuf_local)
{
  int i, j, k, jtag, numbonds;
  int nlocal = atom->nlocal;

  j = 2;
  buf[0] = nlocal;
  for (i = 0; i < nlocal; i++) {
    buf[j-1] = atom->tag[i];
    buf[j+0] = atom->type[i];
    buf[j+1] = sbo[i];
    buf[j+2] = nlp[i];
    buf[j+3] = avq[i];
    buf[j+4] = numneigh[i];
    numbonds = nint(buf[j+4]);

    for (k = 5; k < 5+numbonds; k++) {  
      buf[j+k] = neighid[i][k-5];
    }
    j += (5+numbonds);

    if (atom->molecule == NULL ) buf[j] = 0.0;
    else buf[j] = atom->molecule[i];
    j ++;

    for (k = 0; k < numbonds; k++) {	
      buf[j+k] = abo[i][k];
    }
    j += (1+numbonds);
  }
  nbuf_local = j - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::RecvBuffer( reax_system *system, double *buf, 
		int nbuf, int nbuf_local, int natoms, int maxnum)
{
  int i, j, k, l, itype, jtype, itag, jtag;
  int inode, nlocal_tmp, numbonds, molid;
  int nlocal = atom->nlocal;
  int ntimestep = update->ntimestep;
  double sbotmp, nlptmp, avqtmp, abotmp;

  double cutof3 = reaxc->control->bg_cut;
  MPI_Request irequest, irequest2;
  MPI_Status istatus;

  if (me == 0 ){ 
    fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",ntimestep);
    fprintf(fp,"# \n");
    fprintf(fp,"# Number of particles %d \n",natoms);
    fprintf(fp,"# \n");
    fprintf(fp,"# Max number of bonds per atom %d with "
	    "coarse bond order cutoff %5.3f \n",maxnum,cutof3);
    fprintf(fp,"# Particle connection table and bond orders \n");
    fprintf(fp,"# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q \n");
  }

  j = 2;
  if (me == 0) {
    for (inode = 0; inode < nprocs; inode ++) {
      if (inode == 0) {
	nlocal_tmp = nlocal;
      } else {
	MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Wait(&irequest,&istatus);
	nlocal_tmp = nint(buf[0]);
      }
      j = 2;
      for (i = 0; i < nlocal_tmp; i ++) {
	itag = nint(buf[j-1]);
	itype = nint(buf[j+0]);
	sbotmp = buf[j+1];
	nlptmp = buf[j+2];
	avqtmp = buf[j+3];
	numbonds = nint(buf[j+4]);

	fprintf(fp," %d %d %d",itag,itype,numbonds);

	for (k = 5; k < 5+numbonds; k++) {
          jtag = nint(buf[j+k]);
	  fprintf(fp," %d",jtag);
	}
	j += (5+numbonds);

	fprintf(fp," %d",nint(buf[j]));
	j ++;

	for (k = 0; k < numbonds; k++) {
	  abotmp = buf[j+k];
	  fprintf(fp,"%14.3f",abotmp);
	}
	j += (1+numbonds);
	fprintf(fp,"%14.3f%14.3f%14.3f\n",sbotmp,nlptmp,avqtmp);
      }
    }
  } else {
    MPI_Isend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world,&irequest2);
    MPI_Wait(&irequest2,&istatus);
  }
  if(me ==0) fprintf(fp,"# \n");

}

/* ---------------------------------------------------------------------- */

int FixReaxCBonds::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixReaxCBonds::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::allocate()
{

  irepeat = 0;

  sbo = NULL;
  nlp = NULL;
  avq = NULL;
  numneigh = NULL;
  neighid = NULL;
  tmpid = NULL;
  abo = NULL;
  tmpabo = NULL;

  memory->create(sbo,nmax,"reax/c/bonds:sbo");
  memory->create(nlp,nmax,"reax/c/bonds:nlp");
  memory->create(avq,nmax,"reax/c/bonds:avq");
  memory->create(numneigh,nmax,"reax/c/bonds:numneigh");
  memory->create(abo,nmax,MAXBOND,"reax/c/bonds:abo");
  memory->create(neighid,nmax,MAXBOND,"reax/c/bonds:neighid");
  memory->create(tmpabo,nmax,MAXBOND,"reax/c/bonds:tmpabo");
  memory->create(tmpid,nmax,MAXBOND,"reax/c/bonds:tmpid");

  irepeat = 0;

}

/* ---------------------------------------------------------------------- */

double FixReaxCBonds::memory_usage()
{
  double bytes;

  bytes += 3.0*nmax*sizeof(double);
  bytes += nmax*sizeof(int);
  bytes += 2.0*nmax*MAXBOND*sizeof(double);
  bytes += 2.0*nmax*MAXBOND*sizeof(int);

  return bytes;
}

