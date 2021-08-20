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
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include "fix_reaxff_bonds.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace ReaxFF;

/* ---------------------------------------------------------------------- */

FixReaxFFBonds::FixReaxFFBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix reaxff/bonds command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;
  nmax = atom->nmax;

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);

  if (nevery <= 0)
    error->all(FLERR,"Illegal fix reaxff/bonds command");

  if (me == 0) {
    char *suffix = strrchr(arg[4],'.');
    if (suffix && strcmp(suffix,".gz") == 0) {
#ifdef LAMMPS_GZIP
      auto gzip = fmt::format("gzip -6 > {}",arg[4]);
#ifdef _WIN32
      fp = _popen(gzip.c_str(),"wb");
#else
      fp = popen(gzip.c_str(),"w");
#endif
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else fp = fopen(arg[4],"w");

    if (!fp)
      error->one(FLERR,fmt::format("Cannot open fix reaxff/bonds file {}: "
                                   "{}",arg[4],utils::getsyserror()));
  }

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for fix reaxff bonds");

  abo = nullptr;
  neighid = nullptr;
  numneigh = nullptr;

  allocate();
}

/* ---------------------------------------------------------------------- */

FixReaxFFBonds::~FixReaxFFBonds()
{
  MPI_Comm_rank(world,&me);

  destroy();

  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::init()
{
  reaxff = (PairReaxFF *) force->pair_match("^reax..",0);
  if (reaxff == nullptr) error->all(FLERR,"Cannot use fix reaxff/bonds without "
                                "pair_style reaxff, reaxff/kk, or reaxff/omp");
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::end_of_step()
{
  Output_ReaxFF_Bonds();
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::Output_ReaxFF_Bonds()

{
  int i, j;
  int nbuf, nbuf_local;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);

  if (atom->nmax > nmax) {
    destroy();
    nmax = atom->nmax;
    allocate();
  }

  for (i = 0; i < nmax; i++) {
    numneigh[i] = 0;
    for (j = 0; j < MAXREAXBOND; j++) {
      neighid[i][j] = 0;
      abo[i][j] = 0.0;
    }
  }
  numbonds = FindBond();

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);

  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memory->create(buf,nbuf,"reaxff/bonds:buf");
  for (i = 0; i < nbuf; i ++) buf[i] = 0.0;

  // Pass information to buffer
  PassBuffer(buf, nbuf_local);

  // Receive information from buffer for output
  RecvBuffer(buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memory->destroy(buf);

}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::FindBond()
{
  int *ilist, i, ii, inum;
  int j, pj, nj;
  tagint jtag;
  double bo_tmp,bo_cut;

  inum = reaxff->list->inum;
  ilist = reaxff->list->ilist;
  bond_data *bo_ij;
  bo_cut = reaxff->api->control->bg_cut;

  tagint *tag = atom->tag;
  int numbonds = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;

    for (pj = Start_Index(i, reaxff->api->lists); pj < End_Index(i, reaxff->api->lists); ++pj) {
      bo_ij = &(reaxff->api->lists->select.bond_list[pj]);
      j = bo_ij->nbr;
      jtag = tag[j];
      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp > bo_cut) {
        neighid[i][nj] = jtag;
        abo[i][nj] = bo_tmp;
        nj ++;
      }
    }
    numneigh[i] = nj;
    if (nj > numbonds) numbonds = nj;
  }
  return numbonds;
}
/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::PassBuffer(double *buf, int &nbuf_local)
{
  int i, j, k, numbonds;
  int nlocal = atom->nlocal;

  j = 2;
  buf[0] = nlocal;
  for (i = 0; i < nlocal; i++) {
    buf[j-1] = atom->tag[i];
    buf[j+0] = atom->type[i];
    buf[j+1] = reaxff->api->workspace->total_bond_order[i];
    buf[j+2] = reaxff->api->workspace->nlp[i];
    buf[j+3] = atom->q[i];
    buf[j+4] = numneigh[i];
    numbonds = nint(buf[j+4]);

    for (k = 5; k < 5+numbonds; k++) {
      buf[j+k] = neighid[i][k-5];
    }
    j += (5+numbonds);

    if (atom->molecule == nullptr) buf[j] = 0.0;
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

void FixReaxFFBonds::RecvBuffer(double *buf, int nbuf, int nbuf_local,
                               int natoms, int maxnum)
{
  int i, j, k, itype;
  int inode, nlocal_tmp, numbonds;
  tagint itag,jtag;
  int nlocal = atom->nlocal;
  bigint ntimestep = update->ntimestep;
  double sbotmp, nlptmp, avqtmp, abotmp;

  double cutof3 = reaxff->api->control->bg_cut;
  MPI_Request irequest, irequest2;

  if (me == 0) {
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
        MPI_Wait(&irequest,MPI_STATUS_IGNORE);
        nlocal_tmp = nint(buf[0]);
      }
      j = 2;
      for (i = 0; i < nlocal_tmp; i ++) {
        itag = static_cast<tagint> (buf[j-1]);
        itype = nint(buf[j+0]);
        sbotmp = buf[j+1];
        nlptmp = buf[j+2];
        avqtmp = buf[j+3];
        numbonds = nint(buf[j+4]);

        fprintf(fp," " TAGINT_FORMAT " %d %d",itag,itype,numbonds);

        for (k = 5; k < 5+numbonds; k++) {
          jtag = static_cast<tagint> (buf[j+k]);
          fprintf(fp," " TAGINT_FORMAT,jtag);
        }
        j += (5+numbonds);

        fprintf(fp," " TAGINT_FORMAT,static_cast<tagint> (buf[j]));
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
    MPI_Wait(&irequest2,MPI_STATUS_IGNORE);
  }
  if (me ==0) fprintf(fp,"# \n");

}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::destroy()
{
  memory->destroy(abo);
  memory->destroy(neighid);
  memory->destroy(numneigh);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::allocate()
{
  memory->create(abo,nmax,MAXREAXBOND,"reaxff/bonds:abo");
  memory->create(neighid,nmax,MAXREAXBOND,"reaxff/bonds:neighid");
  memory->create(numneigh,nmax,"reaxff/bonds:numneigh");
}

/* ---------------------------------------------------------------------- */

double FixReaxFFBonds::memory_usage()
{
  double bytes;

  bytes = 3.0*nmax*sizeof(double);
  bytes += (double)nmax*sizeof(int);
  bytes += (double)1.0*nmax*MAXREAXBOND*sizeof(double);
  bytes += (double)1.0*nmax*MAXREAXBOND*sizeof(int);

  return bytes;
}
