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
   Contributing author: Aidan Thompson (Sandia)
------------------------------------------------------------------------- */

#ifdef LAMMPS_BIGBIG
#error LAMMPS_BIGBIG not supported by this file
#endif

#include "stdlib.h"
#include "string.h"
#include "fix_reax_bonds.h"
#include "pair_reax_fortran.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxBonds::FixReaxBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix reax/bonds command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 1) error->all(FLERR,"Illegal fix reax/bonds command");

  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/bonds file %s",arg[4]);
      error->one(FLERR,str);
    }
  }
}

/* ---------------------------------------------------------------------- */

FixReaxBonds::~FixReaxBonds()
{
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   perform initial write
------------------------------------------------------------------------- */

void FixReaxBonds::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxBonds::init()
{
  // insure ReaxFF is defined

  if (force->pair_match("reax",1) == NULL)
    error->all(FLERR,"Cannot use fix reax/bonds without pair_style reax");
}

/* ---------------------------------------------------------------------- */

void FixReaxBonds::end_of_step()
{
  OutputReaxBonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxBonds::OutputReaxBonds(bigint ntimestep, FILE *fp)
{
  int nparticles,nparticles_tot,nbuf,nbuf_local,most,j;
  int ii,jn,mbond,numbonds,nsbmax,nsbmax_most;
  int nprocs,nlocal_tmp,itmp;
  int k,kk,jj,jbufknum;
  double cutof3;
  double *buf;
  MPI_Request irequest;
  MPI_Status istatus;

  MPI_Comm_size(world,&nprocs);

  nparticles = atom->nlocal;
  nparticles_tot = static_cast<int> (atom->natoms);

  jn = ReaxParams::nat;
  mbond = ReaxParams::mbond;
  FORTRAN(getnsbmax,GETNSBMAX)(&nsbmax);
  FORTRAN(getcutof3,GETCUTOF3)(&cutof3);
  MPI_Allreduce(&nparticles,&most,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nsbmax,&nsbmax_most,1,MPI_INT,MPI_MAX,world);

  if (me == 0) {
    fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",ntimestep);
    fprintf(fp,"# \n");
    fprintf(fp,"# Number of particles %d \n",nparticles_tot);
    fprintf(fp,"# \n");
    fprintf(fp,"# Max number of bonds per atom %d with "
            "coarse bond order cutoff %5.3f \n",
            nsbmax_most,cutof3);
    fprintf(fp,"# Particle connection table and bond orders \n");
    fprintf(fp,"# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q \n");
  }

  // allocate a temporary buffer for the snapshot info
  // big enough for largest number of atoms on any one proc
  // nbuf_local = size of local buffer for table of atom bonds

  nbuf = 1+(2*nsbmax_most+7)*most;
  memory->create(buf,nbuf,"reax/bonds:buf");

  j = 0;
  buf[j++] = nparticles;
  for (int iparticle=0;iparticle<nparticles;iparticle++) {
    buf[j++] = atom->tag[iparticle];                   //atom tag
    buf[j++] = FORTRAN(cbkia,CBKIA).iag[iparticle];    //atom type
    jbufknum = j++;
    numbonds = FORTRAN(cbkia,CBKIA).iag[iparticle+jn];

    // connection table based on coarse bond order cutoff (> cutof3)

    kk = 0;
    for (k=0;k<numbonds;k++) {
      ii = FORTRAN(cbknubon2,CBKNUBON2).nubon1[iparticle+jn*k];
      if (FORTRAN(cbkbo,CBKBO).bo[ii-1] > cutof3) {
        kk++;
        jj = FORTRAN(cbkia,CBKIA).iag[iparticle+jn*(k+2)];
        buf[j++] = FORTRAN(cbkc,CBKC).itag[jj-1];
      }
    }
    buf[jbufknum] = kk; //no.bonds
    buf[j++]=FORTRAN(cbkia,CBKIA).iag[iparticle+jn*(mbond+2)]; //molec.id

    // bond orders (> cutof3)

    kk = 0;
    for (k=0;k<numbonds;k++) {
      ii = FORTRAN(cbknubon2,CBKNUBON2).nubon1[iparticle+jn*k];
      if (FORTRAN(cbkbo,CBKBO).bo[ii-1] > cutof3) {
        kk++;
        buf[j++] = FORTRAN(cbkbo,CBKBO).bo[ii-1];
      }
    }

    // atom bond order (abo), no. of lone pairs (vlp), charge (ch)

    buf[j++] = FORTRAN(cbkabo,CBKABO).abo[iparticle];
    buf[j++] = FORTRAN(cbklonpar,CBKLONPAR).vlp[iparticle];
    buf[j++] = atom->q[iparticle];
  }
  nbuf_local = j-1;

  // node 0 pings each node, receives their buffer, writes to file
  // all other nodes wait for ping, send buffer to node 0

  if (me == 0) {
    for (int inode = 0; inode<nprocs; inode++) {
      j = 0;
      if (inode == 0) {
        nlocal_tmp = nparticles;
        j++;
      } else {
        MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Send(&itmp,0,MPI_INT,inode,0,world);
        MPI_Wait(&irequest,&istatus);
        nlocal_tmp = nint(buf[j++]);
      }

      for (int iparticle=0;iparticle<nlocal_tmp;iparticle++) {

        // print atom tag, atom type, no.bonds

        numbonds = nint(buf[j+2]);
        fprintf(fp," %d %d %d",nint(buf[j]),nint(buf[j+1]),numbonds);
        j += 3;
        if (numbonds > nsbmax_most) {
          char str[128];
          sprintf(str,"Fix reax/bonds numbonds > nsbmax_most");
          error->one(FLERR,str);
        }

        // print connection table

        for (k=0;k<numbonds;k++)
          fprintf(fp," %d",nint(buf[j++]));

        // print molecule id

        fprintf(fp," %d",nint(buf[j++]));

        // print bond orders

        for (k=0;k<numbonds;k++)
          fprintf(fp,"%14.3f",buf[j++]);

        // print sum of bond orders, no. of lone pairs, charge

        fprintf(fp,"%14.3f%14.3f%14.3f\n",buf[j],buf[j+1],buf[j+2]);
        j+=3;
      }
    }

  } else {
    MPI_Recv(&itmp,0,MPI_INT,0,0,world,&istatus);
    MPI_Rsend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world);
  }

  if (me == 0) fprintf(fp,"# \n");

  memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

int FixReaxBonds::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}
