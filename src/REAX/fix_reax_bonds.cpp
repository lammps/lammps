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

/* ---------------------------------------------------------------------- */

FixReaxBonds::FixReaxBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix reax/bonds command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  if (nevery < 1) error->all("Illegal fix reax/bonds command");

  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/bonds file %s",arg[4]);
      error->one(str);
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
    error->all("Cannot use fix reax/bonds without pair_style reax");
}

/* ---------------------------------------------------------------------- */

void FixReaxBonds::end_of_step()
{
  OutputReaxBonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxBonds::OutputReaxBonds(int timestep, FILE *fp) 
{
  int nparticles,nparticles_tot,nbuf,nbuf_local,most,j;
  int ii,jn,mbond,numbonds,nsbmax,nsbmax_most;
  int nprocs,nlocal_tmp,itmp;
  double cutof3;
  double *buf;
  MPI_Request irequest;
  MPI_Status istatus;

  MPI_Comm_size(world,&nprocs);
 
  nparticles = atom->nlocal;
  nparticles_tot = static_cast<int> (atom->natoms);
 
  mbond = ReaxParams::mbond;
  FORTRAN(getnsbmax,GETNSBMAX)(&nsbmax);
  FORTRAN(getcutof3,GETCUTOF3)(&cutof3);
  MPI_Allreduce(&nparticles,&most,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nsbmax,&nsbmax_most,1,MPI_INT,MPI_MAX,world);
 
  if (me == 0) {
    fprintf(fp,"# Timestep %d \n",timestep);
    fprintf(fp,"# \n");
    fprintf(fp,"# Number of particles %d \n",nparticles_tot);
    fprintf(fp,"# \n");
    fprintf(fp,"# Max.number of bonds per atom %d with "
	    "coarse bond order cutoff %5.3f \n",
	    nsbmax_most,cutof3);
    fprintf(fp,"# Particle connection table and bond orders \n");
    fprintf(fp,"# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q \n");
  }
 
  // allocate a temporary buffer for the snapshot info
  // big enough for largest number of atoms on any one proc
  // nbuf_local = size of local buffer for table of atom bonds
 
  nbuf = 1+(2*nsbmax_most+7)*most;
  buf = (double *) memory->smalloc(nbuf*sizeof(double),"reax/bonds:buf");
 
  j = 2;
  jn = ReaxParams::nat;
  buf[0] = nparticles;
  for (int iparticle=0;iparticle<nparticles;iparticle++) {
    buf[j-1] = atom->tag[iparticle];                   //atom tag
    buf[j+0] = FORTRAN(cbkia,CBKIA).iag[iparticle];    //atom type
    buf[j+1] = FORTRAN(cbkia,CBKIA).iag[iparticle+jn]; //no.bonds
    int k;
    numbonds = nint(buf[j+1]);

    // connection table based on coarse bond order cutoff (> cutof3)
    for (k=2;k<2+numbonds;k++) {
      ii = FORTRAN(cbkia,CBKIA).iag[iparticle+jn*k];
      buf[j+k] = FORTRAN(cbkc,CBKC).itag[ii-1];
    }
    buf[j+k]=FORTRAN(cbkia,CBKIA).iag[iparticle+jn*(mbond+2)]; //molec.id
    j+=(3+numbonds);

    // bond orders (> cutof3)
    for (k=0;k<numbonds;k++) {
      ii = FORTRAN(cbknubon2,CBKNUBON2).nubon1[iparticle+jn*k];
      buf[j+k] = FORTRAN(cbkbo,CBKBO).bo[ii-1];
    }
    // sum of bond orders (abo), no. of lone pairs (vlp), charge (ch)
    buf[j+k] = FORTRAN(cbkabo,CBKABO).abo[iparticle];
    buf[j+k+1] = FORTRAN(cbklonpar,CBKLONPAR).vlp[iparticle];
    //    buf[j+k+2] = FORTRAN(cbkch,CBKCH).ch[iparticle];
    buf[j+k+2] = atom->q[iparticle];
    j+=(4+numbonds);
  }
  nbuf_local = j-1;
 
  // node 0 pings each node, receives their buffer, writes to file
  // all other nodes wait for ping, send buffer to node 0
 
  if (me == 0) {
    for (int inode = 0; inode<nprocs; inode++) {
      if (inode == 0) {
	nlocal_tmp = nparticles;
      } else {
	MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Send(&itmp,0,MPI_INT,inode,0,world);
        MPI_Wait(&irequest,&istatus);
	nlocal_tmp = nint(buf[0]);
      }
 
      j = 2;
      for (int iparticle=0;iparticle<nlocal_tmp;iparticle++) {
	// print atom tag, atom type, no.bonds

	fprintf(fp," %d %d %d",nint(buf[j-1]),nint(buf[j+0]),nint(buf[j+1]));
	int k;
	numbonds = nint(buf[j+1]);
	if (numbonds > nsbmax_most) {
	  char str[128];
	  sprintf(str,"Fix reax/bonds numbonds > nsbmax_most");
	  error->one(str);
	}

	// print connection table

	for (k=2;k<2+numbonds;k++)
	  fprintf(fp," %d",nint(buf[j+k]));
	fprintf(fp," %d",nint(buf[j+k]));
	j+=(3+numbonds);

	// print bond orders

	for (k=0;k<numbonds;k++)
	  fprintf(fp,"%14.3f",buf[j+k]);

	// print sum of bond orders, no. of lone pairs, charge

	fprintf(fp,"%14.3f%14.3f%14.3f\n",buf[j+k],buf[j+k+1],buf[j+k+2]);
	j+=(4+numbonds);
      }
    }

  } else {
    MPI_Recv(&itmp,0,MPI_INT,0,0,world,&istatus);
    MPI_Rsend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world);
  }

  if (me == 0) fprintf(fp,"# \n");

  memory->sfree(buf);
}

/* ---------------------------------------------------------------------- */

int FixReaxBonds::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}
