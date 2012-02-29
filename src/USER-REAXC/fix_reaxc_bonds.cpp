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
   Contributing author: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_reaxc_bonds.h"
#include "atom.h"
#include "pair_reax_c.h"
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

FixReaxCBonds::FixReaxCBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix reax/c/bonds command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  if (nevery < 1) error->all(FLERR,"Illegal fix reax/c/bonds command");

  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/c/bonds file %s",arg[4]);
      error->one(FLERR,str);
    }
  }
}

/* ---------------------------------------------------------------------- */

FixReaxCBonds::~FixReaxCBonds()
{
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   perform initial write
------------------------------------------------------------------------- */

void FixReaxCBonds::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::init()
{
  // insure ReaxFF/C is defined
  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/bonds without "
		  "pair_style reax/c");

  // Notify pair_reax_c to calculation bonding information
  reaxc->fixbond_flag = 1;

}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::end_of_step()
{
  OutputReaxCBonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCBonds::OutputReaxCBonds(bigint ntimestep, FILE *fp)

{
  int npart,npart_tot,nbuf,nbuf_local,most,j;
  int ii,jn,mbond,numbonds,nsbmax,nsbmax_most;
  int nprocs,nlocal_tmp;
  double cutof3;
  double *buf;
  MPI_Request irequest, irequest2;
  MPI_Status istatus;
  MPI_Comm_size(world,&nprocs);
 
  npart = atom->nlocal;
  npart_tot = static_cast<int> (atom->natoms);

  mbond = MAX_BOND; 		 	// max bond per atom allowed
  nsbmax = reaxc->system->my_bonds;	// max bond for each atom
  cutof3 = reaxc->control->bg_cut;  	// bond order cutoff for determining bonds

  // get maxval from all nodes
  MPI_Allreduce(&npart,&most,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nsbmax,&nsbmax_most,1,MPI_INT,MPI_MAX,world);
 
  if (me == 0) {
    fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",ntimestep);
    fprintf(fp,"# \n");
    fprintf(fp,"# Number of particles %d \n",npart_tot);
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
  memory->create(buf,nbuf,"reax/c/bonds:buf");
 
  // put bonding information in buffers
  j = 2;
  buf[0] = npart;
  for (int ipart = 0; ipart < npart; ipart++) {
    buf[j-1] = atom->tag[ipart];  		//atom tag
    buf[j+0] = atom->type[ipart];  		//atom type
    buf[j+1] = reaxc->system->my_atoms[ipart].numbonds;  	// no of bonds around i
    numbonds = nint(buf[j+1]);

    int k;
    // connection table based on coarse bond order cutoff (> cutof3)
    for (k=2;k<2+numbonds;k++) 
      buf[j+k] = reaxc->system->my_atoms[ipart].nbr_id[k-1];	// calculated in reaxc_traj

    // molecule id
    if (atom->molecule == NULL )
      buf[j+k] = 0;
    else
      buf[j+k] = atom->molecule[ipart];  	

    j+=(3+numbonds);

    // get bond order values for the connection table
    for (k=0;k<numbonds;k++) 
      buf[j+k] = reaxc->system->my_atoms[ipart].nbr_bo[k+1];

    // sum of bond orders (abo), no. of lone pairs (vlp), charge (ch)
    buf[j+k] = reaxc->workspace->total_bond_order[ipart];
    buf[j+k+1] = reaxc->workspace->nlp[ipart]; 
    buf[j+k+2] = atom->q[ipart];
    j+=(4+numbonds);
  }
  nbuf_local = j-1;
 
  // node 0 pings each node, receives their buffer, writes to file
  // all other nodes wait for ping, send buffer to node 0
 
  if (me == 0) {
    for (int inode = 0; inode<nprocs; inode++) {
      if (inode == 0) {
	nlocal_tmp = npart;
      } else {
	MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Wait(&irequest,&istatus);
	nlocal_tmp = nint(buf[0]);
      }
 
      j = 2;
      for (int ipart=0;ipart<nlocal_tmp;ipart++) {
	// print atom tag, atom type, no.bonds

	fprintf(fp," %d %d %d",nint(buf[j-1]),nint(buf[j+0]),nint(buf[j+1]));
	int k;
	numbonds = nint(buf[j+1]);
	if (numbonds > nsbmax_most) {
	  char str[128];
	  sprintf(str,"Fix reax/c/bonds numbonds > nsbmax_most");
	  error->one(FLERR,str);
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
    MPI_Isend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world,&irequest2);
    MPI_Wait(&irequest2,&istatus);
  }

  if (me == 0) fprintf(fp,"# \n");

  memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

int FixReaxCBonds::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}
