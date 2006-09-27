/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tim Lau (MIT)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eam_fs.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairEAMFS::PairEAMFS()
{
  one_coeff = 1;
}

/* ---------------------------------------------------------------------- */

PairEAMFS::~PairEAMFS() {

  // deallocate array unique to derived class, parent will deallocate the rest

  if (frho) memory->destroy_3d_double_array(rhor_fs);

  // insure derived-class's deallocate is called
  // set ptr to NULL to prevent parent class from calling it's deallocate

  if (frho_0) interpolate_deallocate();
  frho_0 = NULL;
}

/* ---------------------------------------------------------------------- */

void PairEAMFS::compute(int eflag, int vflag)
{
  int i,j,k,m,numneigh,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,p,fforce,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  int *neighs;
  double **f;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(rho);
    memory->sfree(fp);
    nmax = atom->nmax;
    rho = (double *) memory->smalloc(nmax*sizeof(double),"eam:rho");
    fp = (double *) memory->smalloc(nmax*sizeof(double),"eam:fp");
  }

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  if (vflag == 2) f = update->f_pair;
  else f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms
  // FS has type-specific rho functional

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	p = sqrt(rsq)*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);
	rho[i] += 
	  ((rhor_fs_3[jtype][itype][m]*p + rhor_fs_2[jtype][itype][m])*p + 
	   rhor_fs_1[jtype][itype][m])*p + rhor_fs_0[jtype][itype][m];
	if (newton_pair || j < nlocal)
	  rho[j] += 
	    ((rhor_fs_3[itype][jtype][m]*p + rhor_fs_2[itype][jtype][m])*p + 
	     rhor_fs_1[itype][jtype][m])*p + rhor_fs_0[itype][jtype][m];
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // FS is same as standard EAM

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    fp[i] = (frho_6[itype][m]*p + frho_5[itype][m])*p + frho_4[itype][m];
    if (eflag) {
      phi = ((frho_3[itype][m]*p + frho_2[itype][m])*p + 
	     frho_1[itype][m])*p + frho_0[itype][m];
      eng_vdwl += phi;
    }
  }

  // communicate derivative of embedding function

  comm->comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms
  // FS has type-specific rhoip and rhojp

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
	jtype = type[j];
	r = sqrt(rsq);
	p = r*rdr + 1.0;
	m = static_cast<int> (p);
	m = MIN(m,nr-1);
	p -= m;
	p = MIN(p,1.0);

	// rhoip = derivative of (density at atom j due to atom i)
	// rhojp = derivative of (density at atom i due to atom j)
	// phi = pair potential energy
	// phip = phi'
	// z2 = phi * r
	// z2p = (phi * r)' = (phi' r) + phi
	// psip needs both fp[i] and fp[j] terms since r_ij appears in two
	//   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
	//   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

	rhoip =
	  (rhor_fs_6[itype][jtype][m]*p + rhor_fs_5[itype][jtype][m])*p + 
	  rhor_fs_4[itype][jtype][m];
	rhojp =
	  (rhor_fs_6[jtype][itype][m]*p + rhor_fs_5[jtype][itype][m])*p +
	  rhor_fs_4[jtype][itype][m];
	z2 = ((z2r_3[itype][jtype][m]*p + z2r_2[itype][jtype][m])*p + 
	      z2r_1[itype][jtype][m])*p + z2r_0[itype][jtype][m];
	z2p = (z2r_6[itype][jtype][m]*p + z2r_5[itype][jtype][m])*p + 
	  z2r_4[itype][jtype][m];
	
	recip = 1.0/r;
	phi = z2*recip;
	phip = z2p*recip - phi*recip;
	psip = fp[i]*rhojp + fp[j]*rhoip + phip;
	fforce = psip*recip;
	f[i][0] -= delx*fforce;
	f[i][1] -= dely*fforce;
	f[i][2] -= delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] += delx*fforce;
	  f[j][1] += dely*fforce;
	  f[j][2] += delz*fforce;
	}

	if (eflag) {
	  if (newton_pair || j < nlocal) eng_vdwl += phi;
	  else eng_vdwl += 0.5*phi;
	}

	if (vflag == 1) {
	  if (newton_pair || j < nlocal) {
	    virial[0] -= delx*delx*fforce;
	    virial[1] -= dely*dely*fforce;
	    virial[2] -= delz*delz*fforce;
	    virial[3] -= delx*dely*fforce;
	    virial[4] -= delx*delz*fforce;
	    virial[5] -= dely*delz*fforce;
	  } else {
	    virial[0] -= 0.5*delx*delx*fforce;
	    virial[1] -= 0.5*dely*dely*fforce;
	    virial[2] -= 0.5*delz*delz*fforce;
	    virial[3] -= 0.5*delx*dely*fforce;
	    virial[4] -= 0.5*delx*delz*fforce;
	    virial[5] -= 0.5*dely*delz*fforce;
	  }
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEAMFS::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all("Incorrect args for pair coefficients");

  // insure I,J args are * *
  // can only have one file, so clear all setflag settings

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // read Finnis/Sinclair modified setfl file
  // first clear setflag since are doing this once for I,J = *,*
  // read for all i,j pairs where ith,jth mapping is non-zero
  // set setflag i,j for non-zero pairs
  // set mass of atom type if i = j

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 0;

  int itable,ith,jth;
  int ilo,ihi,jlo,jhi;
  ilo = jlo = 1;
  ihi = jhi = n;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = jlo; j <= jhi; j++) {
      ith = atoi(arg[2+i]);
      jth = atoi(arg[2+j]);
      if (ith > 0 && jth > 0) {
	itable = read_setfl(arg[2],ith,jth);
	if (i == j) atom->set_mass(i,tables[itable].mass);
	tabindex[i][j] = itable;
	setflag[i][j] = 1;
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAMFS::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all("All EAM pair coeffs are not set");

  // EAM has only one cutoff = max of all pairwise cutoffs
  // determine max by checking table assigned to all type pairs
  // only setflag[i][j] = 1 is relevant (if hybrid, some may not be set)

  cutmax = 0.0;
  for (int ii = 1; ii <= atom->ntypes; ii++) {
    for (int jj = ii; jj <= atom->ntypes; jj++) {
      if (setflag[ii][jj] == 0) continue;
      cutmax = MAX(cutmax,tables[tabindex[ii][jj]].cut);
    }
  }

  return cutmax;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMFS::init_style()
{
  // set communication sizes in comm class

  comm->maxforward_pair = MAX(comm->maxforward_pair,1);
  comm->maxreverse_pair = MAX(comm->maxreverse_pair,1);

  // copy read-in-tables to multi-type setfl format
  // interpolate final spline coeffs
  
  store_setfl();
  interpolate();
  
  cutforcesq = cutmax*cutmax;
}

/* ----------------------------------------------------------------------
   read ith,jth potential values from a multi-element alloy EAM/FS file
   read values into table and bcast values
   this file has different format than standard EAM setfl file
------------------------------------------------------------------------- */

int PairEAMFS::read_setfl(char *file, int ith, int jth)
{
  // check if ith,jth portion of same file has already been read
  // if yes, return index of table entry
  // if no, extend table list

  for (int i = 0; i < ntables; i++)
    if (ith == tables[i].ith && jth == tables[i].jth) return i;

  tables = (Table *) 
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");

  Table *tb = &tables[ntables];
  int n = strlen(file) + 1;
  tb->filename = new char[n];
  strcpy(tb->filename,file);
  tb->ith = ith;
  tb->jth = jth;

  // open potential file

  int me = comm->me;
  FILE *fp;
  char line[MAXLINE];

  if (me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",file);
      error->one(str);
    }
  }

  // read and broadcast header

  int ntypes;
  if (me == 0) {
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d",&ntypes);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %d %lg %lg",
	   &tb->nrho,&tb->drho,&tb->nr,&tb->dr,&tb->cut);
  }

  MPI_Bcast(&ntypes,1,MPI_INT,0,world);
  MPI_Bcast(&tb->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&tb->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->nr,1,MPI_INT,0,world);
  MPI_Bcast(&tb->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->cut,1,MPI_DOUBLE,0,world);

  // check if ith,jth are consistent with ntypes

  if (ith > ntypes || jth > ntypes)
    error->all("Requested atom types in EAM setfl file do not exist");

  // allocate potential arrays and read/bcast them
  // skip sections of file that don't correspond to ith,jth
  // extract mass, frho for i,i from ith element section
  // extract rhor for i,j from jth array of ith element section
  // extract z2r for i,j from ith,jth array of z2r section
  // note that ith can be < or > than jth
  // set zr to NULL (funcl array) so it can be deallocated

  tb->frho = new double[tb->nrho+1];
  tb->rhor = new double[tb->nr+1];
  tb->z2r = new double[tb->nr+1];
  tb->zr = NULL;

  int i,j,tmp;
  double mass;

  for (i = 1; i <= ntypes; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fp);
      sscanf(line,"%d %lg",&tmp,&mass);
    }
    MPI_Bcast(&mass,1,MPI_DOUBLE,0,world);

    if (i == ith && ith == jth) {
      tb->mass = mass;
      if (me == 0) grab(fp,tb->nrho,&tb->frho[1]);
      MPI_Bcast(&tb->frho[1],tb->nrho,MPI_DOUBLE,0,world);

      for (j = 1; j <= ntypes; j++) {
	if (j == jth) {
	  if (me == 0) grab(fp,tb->nr,&tb->rhor[1]);
	  MPI_Bcast(&tb->rhor[1],tb->nr,MPI_DOUBLE,0,world);
	} else if (me == 0) skip(fp,tb->nr);
      }

    } else if (i == ith && ith != jth) {
      if (me == 0) skip(fp,tb->nrho);
      for (j = 1; j <= ntypes; j++) {
	if (j == jth) {
	  if (me == 0) grab(fp,tb->nr,&tb->rhor[1]);
	  MPI_Bcast(&tb->rhor[1],tb->nr,MPI_DOUBLE,0,world);
	} else if (me == 0) skip(fp,tb->nr);
      }

    } else {
      if (me == 0) skip(fp,tb->nrho);
      if (me == 0) for (j = 1; j <= ntypes; j = j + 1) skip(fp,tb->nr);
    }
  }

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= i; j++) {
      if ((i == ith && j == jth) || (j == ith && i == jth)) {
	if (me == 0) grab(fp,tb->nr,&tb->z2r[1]);
	MPI_Bcast(&tb->z2r[1],tb->nr,MPI_DOUBLE,0,world);
      } else if (me == 0) skip(fp,tb->nr);
    }
  }

  // close the potential file

  if (me == 0) fclose(fp);

  ntables++;
  return ntables-1;
}

/* ----------------------------------------------------------------------
   store read-in setfl values in multi-type setfl format
------------------------------------------------------------------------- */

void PairEAMFS::store_setfl()
{
  int i,j,m;

  int ntypes = atom->ntypes;
  
  // set nrho,nr,drho,dr from any i,i table entry since all the same

  for (i = 1; i <= ntypes; i++)
    if (setflag[i][i]) break;

  nrho = tables[tabindex[i][i]].nrho;
  nr = tables[tabindex[i][i]].nr;
  drho = tables[tabindex[i][i]].drho;
  dr = tables[tabindex[i][i]].dr;

  // allocate multi-type setfl arrays

  if (frho) {
    memory->destroy_2d_double_array(frho);
    memory->destroy_3d_double_array(rhor_fs);
    memory->destroy_2d_double_array(zrtmp);
    memory->destroy_3d_double_array(z2r);
  }

  frho = (double **) 
    memory->create_2d_double_array(ntypes+1,nrho+1,"eam/fs:frho");
  rhor_fs = (double ***)
    memory->create_3d_double_array(ntypes+1,ntypes+1,nr+1,"eam/fs:rhor_fs");
  zrtmp = (double **)
    memory->create_2d_double_array(ntypes+1,nr+1,"eam/fs:zrtmp");
  z2r = (double ***)
    memory->create_3d_double_array(ntypes+1,ntypes+1,nr+1,"eam/fs:frho");

  // copy from read-in tables to multi-type setfl arrays
  // frho,rhor are 1:ntypes, z2r is 1:ntypes,1:ntypes
  // skip if setflag i,j = 0 (if hybrid, some may not be set)

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) {
      if (setflag[i][j] == 0) continue;
      Table *tb = &tables[tabindex[i][j]];
      if (i == j) for (m = 1; m <= nrho; m++) frho[i][m] = tb->frho[m];
      for (m = 1; m <= nr; m++) {
	rhor_fs[i][j][m] = tb->rhor[m];
	z2r[i][j][m] = tb->z2r[m];
      }
    }
}

/* ----------------------------------------------------------------------
   interpolate EAM potentials
------------------------------------------------------------------------- */

void PairEAMFS::interpolate()
{
  // free memory from previous interpolation

  if (frho_0) interpolate_deallocate();

  // interpolation spacings

  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  // allocate coeff arrays

  int n = atom->ntypes;

  frho_0 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_0");
  frho_1 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_1");
  frho_2 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_2");
  frho_3 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_3");
  frho_4 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_4");
  frho_5 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_5");
  frho_6 = memory->create_2d_double_array(n+1,nrho+1,"eam/fs:frho_6");

  rhor_fs_0 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_0");
  rhor_fs_1 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_1");
  rhor_fs_2 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_2");
  rhor_fs_3 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_3");
  rhor_fs_4 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_4");
  rhor_fs_5 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_5");
  rhor_fs_6 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:rhor_fs_6");

  z2r_0 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_0");
  z2r_1 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_1");
  z2r_2 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_2");
  z2r_3 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_3");
  z2r_4 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_4");
  z2r_5 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_5");
  z2r_6 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam/fs:z2r_6");

  // frho interpolation for 1:ntypes
  // skip if setflag = 0 (if hybrid, some may not be set)
  // if skip, set frho arrays to 0.0, since they will still be accessed
  //   for non-EAM atoms when compute() calculates embedding function

  int i,j,m;

  for (i = 1; i <= atom->ntypes; i++) {
    if (setflag[i][i] == 0) {
      for (j = 1; j <= n; j++)
      	for (m = 1; m <= nrho; m++)
      	  frho_0[j][m] = frho_1[j][m] = frho_2[j][m] =  frho_3[j][m] =
      	    frho_4[j][m] = frho_5[j][m] = frho_6[j][m] = 0.0;
      continue;
    }

    for (m = 1; m <= nrho; m++) frho_0[i][m] = frho[i][m];

    frho_1[i][1] = frho_0[i][2]-frho_0[i][1];
    frho_1[i][2] = 0.5*(frho_0[i][3]-frho_0[i][1]);
    frho_1[i][nrho-1] = 0.5*(frho_0[i][nrho]-frho_0[i][nrho-2]);
    frho_1[i][nrho] = frho_0[i][nrho]-frho_0[i][nrho-1];

    for (m = 3; m <= nrho-2; m++)
      frho_1[i][m] = ((frho_0[i][m-2]-frho_0[i][m+2]) + 
		       8.0*(frho_0[i][m+1]-frho_0[i][m-1]))/12.0;

    for (m = 1; m <= nrho-1; m++) {
      frho_2[i][m] = 3.*(frho_0[i][m+1]-frho_0[i][m]) - 
	2.0*frho_1[i][m] - frho_1[i][m+1];
      frho_3[i][m] = frho_1[i][m] + frho_1[i][m+1] - 
	2.0*(frho_0[i][m+1]-frho_0[i][m]);
    }

    frho_2[i][nrho] = 0.0;
    frho_3[i][nrho] = 0.0;

    for (m = 1; m <= nrho; m++) {
      frho_4[i][m] = frho_1[i][m]/drho;
      frho_5[i][m] = 2.0*frho_2[i][m]/drho;
      frho_6[i][m] = 3.0*frho_3[i][m]/drho;
    }
  }

  // rhor interpolation for 1:ntypes
  // skip if setflag = 0 (if hybrid, some may not be set)

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = 1; j <= atom->ntypes; j++) {
      if (setflag[i][j] == 0) continue;

      for (m = 1; m <= nr; m++) rhor_fs_0[i][j][m] = rhor_fs[i][j][m];

      rhor_fs_1[i][j][1] = rhor_fs_0[i][j][2]-rhor_fs_0[i][j][1];
      rhor_fs_1[i][j][2] = 0.5*(rhor_fs_0[i][j][3]-rhor_fs_0[i][j][1]);
      rhor_fs_1[i][j][nr-1] = 0.5*(rhor_fs_0[i][j][nr]-rhor_fs_0[i][j][nr-2]);
      rhor_fs_1[i][j][nr] = 0.0;

      for (m = 3; m <= nr-2; m++)
	rhor_fs_1[i][j][m] = ((rhor_fs_0[i][j][m-2]-rhor_fs_0[i][j][m+2]) + 
			8.0*(rhor_fs_0[i][j][m+1]-rhor_fs_0[i][j][m-1]))/12.;

      for (m = 1; m <= nr-1; m++) {
	rhor_fs_2[i][j][m] = 3.0*(rhor_fs_0[i][j][m+1]-rhor_fs_0[i][j][m]) - 
	  2.0*rhor_fs_1[i][j][m] - rhor_fs_1[i][j][m+1];
	rhor_fs_3[i][j][m] = rhor_fs_1[i][j][m] + rhor_fs_1[i][j][m+1] - 
	  2.0*(rhor_fs_0[i][j][m+1]-rhor_fs_0[i][j][m]);
      }

      rhor_fs_2[i][j][nr] = 0.0;
      rhor_fs_3[i][j][nr] = 0.0;

      for (m = 1; m <= nr; m++) {
	rhor_fs_4[i][j][m] = rhor_fs_1[i][j][m]/dr;
	rhor_fs_5[i][j][m] = 2.0*rhor_fs_2[i][j][m]/dr;
	rhor_fs_6[i][j][m] = 3.0*rhor_fs_3[i][j][m]/dr;
      }
    }
  }

  // z2r interpolation for 1:ntypes,1:ntypes
  // skip if setflag i,j = 0 (if hybrid, some may not be set)
  // set j,i coeffs = i,j coeffs

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] == 0) continue;

      for (m = 1; m <= nr; m++) z2r_0[i][j][m] = z2r[i][j][m];

      z2r_1[i][j][1] = z2r_0[i][j][2]-z2r_0[i][j][1];
      z2r_1[i][j][2] = 0.5*(z2r_0[i][j][3]-z2r_0[i][j][1]);
      z2r_1[i][j][nr-1] = 0.5*(z2r_0[i][j][nr]-z2r_0[i][j][nr-2]);
      z2r_1[i][j][nr] = 0.0;

      for (m = 3; m <= nr-2; m++) 
	z2r_1[i][j][m] = ((z2r_0[i][j][m-2]-z2r_0[i][j][m+2]) + 
			   8.0*(z2r_0[i][j][m+1]-z2r_0[i][j][m-1]))/12.;

      for (m = 1; m <= nr-1; m++) {
	z2r_2[i][j][m] = 3.0*(z2r_0[i][j][m+1]-z2r_0[i][j][m]) - 
	  2.0*z2r_1[i][j][m] - z2r_1[i][j][m+1];
	z2r_3[i][j][m] = z2r_1[i][j][m] + z2r_1[i][j][m+1] - 
	  2.0*(z2r_0[i][j][m+1]-z2r_0[i][j][m]);
      }

      z2r_2[i][j][nr] = 0.0;
      z2r_3[i][j][nr] = 0.0;

      for (m = 1; m <= nr; m++) {
	z2r_4[i][j][m] = z2r_1[i][j][m]/dr;
	z2r_5[i][j][m] = 2.0*z2r_2[i][j][m]/dr;
	z2r_6[i][j][m] = 3.0*z2r_3[i][j][m]/dr;
      }

      for (m = 1; m <= nr; m++) {
	z2r_0[j][i][m] = z2r_0[i][j][m];
	z2r_1[j][i][m] = z2r_1[i][j][m];
	z2r_2[j][i][m] = z2r_2[i][j][m];
	z2r_3[j][i][m] = z2r_3[i][j][m];
	z2r_4[j][i][m] = z2r_4[i][j][m];
	z2r_5[j][i][m] = z2r_5[i][j][m];
	z2r_6[j][i][m] = z2r_6[i][j][m];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   deallocate spline interpolation arrays
------------------------------------------------------------------------- */

void PairEAMFS::interpolate_deallocate()
{
  memory->destroy_2d_double_array(frho_0);
  memory->destroy_2d_double_array(frho_1);
  memory->destroy_2d_double_array(frho_2);
  memory->destroy_2d_double_array(frho_3);
  memory->destroy_2d_double_array(frho_4);
  memory->destroy_2d_double_array(frho_5);
  memory->destroy_2d_double_array(frho_6);

  memory->destroy_3d_double_array(rhor_fs_0);
  memory->destroy_3d_double_array(rhor_fs_1);
  memory->destroy_3d_double_array(rhor_fs_2);
  memory->destroy_3d_double_array(rhor_fs_3);
  memory->destroy_3d_double_array(rhor_fs_4);
  memory->destroy_3d_double_array(rhor_fs_5);
  memory->destroy_3d_double_array(rhor_fs_6);

  memory->destroy_3d_double_array(z2r_0);
  memory->destroy_3d_double_array(z2r_1);
  memory->destroy_3d_double_array(z2r_2);
  memory->destroy_3d_double_array(z2r_3);
  memory->destroy_3d_double_array(z2r_4);
  memory->destroy_3d_double_array(z2r_5);
  memory->destroy_3d_double_array(z2r_6);
}

/* ---------------------------------------------------------------------- */

void PairEAMFS::single(int i, int j, int itype, int jtype,
		     double rsq, double factor_coul, double factor_lj,
		     int eflag, One &one)
{
  double r,p,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  int m;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);
  
  rhoip = (rhor_fs_6[itype][jtype][m]*p + rhor_fs_5[itype][jtype][m])*p +
    rhor_fs_4[itype][jtype][m];
  rhojp = (rhor_fs_6[jtype][itype][m]*p + rhor_fs_5[jtype][itype][m])*p +
    rhor_fs_4[jtype][itype][m];
  z2 = ((z2r_3[itype][jtype][m]*p + z2r_2[itype][jtype][m])*p +
	z2r_1[itype][jtype][m])*p + z2r_0[itype][jtype][m];
  z2p = (z2r_6[itype][jtype][m]*p + z2r_5[itype][jtype][m])*p + 
    z2r_4[itype][jtype][m];

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp[i]*rhojp + fp[j]*rhoip + phip;
  one.fforce = -psip*recip;

  if (eflag) {
    one.eng_vdwl = phi;
    one.eng_coul = 0.0;
  }
}
