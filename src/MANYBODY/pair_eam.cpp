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
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eam.h"
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

PairEAM::PairEAM()
{
  nmax = 0;
  rho = NULL;
  fp = NULL;
  
  ntables = 0;
  tables = NULL;
  frho = NULL;
  frho_0 = NULL;

  // set rhor to NULL so memory deallocation will work
  // even from derived classes that don't use rhor

  rhor = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAM::~PairEAM()
{
  memory->sfree(rho);
  memory->sfree(fp);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    memory->destroy_2d_int_array(tabindex);
  }

  for (int m = 0; m < ntables; m++) {
    delete [] tables[m].filename;
    delete [] tables[m].frho;
    delete [] tables[m].rhor;
    delete [] tables[m].zr;
    delete [] tables[m].z2r;
  }
  memory->sfree(tables);

  if (frho) {
    memory->destroy_2d_double_array(frho);
    memory->destroy_2d_double_array(rhor);
    memory->destroy_2d_double_array(zrtmp);
    memory->destroy_3d_double_array(z2r);
  }

  if (frho_0) interpolate_deallocate();
}

/* ---------------------------------------------------------------------- */

void PairEAM::compute(int eflag, int vflag)
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
	rho[i] += ((rhor_3[jtype][m]*p + rhor_2[jtype][m])*p + 
		   rhor_1[jtype][m])*p + rhor_0[jtype][m];
	if (newton_pair || j < nlocal)
	  rho[j] += ((rhor_3[itype][m]*p + rhor_2[itype][m])*p + 
		     rhor_1[itype][m])*p + rhor_0[itype][m];
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

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

	rhoip = (rhor_6[itype][m]*p + rhor_5[itype][m])*p + 
	  rhor_4[itype][m];
	rhojp = (rhor_6[jtype][m]*p + rhor_5[jtype][m])*p + 
	  rhor_4[jtype][m];
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
   allocate all arrays
------------------------------------------------------------------------- */

void PairEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");
  tabindex = memory->create_2d_int_array(n+1,n+1,"pair:tabindex");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEAM::settings(int narg, char **arg)
{
  if (narg > 0) error->all("Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   reading multiple funcfl files defines a funcfl alloy simulation
------------------------------------------------------------------------- */

void PairEAM::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all("Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // read funcfl file only for i,i pairs
  // only setflag i,i will be set
  // set mass of each atom type

  int itable;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (i == j) {
	itable = read_funcfl(arg[2]);
	atom->set_mass(i,tables[itable].mass);
	tabindex[i][i] = itable;
	setflag[i][i] = 1;
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAM::init_one(int i, int j)
{
  // only setflag I,I was set by coeff
  // mixing will occur in init_style if both I,I and J,J were set

  if (setflag[i][i] == 0 || setflag[j][j] == 0)
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

void PairEAM::init_style()
{
  // set communication sizes in comm class

  comm->maxforward_pair = MAX(comm->maxforward_pair,1);
  comm->maxreverse_pair = MAX(comm->maxreverse_pair,1);

  // convert read-in funcfl tables to multi-type setfl format and mix I,J
  // interpolate final spline coeffs
  
  convert_funcfl();
  interpolate();
  
  cutforcesq = cutmax*cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a single element EAM file
   read values into table and bcast values
------------------------------------------------------------------------- */

int PairEAM::read_funcfl(char *file)
{
  // check if same file has already been read
  // if yes, return index of table entry
  // if no, extend table list

  for (int i = 0; i < ntables; i++)
    if (strcmp(file,tables->filename) == 0) return i;

  tables = (Table *) 
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");

  Table *tb = &tables[ntables];
  int n = strlen(file) + 1;
  tb->filename = new char[n];
  strcpy(tb->filename,file);
  tb->ith = tb->jth = 0;

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

  int tmp;
  if (me == 0) {
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg",&tmp,&tb->mass);
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %d %lg %lg",
	   &tb->nrho,&tb->drho,&tb->nr,&tb->dr,&tb->cut);
  }

  MPI_Bcast(&tb->mass,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&tb->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->nr,1,MPI_INT,0,world);
  MPI_Bcast(&tb->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->cut,1,MPI_DOUBLE,0,world);

  // allocate potential arrays and read/bcast them
  // set z2r to NULL (setfl array) so it can be deallocated

  tb->frho = new double[tb->nrho+1];
  tb->zr = new double[tb->nr+1];
  tb->rhor = new double[tb->nr+1];
  tb->z2r = NULL;

  if (me == 0) grab(fp,tb->nrho,&tb->frho[1]);
  MPI_Bcast(&tb->frho[1],tb->nrho,MPI_DOUBLE,0,world);

  if (me == 0) grab(fp,tb->nr,&tb->zr[1]);
  MPI_Bcast(&tb->zr[1],tb->nr,MPI_DOUBLE,0,world);

  if (me == 0) grab(fp,tb->nr,&tb->rhor[1]);
  MPI_Bcast(&tb->rhor[1],tb->nr,MPI_DOUBLE,0,world);

  // close the potential file

  if (me == 0) fclose(fp);

  ntables++;
  return ntables-1;
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potentials to multi-type setfl format
------------------------------------------------------------------------- */

void PairEAM::convert_funcfl()
{
  int i,j,k,m;

  int ntypes = atom->ntypes;

  // determine max values for all i,i type pairs
  // skip if setflag = 0 (if hybrid, some may not be set)

  double rmax,rhomax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 1; i <= ntypes; i++) {
    if (setflag[i][i] == 0) continue;
    Table *tb = &tables[tabindex[i][i]];
    dr = MAX(dr,tb->dr);
    drho = MAX(drho,tb->drho);
    rmax = MAX(rmax,(tb->nr-1) * tb->dr);
    rhomax = MAX(rhomax,(tb->nrho-1) * tb->drho);
  }

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int> (rmax/dr + 0.5);
  nrho = static_cast<int> (rhomax/drho + 0.5);

  // allocate multi-type setfl arrays

  if (frho) {
    memory->destroy_2d_double_array(frho);
    memory->destroy_2d_double_array(rhor);
    memory->destroy_2d_double_array(zrtmp);
    memory->destroy_3d_double_array(z2r);
  }

  frho = (double **) 
    memory->create_2d_double_array(ntypes+1,nrho+1,"eam:frho");
  rhor = (double **)
    memory->create_2d_double_array(ntypes+1,nr+1,"eam:rhor");
  zrtmp = (double **)
    memory->create_2d_double_array(ntypes+1,nr+1,"eam:zrtmp");
  z2r = (double ***)
    memory->create_3d_double_array(ntypes+1,ntypes+1,nr+1,"eam:frho");

  // interpolate all potentials to a single grid and cutoff for all atom types
  // frho,rhor are 1:ntypes, z2r is 1:ntypes,1:ntypes
  // skip if setflag i,i or j,j = 0 (if hybrid, some may not be set)

  double r,p,cof1,cof2,cof3,cof4;
  
  for (i = 1; i <= ntypes; i++) {
    if (setflag[i][i] == 0) continue;
    Table *tb = &tables[tabindex[i][i]];
    for (m = 1; m <= nrho; m++) {
      r = (m-1)*drho;
      p = r/tb->drho + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,tb->nrho-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = 0.166666667*p*(p*p-1.0);
      frho[i][m] = cof1*tb->frho[k-1] + cof2*tb->frho[k] + 
	cof3*tb->frho[k+1] + cof4*tb->frho[k+2];
    }
  }

  for (i = 1; i <= ntypes; i++) {
    if (setflag[i][i] == 0) continue;
    Table *tb = &tables[tabindex[i][i]];
    for (m = 1; m <= nr; m++) {
      r = (m-1)*dr;
      p = r/tb->dr + 1.0;
      k = static_cast<int> (p);
      k = MIN(k,tb->nr-2);
      k = MAX(k,2);
      p -= k;
      p = MIN(p,2.0);
      cof1 = -0.166666667*p*(p-1.0)*(p-2.0);
      cof2 = 0.5*(p*p-1.0)*(p-2.0);
      cof3 = -0.5*p*(p+1.0)*(p-2.0);
      cof4 = 0.166666667*p*(p*p-1.0);
      rhor[i][m] = cof1*tb->rhor[k-1] + cof2*tb->rhor[k] +
	cof3*tb->rhor[k+1] + cof4*tb->rhor[k+2];
      zrtmp[i][m] = cof1*tb->zr[k-1] + cof2*tb->zr[k] +
	cof3*tb->zr[k+1] + cof4*tb->zr[k+2];
    }
  }

  for (i = 1; i <= ntypes; i++)
    for (j = i; j <= ntypes; j++) {
      if (setflag[i][i] == 0 || setflag[j][j] == 0) continue;
      for (m = 1; m <= nr; m++)
	z2r[i][j][m] = 27.2*0.529 * zrtmp[i][m]*zrtmp[j][m];
    }
}

/* ----------------------------------------------------------------------
   interpolate EAM potentials
------------------------------------------------------------------------- */

void PairEAM::interpolate()
{
  // free memory from previous interpolation

  if (frho_0) interpolate_deallocate();

  // interpolation spacings

  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  // allocate coeff arrays

  int n = atom->ntypes;

  frho_0 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_0");
  frho_1 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_1");
  frho_2 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_2");
  frho_3 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_3");
  frho_4 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_4");
  frho_5 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_5");
  frho_6 = memory->create_2d_double_array(n+1,nrho+1,"eam:frho_6");

  rhor_0 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_0");
  rhor_1 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_1");
  rhor_2 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_2");
  rhor_3 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_3");
  rhor_4 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_4");
  rhor_5 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_5");
  rhor_6 = memory->create_2d_double_array(n+1,nr+1,"eam:rhor_6");

  z2r_0 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_0");
  z2r_1 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_1");
  z2r_2 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_2");
  z2r_3 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_3");
  z2r_4 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_4");
  z2r_5 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_5");
  z2r_6 = memory->create_3d_double_array(n+1,n+1,nr+1,"eam:z2r_6");

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
    if (setflag[i][i] == 0) continue;

    for (m = 1; m <= nr; m++) rhor_0[i][m] = rhor[i][m];

    rhor_1[i][1] = rhor_0[i][2]-rhor_0[i][1];
    rhor_1[i][2] = 0.5*(rhor_0[i][3]-rhor_0[i][1]);
    rhor_1[i][nr-1] = 0.5*(rhor_0[i][nr]-rhor_0[i][nr-2]);
    rhor_1[i][nr] = 0.0;

    for (m = 3; m <= nr-2; m++)
      rhor_1[i][m] = ((rhor_0[i][m-2]-rhor_0[i][m+2]) + 
		       8.0*(rhor_0[i][m+1]-rhor_0[i][m-1]))/12.;

    for (m = 1; m <= nr-1; m++) {
      rhor_2[i][m] = 3.0*(rhor_0[i][m+1]-rhor_0[i][m]) - 
	2.0*rhor_1[i][m] - rhor_1[i][m+1];
      rhor_3[i][m] = rhor_1[i][m] + rhor_1[i][m+1] - 
	2.0*(rhor_0[i][m+1]-rhor_0[i][m]);
    }

    rhor_2[i][nr] = 0.0;
    rhor_3[i][nr] = 0.0;

    for (m = 1; m <= nr; m++) {
      rhor_4[i][m] = rhor_1[i][m]/dr;
      rhor_5[i][m] = 2.0*rhor_2[i][m]/dr;
      rhor_6[i][m] = 3.0*rhor_3[i][m]/dr;
    }
  }

  // z2r interpolation for 1:ntypes,1:ntypes
  // skip if setflag i,i or j,j = 0 (if hybrid, some may not be set)
  // set j,i coeffs = i,j coeffs

  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (setflag[i][i] == 0 || setflag[j][j] == 0) continue;

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
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairEAM::grab(FILE *fp, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fp);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while (ptr = strtok(NULL," \t\n\r\f")) list[i++] = atof(ptr);
  }
}

/* ----------------------------------------------------------------------
   skip n values from file fp
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairEAM::skip(FILE *fp, int n)
{
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fp);
    strtok(line," \t\n\r\f");
    i++;
    while (strtok(NULL," \t\n\r\f")) i++;
  }
}

/* ----------------------------------------------------------------------
   deallocate spline interpolation arrays
------------------------------------------------------------------------- */

void PairEAM::interpolate_deallocate()
{
  memory->destroy_2d_double_array(frho_0);
  memory->destroy_2d_double_array(frho_1);
  memory->destroy_2d_double_array(frho_2);
  memory->destroy_2d_double_array(frho_3);
  memory->destroy_2d_double_array(frho_4);
  memory->destroy_2d_double_array(frho_5);
  memory->destroy_2d_double_array(frho_6);

  memory->destroy_2d_double_array(rhor_0);
  memory->destroy_2d_double_array(rhor_1);
  memory->destroy_2d_double_array(rhor_2);
  memory->destroy_2d_double_array(rhor_3);
  memory->destroy_2d_double_array(rhor_4);
  memory->destroy_2d_double_array(rhor_5);
  memory->destroy_2d_double_array(rhor_6);

  memory->destroy_3d_double_array(z2r_0);
  memory->destroy_3d_double_array(z2r_1);
  memory->destroy_3d_double_array(z2r_2);
  memory->destroy_3d_double_array(z2r_3);
  memory->destroy_3d_double_array(z2r_4);
  memory->destroy_3d_double_array(z2r_5);
  memory->destroy_3d_double_array(z2r_6);
}

/* ---------------------------------------------------------------------- */

void PairEAM::single(int i, int j, int itype, int jtype,
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

  rhoip = (rhor_6[itype][m]*p + rhor_5[itype][m])*p + 
    rhor_4[itype][m];
  rhojp = (rhor_6[jtype][m]*p + rhor_5[jtype][m])*p + 
    rhor_4[jtype][m];
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

/* ---------------------------------------------------------------------- */

void PairEAM::single_embed(int i, int itype, double &fpi,
			   int eflag, double &phi)
{
  double p = rho[i]*rdrho + 1.0;
  int m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;

  fpi = (frho_6[itype][m]*p + frho_5[itype][m])*p + frho_4[itype][m];
  if (eflag)
    phi = ((frho_3[itype][m]*p + frho_2[itype][m])*p + 
	   frho_1[itype][m])*p + frho_0[itype][m];
}

/* ---------------------------------------------------------------------- */

int PairEAM::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAM::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairEAM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

int PairEAM::memory_usage()
{
  int bytes = 2 * nmax * sizeof(double);
  return bytes;
}
