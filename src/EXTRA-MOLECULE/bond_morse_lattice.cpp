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
   Contributing author (bond_morse): Jeff Greathouse (SNL)
   Modification to morse/lattice : Thomas D Swinburne (CNRS / CINaM)
------------------------------------------------------------------------- */

#include "bond_morse_lattice.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondMorseLattice::BondMorseLattice(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondMorseLattice::~BondMorseLattice()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(b0);
    memory->destroy(bx);
    memory->destroy(by);
    memory->destroy(bz);
    memory->destroy(kappa);
  }
}

/* ---------------------------------------------------------------------- */

void BondMorseLattice::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond,kharm;
  double delsq,nbx,nby,nbz,db,balpha,bsq,bproj,bperpsq,fharmx,fharmy,fharmz;

  ebond = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    // displacement between i2 and i1. Note that the parts below are sign sensitive, so we take x[i2] - x[i1]
	delx = x[i2][0] - x[i1][0];
    dely = x[i2][1] - x[i1][1];
    delz = x[i2][2] - x[i1][2];
	delsq = delx*delx + dely*dely + delz*delz;

    // scalar project [delx,dely,delz] onto the input vector [bx,by,bz]
	bsq = bx[type]*bx[type] + by[type]*by[type] + bz[type]*bz[type];
    bproj = (delx*bx[type] + dely*by[type] + delz*bz[type]) / sqrt(bsq);
	
	// normalize the input vector
	nbx = bx[type]/bsq;
	nby = by[type]/bsq;
	nbz = bz[type]/bsq;
	
	// displacement of the projected bond from the equilibrium bond
	db = bproj - b0[type];

	// square of the perpendicular component of the bond obtained from the Pythagoras theorem
	bperpsq = delsq - bproj*bproj;

    // spring constant for the harmonic potential along the perpendicular
    kharm = kappa[type]*d0[type]*alpha[type]*alpha[type];
 
	// bond energy
	balpha = exp(-alpha[type]*db);
    if (eflag) {
		// along the bond
		ebond = d0[type]*(1-balpha)*(1-balpha);
		// perpendicular to the bond
		ebond += 0.5*kharm*bperpsq;
	}
	
	// force along the bond
	if (bproj > 0.0) fbond = -2.0*d0[type]*alpha[type]*(1-balpha)*balpha;
    else fbond = 0.0;
	
	// harmonic force
	fharmx = -kharm*(delx - nbx*bproj);
	fharmy = -kharm*(dely - nby*bproj);
	fharmz = -kharm*(delz - nbz*bproj);

    // apply the forces to the atoms. i2 gets positive force 
    if (newton_bond || i2 < nlocal) {
      f[i2][0] += fbond*nbx + fharmx;
      f[i2][1] += fbond*nby + fharmy;
      f[i2][2] += fbond*nbz + fharmz;
    }
	// i1 gets negative force
    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= fbond*nbx + fharmx;
      f[i1][1] -= fbond*nby + fharmy;
      f[i1][2] -= fbond*nbz + fharmz;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondMorseLattice::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(d0,n+1,"bond:d0");
  memory->create(alpha,n+1,"bond:alpha");
  memory->create(kappa,n+1,"bond:kappa");
  memory->create(bx,n+1,"bond:bx");
  memory->create(by,n+1,"bond:by");
  memory->create(bz,n+1,"bond:bz");
  memory->create(b0,n+1,"bond:b0");
  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondMorseLattice::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double d0_one = utils::numeric(FLERR,arg[1],false,lmp);
  double alpha_one = utils::numeric(FLERR,arg[2],false,lmp);
  double b0_one = utils::numeric(FLERR,arg[3],false,lmp);
  double kappa_one = utils::numeric(FLERR,arg[4],false,lmp);
  double bx_one = utils::numeric(FLERR,arg[5],false,lmp);
  double by_one = utils::numeric(FLERR,arg[6],false,lmp);
  double bz_one = utils::numeric(FLERR,arg[7],false,lmp);

  double bm_one = bx_one*bx_one + by_one*by_one + bz_one*bz_one;

  if (bm_one < 1.0e-6) error->all(FLERR,"Reference bond vector too small");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    d0[i] = d0_one;
    alpha[i] = alpha_one;
    b0[i] = b0_one;
    kappa[i] = kappa_one;
    bx[i] = bx_one;
    by[i] = by_one;
    bz[i] = bz_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondMorseLattice::equilibrium_distance(int i)
{
  return b0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondMorseLattice::write_restart(FILE *fp)
{
  fwrite(&d0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&alpha[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&b0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondMorseLattice::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&d0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&alpha[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&b0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
  }
  MPI_Bcast(&d0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&b0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondMorseLattice::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,d0[i],alpha[i],b0[i]);
}

/* ---------------------------------------------------------------------- */

double BondMorseLattice::single(int type, double rsq, int /*i*/, int /*j*/,
                         double &fforce)
{
  error->all(FLERR,"BondMorseLattice breaks single() assumptions- no rsq!");
  double r = sqrt(rsq);
  double dr = r - b0[type];
  double ralpha = exp(-alpha[type]*dr);
  fforce = 0;
  if (r > 0.0) fforce = -2.0*d0[type]*alpha[type]*(1-ralpha)*ralpha/r;
  return d0[type]*(1-ralpha)*(1-ralpha);
}

/* ---------------------------------------------------------------------- */

void *BondMorseLattice::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"b0")==0) return (void*) b0;
  return nullptr;
}
