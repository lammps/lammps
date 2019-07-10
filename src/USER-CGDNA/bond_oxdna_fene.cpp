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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "bond_oxdna_fene.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondOxdnaFene::BondOxdnaFene(LAMMPS *lmp) : Bond(lmp)
{

}

/* ---------------------------------------------------------------------- */

BondOxdnaFene::~BondOxdnaFene()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(Delta);
    memory->destroy(r0);

  }
}


/* ----------------------------------------------------------------------
    compute vector COM-sugar-phosphate backbone interaction site in oxDNA
------------------------------------------------------------------------- */
void BondOxdnaFene::compute_interaction_sites(double e1[3],
  double /*e2*/[3], double r[3])
{
  double d_cs=-0.4;

  r[0] = d_cs*e1[0];
  r[1] = d_cs*e1[1];
  r[2] = d_cs*e1[2];

}

/* ----------------------------------------------------------------------
   compute function for oxDNA FENE-bond interaction
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */
void BondOxdnaFene::compute(int eflag, int vflag)
{
  int a,b,in,type;
  double delf[3],delta[3],deltb[3]; // force, torque increment;;
  double delr[3],ebond,fbond;
  double rsq,Deltasq,rlogarg;
  double r,rr0,rr0sq;
  // vectors COM-backbone site in lab frame
  double ra_cs[3],rb_cs[3];

  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  ebond = 0.0;
  ev_init(eflag,vflag);

  // loop over FENE bonds

  for (in = 0; in < nbondlist; in++) {

    a = bondlist[in][1];
    b = bondlist[in][0];
    type = bondlist[in][2];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);
    qb=bonus[b].quat;
    MathExtra::q_to_exyz(qb,bx,by,bz);

    // vector COM-backbone site a and b
    compute_interaction_sites(ax,ay,ra_cs);
    compute_interaction_sites(bx,by,rb_cs);

    // vector backbone site b to a
    delr[0] = x[a][0] + ra_cs[0] - x[b][0] - rb_cs[0];
    delr[1] = x[a][1] + ra_cs[1] - x[b][1] - rb_cs[1];
    delr[2] = x[a][2] + ra_cs[2] - x[b][2] - rb_cs[2];
    rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];
    r = sqrt(rsq);

    rr0 = r - r0[type];
    rr0sq = rr0*rr0;
    Deltasq = Delta[type] * Delta[type];
    rlogarg = 1.0 - rr0sq/Deltasq;

    // if r -> Delta, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*Delta something serious is wrong, abort

    if (rlogarg < 0.1) {
      char str[128];
      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " "
              TAGINT_FORMAT " " TAGINT_FORMAT " %g",
              update->ntimestep,atom->tag[a],atom->tag[b],r);
      error->warning(FLERR,str,0);
      if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
    }

    fbond = -k[type]*rr0/rlogarg/Deltasq/r;
    delf[0] = delr[0]*fbond;
    delf[1] = delr[1]*fbond;
    delf[2] = delr[2]*fbond;

    // energy

    if (eflag) {
      ebond = -0.5 * k[type]*log(rlogarg);
    }

    // apply force and torque to each of 2 atoms

    if (newton_bond || a < nlocal) {

      f[a][0] += delf[0];
      f[a][1] += delf[1];
      f[a][2] += delf[2];

      MathExtra::cross3(ra_cs,delf,delta);

      torque[a][0] += delta[0];
      torque[a][1] += delta[1];
      torque[a][2] += delta[2];

    }

    if (newton_bond || b < nlocal) {

      f[b][0] -= delf[0];
      f[b][1] -= delf[1];
      f[b][2] -= delf[2];

      MathExtra::cross3(rb_cs,delf,deltb);

      torque[b][0] -= deltb[0];
      torque[b][1] -= deltb[1];
      torque[b][2] -= deltb[2];

    }

    // increment energy and virial
    if (evflag) ev_tally(a,b,nlocal,newton_bond,ebond,fbond,delr[0],delr[1],delr[2]);

  }

}

/* ---------------------------------------------------------------------- */

void BondOxdnaFene::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(Delta,n+1,"bond:Delta");
  memory->create(r0,n+1,"bond:r0");
  memory->create(setflag,n+1,"bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondOxdnaFene::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for bond coefficients in oxdna/fene");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double Delta_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);

  int count = 0;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    Delta[i] = Delta_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients in oxdna/fene");

}

/* ----------------------------------------------------------------------
   set special_bond settings and check if valid
------------------------------------------------------------------------- */

void BondOxdnaFene::init_style()
{
  /* special bonds have to be lj = 0 1 1 and coul = 1 1 1 to exclude
     the ss excluded volume interaction between nearest neighbors   */

  force->special_lj[1] = 0.0;
  force->special_lj[2] = 1.0;
  force->special_lj[3] = 1.0;
  force->special_coul[1] = 1.0;
  force->special_coul[2] = 1.0;
  force->special_coul[3] = 1.0;

  fprintf(screen,"Finding 1-2 1-3 1-4 neighbors ...\n"
                 " Special bond factors lj:   %-10g %-10g %-10g\n"
                 " Special bond factors coul: %-10g %-10g %-10g\n",
                 force->special_lj[1],force->special_lj[2],force->special_lj[3],
                 force->special_coul[1],force->special_coul[2],force->special_coul[3]);

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0 ||
      force->special_coul[1] != 1.0 || force->special_coul[2] != 1.0 || force->special_coul[3] != 1.0)
  {
    if (comm->me == 0)
      error->warning(FLERR,"Use special bonds lj = 0,1,1 and coul = 1,1,1 with bond style oxdna/fene");
  }

}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondOxdnaFene::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Delta[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondOxdnaFene::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Delta[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Delta[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondOxdnaFene::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,k[i],r0[i],Delta[i]);
}

/* ---------------------------------------------------------------------- */

double BondOxdnaFene::single(int type, double rsq, int /*i*/, int /*j*/,
                        double &fforce)
{
  double r = sqrt(rsq);
  double rr0 = r - r0[type];
  double rr0sq = rr0*rr0;
  double Deltasq = Delta[type] * Delta[type];
  double rlogarg = 1.0 - rr0sq/Deltasq;

  // if r -> Delta, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*Delta something serious is wrong, abort

  if (rlogarg < 0.1) {
    char str[128];
    sprintf(str,"FENE bond too long: " BIGINT_FORMAT " %g",
            update->ntimestep,sqrt(rsq));
    error->warning(FLERR,str,0);
    if (rlogarg <= -3.0) error->one(FLERR,"Bad FENE bond");
  }

  double eng = -0.5 * k[type]*log(rlogarg);
  fforce = -k[type]*rr0/rlogarg/Deltasq/r;

  return eng;
}
