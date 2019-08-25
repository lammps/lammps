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
   Contributing authors: Steven E. Strong and Joel D. Eaves
   Joel.Eaves@Colorado.edu
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_flow_gauss.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "error.h"
#include "citeme.h"
#include "respa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_flow_gauss[] =
  "Gaussian dynamics package:\n\n"
  "@Article{strong_water_2017,\n"
  "title = {The Dynamics of Water in Porous Two-Dimensional Crystals},\n"
  "volume = {121},\n"
  "number = {1},\n"
  "url = {http://dx.doi.org/10.1021/acs.jpcb.6b09387},\n"
  "doi = {10.1021/acs.jpcb.6b09387},\n"
  "urldate = {2016-12-07},\n"
  "journal = {J. Phys. Chem. B},\n"
  "author = {Strong, Steven E. and Eaves, Joel D.},\n"
  "year = {2017},\n"
  "pages = {189--207}\n"
  "}\n\n";

FixFlowGauss::FixFlowGauss(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_flow_gauss);

  if (narg < 6) error->all(FLERR,"Not enough input arguments");

  // a group which conserves momentum must also conserve particle number
  dynamic_group_allow = 0;

  scalar_flag = 1;
  vector_flag = 1;
  extscalar = 1;
  extvector = 1;
  size_vector = 3;
  global_freq = 1;    //data available every timestep
  respa_level_support = 1;
  //default respa level=outermost level is set in init()

  dimension = domain->dimension;

  //get inputs
  int tmpFlag;
  for (int ii=0; ii<3; ii++)
  {
    tmpFlag=force->inumeric(FLERR,arg[3+ii]);
    if (tmpFlag==1 || tmpFlag==0)
      flow[ii]=tmpFlag;
    else
      error->all(FLERR,"Constraint flags must be 1 or 0");
  }

  // by default, do not compute work done
  workflag=0;

  // process optional keyword
  int iarg = 6;
  while (iarg < narg) {
    if ( strcmp(arg[iarg],"energy") == 0 ) {
      if ( iarg+2 > narg ) error->all(FLERR,"Illegal energy keyword");
      if ( strcmp(arg[iarg+1],"yes") == 0 ) workflag = 1;
      else if ( strcmp(arg[iarg+1],"no") != 0 )
        error->all(FLERR,"Illegal energy keyword");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix flow/gauss command");
  }

  //error checking
  if (dimension == 2) {
    if (flow[2])
      error->all(FLERR,"Can't constrain z flow in 2d simulation");
  }

  dt=update->dt;
  pe_tot=0.0;
}

/* ---------------------------------------------------------------------- */

int FixFlowGauss::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFlowGauss::init()
{
  //if respa level specified by fix_modify, then override default (outermost)
  //if specified level too high, set to max level
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0)
      ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ----------------------------------------------------------------------
   setup is called after the initial evaluation of forces before a run, so we
   must remove the total force here too
   ------------------------------------------------------------------------- */
void FixFlowGauss::setup(int vflag)
{
  //need to compute work done if set fix_modify energy yes
  if (thermo_energy)
    workflag=1;

  //get total mass of group
  mTot=group->mass(igroup);
  if (mTot <= 0.0)
    error->all(FLERR,"Invalid group mass in fix flow/gauss");

  if (strstr(update->integrate_style,"respa")) {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
  else
    post_force(vflag);
}

/* ----------------------------------------------------------------------
   this is where Gaussian dynamics constraint is applied
   ------------------------------------------------------------------------- */
void FixFlowGauss::post_force(int /*vflag*/)
{
  double **f   = atom->f;
  double **v   = atom->v;

  int *mask    = atom->mask;
  int *type    = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  int nlocal   = atom->nlocal;

  int ii,jj;

  //find the total force on all atoms
  //initialize to zero
  double f_thisProc[3];
  for (ii=0; ii<3; ii++)
    f_thisProc[ii]=0.0;

  //add all forces on each processor
  for(ii=0; ii<nlocal; ii++)
    if (mask[ii] & groupbit)
      for (jj=0; jj<3; jj++)
        if (flow[jj])
          f_thisProc[jj] += f[ii][jj];

  //add the processor sums together
  MPI_Allreduce(f_thisProc, f_tot, 3, MPI_DOUBLE, MPI_SUM, world);

  //compute applied acceleration
  for (ii=0; ii<3; ii++)
    a_app[ii] = -f_tot[ii] / mTot;

  //apply added accelleration to each atom
  double f_app[3];
  double peAdded=0.0;
  for( ii = 0; ii<nlocal; ii++)
    if (mask[ii] & groupbit) {
      if (rmass) {
        f_app[0] = a_app[0]*rmass[ii];
        f_app[1] = a_app[1]*rmass[ii];
        f_app[2] = a_app[2]*rmass[ii];
      } else {
        f_app[0] = a_app[0]*mass[type[ii]];
        f_app[1] = a_app[1]*mass[type[ii]];
        f_app[2] = a_app[2]*mass[type[ii]];
      }

      f[ii][0] += f_app[0]; //f_app[jj] is 0 if flow[jj] is false
      f[ii][1] += f_app[1];
      f[ii][2] += f_app[2];

      //calculate added energy, since more costly, only do this if requested
      if (workflag)
        peAdded += f_app[0]*v[ii][0] + f_app[1]*v[ii][1] + f_app[2]*v[ii][2];
    }

  //finish calculation of work done, sum over all procs
  if (workflag) {
    double pe_tmp=0.0;
    MPI_Allreduce(&peAdded,&pe_tmp,1,MPI_DOUBLE,MPI_SUM,world);
    pe_tot += pe_tmp;
  }

}

void FixFlowGauss::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ----------------------------------------------------------------------
   negative of work done by this fix
   This is only computed if requested, either with fix_modify energy yes, or with the energy keyword. Otherwise returns 0.
   ------------------------------------------------------------------------- */
double FixFlowGauss::compute_scalar()
{
  return -pe_tot*dt;
}

/* ----------------------------------------------------------------------
   return components of applied force
   ------------------------------------------------------------------------- */
double FixFlowGauss::compute_vector(int n)
{
  return -f_tot[n];
}
