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
   Contributing authors: Ezwanur Rahman, J.T. Foster (UTSA)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_peri_ves.h"
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_peri_neigh.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPeriVES::PairPeriVES(LAMMPS *lmp) : Pair(lmp)
{
  for (int i = 0; i < 6; i++) virial[i] = 0.0;
  no_virial_fdotr_compute = 1;
  single_enable = 0;

  ifix_peri = -1;

  nmax = 0;
  s0_new = NULL;
  theta = NULL;

  bulkmodulus = NULL;
  shearmodulus = NULL;
  s00 = alpha = NULL;
  cut = NULL;
  m_lambdai = NULL;
  m_taubi = NULL;

  // set comm size needed by this Pair
  // comm_reverse not needed

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

PairPeriVES::~PairPeriVES()
{
  if (ifix_peri >= 0) modify->delete_fix("PERI_NEIGH");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(bulkmodulus);
    memory->destroy(shearmodulus);
    memory->destroy(s00);
    memory->destroy(alpha);
    memory->destroy(cut);
    memory->destroy(m_lambdai);
    memory->destroy(m_taubi);
    memory->destroy(theta);
    memory->destroy(s0_new);
  }
}

/* ---------------------------------------------------------------------- */

void PairPeriVES::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp0,ytmp0,ztmp0,delx0,dely0,delz0,rsq0;
  double rsq,r,dr,rk,evdwl,fpair,fbond;
  double deltaed,fbondViscoElastic,fbondFinal;
  double decay,betai,lambdai,edbNp1,rkNew;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double d_ij,delta,stretch;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double timestepsize = update->dt;
  double *vfrac = atom->vfrac;
  double *s0 = atom->s0;
  double **x0 = atom->x0;
  double **r0 = ((FixPeriNeigh *) modify->fix[ifix_peri])->r0;
  double **deviatorextention = 
    ((FixPeriNeigh *) modify->fix[ifix_peri])->deviatorextention;
  double **deviatorBackextention = 
    ((FixPeriNeigh *) modify->fix[ifix_peri])->deviatorBackextention;
  tagint **partner = ((FixPeriNeigh *) modify->fix[ifix_peri])->partner;
  int *npartner = ((FixPeriNeigh *) modify->fix[ifix_peri])->npartner;
  double *wvolume = ((FixPeriNeigh *) modify->fix[ifix_peri])->wvolume;

  // lc = lattice constant
  // init_style guarantees it's the same in x, y, and z

  double lc = domain->lattice->xlattice;
  double half_lc = 0.5*lc;
  double vfrac_scale = 1.0;

  // short-range forces

  int newton_pair = force->newton_pair;
  int periodic = domain->xperiodic || domain->yperiodic || domain->zperiodic;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // need minimg() for x0 difference since not ghosted

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0,dely0,delz0);
      rsq0 = delx0*delx0 + dely0*dely0 + delz0*delz0;
      jtype = type[j];

      r = sqrt(rsq);

      // short-range interaction distance based on initial particle position
      // 0.9 and 1.35 are constants

      d_ij = MIN(0.9*sqrt(rsq0),1.35*lc);

      // short-range contact forces
      // 15 is constant taken from the EMU Theory Manual
      // Silling, 12 May 2005, p 18

      if (r < d_ij) {
        dr = r - d_ij;

        // kshort based upon short-range force constant
        // of the bond-based theory used in PMB model

        double kshort = (15.0 * 18.0 * bulkmodulus[itype][itype]) /
          (3.141592653589793 * cutsq[itype][jtype] * cutsq[itype][jtype]);
        rk = (kshort * vfrac[j]) * (dr / cut[itype][jtype]);

        if (r > 0.0) fpair = -(rk/r);
        else fpair = 0.0;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) evdwl = 0.5*rk*dr;
        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,
                             fpair*vfrac[i],delx,dely,delz);
      }
    }
  }

  // grow bond forces array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(s0_new);
    memory->destroy(theta);
    nmax = atom->nmax;
    memory->create(s0_new,nmax,"pair:s0_new");
    memory->create(theta,nmax,"pair:theta");
  }

  // Compute the dilatation on each particle
  compute_dilatation();

  // communicate dilatation (theta) of each particle
  comm->forward_comm_pair(this);

  // communicate weighted volume (wvolume) upon every reneighbor

  if (neighbor->ago == 0)
    comm->forward_comm_fix(modify->fix[ifix_peri]);

  // volume-dependent part of the energy

  if (eflag) {
    for (i = 0; i < nlocal; i++) {
      itype = type[i];
      if (eflag_global)
        eng_vdwl += 0.5 * bulkmodulus[itype][itype] * (theta[i] * theta[i]);
      if (eflag_atom)
        eatom[i] += 0.5 * bulkmodulus[itype][itype] * (theta[i] * theta[i]);
    }
  }

  // loop over my particles and their partners
  // partner list contains all bond partners, so I-J appears twice
  // if bond already broken, skip this partner
  // first = true if this is first neighbor of particle i

  bool first;
  double omega_minus, omega_plus;

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    itype = type[i];
    jnum = npartner[i];
    first = true;
    
    for (jj = 0; jj < jnum; jj++) {
      if (partner[i][jj] == 0) continue;
      j = atom->map(partner[i][jj]);

      // check if lost a partner without first breaking bond

      if (j < 0) {
        partner[i][jj] = 0;
        continue;
      }

      // compute force density, add to PD equation of motion

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (periodic) domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0,dely0,delz0);
      jtype = type[j];
      delta = cut[itype][jtype];
      r = sqrt(rsq);
      dr = r - r0[i][jj];

      // avoid roundoff errors

      if (fabs(dr) < 2.2204e-016) dr = 0.0;

      // scale vfrac[j] if particle j near the horizon

      if ((fabs(r0[i][jj] - delta)) <= half_lc)
        vfrac_scale = (-1.0/(2*half_lc))*(r0[i][jj]) +
          (1.0 + ((delta - half_lc)/(2*half_lc) ) );
      else vfrac_scale = 1.0;

      omega_plus  = influence_function(-1.0*delx0,-1.0*dely0,-1.0*delz0);
      omega_minus = influence_function(delx0,dely0,delz0);
        
      rk = ( (3.0 * bulkmodulus[itype][itype]) * vfrac[j] * vfrac_scale *
        ( (omega_plus * theta[i] / wvolume[i]) +
          ( omega_minus * theta[j] / wvolume[j] ) ) ) * r0[i][jj];

      if (r > 0.0) fbond = -(rk/r);
      else fbond = 0.0;

      // for viscoelasticity
      lambdai=m_lambdai[itype][itype];
      double taui = m_taubi[itype][itype];  
      double c1 = taui/timestepsize;
      decay=exp(-1.0/c1);
      betai=1.-c1*(1.-decay);

      double deviatoric_extension = 
        dr - (theta[i]* r0[i][jj] / 3.0);
      deltaed = deviatoric_extension-deviatorextention[i][jj];
 
      // back extention at current step

      edbNp1 = deviatorextention[i][jj]*(1-decay) + 
        deviatorBackextention[i][jj]*decay+betai*deltaed;

      rkNew = ((1-lambdai)*15.0) * 
        ( shearmodulus[itype][itype] * vfrac[j] * vfrac_scale ) *
        ( (omega_plus / wvolume[i]) + (omega_minus / wvolume[j]) ) * 
        deviatoric_extension;
      rkNew += (lambdai*15.0) * 
        ( shearmodulus[itype][itype] * vfrac[j] * vfrac_scale ) *
        ( (omega_plus / wvolume[i]) + (omega_minus / wvolume[j]) ) * 
        (deviatoric_extension-edbNp1);

      if (r > 0.0) fbondViscoElastic = -(rkNew/r);
      else fbondViscoElastic = 0.0;

      // total Force: elastic + viscoelastic 

      fbondFinal=fbond+fbondViscoElastic;
      fbond=fbondFinal;
         
      f[i][0] += delx*fbond;
      f[i][1] += dely*fbond;
      f[i][2] += delz*fbond;

      // since I-J is double counted, set newton off & use 1/2 factor and I,I

      if (eflag) evdwl =  ((0.5 * 15 * (1 - lambdai) * shearmodulus[itype][itype]/wvolume[i] *
                    omega_plus * deviatoric_extension * 
                    deviatoric_extension) + 
                    (0.5 * 15 * lambdai * shearmodulus[itype][itype]/wvolume[i] *
                    omega_plus * (deviatoric_extension-edbNp1) * 
                    (deviatoric_extension-edbNp1))) * vfrac[j] * vfrac_scale;
      if (evflag) ev_tally(i,i,nlocal,0,0.5*evdwl,0.0,
                           0.5*fbond*vfrac[i],delx,dely,delz);

      // find stretch in bond I-J and break if necessary
      // use s0 from previous timestep

      // store current deviatoric extention

      deviatorextention[i][jj]=deviatoric_extension;
      deviatorBackextention[i][jj]=edbNp1;

      stretch = dr / r0[i][jj];
      if (stretch > MIN(s0[i],s0[j])) partner[i][jj] = 0;

      // update s0 for next timestep

      if (first)
         s0_new[i] = s00[itype][jtype] - (alpha[itype][jtype] * stretch);
      else
         s0_new[i] = MAX(s0_new[i],s00[itype][jtype] -
                         (alpha[itype][jtype] * stretch));

      first = false;
    }  
  }

  // store new s0

  for (i = 0; i < nlocal; i++) s0[i] = s0_new[i];
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPeriVES::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(bulkmodulus,n+1,n+1,"pair:bulkmodulus");
  memory->create(shearmodulus,n+1,n+1,"pair:shearmodulus");
  memory->create(s00,n+1,n+1,"pair:s00");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(m_lambdai,n+1,n+1,"pair:m_lambdai");
  memory->create(m_taubi,n+1,n+1,"pair:m_taubi");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPeriVES::settings(int narg, char **arg)
{
  if (narg) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPeriVES::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double bulkmodulus_one = atof(arg[2]);
  double shearmodulus_one = atof(arg[3]);
  double cut_one = atof(arg[4]);
  double s00_one = atof(arg[5]);
  double alpha_one = atof(arg[6]);
  double mlambdai_one = atof(arg[7]);
  double mtaui_one = atof(arg[8]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      bulkmodulus[i][j] = bulkmodulus_one;
      shearmodulus[i][j] = shearmodulus_one;
      cut[i][j] = cut_one;
      s00[i][j] = s00_one;
      alpha[i][j] = alpha_one;
      m_lambdai[i][j] = mlambdai_one;
      m_taubi[i][j] = mtaui_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPeriVES::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  bulkmodulus[j][i] = bulkmodulus[i][j];
  shearmodulus[j][i] = shearmodulus[i][j];
  s00[j][i] = s00[i][j];
  alpha[j][i] = alpha[i][j];
  cut[j][i] = cut[i][j];
  m_lambdai[j][i] = m_lambdai[i][j];
  m_taubi[j][i] = m_taubi[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPeriVES::init_style()
{
  // error checks

  if (!atom->peri_flag) 
    error->all(FLERR,"Pair style peri requires atom style peri");
  if (atom->map_style == 0)
    error->all(FLERR,"Pair peri requires an atom map, see atom_modify");

  if (domain->lattice == NULL)
    error->all(FLERR,"Pair peri requires a lattice be defined");
  if (domain->lattice->xlattice != domain->lattice->ylattice ||
      domain->lattice->xlattice != domain->lattice->zlattice ||
      domain->lattice->ylattice != domain->lattice->zlattice)
    error->all(FLERR,"Pair peri lattice is not identical in x, y, and z");

  // if first init, create Fix needed for storing fixed neighbors

  if (ifix_peri == -1) {
    char **fixarg = new char*[3];
    fixarg[0] = (char *) "PERI_NEIGH";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "PERI_NEIGH";
    modify->add_fix(3,fixarg);
    delete [] fixarg;
  }

  // find associated PERI_NEIGH fix that must exist
  // could have changed locations in fix list since created

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"PERI_NEIGH") == 0) ifix_peri = i;
  if (ifix_peri == -1) error->all(FLERR,"Fix peri neigh does not exist");

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPeriVES::write_restart(FILE *fp)
{
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&bulkmodulus[i][j],sizeof(double),1,fp);
        fwrite(&shearmodulus[i][j],sizeof(double),1,fp);
        fwrite(&s00[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&m_lambdai[i][j],sizeof(double),1,fp);
        fwrite(&m_taubi[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPeriVES::read_restart(FILE *fp)
{
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&bulkmodulus[i][j],sizeof(double),1,fp);
          fread(&shearmodulus[i][j],sizeof(double),1,fp);
          fread(&s00[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          fread(&m_lambdai[i][j],sizeof(double),1,fp);
          fread(&m_taubi[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&bulkmodulus[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&shearmodulus[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&s00[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&m_lambdai[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&m_taubi[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairPeriVES::memory_usage()
{
  double bytes = 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   influence function definition
------------------------------------------------------------------------- */

double PairPeriVES::influence_function(double xi_x, double xi_y, double xi_z)
{
  double r = sqrt(xi_x*xi_x + xi_y*xi_y + xi_z*xi_z);
  double omega;

  if (fabs(r) < 2.2204e-016)
    error->one(FLERR,"Divide by 0 in influence function of pair peri/lps");
  omega = 1.0/r;
  return omega;
}

/* ---------------------------------------------------------------------- */

void PairPeriVES::compute_dilatation()
{
  int i,j,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double xtmp0,ytmp0,ztmp0,delx0,dely0,delz0;
  double rsq,r,dr;
  double delta;

  double **x = atom->x;
  int *type = atom->type;
  double **x0 = atom->x0;
  int nlocal = atom->nlocal;
  double *vfrac = atom->vfrac;
  double vfrac_scale = 1.0;

  double lc = domain->lattice->xlattice;
  double half_lc = 0.5*lc;


  double **r0   = ((FixPeriNeigh *) modify->fix[ifix_peri])->r0;
  tagint **partner = ((FixPeriNeigh *) modify->fix[ifix_peri])->partner;
  int *npartner = ((FixPeriNeigh *) modify->fix[ifix_peri])->npartner;
  double *wvolume = ((FixPeriNeigh *) modify->fix[ifix_peri])->wvolume;

  int periodic = domain->xperiodic || domain->yperiodic || domain->zperiodic;

  // compute the dilatation theta

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    xtmp0 = x0[i][0];
    ytmp0 = x0[i][1];
    ztmp0 = x0[i][2];
    jnum = npartner[i];
    theta[i] = 0.0;
    itype = type[i];

    for (jj = 0; jj < jnum; jj++) {

      // if bond already broken, skip this partner

      if (partner[i][jj] == 0) continue;

      // look up local index of this partner particle

      j = atom->map(partner[i][jj]);

      // skip if particle is "lost"

      if (j < 0) continue;

      // compute force density and add to PD equation of motion

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      if (periodic) domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      delx0 = xtmp0 - x0[j][0];
      dely0 = ytmp0 - x0[j][1];
      delz0 = ztmp0 - x0[j][2];
      if (periodic) domain->minimum_image(delx0,dely0,delz0);

      r = sqrt(rsq);
      dr = r - r0[i][jj];
      if (fabs(dr) < 2.2204e-016) dr = 0.0;

      jtype = type[j];
      delta = cut[itype][jtype];

      // scale vfrac[j] if particle j near the horizon

      if ((fabs(r0[i][jj] - delta)) <= half_lc)
        vfrac_scale = (-1.0/(2*half_lc))*(r0[i][jj]) +
          (1.0 + ((delta - half_lc)/(2*half_lc) ) );
      else vfrac_scale = 1.0;

      theta[i] += influence_function(delx0, dely0, delz0) * r0[i][jj] * dr *
        vfrac[j] * vfrac_scale;

    }

    // if wvolume[i] is zero, then particle i has no bonds
    // therefore, the dilatation is set to

    if (wvolume[i] != 0.0) theta[i] = (3.0/wvolume[i]) * theta[i];
    else theta[i] = 0;
  }
}


/* ----------------------------------------------------------------------
   communication routines
---------------------------------------------------------------------- */

int PairPeriVES::pack_comm(int n, int *list, double *buf,
                           int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = theta[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairPeriVES::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    theta[i] = buf[m++];
  }
}
