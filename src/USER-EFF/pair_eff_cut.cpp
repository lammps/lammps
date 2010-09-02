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
   Contributing authors: Andres Jaramillo-Botero and Julius Su (Caltech)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eff_cut.h"
#include "pair_eff_inline.h"
#include "atom.h"
#include "update.h"
#include "min.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairEffCut::PairEffCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  nmax = 0;
  min_eradius = NULL;
  min_erforce = NULL;
}

/* ---------------------------------------------------------------------- */

PairEffCut::~PairEffCut()
{
  memory->sfree(min_eradius);
  memory->sfree(min_erforce);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    memory->destroy_2d_double_array(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairEffCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,energy;
  double fpair,fx,fy,fz,e1rforce,e2rforce,e1rvirial,e2rvirial;
  double rsq,rc,forcecoul,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *erforce = atom->erforce;
  double *eradius = atom->eradius;
  int *spin = atom->spin;	
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
  
    // add electron kinetic energy

    if (spin[i] != 0) {
      e1rforce = energy = 0;
      
      energy = 1.5 / (eradius[i] * eradius[i]);
      e1rforce = 3.0 / (eradius[i] * eradius[i] * eradius[i]);
      
      erforce[i] += e1rforce;
    
      // electronic ke accumulates into ecoul (pot)

      if (eflag) ecoul = energy;     // KE e-wavefunction
      if (evflag) {
        ev_tally_eff(i,i,nlocal,newton_pair,ecoul,0.0);
        if (flexible_pressure_flag)  // only on electron
          ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rforce*eradius[i]);
      }
    }
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      rc = sqrt(rsq);
      
      if (j < nall) factor_coul = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	j %= nall;
      }
      
      jtype = type[j];
      double taper = sqrt(cutsq[itype][jtype]);
      if (rsq < cutsq[itype][jtype]) {
	
	// nuclei-nuclei interaction

	if (spin[i] == 0 && spin[j] == 0) {
	  energy = fx = fy = fz = 0;
	  double qxq = qqrd2e*qtmp*q[j];
	  forcecoul = qxq/rsq;
	  
	  double dist = rc / taper;
	  double spline = cutoff(dist);
	  double dspline = dcutoff(dist) / taper;
	  
	  energy = factor_coul*qxq/rc;
	  fpair = forcecoul*spline-energy*dspline;
	  fpair = qqrd2e*fpair/rc;
	  energy = spline*energy;
	  
	  fx = delx*fpair;
	  fy = dely*fpair;
	  fz = delz*fpair;
	  
	  f[i][0] += fx;	
	  f[i][1] += fy;
	  f[i][2] += fz;
	  if (newton_pair || j < nlocal) {
	    f[j][0] -= fx;
	    f[j][1] -= fy;
	    f[j][2] -= fz;
	  }
	  
	  if (eflag) ecoul = energy;	// Electrostatics:N-N
          if (evflag)
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                                   fx,fy,fz,delx,dely,delz);
	}
	
	// I is nucleus, J is electron

	if (spin[i] == 0 && spin[j] != 0) {
	  energy = fpair = e1rforce = fx = fy = fz = e1rvirial = 0;
	  ElecNucElec(-q[i],rc,eradius[j],&energy,&fpair,&e1rforce,i,j);
	  
	  double dist = rc / taper;
	  double spline = cutoff(dist);
	  double dspline = dcutoff(dist) / taper;
	  
	  fpair = qqrd2e * (fpair * spline - energy * dspline);
	  energy = qqrd2e * spline * energy;
	  
	  e1rforce = qqrd2e * spline * e1rforce;
	  erforce[j] += e1rforce;
          e1rvirial = eradius[j] * e1rforce;
	  
	  SmallRForce(delx,dely,delz,rc,fpair,&fx,&fy,&fz);
	  f[i][0] += fx;
	  f[i][1] += fy;
	  f[i][2] += fz;
	  if (newton_pair || j < nlocal) {
	    f[j][0] -= fx;
	    f[j][1] -= fy;
	    f[j][2] -= fz;
	  }
	  
	  if (eflag) ecoul = energy;   // Electrostatics:N-e
	  if (evflag) {
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                                 fx,fy,fz,delx,dely,delz);
            if (flexible_pressure_flag) // only on electron
              ev_tally_eff(j,j,nlocal,newton_pair,0.0,e1rvirial);
          }
	}

	// I is electon, J is nucleus

	if (spin[i] != 0 && spin[j] == 0) {
	  energy = fpair = e1rforce = fx = fy = fz = e1rvirial = 0;
	  ElecNucElec(-q[j],rc,eradius[i],&energy,&fpair,&e1rforce,j,i);
	  
	  double dist = rc / taper;
	  double spline = cutoff(dist);
	  double dspline = dcutoff(dist) / taper;
	  
	  fpair = qqrd2e * (fpair * spline - energy * dspline);
	  energy = qqrd2e * spline * energy;
	  
	  e1rforce = qqrd2e * spline * e1rforce;
	  erforce[i] += e1rforce;
          e1rvirial = eradius[i] * e1rforce;
	  
	  SmallRForce(delx,dely,delz,rc,fpair,&fx,&fy,&fz);
	  f[i][0] += fx;
	  f[i][1] += fy;
	  f[i][2] += fz;
	  if (newton_pair || j < nlocal) {
	    f[j][0] -= fx;
	    f[j][1] -= fy;
	    f[j][2] -= fz;
	  }
	  
	  if (eflag) ecoul = energy;	//Electrostatics-e-N
	  if (evflag) {
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,ecoul,
                                   fx,fy,fz,delx,dely,delz);
            if (flexible_pressure_flag)  // only on electron
              ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rvirial);
          }
	}
	
	// electron-electron interaction

	if (spin[i] && spin[j]) {
	  energy = fpair = fx = fy= fz = 
	    e1rforce = e2rforce = e1rvirial = e2rvirial = 0.0;
	  ElecElecElec(rc,eradius[i],eradius[j],&energy,&fpair,
		       &e1rforce,&e2rforce,i,j);
	  
	  double s_energy, s_fpair, s_e1rforce, s_e2rforce;
	  s_energy = s_fpair = s_e1rforce = s_e2rforce = 0.0;
	  
          // as with the electron ke,
	  // the Pauli term is also accumulated into ecoul (pot)

	  PauliElecElec(spin[j] == spin[i],rc,eradius[i],eradius[j],
			&s_energy,&s_fpair,&s_e1rforce,&s_e2rforce,i,j);
	  
	  double dist = rc / taper;
	  double spline = cutoff(dist);
	  double dspline = dcutoff(dist) / taper;
	  
	  // apply spline cutoff

	  s_fpair = qqrd2e * (s_fpair * spline - s_energy * dspline);
	  s_energy = qqrd2e * spline * s_energy;
	  
	  fpair = qqrd2e * (fpair * spline - energy * dspline);
	  energy = qqrd2e * spline * energy;
	  
	  e1rforce = qqrd2e * spline * (e1rforce + s_e1rforce); 
	  e2rforce = qqrd2e * spline * (e2rforce + s_e2rforce); 
	  
	  // Cartesian and radial forces

	  SmallRForce(delx, dely, delz, rc, fpair + s_fpair, &fx, &fy, &fz);
	  erforce[i] += e1rforce;		
	  erforce[j] += e2rforce;	

          // radial virials

          e1rvirial = eradius[i] * e1rforce;
          e2rvirial = eradius[j] * e2rforce;
	  
	  f[i][0] += fx;
	  f[i][1] += fy;
	  f[i][2] += fz;
	  if (newton_pair || j < nlocal) {
	    f[j][0] -= fx;
	    f[j][1] -= fy;
	    f[j][2] -= fz;
	  }
	  
	  if (eflag) ecoul = energy + s_energy;  // Electrostatics+Pauli: e-e
          if (evflag) {
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,
                       ecoul,fx,fy,fz,delx,dely,delz);
            if (flexible_pressure_flag)          // on both electrons 
              ev_tally_eff(i,j,nlocal,newton_pair,0.0,e1rvirial+e2rvirial);
          }
	}
      }
    }
    
    // limit the electron size for periodic systems, to max=half-box-size
    // limit_size_stiffness for electrons

    if (spin[i] && limit_size_flag) {
      double half_box_length=0, dr, k=1.0; 
      e1rforce = energy = 0;
      
      if (domain->xperiodic == 1 || domain->yperiodic == 1 ||
	  domain->zperiodic == 1) {
	delx = domain->boxhi[0]-domain->boxlo[0];
	dely = domain->boxhi[1]-domain->boxlo[1];
	delz = domain->boxhi[2]-domain->boxlo[2];
	half_box_length = 0.5 * MIN(delx, MIN(dely, delz));
	if (eradius[i] > half_box_length) {
	  dr = eradius[i]-half_box_length;
	  energy=0.5*k*dr*dr;
	  e1rforce=-k*dr;
	}				
      }
      
      erforce[i] += e1rforce;

      // constraint radial energy accumulated as ecoul

      if (eflag) ecoul = energy;   // Radial constraint energy
      if (evflag) {
        ev_tally_eff(i,i,nlocal,newton_pair,ecoul,0.0);
        if (flexible_pressure_flag)  // only on electron
          ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
      }
    }		
  }
  
  if (vflag_fdotr) {
    virial_compute();
    if (flexible_pressure_flag) virial_eff_compute();
  }
}	

/* ----------------------------------------------------------------------
   eff-specific contribution to global virial
------------------------------------------------------------------------- */

void PairEffCut::virial_eff_compute()
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double e_virial;
  int *spin = atom->spin;

  // sum over force on all particles including ghosts
  
  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
    
  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts
    
  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
    
    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per-atom accumulators
   for virial radial electronic contributions
------------------------------------------------------------------------- */

void PairEffCut::ev_tally_eff(int i, int j, int nlocal, int newton_pair,
			      double ecoul, double e_virial)
{
  double ecoulhalf,epairhalf;
  double partial_evirial = e_virial/3.0;

  int *spin = atom->spin;

  // accumulate electronic wavefunction ke and radial constraint as ecoul

  if (eflag_either) {
    if (eflag_global) {
      ecoulhalf = 0.5*ecoul;
      if (i < nlocal)
        eng_coul += ecoulhalf;
      if (j < nlocal) 
        eng_coul += ecoulhalf;
    }
    if (eflag_atom) {
      epairhalf = 0.5 *  ecoul;
      if (i < nlocal) eatom[i] += epairhalf;
      if (j < nlocal) eatom[j] += epairhalf;
    }
  }
  
  if (vflag_either) {
    if (vflag_global) {
      if (spin[i] && i < nlocal) {
        virial[0] += 0.5*partial_evirial;
        virial[1] += 0.5*partial_evirial;
        virial[2] += 0.5*partial_evirial;
      }
      if (spin[j] && j < nlocal) {
        virial[0] += 0.5*partial_evirial;
        virial[1] += 0.5*partial_evirial;
        virial[2] += 0.5*partial_evirial;
      }
    }
    if (vflag_atom) {
      if (spin[i]) {
        if (newton_pair || i < nlocal) {
          vatom[i][0] += 0.5*partial_evirial;
          vatom[i][1] += 0.5*partial_evirial;
          vatom[i][2] += 0.5*partial_evirial;
        }
      }
      if (spin[j]) { 
        if (newton_pair || j < nlocal) {
          vatom[j][0] += 0.5*partial_evirial;
          vatom[j][1] += 0.5*partial_evirial;
          vatom[j][2] += 0.5*partial_evirial;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEffCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  
  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");
  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
}

/* ---------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEffCut::settings(int narg, char **arg)
{
  if (narg != 1 && narg != 3) error->all("Illegal pair_style command");

  if (narg == 1) {
    cut_global = force->numeric(arg[0]);
    limit_size_flag = 0;
    flexible_pressure_flag = 0;
  } else if (narg == 3) {
    cut_global = force->numeric(arg[0]);
    limit_size_flag = force->inumeric(arg[1]);
    flexible_pressure_flag = force->inumeric(arg[2]);
  }

  // reset cutoffs that have been explicitly set
  
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairEffCut::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();
  
  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);
  
  double cut_one = cut_global;
  if (narg == 3) cut_one = atof(arg[2]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  
  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEffCut::init_style()
{
  // error and warning checks

  if (!atom->q_flag || !atom->spin_flag || 
      !atom->eradius_flag || !atom->erforce_flag)
    error->all("Pair eff/cut requires atom attributes "
	       "q, spin, eradius, erforce");
  if (comm->ghost_velocity == 0)
    error->all("Pair eff/cut requires ghost atoms store velocity");

  // add hook to minimizer for eradius and erforce

  if (update->whichflag == 2)
    int ignore = update->minimize->request(this,1,0.01);
 
  // need a half neigh list and optionally a granular history neigh list
 
  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEffCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  
  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) fread(&cut[i][j],sizeof(double),1,fp);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   returns pointers to the log() of electron radius and corresponding force
   minimizer operates on log(radius) so radius never goes negative
   these arrays are stored locally by pair style
------------------------------------------------------------------------- */

void PairEffCut::min_xf_pointers(int ignore, double **xextra, double **fextra)
{
  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->sfree(min_eradius);
    memory->sfree(min_erforce);
    nmax = atom->nmax;
    min_eradius = (double *) memory->smalloc(nmax*sizeof(double),
					     "pair:min_eradius");
    min_erforce = (double *) memory->smalloc(nmax*sizeof(double),
					     "pair:min_erforce");
  }

  *xextra = min_eradius;
  *fextra = min_erforce;
}

/* ----------------------------------------------------------------------
   minimizer requests the log() of electron radius and corresponding force
   calculate and store in min_eradius and min_erforce
------------------------------------------------------------------------- */

void PairEffCut::min_xf_get(int ignore)
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) {
      min_eradius[i] = log(eradius[i]);
      min_erforce[i] = eradius[i]*erforce[i];
    } else min_eradius[i] = min_erforce[i] = 0.0;
}

/* ----------------------------------------------------------------------
   minimizer has changed the log() of electron radius
   propagate the change back to eradius
------------------------------------------------------------------------- */

void PairEffCut::min_x_set(int ignore)
{
  double *eradius = atom->eradius;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) 
    if (spin[i]) eradius[i] = exp(min_eradius[i]);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairEffCut::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
