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
       Contributing Author: Harsh Hemani and Manoj Warrier
       (Bhabha Atomc Research Center, Computational Analysis Division,
        Visakhapatnam, India)
       harshscience777@gmail.com, Manoj.Warrier@gmail.com
---------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstring>
#include "fix_esfriction.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "error.h"
#include "force.h"
#include "respa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixESFriction::FixESFriction(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix es-friction command");
  // args -- 0: fix, 1: fix-id, 2: group-id
  if (strcmp(update->unit_style,"metal") != 0) {
    error->all(FLERR, "fix es-friction supports only metal units");
  }
  Z1 = utils::numeric(FLERR, arg[3], false, lmp); // atomic number of PKA
  Z2 = utils::numeric(FLERR, arg[4], false, lmp); // atomic number of target
  A1 = utils::numeric(FLERR, arg[5], false, lmp); // mass number of PKA
  A2 = utils::numeric(FLERR, arg[6], false, lmp); // mass number of target
  Ne_val = utils::numeric(FLERR, arg[7], false, lmp);  // no. of valance shell electrons
  E_pka = utils::numeric(FLERR, arg[8], false, lmp);   // PKA energy
  /* ES_cutoff is the energy cutoff for applying ES.
     That is, atoms with energy below this will not experience ES */
  ES_cutoff = utils::numeric(FLERR, arg[9], false, lmp);
  // Compute total energy lost due to ES as per NRT model
  double E_hat = get_ehat();
  ES_loss_NRT = E_pka - E_hat;
  // following computes coefficient of friction (beta) using LS model
  double number_density = atom->natoms / 
    (domain->xprd * domain->yprd * domain->zprd) * 1e30; // number-per-m3
  double E_0 = h_cross*h_cross / (2*m_e) * 
    pow( 3*pi*pi*Ne_val*(number_density), 2.0/3.0) * joules_to_ev; // eV
  double v_0 = 5.931e5 * sqrt(E_0) * 1e-2; // ang/ps
  // lambda is in eV-Angstrom-ps units
  double lambda = 8.0 * pi * (e_charge * e_charge * 10) * bohr_radius * 
    pow(Z1, 1.0/6.0) * Z1*Z2 / ( v_0 * pow( pow(Z1, 2.0/3.0) + 
    pow(Z2, 2.0/3.0) , 3.0/2.0) ); 
  beta = number_density * lambda / 1.0e30; //eV-ps-angstrom^-2 units
  // Init energy lost due to ES as per LS model (calculated numerically)
  ES_loss_LS = 0.0;
  MPI_Comm_rank(world, &my_rank);
  if (my_rank==0){
    if(screen){
      fprintf(screen, "E_pka = %g eV\n", E_pka);
      fprintf(screen, "Energy loss due to ES (from NRT model) = %g\n", 
        ES_loss_NRT);
    }
    if (logfile){
      fprintf(logfile, "E_pka = %g eV\n", E_pka);
      fprintf(logfile, "Energy loss due to ES (from NRT model) = %g\n", 
        ES_loss_NRT);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixESFriction::get_kenergy(double A, double vx, double vy, double vz)
{
  double gms_to_kg = 1e-3;  // grams to kg
  double apps_to_mps = 100; // amstrong per picoseconds -to- m/s
  return (0.5 * A * ( vx*vx + vy*vy + vz*vz ) * gms_to_kg * apps_to_mps * 
    apps_to_mps * joules_to_ev / Nav); // in eV
}

/* ----------------------------------------------------------------------
   Ehat is total energy lost due to ES based on NRT model 
----------------------------------------------------------------------- */

double FixESFriction::get_ehat()
{
  double Emax = 25.0 * pow(Z1, 4.0/3.0) * A1;
  if (E_pka >= Emax){
    if (screen){
      fprintf(screen, "E_pka >= Emax (%g > %g) in esfriction fix.\n", E_pka, 
        Emax);
      fprintf(screen, "Hence E_hat calculation (NRT) might be invalid.\n");
    }
    if (logfile){
      fprintf(logfile, "E_pka >= Emax (%g > %g) in esfriction fix.\n", E_pka, 
        Emax);
      fprintf(logfile, "Hence E_hat calculation (NRT) might be invalid.\n");
    }
  }
  double eps, gofeps, kay;
  double const1 = pow( (9.0 * pi*pi / 128.0), (1.0/3.0) );
  double a = ( const1 * (bohr_radius /  10.0) ) / sqrt( pow(Z1, 2/3.0) + 
    pow(Z2, 2/3.0) ) ;
  eps = ( A2 / (A1+A2) ) * ( (a / ( Z1 * Z2 * e_charge * e_charge ) ) ) * 
    E_pka;
  gofeps = (3.4008 * pow(eps, 1/6.0)) + (0.40244 * pow(eps, 3.0/4.0)) + eps;
  kay = 0.1337 * pow(Z1, 1/6.0) * sqrt(Z1/A1);
  double Ehat = E_pka / (1.0 + kay * gofeps);
  return Ehat;
}

/* ---------------------------------------------------------------------- */
FixESFriction::~FixESFriction()
{
}
/* ---------------------------------------------------------------------- */

int FixESFriction::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixESFriction::init()
{
}
/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
Apply frictional force due to electronic stopping on each atom 
proportional to its velocity.
---------------------------------------------------------------------- */

void FixESFriction::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double E_i = 0.0; // energy of ith local atom
  double ES_loss_LS_local = 0.0; // ES loss (LS model) for this process rank
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      E_i = get_kenergy(A2, v[i][0], v[i][1], v[i][2]);
      if(E_i > ES_cutoff){ 
        f[i][0] -= beta * v[i][0];
        f[i][1] -= beta * v[i][1];
        f[i][2] -= beta * v[i][2];
        ES_loss_LS_local += beta * ( pow(v[i][0], 2) + pow(v[i][1], 2) + 
          pow(v[i][2], 2) ) * update->dt;
      }
    }
  }
  double ES_loss_step = 0; // energy lost in this time step
  MPI_Allreduce(&ES_loss_LS_local, &ES_loss_step, 1, MPI_DOUBLE, MPI_SUM, 
    world);
  ES_loss_LS += ES_loss_step; // energy lost till now
}

/* ---------------------------------------------------------------------- */


void FixESFriction::post_run()
{
  if(my_rank==0){
    if (screen){
      fprintf(screen, 
        "Energy Lost due to ES (based on NRT model): %g\n", ES_loss_NRT);
      fprintf(screen, 
        "Energy Lost due to ES (numerically computed using LS model): %g\n", 
        ES_loss_LS);
    }
  if (logfile){
      fprintf(logfile, 
        "Energy Lost due to ES (based on NRT model): %g\n", ES_loss_NRT);
      fprintf(logfile, 
        "Energy Lost due to ES (numerically computed using LS model): %g\n", 
        ES_loss_LS);
    }

  }
}
/* ---------------------------------------------------------------------- */
