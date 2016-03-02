/* ----------------------------------------------------------------------
   Lammps - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the USER-MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License v.2.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterised by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "math.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"


#include "manifold_thermostat_nh.h"
#include "manifold_rattle_integrator.h"

using namespace LAMMPS_NS;
using namespace FixConst;


enum {CONSTANT,EQUAL};
enum {NOBIAS,BIAS};





ManifoldThermostatNH::ManifoldThermostatNH( LAMMPS *lmp, int &narg,
                                            char ***arg_ptr, const char *f_id,
                                            int argi, int igroup, int groupbit )
  : Pointers(lmp), id_temp(NULL), temperature(NULL), eta(NULL), eta_dot(NULL),
    eta_dotdot(NULL), eta_mass(NULL), fix_id(f_id),
    igroup(igroup), groupbit(groupbit)
{
  dthalf = dt4 = dt8 = 0;
	
  t_start = t_stop = t_period = t_current = t_target = ke_target = 0.0;
  t_freq = drag = tdrag_factor = 0;
  
  boltz = force->boltz, nktv2p = force->nktv2p;
  tdof = 0;
  mtchain = 3;
  factor_eta = 0.0;
  which = got_temp = 0;

  // argi is passed from somewhere else.
  while( argi < narg ){
    char **arg = *arg_ptr;
    if( comm->me == 0 ){
      fprintf(screen,"Parsing arg \"%s\"\n",arg[argi]);
    }
    if( is_legal_keyword( arg[argi] ) ){
      int status = parse_arg(argi, narg, arg_ptr);
      if( status ){
        char msg[256];
        sprintf(msg,"Error parsing arg \"%s\".\n", arg[argi]);
        error->all(FLERR, msg);
      }
    }else{
      ++argi;
    }
  }

	
  if( !got_temp ) error->all(FLERR,"Fix nvt/manifold/rattle needs 'temp'!");

  if( t_period < 0.0 ){
    error->all(FLERR,"Fix nvt/manifold/rattle damping parameter must be > 0.0");
  }

  // Create temperature compute:
  int n = strlen(fix_id)+6;
  id_temp = new char[n];
  strcpy(id_temp,fix_id);
  strcat(id_temp,"_temp");
  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char*) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  int icompute = modify->find_compute(id_temp);
  if( icompute < 0 ){
    error->all(FLERR,"Temperature ID for fix nvt/manifold/rattle "
               "does not exist");
  }
  temperature = modify->compute[icompute];
  if( temperature->tempbias ) which = BIAS;
  else                        which = NOBIAS;
  
  // Set t_freq from t_period
  t_freq = 1.0 / t_period;

  // Init Nosé-Hoover chain:
  eta        = new double[mtchain];
  eta_dot    = new double[mtchain+1];
  eta_dotdot = new double[mtchain];
  eta_mass   = new double[mtchain];
  eta_dot[mtchain] = 0.0;

  eta_dot[mtchain] = 0.0;
  for( int ich = 0; ich < mtchain; ++ich ){
    eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
  }
}



ManifoldThermostatNH::~ManifoldThermostatNH()
{
  // Deallocate heap-allocated objects.
  if( eta )        delete[] eta;
  if( eta_dot )    delete[] eta_dot;
  if( eta_dotdot ) delete[] eta_dotdot;
  if( eta_mass )   delete[] eta_mass;

  modify->delete_compute(id_temp);
  if( id_temp )    delete[] id_temp;

}

int ManifoldThermostatNH::parse_arg( int argi, int &narg, char ***arg_ptr )
{
  if( strcmp( (*arg_ptr)[argi], "temp") == 0 ){
    if( argi+3 >= narg )
      error->all(FLERR,"Keyword 'temp' needs 3 arguments");

    t_start  = force->numeric(FLERR,(*arg_ptr)[argi+1]);
    t_stop   = force->numeric(FLERR,(*arg_ptr)[argi+2]);
    t_period = force->numeric(FLERR,(*arg_ptr)[argi+3]);
    t_target = t_start;
    got_temp = 1;

    shift_args_back( argi, arg_ptr, 4, narg );
    
    return 0;
    
  }else if( strcmp( (*arg_ptr)[argi], "tchain" ) == 0 ){
    if( argi+1 >= narg )
      error->all(FLERR,"Keyword 'tchain' needs 1 argument");
    
    mtchain = force->inumeric(FLERR, (*arg_ptr)[argi+1]);
    shift_args_back( argi, arg_ptr, 2, narg );

    return 0;
  }
  return -1; 
}



int ManifoldThermostatNH::is_legal_keyword( const char *arg )
{
  return (strcmp(arg, "temp") == 0) || (strcmp(arg, "tchain") == 0 );
}


void ManifoldThermostatNH::init()
{
  int icompute = modify->find_compute(id_temp);
  if( icompute < 0 ){
    error->all(FLERR,"Temperature ID for fix nvt/manifold/rattle "
               "does not exist");
  }
  temperature = modify->compute[icompute];
  if( temperature->tempbias ) which = BIAS;
  else                        which = NOBIAS;
}



double ManifoldThermostatNH::memory_usage()
{
  double bytes = 0.0;

  return bytes;
}


// All of this is copied/stolen from fix_nh:
void ManifoldThermostatNH::setup(int vflag)
{
  compute_temp_target();

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  // Compute/set eta-masses:
  double inv_t_freq2 = 1.0 / (t_freq*t_freq);
  eta_mass[0] = tdof * boltz * t_target * inv_t_freq2;
  for( int ich = 1; ich < mtchain; ++ich ){
    eta_mass[ich] = boltz * t_target * inv_t_freq2;
  }
  for( int ich = 1; ich < mtchain; ++ich ){
    eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1] - 
                       boltz * t_target ) / eta_mass[ich];
  }

}


void ManifoldThermostatNH::compute_temp_target()
{
  t_current = temperature->compute_scalar();
  tdof      = temperature->dof;	

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0){
    delta /= update->endstep - update->beginstep;
  }

  tdof = temperature->dof;
  t_target = t_start + delta * (t_stop-t_start);
  ke_target = tdof * boltz * t_target;

}




void ManifoldThermostatNH::nhc_temp_integrate()
{
  int ich;
  // t_current = temperature->compute_scalar();
  // tdof = temperature->dof;
  compute_temp_target();
  
  double expfac, kecurrent = tdof * boltz * t_current;
  double inv_t_freq2 = 1.0 / (t_freq*t_freq);
  eta_mass[0] = tdof * boltz * t_target * inv_t_freq2;
  for( int ich = 1; ich < mtchain; ++ich ){
    eta_mass[ich] = boltz * t_target * inv_t_freq2;
  }

  if( eta_mass[0] > 0.0 ){
    eta_dotdot[0] = (kecurrent - ke_target)/eta_mass[0];
  }else{
    eta_dotdot[0] = 0;
  }

  for( ich = mtchain-1; ich > 0; --ich ){
    expfac = exp(-dt8*eta_dot[ich+1]);
    eta_dot[ich] *= expfac;
    eta_dot[ich] += eta_dotdot[ich] * dt4;
    eta_dot[ich] *= tdrag_factor * expfac;

  }

  expfac = exp(-dt8*eta_dot[1]);
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0] * dt4;
  eta_dot[0] *= tdrag_factor * expfac;

  factor_eta = exp(-dthalf*eta_dot[0]);

  
  if( factor_eta == 0 ){
    char msg[2048];
    sprintf(msg, "WTF, factor_eta is 0! dthalf = %f, eta_dot[0] = %f",
            dthalf, eta_dot[0]);
    error->all(FLERR,msg);
  }

  nh_v_temp();

  t_current *= factor_eta*factor_eta;
  kecurrent = tdof * boltz * t_current;
  if( eta_mass[0] > 0.0 ){
    eta_dotdot[0] = (kecurrent - ke_target) / eta_mass[0];
  }else{
    eta_dotdot[0] = 0.0;
  }

  for( int ich = 1; ich < mtchain; ++ich ){
    eta[ich] += dthalf*eta_dot[ich];
  }
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0]*dt4;
  eta_dot[0] *= expfac;

  for( int ich = 1; ich < mtchain; ++ich ){
    expfac = exp(-dt8*eta_dot[ich+1]);
    eta_dot[ich] *= expfac;
    eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1]
                       - boltz*t_target) / eta_mass[ich];
    eta_dot[ich] *= eta_dotdot[ich] * dt4;
    eta_dot[ich] *= expfac;
  }
}

void ManifoldThermostatNH::nh_v_temp()
{
  double **v = atom->v;
  int *mask  = atom->mask;
  int nlocal = atom->nlocal;
  if( igroup == atom->firstgroup) nlocal = atom->nfirst;

  if( which == NOBIAS ){
    for( int i = 0; i < nlocal; ++i ){
      if( mask[i] & groupbit ){
        v[i][0] *= factor_eta;
        v[i][1] *= factor_eta;
        v[i][2] *= factor_eta;
      }
    }
  }else if( which == BIAS ){
    for( int i = 0; i < nlocal; ++i ){
      if( mask[i] & groupbit ){
        temperature->remove_bias(i,v[i]);
        v[i][0] *= factor_eta;
        v[i][1] *= factor_eta;
        v[i][2] *= factor_eta;
        temperature->restore_bias(i,v[i]);
      }
    }
  }

}


void ManifoldThermostatNH::print_stats( const char *header )
{
  
}


void ManifoldThermostatNH::reset_dt()
{
  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  tdrag_factor = 1.0 - (update->dt * t_freq * drag);
}


void ManifoldThermostatNH::print_stuff()
{
  if( comm->me != 0 )  return;
  fprintf(screen, "\n\n********* Thermostat settings ************\n");
  fprintf(screen, "Target temperature is %f\n", t_target);
  
  fprintf(screen, "\n\n******************************************\n\n\n");
}

void ManifoldThermostatNH::update_t_current()
{
  if (which == BIAS && neighbor->ago == 0)
    t_current = temperature->compute_scalar();
}
