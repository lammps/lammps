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
   License: GNU General Public License.

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
#include "force.h"
#include "input.h"
#include "manifold_factory.h"
#include "update.h"
#include "variable.h"


#include "manifold_rattle_integrator.h"

using namespace LAMMPS_NS;

enum { CONST, EQUAL }; // For treating the variables.

inline int was_var( const char *str )
{
  if( strlen(str) > 2 ){
    return (str[0] == 'v') && (str[1] == '_');
  }else{
    return 0;
  }
}

void LAMMPS_NS::shift_args_back( int argi, char ***arg_ptr, int delta, int &narg )
{
  for( int i = argi; i < narg-delta; ++i ){
    (*arg_ptr)[i] = (*arg_ptr[i+delta]);
  }
  narg -= delta;
}



/* -------------------------------------------------------------------------
   Integrator stuff. Constructor:
   ------------------------------------------------------------------------ */
ManifoldRattleIntegratorBase::ManifoldRattleIntegratorBase(LAMMPS *lmp, int &narg,
                                                           char ***arg_ptr )
  : Pointers(lmp), tstrs(NULL), tvars(NULL), tstyle(NULL), is_var(NULL)
{
  char **arg = *arg_ptr;

  // Defaults for member variables:
  nevery = 0;
  dtv = dtf = 0;
  
  tolerance = force->numeric( FLERR, arg[3] );
  max_iter  = force->numeric( FLERR, arg[4] );

  update_style = 1;
  
  int argi = 6; // 6 for general args, + nvars for manifold vars.

	  
  while( argi < narg ){
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
}

/* -------------------------------------------------------------------------
   Destructor:
   ------------------------------------------------------------------------ */
ManifoldRattleIntegratorBase::~ManifoldRattleIntegratorBase()
{
  if( tstrs ){
    for( int i = 0; i < nvars; ++i ){
      delete [] tstrs[i];
    }
    delete [] tstrs;
  }

  if( tvars  ) delete [] tvars;
  if( tstyle ) delete [] tstyle;
  if( is_var ) delete [] is_var;
  
}


void ManifoldRattleIntegratorBase::init_manifold( int narg, char **arg,
                                                  manifold *man_ptr )
{
  if( !man_ptr ) error->all(FLERR,"manifold ptr was NULL");

  tstrs  = new char*[nvars];
  tvars  = new int[nvars];
  tstyle = new int[nvars];
  is_var = new int[nvars];

  if( !tstrs || !tvars || !tstyle || !is_var ){
    error->all(FLERR, "Error creating manifold arg arrays");
  }


  for( int i = 0; i < nvars; ++i ){
    int len = 0, offset = 0;
    if( was_var( arg[i+6] ) ){
      len = strlen(arg[i+6]) - 1; // -1 because -2 for v_, +1 for \0.
      is_var[i] = 1;
      offset = 2;
    }else{
      force->numeric(FLERR,arg[i+6]); // Check if legal number.
      len = strlen( arg[i+6] ) + 1; // +1 for \0.
      is_var[i] = 0;
    }
    tstrs[i] = new char[len];
    if( tstrs[i] == NULL ) error->all(FLERR,"Error allocating space for args.");
    strcpy( tstrs[i], arg[i+6] + offset );
  }


  *(man_ptr->get_params()) = new double[nvars];
  if( !(*man_ptr->get_params()) ) error->all(FLERR,"Failed to allocate params!");

  
  for( int i = 0; i < nvars; ++i ){
    // If param i was variable type, it will be set later...
    (*man_ptr->get_params())[i] = is_var[i] ? 0.0 : force->numeric( FLERR, arg[i+6] );
  }
  man_ptr->post_param_init();
}




void ManifoldRattleIntegratorBase::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}




// Identifies keyword in (*arg_ptr)[argi], handles it properly, and strips
// the relevant args from (*arg_ptr) and subtracts the number of args handled
// from narg (so that (*arg_ptr) contains narg elements again.
// Returns 0 on success, -1 otherwise.
int ManifoldRattleIntegratorBase::parse_arg( int argi, int &narg,
                                             char ***arg_ptr )
{
  if( strcmp((*arg_ptr)[argi], "every") == 0 ){
    if( argi + 1 >= narg ){
      error->all(FLERR,"keyword every needs one argument");
    }
    nevery = force->inumeric(FLERR,(*arg_ptr)[argi+1]);
    shift_args_back( argi, arg_ptr, 2, narg );
    
    return 0;
  }
  if( strcmp((*arg_ptr)[argi], "style") == 0 ){
    if( argi + 1 >= narg ){
      error->all(FLERR,"keyword every needs one argument");
    }
    update_style = force->inumeric(FLERR, (*arg_ptr)[argi+1]);
    if( update_style < 1 || update_style > 2 ){
      error->warning(FLERR,"update style not known, reverting to standard");
    }
    
    shift_args_back( argi, arg_ptr, 2, narg );

    return 0;
  }
  
  return -1;
}



int ManifoldRattleIntegratorBase::is_legal_keyword( const char *arg )
{
  return strcmp(arg, "every") == 0;
}









void ManifoldRattleIntegratorBase::print_stats(const char *header)
{
  double n = stats.natoms;
  if( n > 0 ){
    stats.x_iters_per_atom += stats.x_iters / n;
    stats.v_iters_per_atom += stats.v_iters / n;
  }

  double x_iters = 0, v_iters = 0;
  bigint ntimestep = update->ntimestep;
  int me = -1;


  MPI_Comm_rank(world,&me);
  MPI_Allreduce(&stats.x_iters_per_atom,&x_iters,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&stats.v_iters_per_atom,&v_iters,1,MPI_DOUBLE,MPI_SUM,world);

  // Set iters back to zero:
  stats.x_iters_per_atom = stats.x_iters = 0;
  stats.v_iters_per_atom = stats.v_iters = 0;


  if( me == 0 ){
    double inv_tdiff = 1.0/( static_cast<double>(ntimestep) - stats.last_out );
    stats.last_out = ntimestep;

    fprintf(screen, "%s stats for time step " BIGINT_FORMAT " on %d atoms:\n",
            header, ntimestep, stats.natoms);
    fprintf(screen, "  iters/atom: x = %f, v = %f, dofs removed %d",
            x_iters * inv_tdiff, v_iters * inv_tdiff, stats.dofs_removed);
    fprintf(screen,"\n");
  }
}



/* ----------------------------------------------------------------------------
   Modify the force so that the new positions remain approximately on manifold.
   Time step velocity (half-way) and position accordingly.
-------------------------------------------------------------------------------- */
void ManifoldRattleIntegratorBase::nve_x_rattle(int igroup, int groupbit)
{
  double dtfm;
  // update v and x of atoms in group
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms = 0;
  
  if (igroup == atom->firstgroup){
    nlocal = atom->nfirst;
  }
  

  if (rmass) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        natoms++;
        dtfm = dtf / rmass[i];
        rattle_manifold_x( x[i], v[i], f[i], dtv, dtfm, atom->tag[i] );
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        natoms++;
        dtfm = dtf / mass[type[i]];
        rattle_manifold_x( x[i], v[i], f[i], dtv, dtfm, atom->tag[i] );
      }
    }
  }

  if( nevery > 0 ){
    // Count ALL atoms this fix works on:
    MPI_Allreduce(&natoms,&stats.natoms,1,MPI_INT,MPI_SUM,world);
  }
}


/* ----------------------------------------------------------------------------
   Time step (other half) the velocities so that they remain in manifold.
-------------------------------------------------------------------------------- */
void ManifoldRattleIntegratorBase::nve_v_rattle(int igroup, int groupbit)
{
  double dtfm;

  // update v of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        rattle_manifold_v( v[i], f[i], x[i], dtfm );
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        rattle_manifold_v( v[i], f[i], x[i], dtfm );
      }
    }
  }
}




/******************* ManifoldRattleIntegrator stuff. *****************/



/* ----------------------------------------------------------------------
   Performs a rattle update for x
------------------------------------------------------------------------- */
void ManifoldRattleIntegrator::rattle_manifold_x(double *x, double *v,
                                                 double *f, double dtv,
                                                 double dtfm, int tagi )
{
  /*
    A RATTLE update for the position constraint.
    Original update is x += dtv * v_1/2
    Now you do
    v_1/2(lambda) = v_0 + dtfm * ( f + lambda*n_old )
    and solve
    xold - xnew + dtv * v_1/2(lambda) = 0
    g(xnew) = 0
    for x and lambda. The lambda you find then gives v_1/2 as well.
  */
  double xo[3];      // Previous position to update from.
  double vo[3];      // Previous velocity to update from.
  double l = 0;      // Lagrangian multiplier for constraint forces.
  double R[4];       // System that is 0.
  double dx[4];      // Update that follows from Newton iteration.
  double no[3];      // Normal at xo.
  double nn[3];      // Normal at x, the new position.
  double res;        // Residual.
  int iters = 0;     // Iterations used

  double c  = dtfm*dtv; // Used for iterating in the Newton loop:
  double no_nn, nn_R;

  vo[0] = v[0];
  vo[1] = v[1];
  vo[2] = v[2];
  
  xo[0] = x[0];
  xo[1] = x[1];
  xo[2] = x[2];

  double gg = ptr_m->g_and_n(x,no);
  nn[0] = no[0];
  nn[1] = no[1];
  nn[2] = no[2];

  double vt[3];
  vt[0] = vo[0] + dtfm*f[0];
  vt[1] = vo[1] + dtfm*f[1];
  vt[2] = vo[2] + dtfm*f[2];
  double no_dt[3];
  no_dt[0] = dtfm*no[0];
  no_dt[1] = dtfm*no[1];
  no_dt[2] = dtfm*no[2];

  // Assume that no_nn is roughly constant during iteration:
  
  const double c_inv = 1.0 / c;

  
  while ( 1 ) {
    v[0] = vt[0] - l*no_dt[0];
    v[1] = vt[1] - l*no_dt[1];
    v[2] = vt[2] - l*no_dt[2];

    R[0] = xo[0] - x[0] + dtv * v[0];
    R[1] = xo[1] - x[1] + dtv * v[1];
    R[2] = xo[2] - x[2] + dtv * v[2];
    R[3] = gg;

    // Analytic solution to system J*(dx,dy,dz,dl)^T = R
    // no_nn = no[0]*nn[0] + no[1]*nn[1] + no[2]*nn[2];
    nn_R  = nn[0]*R[0]  + nn[1]*R[1]  + nn[2]*R[2];
    no_nn = no[0]*nn[0] + no[1]*nn[1] + no[2]*nn[2];
    double n_inv = 1.0 / no_nn;

    // fprintf( screen, "nn_R = %f, no_nn = %f\n", nn_R, no_nn );

    dx[3] = -nn_R - R[3];
    dx[3] *= n_inv;
    dx[0] = -R[0] - no[0]*dx[3];
    dx[1] = -R[1] - no[1]*dx[3];
    dx[2] = -R[2] - no[2]*dx[3];

    dx[3] *= c_inv;


    double smooth = 1.0; // / sqrt( 1.0 + res*res );
    x[0] -= smooth*dx[0];
    x[1] -= smooth*dx[1];
    x[2] -= smooth*dx[2];
    l    -= smooth*dx[3];

    /*
    fprintf( screen, "dx = (%e, %e, %e, %e), R = (%e, %e, %e, %e)\n",
             dx[0], dx[1], dx[2], dx[3], R[0], R[1], R[2], R[3] );
    fprintf( screen, "x  = (%f, %f, %f, %f), n = (%f, %f, %f)\n",
             x[0], x[1], x[2], l, nn[0], nn[1], nn[2] );
    */
    
    res = infnorm<4>(R);
    ++iters;

    if( (res < tolerance) || (iters >= max_iter) ) break;
    
    // Update nn and g.
    gg = ptr_m->g(x);
    ptr_m->n(x,nn);
    // gg = ptr_m->g(x);
  }

  if( iters >= max_iter && res > tolerance ){
    char msg[2048];
    sprintf(msg,"Failed to constrain atom %d (x = (%f, %f, %f)! res = %e, iters = %d\n",
            tagi, x[0], x[1], x[2], res, iters);
    error->one(FLERR,msg);
  }
  
  // "sync" x and v:
  v[0] = vt[0] - l*no_dt[0];
  v[1] = vt[1] - l*no_dt[1];
  v[2] = vt[2] - l*no_dt[2];

  stats.x_iters += iters;
}


/* ----------------------------------------------------------------------
   Performs a rattle update for v
------------------------------------------------------------------------- */
void ManifoldRattleIntegrator::rattle_manifold_v(double *v, double *f,
                                                 double *x, double dtfm,
                                                 int tagi )
{
  /*
    The original update was
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
    
    Now you add the rattle-like update:
    vold - vnew + dtfm * F + mu * n_new = 0
    dot( vnew, n_new ) = 0
  */
  double vo[3];      // V at t + 1/2 dt
  double l = 0;      // Lagrangian multiplier for constraint forces.
  double R[4];       // System that is 0.
  double dv[4];      // Update that follows from Newton iteration.
  double n[3];       // Normal.
  double res;        // Residual.
  int iters = 0;     // Iterations used

  double c  = dtfm; // Used for iterating in the Newton loop:
  double nn2, nn_R;

  vo[0] = v[0];
  vo[1] = v[1];
  vo[2] = v[2];

  // Initial guess is unconstrained update:
  v[0] += dtfm*f[0];
  v[1] += dtfm*f[1];
  v[2] += dtfm*f[2];

  ptr_m->n(x,n);

  double vt[3];
  vt[0] = vo[0] + dtfm*f[0];
  vt[1] = vo[1] + dtfm*f[1];
  vt[2] = vo[2] + dtfm*f[2];
  double no_dt[3];
  no_dt[0] = dtfm*n[0];
  no_dt[1] = dtfm*n[1];
  no_dt[2] = dtfm*n[2];

  nn2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  const double n_inv = 1.0 / nn2;
  const double c_inv = 1.0 / c;
  
  switch( update_style ){
    default:
    case 1:
      do{
        R[0] = vt[0] - v[0]  - l * no_dt[0];
        R[1] = vt[1] - v[1]  - l * no_dt[1];
        R[2] = vt[2] - v[2]  - l * no_dt[2];
        R[3] = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];

        // Analytic solution to system J*(dx,dy,dz,dl)^T = R
        nn_R = n[0]*R[0] + n[1]*R[1] + n[2]*R[2];
    
        dv[3] = -nn_R - R[3];
        dv[3] *= n_inv;
        dv[0] = -n[0]*dv[3] - R[0];
        dv[1] = -n[1]*dv[3] - R[1];
        dv[2] = -n[2]*dv[3] - R[2];
        dv[3] *= c_inv;
        
        v[0] -= dv[0];
        v[1] -= dv[1];
        v[2] -= dv[2];
        l    -= dv[3];
    
        res = infnorm<4>(R);
        ++iters;
      }while( (res > tolerance) && (iters < max_iter) );

      if( iters >= max_iter && res >= tolerance ){
        char msg[2048];
        sprintf(msg,"Failed to constrain atom %d (x = (%f, %f, %f)! res = %e, iters = %d\n",
                tagi, x[0], x[1], x[2], res, iters);
        error->all(FLERR,msg);
      }
      break;

    case 2:{
      double mass = dtfm / update->dt;
      double vhnp = vo[0]*n[0] + vo[1]*n[1] + vo[2]*n[2];
      double fpnp = f[0]*n[0] + f[1]*n[1] + f[2]*n[2];
      double phnp = dtfm*vhnp / update->dt;

      l = (2.0*mass/update->dt)*phnp + fpnp;
      v[0] += 0.5*dtv*f[0] - n[0]*(phnp + fpnp);
      v[1] += 0.5*dtv*f[1] - n[1]*(phnp + fpnp);
      v[2] += 0.5*dtv*f[2] - n[2]*(phnp + fpnp);
      break;
    }
  }

  stats.v_iters += iters;
}





int ManifoldRattleIntegratorBase::dof( int igroup, int groupbit )
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms = 0;
  for( int i = 0; i < nlocal; ++i ){
    if(mask[i] & groupbit) ++natoms;
  }

  int dofs;
  MPI_Allreduce( &natoms, &dofs, 1, MPI_INT, MPI_SUM, world );

  // Make sure that, if there is just no or one atom, no dofs are subtracted,
  // since for the first atom already 3 dofs are subtracted because of the
  // centre of mass corrections:
  if( dofs <= 1 ) dofs = 0;
  stats.dofs_removed = dofs;
  
  return dofs;
}



ManifoldRattleIntegrator::ManifoldRattleIntegrator( LAMMPS *lmp, int &narg,
                                                    char ***arg_ptr ) :
  ManifoldRattleIntegratorBase( lmp, narg, arg_ptr ), ptr_m(NULL)
{
  char **arg = *arg_ptr;
  ptr_m = create_manifold(arg[5], lmp, narg, arg);
  if( !ptr_m ){
    char msg[2048];
    sprintf(msg,"Could not create manifold of type \"%s\"", arg[5]);
    error->all(FLERR,msg);
  }
  
  nvars = ptr_m->nparams();
  
  init_manifold( narg, arg, ptr_m );
}

ManifoldRattleIntegrator::~ManifoldRattleIntegrator()
{
  if( ptr_m ) delete ptr_m;
}





void ManifoldRattleIntegrator::init()
{
  // Makes sure the manifold params are set initially.
  update_var_params();

  reset_dt();
}

void ManifoldRattleIntegrator::update_params()
{
  update_var_params();
}

void ManifoldRattleIntegrator::update_var_params()
{
  if( nevery > 0 ){
    stats.x_iters = 0;
    stats.v_iters = 0;
    stats.natoms  = 0;
    stats.x_iters_per_atom = 0.0;
    stats.v_iters_per_atom = 0.0;
  }
	
  double **ptr_params = ptr_m->get_params();
  for( int i = 0; i < nvars; ++i ){
    if( is_var[i] ){
      tvars[i] = input->variable->find(tstrs[i]);
      if( tvars[i] < 0 ){
        error->all(FLERR,
                   "Variable name for fix nve/manifold/rattle does not exist");
      }
      if( input->variable->equalstyle(tvars[i]) ){
        tstyle[i] = EQUAL;
        double new_val = input->variable->compute_equal(tvars[i]);
        // fprintf( stdout, "New value of var %d is now %f\n", i+1, new_val );
        *(ptr_params[i]) = new_val;
      }else{
        error->all(FLERR,
                   "Variable for fix nve/manifold/rattle is invalid style");
      }
    }
  }
}


double ManifoldRattleIntegrator::memory_usage()
{
  double bytes = 0.0;

  bytes += sizeof(statistics);
  bytes += sizeof(*ptr_m) + sizeof(ptr_m);
  bytes += nvars*sizeof(double) + sizeof(double*);
  bytes += nvars*( sizeof(char*) + 3*sizeof(int) );
  return bytes;
}


