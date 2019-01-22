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


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "group.h"
#include <cmath>
#include "input.h"
#include "variable.h"
#include "citeme.h"
#include "memory.h"
#include "comm.h"


#include "fix_nve_manifold_rattle.h"
#include "manifold_factory.h"
#include "manifold.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace user_manifold;


enum { CONST, EQUAL }; // For treating the variables.


static const char* cite_fix_nve_manifold_rattle =
  "fix nve/manifold/rattle command:\n\n"
  "@article{paquay-2016,\n"
  "   author        = {Paquay, Stefan and Kusters, Remy},\n"
  "   doi           = {10.1016/j.bpj.2016.02.017},\n"
  "   issn          = {0006-3495},\n"
  "   journal       = {Biophysical Journal},\n"
  "   month         = apr,\n"
  "   number        = {6},\n"
  "   pages         = {1226--1233},\n"
  "   title         = {{A Method for Molecular Dynamics on Curved Surfaces}},\n"
  "   volume        = {110},\n"
  "   year          = {2016}\n"
  "}\n\n";


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
FixNVEManifoldRattle::FixNVEManifoldRattle( LAMMPS *lmp, int &narg, char **arg,
                                            int error_on_unknown_keyword )
  : Fix(lmp,narg,arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_nve_manifold_rattle);
  if (narg < 6 ) error->all(FLERR, "Illegal fix nve/manifold/rattle command");

  // Set all bits/settings:
  time_integrate = 1;
  dynamic_group_allow = 1;
  size_vector = 0;
  dof_flag = 1;

  nevery = 0;
  next_output = 0;
  dtv = dtf = 0;

  tolerance = force->numeric( FLERR, arg[3] );
  max_iter  = force->numeric( FLERR, arg[4] );

  ptr_m = create_manifold(arg[5], lmp, narg, arg);
  if (!ptr_m) {
    error->all(FLERR,"Error creating manifold pointer");
  }

  nvars = ptr_m->nparams();
  tstrs  = new char*[nvars];
  tvars  = new int[nvars];
  tstyle = new int[nvars];
  is_var = new int[nvars];

  if (!tstrs || !tvars || !tstyle || !is_var) {
    error->all(FLERR, "Error creating manifold arg arrays");
  }

  // Check if you have enough args:
  if (6 + nvars > narg) {
    char msg[2048];
    sprintf(msg, "Not enough args for manifold %s, %d expected but got %d\n",
            ptr_m->id(), nvars, narg - 6);
    error->all(FLERR, msg);
  }
  // Loop over manifold args:
  for( int i = 0; i < nvars; ++i ){
    int len = 0, offset = 0;
    if (was_var( arg[i+6] )) {
      len = strlen(arg[i+6]) - 1; // -1 because -2 for v_, +1 for \0.
      is_var[i] = 1;
      offset = 2;
    } else {
      force->numeric(FLERR,arg[i+6]); // Check if legal number.
      len = strlen( arg[i+6] ) + 1; // +1 for \0.
      is_var[i] = 0;
    }
    tstrs[i] = new char[len];
    if (tstrs[i] == NULL ) error->all(FLERR,"Error allocating space for args.");
    strcpy( tstrs[i], arg[i+6] + offset );
  }

  ptr_m->params = new double[nvars];
  if (!ptr_m->params ) error->all(FLERR,"Failed to allocate params!");
  for( int i = 0; i < nvars; ++i ){
    // If param i was variable type, it will be set later...
    ptr_m->params[i] = is_var[i] ? 0.0 : force->numeric( FLERR, arg[i+6] );
  }
  ptr_m->post_param_init();


  // Loop over rest of args:
  int argi = 6 + nvars;
  while( argi < narg ){
    if (strcmp(arg[argi], "every") == 0) {
      nevery = force->inumeric(FLERR,arg[argi+1]);
      next_output = update->ntimestep + nevery;
      if (comm->me == 0) {
        fprintf(screen,"Outputing every %d steps, next is %d\n",
                        nevery, next_output);
      }
      argi += 2;
    } else if (error_on_unknown_keyword) {
      char msg[2048];
      sprintf(msg,"Error parsing arg \"%s\".\n", arg[argi]);
      error->all(FLERR, msg);
    } else {
      argi += 1;
    }
  }

}


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
FixNVEManifoldRattle::~FixNVEManifoldRattle()
{
  if (tstrs) {
    for( int i = 0; i < nvars; ++i ){
      delete [] tstrs[i];
    }
    delete [] tstrs;
  }

  if (tvars ) delete [] tvars;
  if (tstyle) delete [] tstyle;
  if (is_var) delete [] is_var;
}




/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;



}

void FixNVEManifoldRattle::print_stats( const char *header )
{
  double n = stats.natoms;
  if (n > 0) {
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


  if (me == 0) {
    double inv_tdiff = 1.0/( static_cast<double>(ntimestep) - stats.last_out );
    stats.last_out = ntimestep;

    fprintf(screen, "%s stats for time step " BIGINT_FORMAT " on %d atoms:\n",
            header, ntimestep, stats.natoms);
    fprintf(screen, "  iters/atom: x = %f, v = %f, dofs removed %d",
            x_iters * inv_tdiff, v_iters * inv_tdiff, stats.dofs_removed);
    fprintf(screen,"\n");
  }

  stats.x_iters_per_atom = 0;
  stats.v_iters_per_atom = 0;
  stats.x_iters = 0;
  stats.v_iters = 0;
}


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
int FixNVEManifoldRattle::was_var( const char *str )
{
  if (strlen(str) > 2) {
    return (str[0] == 'v') && (str[1] == '_');
  } else {
    return 0;
  }
}



/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
int FixNVEManifoldRattle::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  if (nevery > 0) mask |= END_OF_STEP;

  return mask;
}


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::init()
{
  // Makes sure the manifold params are set initially.

  update_var_params();
  reset_dt();
}


void FixNVEManifoldRattle::update_var_params()
{

  double *ptr_params = ptr_m->params;

  for( int i = 0; i < nvars; ++i ){
    if (is_var[i]) {
      tvars[i] = input->variable->find(tstrs[i]);
      if (tvars[i] < 0) {
        error->all(FLERR,
                   "Variable name for fix nve/manifold/rattle does not exist");
      }
      if (input->variable->equalstyle(tvars[i])) {
        tstyle[i] = EQUAL;
        double new_val = input->variable->compute_equal(tvars[i]);

        ptr_params[i] = new_val;
      } else {
        error->all(FLERR,
                   "Variable for fix nve/manifold/rattle is invalid style");
      }
    }
  }
}



/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
int FixNVEManifoldRattle::dof(int /*igroup*/)
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
  if (dofs <= 1) dofs = 0;
  stats.dofs_removed = dofs;

  return dofs;
}



/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
double FixNVEManifoldRattle::memory_usage()
{
  double bytes = 0.0;

  bytes += sizeof(statistics);
  bytes += sizeof(*ptr_m) + sizeof(ptr_m);
  bytes += nvars*sizeof(double) + sizeof(double*);
  bytes += nvars*( sizeof(char*) + 3*sizeof(int) );
  return bytes;
}




/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::initial_integrate(int /*vflag*/)
{
  update_var_params();
  nve_x_rattle(igroup, groupbit);
}


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::final_integrate()
{
  nve_v_rattle(igroup, groupbit);
}



/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::end_of_step()
{
  if (nevery && (update->ntimestep == next_output)){
    if (comm->me == 0) {
      print_stats( "nve/manifold/rattle" );
      next_output += nevery;
    }
  }
}

/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::nve_x_rattle(int igroup, int groupbit)
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

  if (nevery > 0) {
    // Count ALL atoms this fix works on:
    MPI_Allreduce(&natoms,&stats.natoms,1,MPI_INT,MPI_SUM,world);
  }
}

/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::nve_v_rattle(int igroup, int groupbit)
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
        rattle_manifold_v( v[i], f[i], x[i], dtfm, atom->tag[i] );
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        rattle_manifold_v( v[i], f[i], x[i], dtfm, atom->tag[i] );
      }
    }
  }
}


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::rattle_manifold_x(double *x, double *v,
                                             double *f, double dtv,
                                             double dtfm, tagint tagi )
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


    x[0] -= dx[0];
    x[1] -= dx[1];
    x[2] -= dx[2];
    l    -= dx[3];

    res = infnorm<4>(R);
    ++iters;

    if ((res < tolerance) || (iters >= max_iter)) break;

    // Update nn and g.
    gg = ptr_m->g(x);
    ptr_m->n(x,nn);
    // gg = ptr_m->g(x);
  }

  if (iters >= max_iter && res > tolerance) {
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


/* -----------------------------------------------------------------------------
   ---------------------------------------------------------------------------*/
void FixNVEManifoldRattle::rattle_manifold_v(double *v, double *f,
                                             double *x, double dtfm,
                                             tagint tagi )
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

  if (iters >= max_iter && res >= tolerance) {
          char msg[2048];
          sprintf(msg,"Failed to constrain atom %d (x = (%f, %f, %f)! res = %e, iters = %d\n",
                  tagi, x[0], x[1], x[2], res, iters);
          error->all(FLERR,msg);
  }

  stats.v_iters += iters;
}
