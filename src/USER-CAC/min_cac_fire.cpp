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

#include <cmath>
#include "min_cac_fire.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5
#define DT_GROW 1.1
#define DT_SHRINK 0.5
#define ALPHA0 0.1
#define ALPHA_SHRINK 0.99
#define TMAX 10.0

/* ---------------------------------------------------------------------- */

CACMinFire::CACMinFire(LAMMPS *lmp) : Min(lmp) {
  copy_flag=1;
}

/* ---------------------------------------------------------------------- */

void CACMinFire::init()
{
  Min::init();
  if (!atom->CAC_flag) error->all(FLERR,"CAC min styles require a CAC atom style");
  if (!atom->CAC_pair_flag) error->all(FLERR,"CAC min styles require a CAC pair style");
  dt = update->dt;
  dtmax = TMAX * dt;
  alpha = ALPHA0;
  last_negative = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void CACMinFire::setup_style()
{

  int *npoly = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int nodes_per_element;

  double *min_v = atom->min_v;
  nvec=atom->dense_count;

  for (int i=0; i < nvec; i ++) min_v[i] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void CACMinFire::reset_vectors()
{
  double *min_x = atom->min_x;
  double *min_f = atom->min_f;
  nvec=atom->dense_count;

  if (nvec) xvec = min_x;
  if (nvec) fvec = min_f;

}

/* ---------------------------------------------------------------------- */

int CACMinFire::iterate(int maxiter)
{
  bigint ntimestep;
  double vmax,vdotf,vdotfall,vdotv,vdotvall,fdotf,fdotfall;
  double scale1,scale2;
  double dtvone,dtv,dtf,dtfm;
  int i, flag,flagall;

  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int *npoly = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;

  alpha_final = 0.0;

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // vdotfall = v dot f

    double *f = atom->min_f;
    double *v = atom->min_v;
    nvec=atom->dense_count;

    vdotf = 0.0;
    for (i = 0; i < nvec; i+=3) {
      vdotf += v[i]*f[i] + v[i+1]*f[i+1] + v[i+2]*f[i+2];
    }
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // sum vdotf over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      vdotf = vdotfall;
      MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // if more than DELAYSTEP since v dot f was negative:
    // increase timestep and decrease alpha

    if (vdotfall > 0.0) {
      vdotv = 0.0;
      for (i = 0; i < nvec; i+=3) {
        vdotv += v[i]*v[i] + v[i+1]*v[i+1] + v[i+2]*v[i+2];
      }
      MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum vdotv over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        vdotv = vdotvall;
        MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      fdotf = 0.0;
      for (i = 0; i < nvec; i+=3) {
        fdotf += f[i]*f[i] + f[i+1]*f[i+1] + f[i+2]*f[i+2];
      }
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum fdotf over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        fdotf = fdotfall;
        MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      scale1 = 1.0 - alpha;
      if (fdotfall == 0.0) scale2 = 0.0;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);
      for (i = 0; i < nvec; i+=3) {
        v[i] =  scale1*v[i] + scale2*f[i];
        v[i+1] =  scale1*v[i+1] + scale2*f[i+1];
        v[i+2] =  scale1*v[i+2] + scale2*f[i+2];
      }

      if (ntimestep - last_negative > DELAYSTEP) {
        dt = MIN(dt*DT_GROW,dtmax);
        alpha *= ALPHA_SHRINK;
      }

    // else (v dot f) <= 0:
    // decrease timestep, reset alpha, set v = 0

    } else {
      last_negative = ntimestep;
      dt *= DT_SHRINK;
      alpha = ALPHA0;
      for (i = 0; i < nvec; i+=3) {
        v[i] =  v[i+1] = v[i+2] = 0.0;
      }
    }

    // limit timestep so no particle moves further than dmax

    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;

    dtvone = dt;

    for (int i = 0; i < nvec; i+=3) {
      vmax = MAX(fabs(v[i]),fabs(v[i+1]));
      vmax = MAX(vmax,fabs(v[i+2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    // min dtv over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      dtvone = dtv;
      MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,universe->uworld);
    }

    dtf = dtv * force->ftm2v;

    // Euler integration step

    double *x = atom->min_x;

    if (rmass) {
      for (int i = 0; i < nvec; i++) {
        dtfm = dtf / rmass[i];
        x[i] += dtv * v[i];
        v[i] += dtfm * f[i];

      }
    } else {

        int dense = 0;
        // nodal loops required to get mass
        for(int element_counter=0; element_counter < atom->nlocal; element_counter++) {
          for (int poly_counter = 0; poly_counter < npoly[element_counter]; poly_counter++) {
            dtfm = dtf / mass[node_types[element_counter][poly_counter]];
            for(int node_counter=0; node_counter < nodes_per_element_list[element_type[element_counter]]; node_counter++){

              x[dense+0] += dtv * v[dense+0];
              x[dense+1] += dtv * v[dense+1];
              x[dense+2] += dtv * v[dense+2];
              v[dense+0] += dtfm * f[dense+0];
              v[dense+1] += dtfm * f[dense+1];
              v[dense+2] += dtfm * f[dense+2];

              dense+=3;
            }
          }
        }
    }


    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    // energy tolerance criterion
    // only check after DELAYSTEP elapsed since velocties reset to 0
    // sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          return ETOL;
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    // force tolerance criterion
    // sync across replicas if running multi-replica minimization

    if (update->ftol > 0.0) {
      fdotf = fnorm_sqr();
      if (update->multireplica == 0) {
        if (fdotf < update->ftol*update->ftol) return FTOL;
      } else {
        if (fdotf < update->ftol*update->ftol) flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}


/* ----------------------------------------------------------------------
   copy dense arrays to atomvec arrays for energy_force evaluation
------------------------------------------------------------------------- */

void CACMinFire::copy_vectors(){
int *npoly = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_forces = atom->nodal_forces;
  double ****nodal_velocities = atom->nodal_velocities;
  double *min_x = atom->min_x;
  double *min_f = atom->min_f;
  double *min_v = atom->min_v;
  double **x = atom->x;
  int nodes_per_element;



  //copy contents to these vectors
  int dense_count_x=0;
  int dense_count_f=0;
  int dense_count_v=0;
  for(int element_counter=0; element_counter < atom->nlocal; element_counter++){
    for(int poly_counter=0; poly_counter < npoly[element_counter]; poly_counter++){
      for(int node_counter=0; node_counter < nodes_per_element_list[element_type[element_counter]]; node_counter++){
         nodal_positions[element_counter][poly_counter][node_counter][0] = min_x[dense_count_x++];
         nodal_positions[element_counter][poly_counter][node_counter][1] = min_x[dense_count_x++];
         nodal_positions[element_counter][poly_counter][node_counter][2] = min_x[dense_count_x++];
         nodal_forces[element_counter][poly_counter][node_counter][0] = min_f[dense_count_f++];
         nodal_forces[element_counter][poly_counter][node_counter][1] = min_f[dense_count_f++];
         nodal_forces[element_counter][poly_counter][node_counter][2] = min_f[dense_count_f++];
         nodal_velocities[element_counter][poly_counter][node_counter][0] = min_v[dense_count_v++];
         nodal_velocities[element_counter][poly_counter][node_counter][1] = min_v[dense_count_v++];
         nodal_velocities[element_counter][poly_counter][node_counter][2] = min_v[dense_count_v++];
       }
     }
  }

    // update x for elements and atoms using nodal variables
  for (int i = 0; i < atom->nlocal; i++){
    //determine element type

    nodes_per_element=nodes_per_element_list[element_type[i]];
    x[i][0] = 0;
    x[i][1] = 0;
    x[i][2] = 0;

    for (int poly_counter = 0; poly_counter < npoly[i];poly_counter++) {
      for(int k=0; k<nodes_per_element; k++){
        x[i][0] += nodal_positions[i][poly_counter][k][0];
        x[i][1] += nodal_positions[i][poly_counter][k][1];
        x[i][2] += nodal_positions[i][poly_counter][k][2];
      }
    }
  x[i][0] = x[i][0] / nodes_per_element / npoly[i];
  x[i][1] = x[i][1] / nodes_per_element / npoly[i];
  x[i][2] = x[i][2] / nodes_per_element / npoly[i];
  }

}

