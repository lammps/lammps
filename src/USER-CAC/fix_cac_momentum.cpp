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

#include <cstdlib>
#include <cstring>
#include "fix_cac_momentum.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCACMomentum::FixCACMomentum(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix cac/momentum command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix cac/momentum command");

  dynamic = linear = angular = rescale = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix cac/momentum command");
      linear = 1;
      xflag = force->inumeric(FLERR,arg[iarg+1]);
      yflag = force->inumeric(FLERR,arg[iarg+2]);
      zflag = force->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angular") == 0) {
      angular = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"rescale") == 0) {
      rescale = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix cac/momentum command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR,"Illegal fix cac/momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR,"Illegal fix cac/momentum command");

  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

int FixCACMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCACMomentum::init()
{
  if (!atom->CAC_flag) error->all(FLERR,"fix cac/momentum requires a CAC atom style");
  if (group->dynamic[igroup]) {
    dynamic = 1;
  } else {
   if (group->count(igroup) == 0)
     error->all(FLERR,"Fix cac/momentum group has no atoms or elements");
  }
  quadrature_init();
}

/* ---------------------------------------------------------------------- */

void FixCACMomentum::end_of_step()
{
  double **x = atom->x;
  double **v = atom->v;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ekin_old,ekin_new;
  ekin_old = ekin_new = 0.0;
  
  //compute system mass
  masstotal = 0;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int **node_types = atom->node_types;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **element_scale = atom->element_scale; 
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;	
  double iso_volume;
  double unit_cell_mapped[3];

  double one = 0.0;

  for (int i = 0; i < nlocal; i++)
  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
    if (mask[i] & groupbit){
      if(element_type[i]==1)
      one += mass[node_types[i][ipoly]]*element_scale[i][0]*element_scale[i][1]*element_scale[i][2];
      if(element_type[i]==0)
      one += mass[node_types[i][ipoly]];
    }

  MPI_Allreduce(&one,&masstotal,1,MPI_DOUBLE,MPI_SUM,world);
  // do nothing is group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;
  
  if (linear) {
    double vcm[3];
    //compute center of mass velocity
    double p[3],massone;
    p[0] = p[1] = p[2] = 0.0;
    double shape_func;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(element_type[i]==1){
        nodes_per_element = nodes_count_list[element_type[i]];
        unit_cell_mapped[0] = 2 / double(element_scale[i][0]);
        unit_cell_mapped[1] = 2 / double(element_scale[i][1]);
        unit_cell_mapped[2] = 2 / double(element_scale[i][2]);
        iso_volume=unit_cell_mapped[0]*unit_cell_mapped[1]*unit_cell_mapped[2];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
        for(int ii=0;ii<nodes_per_element;ii++){
          for (int qi=0; qi< quadrature_node_count;qi++){
            for (int qj=0; qj< quadrature_node_count;qj++){
              for (int qk=0; qk< quadrature_node_count;qk++){
              shape_func= mass[node_types[i][ipoly]]*
              quadrature_weights[qi] * quadrature_weights[qj] * quadrature_weights[qk] *
              shape_function(quadrature_abcissae[qi],quadrature_abcissae[qj],quadrature_abcissae[qk],2,ii+1)/iso_volume;
              p[0]+=shape_func*nodal_velocities[i][ipoly][ii][0];
              p[1]+=shape_func*nodal_velocities[i][ipoly][ii][1];
              p[2]+=shape_func*nodal_velocities[i][ipoly][ii][2];
              }
            }
          }
        }
        }
        else if(element_type[i]==0){
        massone = mass[node_types[i][0]];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
        }
      }
  

  MPI_Allreduce(p,vcm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    vcm[0] /= masstotal;
    vcm[1] /= masstotal;
    vcm[2] /= masstotal;
  }

    //group->vcm(igroup,masstotal,vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
        for(int ii=0;ii<nodes_per_element;ii++){
        if (xflag) nodal_velocities[i][ipoly][ii][0] -= vcm[0];
        if (yflag) nodal_velocities[i][ipoly][ii][1] -= vcm[1];
        if (zflag) nodal_velocities[i][ipoly][ii][2] -= vcm[2];
        }
      }
  }

  if (angular) {
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    //group->xcm(igroup,masstotal,xcm);
    //compute center of mass
    double c[3],massone,nodesum[3],nodevel[3];
    double dx,dy,dz;
    c[0] = c[1] = c[2] = 0.0;
    double shape_func;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(element_type[i]==1){
        nodes_per_element = nodes_count_list[element_type[i]];
        unit_cell_mapped[0] = 2 / double(element_scale[i][0]);
        unit_cell_mapped[1] = 2 / double(element_scale[i][1]);
        unit_cell_mapped[2] = 2 / double(element_scale[i][2]);
        iso_volume=unit_cell_mapped[0]*unit_cell_mapped[1]*unit_cell_mapped[2];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
        for(int ii=0;ii<nodes_per_element;ii++){
          for (int qi=0; qi< quadrature_node_count;qi++){
            for (int qj=0; qj< quadrature_node_count;qj++){
              for (int qk=0; qk< quadrature_node_count;qk++){
              shape_func= mass[node_types[i][ipoly]]*
              quadrature_weights[qi] * quadrature_weights[qj] * quadrature_weights[qk] *
              shape_function(quadrature_abcissae[qi],quadrature_abcissae[qj],quadrature_abcissae[qk],2,ii+1)/iso_volume;
              c[0]+=shape_func*nodal_positions[i][ipoly][ii][0];
              c[1]+=shape_func*nodal_positions[i][ipoly][ii][1];
              c[2]+=shape_func*nodal_positions[i][ipoly][ii][2];
              }
            }
          }
        }
        }
        else if(element_type[i]==0){
        massone = mass[node_types[i][0]];
        c[0] += x[i][0]*massone;
        c[1] += x[i][1]*massone;
        c[2] += x[i][2]*massone;
        }
      }
  

  MPI_Allreduce(c,xcm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    xcm[0] /= masstotal;
    xcm[1] /= masstotal;
    xcm[2] /= masstotal;
  }
    //compute angular momentum
    c[0] = c[1] = c[2] = 0.0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(element_type[i]==1){
        nodes_per_element = nodes_count_list[element_type[i]];
        unit_cell_mapped[0] = 2 / double(element_scale[i][0]);
        unit_cell_mapped[1] = 2 / double(element_scale[i][1]);
        unit_cell_mapped[2] = 2 / double(element_scale[i][2]);
        iso_volume=unit_cell_mapped[0]*unit_cell_mapped[1]*unit_cell_mapped[2];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
          for (int qi=0; qi< quadrature_node_count;qi++){
            for (int qj=0; qj< quadrature_node_count;qj++){
              for (int qk=0; qk< quadrature_node_count;qk++){
              nodesum[0] = nodesum[1] = nodesum[2] = nodevel[0] = nodevel[1] = nodevel[2] =0;
              for(int ii=0;ii<nodes_per_element;ii++){
                shape_func=shape_function(quadrature_abcissae[qi],quadrature_abcissae[qj],quadrature_abcissae[qk],2,ii+1);
                nodesum[0]+=shape_func*nodal_positions[i][ipoly][ii][0];
                nodesum[1]+=shape_func*nodal_positions[i][ipoly][ii][1];
                nodesum[2]+=shape_func*nodal_positions[i][ipoly][ii][2];
              }  
              for(int ii=0;ii<nodes_per_element;ii++){
                shape_func=shape_function(quadrature_abcissae[qi],quadrature_abcissae[qj],quadrature_abcissae[qk],2,ii+1);
                nodevel[0]+=shape_func*nodal_velocities[i][ipoly][ii][0];
                nodevel[1]+=shape_func*nodal_velocities[i][ipoly][ii][1];
                nodevel[2]+=shape_func*nodal_velocities[i][ipoly][ii][2];
              }
              dx = nodesum[0]-xcm[0];
              dy = nodesum[1]-xcm[1];
              dz = nodesum[2]-xcm[2];
              shape_func= mass[node_types[i][ipoly]]*
              quadrature_weights[qi] * quadrature_weights[qj] * quadrature_weights[qk]/iso_volume;
              c[0] += shape_func*(dy*nodevel[2] - dz*nodevel[1]);
              c[1] += shape_func*(dz*nodevel[0] - dx*nodevel[2]);
              c[2] += shape_func*(dx*nodevel[1] - dy*nodevel[0]);
              }
            }
          }
        }
        else if(element_type[i]==0){
        massone = mass[node_types[i][0]];
        dx = x[i][0]-xcm[0];
        dy = x[i][1]-xcm[1];
        dz = x[i][2]-xcm[2];
        c[0] += massone * (dy*v[i][2] - dz*v[i][1]);
        c[1] += massone * (dz*v[i][0] - dx*v[i][2]);
        c[2] += massone * (dx*v[i][1] - dy*v[i][0]);
        }
      }
  

  MPI_Allreduce(c,angmom,3,MPI_DOUBLE,MPI_SUM,world);
    //group->inertia(igroup,xcm,inertia);
    //compute moment of inertia matrix
  double ione[3][3];
  for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
      ione[ii][jj] = 0.0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(element_type[i]==1){
        nodes_per_element = nodes_count_list[element_type[i]];
        unit_cell_mapped[0] = 2 / double(element_scale[i][0]);
        unit_cell_mapped[1] = 2 / double(element_scale[i][1]);
        unit_cell_mapped[2] = 2 / double(element_scale[i][2]);
        iso_volume=unit_cell_mapped[0]*unit_cell_mapped[1]*unit_cell_mapped[2];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
          for (int qi=0; qi< quadrature_node_count;qi++){
            for (int qj=0; qj< quadrature_node_count;qj++){
              for (int qk=0; qk< quadrature_node_count;qk++){
              nodesum[0] = nodesum[1] = nodesum[2] = 0;
              for(int ii=0;ii<nodes_per_element;ii++){
                shape_func=shape_function(quadrature_abcissae[qi],quadrature_abcissae[qj],quadrature_abcissae[qk],2,ii+1);
                nodesum[0]+=shape_func*nodal_positions[i][ipoly][ii][0];
                nodesum[1]+=shape_func*nodal_positions[i][ipoly][ii][1];
                nodesum[2]+=shape_func*nodal_positions[i][ipoly][ii][2];
              } 
              dx = nodesum[0]-xcm[0];
              dy = nodesum[1]-xcm[1];
              dz = nodesum[2]-xcm[2];
              shape_func= mass[node_types[i][ipoly]]*
              quadrature_weights[qi] * quadrature_weights[qj] * quadrature_weights[qk]/iso_volume;
              ione[0][0] += shape_func * (dy*dy + dz*dz);
              ione[1][1] += shape_func * (dx*dx + dz*dz);
              ione[2][2] += shape_func * (dx*dx + dy*dy);
              ione[0][1] -= shape_func * dx*dy;
              ione[1][2] -= shape_func * dy*dz;
              ione[0][2] -= shape_func * dx*dz;
              }
            }
          }
        }
        else if(element_type[i]==0){
        massone = mass[node_types[i][0]];
        dx = x[i][0]-xcm[0];
        dy = x[i][1]-xcm[1];
        dz = x[i][2]-xcm[2];
        ione[0][0] += massone * (dy*dy + dz*dz);
        ione[1][1] += massone * (dx*dx + dz*dz);
        ione[2][2] += massone * (dx*dx + dy*dy);
        ione[0][1] -= massone * dx*dy;
        ione[1][2] -= massone * dy*dz;
        ione[0][2] -= massone * dx*dz;
        }
      }

    ione[1][0] = ione[0][1];
    ione[2][1] = ione[1][2];
    ione[2][0] = ione[0][2];

    MPI_Allreduce(&ione[0][0],&inertia[0][0],9,MPI_DOUBLE,MPI_SUM,world);
    //compute angular rates
    group->omega(angmom,inertia,omega);
    

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
        for(int ipoly=0; ipoly < poly_count[i]; ipoly++)  
        for(int ii=0;ii<nodes_per_element;ii++){
        dx = nodal_positions[i][ipoly][ii][0] - xcm[0];
        dy = nodal_positions[i][ipoly][ii][1] - xcm[1];
        dz = nodal_positions[i][ipoly][ii][2] - xcm[2];
        nodal_velocities[i][ipoly][ii][0] -= omega[1]*dz - omega[2]*dy;
        nodal_velocities[i][ipoly][ii][1] -= omega[2]*dx - omega[0]*dz;
        nodal_velocities[i][ipoly][ii][2] -= omega[0]*dy - omega[1]*dx;
        }
      }
  }
}

/* ---------------------------------------------------------------------- */
void FixCACMomentum::quadrature_init(){

  quadrature_node_count=2;
  memory->create(quadrature_weights,quadrature_node_count,"pairCAC:quadrature_weights");
  memory->create(quadrature_abcissae,quadrature_node_count,"pairCAC:quadrature_abcissae");
  quadrature_weights[0]=1;
  quadrature_weights[1]=1;
  quadrature_abcissae[0]=-0.5773502691896258;
  quadrature_abcissae[1]=0.5773502691896258;
  
}

/* ---------------------------------------------------------------------- */
double FixCACMomentum::shape_function(double s, double t, double w, int flag, int index){
double shape_function;
if(flag==2){
    if(index==1){
    shape_function=(1-s)*(1-t)*(1-w)/8;
    }
    else if(index==2){
    shape_function=(1+s)*(1-t)*(1-w)/8;
    }
    else if(index==3){
    shape_function=(1+s)*(1+t)*(1-w)/8;
    }
    else if(index==4){
    shape_function=(1-s)*(1+t)*(1-w)/8;
    }
    else if(index==5){
    shape_function=(1-s)*(1-t)*(1+w)/8;
    }
    else if(index==6){
    shape_function=(1+s)*(1-t)*(1+w)/8;
    }
    else if(index==7){
    shape_function=(1+s)*(1+t)*(1+w)/8;
    }
    else if(index==8){
    shape_function=(1-s)*(1+t)*(1+w)/8;
    }
}
return shape_function;
}