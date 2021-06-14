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

#include "displace_nodes.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "irregular.h"
#include "group.h"
#include "math_const.h"
#include "random_park.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "atom_vec_body.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{MOVE,RAMP,RANDOM,ROTATE};

/* ---------------------------------------------------------------------- */

DisplaceNodes::DisplaceNodes(LAMMPS *lmp) : Command(lmp)
{
  mvec = NULL;
}

/* ---------------------------------------------------------------------- */

DisplaceNodes::~DisplaceNodes()
{
  memory->destroy(mvec);
}

/* ---------------------------------------------------------------------- */

void DisplaceNodes::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"displace_nodes command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace_nodes command");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot displace_nodes after "
               "reading restart file with per-atom info");

  if (comm->me == 0 && screen) fprintf(screen,"Displacing nodes ...\n");

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_nodes group ID");
  groupbit = group->bitmask[igroup];

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Attempting to displace atoms in rigid bodies");

  int style = -1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else error->all(FLERR,"Illegal displace_nodes command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);


  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector or specified variable(s)

  if (style == MOVE) {
    move(0,arg[2],xscale);
    move(1,arg[3],yscale);
    move(2,arg[4],zscale);
  }

  // move atoms in ramped fashion

  if (style == RAMP) {

    int d_dim = 0;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = xscale*utils::numeric(FLERR,arg[4],false,lmp);
    } else if (d_dim == 1) {
      d_lo = yscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = yscale*utils::numeric(FLERR,arg[4],false,lmp);
    } else if (d_dim == 2) {
      d_lo = zscale*utils::numeric(FLERR,arg[3],false,lmp);
      d_hi = zscale*utils::numeric(FLERR,arg[4],false,lmp);
    }

    int coord_dim = 0;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = xscale*utils::numeric(FLERR,arg[7],false,lmp);
    } else if (coord_dim == 1) {
      coord_lo = yscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = yscale*utils::numeric(FLERR,arg[7],false,lmp);
    } else if (coord_dim == 2) {
      coord_lo = zscale*utils::numeric(FLERR,arg[6],false,lmp);
      coord_hi = zscale*utils::numeric(FLERR,arg[7],false,lmp);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double ****nodal_positions = atom->nodal_positions;
    int *element_type = atom->element_type;
    int *npoly = atom->poly_count;
    int **node_types = atom->node_types;
    int *nodes_count_list = atom->nodes_per_element_list;
    double fraction,dramp;
    int nodes_per_element;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][d_dim] = 0;
          nodes_per_element = nodes_count_list[element_type[i]];
          for (int ip = 0; ip < npoly[i]; ip++) {	
            for(int in=0; in<nodes_per_element; in++){
              fraction = (nodal_positions[i][ip][in][coord_dim] - coord_lo) / (coord_hi - coord_lo);
              fraction = MAX(fraction,0.0);
              fraction = MIN(fraction,1.0);
              dramp = d_lo + fraction*(d_hi - d_lo);
              nodal_positions[i][ip][in][d_dim] += dramp;
              x[i][d_dim] += nodal_positions[i][ip][in][d_dim];
            }
          }
          x[i][d_dim] = x[i][d_dim] / nodes_per_element / npoly[i];
        
      }
    }
  }


  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0) {
    char str[128];
    sprintf(str,"Lost atoms via displace_nodes: original " BIGINT_FORMAT
            " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   move atoms either by specified numeric displacement or variable evaluation
------------------------------------------------------------------------- */

void DisplaceNodes::move(int idim, char *arg, double scale)
{
  double **x = atom->x;
  double ****nodal_positions = atom->nodal_positions;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *element_type = atom->element_type;
  int *npoly = atom->poly_count;
  int **node_types = atom->node_types;
  int *nodes_count_list = atom->nodes_per_element_list;	
  
  int nodes_per_element;

  if (strstr(arg,"v_") != arg) {
    double delta = scale*utils::numeric(FLERR,arg,false,lmp);
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        x[i][idim] = 0;
        nodes_per_element = nodes_count_list[element_type[i]];
        for (int ip = 0; ip < npoly[i]; ip++) {	
          for(int in=0; in<nodes_per_element; in++){	
            nodal_positions[i][ip][in][idim] += delta;
            x[i][idim] += nodal_positions[i][ip][in][idim];
          }
        }
        x[i][idim] = x[i][idim] / nodes_per_element / npoly[i];
      }
    } 
  } else {
    int ivar = input->variable->find(arg+2);
    if (ivar < 0)
      error->all(FLERR,"Variable name for displace_nodes does not exist");

    if (input->variable->equalstyle(ivar)) {
      double delta = scale * input->variable->compute_equal(ivar);
      for (int i = 0; i < nlocal; i++){
        if (mask[i] & groupbit){
          x[i][idim] = 0;
          nodes_per_element = nodes_count_list[element_type[i]];
          for (int ip = 0; ip < npoly[i]; ip++) {	
            for(int in=0; in<nodes_per_element; in++){
              nodal_positions[i][ip][in][idim] += delta;
              x[i][idim] += nodal_positions[i][ip][in][idim];
            }
          }
          x[i][idim] = x[i][idim] / nodes_per_element / npoly[i];
        }
      } 
    } else if (input->variable->atomstyle(ivar)) {
      if (mvec == NULL) memory->create(mvec,nlocal,"displace_nodes:mvec");
      input->variable->compute_atom(ivar,igroup,mvec,1,0);
      for (int i = 0; i < nlocal; i++){
        if (mask[i] & groupbit){
          x[i][idim] = 0;
          nodes_per_element = nodes_count_list[element_type[i]];
          for (int ip = 0; ip < npoly[i]; ip++) {	
            for(int in=0; in<nodes_per_element; in++){
              nodal_positions[i][ip][in][idim] += scale*mvec[i];
              x[i][idim] += nodal_positions[i][ip][in][idim];
            }
          }
          x[i][idim] = x[i][idim] / nodes_per_element / npoly[i];
        }
      } 
    } else error->all(FLERR,"Variable for displace_nodes is invalid style");
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_nodes input line
------------------------------------------------------------------------- */

void DisplaceNodes::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace_nodes command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace_nodes command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace_nodes command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace_nodes command");
  }
}