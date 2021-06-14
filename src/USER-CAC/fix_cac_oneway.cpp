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
   Contributing author: Axel Kohlmeyer (ICTP, Italy)
------------------------------------------------------------------------- */

#include <cstring>
#include "fix_cac_oneway.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "region.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE=-1,X=0,Y=1,Z=2,XYZMASK=3,MINUS=4,PLUS=0};

/* ---------------------------------------------------------------------- */

CACFixOneWay::CACFixOneWay(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  direction = NONE;
  regionidx = 0;
  regionstr = NULL;

  if (narg < 6) error->all(FLERR,"Illegal fix oneway command");

  nevery =   nevery = utils::numeric(FLERR,arg[3],false,lmp);
  if (nevery < 1) error->all(FLERR,"Illegal fix oneway command");

  int len = strlen(arg[4]);
  regionstr = new char[len+1];
  strcpy(regionstr,arg[4]);

  if (strcmp(arg[5], "x") == 0) direction = X|PLUS;
  if (strcmp(arg[5], "X") == 0) direction = X|PLUS;
  if (strcmp(arg[5], "y") == 0) direction = Y|PLUS;
  if (strcmp(arg[5], "Y") == 0) direction = Y|PLUS;
  if (strcmp(arg[5], "z") == 0) direction = Z|PLUS;
  if (strcmp(arg[5], "Z") == 0) direction = Z|PLUS;
  if (strcmp(arg[5],"-x") == 0) direction = X|MINUS;
  if (strcmp(arg[5],"-X") == 0) direction = X|MINUS;
  if (strcmp(arg[5],"-y") == 0) direction = Y|MINUS;
  if (strcmp(arg[5],"-Y") == 0) direction = Y|MINUS;
  if (strcmp(arg[5],"-z") == 0) direction = Z|MINUS;
  if (strcmp(arg[5],"-Z") == 0) direction = Z|MINUS;

  global_freq = nevery;
}

/* ---------------------------------------------------------------------- */

CACFixOneWay::~CACFixOneWay()
{
  if (regionstr) delete[] regionstr;
}

/* ---------------------------------------------------------------------- */

int CACFixOneWay::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void CACFixOneWay::init()
{
  //check if CAC atom style is active
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
  //find region id provided as an argument
  regionidx = domain->find_region(regionstr);
  if (regionidx < 0)
    error->all(FLERR,"Region for fix oneway does not exist");
}

/* ---------------------------------------------------------------------- */

void CACFixOneWay::end_of_step()
{
  Region *region = domain->regions[regionidx];
  region->prematch();

  const int idx = direction & XYZMASK;
  const double * const * const x = atom->x;
  double * const * const v = atom->v;
  const int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;	

  for (int i = 0; i < nlocal; ++i) {
    if ((mask[i] & groupbit) && region->match(x[i][0],x[i][1],x[i][2])) {
      nodes_per_element = nodes_count_list[element_type[i]];
      v[i][idx] = 0;
      for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {	
        for(int k=0; k<nodes_per_element; k++){	
          if (direction & MINUS) {
            if (nodal_velocities[i][poly_counter][k][idx] > 0.0) 
              nodal_velocities[i][poly_counter][k][idx] = -nodal_velocities[i][poly_counter][k][idx];
          } 
          else {
            if (nodal_velocities[i][poly_counter][k][idx] < 0.0) 
              nodal_velocities[i][poly_counter][k][idx] = -nodal_velocities[i][poly_counter][k][idx];
          }
          v[i][idx]+=nodal_velocities[i][poly_counter][k][idx];
        }
      }
      v[i][idx] = v[i][idx] / nodes_per_element / poly_count[i];
    }
  }
}

