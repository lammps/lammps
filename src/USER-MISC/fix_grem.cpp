/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Force scaling fix for gREM.
   Cite: http://dx.doi.org/10.1063/1.3432176
   Cite: http://dx.doi.org/10.1021/acs.jpcb.5b07614
   
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Edyta Malolepsza (Broad Institute)
   Contributing author: David Stelter (Boston University)
   Contributing author: Tom Keyes (Boston University)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "comm.h"
#include "fix_grem.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "input.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixGrem::FixGrem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all(FLERR,"Illegal fix grem command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // tbath - temp of bath, the same as defined in thermostat

  lambda = force->numeric(FLERR,arg[3]);
  eta = force->numeric(FLERR,arg[4]);
  h0 = force->numeric(FLERR,arg[5]);
  tbath = force->numeric(FLERR,arg[6]);
  pressref = force->numeric(FLERR,arg[7]);

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[5];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure/grem";
  newarg[3] = id_temp;
  newarg[4] = id;
  modify->add_compute(5,newarg);
  delete [] newarg;
  
  // create a new compute ke style
  // id = fix-ID + ke

  n = strlen(id) + 8;
  id_ke = new char[n];
  strcpy(id_ke,id);
  strcat(id_ke,"_ke");

  newarg = new char*[3];
  newarg[0] = id_ke;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "ke";
  modify->add_compute(3,newarg);
  delete [] newarg;
  keflag = 1;

  // create a new compute pe style
  // id = fix-ID + pe

  n = strlen(id) + 9;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;
  peflag = 1;

}

/* ---------------------------------------------------------------------- */

FixGrem::~FixGrem()
{
  // delete temperature, pressure and energies if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  if (keflag) modify->delete_compute(id_ke);
  if (peflag) modify->delete_compute(id_pe);
  delete [] id_temp;
  delete [] id_press;
  delete [] id_ke;
  delete [] id_pe;

}

/* ---------------------------------------------------------------------- */

int FixGrem::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGrem::init()
{

  // add check of sign of eta

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix scaledforce does not exist");
  temperature = modify->compute[icompute];

  icompute = modify->find_compute(id_ke);
  if (icompute < 0)
    error->all(FLERR,"Ke ID for fix scaledforce does not exist");
  ke = modify->compute[icompute];

  icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Pe ID for fix scaledforce does not exist");
  pe = modify->compute[icompute];  
}

/* ---------------------------------------------------------------------- */

void FixGrem::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGrem::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGrem::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double tmpvolume = domain->xprd * domain->yprd * domain->zprd;
  double tmppe = pe->compute_scalar();
  // potential energy
  double tmpenthalpy = tmppe+pressref*tmpvolume/(force->nktv2p);

  double teffective = lambda+eta*(tmpenthalpy-h0);
  scale_grem = tbath/teffective;
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] *= scale_grem;
      f[i][1] *= scale_grem;
      f[i][2] *= scale_grem;
    }
  pe->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   extract scale factor
------------------------------------------------------------------------- */

void *FixGrem::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"scale_grem") == 0) {
    return &scale_grem;
  }
  return NULL;
}
