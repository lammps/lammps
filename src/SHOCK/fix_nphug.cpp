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

/*

  This fix applies the NPHug Hugoniostat method of Ravelo et al.

(Ravelo, Holian, Germann, and Lomdahl, PRB 70 014103 (2004))

It uses the Nose-Hoover thermostat and barostat (fix_nh.html).
The Nose-Hoover barostat is used to compress the system 
to a specified final stress state. This is done either
hydrostatically (using keyword iso, aniso, or tri) or uniaxially
(using keywords x, y, or z).  In the hydrostatic case,
the cell dimensions change dynamically so that the average axial stress
in all three directions converges towards the specified target value. 
In the uniaxial case, the chosen cell dimension changes dynamically 
so that the average
axial stress in that direction converges towards the target value. The
other two cell dimensions are kept fixed (zero lateral strain).

This leads to the following restrictions on the keywords:

- The specified initial and
final target pressures must be the same.

- Only one of the following keywords may be used:
iso, aniso, tri, x, y, z, 

- The keywords xy, xz, yz may not be used.

- The only admissible value for the couple keyword is xyz,
which has the same effect as keyword iso

- The drag parameter is proportional to the beta_H and beta_p
damping coefficients in the Ravelo paper.

- The temp keyword serves only to set the value of tdamp. The initial
and final target temperatures are ignored. 

- The values of tdamp and pdamp are inversely proportional to the
coupling rate nu_H and nu_p in the Ravelo paper

- All other keywords function in the same way. 

*/





#include "string.h"
#include "fix_nphug.h"
#include "modify.h"
#include "error.h"
#include "update.h"
#include "compute.h"
#include "force.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixNPHug::FixNPHug(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  // Hard-code to use initial state of system

  v0_set = p0_set = e0_set = 0;

  // Hard-code to use z direction

  direction = 2;
  if (p_flag[0] == 1 || p_flag[1] == 1 ||
      p_flag[3] == 1 || p_flag[4] == 1 || p_flag[5] == 1)
    error->all("Only pressure control in z direction to be used with fix nphug");
  if (p_flag[2] == 0)
    error->all("Pressure control in z direction must be used with fix nphug");

  if (!tstat_flag)
    error->all("Temperature control must be used with fix nphug");
  if (!pstat_flag)
    error->all("Pressure control must be used with fix nphug");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  // and thus its KE/temperature contribution should use group all

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
  
  newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  modify->add_compute(4,newarg);
  delete [] newarg;
  pflag = 1;

  // create a new compute potential energy compute

  n = strlen(id) + 3;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char*)"all";
  newarg[2] = (char*)"pe";
  modify->add_compute(3,newarg);
  delete [] newarg;
  peflag = 1;
}

/* ---------------------------------------------------------------------- */

FixNPHug::~FixNPHug()
{

  // temp and press computes handled by base class
  // delete pe compute

  if (peflag) modify->delete_compute(id_pe);
  delete [] id_pe;

}

/* ---------------------------------------------------------------------- */

void FixNPHug::init()
{
  // Call base class init()

  FixNH::init();

  // set pe ptr

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0) 
    error->all("Potential energy ID for fix nvt/nph/npt does not exist");
  pe = modify->compute[icompute];
}


/* ----------------------------------------------------------------------
   compute initial state before integrator starts 
------------------------------------------------------------------------- */

void FixNPHug::setup(int vflag)
{
  FixNH::setup(vflag);

  if ( v0_set == 0 ) {
    v0 = compute_vol();
    v0_set = 1;
  } 

  if ( p0_set == 0 ) {
    p0 = p_current[direction];
    p0_set = 1;
  }

  if ( e0_set == 0 ) {
    e0 = compute_etotal();
    e0_set = 1;
  }

}

/* ----------------------------------------------------------------------
   compute target temperature and kinetic energy
-----------------------------------------------------------------------*/

void FixNPHug::compute_temp_target()
{
  t_target = t_current + compute_hugoniot();
  ke_target = tdof * boltz * t_target;
  if (ke_target < 0.0) ke_target = 0.0;
  
  // If t_target is very small, need to choose 
  // more reasonable value for use by barostat and 
  // thermostat masses. ke_target is left as is.

  if (t_target <= 1.0e-6) {
    if (strcmp(update->unit_style,"lj") == 0) t0 = 1.0;
    else t0 = 300.0;
  }
  pressure->addstep(update->ntimestep+1);
  pe->addstep(update->ntimestep+1);
}


/* ---------------------------------------------------------------------- */

double FixNPHug::compute_etotal()
{
  double epot,ekin,etot;
  epot = pe->compute_scalar();
  if (thermo_energy) epot -= compute_scalar();
  ekin = temperature->compute_scalar();
  ekin *= 0.5 * temperature->dof * force->boltz;
  etot = epot+ekin;
  return etot;
}

/* ---------------------------------------------------------------------- */

double FixNPHug::compute_vol()
{
  if (domain->dimension == 3)
    return domain->xprd * domain->yprd * domain->zprd;
  else
    return domain->xprd * domain->yprd;
}

/* ----------------------------------------------------------------------
   Computes the deviation of the current point 
   from the Hugoniot in energy units.
------------------------------------------------------------------------- */

double FixNPHug::compute_hugoniot()
{
  double v, e, p;
  double dhugo;
  
  e = compute_etotal();
  
  temperature->compute_vector();
  pressure->compute_vector();
  p = pressure->vector[direction];
  
  v = compute_vol();
  
  dhugo = (0.5 * (p + p0 ) * ( v0 - v)) / 
    force->nktv2p + e0 - e;

  return dhugo;
}
