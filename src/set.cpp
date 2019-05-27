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

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "set.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "atom_vec_body.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{ATOM_SELECT,MOL_SELECT,TYPE_SELECT,GROUP_SELECT,REGION_SELECT};

enum{TYPE,TYPE_FRACTION,MOLECULE,X,Y,Z,CHARGE,MASS,SHAPE,LENGTH,TRI,
     DIPOLE,DIPOLE_RANDOM,SPIN,SPIN_RANDOM,QUAT,QUAT_RANDOM,
     THETA,THETA_RANDOM,ANGMOM,OMEGA,
     DIAMETER,DENSITY,VOLUME,IMAGE,BOND,ANGLE,DIHEDRAL,IMPROPER,
     MESO_E,MESO_CV,MESO_RHO,EDPD_TEMP,EDPD_CV,CC,SMD_MASS_DENSITY,
     SMD_CONTACT_RADIUS,DPDTHETA,INAME,DNAME,VX,VY,VZ};

#define BIG INT_MAX

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Set command with no atoms existing");
  if (narg < 3) error->all(FLERR,"Illegal set command");

  // style and ID info

  if (strcmp(arg[0],"atom") == 0) style = ATOM_SELECT;
  else if (strcmp(arg[0],"mol") == 0) style = MOL_SELECT;
  else if (strcmp(arg[0],"type") == 0) style = TYPE_SELECT;
  else if (strcmp(arg[0],"group") == 0) style = GROUP_SELECT;
  else if (strcmp(arg[0],"region") == 0) style = REGION_SELECT;
  else error->all(FLERR,"Illegal set command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id,arg[1]);
  select = NULL;
  selection(atom->nlocal);

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen,"Setting atom values ...\n");

  int allcount,origarg;

  int iarg = 2;
  while (iarg < narg) {
    varflag = varflag1 = varflag2 = varflag3 = varflag4 = 0;
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      set(TYPE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"type/fraction") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      newtype = force->inumeric(FLERR,arg[iarg+1]);
      fraction = force->numeric(FLERR,arg[iarg+2]);
      ivalue = force->inumeric(FLERR,arg[iarg+3]);
      if (newtype <= 0 || newtype > atom->ntypes)
        error->all(FLERR,"Invalid value in set command");
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR,"Invalid value in set command");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      setrandom(TYPE_FRACTION);
      iarg += 4;

    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (!atom->molecule_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MOLECULE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(X);
      iarg += 2;

    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(Y);
      iarg += 2;

    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(Z);
      iarg += 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VX);
      iarg += 2;

    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VY);
      iarg += 2;

    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      set(VZ);
      iarg += 2;

    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->q_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(CHARGE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rmass_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(MASS);
      iarg += 2;

    } else if (strcmp(arg[iarg],"shape") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->ellipsoid_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(SHAPE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"length") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->line_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(LENGTH);
      iarg += 2;

    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->tri_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(TRI);
      iarg += 2;

    } else if (strcmp(arg[iarg],"dipole") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->mu_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(DIPOLE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"dipole/random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      dvalue = force->numeric(FLERR,arg[iarg+2]);
      if (!atom->mu_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      if (dvalue <= 0.0)
        error->all(FLERR,"Invalid dipole length in set command");
      setrandom(DIPOLE_RANDOM);
      iarg += 3;

    } else if (strcmp(arg[iarg],"spin") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else xvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else yvalue = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) varparse(arg[iarg+4],4);
      else zvalue = force->numeric(FLERR,arg[iarg+4]);
      if (!atom->sp_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(SPIN);
      iarg += 5;

    } else if (strcmp(arg[iarg],"spin/random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      dvalue = force->numeric(FLERR,arg[iarg+2]);
      if (!atom->sp_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      if (dvalue <= 0.0)
        error->all(FLERR,"Invalid dipole length in set command");
      setrandom(SPIN_RANDOM);
      iarg += 3;

    } else if (strcmp(arg[iarg],"quat") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) varparse(arg[iarg+4],4);
      else wvalue = force->numeric(FLERR,arg[iarg+4]);
      if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(QUAT);
      iarg += 5;

    } else if (strcmp(arg[iarg],"quat/random") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (!atom->ellipsoid_flag && !atom->tri_flag && !atom->body_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      setrandom(QUAT_RANDOM);
      iarg += 2;

    } else if (strcmp(arg[iarg],"theta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        dvalue = force->numeric(FLERR,arg[iarg+1]);
        dvalue *= MY_PI/180.0;
      }
      if (!atom->line_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(THETA);
      iarg += 2;

    } else if (strcmp(arg[iarg],"theta/random") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (!atom->line_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0)
        error->all(FLERR,"Invalid random number seed in set command");
      set(THETA_RANDOM);
      iarg += 2;

    } else if (strcmp(arg[iarg],"angmom") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->angmom_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(ANGMOM);
      iarg += 4;

    } else if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else xvalue = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
      else yvalue = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
      else zvalue = force->numeric(FLERR,arg[iarg+3]);
      if (!atom->omega_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(OMEGA);
      iarg += 4;

    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->radius_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      set(DIAMETER);
      iarg += 2;

    } else if (strcmp(arg[iarg],"density") == 0 ||
               (strcmp(arg[iarg],"density/disc") == 0)) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rmass_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (dvalue <= 0.0) error->all(FLERR,"Invalid density in set command");
      discflag = 0;
      if (strcmp(arg[iarg],"density/disc") == 0) {
        discflag = 1;
        if (domain->dimension != 2)
          error->all(FLERR,"Density/disc option requires 2d simulation");
      }
      set(DENSITY);
      iarg += 2;

    } else if (strcmp(arg[iarg],"volume") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->vfrac_flag)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (dvalue <= 0.0) error->all(FLERR,"Invalid volume in set command");
      set(VOLUME);
      iarg += 2;

    } else if (strcmp(arg[iarg],"image") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal set command");
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1],"NULL") != 0) {
        ximageflag = 1;
        if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
        else ximage = force->inumeric(FLERR,arg[iarg+1]);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0) {
        yimageflag = 1;
        if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) varparse(arg[iarg+2],2);
        else yimage = force->inumeric(FLERR,arg[iarg+2]);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0) {
        zimageflag = 1;
        if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) varparse(arg[iarg+3],3);
        else zimage = force->inumeric(FLERR,arg[iarg+3]);
      }
      if (ximageflag && ximage && !domain->xperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (yimageflag && yimage && !domain->yperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (zimageflag && zimage && !domain->zperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      set(IMAGE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->bonds_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nbondtypes)
        error->all(FLERR,"Invalid value in set command");
      topology(BOND);
      iarg += 2;

    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->angles_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nangletypes)
        error->all(FLERR,"Invalid value in set command");
      topology(ANGLE);
      iarg += 2;

    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->dihedrals_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->ndihedraltypes)
        error->all(FLERR,"Invalid value in set command");
      topology(DIHEDRAL);
      iarg += 2;

    } else if (strcmp(arg[iarg],"improper") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      ivalue = force->inumeric(FLERR,arg[iarg+1]);
      if (atom->avec->impropers_allow == 0)
        error->all(FLERR,"Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nimpropertypes)
        error->all(FLERR,"Invalid value in set command");
      topology(IMPROPER);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso/e") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->e_flag)
        error->all(FLERR,"Cannot set meso/e for this atom style");
      set(MESO_E);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso/cv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->cv_flag)
            error->all(FLERR,"Cannot set meso/cv for this atom style");
      set(MESO_CV);
      iarg += 2;

    } else if (strcmp(arg[iarg],"meso/rho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      if (!atom->rho_flag)
        error->all(FLERR,"Cannot set meso/rho for this atom style");
      set(MESO_RHO);
      iarg += 2;

    } else if (strcmp(arg[iarg],"edpd/temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        dvalue = force->numeric(FLERR,arg[iarg+1]);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->edpd_flag)
        error->all(FLERR,"Cannot set edpd/temp for this atom style");
      set(EDPD_TEMP);
      iarg += 2;

    } else if (strcmp(arg[iarg],"edpd/cv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        dvalue = force->numeric(FLERR,arg[iarg+1]);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->edpd_flag)
        error->all(FLERR,"Cannot set edpd/cv for this atom style");
      set(EDPD_CV);
      iarg += 2;

    } else if (strcmp(arg[iarg],"cc") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal set command");
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        cc_index = force->inumeric(FLERR,arg[iarg+1]);
        dvalue = force->numeric(FLERR,arg[iarg+2]);
        if (cc_index < 1) error->all(FLERR,"Illegal set command");
      }
      if (!atom->tdpd_flag)
        error->all(FLERR,"Cannot set cc for this atom style");
      set(CC);
      iarg += 3;

    } else if (strcmp(arg[iarg],"smd/mass/density") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
          if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
          else dvalue = force->numeric(FLERR,arg[iarg+1]);
          if (!atom->smd_flag)
            error->all(FLERR,"Cannot set smd/mass/density for this atom style");
          set(SMD_MASS_DENSITY);
          iarg += 2;

    } else if (strcmp(arg[iarg],"smd/contact/radius") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
          if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
          else dvalue = force->numeric(FLERR,arg[iarg+1]);
          if (!atom->smd_flag)
            error->all(FLERR,"Cannot set smd/contact/radius "
                       "for this atom style");
          set(SMD_CONTACT_RADIUS);
          iarg += 2;

    } else if (strcmp(arg[iarg],"dpd/theta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else {
        dvalue = force->numeric(FLERR,arg[iarg+1]);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->dpd_flag)
        error->all(FLERR,"Cannot set dpd/theta for this atom style");
      set(DPDTHETA);
      iarg += 2;

    } else if (strstr(arg[iarg],"i_") == arg[iarg]) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else ivalue = force->inumeric(FLERR,arg[iarg+1]);
      int flag;
      index_custom = atom->find_custom(&arg[iarg][2],flag);
      if (index_custom < 0 || flag != 0)
        error->all(FLERR,"Set command integer vector does not exist");
      set(INAME);
      iarg += 2;

    } else if (strstr(arg[iarg],"d_") == arg[iarg]) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal set command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) varparse(arg[iarg+1],1);
      else dvalue = force->numeric(FLERR,arg[iarg+1]);
      int flag;
      index_custom = atom->find_custom(&arg[iarg][2],flag);
      if (index_custom < 0 || flag != 1)
        error->all(FLERR,"Set command floating point vector does not exist");
      set(DNAME);
      iarg += 2;

    } else error->all(FLERR,"Illegal set command");

    // statistics
    // for CC option, include species index

    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);

    if (comm->me == 0) {

      if (screen) {
        if (strcmp(arg[origarg],"cc") == 0)
          fprintf(screen,"  %d settings made for %s index %s\n",
                  allcount,arg[origarg],arg[origarg+1]);
        else
          fprintf(screen,"  %d settings made for %s\n",
                  allcount,arg[origarg]);
      }
      if (logfile) {
        if (strcmp(arg[origarg],"cc") == 0)
          fprintf(logfile,"  %d settings made for %s index %s\n",
                  allcount,arg[origarg],arg[origarg+1]);
        else
          fprintf(logfile,"  %d settings made for %s\n",
                  allcount,arg[origarg]);
      }
    }
  }

  // free local memory

  delete [] id;
  delete [] select;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, MOLECULE, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set::selection(int n)
{
  delete [] select;
  select = new int[n];
  int nlo,nhi;

  if (style == ATOM_SELECT) {
    if (atom->tag_enable == 0)
      error->all(FLERR,"Cannot use set atom with no atom IDs defined");
    bigint nlobig,nhibig;
    force->boundsbig(FLERR,id,MAXTAGINT,nlobig,nhibig);

    tagint *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == MOL_SELECT) {
    if (atom->molecule_flag == 0)
      error->all(FLERR,"Cannot use set mol with no molecule IDs defined");
    bigint nlobig,nhibig;
    force->boundsbig(FLERR,id,MAXTAGINT,nlobig,nhibig);

    tagint *molecule = atom->molecule;
    for (int i = 0; i < n; i++)
      if (molecule[i] >= nlobig && molecule[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    force->bounds(FLERR,id,atom->ntypes,nlo,nhi);

    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR,"Could not find set group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all(FLERR,"Set region ID does not exist");
    domain->regions[iregion]->prematch();

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
        select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set owned atom properties directly
   either scalar or per-atom values from atom-style variable(s)
------------------------------------------------------------------------- */

void Set::set(int keyword)
{
  // evaluate atom-style variable(s) if necessary

  vec1 = vec2 = vec3 = vec4 = NULL;

  if (varflag) {
    int nlocal = atom->nlocal;
    if (varflag1) {
      memory->create(vec1,nlocal,"set:vec1");
      input->variable->compute_atom(ivar1,0,vec1,1,0);
    }
    if (varflag2) {
      memory->create(vec2,nlocal,"set:vec2");
      input->variable->compute_atom(ivar2,0,vec2,1,0);
    }
    if (varflag3) {
      memory->create(vec3,nlocal,"set:vec3");
      input->variable->compute_atom(ivar3,0,vec3,1,0);
    }
    if (varflag4) {
      memory->create(vec4,nlocal,"set:vec4");
      input->variable->compute_atom(ivar4,0,vec4,1,0);
    }
  }

  // check if properties of atoms in rigid bodies are updated
  // that are cached as per-body data.
  switch (keyword) {
  case X:
  case Y:
  case Z:
  case MOLECULE:
  case MASS:
  case ANGMOM:
  case SHAPE:
  case DIAMETER:
  case DENSITY:
  case QUAT:
  case IMAGE:
    if (modify->check_rigid_list_overlap(select))
      error->warning(FLERR,"Changing a property of atoms in rigid bodies "
                     "that has no effect unless rigid bodies are rebuild");
    break;
  default: // assume no conflict for all other properties
    break;
  }

  // loop over selected atoms

  AtomVecEllipsoid *avec_ellipsoid =
    (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecLine *avec_line = (AtomVecLine *) atom->style_match("line");
  AtomVecTri *avec_tri = (AtomVecTri *) atom->style_match("tri");
  AtomVecBody *avec_body = (AtomVecBody *) atom->style_match("body");

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    // overwrite dvalue, ivalue, xyzw value if variables defined
    // else the input script scalar value remains in place

    if (varflag) {
      if (varflag1) {
        dvalue = xvalue = vec1[i];
        ivalue = static_cast<int> (dvalue);
      }
      if (varflag2) yvalue = vec2[i];
      if (varflag3) zvalue = vec3[i];
      if (varflag4) wvalue = vec4[i];
    }

    // set values in per-atom arrays
    // error check here in case atom-style variables generated bogus value

    if (keyword == TYPE) {
      if (ivalue <= 0 || ivalue > atom->ntypes)
        error->one(FLERR,"Invalid value in set command");
      atom->type[i] = ivalue;
    }
    else if (keyword == MOLECULE) atom->molecule[i] = ivalue;
    else if (keyword == X) atom->x[i][0] = dvalue;
    else if (keyword == Y) atom->x[i][1] = dvalue;
    else if (keyword == Z) atom->x[i][2] = dvalue;
    else if (keyword == VX) atom->v[i][0] = dvalue;
    else if (keyword == VY) atom->v[i][1] = dvalue;
    else if (keyword == VZ) atom->v[i][2] = dvalue;
    else if (keyword == CHARGE) atom->q[i] = dvalue;
    else if (keyword == MASS) {
      if (dvalue <= 0.0) error->one(FLERR,"Invalid mass in set command");
      atom->rmass[i] = dvalue;
    }
    else if (keyword == DIAMETER) {
      if (dvalue < 0.0) error->one(FLERR,"Invalid diameter in set command");
      atom->radius[i] = 0.5 * dvalue;
    }
    else if (keyword == VOLUME) {
      if (dvalue <= 0.0) error->one(FLERR,"Invalid volume in set command");
      atom->vfrac[i] = dvalue;
    }
    else if (keyword == MESO_E) atom->e[i] = dvalue;
    else if (keyword == MESO_CV) atom->cv[i] = dvalue;
    else if (keyword == MESO_RHO) atom->rho[i] = dvalue;

    else if (keyword == EDPD_TEMP) atom->edpd_temp[i] = dvalue;
    else if (keyword == EDPD_CV) atom->edpd_cv[i] = dvalue;
    else if (keyword == CC) atom->cc[i][cc_index-1] = dvalue;

    else if (keyword == SMD_MASS_DENSITY) {
      // set mass from volume and supplied mass density
      atom->rmass[i] = atom->vfrac[i] * dvalue;
    }
    else if (keyword == SMD_CONTACT_RADIUS) atom->contact_radius[i] = dvalue;

    else if (keyword == DPDTHETA) {
      if (dvalue >= 0.0) atom->dpdTheta[i] = dvalue;
      else {
        double onemass;
        if (atom->rmass) onemass = atom->rmass[i];
        else onemass = atom->mass[atom->type[i]];
        double vx = atom->v[i][0];
        double vy = atom->v[i][1];
        double vz = atom->v[i][2];
        double tfactor = force->mvv2e / (domain->dimension * force->boltz);
        atom->dpdTheta[i] = tfactor * onemass * (vx*vx + vy*vy + vz*vz);
      }
    }

    // set shape of ellipsoidal particle

    else if (keyword == SHAPE) {
      if (xvalue < 0.0 || yvalue < 0.0 || zvalue < 0.0)
        error->one(FLERR,"Invalid shape in set command");
      if (xvalue > 0.0 || yvalue > 0.0 || zvalue > 0.0) {
        if (xvalue == 0.0 || yvalue == 0.0 || zvalue == 0.0)
          error->one(FLERR,"Invalid shape in set command");
      }
      avec_ellipsoid->set_shape(i,0.5*xvalue,0.5*yvalue,0.5*zvalue);
    }

    // set length of line particle

    else if (keyword == LENGTH) {
      if (dvalue < 0.0) error->one(FLERR,"Invalid length in set command");
      avec_line->set_length(i,dvalue);
    }

    // set corners of tri particle

    else if (keyword == TRI) {
      if (dvalue < 0.0) error->one(FLERR,"Invalid length in set command");
      avec_tri->set_equilateral(i,dvalue);
    }

    // set rmass via density
    // if radius > 0.0, treat as sphere or disc
    // if shape > 0.0, treat as ellipsoid (or ellipse, when uncomment below)
    // if length > 0.0, treat as line
    // if area > 0.0, treat as tri
    // else set rmass to density directly

    else if (keyword == DENSITY) {
      if (dvalue <= 0.0) error->one(FLERR,"Invalid density in set command");
      if (atom->radius_flag && atom->radius[i] > 0.0)
        if (discflag)
          atom->rmass[i] = MY_PI*atom->radius[i]*atom->radius[i] * dvalue;
        else
          atom->rmass[i] = 4.0*MY_PI/3.0 *
            atom->radius[i]*atom->radius[i]*atom->radius[i] * dvalue;
      else if (atom->ellipsoid_flag && atom->ellipsoid[i] >= 0) {
        double *shape = avec_ellipsoid->bonus[atom->ellipsoid[i]].shape;
        // enable 2d ellipse (versus 3d ellipsoid) when time integration
        //   options (fix nve/asphere, fix nh/asphere) are also implemented
        // if (discflag)
        // atom->rmass[i] = MY_PI*shape[0]*shape[1] * dvalue;
        // else
        atom->rmass[i] = 4.0*MY_PI/3.0 * shape[0]*shape[1]*shape[2] * dvalue;
      } else if (atom->line_flag && atom->line[i] >= 0) {
        double length = avec_line->bonus[atom->line[i]].length;
        atom->rmass[i] = length * dvalue;
      } else if (atom->tri_flag && atom->tri[i] >= 0) {
        double *c1 = avec_tri->bonus[atom->tri[i]].c1;
        double *c2 = avec_tri->bonus[atom->tri[i]].c2;
        double *c3 = avec_tri->bonus[atom->tri[i]].c3;
        double c2mc1[3],c3mc1[3];
        MathExtra::sub3(c2,c1,c2mc1);
        MathExtra::sub3(c3,c1,c3mc1);
        double norm[3];
        MathExtra::cross3(c2mc1,c3mc1,norm);
        double area = 0.5 * MathExtra::len3(norm);
        atom->rmass[i] = area * dvalue;
      } else atom->rmass[i] = dvalue;
    }

    // set dipole moment

    else if (keyword == DIPOLE) {
      double **mu = atom->mu;
      mu[i][0] = xvalue;
      mu[i][1] = yvalue;
      mu[i][2] = zvalue;
      mu[i][3] = sqrt(mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] +
                      mu[i][2]*mu[i][2]);
    }

    // set magnetic moments

    else if (keyword == SPIN) {
      double **sp = atom->sp;
      double inorm = 1.0/sqrt(xvalue*xvalue+yvalue*yvalue+zvalue*zvalue);
      sp[i][0] = inorm*xvalue;
      sp[i][1] = inorm*yvalue;
      sp[i][2] = inorm*zvalue;
      sp[i][3] = dvalue;
    }

    // set quaternion orientation of ellipsoid or tri or body particle
    // set quaternion orientation of ellipsoid or tri or body particle
    // enforce quat rotation vector in z dir for 2d systems

    else if (keyword == QUAT) {
      double *quat = NULL;
      if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
        quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
      else if (avec_tri && atom->tri[i] >= 0)
        quat = avec_tri->bonus[atom->tri[i]].quat;
      else if (avec_body && atom->body[i] >= 0)
        quat = avec_body->bonus[atom->body[i]].quat;
      else
        error->one(FLERR,"Cannot set quaternion for atom that has none");
      if (domain->dimension == 2 && (xvalue != 0.0 || yvalue != 0.0))
        error->one(FLERR,"Cannot set quaternion with xy components "
                   "for 2d system");

      double theta2 = MY_PI2 * wvalue/180.0;
      double sintheta2 = sin(theta2);
      quat[0] = cos(theta2);
      quat[1] = xvalue * sintheta2;
      quat[2] = yvalue * sintheta2;
      quat[3] = zvalue * sintheta2;
      MathExtra::qnormalize(quat);
    }

    // set theta of line particle

    else if (keyword == THETA) {
      if (atom->line[i] < 0)
        error->one(FLERR,"Cannot set theta for atom that is not a line");
      avec_line->bonus[atom->line[i]].theta = dvalue;
    }

    // set angmom or omega of particle

    else if (keyword == ANGMOM) {
      atom->angmom[i][0] = xvalue;
      atom->angmom[i][1] = yvalue;
      atom->angmom[i][2] = zvalue;
    }

    else if (keyword == OMEGA) {
      atom->omega[i][0] = xvalue;
      atom->omega[i][1] = yvalue;
      atom->omega[i][2] = zvalue;
    }

    // reset any or all of 3 image flags

    else if (keyword == IMAGE) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if (varflag1) ximage = static_cast<int>(xvalue);
      if (varflag2) yimage = static_cast<int>(yvalue);
      if (varflag3) zimage = static_cast<int>(zvalue);
      if (ximageflag) xbox = ximage;
      if (yimageflag) ybox = yimage;
      if (zimageflag) zbox = zimage;
      atom->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }

    // set value for custom integer or double vector

    else if (keyword == INAME) {
      atom->ivector[index_custom][i] = ivalue;
    }

    else if (keyword == DNAME) {
      atom->dvector[index_custom][i] = dvalue;
    }

    count++;
  }

  // update bonus data numbers
  if (keyword == SHAPE) {
    bigint nlocal_bonus = avec_ellipsoid->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus,&atom->nellipsoids,1,
                  MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (keyword == LENGTH) {
    bigint nlocal_bonus = avec_line->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus,&atom->nlines,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (keyword == TRI) {
    bigint nlocal_bonus = avec_tri->nlocal_bonus;
    MPI_Allreduce(&nlocal_bonus,&atom->ntris,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // clear up per-atom memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
  memory->destroy(vec4);
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom coordinates
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set::setrandom(int keyword)
{
  int i;

  AtomVecEllipsoid *avec_ellipsoid =
    (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecLine *avec_line = (AtomVecLine *) atom->style_match("line");
  AtomVecTri *avec_tri = (AtomVecTri *) atom->style_match("tri");
  AtomVecBody *avec_body = (AtomVecBody *) atom->style_match("body");

  RanPark *random = new RanPark(lmp,1);
  double **x = atom->x;
  int seed = ivalue;

  // set fraction of atom types to newtype

  if (keyword == TYPE_FRACTION) {
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (select[i]) {
        random->reset(seed,x[i]);
        if (random->uniform() > fraction) continue;
        atom->type[i] = newtype;
        count++;
      }

  // set dipole moments to random orientations in 3d or 2d
  // dipole length is determined by dipole type array

  } else if (keyword == DIPOLE_RANDOM) {
    double **mu = atom->mu;
    int nlocal = atom->nlocal;

    double msq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          mu[i][0] = random->uniform() - 0.5;
          mu[i][1] = random->uniform() - 0.5;
          mu[i][2] = random->uniform() - 0.5;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][2] *= scale;
          mu[i][3] = dvalue;
          count++;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          mu[i][0] = random->uniform() - 0.5;
          mu[i][1] = random->uniform() - 0.5;
          mu[i][2] = 0.0;
          msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
          scale = dvalue/sqrt(msq);
          mu[i][0] *= scale;
          mu[i][1] *= scale;
          mu[i][3] = dvalue;
          count++;
        }
    }


  // set spin moments to random orientations in 3d or 2d
  // spin length is fixed to unity

  } else if (keyword == SPIN_RANDOM) {
    double **sp = atom->sp;
    int nlocal = atom->nlocal;

    double sp_sq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          sp[i][0] = random->uniform() - 0.5;
          sp[i][1] = random->uniform() - 0.5;
          sp[i][2] = random->uniform() - 0.5;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1] + sp[i][2]*sp[i][2];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][2] *= scale;
          sp[i][3] = dvalue;
          count++;
        }

    } else {
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          random->reset(seed,x[i]);
          sp[i][0] = random->uniform() - 0.5;
          sp[i][1] = random->uniform() - 0.5;
          sp[i][2] = 0.0;
          sp_sq = sp[i][0]*sp[i][0] + sp[i][1]*sp[i][1];
          scale = 1.0/sqrt(sp_sq);
          sp[i][0] *= scale;
          sp[i][1] *= scale;
          sp[i][3] = dvalue;
          count++;
        }
    }

  // set quaternions to random orientations in 3d and 2d

  } else if (keyword == QUAT_RANDOM) {
    int nlocal = atom->nlocal;
    double *quat;

    if (domain->dimension == 3) {
      double s,t1,t2,theta1,theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          else if (avec_tri && atom->tri[i] >= 0)
            quat = avec_tri->bonus[atom->tri[i]].quat;
          else if (avec_body && atom->body[i] >= 0)
            quat = avec_body->bonus[atom->body[i]].quat;
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          random->reset(seed,x[i]);
          s = random->uniform();
          t1 = sqrt(1.0-s);
          t2 = sqrt(s);
          theta1 = 2.0*MY_PI*random->uniform();
          theta2 = 2.0*MY_PI*random->uniform();
          quat[0] = cos(theta2)*t2;
          quat[1] = sin(theta1)*t1;
          quat[2] = cos(theta1)*t1;
          quat[3] = sin(theta2)*t2;
          count++;
        }

    } else {
      double theta2;
      for (i = 0; i < nlocal; i++)
        if (select[i]) {
          if (avec_ellipsoid && atom->ellipsoid[i] >= 0)
            quat = avec_ellipsoid->bonus[atom->ellipsoid[i]].quat;
          else if (avec_body && atom->body[i] >= 0)
            quat = avec_body->bonus[atom->body[i]].quat;
          else
            error->one(FLERR,"Cannot set quaternion for atom that has none");

          random->reset(seed,x[i]);
          theta2 = MY_PI*random->uniform();
          quat[0] = cos(theta2);
          quat[1] = 0.0;
          quat[2] = 0.0;
          quat[3] = sin(theta2);
          count++;
        }
    }

  // set theta to random orientation in 2d

  } else if (keyword == THETA_RANDOM) {
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (select[i]) {
        if (atom->line[i] < 0)
          error->one(FLERR,"Cannot set theta for atom that is not a line");
        random->reset(seed,x[i]);
        avec_line->bonus[atom->line[i]].theta = MY_2PI*random->uniform();
        count++;
      }
    }
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

void Set::topology(int keyword)
{
  int m,atom1,atom2,atom3,atom4;

  // error check

  if (atom->molecular == 2)
    error->all(FLERR,"Cannot set bond topology types for atom style template");

  // border swap to acquire ghost atom info
  // enforce PBC before in case atoms are outside box
  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc

  if (comm->me == 0 && screen) fprintf(screen,"  system init for set ...\n");
  lmp->init();

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // select both owned and ghost atoms

  selection(atom->nlocal + atom->nghost);

  // for BOND, each of 2 atoms must be in group

  if (keyword == BOND) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_bond[i]; m++) {
        atom1 = atom->map(atom->bond_atom[i][m]);
        if (atom1 == -1) error->one(FLERR,"Bond atom missing in set command");
        if (select[i] && select[atom1]) {
          atom->bond_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for ANGLE, each of 3 atoms must be in group

  if (keyword == ANGLE) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_angle[i]; m++) {
        atom1 = atom->map(atom->angle_atom1[i][m]);
        atom2 = atom->map(atom->angle_atom2[i][m]);
        atom3 = atom->map(atom->angle_atom3[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Angle atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3]) {
          atom->angle_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for DIHEDRAL, each of 4 atoms must be in group

  if (keyword == DIHEDRAL) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_dihedral[i]; m++) {
        atom1 = atom->map(atom->dihedral_atom1[i][m]);
        atom2 = atom->map(atom->dihedral_atom2[i][m]);
        atom3 = atom->map(atom->dihedral_atom3[i][m]);
        atom4 = atom->map(atom->dihedral_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Dihedral atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->dihedral_type[i][m] = ivalue;
          count++;
        }
      }
  }

  // for IMPROPER, each of 4 atoms must be in group

  if (keyword == IMPROPER) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (m = 0; m < atom->num_improper[i]; m++) {
        atom1 = atom->map(atom->improper_atom1[i][m]);
        atom2 = atom->map(atom->improper_atom2[i][m]);
        atom3 = atom->map(atom->improper_atom3[i][m]);
        atom4 = atom->map(atom->improper_atom4[i][m]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Improper atom missing in set command");
        if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
          atom->improper_type[i][m] = ivalue;
          count++;
        }
      }
  }
}

/* ---------------------------------------------------------------------- */

void Set::varparse(char *name, int m)
{
  varflag = 1;

  name = &name[2];
  int n = strlen(name) + 1;
  char *str = new char[n];
  strcpy(str,name);

  int ivar = input->variable->find(str);
  delete [] str;

  if (ivar < 0)
    error->all(FLERR,"Variable name for set command does not exist");
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable for set command is invalid style");

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  } else if (m == 4) {
    varflag4 = 1; ivar4 = ivar;
  }
}
