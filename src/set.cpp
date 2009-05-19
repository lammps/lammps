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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "set.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "random_park.h"
#include "math_extra.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{ATOM,GROUP,REGION};
enum{TYPE,TYPE_FRACTION,MOLECULE,
       X,Y,Z,VX,VY,VZ,CHARGE,
       DIPOLE,DIPOLE_RANDOM,QUAT,QUAT_RANDOM,
       DIAMETER,DENSITY,VOLUME,IMAGE,
       BOND,ANGLE,DIHEDRAL,IMPROPER};

/* ---------------------------------------------------------------------- */

Set::Set(LAMMPS *lmp) : Pointers(lmp)
{
  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all("Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all("Set command with no atoms existing");
  if (narg < 3) error->all("Illegal set command");

  // style and ID

  if (strcmp(arg[0],"atom") == 0) style = ATOM;
  else if (strcmp(arg[0],"group") == 0) style = GROUP;
  else if (strcmp(arg[0],"region") == 0) style = REGION;
  else error->all("Illegal set command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id,arg[1]);
  select = NULL;

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen,"Setting atom values ...\n");

  int allcount,origarg;

  int iarg = 2;
  while (iarg < narg) {
    count = 0;
    origarg = iarg;
    
    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (ivalue <= 0 || ivalue > atom->ntypes)
	error->all("Invalid value in set command");
      set(TYPE);
      iarg += 2;
    } else if (strcmp(arg[iarg],"type/fraction") == 0) {
      if (iarg+4 > narg) error->all("Illegal set command");
      newtype = atoi(arg[iarg+1]);
      fraction = atof(arg[iarg+2]);
      ivalue = atoi(arg[iarg+3]);
      if (fraction < 0.0 || fraction > 1.0) 
	error->all("Invalid value in set command");
      if (ivalue <= 0) error->all("Invalid random number seed in set command");
      setrandom(TYPE_FRACTION);
      iarg += 4;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (!atom->molecule_flag)
	error->all("Cannot set this attribute for this atom style");
      set(MOLECULE);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(X);
      iarg += 2;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(Y);
      iarg += 2;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(Z);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VX);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VZ);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->q_flag)
	error->all("Cannot set this attribute for this atom style");
      set(CHARGE);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dipole") == 0) {
      if (iarg+4 > narg) error->all("Illegal set command");
      xvalue = atof(arg[iarg+1]);
      yvalue = atof(arg[iarg+2]);
      zvalue = atof(arg[iarg+3]);
      if (!atom->mu_flag || atom->dipole == NULL)
	error->all("Cannot set this attribute for this atom style");
      set(DIPOLE);
      iarg += 4;
    } else if (strcmp(arg[iarg],"dipole/random") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (!atom->mu_flag || atom->shape == NULL)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0) error->all("Invalid random number seed in set command");
      setrandom(DIPOLE_RANDOM);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quat") == 0) {
      if (iarg+5 > narg) error->all("Illegal set command");
      xvalue = atof(arg[iarg+1]);
      yvalue = atof(arg[iarg+2]);
      zvalue = atof(arg[iarg+3]);
      wvalue = atof(arg[iarg+4]);
      if (!atom->quat_flag)
	error->all("Cannot set this attribute for this atom style");
      set(QUAT);
      iarg += 5;
    } else if (strcmp(arg[iarg],"quat/random") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (!atom->quat_flag)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0) error->all("Invalid random number seed in set command");
      setrandom(QUAT_RANDOM);
      iarg += 2;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->radius_flag)
        error->all("Cannot set this attribute for this atom style");
      set(DIAMETER);
      iarg += 2;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->density_flag)
        error->all("Cannot set this attribute for this atom style");
      set(DENSITY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"volume") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->vfrac_flag)
        error->all("Cannot set this attribute for this atom style");
      set(VOLUME);
      iarg += 2;
    } else if (strcmp(arg[iarg],"image") == 0) {
      if (iarg+4 > narg) error->all("Illegal set command");
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1],"NULL") != 0 && domain->xperiodic) {
	ximageflag = 1;
	ximage = atoi(arg[iarg+1]);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0 && domain->yperiodic) {
	yimageflag = 1;
	yimage = atoi(arg[iarg+2]);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0 && domain->zperiodic) {
	zimageflag = 1;
	zimage = atoi(arg[iarg+3]);
      }
      set(IMAGE);
      iarg += 4;
    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (atom->avec->bonds_allow == 0)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nbondtypes)
	error->all("Invalid value in set command");
      topology(BOND);
      iarg += 2;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (atom->avec->angles_allow == 0)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nangletypes)
	error->all("Invalid value in set command");
      topology(ANGLE);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (atom->avec->dihedrals_allow == 0)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->ndihedraltypes)
	error->all("Invalid value in set command");
      topology(DIHEDRAL);
      iarg += 2;
    } else if (strcmp(arg[iarg],"improper") == 0) {
      if (iarg+2 > narg) error->all("Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (atom->avec->impropers_allow == 0)
	error->all("Cannot set this attribute for this atom style");
      if (ivalue <= 0 || ivalue > atom->nimpropertypes)
	error->all("Invalid value in set command");
      topology(IMPROPER);
      iarg += 2;
    } else error->all("Illegal set command");    

    // statistics

    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
    
    if (comm->me == 0) {
      if (screen) fprintf(screen,"  %d settings made for %s\n",
			  allcount,arg[origarg]);
      if (logfile) fprintf(logfile,"  %d settings made for %s\n",
			   allcount,arg[origarg]);
    }
  }

  // free local memory

  delete [] id;
  delete [] select;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set::selection(int n)
{
  delete [] select;
  select = new int[n];

  if (style == ATOM) {
    if (atom->tag_enable == 0)
      error->all("Cannot use set atom with no atom IDs defined");
    int idatom = atoi(id);

    int *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (idatom == tag[i]) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP) {
    int igroup = group->find(id);
    if (igroup == -1) error->all("Could not find set group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all("Set region ID does not exist");

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set an owned atom property directly
------------------------------------------------------------------------- */

void Set::set(int keyword)
{
  if (keyword == DIPOLE) atom->check_dipole();

  selection(atom->nlocal);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    if (select[i]) {
      if (keyword == TYPE) atom->type[i] = ivalue;
      else if (keyword == MOLECULE) atom->molecule[i] = ivalue;
      else if (keyword == X) atom->x[i][0] = dvalue;
      else if (keyword == Y) atom->x[i][1] = dvalue;
      else if (keyword == Z) atom->x[i][2] = dvalue;
      else if (keyword == VX) atom->v[i][0] = dvalue;
      else if (keyword == VY) atom->v[i][1] = dvalue;
      else if (keyword == VZ) atom->v[i][2] = dvalue;
      else if (keyword == CHARGE) atom->q[i] = dvalue;

      // set radius from diameter
      // set rmass if both rmass and density are defined

      else if (keyword == DIAMETER) {
	atom->radius[i] = 0.5 * dvalue;
	if (atom->rmass_flag && atom->density_flag) 
	  atom->rmass[i] = 4.0*PI/3.0 * 
	    atom->radius[i]*atom->radius[i]*atom->radius[i] * atom->density[i];

      // set density
      // set rmass (granular) if both rmass and radius are defined
      // set rmass (peri) if both rmass and vfrac are defined
	
      } else if (keyword == DENSITY) {
	atom->density[i] = dvalue;
	if (atom->rmass_flag && atom->radius_flag) 
	  atom->rmass[i] = 4.0*PI/3.0 * 
	    atom->radius[i]*atom->radius[i]*atom->radius[i] * 
	    atom->density[i];
	else if (atom->rmass_flag && atom->vfrac_flag)
	  atom->rmass[i] = dvalue;

      } else if (keyword == VOLUME) atom->vfrac[i] = dvalue;

      // reset any or all of 3 image flags

      else if (keyword == IMAGE) {
	int xbox = (atom->image[i] & 1023) - 512;
	int ybox = (atom->image[i] >> 10 & 1023) - 512;
	int zbox = (atom->image[i] >> 20) - 512;
	if (ximageflag) xbox = ximage;
	if (yimageflag) ybox = yimage;
	if (zimageflag) zbox = zimage;
	atom->image[i] = ((zbox + 512 & 1023) << 20) |
	  ((ybox + 512 & 1023) << 10) |	(xbox + 512 & 1023);

      } else if (keyword == DIPOLE) {
	if (atom->dipole[atom->type[i]] > 0.0) {
	  double **mu = atom->mu;
	  mu[i][0] = xvalue;
	  mu[i][1] = yvalue;
	  mu[i][2] = zvalue;
	  double lensq = 
	    mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
	  double scale = atom->dipole[atom->type[i]]/sqrt(lensq);
	  mu[i][0] *= scale;
	  mu[i][1] *= scale;
	  mu[i][2] *= scale;
	}

      } else if (keyword == QUAT) {
	double PI = 4.0*atan(1.0);
	double theta2 = 0.5 * PI * wvalue/180.0;
	double sintheta2 = sin(theta2);
	double **quat = atom->quat;
	quat[i][0] = cos(theta2);
	quat[i][1] = xvalue * sintheta2;
	quat[i][2] = yvalue * sintheta2;
	quat[i][3] = zvalue * sintheta2;
	MathExtra::normalize4(quat[i]);
      }
      count++;
    }
}

/* ----------------------------------------------------------------------
   set an owned atom property randomly
   set seed based on atom tag
   make atom result independent of what proc owns it
------------------------------------------------------------------------- */

void Set::setrandom(int keyword)
{
  int i;

  if (keyword == DIPOLE_RANDOM) atom->check_dipole();

  selection(atom->nlocal);
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
    int *type = atom->type;
    double *dipole = atom->dipole;
    double **mu = atom->mu;
    int nlocal = atom->nlocal;

    double msq,scale;

    if (domain->dimension == 3) {
      for (i = 0; i < nlocal; i++)
	if (select[i]) {
	  if (dipole[type[i]] > 0.0) {
	    random->reset(seed,x[i]);
	    mu[i][0] = random->uniform() - 0.5;
	    mu[i][1] = random->uniform() - 0.5;
	    mu[i][2] = random->uniform() - 0.5;
	    msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
	    scale = dipole[type[i]]/sqrt(msq);
	    mu[i][0] *= scale;
	    mu[i][1] *= scale;
	    mu[i][2] *= scale;
	    count++;
	  }
	}

    } else {
      for (i = 0; i < nlocal; i++)
	if (select[i]) {
	  if (dipole[type[i]] > 0.0) {
	    random->reset(seed,x[i]);
	    mu[i][0] = random->uniform() - 0.5;
	    mu[i][1] = random->uniform() - 0.5;
	    mu[i][2] = 0.0;
	    msq = mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1];
	    scale = dipole[type[i]]/sqrt(msq);
	    mu[i][0] *= scale;
	    mu[i][1] *= scale;
	    count++;
	  }
	}
    }

  // set quaternions to random orientations in 3d or 2d
  // no need to normalize quats since creations algorithms already do

  } else if (keyword == QUAT_RANDOM) {
    double **quat = atom->quat;
    int nlocal = atom->nlocal;

    if (domain->dimension == 3) {
      double s,t1,t2,theta1,theta2;
      double PI = 4.0*atan(1.0);
      for (i = 0; i < nlocal; i++)
	if (select[i]) {
	  random->reset(seed,x[i]);
	  s = random->uniform();
	  t1 = sqrt(1.0-s);
	  t2 = sqrt(s);
	  theta1 = 2.0*PI*random->uniform();
	  theta2 = 2.0*PI*random->uniform();
	  quat[i][0] = cos(theta2)*t2;
	  quat[i][1] = sin(theta1)*t1;
	  quat[i][2] = cos(theta1)*t1;
	  quat[i][3] = sin(theta2)*t2;
	  count++;
	}

    } else {
      double theta2;
      double PI = 4.0*atan(1.0);
      for (i = 0; i < nlocal; i++)
	if (select[i]) {
	  random->reset(seed,x[i]);
	  theta2 = PI*random->uniform();
	  quat[i][0] = cos(theta2);
	  quat[i][1] = 0.0;
	  quat[i][2] = 0.0;
	  quat[i][3] = sin(theta2);
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
	if (atom1 == -1) error->one("Bond atom missing in set command");
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
	  error->one("Angle atom missing in set command");
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
	  error->one("Dihedral atom missing in set command");
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
	  error->one("Improper atom missing in set command");
	if (select[atom1] && select[atom2] && select[atom3] && select[atom4]) {
	  atom->improper_type[i][m] = ivalue;
	  count++;
	}
      }
  }
}
