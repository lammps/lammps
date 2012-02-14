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
   Contributing author: Craig Tenney (University of Notre Dame)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_restrain.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{DIHEDRAL};

#define TOLERANCE 0.05

/* ---------------------------------------------------------------------- */

FixRestrain::FixRestrain(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int iarg = 6;
  if (narg < iarg) error->all(FLERR,"Illegal fix restrain command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  time_depend = 1;

  // parse standard arguments

  int n_atoms;
  k_start = force->numeric(arg[3]);
  k_stop = force->numeric(arg[4]);
  if (strcmp(arg[5], "dihedral") == 0) {
    rstyle = DIHEDRAL;
    n_atoms = 4;
  } else error->all(FLERR,"Illegal fix restrain command");

  n_bonds = (narg - iarg) / (n_atoms + 1);
  if (narg != iarg + n_bonds * (n_atoms + 1))
    error->all(FLERR,"Illegal fix restrain command");

  // allocate arrays

  memory->create(atom_id,n_bonds,n_atoms,"restrain:atom_id");
  memory->create(target,n_bonds,"restrain:taret");

  // grab atom_ids and restraint target values

  int ibond = 0;
  while (iarg < narg) {
    for (int i = 0; i < n_atoms; i++) {
      atom_id[ibond][i] = force->inumeric(arg[iarg]);
      iarg++;
    }
    target[ibond] = force->numeric(arg[iarg]);
    iarg++;
    ibond++;
  }

  // special treatment for dihedral restraints

  if (rstyle == DIHEDRAL) {
    cos_shift = (double *)
      memory->smalloc(n_bonds * sizeof(double),"restrain:cos_shift");
    sin_shift = (double *)
      memory->smalloc(n_bonds * sizeof(double),"restrain:sin_shift");
    for (ibond = 0; ibond < n_bonds; ibond++) {
      double my_arg = MY_PI * (180.0 + target[ibond]) / 180.0;
      cos_shift[ibond] = cos(my_arg);
      sin_shift[ibond] = sin(my_arg);
    }
  }

  // require atom map to lookup atom IDs

  if (atom->map_style == 0) 
    error->all(FLERR,"Fix restrain requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

FixRestrain::~FixRestrain()
{
  memory->destroy(atom_id);
  memory->destroy(target);

  if (rstyle == DIHEDRAL) {
    memory->sfree(cos_shift);
    memory->sfree(sin_shift);
  }
}

/* ---------------------------------------------------------------------- */

int FixRestrain::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRestrain::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRestrain::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixRestrain::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrain::post_force(int vflag)
{
  energy = 0.0;
  if (rstyle == DIHEDRAL) restrain_dihedral();
}

/* ---------------------------------------------------------------------- */

void FixRestrain::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrain::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   apply dihedral restraints
   adopted from dihedral_charmm
---------------------------------------------------------------------- */

void FixRestrain::restrain_dihedral()
{
  int i1,i2,i3,i4,i,m,n;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;  
  double c,s,p,sx2,sy2,sz2;

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  double k_step = k_start + delta * (k_stop - k_start);

  for (n = 0; n < n_bonds; n++) {
    i1 = atom->map(atom_id[n][0]);
    i2 = atom->map(atom_id[n][1]);
    i3 = atom->map(atom_id[n][2]);
    i4 = atom->map(atom_id[n][3]);

    // insure exactly one processor computes restraint

    if (newton_bond) {
      if (i2 == -1 || i2 >= nlocal) continue;
      if (i1 == -1 || i3 == -1 || i4 == -1) {
	char str[128];
	sprintf(str,
		"Restrain atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		atom_id[n][0],atom_id[n][1],atom_id[n][2],atom_id[n][3],
		comm->me,update->ntimestep);
	error->one(FLERR,str);
      }
    } else {
      if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal) &&
	  (i3 == -1 || i3 >= nlocal) && (i4 == -1 || i3 >= nlocal)) continue;
      if (i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1) {
	char str[128];
	sprintf(str,
		"Restrain atoms %d %d %d %d missing on proc %d at step " 
		BIGINT_FORMAT,
		atom_id[n][0],atom_id[n][1],atom_id[n][2],atom_id[n][3],
		comm->me,update->ntimestep);
	error->one(FLERR,str);
      }
    }

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
    domain->minimum_image(vb1x,vb1y,vb1z);

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    domain->minimum_image(vb2x,vb2y,vb2z);

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;
    domain->minimum_image(vb2xm,vb2ym,vb2zm);

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    domain->minimum_image(vb3x,vb3y,vb3z);
    
    ax = vb1y*vb2zm - vb1z*vb2ym;
    ay = vb1z*vb2xm - vb1x*vb2zm;
    az = vb1x*vb2ym - vb1y*vb2xm;
    bx = vb3y*vb2zm - vb3z*vb2ym;
    by = vb3z*vb2xm - vb3x*vb2zm;
    bz = vb3x*vb2ym - vb3y*vb2xm;

    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
    rg = sqrt(rgsq);
    
    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
	char str[128];
	sprintf(str,"Restrain problem: %d " BIGINT_FORMAT " %d %d %d %d",
		me,update->ntimestep,
		atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	error->warning(FLERR,str);
	fprintf(screen,"  1st atom: %d %g %g %g\n",
		me,x[i1][0],x[i1][1],x[i1][2]);
	fprintf(screen,"  2nd atom: %d %g %g %g\n",
		me,x[i2][0],x[i2][1],x[i2][2]);
	fprintf(screen,"  3rd atom: %d %g %g %g\n",
		me,x[i3][0],x[i3][1],x[i3][2]);
	fprintf(screen,"  4th atom: %d %g %g %g\n",
		me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
         
    m = 1;  //multiplicity
    p = 1.0;
    df1 = 0.0;
    
    for (i = 0; i < m; i++) {
      ddf1 = p*c - df1*s;
      df1 = p*s + df1*c;
      p = ddf1;
    }

    p = p*cos_shift[n] + df1*sin_shift[n];
    df1 = df1*cos_shift[n] - ddf1*sin_shift[n];
    df1 *= -m;
    p += 1.0;
 
    energy = k_step * p; 
       
    fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
    hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
    fga = fg*ra2inv*rginv;
    hgb = hg*rb2inv*rginv;
    gaa = -ra2inv*rg;
    gbb = rb2inv*rg;
    
    dtfx = gaa*ax;
    dtfy = gaa*ay;
    dtfz = gaa*az;
    dtgx = fga*ax - hgb*bx;
    dtgy = fga*ay - hgb*by;
    dtgz = fga*az - hgb*bz;
    dthx = gbb*bx;
    dthy = gbb*by;
    dthz = gbb*bz;
    
    df = -k_step * df1;
 
    sx2 = df*dtgx;
    sy2 = df*dtgy;
    sz2 = df*dtgz;

    f1[0] = df*dtfx;
    f1[1] = df*dtfy;
    f1[2] = df*dtfz;

    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];

    f4[0] = df*dthx;
    f4[1] = df*dthy;
    f4[2] = df*dthz;

    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }
  }
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixRestrain::compute_scalar()
{
  MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
  return energy_all;
}
