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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_rigid.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-7
#define MAXJACOBI 50

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixRigid::FixRigid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int i,ibody;

  time_integrate = 1;
  rigid_flag = 1;
  virial_flag = 1;
  create_attribute = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  extended = dorientflag = qorientflag = 0;
  body = NULL;
  displace = NULL;
  eflags = NULL;
  dorient = NULL;
  qorient = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // parse command-line args
  // set nbody and body[i] for each atom

  if (narg < 4) error->all("Illegal fix rigid command");
  
  // single rigid body
  // nbody = 1
  // all atoms in fix group are part of body

  int iarg;

  if (strcmp(arg[3],"single") == 0) {
    iarg = 4;
    nbody = 1;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = 0;
    }

  // each molecule in fix group is a rigid body
  // maxmol = largest molecule #
  // ncount = # of atoms in each molecule (have to sum across procs)
  // nbody = # of non-zero ncount values
  // use nall as incremented ptr to set body[] values for each atom

  } else if (strcmp(arg[3],"molecule") == 0) {
    iarg = 4;
    if (atom->molecular == 0)
      error->all("Must use a molecular atom style with fix rigid molecule");

    int *mask = atom->mask;
    int *molecule = atom->molecule;
    int nlocal = atom->nlocal;

    int maxmol = -1;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) maxmol = MAX(maxmol,molecule[i]);

    int itmp;
    MPI_Allreduce(&maxmol,&itmp,1,MPI_INT,MPI_MAX,world);
    maxmol = itmp + 1;

    int *ncount = new int[maxmol];
    for (i = 0; i < maxmol; i++) ncount[i] = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) ncount[molecule[i]]++;

    int *nall = new int[maxmol];
    MPI_Allreduce(ncount,nall,maxmol,MPI_INT,MPI_SUM,world);

    nbody = 0;
    for (i = 0; i < maxmol; i++)
      if (nall[i]) nall[i] = nbody++;
      else nall[i] = -1;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = nall[molecule[i]];
    }

    delete [] ncount;
    delete [] nall;

  // each listed group is a rigid body
  // check if all listed groups exist
  // an atom must belong to fix group and listed group to be in rigid body
  // error if atom belongs to more than 1 rigid body

  } else if (strcmp(arg[3],"group") == 0) {
    if (narg < 5) error->all("Illegal fix rigid command");
    nbody = atoi(arg[4]);
    if (nbody <= 0) error->all("Illegal fix rigid command");
    if (narg < 5+nbody) error->all("Illegal fix rigid command");
    iarg = 5 + nbody;

    int *igroups = new int[nbody];
    for (ibody = 0; ibody < nbody; ibody++) {
      igroups[ibody] = group->find(arg[5+ibody]);
      if (igroups[ibody] == -1) 
	error->all("Could not find fix rigid group ID");
    }

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit)
	for (ibody = 0; ibody < nbody; ibody++)
	  if (mask[i] & group->bitmask[igroups[ibody]]) {
	    if (body[i] >= 0) flag = 1;
	    body[i] = ibody;
	  }
    }

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) 
      error->all("One or more atoms belong to multiple rigid bodies");

    delete [] igroups;

  } else error->all("Illegal fix rigid command");

  // error check on nbody

  if (nbody == 0) error->all("No rigid bodies defined");

  // create all nbody-length arrays

  nrigid = (int *) memory->smalloc(nbody*sizeof(int),"rigid:nrigid");
  masstotal = (double *)
    memory->smalloc(nbody*sizeof(double),"rigid:masstotal");
  xcm = memory->create_2d_double_array(nbody,3,"rigid:xcm");
  vcm = memory->create_2d_double_array(nbody,3,"rigid:vcm");
  fcm = memory->create_2d_double_array(nbody,3,"rigid:fcm");
  inertia = memory->create_2d_double_array(nbody,3,"rigid:inertia");
  ex_space = memory->create_2d_double_array(nbody,3,"rigid:ex_space");
  ey_space = memory->create_2d_double_array(nbody,3,"rigid:ey_space");
  ez_space = memory->create_2d_double_array(nbody,3,"rigid:ez_space");
  angmom = memory->create_2d_double_array(nbody,3,"rigid:angmom");
  omega = memory->create_2d_double_array(nbody,3,"rigid:omega");
  torque = memory->create_2d_double_array(nbody,3,"rigid:torque");
  quat = memory->create_2d_double_array(nbody,4,"rigid:quat");
  image = (int *) memory->smalloc(nbody*sizeof(int),"rigid:image");
  fflag = memory->create_2d_double_array(nbody,3,"rigid:fflag");
  tflag = memory->create_2d_double_array(nbody,3,"rigid:tflag");

  sum = memory->create_2d_double_array(nbody,6,"rigid:sum");
  all = memory->create_2d_double_array(nbody,6,"rigid:all");
  remapflag = memory->create_2d_int_array(nbody,4,"rigid:remapflag");

  // initialize force/torque flags to default = 1.0

  vector_flag = 1;
  size_vector = 12*nbody;
  scalar_vector_freq = 1;
  extvector = 0;

  for (i = 0; i < nbody; i++) {
    fflag[i][0] = fflag[i][1] = fflag[i][2] = 1.0;
    tflag[i][0] = tflag[i][1] = tflag[i][2] = 1.0;
  }

  // parse optional args that set fflag and tflag

  while (iarg < narg) {
    if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix rigid command");

      int mlo,mhi;
      force->bounds(arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all("Illegal fix rigid command");
      if (strcmp(arg[iarg+3],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all("Illegal fix rigid command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all("Illegal fix rigid command");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
	fflag[m-1][0] = xflag;
	fflag[m-1][1] = yflag;
	fflag[m-1][2] = zflag;
	count++;
      }
      if (count == 0) error->all("Illegal fix rigid command");

      iarg += 5;
    } else if (strcmp(arg[iarg],"torque") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix rigid command");

      int mlo,mhi;
      force->bounds(arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all("Illegal fix rigid command");
      if (strcmp(arg[iarg+3],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all("Illegal fix rigid command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all("Illegal fix rigid command");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
	tflag[m-1][0] = xflag;
	tflag[m-1][1] = yflag;
	tflag[m-1][2] = zflag;
	count++;
      }
      if (count == 0) error->all("Illegal fix rigid command");

      iarg += 5;
    } else error->all("Illegal fix rigid command");
  }

  // initialize vector output quantities in case accessed before run

  for (i = 0; i < nbody; i++) {
    xcm[i][0] = xcm[i][1] = xcm[i][2] = 0.0;
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    fcm[i][0] = fcm[i][1] = fcm[i][2] = 0.0;
    torque[i][0] = torque[i][1] = torque[i][2] = 0.0;
  }

  // nrigid[n] = # of atoms in Nth rigid body
  // error if one or zero atoms
  
  int *ncount = new int[nbody];
  for (ibody = 0; ibody < nbody; ibody++) ncount[ibody] = 0;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) ncount[body[i]]++;
  
  MPI_Allreduce(ncount,nrigid,nbody,MPI_INT,MPI_SUM,world);
  delete [] ncount;

  for (ibody = 0; ibody < nbody; ibody++)
    if (nrigid[ibody] <= 1) error->all("One or zero atoms in rigid body");

  // set image flags for each rigid body to default values
  // will be reset during init() based on xcm and then by pre_neighbor()
  // set here, so image value will persist from run to run

  for (ibody = 0; ibody < nbody; ibody++)
    image[ibody] = (512 << 20) | (512 << 10) | 512;

  // bitmasks for properties of extended particles

  INERTIA_SPHERE_RADIUS = 1;
  INERTIA_SPHERE_SHAPE = 2;
  INERTIA_ELLIPSOID = 4;
  ORIENT_DIPOLE = 8;
  ORIENT_QUAT = 16;
  OMEGA = 32;
  ANGMOM = 64;
  TORQUE = 128;

  // print statistics

  int nsum = 0;
  for (ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];
  
  if (comm->me == 0) {
    if (screen) fprintf(screen,"%d rigid bodies with %d atoms\n",nbody,nsum);
    if (logfile) fprintf(logfile,"%d rigid bodies with %d atoms\n",nbody,nsum);
  }
}

/* ---------------------------------------------------------------------- */

FixRigid::~FixRigid()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  
  // delete locally stored arrays
  
  memory->sfree(body);
  memory->destroy_2d_double_array(displace);
  memory->sfree(eflags);
  memory->destroy_2d_double_array(dorient);
  memory->destroy_2d_double_array(qorient);
  
  // delete nbody-length arrays

  memory->sfree(nrigid);
  memory->sfree(masstotal);
  memory->destroy_2d_double_array(xcm);
  memory->destroy_2d_double_array(vcm);
  memory->destroy_2d_double_array(fcm);
  memory->destroy_2d_double_array(inertia);
  memory->destroy_2d_double_array(ex_space);
  memory->destroy_2d_double_array(ey_space);
  memory->destroy_2d_double_array(ez_space);
  memory->destroy_2d_double_array(angmom);
  memory->destroy_2d_double_array(omega);
  memory->destroy_2d_double_array(torque);
  memory->destroy_2d_double_array(quat);
  memory->sfree(image);
  memory->destroy_2d_double_array(fflag);
  memory->destroy_2d_double_array(tflag);

  memory->destroy_2d_double_array(sum);
  memory->destroy_2d_double_array(all);
  memory->destroy_2d_int_array(remapflag);
}

/* ---------------------------------------------------------------------- */

int FixRigid::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigid::init()
{
  int i,itype,ibody;

  triclinic = domain->triclinic;

  // warn if more than one rigid fix

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"rigid") == 0) count++;
  if (count > 1 && comm->me == 0) error->warning("More than one fix rigid");

  // error if npt,nph fix comes before rigid fix

  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) break;
    if (strcmp(modify->fix[i]->style,"nph") == 0) break;
  }
  if (i < modify->nfix) {
    for (int j = i; j < modify->nfix; j++)
      if (strcmp(modify->fix[j]->style,"rigid") == 0)
	error->all("Rigid fix must come before NPT/NPH fix");
  }

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    step_respa = ((Respa *) update->integrate)->step;

  // extended = 1 if any particle in a rigid body is finite size

  extended = dorientflag = qorientflag = 0;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **shape = atom->shape;
  double *dipole = atom->dipole;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (atom->radius_flag || atom->avec->shape_type) {
    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      else if (shape && shape[type[i]][0] > 0.0) flag = 1;
    }

    MPI_Allreduce(&flag,&extended,1,MPI_INT,MPI_MAX,world);
  }

  // grow extended arrays and set extended flags for each particle
  // vorientflag = 1 if any particles store dipole orientation
  // qorientflag = 1 if any particles store quat orientation

  if (extended) {
    if (atom->mu_flag) dorientflag = 1;
    if (atom->quat_flag) qorientflag = 1;
    grow_arrays(atom->nmax);

    for (i = 0; i < nlocal; i++) {
      eflags[i] = 0;
      if (body[i] < 0) continue;

      // set INERTIA if radius or shape > 0.0

      if (radius) {
	if (radius[i] > 0.0) eflags[i] |= INERTIA_SPHERE_RADIUS;
      } else if (shape) {
	if (shape[type[i]][0] > 0.0) {
	  if (shape[type[i]][0] == shape[type[i]][1] &&
	      shape[type[i]][0] == shape[type[i]][2])
	    eflags[i] |= INERTIA_SPHERE_SHAPE;
	  else eflags[i] |= INERTIA_ELLIPSOID;
	}
      }

      // other flags only set if particle is finite size
      // set DIPOLE if atom->mu and atom->dipole exist and dipole[itype] > 0.0
      // set QUAT if atom->quat exists (could be ellipsoid or sphere)
      // set TORQUE if atom->torque exists
      // set exactly one of OMEGA and ANGMOM so particle contributes once
      // set OMEGA if either radius or rmass exists
      // set ANGMOM if shape and mass exist
      // set OMEGE if atom->angmom doesn't exist

      if (eflags[i] == 0) continue;

      if (atom->mu_flag && dipole && dipole[type[i]] > 0.0)
	eflags[i] |= ORIENT_DIPOLE;
      if (atom->quat_flag) eflags[i] |= ORIENT_QUAT;
      if (atom->torque_flag) eflags[i] |= TORQUE;
      if ((radius || rmass) && atom->omega_flag) eflags[i] |= OMEGA;
      else if (shape && mass && atom->angmom_flag) eflags[i] |= ANGMOM;
      else if (atom->omega_flag) eflags[i] != OMEGA;
      else error->one("Could not set finite-size particle attribute "
		      "in fix rigid");
    }
  }

  // compute masstotal & center-of-mass of each rigid body
  
  double **x = atom->x;
  int *image = atom->image;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double massone,xunwrap,yunwrap,zunwrap;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    sum[ibody][0] += xunwrap * massone;
    sum[ibody][1] += yunwrap * massone;
    sum[ibody][2] += zunwrap * massone;
    sum[ibody][3] += massone;
  }
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
  
  for (ibody = 0; ibody < nbody; ibody++) {
    masstotal[ibody] = all[ibody][3];
    xcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    xcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    xcm[ibody][2] = all[ibody][2]/masstotal[ibody];
  }
  
  // remap the xcm of each body back into simulation box if needed
  // only really necessary the 1st time a run is performed

  pre_neighbor();

  // compute 6 moments of inertia of each body
  // dx,dy,dz = coords relative to center-of-mass
  
  double dx,dy,dz,rad;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += massone * (dy*dy + dz*dz);
    sum[ibody][1] += massone * (dx*dx + dz*dz);
    sum[ibody][2] += massone * (dx*dx + dy*dy);
    sum[ibody][3] -= massone * dx*dy;
    sum[ibody][4] -= massone * dy*dz;
    sum[ibody][5] -= massone * dx*dz;
  }

  // extended particles contribute extra terms to moments of inertia

  if (extended) {
    double ex[3],ey[3],ez[3],idiag[3];
    double p[3][3],ptrans[3][3],ispace[3][3],itemp[3][3];

    double *radius = atom->radius;
    double **shape = atom->shape;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      itype = type[i];
      if (rmass) massone = rmass[i];
      else massone = mass[itype];

      if (eflags[i] & INERTIA_SPHERE_RADIUS) {
	sum[ibody][0] += 0.4 * massone * radius[i]*radius[i];
	sum[ibody][1] += 0.4 * massone * radius[i]*radius[i];
	sum[ibody][2] += 0.4 * massone * radius[i]*radius[i];
      }
      if (eflags[i] & INERTIA_SPHERE_SHAPE) {
	rad = shape[type[i]][0];
	sum[ibody][0] += 0.4 * massone * rad*rad;
	sum[ibody][1] += 0.4 * massone * rad*rad;
	sum[ibody][2] += 0.4 * massone * rad*rad;
      }
      if (eflags[i] & INERTIA_ELLIPSOID) {
	MathExtra::quat_to_mat(atom->quat[i],p);
	MathExtra::quat_to_mat_trans(atom->quat[i],ptrans);
	idiag[0] = 0.2*massone *
	  (shape[itype][1]*shape[itype][1]+shape[itype][2]*shape[itype][2]);
	idiag[1] = 0.2*massone *
	  (shape[itype][0]*shape[itype][0]+shape[itype][2]*shape[itype][2]);
	idiag[2] = 0.2*massone *
	  (shape[itype][0]*shape[itype][0]+shape[itype][1]*shape[itype][1]);
	MathExtra::diag_times3(idiag,ptrans,itemp);
	MathExtra::times3(p,itemp,ispace);
	sum[ibody][0] += ispace[0][0];
	sum[ibody][1] += ispace[1][1];
	sum[ibody][2] += ispace[2][2];
	sum[ibody][3] += ispace[0][1];
	sum[ibody][4] += ispace[1][2];
	sum[ibody][5] += ispace[0][2];
      }
    }
  }
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
  
  // inertia = 3 eigenvalues = principal moments of inertia
  // ex_space,ey_space,ez_space = 3 eigenvectors = principal axes of rigid body
  
  double **tensor = memory->create_2d_double_array(3,3,"fix_rigid:tensor");
  double **evectors = memory->create_2d_double_array(3,3,"fix_rigid:evectors");

  int ierror;
  double ez0,ez1,ez2;

  for (ibody = 0; ibody < nbody; ibody++) {
    tensor[0][0] = all[ibody][0];
    tensor[1][1] = all[ibody][1];
    tensor[2][2] = all[ibody][2];
    tensor[0][1] = tensor[1][0] = all[ibody][3];
    tensor[1][2] = tensor[2][1] = all[ibody][4];
    tensor[0][2] = tensor[2][0] = all[ibody][5];

    ierror = jacobi(tensor,inertia[ibody],evectors);
    if (ierror) error->all("Insufficient Jacobi rotations for rigid body");

    ex_space[ibody][0] = evectors[0][0];
    ex_space[ibody][1] = evectors[1][0];
    ex_space[ibody][2] = evectors[2][0];
    
    ey_space[ibody][0] = evectors[0][1];
    ey_space[ibody][1] = evectors[1][1];
    ey_space[ibody][2] = evectors[2][1];
    
    ez_space[ibody][0] = evectors[0][2];
    ez_space[ibody][1] = evectors[1][2];
    ez_space[ibody][2] = evectors[2][2];

    // if any principal moment < scaled EPSILON, set to 0.0
  
    double max;
    max = MAX(inertia[ibody][0],inertia[ibody][1]);
    max = MAX(max,inertia[ibody][2]);
  
    if (inertia[ibody][0] < EPSILON*max) inertia[ibody][0] = 0.0;
    if (inertia[ibody][1] < EPSILON*max) inertia[ibody][1] = 0.0;
    if (inertia[ibody][2] < EPSILON*max) inertia[ibody][2] = 0.0;

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd evector if needed
  
    ez0 = ex_space[ibody][1]*ey_space[ibody][2] -
      ex_space[ibody][2]*ey_space[ibody][1];
    ez1 = ex_space[ibody][2]*ey_space[ibody][0] -
      ex_space[ibody][0]*ey_space[ibody][2];
    ez2 = ex_space[ibody][0]*ey_space[ibody][1] -
      ex_space[ibody][1]*ey_space[ibody][0];
  
    if (ez0*ez_space[ibody][0] + ez1*ez_space[ibody][1] + 
	ez2*ez_space[ibody][2] < 0.0) {
      ez_space[ibody][0] = -ez_space[ibody][0];
      ez_space[ibody][1] = -ez_space[ibody][1];
      ez_space[ibody][2] = -ez_space[ibody][2];
    }

    // create initial quaternion
  
    q_from_exyz(ex_space[ibody],ey_space[ibody],ez_space[ibody],quat[ibody]);
  }

  // free temporary memory
  
  memory->destroy_2d_double_array(tensor);
  memory->destroy_2d_double_array(evectors);
  
  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  double qc[4];

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) {
      displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
      continue;
    }

    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;
    
    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }
    
    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];
    
    displace[i][0] = dx*ex_space[ibody][0] + dy*ex_space[ibody][1] +
      dz*ex_space[ibody][2];
    displace[i][1] = dx*ey_space[ibody][0] + dy*ey_space[ibody][1] +
      dz*ey_space[ibody][2];
    displace[i][2] = dx*ez_space[ibody][0] + dy*ez_space[ibody][1] +
      dz*ez_space[ibody][2];
    
    if (extended) {
      double **mu = atom->mu;
      double *dipole = atom->dipole;
      int *type = atom->type;

      if (eflags[i] & ORIENT_DIPOLE) {
	dorient[i][0] = mu[i][0]*ex_space[ibody][0] + 
	  mu[i][1]*ex_space[ibody][1] + mu[i][2]*ex_space[ibody][2];
	dorient[i][1] = mu[i][0]*ey_space[ibody][0] + 
	  mu[i][1]*ey_space[ibody][1] + mu[i][2]*ey_space[ibody][2];
	dorient[i][2] = mu[i][0]*ez_space[ibody][0] + 
	  mu[i][1]*ez_space[ibody][1] + mu[i][2]*ez_space[ibody][2];
	MathExtra::snormalize3(dipole[type[i]],dorient[i],dorient[i]);
      } else if (dorientflag)
	dorient[i][0] = dorient[i][1] = dorient[i][2] = 0.0;

      if (eflags[i] & ORIENT_QUAT) {
	qconjugate(quat[ibody],qc);
	quatquat(qc,atom->quat[i],qorient[i]);
	qnormalize(qorient[i]);
      } else if (qorientflag)
	qorient[i][0] = qorient[i][1] = qorient[i][2] = qorient[i][3] = 0.0;
    }
  }
  
  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles contribute extra terms to moments of inertia
  
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += massone * 
      (displace[i][1]*displace[i][1] + displace[i][2]*displace[i][2]);
    sum[ibody][1] += massone * 
      (displace[i][0]*displace[i][0] + displace[i][2]*displace[i][2]);
    sum[ibody][2] += massone * 
      (displace[i][0]*displace[i][0] + displace[i][1]*displace[i][1]);
    sum[ibody][3] -= massone * displace[i][0]*displace[i][1];
    sum[ibody][4] -= massone * displace[i][1]*displace[i][2];
    sum[ibody][5] -= massone * displace[i][0]*displace[i][2];
  }

  if (extended) {
    double ex[3],ey[3],ez[3],idiag[3];
    double p[3][3],ptrans[3][3],ispace[3][3],itemp[3][3];

    double *radius = atom->radius;
    double **shape = atom->shape;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      itype = type[i];
      if (rmass) massone = rmass[i];
      else massone = mass[itype];

      if (eflags[i] & INERTIA_SPHERE_RADIUS) {
	sum[ibody][0] += 0.4 * massone * radius[i]*radius[i];
	sum[ibody][1] += 0.4 * massone * radius[i]*radius[i];
	sum[ibody][2] += 0.4 * massone * radius[i]*radius[i];
      }
      if (eflags[i] & INERTIA_SPHERE_SHAPE) {
	rad = shape[type[i]][0];
	sum[ibody][0] += 0.4 * massone * rad*rad;
	sum[ibody][1] += 0.4 * massone * rad*rad;
	sum[ibody][2] += 0.4 * massone * rad*rad;
      }
      if (eflags[i] & INERTIA_ELLIPSOID) {
	MathExtra::quat_to_mat(qorient[i],p);
	MathExtra::quat_to_mat_trans(qorient[i],ptrans);
	idiag[0] = 0.2*massone *
	  (shape[itype][1]*shape[itype][1]+shape[itype][2]*shape[itype][2]);
	idiag[1] = 0.2*massone *
	  (shape[itype][0]*shape[itype][0]+shape[itype][2]*shape[itype][2]);
	idiag[2] = 0.2*massone *
	  (shape[itype][0]*shape[itype][0]+shape[itype][1]*shape[itype][1]);
	MathExtra::diag_times3(idiag,ptrans,itemp);
	MathExtra::times3(p,itemp,ispace);
	sum[ibody][0] += ispace[0][0];
	sum[ibody][1] += ispace[1][1];
	sum[ibody][2] += ispace[2][2];
	sum[ibody][3] += ispace[0][1];
	sum[ibody][4] += ispace[1][2];
	sum[ibody][5] += ispace[0][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    if (inertia[ibody][0] == 0.0) {
      if (fabs(all[ibody][0]) > TOLERANCE)
	error->all("Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][0]-inertia[ibody][0])/inertia[ibody][0]) > 
	  TOLERANCE) error->all("Fix rigid: Bad principal moments");
    }
    if (inertia[ibody][1] == 0.0) {
      if (fabs(all[ibody][1]) > TOLERANCE)
	error->all("Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][1]-inertia[ibody][1])/inertia[ibody][1]) > 
	  TOLERANCE) error->all("Fix rigid: Bad principal moments");
    }
    if (inertia[ibody][2] == 0.0) {
      if (fabs(all[ibody][2]) > TOLERANCE)
	error->all("Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][2]-inertia[ibody][2])/inertia[ibody][2]) > 
	  TOLERANCE) error->all("Fix rigid: Bad principal moments");
    }
    if (fabs(all[ibody][3]) > TOLERANCE || 
	fabs(all[ibody][4]) > TOLERANCE ||
	fabs(all[ibody][5]) > TOLERANCE)
      error->all("Fix rigid: Bad principal moments");
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::setup(int vflag)
{
  int i,n,ibody;
  double massone,radone;
  
  // vcm = velocity of center-of-mass of each rigid body
  // fcm = force on center-of-mass of each rigid body

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += v[i][0] * massone;
    sum[ibody][1] += v[i][1] * massone;
    sum[ibody][2] += v[i][2] * massone;
    sum[ibody][3] += f[i][0];
    sum[ibody][4] += f[i][1];
    sum[ibody][5] += f[i][2];
  }
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    vcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    vcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    vcm[ibody][2] = all[ibody][2]/masstotal[ibody];
    fcm[ibody][0] = all[ibody][3];
    fcm[ibody][1] = all[ibody][4];
    fcm[ibody][2] = all[ibody][5];
  }

  // angmom = angular momentum of each rigid body
  // torque = torque on each rigid body
  
  int *image = atom->image;
  double **x = atom->x;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double xunwrap,yunwrap,zunwrap,dx,dy,dz;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];
    
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += dy * massone*v[i][2] - dz * massone*v[i][1];
    sum[ibody][1] += dz * massone*v[i][0] - dx * massone*v[i][2];
    sum[ibody][2] += dx * massone*v[i][1] - dy * massone*v[i][0];
    sum[ibody][3] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][4] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][5] += dx * f[i][1] - dy * f[i][0];
  }

  // extended particles add their rotation/torque to angmom/torque of body

  if (extended) {
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double **torque_one = atom->torque;
    double *radius = atom->radius;
    double **shape = atom->shape;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & OMEGA) {
	if (rmass) massone = rmass[i];
	else massone = mass[type[i]];
	if (radius) radone = radius[i];
	else radone = shape[type[i]][0];
	sum[ibody][0] += 0.4 * massone * radone*radone * omega_one[i][0];
	sum[ibody][1] += 0.4 * massone * radone*radone * omega_one[i][1];
	sum[ibody][2] += 0.4 * massone * radone*radone * omega_one[i][2];
      }

      if (eflags[i] & ANGMOM) {
	sum[ibody][0] += angmom_one[i][0];
	sum[ibody][1] += angmom_one[i][1];
	sum[ibody][2] += angmom_one[i][2];
      }
      if (eflags[i] & TORQUE) {
	sum[ibody][3] += torque_one[i][0];
	sum[ibody][4] += torque_one[i][1];
	sum[ibody][5] += torque_one[i][2];
      }
    }
  }
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    angmom[ibody][0] = all[ibody][0];
    angmom[ibody][1] = all[ibody][1];
    angmom[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];
  }

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set velocities from angmom & omega

  for (ibody = 0; ibody < nbody; ibody++)
    omega_from_angmom(angmom[ibody],ex_space[ibody],ey_space[ibody],
		      ez_space[ibody],inertia[ibody],omega[ibody]);
  set_v();

  // guesstimate virial as 2x the set_v contribution

  if (vflag_global)
    for (n = 0; n < 6; n++) virial[n] *= 2.0;
  if (vflag_atom) {
    for (i = 0; i < nlocal; i++)
      for (n = 0; n < 6; n++)
	vatom[i][n] *= 2.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate(int vflag)
{
  double dtfm;

  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step
  
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
  
    // update xcm by full step
  
    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];
  
    // update angular momentum by 1/2 step
    
    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];

    // update quaternion a full step
    // returns new normalized quat
    // returns ex_space,ey_space,ez_space for new quat
    // returns omega at 1/2 step (depends on angmom and quat)

    richardson(quat[ibody],omega[ibody],angmom[ibody],inertia[ibody],
	       ex_space[ibody],ey_space[ibody],ez_space[ibody]);
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega
  
  set_xv();
}

/* ---------------------------------------------------------------------- */

void FixRigid::richardson(double *q, double *w,
			  double *m, double *moments,
			  double *ex, double *ey, double *ez)
{
  // compute omega at 1/2 step from m at 1/2 step and q at 0
  
  omega_from_angmom(m,ex,ey,ez,moments,w);

  // full update from dq/dt = 1/2 w q

  double wq[4];
  vecquat(w,q,wq);

  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];

  qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];

  qnormalize(qhalf);

  // udpate ex,ey,ez from qhalf
  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  exyz_from_q(qhalf,ex,ey,ez);
  omega_from_angmom(m,ex,ey,ez,moments,w);
  vecquat(w,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];

  qnormalize(qhalf);

  // corrected Richardson update

  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];

  qnormalize(q);
  exyz_from_q(q,ex,ey,ez);
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate()
{
  int i,ibody;
  double dtfm;
  double xy,xz,yz;

  // sum over atoms to get force and torque on rigid body
  
  int *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  int xbox,ybox,zbox;
  double xunwrap,yunwrap,zunwrap,dx,dy,dz;
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];
      
    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];
    
    sum[ibody][3] += dy*f[i][2] - dz*f[i][1];
    sum[ibody][4] += dz*f[i][0] - dx*f[i][2];
    sum[ibody][5] += dx*f[i][1] - dy*f[i][0];
  }

  // extended particles add their torque to torque of body

  if (extended) {
    double **torque_one = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & TORQUE) {
	sum[ibody][3] += torque_one[i][0];
	sum[ibody][4] += torque_one[i][1];
	sum[ibody][5] += torque_one[i][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
  
  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];

    // update vcm by 1/2 step
  
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
  
    // update angular momentum by 1/2 step
    
    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  for (ibody = 0; ibody < nbody; ibody++) 
    omega_from_angmom(angmom[ibody],ex_space[ibody],ey_space[ibody],
		      ez_space[ibody],inertia[ibody],omega[ibody]);

  set_v();
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtq = 0.5 * step_respa[ilevel];

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate_respa(int ilevel)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ----------------------------------------------------------------------
   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to 1st definition of rigid body or due to box flip
   if don't do this, then atoms of a body which drifts far away
     from a triclinic box will be remapped back into box
     with huge displacements when the box tilt changes via set_x() 
------------------------------------------------------------------------- */

void FixRigid::pre_neighbor()
{
  int original,oldimage,newimage;

  for (int ibody = 0; ibody < nbody; ibody++) {
    original = image[ibody];
    domain->remap(xcm[ibody],image[ibody]);
    
    if (original == image[ibody]) remapflag[ibody][3] = 0;
    else {
      oldimage = original & 1023;
      newimage = image[ibody] & 1023;
      remapflag[ibody][0] = newimage - oldimage;
      oldimage = (original >> 10) & 1023;
      newimage = (image[ibody] >> 10) & 1023;
      remapflag[ibody][1] = newimage - oldimage;
      oldimage = original >> 20;
      newimage = image[ibody] >> 20;
      remapflag[ibody][2] = newimage - oldimage;
      remapflag[ibody][3] = 1;
    }
  }

  // adjust image flags of any atom in a rigid body whose xcm was remapped

  int *atomimage = atom->image;
  int nlocal = atom->nlocal;

  int ibody,idim,otherdims;

  for (int i = 0; i < nlocal; i++) {
    if (body[i] == -1) continue;
    if (remapflag[body[i]][3] == 0) continue;
    ibody = body[i];

    if (remapflag[ibody][0]) {
      idim = atomimage[i] & 1023;
      otherdims = atomimage[i] ^ idim;
      idim -= remapflag[ibody][0];
      idim &= 1023;
      atomimage[i] = otherdims | idim;
    }
    if (remapflag[ibody][1]) {
      idim = (atomimage[i] >> 10) & 1023;
      otherdims = atomimage[i] ^ (idim << 10);
      idim -= remapflag[ibody][1];
      idim &= 1023;
      atomimage[i] = otherdims | (idim << 10);
    }
    if (remapflag[ibody][2]) {
      idim = atomimage[i] >> 20;
      otherdims = atomimage[i] ^ (idim << 20);
      idim -= remapflag[ibody][2];
      idim &= 1023;
      atomimage[i] = otherdims | (idim << 20);
    }
  }
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_rigid for atoms in igroup 
------------------------------------------------------------------------- */

int FixRigid::dof(int igroup)
{
  int groupbit = group->bitmask[igroup];

  // nall = # of point particles in each rigid body
  // mall = # of finite-size particles in each rigid body
  // particles must also be in temperature group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *ncount = new int[nbody];
  int *mcount = new int[nbody];
  for (int ibody = 0; ibody < nbody; ibody++)
    ncount[ibody] = mcount[ibody] = 0;

  for (int i = 0; i < nlocal; i++)
    if (body[i] >= 0 && mask[i] & groupbit) {
      if (extended && eflags[i]) mcount[body[i]]++;
      else ncount[body[i]]++;
    }

  int *nall = new int[nbody];
  int *mall = new int[nbody];
  MPI_Allreduce(ncount,nall,nbody,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(mcount,mall,nbody,MPI_INT,MPI_SUM,world);

  // warn if nall+mall != nrigid for any body included in temperature group

  int flag = 0;
  for (int ibody = 0; ibody < nbody; ibody++) {
    if (nall[ibody]+mall[ibody] > 0 && 
	nall[ibody]+mall[ibody] != nrigid[ibody]) flag = 1;
  }
  if (flag && comm->me == 0)
    error->warning("Computing temperature of portions of rigid bodies");

  // remove appropriate DOFs for each rigid body wholly in temperature group
  // N = # of point particles in body
  // M = # of finite-size particles in body
  // 3d body has 3N + 6M dof to start with
  // 2d body has 2N + 3M dof to start with
  // 3d point-particle body with all non-zero I should have 6 dof, remove 3N-6
  // 3d point-particle body (linear) with a 0 I should have 5 dof, remove 3N-5
  // 2d point-particle body should have 3 dof, remove 2N-3
  // 3d body with any finite-size M should have 6 dof, remove (3N+6M) - 6
  // 2d body with any finite-size M should have 3 dof, remove (2N+3M) - 3

  int n = 0;
  if (domain->dimension == 3) {
    for (int ibody = 0; ibody < nbody; ibody++)
      if (nall[ibody]+mall[ibody] == nrigid[ibody]) {
	n += 3*nall[ibody] + 6*mall[ibody] - 6;
	if (inertia[ibody][0] == 0.0 || inertia[ibody][1] == 0.0 || 
	    inertia[ibody][2] == 0.0) n++;
      }
  } else if (domain->dimension == 2) {
    for (int ibody = 0; ibody < nbody; ibody++)
      if (nall[ibody]+mall[ibody] == nrigid[ibody])
	n += 2*nall[ibody] + 3*mall[ibody] - 3;
  }

  delete [] ncount;
  delete [] mcount;
  delete [] nall;
  delete [] mall;

  return n;
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

void FixRigid::deform(int flag)
{
  if (flag == 0) 
    for (int ibody = 0; ibody < nbody; ibody++)
      domain->x2lamda(xcm[ibody],xcm[ibody]);
  else
    for (int ibody = 0; ibody < nbody; ibody++)
      domain->lamda2x(xcm[ibody],xcm[ibody]);
}

/* ----------------------------------------------------------------------
   compute evalues and evectors of 3x3 real symmetric matrix
   based on Jacobi rotations
   adapted from Numerical Recipes jacobi() function
------------------------------------------------------------------------- */

int FixRigid::jacobi(double **matrix, double *evalues, double **evectors)
{
  int i,j,k;
  double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) evectors[i][j] = 0.0;
    evectors[i][i] = 1.0;
  }
  for (i = 0; i < 3; i++) {
    b[i] = evalues[i] = matrix[i][i];
    z[i] = 0.0;
  }
  
  for (int iter = 1; iter <= MAXJACOBI; iter++) {
    sm = 0.0;
    for (i = 0; i < 2; i++)
      for (j = i+1; j < 3; j++)
	sm += fabs(matrix[i][j]);
    if (sm == 0.0) return 0;
    
    if (iter < 4) tresh = 0.2*sm/(3*3);
    else tresh = 0.0;
    
    for (i = 0; i < 2; i++) {
      for (j = i+1; j < 3; j++) {
	g = 100.0*fabs(matrix[i][j]);
	if (iter > 4 && fabs(evalues[i])+g == fabs(evalues[i])
	    && fabs(evalues[j])+g == fabs(evalues[j]))
	  matrix[i][j] = 0.0;
	else if (fabs(matrix[i][j]) > tresh) {
	  h = evalues[j]-evalues[i];
	  if (fabs(h)+g == fabs(h)) t = (matrix[i][j])/h;
	  else {
	    theta = 0.5*h/(matrix[i][j]);
	    t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt(1.0+t*t);
	  s = t*c;
	  tau = s/(1.0+c);
	  h = t*matrix[i][j];
	  z[i] -= h;
	  z[j] += h;
	  evalues[i] -= h;
	  evalues[j] += h;
	  matrix[i][j] = 0.0;
	  for (k = 0; k < i; k++) rotate(matrix,k,i,k,j,s,tau);
	  for (k = i+1; k < j; k++) rotate(matrix,i,k,k,j,s,tau);
	  for (k = j+1; k < 3; k++) rotate(matrix,i,k,j,k,s,tau);
	  for (k = 0; k < 3; k++) rotate(evectors,k,i,k,j,s,tau);
	}
      }
    }
    
    for (i = 0; i < 3; i++) {
      evalues[i] = b[i] += z[i];
      z[i] = 0.0;
    }
  }
  return 1;
}

/* ----------------------------------------------------------------------
   perform a single Jacobi rotation
------------------------------------------------------------------------- */

void FixRigid::rotate(double **matrix, int i, int j, int k, int l,
		      double s, double tau)
{
  double g = matrix[i][j];
  double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}

/* ----------------------------------------------------------------------
   create unit quaternion from space-frame ex,ey,ez
   ex,ey,ez are columns of a rotation matrix
------------------------------------------------------------------------- */

void FixRigid::q_from_exyz(double *ex, double *ey, double *ez, double *q)
{
  // squares of quaternion components
  
  double q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0);
  double q1sq = q0sq - 0.5 * (ey[1] + ez[2]);
  double q2sq = q0sq - 0.5 * (ex[0] + ez[2]);
  double q3sq = q0sq - 0.5 * (ex[0] + ey[1]);
  
  // some component must be greater than 1/4 since they sum to 1
  // compute other components from it
  
  if (q0sq >= 0.25) {
    q[0] = sqrt(q0sq);
    q[1] = (ey[2] - ez[1]) / (4.0*q[0]);
    q[2] = (ez[0] - ex[2]) / (4.0*q[0]);
    q[3] = (ex[1] - ey[0]) / (4.0*q[0]);
  } else if (q1sq >= 0.25) {
    q[1] = sqrt(q1sq);
    q[0] = (ey[2] - ez[1]) / (4.0*q[1]);
    q[2] = (ey[0] + ex[1]) / (4.0*q[1]);
    q[3] = (ex[2] + ez[0]) / (4.0*q[1]);
  } else if (q2sq >= 0.25) {
    q[2] = sqrt(q2sq);
    q[0] = (ez[0] - ex[2]) / (4.0*q[2]);
    q[1] = (ey[0] + ex[1]) / (4.0*q[2]);
    q[3] = (ez[1] + ey[2]) / (4.0*q[2]);
  } else if (q3sq >= 0.25) {
    q[3] = sqrt(q3sq);
    q[0] = (ex[1] - ey[0]) / (4.0*q[3]);
    q[1] = (ez[0] + ex[2]) / (4.0*q[3]);
    q[2] = (ez[1] + ey[2]) / (4.0*q[3]);
  } else
    error->all("Quaternion creation numeric error");

  qnormalize(q);
}

/* ----------------------------------------------------------------------
   compute space-frame ex,ey,ez from current quaternion q
   ex,ey,ez = space-frame coords of 1st,2nd,3rd principal axis
   operation is ex = q' d q = Q d, where d is (1,0,0) = 1st axis in body frame
------------------------------------------------------------------------- */

void FixRigid::exyz_from_q(double *q, double *ex, double *ey, double *ez)
{
  ex[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  ex[1] = 2.0 * (q[1]*q[2] + q[0]*q[3]);
  ex[2] = 2.0 * (q[1]*q[3] - q[0]*q[2]);
  
  ey[0] = 2.0 * (q[1]*q[2] - q[0]*q[3]);
  ey[1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  ey[2] = 2.0 * (q[2]*q[3] + q[0]*q[1]);
  
  ez[0] = 2.0 * (q[1]*q[3] + q[0]*q[2]);
  ez[1] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
  ez[2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

/* ----------------------------------------------------------------------
   vector-quaternion multiply: c = a*b, where a = (0,a)
------------------------------------------------------------------------- */

void FixRigid::vecquat(double *a, double *b, double *c)
{
  c[0] = -(a[0]*b[1] + a[1]*b[2] + a[2]*b[3]);
  c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
  c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
  c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
}

/* ----------------------------------------------------------------------
   quaternion-vector multiply: c = a*b, where b = (0,b)
------------------------------------------------------------------------- */

void FixRigid::quatvec(double *a, double *b, double *c)
{
  c[0] = - (a[1]*b[0] + a[2]*b[1] + a[3]*b[2]);
  c[1] = a[0]*b[0] + a[2]*b[2] - a[3]*b[1];
  c[2] = a[0]*b[1] + a[3]*b[0] - a[1]*b[2];
  c[3] = a[0]*b[2] + a[1]*b[1] - a[2]*b[0];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */

void FixRigid::quatquat(double *a, double *b, double *c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  c[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

void FixRigid::qconjugate(double *q, double *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

void FixRigid::qnormalize(double *q)
{
  double norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum, both in space frame
   only know Idiag so need to do M = Iw in body frame
   ex,ey,ez are column vectors of rotation matrix P
   Mbody = P_transpose Mspace
   wbody = Mbody / Idiag
   wspace = P wbody
   set wbody component to 0.0 if inertia component is 0.0
     otherwise body can spin easily around that axis
------------------------------------------------------------------------- */

void FixRigid::omega_from_angmom(double *m, double *ex, double *ey, double *ez,
				 double *idiag, double *w)
{
  double wbody[3];

  if (idiag[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / idiag[0];
  if (idiag[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / idiag[1];
  if (idiag[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / idiag[2];

  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   compute angular momentum from omega, both in space frame
   only know Idiag so need to do M = Iw in body frame
   ex,ey,ez are column vectors of rotation matrix P
   wbody = P_transpose wspace
   Mbody = Idiag wbody
   Mspace = P Mbody
------------------------------------------------------------------------- */

void FixRigid::angmom_from_omega(double *w, 
				  double *ex, double *ey, double *ez,
				  double *idiag, double *m)
{
  double mbody[3];

  mbody[0] = (w[0]*ex[0] + w[1]*ex[1] + w[2]*ex[2]) * idiag[0];
  mbody[1] = (w[0]*ey[0] + w[1]*ey[1] + w[2]*ey[2]) * idiag[1];
  mbody[2] = (w[0]*ez[0] + w[1]*ez[1] + w[2]*ez[2]) * idiag[2];

  m[0] = mbody[0]*ex[0] + mbody[1]*ey[0] + mbody[2]*ez[0];
  m[1] = mbody[0]*ex[1] + mbody[1]*ey[1] + mbody[2]*ez[1];
  m[2] = mbody[0]*ex[2] + mbody[1]*ey[2] + mbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_xv()
{
  int ibody,itype;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double ione[3],exone[3],eyone[3],ezone[3],vr[6],p[3][3];

  int *image = atom->image;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass; 
  double *mass = atom->mass; 
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set x and v of each atom
  
  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    // save old positions and velocities for virial

    if (evflag) {
      if (triclinic == 0) {
	x0 = x[i][0] + xbox*xprd;
	x1 = x[i][1] + ybox*yprd;
	x2 = x[i][2] + zbox*zprd;
      } else {
	x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
	x1 = x[i][1] + ybox*yprd + zbox*yz;
	x2 = x[i][2] + zbox*zprd;
      }
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    x[i][0] = ex_space[ibody][0]*displace[i][0] +
      ey_space[ibody][0]*displace[i][1] + 
      ez_space[ibody][0]*displace[i][2];
    x[i][1] = ex_space[ibody][1]*displace[i][0] +
      ey_space[ibody][1]*displace[i][1] + 
      ez_space[ibody][1]*displace[i][2];
    x[i][2] = ex_space[ibody][2]*displace[i][0] +
      ey_space[ibody][2]*displace[i][1] + 
      ez_space[ibody][2]*displace[i][2];

    v[i][0] = omega[ibody][1]*x[i][2] - omega[ibody][2]*x[i][1] +
      vcm[ibody][0];
    v[i][1] = omega[ibody][2]*x[i][0] - omega[ibody][0]*x[i][2] +
      vcm[ibody][1];
    v[i][2] = omega[ibody][0]*x[i][1] - omega[ibody][1]*x[i][0] +
      vcm[ibody][2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    if (triclinic == 0) {
      x[i][0] += xcm[ibody][0] - xbox*xprd;
      x[i][1] += xcm[ibody][1] - ybox*yprd;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    } else {
      x[i][0] += xcm[ibody][0] - xbox*xprd - ybox*xy - zbox*xz;
      x[i][1] += xcm[ibody][1] - ybox*yprd - zbox*yz;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    }

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2]; 

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set orientation, omega, angmom of each extended particle
  
  if (extended) {
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double *dipole = atom->dipole;
    double **shape = atom->shape;

    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];
      
      if (eflags[i] & ORIENT_DIPOLE) {
	MathExtra::quat_to_mat(quat[ibody],p);
	MathExtra::times_column3(p,dorient[i],atom->mu[i]);
	MathExtra::snormalize3(dipole[type[i]],atom->mu[i],atom->mu[i]);
      }
      if (eflags[i] & ORIENT_QUAT) {
	quatquat(quat[ibody],qorient[i],atom->quat[i]);
	qnormalize(atom->quat[i]);
      }
      if (eflags[i] & OMEGA) {
	omega_one[i][0] = omega[ibody][0];
	omega_one[i][1] = omega[ibody][1];
	omega_one[i][2] = omega[ibody][2];
      }
      if (eflags[i] & ANGMOM) {
	itype = type[i];
	ione[0] = 0.2*mass[itype] *
	  (shape[itype][1]*shape[itype][1] + shape[itype][2]*shape[itype][2]);
	ione[1] = 0.2*mass[itype] *
	  (shape[itype][0]*shape[itype][0] + shape[itype][2]*shape[itype][2]);
	ione[2] = 0.2*mass[itype] * 
	  (shape[itype][0]*shape[itype][0] + shape[itype][1]*shape[itype][1]);
	exyz_from_q(atom->quat[i],exone,eyone,ezone);
	angmom_from_omega(omega[ibody],exone,eyone,ezone,ione,angmom_one[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_v()
{
  int ibody,itype;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double ione[3],exone[3],eyone[3],ezone[3],vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass; 
  double *mass = atom->mass; 
  int *type = atom->type;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    dx = ex_space[ibody][0]*displace[i][0] +
      ey_space[ibody][0]*displace[i][1] + 
      ez_space[ibody][0]*displace[i][2];
    dy = ex_space[ibody][1]*displace[i][0] +
      ey_space[ibody][1]*displace[i][1] + 
      ez_space[ibody][1]*displace[i][2];
    dz = ex_space[ibody][2]*displace[i][0] +
      ey_space[ibody][2]*displace[i][1] + 
      ez_space[ibody][2]*displace[i][2];

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[ibody][1]*dz - omega[ibody][2]*dy + vcm[ibody][0];
    v[i][1] = omega[ibody][2]*dx - omega[ibody][0]*dz + vcm[ibody][1];
    v[i][2] = omega[ibody][0]*dy - omega[ibody][1]*dx + vcm[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2]; 

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;

      if (triclinic == 0) {
	x0 = x[i][0] + xbox*xprd;
	x1 = x[i][1] + ybox*yprd;
	x2 = x[i][2] + zbox*zprd;
      } else {
	x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
	x1 = x[i][1] + ybox*yprd + zbox*yz;
	x2 = x[i][2] + zbox*zprd;
      }

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set omega, angmom of each extended particle
  
  if (extended) {
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double **shape = atom->shape;

    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & OMEGA) {
	omega_one[i][0] = omega[ibody][0];
	omega_one[i][1] = omega[ibody][1];
	omega_one[i][2] = omega[ibody][2];
      }
      if (eflags[i] & ANGMOM) {
	itype = type[i];
	ione[0] = 0.2*mass[itype] *
	  (shape[itype][1]*shape[itype][1] + shape[itype][2]*shape[itype][2]);
	ione[1] = 0.2*mass[itype] *
	  (shape[itype][0]*shape[itype][0] + shape[itype][2]*shape[itype][2]);
	ione[2] = 0.2*mass[itype] * 
	  (shape[itype][0]*shape[itype][0] + shape[itype][1]*shape[itype][1]);
	exyz_from_q(atom->quat[i],exone,eyone,ezone);
	angmom_from_omega(omega[ibody],exone,eyone,ezone,ione,angmom_one[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixRigid::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  if (extended) {
    bytes += nmax * sizeof(int);
    if (dorientflag) bytes = nmax*3 * sizeof(double);
    if (qorientflag) bytes = nmax*4 * sizeof(double);
  }
  return bytes;
} 

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixRigid::grow_arrays(int nmax)
{
  body = (int *) memory->srealloc(body,nmax*sizeof(int),"rigid:body");
  displace = memory->grow_2d_double_array(displace,nmax,3,"rigid:displace");
  if (extended) {
    eflags = (int *)
      memory->srealloc(eflags,nmax*sizeof(int),"rigid:eflags");
    if (dorientflag)
      dorient = memory->grow_2d_double_array(dorient,nmax,3,"rigid:dorient");
    if (qorientflag)
      qorient = memory->grow_2d_double_array(qorient,nmax,4,"rigid:qorient");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixRigid::copy_arrays(int i, int j)
{
  body[j] = body[i];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];
  if (extended) {
    eflags[j] = eflags[i];
    if (dorientflag) {
      dorient[j][0] = dorient[i][0];
      dorient[j][1] = dorient[i][1];
      dorient[j][2] = dorient[i][2];
    }
    if (qorientflag) {
      qorient[j][0] = qorient[i][0];
      qorient[j][1] = qorient[i][1];
      qorient[j][2] = qorient[i][2];
      qorient[j][3] = qorient[i][3];
    }
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixRigid::set_arrays(int i)
{
  body[i] = -1;
  displace[i][0] = 0.0;
  displace[i][1] = 0.0;
  displace[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixRigid::pack_exchange(int i, double *buf)
{
  buf[0] = body[i];
  buf[1] = displace[i][0];
  buf[2] = displace[i][1];
  buf[3] = displace[i][2];
  if (!extended) return 4;

  int m = 4;
  buf[m++] = eflags[i];
  if (dorientflag) {
    buf[m++] = dorient[i][0];
    buf[m++] = dorient[i][1];
    buf[m++] = dorient[i][2];
  }
  if (qorientflag) {
    buf[m++] = qorient[i][0];
    buf[m++] = qorient[i][1];
    buf[m++] = qorient[i][2];
    buf[m++] = qorient[i][3];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixRigid::unpack_exchange(int nlocal, double *buf)
{
  body[nlocal] = static_cast<int> (buf[0]);
  displace[nlocal][0] = buf[1];
  displace[nlocal][1] = buf[2];
  displace[nlocal][2] = buf[3];
  if (!extended) return 4;

  int m = 4;
  eflags[nlocal] = static_cast<int> (buf[m++]);
  if (dorientflag) {
    dorient[nlocal][0] = buf[m++];
    dorient[nlocal][0] = buf[m++];
    dorient[nlocal][0] = buf[m++];
  }
  if (qorientflag) {
    qorient[nlocal][0] = buf[m++];
    qorient[nlocal][0] = buf[m++];
    qorient[nlocal][0] = buf[m++];
    qorient[nlocal][0] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRigid::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
}

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   12 values per body
   xcm = 1,2,3; vcm = 4,5,6; fcm = 7,8,9; torque = 10,11,12
------------------------------------------------------------------------- */

double FixRigid::compute_vector(int n)
{
  int ibody = n/12;
  int iattribute = (n % 12) / 3;
  int index = n % 3;

  if (iattribute == 0) return xcm[ibody][index];
  else if (iattribute == 1) return vcm[ibody][index];
  else if (iattribute == 2) return fcm[ibody][index];
  else return torque[ibody][index];
}
