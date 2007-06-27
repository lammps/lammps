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
#include "string.h"
#include "fix_rigid.h"
#include "atom.h"
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

  rigid_flag = 1;
  virial_flag = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class
  
  body = NULL;
  displace = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  
  // parse command-line args
  // set nbody and body[i] for each atom

  if (narg < 4) error->all("Illegal fix rigid command");
  
  // single rigid body
  // nbody = 1
  // all atoms in fix group are part of body

  if (strcmp(arg[3],"single") == 0) {
    if (narg != 4) error->all("Illegal fix rigid command");

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
    if (narg != 4) error->all("Illegal fix rigid command");
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
    nbody = narg-4;
    if (nbody <= 0) error->all("Illegal fix rigid command");

    int *igroups = new int[nbody];
    for (ibody = 0; ibody < nbody; ibody++) {
      igroups[ibody] = group->find(arg[ibody+4]);
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

  sum = memory->create_2d_double_array(nbody,6,"rigid:sum");
  all = memory->create_2d_double_array(nbody,6,"rigid:all");
  
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

  // print statistics

  int nsum = 0;
  for (ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];
  
  if (comm->me == 0) {
    if (screen) fprintf(screen,"%d rigid bodies with %d atoms\n",nbody,nsum);
    if (logfile) fprintf(logfile,"%d rigid bodies with %d atoms\n",nbody,nsum);
  }

  // zero fix_rigid virial in case pressure uses it before 1st fix_rigid call

  for (int n = 0; n < 6; n++) virial[n] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixRigid::~FixRigid()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  
  // delete locally stored arrays
  
  memory->sfree(body);
  memory->destroy_2d_double_array(displace);

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

  memory->destroy_2d_double_array(sum);
  memory->destroy_2d_double_array(all);
}

/* ---------------------------------------------------------------------- */

int FixRigid::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigid::init()
{
  int i,ibody;

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

  // compute rigid contribution to virial every step if fix NPT,NPH exists

  pressure_flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) pressure_flag = 1;
    if (strcmp(modify->fix[i]->style,"nph") == 0) pressure_flag = 1;
  }

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    step_respa = ((Respa *) update->integrate)->step;

  // compute masstotal & center-of-mass of each rigid body
  
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double massone;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;
    massone = mass[type[i]];
    
    sum[ibody][0] += (x[i][0] + xbox*xprd) * massone;
    sum[ibody][1] += (x[i][1] + ybox*yprd) * massone;
    sum[ibody][2] += (x[i][2] + zbox*zprd) * massone;
    sum[ibody][3] += massone;
  }
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
  
  for (ibody = 0; ibody < nbody; ibody++) {
    masstotal[ibody] = all[ibody][3];
    xcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    xcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    xcm[ibody][2] = all[ibody][2]/masstotal[ibody];
  }
  
  // compute 6 moments of inertia of each body
  // dx,dy,dz = coords relative to center-of-mass
  
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  double dx,dy,dz;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;
    dx = x[i][0] + xbox*xprd - xcm[ibody][0];
    dy = x[i][1] + ybox*yprd - xcm[ibody][1];
    dz = x[i][2] + zbox*zprd - xcm[ibody][2];
    massone = mass[type[i]];
    
    sum[ibody][0] += massone * (dy*dy + dz*dz);
    sum[ibody][1] += massone * (dx*dx + dz*dz);
    sum[ibody][2] += massone * (dx*dx + dy*dy);
    sum[ibody][3] -= massone * dx*dy;
    sum[ibody][4] -= massone * dy*dz;
    sum[ibody][5] -= massone * dx*dz;
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
  
    qcreate(evectors,quat[ibody]);
  }

  // free temporary memory
  
  memory->destroy_2d_double_array(tensor);
  memory->destroy_2d_double_array(evectors);
  
  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
    else {
      ibody = body[i];

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = x[i][0] + xbox*xprd - xcm[ibody][0];
      dy = x[i][1] + ybox*yprd - xcm[ibody][1];
      dz = x[i][2] + zbox*zprd - xcm[ibody][2];
      
      displace[i][0] = dx*ex_space[ibody][0] + dy*ex_space[ibody][1] +
	dz*ex_space[ibody][2];
      displace[i][1] = dx*ey_space[ibody][0] + dy*ey_space[ibody][1] +
	dz*ey_space[ibody][2];
      displace[i][2] = dx*ez_space[ibody][0] + dy*ez_space[ibody][1] +
	dz*ez_space[ibody][2];
    }
  }
  
  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    massone = mass[type[i]];

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
  
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    if (fabs(all[ibody][0]-inertia[ibody][0]) > TOLERANCE || 
	fabs(all[ibody][1]-inertia[ibody][1]) > TOLERANCE ||
	fabs(all[ibody][2]-inertia[ibody][2]) > TOLERANCE)
      error->all("Bad principal moments");
    if (fabs(all[ibody][3]) > TOLERANCE || 
	fabs(all[ibody][4]) > TOLERANCE ||
	fabs(all[ibody][5]) > TOLERANCE)
      error->all("Bad principal moments");
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::setup()
{
  int i,ibody;
  
  // vcm = velocity of center-of-mass of each rigid body
  // fcm = force on center-of-mass of each rigid body

  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  double massone;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    massone = mass[type[i]];
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

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;
    dx = x[i][0] + xbox*xprd - xcm[ibody][0];
    dy = x[i][1] + ybox*yprd - xcm[ibody][1];
    dz = x[i][2] + zbox*zprd - xcm[ibody][2];
    
    massone = mass[type[i]];
    sum[ibody][0] += dy * massone*v[i][2] - dz * massone*v[i][1];
    sum[ibody][1] += dz * massone*v[i][0] - dx * massone*v[i][2];
    sum[ibody][2] += dx * massone*v[i][1] - dy * massone*v[i][0];
    sum[ibody][3] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][4] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][5] += dx * f[i][1] - dy * f[i][0];
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

  // set velocities from angmom & omega
  // guestimate virial as 2x the set_v contribution

  for (ibody = 0; ibody < nbody; ibody++)
    omega_from_mq(angmom[ibody],ex_space[ibody],ey_space[ibody],
		  ez_space[ibody],inertia[ibody],omega[ibody]);

  for (int n = 0; n < 6; n++) virial[n] = 0.0;
  set_v(1);
  for (int n = 0; n < 6; n++) virial[n] *= 2.0;
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate()
{
  double dtfm;

  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step
  
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2];
  
    // update xcm by full step
  
    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];
  
    // update angular momentum by 1/2 step
    
    angmom[ibody][0] += dtf * torque[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2];

    // update quaternion a full step
    // returns new normalized quat
    // returns ex_space,ey_space,ez_space for new quat
    // returns omega at 1/2 step (depends on angmom and quat)

    richardson(quat[ibody],omega[ibody],angmom[ibody],inertia[ibody],
	       ex_space[ibody],ey_space[ibody],ez_space[ibody]);
  }

  // set coords and velocities if atoms in rigid bodies
  // from quarternion and omega
  
  int vflag = 0;
  if (pressure_flag || output->next_thermo == update->ntimestep) vflag = 1;
  set_xv(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRigid::richardson(double *q, double *w,
			  double *m, double *moments,
			  double *ex, double *ey, double *ez)
{
  // compute omega at 1/2 step from m at 1/2 step and q at 0
  
  omega_from_mq(m,ex,ey,ez,moments,w);

  // full update from dq/dt = 1/2 w q

  double wq[4];
  multiply(w,q,wq);

  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];

  normalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];

  normalize(qhalf);

  // udpate ex,ey,ez from qhalf
  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  exyz_from_q(qhalf,ex,ey,ez);
  omega_from_mq(m,ex,ey,ez,moments,w);
  multiply(w,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];

  normalize(qhalf);

  // corrected Richardson update

  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];

  normalize(q);
  exyz_from_q(q,ex,ey,ez);
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate()
{
  int i,ibody;
  double dtfm;

  // sum forces and torques on atoms in rigid body
  
  int *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  int xbox,ybox,zbox;
  double dx,dy,dz;
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
    dx = x[i][0] + xbox*xprd - xcm[ibody][0];
    dy = x[i][1] + ybox*yprd - xcm[ibody][1];
    dz = x[i][2] + zbox*zprd - xcm[ibody][2];
    
    sum[ibody][3] += dy*f[i][2] - dz*f[i][1];
    sum[ibody][4] += dz*f[i][0] - dx*f[i][2];
    sum[ibody][5] += dx*f[i][1] - dy*f[i][0];
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
    vcm[ibody][0] += dtfm * fcm[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2];
  
    // update angular momentum by 1/2 step
    
    angmom[ibody][0] += dtf * torque[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2];
  }

  // set velocities from angmom & omega
  
  for (ibody = 0; ibody < nbody; ibody++) 
    omega_from_mq(angmom[ibody],ex_space[ibody],ey_space[ibody],
		  ez_space[ibody],inertia[ibody],omega[ibody]);

  int vflag = 0;
  if (pressure_flag || output->next_thermo == update->ntimestep) vflag = 1;
  set_v(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate_respa(int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtq = 0.5 * step_respa[ilevel];

  if (ilevel == 0) initial_integrate();
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate_respa(int ilevel)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_rigid for atoms in igroup 
------------------------------------------------------------------------- */

int FixRigid::dof(int igroup)
{
  int groupbit = group->bitmask[igroup];

  // ncount = # of atoms in each rigid body that are also in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *ncount = new int[nbody];
  for (int ibody = 0; ibody < nbody; ibody++) ncount[ibody] = 0;

  for (int i = 0; i < nlocal; i++)
    if (body[i] >= 0 && mask[i] & groupbit) ncount[body[i]]++;

  int *nall = new int[nbody];
  MPI_Allreduce(ncount,nall,nbody,MPI_INT,MPI_SUM,world);

  // remove 3N - 6 dof for each rigid body if more than 2 atoms in igroup
  // remove 3N - 5 dof for each diatomic rigid body in igroup

  int n = 0;
  for (int ibody = 0; ibody < nbody; ibody++) {
    if (nall[ibody] > 2) n += 3*nall[ibody] - 6;
    else if (nall[ibody] == 2) n++;
  }

  delete [] ncount;
  delete [] nall;

  return n;
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box dilation in idim
   called by various fixes that change box
------------------------------------------------------------------------- */

void FixRigid::dilate(int idim,
		      double oldlo, double oldhi, double newlo, double newhi)
{
  double ratio;
  for (int ibody = 0; ibody < nbody; ibody++) {
    ratio = (xcm[ibody][idim] - oldlo) / (oldhi - oldlo);
    xcm[ibody][idim] = newlo + ratio*(newhi - newlo);
  }
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
   create quaternion from rotation matrix (evectors)
------------------------------------------------------------------------- */

void FixRigid::qcreate(double **a, double *q)
{
  // squares of quaternion components
  
  double q0sq = 0.25 * (a[0][0] + a[1][1] + a[2][2] + 1.0);
  double q1sq = q0sq - 0.5 * (a[1][1] + a[2][2]);
  double q2sq = q0sq - 0.5 * (a[0][0] + a[2][2]);
  double q3sq = q0sq - 0.5 * (a[0][0] + a[1][1]);
  
  // some component must be greater than 1/4 since they sum to 1
  // compute other components from it
  
  if (q0sq >= 0.25) {
    q[0] = sqrt(q0sq);
    q[1] = (a[2][1] - a[1][2]) / (4.0*q[0]);
    q[2] = (a[0][2] - a[2][0]) / (4.0*q[0]);
    q[3] = (a[1][0] - a[0][1]) / (4.0*q[0]);
  } else if (q1sq >= 0.25) {
    q[1] = sqrt(q1sq);
    q[0] = (a[2][1] - a[1][2]) / (4.0*q[1]);
    q[2] = (a[0][1] + a[1][0]) / (4.0*q[1]);
    q[3] = (a[2][0] + a[0][2]) / (4.0*q[1]);
  } else if (q2sq >= 0.25) {
    q[2] = sqrt(q2sq);
    q[0] = (a[0][2] - a[2][0]) / (4.0*q[2]);
    q[1] = (a[0][1] + a[1][0]) / (4.0*q[2]);
    q[2] = (a[1][2] + a[2][1]) / (4.0*q[2]);
  } else if (q3sq >= 0.25) {
    q[3] = sqrt(q3sq);
    q[0] = (a[1][0] - a[0][1]) / (4.0*q[3]);
    q[1] = (a[0][2] + a[2][0]) / (4.0*q[3]);
    q[2] = (a[1][2] + a[2][1]) / (4.0*q[3]);
  } else
    error->all("Quaternion creation numeric error");

  normalize(q);
}

/* ----------------------------------------------------------------------
   quaternion multiply: c = a*b where a = (0,a)
------------------------------------------------------------------------- */

void FixRigid::multiply(double *a, double *b, double *c)
{
  c[0] = -(a[0]*b[1] + a[1]*b[2] + a[2]*b[3]);
  c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
  c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
  c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

void FixRigid::normalize(double *q)
{
  double norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in body frame
   set wbody component to 0.0 if inertia component is 0.0
     otherwise body can spin easily around that axis
   project space-frame angular momentum onto body axes
     and divide by principal moments
------------------------------------------------------------------------- */

void FixRigid::omega_from_mq(double *m, double *ex, double *ey, double *ez,
			     double *inertia, double *w)
{
  double wbody[3];

  if (inertia[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / inertia[0];
  if (inertia[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / inertia[1];
  if (inertia[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / inertia[2];

  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
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
   set space-frame coords and velocity of each atom in each rigid body
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_xv(int vflag)
{
  int *image = atom->image;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  int ibody;
  int xbox,ybox,zbox;

  double vold0,vold1,vold2,fc0,fc1,fc2,massone,x0,x1,x2;
  double *mass = atom->mass; 
  double **f = atom->f;
  int *type = atom->type;
  
  // zero out fix_rigid virial

  if (vflag) for (int n = 0; n < 6; n++) virial[n] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    // save old positions and velocities for virial contribution

    if (vflag) {
      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      vold0 = v[i][0];
      vold1 = v[i][1];
      vold2 = v[i][2];
    }

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
    
    x[i][0] += xcm[ibody][0] - xbox*xprd;
    x[i][1] += xcm[ibody][1] - ybox*yprd;
    x[i][2] += xcm[ibody][2] - zbox*zprd;

    // compute body constraint forces for virial

    if (vflag) {
      massone = mass[type[i]];
      fc0 = massone*(v[i][0] - vold0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - vold1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - vold2)/dtf - f[i][2]; 

      virial[0] += 0.5*fc0*x0;
      virial[1] += 0.5*fc1*x1;
      virial[2] += 0.5*fc2*x2;
      virial[3] += 0.5*fc1*x0;
      virial[4] += 0.5*fc2*x0;
      virial[5] += 0.5*fc2*x1;
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in rigid body
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_v(int vflag)
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  int ibody;
  double dx,dy,dz;

  double vold0,vold1,vold2,fc0,fc1,fc2,massone,x0,x1,x2;
  double *mass = atom->mass; 
  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int *image = atom->image;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  int xbox,ybox,zbox;

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

    if (vflag) {
      vold0 = v[i][0];
      vold1 = v[i][1];
      vold2 = v[i][2];
    }

    v[i][0] = omega[ibody][1]*dz - omega[ibody][2]*dy + vcm[ibody][0];
    v[i][1] = omega[ibody][2]*dx - omega[ibody][0]*dz + vcm[ibody][1];
    v[i][2] = omega[ibody][0]*dy - omega[ibody][1]*dx + vcm[ibody][2];

    // compute body constraint forces for virial
    // use unwrapped atom positions

    if (vflag) {
      massone = mass[type[i]];
      fc0 = massone*(v[i][0] - vold0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - vold1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - vold2)/dtf - f[i][2]; 

      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;

      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      virial[0] += 0.5*fc0*x0;
      virial[1] += 0.5*fc1*x1;
      virial[2] += 0.5*fc2*x2;
      virial[3] += 0.5*fc1*x0;
      virial[4] += 0.5*fc2*x0;
      virial[5] += 0.5*fc2*x1;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

int FixRigid::memory_usage()
{
  int nmax = atom->nmax;
  int bytes = nmax * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  return bytes;
} 

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixRigid::grow_arrays(int nmax)
{
  body = (int *) memory->srealloc(body,nmax*sizeof(int),"rigid:body");
  displace = memory->grow_2d_double_array(displace,nmax,3,"rigid:displace");
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
  return 4;
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
  return 4;
}
