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
#include "fix_rigid_small.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "molecule.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "random_mars.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// allocate space for static class variable

FixRigidSmall *FixRigidSmall::frsptr;

#define MAXLINE 256
#define CHUNK 1024
#define ATTRIBUTE_PERBODY 11

#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-7
#define BIG 1.0e20

#define SINERTIA 0.4            // moment of inertia prefactor for sphere
#define EINERTIA 0.4            // moment of inertia prefactor for ellipsoid
#define LINERTIA (1.0/12.0)     // moment of inertia prefactor for line segment

#define DELTA_BODY 10000

enum{FULL_BODY,INITIAL,FINAL,FORCE_TORQUE,VCM_ANGMOM,XCM_MASS,ITENSOR,DOF};

/* ---------------------------------------------------------------------- */

FixRigidSmall::FixRigidSmall(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int i,ibody;

  scalar_flag = 1;
  extscalar = 0;
  global_freq = 1;
  time_integrate = 1;
  rigid_flag = 1;
  virial_flag = 1;
  create_attribute = 1;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  extended = orientflag = dorientflag = 0;
  bodyown = NULL;
  bodytag = NULL;
  atom2body = NULL;
  displace = NULL;
  eflags = NULL;
  orient = NULL;
  dorient = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // parse args for rigid body specification

  if (narg < 4) error->all(FLERR,"Illegal fix rigid/small command");
  if (strcmp(arg[3],"molecule") != 0) 
    error->all(FLERR,"Illegal fix rigid/small command");

  if (atom->molecule_flag == 0)
    error->all(FLERR,"Fix rigid/small requires atom attribute molecule");
  if (atom->map_style == 0)
    error->all(FLERR,"Fix rigid/small requires an atom map, see atom_modify");

  // maxmol = largest molecule #

  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  maxmol = -1;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) maxmol = MAX(maxmol,molecule[i]);

  tagint itmp;
  MPI_Allreduce(&maxmol,&itmp,1,MPI_LMP_TAGINT,MPI_MAX,world);
  maxmol = itmp;

  // parse optional args

  int seed;
  langflag = 0;
  infile = NULL;
  onemol = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"langevin") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix rigid/small command");
      if (strcmp(style,"rigid/small") != 0)
        error->all(FLERR,"Illegal fix rigid/small command");
      langflag = 1;
      t_start = force->numeric(FLERR,arg[iarg+1]);
      t_stop = force->numeric(FLERR,arg[iarg+2]);
      t_period = force->numeric(FLERR,arg[iarg+3]);
      seed = force->inumeric(FLERR,arg[iarg+4]);
      if (t_period <= 0.0)
        error->all(FLERR,"Fix rigid/small langevin period must be > 0.0");
      if (seed <= 0) error->all(FLERR,"Illegal fix rigid/small command");
      iarg += 5;
    } else if (strcmp(arg[iarg],"infile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid/small command");
      delete [] infile;
      int n = strlen(arg[iarg+1]) + 1;
      infile = new char[n];
      strcpy(infile,arg[iarg+1]);
      restart_file = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid/small command");
      int imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule template ID for "
                   "fix rigid/small does not exist");
      if (atom->molecules[imol]->nset > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template for "
                       "fix rigid/small has multiple molecules");
      onemol = atom->molecules[imol];
      iarg += 2;
    } else error->all(FLERR,"Illegal fix rigid/small command");
  }

  // error check and further setup for Molecule template

  if (onemol) {
    if (onemol->xflag == 0)
      error->all(FLERR,"Fix rigid/small molecule must have coordinates");
    if (onemol->typeflag == 0)
      error->all(FLERR,"Fix rigid/small molecule must have atom types");

    // fix rigid/small uses center, masstotal, COM, inertia of molecule

    onemol->compute_center();
    onemol->compute_mass();
    onemol->compute_com();
    onemol->compute_inertia();
  }

  // create rigid bodies based on molecule ID
  // sets bodytag for owned atoms
  // body attributes are computed later by setup_bodies()

  create_bodies();

  // set nlocal_body and allocate bodies I own

  tagint *tag = atom->tag;

  nlocal_body = nghost_body = 0;
  for (i = 0; i < nlocal; i++)
    if (bodytag[i] == tag[i]) nlocal_body++;

  nmax_body = 0;
  while (nmax_body < nlocal_body) nmax_body += DELTA_BODY;
  body = (Body *) memory->smalloc(nmax_body*sizeof(Body),
                                  "rigid/small:body");

  // set bodyown for owned atoms

  nlocal_body = 0;
  for (i = 0; i < nlocal; i++)
    if (bodytag[i] == tag[i]) {
      body[nlocal_body].ilocal = i;
      bodyown[i] = nlocal_body++;
    } else bodyown[i] = -1;


  // bodysize = sizeof(Body) in doubles

  bodysize = sizeof(Body)/sizeof(double);
  if (bodysize*sizeof(double) != sizeof(Body)) bodysize++;

  // set max comm sizes needed by this fix

  comm_forward = 1 + bodysize;
  comm_reverse = 6;

  // bitmasks for properties of extended particles

  POINT = 1;
  SPHERE = 2;
  ELLIPSOID = 4;
  LINE = 8;
  TRIANGLE = 16;
  DIPOLE = 32;
  OMEGA = 64;
  ANGMOM = 128;
  TORQUE = 256;

  MINUSPI = -MY_PI;
  TWOPI = 2.0*MY_PI;

  // atom style pointers to particles that store extra info

  avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  avec_line = (AtomVecLine *) atom->style_match("line");
  avec_tri = (AtomVecTri *) atom->style_match("tri");

  // print statistics

  int one = 0;
  bigint atomone = 0;
  for (int i = 0; i < nlocal; i++) {
    if (bodyown[i] >= 0) one++;
    if (bodytag[i] > 0) atomone++;
  }
  MPI_Allreduce(&one,&nbody,1,MPI_INT,MPI_SUM,world);
  bigint atomall;
  MPI_Allreduce(&atomone,&atomall,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"%d rigid bodies with " BIGINT_FORMAT " atoms\n",
              nbody,atomall);
      fprintf(screen,"  %g = max distance from body owner to body atom\n",
              maxextent);
    }
    if (logfile) {
      fprintf(logfile,"%d rigid bodies with " BIGINT_FORMAT " atoms\n",
              nbody,atomall);
      fprintf(logfile,"  %g = max distance from body owner to body atom\n",
              maxextent);
    }
  }

  // initialize Marsaglia RNG with processor-unique seed

  maxlang = 0;
  langextra = NULL;
  random = NULL;
  if (langflag) random = new RanMars(lmp,seed + comm->me);

  // mass vector for granular pair styles

  mass_body = NULL;
  nmax_mass = 0;

  // firstflag = 1 triggers one-time initialization of rigid body attributes

  firstflag = 1;
}

/* ---------------------------------------------------------------------- */

FixRigidSmall::~FixRigidSmall()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->sfree(body);
  
  memory->destroy(bodyown);
  memory->destroy(bodytag);
  memory->destroy(atom2body);
  memory->destroy(displace);
  memory->destroy(eflags);
  memory->destroy(orient);
  memory->destroy(dorient);

  delete random;
  delete [] infile;

  memory->destroy(langextra);
  memory->destroy(mass_body);
}

/* ---------------------------------------------------------------------- */

int FixRigidSmall::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  if (langflag) mask |= POST_FORCE;
  mask |= PRE_NEIGHBOR;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigidSmall::init()
{
  int i,ibody;

  triclinic = domain->triclinic;

  // warn if more than one rigid fix

  int count = 0;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"rigid") == 0) count++;
  if (count > 1 && me == 0) error->warning(FLERR,"More than one fix rigid");

  // error if npt,nph fix comes before rigid fix

  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) break;
    if (strcmp(modify->fix[i]->style,"nph") == 0) break;
  }
  if (i < modify->nfix) {
    for (int j = i; j < modify->nfix; j++)
      if (strcmp(modify->fix[j]->style,"rigid") == 0)
        error->all(FLERR,"Rigid fix must come before NPT/NPH fix");
  }

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes via local comm
   extended flags, mass, COM, inertia tensor, displacement of each atom
   performed after init() b/c requires communication stencil
     has been setup by comm->borders()
------------------------------------------------------------------------- */

void FixRigidSmall::setup_pre_neighbor()
{
  if (firstflag) {
    firstflag = 0;
    setup_bodies_static();
    setup_bodies_dynamic();
  } else pre_neighbor();
}

/* ----------------------------------------------------------------------
   compute initial fcm and torque on bodies, also initial virial
   reset all particle velocities to be consistent with vcm and omega
------------------------------------------------------------------------- */

void FixRigidSmall::setup(int vflag)
{
  int i,n,ibody;
  double massone,radone;

  //check(1);

  // sum fcm, torque across all rigid bodies
  // fcm = force on COM
  // torque = torque around COM

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *xcm,*fcm,*tcm;
  double dx,dy,dz;
  double unwrap[3];

  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    fcm = body[ibody].fcm;
    fcm[0] = fcm[1] = fcm[2] = 0.0;
    tcm = body[ibody].torque;
    tcm[0] = tcm[1] = tcm[2] = 0.0;
  }

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    fcm = b->fcm;
    fcm[0] += f[i][0];
    fcm[1] += f[i][1];
    fcm[2] += f[i][2];

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    dx = unwrap[0] - xcm[0];
    dy = unwrap[1] - xcm[1];
    dz = unwrap[2] - xcm[2];

    tcm = b->torque;
    tcm[0] += dy * f[i][2] - dz * f[i][1];
    tcm[1] += dz * f[i][0] - dx * f[i][2];
    tcm[2] += dx * f[i][1] - dy * f[i][0];
  }

  // extended particles add their rotation/torque to angmom/torque of body

  if (extended) {
    double **torque = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];
      if (eflags[i] & TORQUE) {
        tcm = b->torque;
        tcm[0] += torque[i][0];
        tcm[1] += torque[i][1];
        tcm[2] += torque[i][2];
      }
    }
  }

  // reverse communicate fcm, torque of all bodies

  commflag = FORCE_TORQUE;
  comm->reverse_comm_variable_fix(this);

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // compute and forward communicate vcm and omega of all bodies

  for (ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];
    MathExtra::angmom_to_omega(b->angmom,b->ex_space,b->ey_space,
                               b->ez_space,b->inertia,b->omega);
  }

  commflag = FINAL;
  comm->forward_comm_variable_fix(this);

  // set velocity/rotation of atoms in rigid bodues

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

void FixRigidSmall::initial_integrate(int vflag)
{
  double dtfm;

  //check(2);

  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];

    // update vcm by 1/2 step

    dtfm = dtf / b->mass;
    b->vcm[0] += dtfm * b->fcm[0];
    b->vcm[1] += dtfm * b->fcm[1];
    b->vcm[2] += dtfm * b->fcm[2];

    // update xcm by full step

    b->xcm[0] += dtv * b->vcm[0];
    b->xcm[1] += dtv * b->vcm[1];
    b->xcm[2] += dtv * b->vcm[2];

    // update angular momentum by 1/2 step

    b->angmom[0] += dtf * b->torque[0];
    b->angmom[1] += dtf * b->torque[1];
    b->angmom[2] += dtf * b->torque[2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(b->angmom,b->ex_space,b->ey_space,
                               b->ez_space,b->inertia,b->omega);
    MathExtra::richardson(b->quat,b->angmom,b->omega,b->inertia,dtq);
    MathExtra::q_to_exyz(b->quat,b->ex_space,b->ey_space,b->ez_space);
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  comm->forward_comm_variable_fix(this);

  // set coords/orient and velocity/rotation of atoms in rigid bodies

  set_xv();
}

/* ----------------------------------------------------------------------
   apply Langevin thermostat to all 6 DOF of rigid bodies I own
   unlike fix langevin, this stores extra force in extra arrays,
     which are added in when final_integrate() calculates a new fcm/torque
------------------------------------------------------------------------- */

void FixRigidSmall::post_force(int vflag)
{
  double gamma1,gamma2;

  // grow langextra if needed

  if (nlocal_body > maxlang) {
    memory->destroy(langextra);
    maxlang = nlocal_body + nghost_body;
    memory->create(langextra,maxlang,6,"rigid/small:langextra");
  }

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);
  double tsqrt = sqrt(t_target);

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double *vcm,*omega,*inertia;

  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    vcm = body[ibody].vcm;
    omega = body[ibody].omega;
    inertia = body[ibody].inertia;

    gamma1 = -body[ibody].mass / t_period / ftm2v;
    gamma2 = sqrt(body[ibody].mass) * tsqrt *
      sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
    langextra[ibody][0] = gamma1*vcm[0] + gamma2*(random->uniform()-0.5);
    langextra[ibody][1] = gamma1*vcm[1] + gamma2*(random->uniform()-0.5);
    langextra[ibody][2] = gamma1*vcm[2] + gamma2*(random->uniform()-0.5);
    
    gamma1 = -1.0 / t_period / ftm2v;
    gamma2 = tsqrt * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
    langextra[ibody][3] = inertia[0]*gamma1*omega[0] +
      sqrt(inertia[0])*gamma2*(random->uniform()-0.5);
    langextra[ibody][4] = inertia[1]*gamma1*omega[1] +
      sqrt(inertia[1])*gamma2*(random->uniform()-0.5);
    langextra[ibody][5] = inertia[2]*gamma1*omega[2] +
      sqrt(inertia[2])*gamma2*(random->uniform()-0.5);
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidSmall::final_integrate()
{
  int i,ibody;
  double dtfm,xy,xz,yz;

  //check(3);

  // sum over atoms to get force and torque on rigid body

  imageint *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];
  double *xcm,*fcm,*tcm;

  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    fcm = body[ibody].fcm;
    fcm[0] = fcm[1] = fcm[2] = 0.0;
    tcm = body[ibody].torque;
    tcm[0] = tcm[1] = tcm[2] = 0.0;
  }

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    fcm = b->fcm;
    fcm[0] += f[i][0];
    fcm[1] += f[i][1];
    fcm[2] += f[i][2];

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    dx = unwrap[0] - xcm[0];
    dy = unwrap[1] - xcm[1];
    dz = unwrap[2] - xcm[2];

    tcm = b->torque;
    tcm[0] += dy*f[i][2] - dz*f[i][1];
    tcm[1] += dz*f[i][0] - dx*f[i][2];
    tcm[2] += dx*f[i][1] - dy*f[i][0];
  }

  // extended particles add their torque to torque of body

  if (extended) {
    double **torque = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;

      if (eflags[i] & TORQUE) {
        tcm = body[atom2body[i]].torque;
        tcm[0] += torque[i][0];
        tcm[1] += torque[i][1];
        tcm[2] += torque[i][2];
      }
    }
  }

  // reverse communicate fcm, torque of all bodies

  commflag = FORCE_TORQUE;
  comm->reverse_comm_variable_fix(this);

  // include Langevin thermostat forces and torques

  if (langflag) {
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      fcm = body[ibody].fcm;
      fcm[0] += langextra[ibody][0];
      fcm[1] += langextra[ibody][1];
      fcm[2] += langextra[ibody][2];
      tcm = body[ibody].torque;
      tcm[0] += langextra[ibody][3];
      tcm[1] += langextra[ibody][4];
      tcm[2] += langextra[ibody][5];
    }
  }

  // update vcm and angmom, recompute omega
  
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];

    // update vcm by 1/2 step

    dtfm = dtf / b->mass;
    b->vcm[0] += dtfm * b->fcm[0];
    b->vcm[1] += dtfm * b->fcm[1];
    b->vcm[2] += dtfm * b->fcm[2];

    // update angular momentum by 1/2 step

    b->angmom[0] += dtf * b->torque[0];
    b->angmom[1] += dtf * b->torque[1];
    b->angmom[2] += dtf * b->torque[2];

    MathExtra::angmom_to_omega(b->angmom,b->ex_space,b->ey_space,
                               b->ez_space,b->inertia,b->omega);
  }

  // forward communicate updated info of all bodies

  commflag = FINAL;
  comm->forward_comm_variable_fix(this);

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}

/* ---------------------------------------------------------------------- */

void FixRigidSmall::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtq = 0.5 * step_respa[ilevel];

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixRigidSmall::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ----------------------------------------------------------------------
   atom to processor assignments have changed,
     so acquire ghost bodies and setup atom2body
   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to 1st definition of rigid body or due to box flip
   if don't do this, then atoms of a body which drifts far away
     from a triclinic box will be remapped back into box
     with huge displacements when the box tilt changes via set_x()
   adjust image flag of body and image flags of all atoms in body
------------------------------------------------------------------------- */

void FixRigidSmall::pre_neighbor()
{
  // remap xcm and image flags of each body as needed

  imageint original,oldimage,newimage;

  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body *b = &body[ibody];

    original = b->image;
    domain->remap(b->xcm,b->image);

    if (original == b->image) b->remapflag[3] = 0;
    else {
      oldimage = original & IMGMASK;
      newimage = b->image & IMGMASK;
      b->remapflag[0] = newimage - oldimage;
      oldimage = (original >> IMGBITS) & IMGMASK;
      newimage = (b->image >> IMGBITS) & IMGMASK;
      b->remapflag[1] = newimage - oldimage;
      oldimage = original >> IMG2BITS;
      newimage = b->image >> IMG2BITS;
      b->remapflag[2] = newimage - oldimage;
      b->remapflag[3] = 1;
    }
  }

  // acquire ghost bodies via forward comm
  // also gets new remapflags needed for adjusting atom image flags
  // reset atom2body for owned atoms
  // forward comm sets atom2body for ghost atoms

  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm_variable_fix(this);
  reset_atom2body();
  //check(4);

  // adjust image flags of any atom in a rigid body whose xcm was remapped

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  int ibody;
  imageint idim,otherdims;

  for (int i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    if (b->remapflag[3] == 0) continue;

    if (b->remapflag[0]) {
      idim = image[i] & IMGMASK;
      otherdims = image[i] ^ idim;
      idim -= b->remapflag[0];
      idim &= IMGMASK;
      image[i] = otherdims | idim;
    }
    if (b->remapflag[1]) {
      idim = (image[i] >> IMGBITS) & IMGMASK;
      otherdims = image[i] ^ (idim << IMGBITS);
      idim -= b->remapflag[1];
      idim &= IMGMASK;
      image[i] = otherdims | (idim << IMGBITS);
    }
    if (b->remapflag[2]) {
      idim = image[i] >> IMG2BITS;
      otherdims = image[i] ^ (idim << IMG2BITS);
      idim -= b->remapflag[2];
      idim &= IMGMASK;
      image[i] = otherdims | (idim << IMG2BITS);
    }
  }
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

int FixRigidSmall::dof(int tgroup)
{
  int i,j;

  // cannot count DOF correctly unless setup_bodies_static() has been called

  if (firstflag) {
    if (comm->me == 0) 
      error->warning(FLERR,"Cannot count rigid body degrees-of-freedom "
                     "before bodies are fully initialized");
    return 0;
  }

  int tgroupbit = group->bitmask[tgroup];

  // counts = 3 values per rigid body I own
  // 0 = # of point particles in rigid body and in temperature group
  // 1 = # of finite-size particles in rigid body and in temperature group
  // 2 = # of particles in rigid body, disregarding temperature group

  memory->create(counts,nlocal_body+nghost_body,3,"rigid/small:counts");
  for (int i = 0; i < nlocal_body+nghost_body; i++)
    counts[i][0] = counts[i][1] = counts[i][2] = 0;

  // tally counts from my owned atoms
  // 0 = # of point particles in rigid body and in temperature group
  // 1 = # of finite-size particles in rigid body and in temperature group
  // 2 = # of particles in rigid body, disregarding temperature group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    j = atom2body[i];
    counts[j][2]++;
    if (mask[i] & tgroupbit) {
      if (extended && eflags[i]) counts[j][1]++;
      else counts[j][0]++;
    }
  }

  commflag = DOF;
  comm->reverse_comm_variable_fix(this);

  // nall = count0 = # of point particles in each rigid body
  // mall = count1 = # of finite-size particles in each rigid body
  // warn if nall+mall != nrigid for any body included in temperature group

  int flag = 0;
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    if (counts[ibody][0]+counts[ibody][1] > 0 &&
        counts[ibody][0]+counts[ibody][1] != counts[ibody][2]) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && me == 0)
    error->warning(FLERR,"Computing temperature of portions of rigid bodies");

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

  double *inertia;

  int n = 0;
  if (domain->dimension == 3) {
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      if (counts[ibody][0]+counts[ibody][1] == counts[ibody][2]) {
        n += 3*counts[ibody][0] + 6*counts[ibody][1] - 6;
        inertia = body[ibody].inertia;
        if (inertia[0] == 0.0 || inertia[1] == 0.0 || inertia[2] == 0.0) n++;
      }
    }
  } else if (domain->dimension == 2) {
    for (int ibody = 0; ibody < nlocal_body; ibody++)
      if (counts[ibody][0]+counts[ibody][1] == counts[ibody][2])
        n += 2*counts[ibody][0] + 3*counts[ibody][1] - 3;
  }

  memory->destroy(counts);

  int nall;
  MPI_Allreduce(&n,&nall,1,MPI_INT,MPI_SUM,world);
  return nall;
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

void FixRigidSmall::deform(int flag)
{
  if (flag == 0)
    for (int ibody = 0; ibody < nlocal_body; ibody++)
      domain->x2lamda(body[ibody].xcm,body[ibody].xcm);
  else
    for (int ibody = 0; ibody < nlocal_body; ibody++)
      domain->lamda2x(body[ibody].xcm,body[ibody].xcm);
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigidSmall::set_xv()
{
  int ibody,itype;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double ione[3],exone[3],eyone[3],ezone[3],vr[6],p[3][3];

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  imageint *image = atom->image;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // set x and v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

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

    MathExtra::matvec(b->ex_space,b->ey_space,b->ez_space,displace[i],x[i]);

    v[i][0] = b->omega[1]*x[i][2] - b->omega[2]*x[i][1] + b->vcm[0];
    v[i][1] = b->omega[2]*x[i][0] - b->omega[0]*x[i][2] + b->vcm[1];
    v[i][2] = b->omega[0]*x[i][1] - b->omega[1]*x[i][0] + b->vcm[2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    if (triclinic == 0) {
      x[i][0] += b->xcm[0] - xbox*xprd;
      x[i][1] += b->xcm[1] - ybox*yprd;
      x[i][2] += b->xcm[2] - zbox*zprd;
    } else {
      x[i][0] += b->xcm[0] - xbox*xprd - ybox*xy - zbox*xz;
      x[i][1] += b->xcm[1] - ybox*yprd - zbox*yz;
      x[i][2] += b->xcm[2] - zbox*zprd;
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
    double theta_body,theta;
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega = atom->omega;
    double **angmom = atom->angmom;
    double **mu = atom->mu;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega[i][0] = b->omega[0];
        omega[i][1] = b->omega[1];
        omega[i][2] = b->omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::quatquat(b->quat,orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b->omega,exone,eyone,ezone,ione,angmom[i]);
      } else if (eflags[i] & LINE) {
        if (b->quat[3] >= 0.0) theta_body = 2.0*acos(b->quat[0]);
        else theta_body = -2.0*acos(b->quat[0]);
        theta = orient[i][0] + theta_body;
        while (theta <= MINUSPI) theta += TWOPI;
        while (theta > MY_PI) theta -= TWOPI;
        lbonus[line[i]].theta = theta;
        omega[i][0] = b->omega[0];
        omega[i][1] = b->omega[1];
        omega[i][2] = b->omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::quatquat(b->quat,orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b->omega,exone,eyone,ezone,
                                   inertiaatom,angmom[i]);
      }
      if (eflags[i] & DIPOLE) {
        MathExtra::quat_to_mat(b->quat,p);
        MathExtra::matvec(p,dorient[i],mu[i]);
        MathExtra::snormalize3(mu[i][3],mu[i],mu[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigidSmall::set_v()
{
  int ibody,itype;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double ione[3],exone[3],eyone[3],ezone[3],delta[3],vr[6];

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    MathExtra::matvec(b->ex_space,b->ey_space,b->ez_space,displace[i],delta);

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = b->omega[1]*delta[2] - b->omega[2]*delta[1] + b->vcm[0];
    v[i][1] = b->omega[2]*delta[0] - b->omega[0]*delta[2] + b->vcm[1];
    v[i][2] = b->omega[0]*delta[1] - b->omega[1]*delta[0] + b->vcm[2];

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

      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;

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
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega = atom->omega;
    double **angmom = atom->angmom;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega[i][0] = b->omega[0];
        omega[i][1] = b->omega[1];
        omega[i][2] = b->omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b->omega,exone,eyone,ezone,ione,
                                   angmom[i]);
      } else if (eflags[i] & LINE) {
        omega[i][0] = b->omega[0];
        omega[i][1] = b->omega[1];
        omega[i][2] = b->omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b->omega,exone,eyone,ezone,
                                   inertiaatom,angmom[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   one-time identification of which atoms are in which rigid bodies
   set bodytag for all owned atoms
------------------------------------------------------------------------- */

void FixRigidSmall::create_bodies()
{
  int i,m,n;
  double unwrap[3];

  // error check on image flags of atoms in rigid bodies

  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *periodicity = domain->periodicity;
  int xbox,ybox,zbox;

  int flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;
    if ((xbox && !periodicity[0]) || (ybox && !periodicity[1]) ||
        (zbox && !periodicity[2])) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"Fix rigid/small atom has non-zero image flag "
                          "in a non-periodic dimension");

  // allocate buffer for passing messages around ring of procs
  // percount = max number of values to put in buffer for each of ncount

  int ncount = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) ncount++;

  int percount = 5;
  double *buf;
  memory->create(buf,ncount*percount,"rigid/small:buf");

  // create map hash for storing unique molecule IDs of my atoms
  // key = molecule ID
  // value = index into per-body data structure
  // n = # of entries in hash

  hash = new std::map<tagint,int>();
  hash->clear();

  // setup hash
  // key = body ID
  // value = index into N-length data structure
  // n = count of unique bodies my atoms are part of

  tagint *molecule = atom->molecule;

  n = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (hash->find(molecule[i]) == hash->end()) (*hash)[molecule[i]] = n++;
  }

  // bbox = bounding box of each rigid body my atoms are part of

  memory->create(bbox,n,6,"rigid/small:bbox");

  for (i = 0; i < n; i++) {
    bbox[i][0] = bbox[i][2] = bbox[i][4] = BIG;
    bbox[i][1] = bbox[i][3] = bbox[i][5] = -BIG;
  }

  // pack my atoms into buffer as molecule ID, unwrapped coords

  double **x = atom->x;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    domain->unmap(x[i],image[i],unwrap);
    buf[m++] = molecule[i];
    buf[m++] = unwrap[0];
    buf[m++] = unwrap[1];
    buf[m++] = unwrap[2];
  }

  // pass buffer around ring of procs
  // func = update bbox with atom coords from every proc
  // when done, have full bbox for every rigid body my atoms are part of

  frsptr = this;
  comm->ring(m,sizeof(double),buf,1,ring_bbox,NULL);

  // ctr = center pt of each rigid body my atoms are part of

  memory->create(ctr,n,6,"rigid/small:bbox");

  for (i = 0; i < n; i++) {
    ctr[i][0] = 0.5 * (bbox[i][0] + bbox[i][1]);
    ctr[i][1] = 0.5 * (bbox[i][2] + bbox[i][3]);
    ctr[i][2] = 0.5 * (bbox[i][4] + bbox[i][5]);
  }

  // idclose = ID of atom in body closest to center pt (smaller ID if tied)
  // rsqclose = distance squared from idclose to center pt

  memory->create(idclose,n,"rigid/small:idclose");
  memory->create(rsqclose,n,"rigid/small:rsqclose");

  for (i = 0; i < n; i++) rsqclose[i] = BIG;

  // pack my atoms into buffer as molecule ID, atom ID, unwrapped coords

  tagint *tag = atom->tag;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    domain->unmap(x[i],image[i],unwrap);
    buf[m++] = molecule[i];
    buf[m++] = ubuf(tag[i]).d;
    buf[m++] = unwrap[0];
    buf[m++] = unwrap[1];
    buf[m++] = unwrap[2];
  }

  // pass buffer around ring of procs
  // func = update idclose,rsqclose with atom IDs from every proc
  // when done, have idclose for every rigid body my atoms are part of

  frsptr = this;
  comm->ring(m,sizeof(double),buf,2,ring_nearest,NULL);

  // set bodytag of all owned atoms, based on idclose
  // find max value of rsqclose across all procs

  double rsqmax = 0.0;
  for (i = 0; i < nlocal; i++) {
    bodytag[i] = 0;
    if (!(mask[i] & groupbit)) continue;
    m = hash->find(molecule[i])->second;
    bodytag[i] = idclose[m];
    rsqmax = MAX(rsqmax,rsqclose[m]);
  }

  // pack my atoms into buffer as bodytag of owning atom, unwrapped coords

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    domain->unmap(x[i],image[i],unwrap);
    buf[m++] = ubuf(bodytag[i]).d;
    buf[m++] = unwrap[0];
    buf[m++] = unwrap[1];
    buf[m++] = unwrap[2];
  }

  // pass buffer around ring of procs
  // func = update rsqfar for atoms belonging to bodies I own
  // when done, have rsqfar for all atoms in bodies I own

  rsqfar = 0.0;
  frsptr = this;
  comm->ring(m,sizeof(double),buf,3,ring_farthest,NULL);

  // find maxextent of rsqfar across all procs
  // if defined, include molecule->maxextent

  MPI_Allreduce(&rsqfar,&maxextent,1,MPI_DOUBLE,MPI_MAX,world);
  maxextent = sqrt(maxextent);
  if (onemol) maxextent = MAX(maxextent,onemol->maxextent);

  // clean up

  delete hash;
  memory->destroy(buf);
  memory->destroy(bbox);
  memory->destroy(ctr);
  memory->destroy(idclose);
  memory->destroy(rsqclose);
}

/* ----------------------------------------------------------------------
   process rigid body atoms from another proc
   update bounding box for rigid bodies my atoms are part of
------------------------------------------------------------------------- */

void FixRigidSmall::ring_bbox(int n, char *cbuf)
{
  std::map<tagint,int> *hash = frsptr->hash;
  double **bbox = frsptr->bbox;

  double *buf = (double *) cbuf;
  int ndatums = n/4;

  int j,imol;
  double *x;

  int m = 0;
  for (int i = 0; i < ndatums; i++, m += 4) {
    imol = static_cast<int> (buf[m]);
    if (hash->find(imol) != hash->end()) {
      j = hash->find(imol)->second;
      x = &buf[m+1];
      bbox[j][0] = MIN(bbox[j][0],x[0]);
      bbox[j][1] = MAX(bbox[j][1],x[0]);
      bbox[j][2] = MIN(bbox[j][2],x[1]);
      bbox[j][3] = MAX(bbox[j][3],x[1]);
      bbox[j][4] = MIN(bbox[j][4],x[2]);
      bbox[j][5] = MAX(bbox[j][5],x[2]);
    }
  }
}

/* ----------------------------------------------------------------------
   process rigid body atoms from another proc
   update nearest atom to body center for rigid bodies my atoms are part of
------------------------------------------------------------------------- */

void FixRigidSmall::ring_nearest(int n, char *cbuf)
{
  std::map<tagint,int> *hash = frsptr->hash;
  double **ctr = frsptr->ctr;
  tagint *idclose = frsptr->idclose;
  double *rsqclose = frsptr->rsqclose;

  double *buf = (double *) cbuf;
  int ndatums = n/5;

  int j,imol;
  tagint tag;
  double delx,dely,delz,rsq;
  double *x;

  int m = 0;
  for (int i = 0; i < ndatums; i++, m += 5) {
    imol = static_cast<int> (buf[m]);
    if (hash->find(imol) != hash->end()) {
      j = hash->find(imol)->second;
      tag = (tagint) ubuf(buf[m+1]).i;
      x = &buf[m+2];
      delx = x[0] - ctr[j][0];
      dely = x[1] - ctr[j][1];
      delz = x[2] - ctr[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq <= rsqclose[j]) {
        if (rsq == rsqclose[j] && tag > idclose[j]) continue;
        idclose[j] = tag;
        rsqclose[j] = rsq;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   process rigid body atoms from another proc
   update rsqfar = distance from owning atom to other atom
------------------------------------------------------------------------- */

void FixRigidSmall::ring_farthest(int n, char *cbuf)
{
  double **x = frsptr->atom->x;
  imageint *image = frsptr->atom->image;
  int nlocal = frsptr->atom->nlocal;

  double *buf = (double *) cbuf;
  int ndatums = n/4;

  int iowner;
  tagint tag;
  double delx,dely,delz,rsq;
  double *xx;
  double unwrap[3];

  int m = 0;
  for (int i = 0; i < ndatums; i++, m += 4) {
    tag = (tagint) ubuf(buf[m]).i;
    iowner = frsptr->atom->map(tag);
    if (iowner < 0 || iowner >= nlocal) continue;
    frsptr->domain->unmap(x[iowner],image[iowner],unwrap);
    xx = &buf[m+1];
    delx = xx[0] - unwrap[0];
    dely = xx[1] - unwrap[1];
    delz = xx[2] - unwrap[2];
    rsq = delx*delx + dely*dely + delz*delz;
    frsptr->rsqfar = MAX(frsptr->rsqfar,rsq);
  }
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes
   extended flags, mass, center-of-mass
   Cartesian and diagonalized inertia tensor
   read per-body attributes from infile if specified
------------------------------------------------------------------------- */

void FixRigidSmall::setup_bodies_static()
{
  int i,itype,ibody;

  // extended = 1 if any particle in a rigid body is finite size
  //              or has a dipole moment

  extended = orientflag = dorientflag = 0;

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double **mu = atom->mu;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (atom->radius_flag || atom->ellipsoid_flag || atom->line_flag ||
      atom->tri_flag || atom->mu_flag) {
    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (bodytag[i] == 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      if (ellipsoid && ellipsoid[i] >= 0) flag = 1;
      if (line && line[i] >= 0) flag = 1;
      if (tri && tri[i] >= 0) flag = 1;
      if (mu && mu[i][3] > 0.0) flag = 1;
    }

    MPI_Allreduce(&flag,&extended,1,MPI_INT,MPI_MAX,world);
  }

  // extended = 1 if using molecule template with finite-size particles

  if (onemol && onemol->radiusflag) extended = 1;

  // grow extended arrays and set extended flags for each particle
  // orientflag = 4 if any particle stores ellipsoid or tri orientation
  // orientflag = 1 if any particle stores line orientation
  // dorientflag = 1 if any particle stores dipole orientation

  if (extended) {
    if (atom->ellipsoid_flag) orientflag = 4;
    if (atom->line_flag) orientflag = 1;
    if (atom->tri_flag) orientflag = 4;
    if (atom->mu_flag) dorientflag = 1;
    grow_arrays(atom->nmax);

    for (i = 0; i < nlocal; i++) {
      eflags[i] = 0;
      if (bodytag[i] == 0) continue;

      // set to POINT or SPHERE or ELLIPSOID or LINE

      if (radius && radius[i] > 0.0) {
        eflags[i] |= SPHERE;
        eflags[i] |= OMEGA;
        eflags[i] |= TORQUE;
      } else if (ellipsoid && ellipsoid[i] >= 0) {
        eflags[i] |= ELLIPSOID;
        eflags[i] |= ANGMOM;
        eflags[i] |= TORQUE;
      } else if (line && line[i] >= 0) {
        eflags[i] |= LINE;
        eflags[i] |= OMEGA;
        eflags[i] |= TORQUE;
      } else if (tri && tri[i] >= 0) {
        eflags[i] |= TRIANGLE;
        eflags[i] |= ANGMOM;
        eflags[i] |= TORQUE;
      } else eflags[i] |= POINT;

      // set DIPOLE if atom->mu and mu[3] > 0.0

      if (atom->mu_flag && mu[i][3] > 0.0)
        eflags[i] |= DIPOLE;
    }
  }

  // acquire ghost bodies via forward comm
  // set atom2body for ghost atoms via forward comm
  // set atom2body for other owned atoms via reset_atom2body()

  nghost_body = 0;
  commflag = FULL_BODY;
  comm->forward_comm_variable_fix(this);
  reset_atom2body();

  // compute mass & center-of-mass of each rigid body

  double **x = atom->x;
  imageint *image = atom->image;

  double *xcm;

  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    xcm = body[ibody].xcm;
    xcm[0] = xcm[1] = xcm[2] = 0.0;
    body[ibody].mass = 0.0;
  }

  double unwrap[3];
  double massone;

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    xcm[0] += unwrap[0] * massone;
    xcm[1] += unwrap[1] * massone;
    xcm[2] += unwrap[2] * massone;
    b->mass += massone;
  }

  // reverse communicate xcm, mass of all bodies

  commflag = XCM_MASS;
  comm->reverse_comm_variable_fix(this);

  for (ibody = 0; ibody < nlocal_body; ibody++) {
    xcm = body[ibody].xcm;
    xcm[0] /= body[ibody].mass;
    xcm[1] /= body[ibody].mass;
    xcm[2] /= body[ibody].mass;
  }

  // overwrite masstotal and center-of-mass with file values
  // inbody[i] = 0/1 if Ith rigid body is initialized by file

  int *inbody;
  if (infile) {
    memory->create(inbody,nlocal_body,"rigid/small:inbody");
    for (ibody = 0; ibody < nlocal_body; ibody++) inbody[ibody] = 0;
    readfile(0,NULL,inbody);
  }

  // set image flags for each rigid body to default values
  // then remap the xcm of each body back into simulation box if needed

  for (ibody = 0; ibody < nlocal_body; ibody++)
    body[ibody].image = ((imageint) IMGMAX << IMG2BITS) | 
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  pre_neighbor();

  // compute 6 moments of inertia of each body in Cartesian reference frame
  // dx,dy,dz = coords relative to center-of-mass
  // symmetric 3x3 inertia tensor stored in Voigt notation as 6-vector

  memory->create(itensor,nlocal_body+nghost_body,6,"rigid/small:itensor");
  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++)
    for (i = 0; i < 6; i++) itensor[ibody][i] = 0.0;

  double dx,dy,dz,rad;
  double *inertia;

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    dx = unwrap[0] - xcm[0];
    dy = unwrap[1] - xcm[1];
    dz = unwrap[2] - xcm[2];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    inertia = itensor[atom2body[i]];
    inertia[0] += massone * (dy*dy + dz*dz);
    inertia[1] += massone * (dx*dx + dz*dz);
    inertia[2] += massone * (dx*dx + dy*dy);
    inertia[3] -= massone * dy*dz;
    inertia[4] -= massone * dx*dz;
    inertia[5] -= massone * dx*dy;
  }

  // extended particles may contribute extra terms to moments of inertia

  if (extended) {
    double ivec[6];
    double *shape,*quatatom,*inertiaatom;
    double length,theta;

    for (i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      inertia = itensor[atom2body[i]];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      if (eflags[i] & SPHERE) {
        inertia[0] += SINERTIA*massone * radius[i]*radius[i];
        inertia[1] += SINERTIA*massone * radius[i]*radius[i];
        inertia[2] += SINERTIA*massone * radius[i]*radius[i];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::inertia_ellipsoid(shape,quatatom,massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      } else if (eflags[i] & LINE) {
        length = lbonus[line[i]].length;
        theta = lbonus[line[i]].theta;
        MathExtra::inertia_line(length,theta,massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::inertia_triangle(inertiaatom,quatatom,massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      }
    }
  }

  // reverse communicate inertia tensor of all bodies

  commflag = ITENSOR;
  comm->reverse_comm_variable_fix(this);

  // overwrite Cartesian inertia tensor with file values

  if (infile) readfile(1,itensor,inbody);

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy_space = 3 evectors = principal axes of rigid body

  int ierror;
  double cross[3];
  double tensor[3][3],evectors[3][3];
  double *ex,*ey,*ez;

  for (ibody = 0; ibody < nlocal_body; ibody++) {
    tensor[0][0] = itensor[ibody][0];
    tensor[1][1] = itensor[ibody][1];
    tensor[2][2] = itensor[ibody][2];
    tensor[1][2] = tensor[2][1] = itensor[ibody][3];
    tensor[0][2] = tensor[2][0] = itensor[ibody][4];
    tensor[0][1] = tensor[1][0] = itensor[ibody][5];

    inertia = body[ibody].inertia;
    ierror = MathExtra::jacobi(tensor,inertia,evectors);
    if (ierror) error->all(FLERR,
                           "Insufficient Jacobi rotations for rigid body");

    ex = body[ibody].ex_space;
    ex[0] = evectors[0][0];
    ex[1] = evectors[1][0];
    ex[2] = evectors[2][0];
    ey = body[ibody].ey_space;
    ey[0] = evectors[0][1];
    ey[1] = evectors[1][1];
    ey[2] = evectors[2][1];
    ez = body[ibody].ez_space;
    ez[0] = evectors[0][2];
    ez[1] = evectors[1][2];
    ez[2] = evectors[2][2];

    // if any principal moment < scaled EPSILON, set to 0.0

    double max;
    max = MAX(inertia[0],inertia[1]);
    max = MAX(max,inertia[2]);

    if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
    if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
    if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed

    MathExtra::cross3(ex,ey,cross);
    if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);

    // create initial quaternion

    MathExtra::exyz_to_q(ex,ey,ez,body[ibody].quat);
  }

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  comm->forward_comm_variable_fix(this);

  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  double qc[4],delta[3];
  double *quatatom;
  double theta_body;

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) {
      displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
      continue;
    }

    Body *b = &body[atom2body[i]];

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    delta[0] = unwrap[0] - xcm[0];
    delta[1] = unwrap[1] - xcm[1];
    delta[2] = unwrap[2] - xcm[2];
    MathExtra::transpose_matvec(b->ex_space,b->ey_space,b->ez_space,
                                delta,displace[i]);

    if (extended) {
      if (eflags[i] & ELLIPSOID) {
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::qconjugate(b->quat,qc);
        MathExtra::quatquat(qc,quatatom,orient[i]);
        MathExtra::qnormalize(orient[i]);
      } else if (eflags[i] & LINE) {
        if (b->quat[3] >= 0.0) theta_body = 2.0*acos(b->quat[0]);
        else theta_body = -2.0*acos(b->quat[0]);
        orient[i][0] = lbonus[line[i]].theta - theta_body;
        while (orient[i][0] <= MINUSPI) orient[i][0] += TWOPI;
        while (orient[i][0] > MY_PI) orient[i][0] -= TWOPI;
        if (orientflag == 4) orient[i][1] = orient[i][2] = orient[i][3] = 0.0;
      } else if (eflags[i] & TRIANGLE) {
        quatatom = tbonus[tri[i]].quat;
        MathExtra::qconjugate(b->quat,qc);
        MathExtra::quatquat(qc,quatatom,orient[i]);
        MathExtra::qnormalize(orient[i]);
      } else if (orientflag == 4) {
        orient[i][0] = orient[i][1] = orient[i][2] = orient[i][3] = 0.0;
      } else if (orientflag == 1)
        orient[i][0] = 0.0;

      if (eflags[i] & DIPOLE) {
        MathExtra::transpose_matvec(b->ex_space,b->ey_space,b->ez_space,
                                    mu[i],dorient[i]);
        MathExtra::snormalize3(mu[i][3],dorient[i],dorient[i]);
      } else if (dorientflag)
        dorient[i][0] = dorient[i][1] = dorient[i][2] = 0.0;
    }
  }

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles may contribute extra terms to moments of inertia

  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++)
    for (i = 0; i < 6; i++) itensor[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    inertia = itensor[atom2body[i]];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    inertia[0] += massone *
      (displace[i][1]*displace[i][1] + displace[i][2]*displace[i][2]);
    inertia[1] += massone *
      (displace[i][0]*displace[i][0] + displace[i][2]*displace[i][2]);
    inertia[2] += massone *
      (displace[i][0]*displace[i][0] + displace[i][1]*displace[i][1]);
    inertia[3] -= massone * displace[i][1]*displace[i][2];
    inertia[4] -= massone * displace[i][0]*displace[i][2];
    inertia[5] -= massone * displace[i][0]*displace[i][1];
  }

  if (extended) {
    double ivec[6];
    double *shape,*inertiaatom;
    double length;

    for (i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      inertia = itensor[atom2body[i]];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      if (eflags[i] & SPHERE) {
        inertia[0] += SINERTIA*massone * radius[i]*radius[i];
        inertia[1] += SINERTIA*massone * radius[i]*radius[i];
        inertia[2] += SINERTIA*massone * radius[i]*radius[i];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        MathExtra::inertia_ellipsoid(shape,orient[i],massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      } else if (eflags[i] & LINE) {
        length = lbonus[line[i]].length;
        MathExtra::inertia_line(length,orient[i][0],massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        MathExtra::inertia_triangle(inertiaatom,orient[i],massone,ivec);
        inertia[0] += ivec[0];
        inertia[1] += ivec[1];
        inertia[2] += ivec[2];
        inertia[3] += ivec[3];
        inertia[4] += ivec[4];
        inertia[5] += ivec[5];
      }
    }
  }

  // reverse communicate inertia tensor of all bodies

  commflag = ITENSOR;
  comm->reverse_comm_variable_fix(this);

  // error check that re-computed momemts of inertia match diagonalized ones
  // do not do test for bodies with params read from infile

  double *inew;

  double norm;
  for (ibody = 0; ibody < nlocal_body; ibody++) {
    if (infile && inbody[ibody]) continue;
    inertia = body[ibody].inertia;

    if (inertia[0] == 0.0) {
      if (fabs(itensor[ibody][0]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((itensor[ibody][0]-inertia[0])/inertia[0]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    if (inertia[1] == 0.0) {
      if (fabs(itensor[ibody][1]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((itensor[ibody][1]-inertia[1])/inertia[1]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    if (inertia[2] == 0.0) {
      if (fabs(itensor[ibody][2]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((itensor[ibody][2]-inertia[2])/inertia[2]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    norm = (inertia[0] + inertia[1] + inertia[2]) / 3.0;
    if (fabs(itensor[ibody][3]/norm) > TOLERANCE ||
        fabs(itensor[ibody][4]/norm) > TOLERANCE ||
        fabs(itensor[ibody][5]/norm) > TOLERANCE)
      error->all(FLERR,"Fix rigid: Bad principal moments");
  }

  // clean up

  memory->destroy(itensor);
  if (infile) memory->destroy(inbody);
}

/* ----------------------------------------------------------------------
   one-time initialization of dynamic rigid body attributes
   Vcm and angmom, computed explicitly from constituent particles
   even if wrong for overlapping particles, is OK,
     since is just setting initial time=0 Vcm and angmom of the body
     which can be estimated value
------------------------------------------------------------------------- */

void FixRigidSmall::setup_bodies_dynamic()
{
  int i,n,ibody;
  double massone,radone;

  // sum vcm, angmom across all rigid bodies
  // vcm = velocity of COM
  // angmom = angular momentum around COM

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *xcm,*vcm,*acm;
  double dx,dy,dz;
  double unwrap[3];

  for (ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    vcm = body[ibody].vcm;
    vcm[0] = vcm[1] = vcm[2] = 0.0;
    acm = body[ibody].angmom;
    acm[0] = acm[1] = acm[2] = 0.0;
  }

  for (i = 0; i < nlocal; i++) {
    if (atom2body[i] < 0) continue;
    Body *b = &body[atom2body[i]];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    vcm = b->vcm;
    vcm[0] += v[i][0] * massone;
    vcm[1] += v[i][1] * massone;
    vcm[2] += v[i][2] * massone;

    domain->unmap(x[i],image[i],unwrap);
    xcm = b->xcm;
    dx = unwrap[0] - xcm[0];
    dy = unwrap[1] - xcm[1];
    dz = unwrap[2] - xcm[2];

    acm = b->angmom;
    acm[0] += dy * massone*v[i][2] - dz * massone*v[i][1];
    acm[1] += dz * massone*v[i][0] - dx * massone*v[i][2];
    acm[2] += dx * massone*v[i][1] - dy * massone*v[i][0];
  }

  // extended particles add their rotation to angmom of body

  if (extended) {
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    double **omega = atom->omega;
    double **angmom = atom->angmom;
    double *radius = atom->radius;
    int *line = atom->line;

    for (i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body *b = &body[atom2body[i]];

      if (eflags[i] & OMEGA) {
        if (eflags[i] & SPHERE) {
          radone = radius[i];
          acm = b->angmom;
          acm[0] += SINERTIA*rmass[i] * radone*radone * omega[i][0];
          acm[1] += SINERTIA*rmass[i] * radone*radone * omega[i][1];
          acm[2] += SINERTIA*rmass[i] * radone*radone * omega[i][2];
        } else if (eflags[i] & LINE) {
          radone = lbonus[line[i]].length;
          b->angmom[2] += LINERTIA*rmass[i] * radone*radone * omega[i][2];
        }
      }
      if (eflags[i] & ANGMOM) {
        acm = b->angmom;
        acm[0] += angmom[i][0];
        acm[1] += angmom[i][1];
        acm[2] += angmom[i][2];
      }
    }
  }

  // reverse communicate vcm, angmom of all bodies

  commflag = VCM_ANGMOM;
  comm->reverse_comm_variable_fix(this);

  // normalize velocity of COM

  for (ibody = 0; ibody < nlocal_body; ibody++) {
    vcm = body[ibody].vcm;
    vcm[0] /= body[ibody].mass;
    vcm[1] /= body[ibody].mass;
    vcm[2] /= body[ibody].mass;
  }
}

/* ----------------------------------------------------------------------
   read per rigid body info from user-provided file
   which = 0 to read total mass and center-of-mass
   which = 1 to read 6 moments of inertia, store in array
   flag inbody = 0 for local bodies whose info is read from file
   nlines = # of lines of rigid body info
   one line = rigid-ID mass xcm ycm zcm ixx iyy izz ixy ixz iyz
   and rigid-ID = mol-ID for fix rigid/small
------------------------------------------------------------------------- */

void FixRigidSmall::readfile(int which, double **array, int *inbody)
{
  int i,j,m,nchunk,eofflag,nlines;
  tagint id;
  FILE *fp;
  char *eof,*start,*next,*buf;
  char line[MAXLINE];
  
  // create local hash with key/value pairs
  // key = mol ID of bodies my atoms own
  // value = index into local body array

  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  hash = new std::map<tagint,int>();
  for (int i = 0; i < nlocal; i++)
    if (bodyown[i] >= 0) (*hash)[atom->molecule[i]] = bodyown[i];

  // open file and read header

  if (me == 0) {
    fp = fopen(infile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix rigid/small infile %s",infile);
      error->one(FLERR,str);
    }

    while (1) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) 
        error->one(FLERR,"Unexpected end of fix rigid/small file");
      start = &line[strspn(line," \t\n\v\f\r")];
      if (*start != '\0' && *start != '#') break;
    }

    sscanf(line,"%d",&nlines);
  }

  MPI_Bcast(&nlines,1,MPI_INT,0,world);
  if (nlines == 0) error->all(FLERR,"Fix rigid file has no lines");

  char *buffer = new char[CHUNK*MAXLINE];
  char **values = new char*[ATTRIBUTE_PERBODY];

  int nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines-nread,CHUNK);
    eofflag = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eofflag) error->all(FLERR,"Unexpected end of fix rigid/small file");

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = atom->count_words(buf);
    *next = '\n';

    if (nwords != ATTRIBUTE_PERBODY)
      error->all(FLERR,"Incorrect rigid body format in fix rigid/small file");
    
    // loop over lines of rigid body attributes
    // tokenize the line into values
    // id = rigid body ID = mol-ID
    // for which = 0, store mass/com in vec/array
    // for which = 1, store interia tensor array, invert 3,4,5 values to Voigt

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      
      values[0] = strtok(buf," \t\n\r\f");
      for (j = 1; j < nwords; j++)
        values[j] = strtok(NULL," \t\n\r\f");

      id = ATOTAGINT(values[0]);
      if (id <= 0 || id > maxmol) 
        error->all(FLERR,"Invalid rigid body ID in fix rigid/small file");
      if (hash->find(id) == hash->end()) {
        buf = next + 1;
        continue;
      }
      m = (*hash)[id];
      inbody[m] = 1;

      if (which == 0) {
        body[m].mass = atof(values[1]);
        body[m].xcm[0] = atof(values[2]);
        body[m].xcm[1] = atof(values[3]);
        body[m].xcm[2] = atof(values[4]);
      } else {
        array[m][0] = atof(values[5]);
        array[m][1] = atof(values[6]);
        array[m][2] = atof(values[7]);
        array[m][3] = atof(values[10]);
        array[m][4] = atof(values[9]);
        array[m][5] = atof(values[8]);
      }

      buf = next + 1;
    }
    
    nread += nchunk;
  }

  if (me == 0) fclose(fp);

  delete [] buffer;
  delete [] values;
  delete hash;
}

/* ----------------------------------------------------------------------
   write out restart info for mass, COM, inertia tensor to file
   identical format to infile option, so info can be read in when restarting
   each proc contributes info for rigid bodies it owns
------------------------------------------------------------------------- */

void FixRigidSmall::write_restart_file(char *file)
{
  FILE *fp;

  // do not write file if bodies have not yet been intialized

  if (firstflag) return;

  // proc 0 opens file and writes header

  if (me == 0) {
    char outfile[128];
    sprintf(outfile,"%s.rigid",file);
    fp = fopen(outfile,"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix rigid restart file %s",outfile);
      error->one(FLERR,str);
    }

    fprintf(fp,"# fix rigid mass, COM, inertia tensor info for "
            "%d bodies on timestep " BIGINT_FORMAT "\n\n",
            nbody,update->ntimestep);
    fprintf(fp,"%d\n",nbody);
  }

  // communication buffer for all my rigid body info
  // max_size = largest buffer needed by any proc

  int ncol = 11;
  int sendrow = nlocal_body;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"rigid/small:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"rigid/small:buf");

  // pack my rigid body info into buf
  // compute I tensor against xyz axes from diagonalized I and current quat
  // Ispace = P Idiag P_transpose
  // P is stored column-wise in exyz_space

  double p[3][3],pdiag[3][3],ispace[3][3];

  for (int i = 0; i < nlocal_body; i++) {
    MathExtra::col2mat(body[i].ex_space,body[i].ey_space,body[i].ez_space,p);
    MathExtra::times3_diag(p,body[i].inertia,pdiag);
    MathExtra::times3_transpose(pdiag,p,ispace);

    buf[i][0] = atom->molecule[body[i].ilocal];
    buf[i][1] = body[i].mass;
    buf[i][2] = body[i].xcm[0];
    buf[i][3] = body[i].xcm[1];
    buf[i][4] = body[i].xcm[2];
    buf[i][5] = ispace[0][0];
    buf[i][6] = ispace[1][1];
    buf[i][7] = ispace[2][2];
    buf[i][8] = ispace[0][1];
    buf[i][9] = ispace[0][2];
    buf[i][10] = ispace[1][2];
  }

  // write one chunk of rigid body info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      for (int i = 0; i < recvrow; i++)
        fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e "
                "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
                static_cast<int> (buf[i][0]),buf[i][1],
                buf[i][2],buf[i][3],buf[i][4],
                buf[i][5],buf[i][6],buf[i][7],buf[i][8],buf[i][9],buf[i][10]);
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  // clean up and close file

  memory->destroy(buf);
  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixRigidSmall::grow_arrays(int nmax)
{
  memory->grow(bodyown,nmax,"rigid/small:bodyown");
  memory->grow(bodytag,nmax,"rigid/small:bodytag");
  memory->grow(atom2body,nmax,"rigid/small:atom2body");
  memory->grow(displace,nmax,3,"rigid/small:displace");
  if (extended) {
    memory->grow(eflags,nmax,"rigid/small:eflags");
    if (orientflag) memory->grow(orient,nmax,orientflag,"rigid/small:orient");
    if (dorientflag) memory->grow(dorient,nmax,3,"rigid/small:dorient");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixRigidSmall::copy_arrays(int i, int j, int delflag)
{
  bodytag[j] = bodytag[i];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];

  if (extended) {
    eflags[j] = eflags[i];
    for (int k = 0; k < orientflag; k++)
      orient[j][k] = orient[i][k];
    if (dorientflag) {
      dorient[j][0] = dorient[i][0];
      dorient[j][1] = dorient[i][1];
      dorient[j][2] = dorient[i][2];
    }
  }

  // if deleting atom J via delflag and J owns a body, then delete it

  if (delflag && bodyown[j] >= 0) {
    bodyown[body[nlocal_body-1].ilocal] = bodyown[j];
    memcpy(&body[bodyown[j]],&body[nlocal_body-1],sizeof(Body));
    nlocal_body--;
  }

  // if atom I owns a body, reset I's body.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's body is already deleted

  if (bodyown[i] >= 0 && i != j) body[bodyown[i]].ilocal = j;
  bodyown[j] = bodyown[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixRigidSmall::set_arrays(int i)
{
  bodyown[i] = -1;
  bodytag[i] = 0;
  atom2body[i] = -1;
  displace[i][0] = 0.0;
  displace[i][1] = 0.0;
  displace[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   initialize a molecule inserted by another fix, e.g. deposit or pour
   called when molecule is created
   nlocalprev = # of atoms on this proc before molecule inserted
   tagprev = atom ID previous to new atoms in the molecule
   xgeom = geometric center of new molecule
   vcm = COM velocity of new molecule
   quat = rotation of new molecule (around geometric center)
          relative to template in Molecule class
------------------------------------------------------------------------- */

void FixRigidSmall::set_molecule(int nlocalprev, tagint tagprev, 
                                 double *xgeom, double *vcm, double *quat)
{
  int m;
  double ctr2com[3],ctr2com_rotate[3];
  double rotmat[3][3];

  int nlocal = atom->nlocal;
  if (nlocalprev == nlocal) return;

  double **x = atom->x;
  tagint *tag = atom->tag;

  for (int i = nlocalprev; i < nlocal; i++) {
    bodytag[i] = tagprev + onemol->comatom;
    if (tag[i]-tagprev == onemol->comatom) bodyown[i] = nlocal_body;

    m = tag[i] - tagprev-1;
    displace[i][0] = onemol->dxbody[m][0];
    displace[i][1] = onemol->dxbody[m][1];
    displace[i][2] = onemol->dxbody[m][2];

    eflags[i] = 0;
    if (onemol->radiusflag) {
      eflags[i] |= SPHERE;
      eflags[i] |= OMEGA;
      eflags[i] |= TORQUE;
    }

    if (bodyown[i] >= 0) {
      if (nlocal_body == nmax_body) grow_body();
      Body *b = &body[nlocal_body];
      b->mass = onemol->masstotal;

      // new COM = Q (onemol->xcm - onemol->center) + xgeom
      // Q = rotation matrix associated with quat

      MathExtra::quat_to_mat(quat,rotmat);
      MathExtra::sub3(onemol->com,onemol->center,ctr2com);
      MathExtra::matvec(rotmat,ctr2com,ctr2com_rotate);
      MathExtra::add3(ctr2com_rotate,xgeom,b->xcm);

      b->vcm[0] = vcm[0];
      b->vcm[1] = vcm[1];
      b->vcm[2] = vcm[2];
      b->inertia[0] = onemol->inertia[0];
      b->inertia[1] = onemol->inertia[1];
      b->inertia[2] = onemol->inertia[2];

      // final quat is product of insertion quat and original quat
      // true even if insertion rotation was not around COM

      MathExtra::quatquat(quat,onemol->quat,b->quat);
      MathExtra::q_to_exyz(b->quat,b->ex_space,b->ey_space,b->ez_space);

      b->angmom[0] = b->angmom[1] = b->angmom[2] = 0.0;
      b->omega[0] = b->omega[1] = b->omega[2] = 0.0;
      b->image = ((imageint) IMGMAX << IMG2BITS) |
        ((imageint) IMGMAX << IMGBITS) | IMGMAX;
      b->ilocal = i;
      nlocal_body++;
    }
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixRigidSmall::pack_exchange(int i, double *buf)
{
  buf[0] = ubuf(bodytag[i]).d;
  buf[1] = displace[i][0];
  buf[2] = displace[i][1];
  buf[3] = displace[i][2];

  // extended attribute info

  int m = 4;
  if (extended) {
    buf[m++] = eflags[i];
    for (int j = 0; j < orientflag; j++)
      buf[m++] = orient[i][j];
    if (dorientflag) {
      buf[m++] = dorient[i][0];
      buf[m++] = dorient[i][1];
      buf[m++] = dorient[i][2];
    }
  }

  // atom not in a rigid body

  if (!bodytag[i]) return m;

  // atom does not own its rigid body

  if (bodyown[i] < 0) {
    buf[m++] = 0;
    return m;
  }

  // body info for atom that owns a rigid body

  buf[m++] = 1;
  memcpy(&buf[m],&body[bodyown[i]],sizeof(Body));
  m += bodysize;
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixRigidSmall::unpack_exchange(int nlocal, double *buf)
{
  bodytag[nlocal] = (tagint) ubuf(buf[0]).i;
  displace[nlocal][0] = buf[1];
  displace[nlocal][1] = buf[2];
  displace[nlocal][2] = buf[3];

  // extended attribute info

  int m = 4;
  if (extended) {
    eflags[nlocal] = static_cast<int> (buf[m++]);
    for (int j = 0; j < orientflag; j++)
      orient[nlocal][j] = buf[m++];
    if (dorientflag) {
      dorient[nlocal][0] = buf[m++];
      dorient[nlocal][1] = buf[m++];
      dorient[nlocal][2] = buf[m++];
    }
  }

  // atom not in a rigid body

  if (!bodytag[nlocal]) {
    bodyown[nlocal] = -1;
    return m;
  }

  // atom does not own its rigid body

  bodyown[nlocal] = static_cast<int> (buf[m++]);
  if (bodyown[nlocal] == 0) {
    bodyown[nlocal] = -1;
    return m;
  }

  // body info for atom that owns a rigid body

  if (nlocal_body == nmax_body) grow_body();
  memcpy(&body[nlocal_body],&buf[m],sizeof(Body));
  m += bodysize;
  body[nlocal_body].ilocal = nlocal;
  bodyown[nlocal] = nlocal_body++;

  return m;
}

/* ----------------------------------------------------------------------
   only pack body info if own or ghost atom owns the body
   for FULL_BODY, send 0/1 flag with every atom
------------------------------------------------------------------------- */

int FixRigidSmall::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j;
  double *xcm,*vcm,*quat,*omega,*ex_space,*ey_space,*ez_space;

  int nlocal = atom->nlocal;

  int m = 0;

  if (commflag == INITIAL) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      xcm = body[bodyown[j]].xcm;
      buf[m++] = xcm[0];
      buf[m++] = xcm[1];
      buf[m++] = xcm[2];
      vcm = body[bodyown[j]].vcm;
      buf[m++] = vcm[0];
      buf[m++] = vcm[1];
      buf[m++] = vcm[2];
      quat = body[bodyown[j]].quat;
      buf[m++] = quat[0];
      buf[m++] = quat[1];
      buf[m++] = quat[2];
      buf[m++] = quat[3];
      omega = body[bodyown[j]].omega;
      buf[m++] = omega[0];
      buf[m++] = omega[1];
      buf[m++] = omega[2];
      ex_space = body[bodyown[j]].ex_space;
      buf[m++] = ex_space[0];
      buf[m++] = ex_space[1];
      buf[m++] = ex_space[2];
      ey_space = body[bodyown[j]].ey_space;
      buf[m++] = ey_space[0];
      buf[m++] = ey_space[1];
      buf[m++] = ey_space[2];
      ez_space = body[bodyown[j]].ez_space;
      buf[m++] = ez_space[0];
      buf[m++] = ez_space[1];
      buf[m++] = ez_space[2];
    }

  } else if (commflag == FINAL) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      vcm = body[bodyown[j]].vcm;
      buf[m++] = vcm[0];
      buf[m++] = vcm[1];
      buf[m++] = vcm[2];
      omega = body[bodyown[j]].omega;
      buf[m++] = omega[0];
      buf[m++] = omega[1];
      buf[m++] = omega[2];
    }

  } else if (commflag == FULL_BODY) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) buf[m++] = 0;
      else {
        buf[m++] = 1;
        memcpy(&buf[m],&body[bodyown[j]],sizeof(Body));
        m += bodysize;
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   only ghost atoms are looped over
   for FULL_BODY, store a new ghost body if this atom owns it
   for other commflag values, only unpack body info if atom owns it
------------------------------------------------------------------------- */

void FixRigidSmall::unpack_comm(int n, int first, double *buf)
{
  int i,j,last,flag;
  double *xcm,*vcm,*quat,*omega,*ex_space,*ey_space,*ez_space;

  int m = 0;
  last = first + n;

  if (commflag == INITIAL) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      xcm = body[bodyown[i]].xcm;
      xcm[0] = buf[m++];
      xcm[1] = buf[m++];
      xcm[2] = buf[m++];
      vcm = body[bodyown[i]].vcm;
      vcm[0] = buf[m++];
      vcm[1] = buf[m++];
      vcm[2] = buf[m++];
      quat = body[bodyown[i]].quat;
      quat[0] = buf[m++];
      quat[1] = buf[m++];
      quat[2] = buf[m++];
      quat[3] = buf[m++];
      omega = body[bodyown[i]].omega;
      omega[0] = buf[m++];
      omega[1] = buf[m++];
      omega[2] = buf[m++];
      ex_space = body[bodyown[i]].ex_space;
      ex_space[0] = buf[m++];
      ex_space[1] = buf[m++];
      ex_space[2] = buf[m++];
      ey_space = body[bodyown[i]].ey_space;
      ey_space[0] = buf[m++];
      ey_space[1] = buf[m++];
      ey_space[2] = buf[m++];
      ez_space = body[bodyown[i]].ez_space;
      ez_space[0] = buf[m++];
      ez_space[1] = buf[m++];
      ez_space[2] = buf[m++];
    }

  } else if (commflag == FINAL) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      vcm = body[bodyown[i]].vcm;
      vcm[0] = buf[m++];
      vcm[1] = buf[m++];
      vcm[2] = buf[m++];
      omega = body[bodyown[i]].omega;
      omega[0] = buf[m++];
      omega[1] = buf[m++];
      omega[2] = buf[m++];
    }

  } else if (commflag == FULL_BODY) {
    for (i = first; i < last; i++) {
      bodyown[i] = static_cast<int> (buf[m++]);
      if (bodyown[i] == 0) bodyown[i] = -1;
      else {
        j = nlocal_body + nghost_body;
        if (j == nmax_body) grow_body();
        memcpy(&body[j],&buf[m],sizeof(Body));
        m += bodysize;
        body[j].ilocal = i;
        bodyown[i] = j;
        nghost_body++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   only ghost atoms are looped over
   only pack body info if atom owns it
------------------------------------------------------------------------- */

int FixRigidSmall::pack_reverse_comm(int n, int first, double *buf)
{
  int i,j,m,last;
  double *fcm,*torque,*vcm,*angmom,*xcm;

  m = 0;
  last = first + n;

  if (commflag == FORCE_TORQUE) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      fcm = body[bodyown[i]].fcm;
      buf[m++] = fcm[0];
      buf[m++] = fcm[1];
      buf[m++] = fcm[2];
      torque = body[bodyown[i]].torque;
      buf[m++] = torque[0];
      buf[m++] = torque[1];
      buf[m++] = torque[2];
    }

  } else if (commflag == VCM_ANGMOM) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      vcm = body[bodyown[i]].vcm;
      buf[m++] = vcm[0];
      buf[m++] = vcm[1];
      buf[m++] = vcm[2];
      angmom = body[bodyown[i]].angmom;
      buf[m++] = angmom[0];
      buf[m++] = angmom[1];
      buf[m++] = angmom[2];
    }

  } else if (commflag == XCM_MASS) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      xcm = body[bodyown[i]].xcm;
      buf[m++] = xcm[0];
      buf[m++] = xcm[1];
      buf[m++] = xcm[2];
      buf[m++] = body[bodyown[i]].mass;
    }

  } else if (commflag == ITENSOR) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      j = bodyown[i];
      buf[m++] = itensor[j][0];
      buf[m++] = itensor[j][1];
      buf[m++] = itensor[j][2];
      buf[m++] = itensor[j][3];
      buf[m++] = itensor[j][4];
      buf[m++] = itensor[j][5];
    }

  } else if (commflag == DOF) {
    for (i = first; i < last; i++) {
      if (bodyown[i] < 0) continue;
      j = bodyown[i];
      buf[m++] = counts[j][0];
      buf[m++] = counts[j][1];
      buf[m++] = counts[j][2];
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   only unpack body info if own or ghost atom owns the body
------------------------------------------------------------------------- */

void FixRigidSmall::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k;
  double *fcm,*torque,*vcm,*angmom,*xcm;

  int nlocal = atom->nlocal;

  int m = 0;

  if (commflag == FORCE_TORQUE) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      fcm = body[bodyown[j]].fcm;
      fcm[0] += buf[m++];
      fcm[1] += buf[m++];
      fcm[2] += buf[m++];
      torque = body[bodyown[j]].torque;
      torque[0] += buf[m++];
      torque[1] += buf[m++];
      torque[2] += buf[m++];
    }

  } else if (commflag == VCM_ANGMOM) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      vcm = body[bodyown[j]].vcm;
      vcm[0] += buf[m++];
      vcm[1] += buf[m++];
      vcm[2] += buf[m++];
      angmom = body[bodyown[j]].angmom;
      angmom[0] += buf[m++];
      angmom[1] += buf[m++];
      angmom[2] += buf[m++];
    }

  } else if (commflag == XCM_MASS) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      xcm = body[bodyown[j]].xcm;
      xcm[0] += buf[m++];
      xcm[1] += buf[m++];
      xcm[2] += buf[m++];
      body[bodyown[j]].mass += buf[m++];
    }

  } else if (commflag == ITENSOR) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      k = bodyown[j];
      itensor[k][0] += buf[m++];
      itensor[k][1] += buf[m++];
      itensor[k][2] += buf[m++];
      itensor[k][3] += buf[m++];
      itensor[k][4] += buf[m++];
      itensor[k][5] += buf[m++];
    }

  } else if (commflag == DOF) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (bodyown[j] < 0) continue;
      k = bodyown[j];
      counts[k][0] += static_cast<int> (buf[m++]);
      counts[k][1] += static_cast<int> (buf[m++]);
      counts[k][2] += static_cast<int> (buf[m++]);
    }
  }
}

/* ----------------------------------------------------------------------
   grow body data structure
------------------------------------------------------------------------- */

void FixRigidSmall::grow_body()
{
  nmax_body += DELTA_BODY;
  body = (Body *) memory->srealloc(body,nmax_body*sizeof(Body),
                                   "rigid/small:body");
}

/* ----------------------------------------------------------------------
   reset atom2body for all owned atoms
   do this via bodyown of atom that owns the body the owned atom is in
   atom2body values can point to original body or any image of the body
------------------------------------------------------------------------- */

void FixRigidSmall::reset_atom2body()
{
  int iowner;

  // iowner = index of atom that owns the body that atom I is in

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    atom2body[i] = -1;
    if (bodytag[i]) {
      iowner = atom->map(bodytag[i]);
      if (iowner == -1) {
        char str[128];
        sprintf(str,
                "Rigid body atoms " TAGINT_FORMAT " " TAGINT_FORMAT 
                " missing on proc %d at step " BIGINT_FORMAT,
                atom->tag[i],bodytag[i],comm->me,update->ntimestep);
        error->one(FLERR,str);
        
      }
      atom2body[i] = bodyown[iowner];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidSmall::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
}

/* ----------------------------------------------------------------------
   zero linear momentum of each rigid body
   set Vcm to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixRigidSmall::zero_momentum()
{
  double *vcm;
  for (int ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    vcm = body[ibody].vcm;
    vcm[0] = vcm[1] = vcm[2] = 0.0;
  }

  // forward communicate of vcm to all ghost copies

  commflag = FINAL;
  comm->forward_comm_variable_fix(this);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_v();
}

/* ----------------------------------------------------------------------
   zero angular momentum of each rigid body
   set angmom/omega to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixRigidSmall::zero_rotation()
{
  double *angmom,*omega;
  for (int ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    angmom = body[ibody].angmom;
    angmom[0] = angmom[1] = angmom[2] = 0.0;
    omega = body[ibody].omega;
    omega[0] = omega[1] = omega[2] = 0.0;
  }

  // forward communicate of omega to all ghost copies

  commflag = FINAL;
  comm->forward_comm_variable_fix(this);

  // set velocity of atoms in rigid bodues

  evflag = 0;
  set_v();
}

/* ---------------------------------------------------------------------- */

void *FixRigidSmall::extract(const char *str, int &dim)
{
  if (strcmp(str,"body") == 0) {
    dim = 1;
    return atom2body;
  }

  if (strcmp(str,"onemol") == 0) {
    dim = 0;
    return onemol;
  }

  // return vector of rigid body masses, for owned+ghost bodies
  // used by granular pair styles, indexed by atom2body

  if (strcmp(str,"masstotal") == 0) {
    dim = 1;

    if (nmax_mass < nmax_body) {
      memory->destroy(mass_body);
      nmax_mass = nmax_body;
      memory->create(mass_body,nmax_mass,"rigid:mass_body");
    }

    int n = nlocal_body + nghost_body;
    for (int i = 0; i < n; i++)
      mass_body[i] = body[i].mass;

    return mass_body;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
   sum local body results across procs
------------------------------------------------------------------------- */

double FixRigidSmall::extract_ke()
{
  double *vcm;

  double ke = 0.0;
  for (int i = 0; i < nlocal_body; i++) {
    vcm = body[i].vcm;
    ke += body[i].mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);
  }

  double keall;
  MPI_Allreduce(&ke,&keall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*keall;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

double FixRigidSmall::extract_erotational()
{
  double wbody[3],rot[3][3];
  double *inertia;

  double erotate = 0.0;
  for (int i = 0; i < nlocal_body; i++) {

    // for Iw^2 rotational term, need wbody = angular velocity in body frame 
    // not omega = angular velocity in space frame

    inertia = body[i].inertia;
    MathExtra::quat_to_mat(body[i].quat,rot);
    MathExtra::transpose_matvec(rot,body[i].angmom,wbody);
    if (inertia[0] == 0.0) wbody[0] = 0.0;
    else wbody[0] /= inertia[0];
    if (inertia[1] == 0.0) wbody[1] = 0.0;
    else wbody[1] /= inertia[1];
    if (inertia[2] == 0.0) wbody[2] = 0.0;
    else wbody[2] /= inertia[2];

    erotate += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
      inertia[2]*wbody[2]*wbody[2];
  }

  double erotateall;
  MPI_Allreduce(&erotate,&erotateall,1,MPI_DOUBLE,MPI_SUM,world);

  return 0.5*erotateall;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
   non-active DOF are removed by fflag/tflag and in tfactor
------------------------------------------------------------------------- */

double FixRigidSmall::compute_scalar()
{
  double wbody[3],rot[3][3];

  double *vcm,*inertia;

  double t = 0.0;

  for (int i = 0; i < nlocal_body; i++) {
    vcm = body[i].vcm;
    t += body[i].mass * (vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);

    // for Iw^2 rotational term, need wbody = angular velocity in body frame 
    // not omega = angular velocity in space frame

    inertia = body[i].inertia;
    MathExtra::quat_to_mat(body[i].quat,rot);
    MathExtra::transpose_matvec(rot,body[i].angmom,wbody);
    if (inertia[0] == 0.0) wbody[0] = 0.0;
    else wbody[0] /= inertia[0];
    if (inertia[1] == 0.0) wbody[1] = 0.0;
    else wbody[1] /= inertia[1];
    if (inertia[2] == 0.0) wbody[2] = 0.0;
    else wbody[2] /= inertia[2];

    t += inertia[0]*wbody[0]*wbody[0] + inertia[1]*wbody[1]*wbody[1] +
      inertia[2]*wbody[2]*wbody[2];
  }

  double tall;
  MPI_Allreduce(&t,&tall,1,MPI_DOUBLE,MPI_SUM,world);

  double tfactor = force->mvv2e / (6.0*nbody * force->boltz);
  tall *= tfactor;
  return tall;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixRigidSmall::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 2 * nmax * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);     // vatom
  if (extended) {
    bytes += nmax * sizeof(int);
    if (orientflag) bytes = nmax*orientflag * sizeof(double);
    if (dorientflag) bytes = nmax*3 * sizeof(double);
  }
  bytes += nmax_body * sizeof(Body);
  return bytes;
}

/* ----------------------------------------------------------------------
   debug method for sanity checking of atom/body data pointers
------------------------------------------------------------------------- */

/*
void FixRigidSmall::check(int flag)
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (bodyown[i] >= 0) {
      if (bodytag[i] != atom->tag[i]) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD AAA");
      }
      if (bodyown[i] < 0 || bodyown[i] >= nlocal_body) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD BBB");
      }
      if (atom2body[i] != bodyown[i]) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD CCC");
      }
      if (body[bodyown[i]].ilocal != i) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD DDD");
      }
    }
  }

  for (int i = 0; i < atom->nlocal; i++) {
    if (bodyown[i] < 0 && bodytag[i] > 0) {
      if (atom2body[i] < 0 || atom2body[i] >= nlocal_body+nghost_body) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD EEE");
      }
      if (bodytag[i] != atom->tag[body[atom2body[i]].ilocal]) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD FFF");
      }
    }
  }

  for (int i = atom->nlocal; i < atom->nlocal + atom->nghost; i++) {
    if (bodyown[i] >= 0) {
      if (bodyown[i] < nlocal_body || 
          bodyown[i] >= nlocal_body+nghost_body) {
        printf("Values %d %d: %d %d %d\n",
               i,atom->tag[i],bodyown[i],nlocal_body,nghost_body);
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD GGG");
      }
      if (body[bodyown[i]].ilocal != i) {
        printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
        errorx->one(FLERR,"BAD HHH");
      }
    }
  }

  for (int i = 0; i < nlocal_body; i++) {
    if (body[i].ilocal < 0 || body[i].ilocal >= atom->nlocal) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD III");
    }
    if (bodytag[body[i].ilocal] != atom->tag[body[i].ilocal] || 
        bodyown[body[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD JJJ");
    }
  }

  for (int i = nlocal_body; i < nlocal_body + nghost_body; i++) {
    if (body[i].ilocal < atom->nlocal || 
        body[i].ilocal >= atom->nlocal + atom->nghost) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD KKK");
    }
    if (bodyown[body[i].ilocal] != i) {
      printf("Proc %d, step %ld, flag %d\n",comm->me,update->ntimestep,flag);
      errorx->one(FLERR,"BAD LLL");
    }
  }
}
*/
