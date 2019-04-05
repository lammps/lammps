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
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)
   Based on fix_rigid (version from 2008).
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "fix_lb_rigid_pc_sphere.h"
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
#include "fix_lb_fluid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* -------------------------------------------------------------------------- */

FixLbRigidPCSphere::FixLbRigidPCSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int i, ibody;

  scalar_flag = 1;
  extscalar = 0;
  time_integrate = 1;
  rigid_flag = 1;
  create_attribute = 1;
  virial_flag = 1;
  thermo_virial = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  body = NULL;
  up = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // by default assume all of the particles interact with the fluid.
  inner_nodes = 0;

  // parse command-line args
  // set nbody and body[i] for each atom

  if (narg < 4) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

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
      error->all(FLERR,"Must use a molecular atom style with "
                 "fix lb/rigid/pc/sphere molecule");

    int *mask = atom->mask;
    tagint *molecule = atom->molecule;
    int nlocal = atom->nlocal;

    tagint maxmol_tag = -1;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) maxmol_tag = MAX(maxmol_tag,molecule[i]);

    tagint itmp;
    MPI_Allreduce(&maxmol_tag,&itmp,1,MPI_LMP_TAGINT,MPI_MAX,world);
    if (itmp+1 > MAXSMALLINT)
      error->all(FLERR,"Too many molecules for fix lb/rigid/pc/sphere");
    int maxmol = (int) itmp;

    int *ncount;
    memory->create(ncount,maxmol+1,"rigid:ncount");
    for (i = 0; i <= maxmol; i++) ncount[i] = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) ncount[molecule[i]]++;

    int *nall;
    memory->create(nall,maxmol+1,"rigid:ncount");
    MPI_Allreduce(ncount,nall,maxmol+1,MPI_LMP_TAGINT,MPI_SUM,world);

    nbody = 0;
    for (i = 0; i <= maxmol; i++)
      if (nall[i]) nall[i] = nbody++;
      else nall[i] = -1;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = nall[molecule[i]];
    }

    memory->destroy(ncount);
    memory->destroy(nall);

  // each listed group is a rigid body
  // check if all listed groups exist
  // an atom must belong to fix group and listed group to be in rigid body
  // error if atom belongs to more than 1 rigid body

  } else if (strcmp(arg[3],"group") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
    nbody = atoi(arg[4]);
    if (nbody <= 0) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
    if (narg < 5+nbody)
      error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
    iarg = 5 + nbody;

    int *igroups = new int[nbody];
    for (ibody = 0; ibody < nbody; ibody++) {
      igroups[ibody] = group->find(arg[5+ibody]);
      if (igroups[ibody] == -1)
        error->all(FLERR,"Could not find fix lb/rigid/pc/sphere group ID");
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
      error->all(FLERR,"One or more atoms belong to multiple rigid bodies");

    delete [] igroups;

  } else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

  // error check on nbody

  if (nbody == 0) error->all(FLERR,"No rigid bodies defined");

  // create all nbody-length arrays

  memory->create(nrigid,nbody,"lb/rigid/pc/sphere:nrigid");
  memory->create(nrigid_shell,nbody,"lb/rigid/pc/sphere:nrigid_shell");
  memory->create(masstotal,nbody,"lb/rigid/pc/sphere:masstotal");
  memory->create(masstotal_shell,nbody,"lb/rigid/pc/sphere:masstotal_shell");
  memory->create(sphereradius,nbody,"lb/rigid/pc/sphere:sphereradius");
  memory->create(xcm,nbody,3,"lb/rigid/pc/sphere:xcm");
  memory->create(xcm_old,nbody,3,"lb/rigid/pc/sphere:xcm_old");
  memory->create(vcm,nbody,3,"lb/rigid/pc/sphere:vcm");
  memory->create(ucm,nbody,3,"lb/rigid/pc/sphere:ucm");
  memory->create(ucm_old,nbody,3,"lb/rigid/pc/sphere:ucm_old");
  memory->create(fcm,nbody,3,"lb/rigid/pc/sphere:fcm");
  memory->create(fcm_old,nbody,3,"lb/rigid/pc/sphere:fcm_old");
  memory->create(fcm_fluid,nbody,3,"lb/rigid/pc/sphere:fcm_fluid");
  memory->create(omega,nbody,3,"lb/rigid/pc/sphere:omega");
  memory->create(torque,nbody,3,"lb/rigid/pc/sphere:torque");
  memory->create(torque_old,nbody,3,"lb/rigid/pc/sphere:torque_old");
  memory->create(torque_fluid,nbody,3,"lb/rigid/pc/sphere:torque_fluid");
  memory->create(torque_fluid_old,nbody,3,"lb/rigid/pc/sphere:torque_fluid_old");
  memory->create(rotate,nbody,3,"lb/rigid/pc/sphere:rotate");
  memory->create(imagebody,nbody,"lb/rigid/pc/sphere:imagebody");
  memory->create(fflag,nbody,3,"lb/rigid/pc/sphere:fflag");
  memory->create(tflag,nbody,3,"lb/rigid/pc/sphere:tflag");

  memory->create(sum,nbody,6,"lb/rigid/pc/sphere:sum");
  memory->create(all,nbody,6,"lb/rigid/pc/sphere:all");
  memory->create(remapflag,nbody,4,"lb/rigid/pc/sphere:remapflag");


  Gamma_MD = new double[nbody];

  // initialize force/torque flags to default = 1.0

  array_flag = 1;
  size_array_rows = nbody;
  size_array_cols = 15;
  global_freq = 1;
  extarray = 0;

  for (i = 0; i < nbody; i++) {
    fflag[i][0] = fflag[i][1] = fflag[i][2] = 1.0;
    tflag[i][0] = tflag[i][1] = tflag[i][2] = 1.0;
  }

  // parse optional args that set fflag and tflag

  while (iarg < narg) {
    if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      int mlo,mhi;
      force->bounds(FLERR,arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
      if (strcmp(arg[iarg+2],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
        fflag[m-1][0] = xflag;
        fflag[m-1][1] = yflag;
        fflag[m-1][2] = zflag;
        count++;
      }
      if (count == 0) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      iarg += 5;
    } else if (strcmp(arg[iarg],"torque") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      int mlo,mhi;
      force->bounds(FLERR,arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
      if (strcmp(arg[iarg+3],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
        tflag[m-1][0] = xflag;
        tflag[m-1][1] = yflag;
        tflag[m-1][2] = zflag;
        count++;
      }
      if (count == 0) error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");

      iarg += 5;
    // specify if certain particles are inside the rigid spherical body,
    // and therefore should not
    } else if(strcmp(arg[iarg],"innerNodes")==0){
      inner_nodes = 1;
      igroupinner = group->find(arg[iarg+1]);
      if(igroupinner == -1)
        error->all(FLERR,"Could not find fix lb/rigid/pc/sphere innerNodes group ID");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lb/rigid/pc/sphere command");
  }

  // initialize vector output quantities in case accessed before run

  for (i = 0; i < nbody; i++) {
    xcm[i][0] = xcm[i][1] = xcm[i][2] = 0.0;
    xcm_old[i][0] = xcm_old[i][1] = xcm_old[i][2] = 0.0;
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    ucm[i][0] = ucm[i][1] = ucm[i][2] = 0.0;
    ucm_old[i][0] = ucm_old[i][1] = ucm_old[i][2] = 0.0;
    fcm[i][0] = fcm[i][1] = fcm[i][2] = 0.0;
    fcm_old[i][0] = fcm_old[i][1] = fcm_old[i][2] = 0.0;
    fcm_fluid[i][0] = fcm_fluid[i][1] = fcm_fluid[i][2] = 0.0;
    torque[i][0] = torque[i][1] = torque[i][2] = 0.0;
    torque_old[i][0] = torque_old[i][1] = torque_old[i][2] = 0.0;
    torque_fluid[i][0] = torque_fluid[i][1] = torque_fluid[i][2] = 0.0;
    torque_fluid_old[i][0] = torque_fluid_old[i][1] = torque_fluid_old[i][2] = 0.0;
  }

  // nrigid[n] = # of atoms in Nth rigid body
  // error if one or zero atoms

  int *ncount = new int[nbody];
  for (ibody = 0; ibody < nbody; ibody++) ncount[ibody] = 0;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) ncount[body[i]]++;

  MPI_Allreduce(ncount,nrigid,nbody,MPI_INT,MPI_SUM,world);

  //count the number of atoms in the shell.
  if (inner_nodes == 1) {
    int *mask = atom->mask;
    for(ibody=0; ibody<nbody; ibody++) ncount[ibody] = 0;
    for(i=0; i<nlocal; i++){
      if(!(mask[i] & group->bitmask[igroupinner])){
        if(body[i] >= 0) ncount[body[i]]++;
      }
    }

    MPI_Allreduce(ncount,nrigid_shell,nbody,MPI_INT,MPI_SUM,world);
  } else {
    for(ibody=0; ibody < nbody; ibody++) nrigid_shell[ibody]=nrigid[ibody];
  }

  delete [] ncount;

  for (ibody = 0; ibody < nbody; ibody++)
    if (nrigid[ibody] <= 1) error->all(FLERR,"One or zero atoms in rigid body");

  // set image flags for each rigid body to default values
  // will be reset during init() based on xcm and then by pre_neighbor()
  // set here, so image value will persist from run to run

  for (ibody = 0; ibody < nbody; ibody++)
    imagebody[ibody] = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // print statistics

  int nsum = 0;
  for (ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];

  if (comm->me == 0) {
    if (screen) fprintf(screen,"%d rigid bodies with %d atoms\n",nbody,nsum);
    if (logfile) fprintf(logfile,"%d rigid bodies with %d atoms\n",nbody,nsum);
  }

  int groupbit_lb_fluid = 0;

  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lb/fluid")==0){
      fix_lb_fluid = (FixLbFluid *)modify->fix[ifix];
      groupbit_lb_fluid = group->bitmask[modify->fix[ifix]->igroup];
    }

   if(groupbit_lb_fluid == 0)
    error->all(FLERR,"the lb/fluid fix must also be used if using the lb/rigid/pc/sphere fix");

   int *mask = atom->mask;
   if(inner_nodes == 1){
     for(int j=0; j<nlocal; j++){
       if((mask[j] & groupbit) && !(mask[j] & group->bitmask[igroupinner]) && !(mask[j] & groupbit_lb_fluid))
         error->one(FLERR,"use the innerNodes keyword in the lb/rigid/pc/sphere fix for atoms which do not interact with the lb/fluid");

  // If inner nodes are present, which should not interact with the fluid, make
  // sure these are not used by the lb/fluid fix to apply a force to the fluid.
       if((mask[j] & groupbit) && (mask[j] & groupbit_lb_fluid) && (mask[j] & group->bitmask[igroupinner]))
         error->one(FLERR,"the inner nodes specified in lb/rigid/pc/sphere should not be included in the lb/fluid fix");
     }
   } else {
     for(int j=0; j<nlocal; j++){
       if((mask[j] & groupbit) && !(mask[j] & groupbit_lb_fluid))
         error->one(FLERR,"use the innerNodes keyword in the lb/rigid/pc/sphere fix for atoms which do not interact with the lb/fluid");
     }
   }

}

/* ---------------------------------------------------------------------- */

FixLbRigidPCSphere::~FixLbRigidPCSphere()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(body);
  memory->destroy(up);

  // delete nbody-length arrays

  memory->destroy(nrigid);
  memory->destroy(nrigid_shell);
  memory->destroy(masstotal);
  memory->destroy(masstotal_shell);
  memory->destroy(sphereradius);
  memory->destroy(xcm);
  memory->destroy(xcm_old);
  memory->destroy(vcm);
  memory->destroy(ucm);
  memory->destroy(ucm_old);
  memory->destroy(fcm);
  memory->destroy(fcm_old);
  memory->destroy(fcm_fluid);
  memory->destroy(omega);
  memory->destroy(torque);
  memory->destroy(torque_old);
  memory->destroy(torque_fluid);
  memory->destroy(torque_fluid_old);
  memory->destroy(rotate);
  memory->destroy(imagebody);
  memory->destroy(fflag);
  memory->destroy(tflag);

  memory->destroy(sum);
  memory->destroy(all);
  memory->destroy(remapflag);

  delete [] Gamma_MD;
}

/* ---------------------------------------------------------------------- */

int FixLbRigidPCSphere::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLbRigidPCSphere::init()
{
  int i,ibody;

  // warn if more than one rigid fix

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"lb/rigid/pc/sphere") == 0) count++;
  if (count > 1 && comm->me == 0) error->warning(FLERR,"More than one fix lb/rigid/pc/sphere");

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *periodicity = domain->periodicity;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double **mu = atom->mu;
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int extended = 0;
  int *mask = atom->mask;

  // Warn if any extended particles are included.
  if (atom->radius_flag || atom->ellipsoid_flag || atom->mu_flag) {
    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      if (ellipsoid && ellipsoid[i] >= 0) flag = 1;
      if (mu && mu[i][3] > 0.0) flag = 1;
    }

    MPI_Allreduce(&flag,&extended,1,MPI_INT,MPI_MAX,world);
  }
  if(extended)
    error->warning(FLERR,"Fix lb/rigid/pc/sphere assumes point particles");

  // compute masstotal & center-of-mass of each rigid body

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double massone,xunwrap,yunwrap,zunwrap;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    if ((xbox && !periodicity[0]) || (ybox && !periodicity[1]) ||
        (zbox && !periodicity[2]))
        error->one(FLERR,"Fix lb/rigid/pc/sphere atom has non-zero image flag "
                   "in a non-periodic dimension");

    xunwrap = x[i][0] + xbox*xprd;
    yunwrap = x[i][1] + ybox*yprd;
    zunwrap = x[i][2] + zbox*zprd;

    sum[ibody][0] += xunwrap * massone;
    sum[ibody][1] += yunwrap * massone;
    sum[ibody][2] += zunwrap * massone;
    sum[ibody][3] += massone;
    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        sum[ibody][4] += massone;
      }
    } else {
      sum[ibody][4] += massone;
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    masstotal[ibody] = all[ibody][3];
    masstotal_shell[ibody] = all[ibody][4];
    xcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    xcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    xcm[ibody][2] = all[ibody][2]/masstotal[ibody];
  }

  // Calculate the radius of the rigid body, and assign the value for gamma:
  double dx,dy,dz;
  double *Gamma = fix_lb_fluid->Gamma;
  double dm_lb = fix_lb_fluid->dm_lb;
  double dt_lb = fix_lb_fluid->dt_lb;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i=0; i<nlocal; i++){
    if(body[i] < 0) continue;
    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        ibody = body[i];

        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;

        xunwrap = x[i][0] + xbox*xprd;
        yunwrap = x[i][1] + ybox*yprd;
        zunwrap = x[i][2] + zbox*zprd;

        dx = xunwrap - xcm[ibody][0];
        dy = yunwrap - xcm[ibody][1];
        dz = zunwrap - xcm[ibody][2];

        sum[ibody][0] += dx*dx + dy*dy + dz*dz;
        sum[ibody][1] += Gamma[type[i]];
      }
    } else {
      ibody = body[i];

      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;

      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;

      dx = xunwrap - xcm[ibody][0];
      dy = yunwrap - xcm[ibody][1];
      dz = zunwrap - xcm[ibody][2];

      sum[ibody][0] += dx*dx + dy*dy + dz*dz;
      sum[ibody][1] += Gamma[type[i]];
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for(ibody=0; ibody < nbody; ibody++){
    sphereradius[ibody] = sqrt(all[ibody][0]/nrigid_shell[ibody]);
    Gamma_MD[ibody] = all[ibody][1]*dm_lb/dt_lb/nrigid_shell[ibody];
  }

  // Check that all atoms in the rigid body have the same value of gamma.
  double eps = 1.0e-7;
  for (i=0; i<nlocal; i++){
    if(body[i] < 0) continue;
    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        ibody = body[i];

        if(Gamma_MD[ibody]*dt_lb/dm_lb - Gamma[type[i]] > eps)
          error->one(FLERR,"All atoms in a rigid body must have the same gamma value");
      }
    } else {
      ibody = body[i];

        if(Gamma_MD[ibody]*dt_lb/dm_lb - Gamma[type[i]] > eps)
          error->one(FLERR,"All atoms in a rigid body must have the same gamma value");
    }
  }


  // remap the xcm of each body back into simulation box if needed
  // only really necessary the 1st time a run is performed

  pre_neighbor();


  // temperature scale factor

  double ndof = 0.0;
  for (ibody = 0; ibody < nbody; ibody++) {
    ndof += fflag[ibody][0] + fflag[ibody][1] + fflag[ibody][2];
    ndof += tflag[ibody][0] + tflag[ibody][1] + tflag[ibody][2];
  }
  if (ndof > 0.0) tfactor = force->mvv2e / (ndof * force->boltz);
  else tfactor = 0.0;

}
/* ---------------------------------------------------------------------- */

void FixLbRigidPCSphere::setup(int vflag)
{
  int i,n,ibody;
  double massone;

  // vcm = velocity of center-of-mass of each rigid body
  // fcm = force on center-of-mass of each rigid body

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  imageint *image = atom->image;

  double unwrap[3];
  double dx,dy,dz;



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

  // omega = angular velocity of each rigid body
  //         Calculated as the average of the angular velocities of the
  //         individual atoms comprising the rigid body.
  // torque = torque on each rigid body
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    domain->unmap(x[i],image[i],unwrap);

    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += (dy * (v[i][2]-vcm[ibody][2]) - dz * (v[i][1]-vcm[ibody][1]))/(dx*dx+dy*dy+dz*dz);
    sum[ibody][1] += (dz * (v[i][0]-vcm[ibody][0]) - dx * (v[i][2]-vcm[ibody][2]))/(dx*dx+dy*dy+dz*dz);
    sum[ibody][2] += (dx * (v[i][1]-vcm[ibody][1]) - dy * (v[i][0]-vcm[ibody][0]))/(dx*dx+dy*dy+dz*dz);
    sum[ibody][3] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][4] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][5] += dx * f[i][1] - dy * f[i][0];
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    omega[ibody][0] = all[ibody][0]/nrigid[ibody];
    omega[ibody][1] = all[ibody][1]/nrigid[ibody];
    omega[ibody][2] = all[ibody][2]/nrigid[ibody];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];
  }

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // Set the velocities
  set_v();

  if (evflag) {
    if (vflag_global)
      for (n = 0; n < 6; n++) virial[n] *= 2.0;
    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (n = 0; n < 6; n++)
          vatom[i][n] *= 2.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLbRigidPCSphere::initial_integrate(int vflag)
{
  double dtfm;

  int i,ibody;

  double massone;
  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  imageint *image = atom->image;

  double unwrap[3];
  double dx,dy,dz;

  int *mask = atom->mask;

  // compute the fluid velocity at the initial particle positions
  compute_up();

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  // Store the fluid velocity at the center of mass

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        sum[ibody][0] += up[i][0]*massone;
        sum[ibody][1] += up[i][1]*massone;
        sum[ibody][2] += up[i][2]*massone;
      }
    } else {
      sum[ibody][0] += up[i][0]*massone;
      sum[ibody][1] += up[i][1]*massone;
      sum[ibody][2] += up[i][2]*massone;
    }
  }
  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    ucm[ibody][0] = all[ibody][0]/masstotal_shell[ibody];
    ucm[ibody][1] = all[ibody][1]/masstotal_shell[ibody];
    ucm[ibody][2] = all[ibody][2]/masstotal_shell[ibody];
  }

  //Store the total torque due to the fluid.
  for (ibody = 0; ibody < nbody; ibody++)
    for(i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for(i = 0; i<nlocal; i++){
    if(body[i] < 0) continue;
    ibody = body[i];

    domain->unmap(x[i],image[i],unwrap);

    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        sum[ibody][0] += Gamma_MD[ibody]*(dy * ((up[i][2]-vcm[ibody][2])) -
                                          dz * ((up[i][1]-vcm[ibody][1])));
        sum[ibody][1] += Gamma_MD[ibody]*(dz * ((up[i][0]-vcm[ibody][0])) -
                                          dx * ((up[i][2]-vcm[ibody][2])));
        sum[ibody][2] += Gamma_MD[ibody]*(dx * ((up[i][1]-vcm[ibody][1])) -
                                          dy * ((up[i][0]-vcm[ibody][0])));
        sum[ibody][3] += -Gamma_MD[ibody]*(v[i][0]-up[i][0]);
        sum[ibody][4] += -Gamma_MD[ibody]*(v[i][1]-up[i][1]);
        sum[ibody][5] += -Gamma_MD[ibody]*(v[i][2]-up[i][2]);
      }
    } else {
      sum[ibody][0] += Gamma_MD[ibody]*(dy * ((up[i][2]-vcm[ibody][2])) -
                                        dz * ((up[i][1]-vcm[ibody][1])));
      sum[ibody][1] += Gamma_MD[ibody]*(dz * ((up[i][0]-vcm[ibody][0])) -
                                        dx * ((up[i][2]-vcm[ibody][2])));
      sum[ibody][2] += Gamma_MD[ibody]*(dx * ((up[i][1]-vcm[ibody][1])) -
                                        dy * ((up[i][0]-vcm[ibody][0])));
      sum[ibody][3] += -Gamma_MD[ibody]*(v[i][0]-up[i][0]);
      sum[ibody][4] += -Gamma_MD[ibody]*(v[i][1]-up[i][1]);
      sum[ibody][5] += -Gamma_MD[ibody]*(v[i][2]-up[i][2]);
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    torque_fluid[ibody][0] = all[ibody][0];
    torque_fluid[ibody][1] = all[ibody][1];
    torque_fluid[ibody][2] = all[ibody][2];
    fcm_fluid[ibody][0] = all[ibody][3];
    fcm_fluid[ibody][1] = all[ibody][4];
    fcm_fluid[ibody][2] = all[ibody][5];
  }

  for (int ibody = 0; ibody < nbody; ibody++) {
    fcm_old[ibody][0] = fcm[ibody][0];
    fcm_old[ibody][1] = fcm[ibody][1];
    fcm_old[ibody][2] = fcm[ibody][2];
    torque_old[ibody][0] = torque[ibody][0];
    torque_old[ibody][1] = torque[ibody][1];
    torque_old[ibody][2] = torque[ibody][2];
    torque_fluid_old[ibody][0] = torque_fluid[ibody][0];
    torque_fluid_old[ibody][1] = torque_fluid[ibody][1];
    torque_fluid_old[ibody][2] = torque_fluid[ibody][2];
    ucm_old[ibody][0] = ucm[ibody][0];
    ucm_old[ibody][1] = ucm[ibody][1];
    ucm_old[ibody][2] = ucm[ibody][2];
    xcm_old[ibody][0] = xcm[ibody][0];
    xcm_old[ibody][1] = xcm[ibody][1];
    xcm_old[ibody][2] = xcm[ibody][2];

    // update xcm by full step

    dtfm = dtf / masstotal[ibody];
    xcm[ibody][0] += dtv * vcm[ibody][0]+(fcm[ibody][0]+fcm_fluid[ibody][0]/force->ftm2v)*fflag[ibody][0]*dtfm*dtv;
    xcm[ibody][1] += dtv * vcm[ibody][1]+(fcm[ibody][1]+fcm_fluid[ibody][1]/force->ftm2v)*fflag[ibody][1]*dtfm*dtv;
    xcm[ibody][2] += dtv * vcm[ibody][2]+(fcm[ibody][2]+fcm_fluid[ibody][2]/force->ftm2v)*fflag[ibody][2]*dtfm*dtv;

    rotate[ibody][0] = omega[ibody][0]*dtv + tflag[ibody][0]*(torque[ibody][0]*force->ftm2v+torque_fluid[ibody][0])*
      dtv*dtv*5.0/(4.0*masstotal[ibody]*sphereradius[ibody]*sphereradius[ibody]);
    rotate[ibody][1] = omega[ibody][1]*dtv + tflag[ibody][1]*(torque[ibody][1]*force->ftm2v+torque_fluid[ibody][1])*
      dtv*dtv*5.0/(4.0*masstotal[ibody]*sphereradius[ibody]*sphereradius[ibody]);
    rotate[ibody][2] = omega[ibody][2]*dtv + tflag[ibody][2]*(torque[ibody][2]*force->ftm2v+torque_fluid[ibody][2])*
      dtv*dtv*5.0/(4.0*masstotal[ibody]*sphereradius[ibody]*sphereradius[ibody]);

    // Approximate vcm
    expminusdttimesgamma = exp(-Gamma_MD[ibody]*dtv*nrigid_shell[ibody]/masstotal[ibody]);
    force_factor = force->ftm2v/Gamma_MD[ibody]/nrigid_shell[ibody];

    if(fflag[ibody][0]==1){
      vcm[ibody][0] = expminusdttimesgamma*(vcm[ibody][0] - ucm[ibody][0] - fcm[ibody][0]*force_factor)
        + ucm[ibody][0] + fcm[ibody][0]*force_factor;
    }
    if(fflag[ibody][1]==1){
      vcm[ibody][1] = expminusdttimesgamma*(vcm[ibody][1] - ucm[ibody][1] - fcm[ibody][1]*force_factor) +
        ucm[ibody][1] + fcm[ibody][1]*force_factor;
    }
    if(fflag[ibody][2]==1){
      vcm[ibody][2] = expminusdttimesgamma*(vcm[ibody][2] - ucm[ibody][2] - fcm[ibody][2]*force_factor) +
        ucm[ibody][2] + fcm[ibody][2]*force_factor;
    }

    // Approximate angmom
    torque_factor = 5.0*Gamma_MD[ibody]*nrigid_shell[ibody]/(3.0*masstotal[ibody]);
    expminusdttimesgamma = exp(-dtv*torque_factor);


    if(tflag[ibody][0]==1){
      omega[ibody][0] = expminusdttimesgamma*(omega[ibody][0] - (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][0] + torque_fluid[ibody][0])) +
                                              (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][0] + torque_fluid[ibody][0]);
    }
    if(tflag[ibody][1]==1){
      omega[ibody][1] = expminusdttimesgamma*(omega[ibody][1] - (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][1] + torque_fluid[ibody][1])) +
                                              (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][1] + torque_fluid[ibody][1]);
    }
    if(tflag[ibody][2]==1){
      omega[ibody][2] = expminusdttimesgamma*(omega[ibody][2] - (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][2] + torque_fluid[ibody][2])) +
                                              (3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
                                              (force->ftm2v*torque[ibody][2] + torque_fluid[ibody][2]);
    }

  }
  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  set_xv();

}


/* ---------------------------------------------------------------------- */

void FixLbRigidPCSphere::final_integrate()
{
  int i,ibody;

  // sum over atoms to get force and torque on rigid body
  double massone;
  imageint *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double unwrap[3];
  double dx,dy,dz;

  int *mask = atom->mask;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];

    domain->unmap(x[i],image[i],unwrap);

    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][3] += dy*f[i][2] - dz*f[i][1];
    sum[ibody][4] += dz*f[i][0] - dx*f[i][2];
    sum[ibody][5] += dx*f[i][1] - dy*f[i][0];
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  //Compute the correction to the velocity and angular momentum due to the non-fluid forces:
  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];

    expminusdttimesgamma = exp(-dtv*Gamma_MD[ibody]*nrigid_shell[ibody]/masstotal[ibody]);
    DMDcoeff= (dtv - (masstotal[ibody]/nrigid_shell[ibody])*(1.0-expminusdttimesgamma)/Gamma_MD[ibody]);

    vcm[ibody][0] += fflag[ibody][0]*DMDcoeff*force->ftm2v*(fcm[ibody][0]-fcm_old[ibody][0])/Gamma_MD[ibody]/dtv/nrigid_shell[ibody];
    vcm[ibody][1] += fflag[ibody][1]*DMDcoeff*force->ftm2v*(fcm[ibody][1]-fcm_old[ibody][1])/Gamma_MD[ibody]/dtv/nrigid_shell[ibody];
    vcm[ibody][2] += fflag[ibody][2]*DMDcoeff*force->ftm2v*(fcm[ibody][2]-fcm_old[ibody][2])/Gamma_MD[ibody]/dtv/nrigid_shell[ibody];

    torque_factor = 5.0*Gamma_MD[ibody]*nrigid_shell[ibody]/(3.0*masstotal[ibody]);
    expminusdttimesgamma = exp(-dtv*torque_factor);
    DMDcoeff = (dtv - (1.0-expminusdttimesgamma)/torque_factor);

    omega[ibody][0] += tflag[ibody][0]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*DMDcoeff*
                                        force->ftm2v*(torque[ibody][0] - torque_old[ibody][0])/dtv;
    omega[ibody][1] += tflag[ibody][1]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*DMDcoeff*
                                        force->ftm2v*(torque[ibody][1] - torque_old[ibody][1])/dtv;
    omega[ibody][2] += tflag[ibody][2]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*DMDcoeff*
                                        force->ftm2v*(torque[ibody][2] - torque_old[ibody][2])/dtv;
  }

  //Next, calculate the correction to the velocity and angular momentum due to the fluid forces:
  //Calculate the fluid velocity at the new particle locations.
  compute_up();

  // store fluid quantities for the total body
   for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
   // Store the fluid velocity at the center of mass, and the total force
   // due to the fluid.
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    domain->unmap(x[i],image[i],unwrap);

    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    if(inner_nodes == 1){
      if(!(mask[i] & group->bitmask[igroupinner])){
        sum[ibody][0] += up[i][0]*massone;
        sum[ibody][1] += up[i][1]*massone;
        sum[ibody][2] += up[i][2]*massone;
        sum[ibody][3] += Gamma_MD[ibody]*(dy * ((up[i][2]-vcm[ibody][2])) -
                                          dz * ((up[i][1]-vcm[ibody][1])));
        sum[ibody][4] += Gamma_MD[ibody]*(dz * ((up[i][0]-vcm[ibody][0])) -
                                          dx * ((up[i][2]-vcm[ibody][2])));
        sum[ibody][5] += Gamma_MD[ibody]*(dx * ((up[i][1]-vcm[ibody][1])) -
                                          dy * ((up[i][0]-vcm[ibody][0])));
      }
    } else {
      sum[ibody][0] += up[i][0]*massone;
      sum[ibody][1] += up[i][1]*massone;
      sum[ibody][2] += up[i][2]*massone;
      sum[ibody][3] += Gamma_MD[ibody]*(dy * ((up[i][2]-vcm[ibody][2])) -
                                        dz * ((up[i][1]-vcm[ibody][1])));
      sum[ibody][4] += Gamma_MD[ibody]*(dz * ((up[i][0]-vcm[ibody][0])) -
                                        dx * ((up[i][2]-vcm[ibody][2])));
      sum[ibody][5] += Gamma_MD[ibody]*(dx * ((up[i][1]-vcm[ibody][1])) -
                                        dy * ((up[i][0]-vcm[ibody][0])));
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    ucm[ibody][0] = all[ibody][0]/masstotal_shell[ibody];
    ucm[ibody][1] = all[ibody][1]/masstotal_shell[ibody];
    ucm[ibody][2] = all[ibody][2]/masstotal_shell[ibody];
    torque_fluid[ibody][0] = all[ibody][3];
    torque_fluid[ibody][1] = all[ibody][4];
    torque_fluid[ibody][2] = all[ibody][5];
  }

  for (ibody = 0; ibody < nbody; ibody++) {

    expminusdttimesgamma = exp(-dtv*Gamma_MD[ibody]*nrigid_shell[ibody]/masstotal[ibody]);
    DMDcoeff= (dtv - (masstotal[ibody]/nrigid_shell[ibody])*(1.0-expminusdttimesgamma)/Gamma_MD[ibody]);

    vcm[ibody][0] += DMDcoeff*fflag[ibody][0]*(ucm[ibody][0]-ucm_old[ibody][0])/dtv;
    vcm[ibody][1] += DMDcoeff*fflag[ibody][1]*(ucm[ibody][1]-ucm_old[ibody][1])/dtv;
    vcm[ibody][2] += DMDcoeff*fflag[ibody][2]*(ucm[ibody][2]-ucm_old[ibody][2])/dtv;

    torque_factor = 5.0*Gamma_MD[ibody]*nrigid_shell[ibody]/(3.0*masstotal[ibody]);
    expminusdttimesgamma = exp(-dtv*torque_factor);
    DMDcoeff = (dtv - (1.0-expminusdttimesgamma)/torque_factor);

    omega[ibody][0] += tflag[ibody][0]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
      DMDcoeff*(torque_fluid[ibody][0] - torque_fluid_old[ibody][0])/dtv;
    omega[ibody][1] += tflag[ibody][1]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
      DMDcoeff*(torque_fluid[ibody][1] - torque_fluid_old[ibody][1])/dtv;
    omega[ibody][2] += tflag[ibody][2]*(3.0/(2.0*nrigid_shell[ibody]*sphereradius[ibody]*sphereradius[ibody]*Gamma_MD[ibody]))*
      DMDcoeff*(torque_fluid[ibody][2] - torque_fluid_old[ibody][2])/dtv;
  }

  set_v();

}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixLbRigidPCSphere::set_v()
{
  int ibody;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xunwrap,yunwrap,zunwrap;

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

    xunwrap = x[i][0] + xbox*xprd;
    yunwrap = x[i][1] + ybox*yprd;
    zunwrap = x[i][2] + zbox*zprd;

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];

    // save old velocities for virial.
    if(evflag){
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = (omega[ibody][1]*dz - omega[ibody][2]*dy) + vcm[ibody][0];
    v[i][1] = (omega[ibody][2]*dx - omega[ibody][0]*dz) + vcm[ibody][1];
    v[i][2] = (omega[ibody][0]*dy - omega[ibody][1]*dx) + vcm[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0] + Gamma_MD[ibody]*(v0-up[i][0]);
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1] + Gamma_MD[ibody]*(v1-up[i][1]);
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2] + Gamma_MD[ibody]*(v2-up[i][2]);

      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;

      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }


  }

}
/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixLbRigidPCSphere::set_xv()
{
  int ibody;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xunwrap,yunwrap,zunwrap;
  double dx,dy,dz;


  // set x and v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

    xunwrap = x[i][0] + xbox*xprd;
    yunwrap = x[i][1] + ybox*yprd;
    zunwrap = x[i][2] + zbox*zprd;

    dx = xunwrap - xcm_old[ibody][0];
    dy = yunwrap - xcm_old[ibody][1];
    dz = zunwrap - xcm_old[ibody][2];

    // save old positions and velocities for virial
    if(evflag){
      x0 = xunwrap;
      x1 = yunwrap;
      x2 = zunwrap;
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    x[i][0] = dx;
    x[i][1] = dy;
    x[i][2] = dz;

    //Perform the rotations:
    dx = x[i][0];
    dy = x[i][1];
    dz = x[i][2];
    x[i][0] = cos(rotate[ibody][2])*dx - sin(rotate[ibody][2])*dy;
    x[i][1] = sin(rotate[ibody][2])*dx + cos(rotate[ibody][2])*dy;

    dx = x[i][0];
    dy = x[i][1];
    dz = x[i][2];
    x[i][0] = cos(rotate[ibody][1])*dx + sin(rotate[ibody][1])*dz;
    x[i][2] = -sin(rotate[ibody][1])*dx + cos(rotate[ibody][1])*dz;

    dx = x[i][0];
    dy = x[i][1];
    dz = x[i][2];
    x[i][1] = cos(rotate[ibody][0])*dy - sin(rotate[ibody][0])*dz;
    x[i][2] = sin(rotate[ibody][0])*dy + cos(rotate[ibody][0])*dz;

    v[i][0] = (omega[ibody][1]*x[i][2] - omega[ibody][2]*x[i][1]) + vcm[ibody][0];
    v[i][1] = (omega[ibody][2]*x[i][0] - omega[ibody][0]*x[i][2]) + vcm[ibody][1];
    v[i][2] = (omega[ibody][0]*x[i][1] - omega[ibody][1]*x[i][0]) + vcm[ibody][2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well
    x[i][0] += xcm[ibody][0] - xbox*xprd;
    x[i][1] += xcm[ibody][1] - ybox*yprd;
    x[i][2] += xcm[ibody][2] - zbox*zprd;

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0] + Gamma_MD[ibody]*(v0-up[i][0]);
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1] + Gamma_MD[ibody]*(v1-up[i][1]);
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2] + Gamma_MD[ibody]*(v2-up[i][2]);

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }

  }


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

void FixLbRigidPCSphere::pre_neighbor()
{
  imageint original,oldimage,newimage;

  for (int ibody = 0; ibody < nbody; ibody++) {
    original = imagebody[ibody];
    domain->remap(xcm[ibody],imagebody[ibody]);

    if (original == imagebody[ibody]) {
      remapflag[ibody][3] = 0;
    } else {
      oldimage = original & IMGMASK;
      newimage = imagebody[ibody] & IMGMASK;
      remapflag[ibody][0] = newimage - oldimage;
      oldimage = (original >> IMGBITS) & IMGMASK;
      newimage = (imagebody[ibody] >> IMGBITS) & IMGMASK;
      remapflag[ibody][1] = newimage - oldimage;
      oldimage = original >> IMG2BITS;
      newimage = imagebody[ibody] >> IMG2BITS;
      remapflag[ibody][2] = newimage - oldimage;
      remapflag[ibody][3] = 1;
    }
  }

  // adjust image flags of any atom in a rigid body whose xcm was remapped

  imageint *atomimage = atom->image;
  int nlocal = atom->nlocal;

  int ibody;
  imageint idim,otherdims;

  for (int i = 0; i < nlocal; i++) {
    if (body[i] == -1) continue;
    if (remapflag[body[i]][3] == 0) continue;
    ibody = body[i];

    if (remapflag[ibody][0]) {
      idim = atomimage[i] & IMGMASK;
      otherdims = atomimage[i] ^ idim;
      idim -= remapflag[ibody][0];
      idim &= IMGMASK;
      atomimage[i] = otherdims | idim;
    }
    if (remapflag[ibody][1]) {
      idim = (atomimage[i] >> IMGBITS) & IMGMASK;
      otherdims = atomimage[i] ^ (idim << IMGBITS);
      idim -= remapflag[ibody][1];
      idim &= IMGMASK;
      atomimage[i] = otherdims | (idim << IMGBITS);
    }
    if (remapflag[ibody][2]) {
      idim = atomimage[i] >> IMG2BITS;
      otherdims = atomimage[i] ^ (idim << IMG2BITS);
      idim -= remapflag[ibody][2];
      idim &= IMGMASK;
      atomimage[i] = otherdims | (idim << IMG2BITS);
    }
  }
}
/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_rigid for atoms in igroup
------------------------------------------------------------------------- */

int FixLbRigidPCSphere::dof(int igroup)
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
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixLbRigidPCSphere::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixLbRigidPCSphere::grow_arrays(int nmax)
{

  memory->grow(body,nmax,"rigid:body");
  memory->grow(up,nmax,3,"rigid:up");

}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixLbRigidPCSphere::copy_arrays(int i, int j, int /*delflag*/)
{
  body[j] = body[i];
  up[j][0] = up[i][0];
  up[j][1] = up[i][1];
  up[j][2] = up[i][2];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixLbRigidPCSphere::set_arrays(int i)
{
  body[i] = -1;
  up[i][0] = 0.0;
  up[i][1] = 0.0;
  up[i][2] = 0.0;

}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixLbRigidPCSphere::pack_exchange(int i, double *buf)
{
  buf[0] = body[i];
  buf[1] = up[i][0];
  buf[2] = up[i][1];
  buf[3] = up[i][2];
  return 4;

}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixLbRigidPCSphere::unpack_exchange(int nlocal, double *buf)
{
  body[nlocal] = static_cast<int> (buf[0]);
  up[nlocal][0] = buf[1];
  up[nlocal][1] = buf[2];
  up[nlocal][2] = buf[3];
  return 4;
}

/* ---------------------------------------------------------------------- */

void FixLbRigidPCSphere::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
   non-active DOF are removed by fflag/tflag and in tfactor
------------------------------------------------------------------------- */

double FixLbRigidPCSphere::compute_scalar()
{

  double inertia;
  double t = 0.0;

  for (int i = 0; i < nbody; i++) {
    t += masstotal[i] * (fflag[i][0]*vcm[i][0]*vcm[i][0] +
                         fflag[i][1]*vcm[i][1]*vcm[i][1] +      \
                         fflag[i][2]*vcm[i][2]*vcm[i][2]);

    // wbody = angular velocity in body frame

    inertia = 2.0*masstotal[i]*sphereradius[i]*sphereradius[i]/5.0;

    t += tflag[i][0]*inertia*omega[i][0]*omega[i][0] +
      tflag[i][1]*inertia*omega[i][1]*omega[i][1] +
      tflag[i][2]*inertia*omega[i][2]*omega[i][2];
  }

  t *= tfactor;
  return t;
}


/* ----------------------------------------------------------------------
   return attributes of a rigid body
   15 values per body
   xcm = 0,1,2; vcm = 3,4,5; fcm = 6,7,8; torque = 9,10,11; image = 12,13,14
------------------------------------------------------------------------- */

double FixLbRigidPCSphere::compute_array(int i, int j)
{
  if (j < 3) return xcm[i][j];
  if (j < 6) return vcm[i][j-3];
  if (j < 9) return (fcm[i][j-6]+fcm_fluid[i][j-6]);
  if (j < 12) return (torque[i][j-9]+torque_fluid[i][j-9]);
  if (j == 12) return (imagebody[i] & IMGMASK) - IMGMAX;
  if (j == 13) return (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
  return (imagebody[i] >> IMG2BITS) - IMGMAX;
}

/* ---------------------------------------------------------------------- */
 void FixLbRigidPCSphere::compute_up(void)
 {
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   double **x = atom->x;
   int i,k;
   int ix,iy,iz;
   int ixp,iyp,izp;
   double dx1,dy1,dz1;
   int isten,ii,jj,kk;
   double r,rsq,weightx,weighty,weightz;
   double ****u_lb = fix_lb_fluid->u_lb;
   int subNbx = fix_lb_fluid->subNbx;
   int subNby = fix_lb_fluid->subNby;
   int subNbz = fix_lb_fluid->subNbz;
   double dx_lb = fix_lb_fluid->dx_lb;
   double dt_lb = fix_lb_fluid->dt_lb;
   int trilinear_stencil = fix_lb_fluid->trilinear_stencil;
   double FfP[64];


  for(i=0; i<nlocal; i++){
    if(mask[i] & groupbit){

      //Calculate nearest leftmost grid point.
      //Since array indices from 1 to subNb-2 correspond to the
      // local subprocessor domain (not indices from 0), use the
      // ceiling value.
      ix = (int)ceil((x[i][0]-domain->sublo[0])/dx_lb);
      iy = (int)ceil((x[i][1]-domain->sublo[1])/dx_lb);
      iz = (int)ceil((x[i][2]-domain->sublo[2])/dx_lb);

      //Calculate distances to the nearest points.
      dx1 = x[i][0] - (domain->sublo[0] + (ix-1)*dx_lb);
      dy1 = x[i][1] - (domain->sublo[1] + (iy-1)*dx_lb);
      dz1 = x[i][2] - (domain->sublo[2] + (iz-1)*dx_lb);

      // Need to convert these to lattice units:
      dx1 = dx1/dx_lb;
      dy1 = dy1/dx_lb;
      dz1 = dz1/dx_lb;


      up[i][0]=0.0; up[i][1]=0.0; up[i][2]=0.0;

      if(trilinear_stencil==0){
        isten=0;
        for(ii=-1; ii<3; ii++){
          rsq=(-dx1+ii)*(-dx1+ii);

          if(rsq>=4) {
            weightx=0.0;
          } else {
            r=sqrt(rsq);
            if(rsq>1){
              weightx=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
            } else {
              weightx=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
            }
          }
          for(jj=-1; jj<3; jj++){
            rsq=(-dy1+jj)*(-dy1+jj);
            if(rsq>=4) {
              weighty=0.0;
            } else {
              r=sqrt(rsq);
              if(rsq>1){
                weighty=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
              } else {
                weighty=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
              }
            }
            for(kk=-1; kk<3; kk++){
              rsq=(-dz1+kk)*(-dz1+kk);
              if(rsq>=4) {
                weightz=0.0;
              } else {
                r=sqrt(rsq);
                if(rsq>1){
                  weightz=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
                } else {
                  weightz=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
                }
              }
              ixp = ix+ii;
              iyp = iy+jj;
              izp = iz+kk;


              if(ixp==-1) ixp=subNbx+2;
              if(iyp==-1) iyp=subNby+2;
              if(izp==-1) izp=subNbz+2;

              FfP[isten] = weightx*weighty*weightz;
              // interpolated velocity based on delta function.
              for(k=0; k<3; k++){
                up[i][k] += u_lb[ixp][iyp][izp][k]*FfP[isten];
              }
            }
          }
        }
      } else {
        FfP[0] = (1.-dx1)*(1.-dy1)*(1.-dz1);
        FfP[1] = (1.-dx1)*(1.-dy1)*dz1;
        FfP[2] = (1.-dx1)*dy1*(1.-dz1);
        FfP[3] = (1.-dx1)*dy1*dz1;
        FfP[4] = dx1*(1.-dy1)*(1.-dz1);
        FfP[5] = dx1*(1.-dy1)*dz1;
        FfP[6] = dx1*dy1*(1.-dz1);
        FfP[7] = dx1*dy1*dz1;

        ixp = (ix+1);
        iyp = (iy+1);
        izp = (iz+1);

        for (k=0; k<3; k++) {   // tri-linearly interpolated velocity at node
          up[i][k] = u_lb[ix][iy][iz][k]*FfP[0]
            + u_lb[ix][iy][izp][k]*FfP[1]
            + u_lb[ix][iyp][iz][k]*FfP[2]
            + u_lb[ix][iyp][izp][k]*FfP[3]
            + u_lb[ixp][iy][iz][k]*FfP[4]
            + u_lb[ixp][iy][izp][k]*FfP[5]
            + u_lb[ixp][iyp][iz][k]*FfP[6]
            + u_lb[ixp][iyp][izp][k]*FfP[7];
        }
      }
      for(k=0; k<3; k++)
        up[i][k] = up[i][k]*dx_lb/dt_lb;

    }
  }

 }
