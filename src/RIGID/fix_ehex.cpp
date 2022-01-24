// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Peter Wirnsberger (University of Cambridge)

   This source file implements the asymmetric version of the enhanced heat
   exchange (eHEX/a) algorithm. The paper is available for download on
   arXiv: https://arxiv.org/pdf/1507.07081.pdf.

   This file is based on fix_heat.cpp written by Paul Crozier (SNL)
   which implements the heat exchange (HEX) algorithm.
------------------------------------------------------------------------- */

#include "fix_ehex.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "fix_shake.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixEHEX::FixEHEX(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
  idregion(nullptr), x(nullptr), f(nullptr), v(nullptr),
  mass(nullptr), rmass(nullptr), type(nullptr), scalingmask(nullptr)
{
  MPI_Comm_rank(world, &me);

  // check
  if (narg < 4) error->all(FLERR,"Illegal fix ehex command: wrong number of parameters ");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  // apply fix every nevery timesteps

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);

  if (nevery <= 0) error->all(FLERR,"Illegal fix ehex command");

  // heat flux into the reservoir

  heat_input = utils::numeric(FLERR,arg[4],false,lmp);

  // optional args

  iregion = -1;

  // NOTE: constraints are deactivated by default

  constraints = 0;

  // NOTE: cluster rescaling is deactivated by default

  cluster = 0;

  // NOTE: hex = 1 means that no coordinate correction is applied in which case eHEX reduces to HEX

  hex = 0;

  int iarg = 5;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ehex command: wrong number of parameters ");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix ehex does not exist");
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;
    }

    // apply constraints (shake/rattle) at the end of the timestep

    else if (strcmp(arg[iarg], "constrain") == 0) {
      constraints = 1;
      iarg += 1;
    }

    // rescale only if the entire molecule is contained within the region

    else if (strcmp(arg[iarg], "com") == 0) {
      cluster = 1;
      iarg += 1;
    }

    // don't apply a coordinate correction if this keyword is specified

    else if (strcmp(arg[iarg], "hex") == 0) {
      hex  = 1;
      iarg+= 1;
    }
    else
      error->all(FLERR, "Illegal fix ehex keyword ");
  }

  // check options

  if (cluster && !constraints)
    error->all(FLERR, "You can only use the keyword 'com' together with the keyword 'constrain' ");

  scale = 1.0;
  scalingmask    = nullptr;
  FixEHEX::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

}


/* ---------------------------------------------------------------------- */

void FixEHEX::grow_arrays(int nmax) {
  memory->grow(scalingmask, nmax,"ehex:scalingmask");
}


/* ---------------------------------------------------------------------- */

FixEHEX::~FixEHEX()
{
  atom->delete_callback(id,Atom::GROW);
  delete [] idregion;
  memory->destroy(scalingmask);

}

/* ---------------------------------------------------------------------- */

int FixEHEX::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEHEX::init()
{
  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix ehex does not exist");
  }

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix ehex group has no atoms");

  fshake = nullptr;
  if (constraints) {

    // check if constraining algorithm is used (FixRattle inherits from FixShake)

    int cnt_shake = 0;
    int id_shake;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp("rattle", modify->fix[i]->style) == 0 ||
          strcmp("shake", modify->fix[i]->style) == 0) {
        cnt_shake++;
        id_shake = i;
      }
    }

    if (cnt_shake > 1)
      error->all(FLERR,"Multiple instances of fix shake/rattle detected (not supported yet)");
    else if (cnt_shake == 1)   {
     fshake = ((FixShake*) modify->fix[id_shake]);
    }
    else if (cnt_shake == 0)
      error->all(FLERR, "Fix ehex was configured with keyword constrain, but shake/rattle was not defined");
  }
}



/* ---------------------------------------------------------------------- */


void FixEHEX::end_of_step() {
  // store local pointers

  x       = atom->x;
  f       = atom->f;
  v       = atom->v;
  mass    = atom->mass;
  rmass   = atom->rmass;
  type    = atom->type;
  nlocal  = atom->nlocal;

  // determine which sites are to be rescaled

  update_scalingmask();

  // rescale velocities

  rescale();

  // if required use shake/rattle to correct coordinates and velocities

  if (constraints && fshake)
    fshake->shake_end_of_step(0);
}



/* ----------------------------------------------------------------------
   Iterate over all atoms, rescale the velocities and apply coordinate
   corrections.
------------------------------------------------------------------------- */

void FixEHEX::rescale() {
  double Kr, Ke, escale;
  double vsub[3],vcm[3], sfr[3];
  double mi;
  double dt;
  double F, mr, epsr_ik, sfvr, eta_ik;

  dt = update->dt;

  // calculate centre of mass properties

  com_properties(vcm, sfr, &sfvr, &Ke, &Kr, &masstotal);

  // heat flux into the reservoir

  F     = heat_input * force->ftm2v * nevery;

  // total mass

  mr    = masstotal;

  // energy scaling factor

  escale = 1. + (F*dt)/Kr;

  // safety check for kinetic energy

  if (escale < 0.0) error->all(FLERR,"Fix ehex kinetic energy went negative");

  scale = sqrt(escale);
  vsub[0] = (scale-1.0) * vcm[0];
  vsub[1] = (scale-1.0) * vcm[1];
  vsub[2] = (scale-1.0) * vcm[2];

  for (int i = 0; i < nlocal; i++) {

    if (scalingmask[i]) {

      mi = (rmass) ? rmass[i] :  mass[type[i]];

      for (int k=0; k<3; k++) {

        // apply coordinate correction unless running in hex mode

        if (!hex) {

            // epsr_ik implements Eq. (20) in the paper

            eta_ik    = mi * F/(2.*Kr) * (v[i][k] - vcm[k]);
            epsr_ik   = eta_ik / (mi*Kr) * (F/48. + sfvr/6.*force->ftm2v) - F/(12.*Kr) * (f[i][k]/mi - sfr[k]/mr)*force->ftm2v;

            x[i][k]  -= dt*dt*dt * epsr_ik;
        }

        // rescale the velocity

        v[i][k]   = scale*v[i][k] - vsub[k];
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

double FixEHEX::compute_scalar()
{
  return scale;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixEHEX::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)atom->nmax * sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
   Update the array scalingmask depending on which individual atoms
   will be rescaled or not.
------------------------------------------------------------------------- */

void FixEHEX::update_scalingmask() {
  int m;
  int lid;
  bool stat;
  int nsites;

  // prematch region

  Region *region = nullptr;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // only rescale molecules whose center of mass if fully contained in the region

  if (cluster) {

    // loop over all clusters

    for (int i=0; i < fshake->nlist; i++) {

      // cluster id

      m    = fshake->list[i];

      // check if the centre of mass of the cluster is inside the region
      // if region == nullptr, just check the group information of all sites

      if      (fshake->shake_flag[m] == 1)      nsites = 3;
      else if (fshake->shake_flag[m] == 2)      nsites = 2;
      else if (fshake->shake_flag[m] == 3)      nsites = 3;
      else if (fshake->shake_flag[m] == 4)      nsites = 4;
      else                                      nsites = 0;

      if (nsites == 0) {
        error->all(FLERR,"Internal error: shake_flag[m] has to be between 1 and 4 for m in nlist");
      }

      stat = check_cluster(fshake->shake_atom[m], nsites, region);

      for (int l=0; l < nsites; l++)  {
        lid = atom->map(fshake->shake_atom[m][l]);
        scalingmask[lid] = stat;
      }
    }

    // check atoms that do not belong to any cluster

    for (int i=0; i<atom->nlocal; i++)  {
      if (fshake->shake_flag[i] == 0)
        scalingmask[i] = rescale_atom(i,region);
    }

  }

  // no clusters, just individual sites (e.g. monatomic system or flexible molecules)

  else {
    for (int i=0; i<atom->nlocal; i++)
      scalingmask[i] =  rescale_atom(i,region);
  }

}


/* ----------------------------------------------------------------------
   Check if the centre of mass of the cluster to be constrained is
   inside the region.
------------------------------------------------------------------------- */

bool FixEHEX::check_cluster(tagint *shake_atom, int n, Region * region) {

  // IMPORTANT NOTE: If any site of the cluster belongs to a group
  //                 which should not be rescaled than all of the sites
  //                 will be ignored!

  double **x     = atom->x;
  double * rmass = atom->rmass;
  double * mass  = atom->mass;
  int    * type  = atom->type;
  int    * mask  = atom->mask;
  double   xcom[3], xtemp[3];
  double   mcluster, mi;
  bool     stat;
  int      lid[4];

  // accumulate mass and centre of mass position

  stat      = true;
  xcom[0]   = 0.;
  xcom[1]   = 0.;
  xcom[2]   = 0.;
  mcluster  = 0;

  for (int i = 0; i < n; i++) {

    // get local id

    lid[i] = atom->map(shake_atom[i]);

    // check if all sites of the cluster belong to the correct group

    stat = stat && (mask[lid[i]] & groupbit);

    if (region && stat)  {

      // check if reduced mass is used

      mi        = (rmass) ? rmass[lid[i]] : mass[type[lid[i]]];
      mcluster += mi;

      // accumulate centre of mass
      // NOTE: you can either use unwrapped coordinates or take site x[lid[0]] as reference,
      //       i.e. reconstruct the molecule around this site and calculate the com.

      for (int k=0; k<3; k++)
        xtemp[k] = x[lid[i]][k] - x[lid[0]][k];

      // take into account pbc

      domain->minimum_image(xtemp);

      for (int k=0; k<3; k++)
        xcom[k] += mi * (x[lid[0]][k] + xtemp[k]) ;
    }
  }

  // check if centre of mass is inside the region (if specified)

  if (region && stat) {

    // check mass

    if (mcluster < 1.e-14) {
      error->all(FLERR, "Fix ehex shake cluster has almost zero mass.");
    }

    // divide by total mass

    for (int k=0; k<3; k++)
      xcom[k] = xcom[k]/mcluster;

    // apply periodic boundary conditions (centre of mass could be outside the box)
    // and check if molecule is inside the region

    domain->remap(xcom);
    stat = stat && region->match(xcom[0], xcom[1], xcom[2]);
  }

  return stat;
}


/* ----------------------------------------------------------------------
   Check if atom i has the correct group and is inside the region.
------------------------------------------------------------------------- */

bool FixEHEX::rescale_atom(int i, Region*region) {
  bool stat;
  double x_r[3];

  // check mask and group

  stat = (atom->mask[i] & groupbit);

  if (region) {

    x_r[0] = atom->x[i][0];
    x_r[1] = atom->x[i][1];
    x_r[2] = atom->x[i][2];

    // apply periodic boundary conditions

    domain->remap(x_r);

    // check if the atom is in the group/region

    stat = stat && region->match(x_r[0],x_r[1],x_r[2]);
  }

  return stat;
}

/* ----------------------------------------------------------------------
   Calculate global properties of the atoms inside the reservoir.
   (e.g. com velocity, kinetic energy, total mass,...)
------------------------------------------------------------------------- */

void FixEHEX::com_properties(double * vr, double * sfr, double *sfvr, double *K, double *Kr, double *mr) {
   double ** f  = atom->f;
   double ** v  = atom->v;
   int nlocal   = atom->nlocal;
   double *rmass= atom->rmass;
   double *mass = atom->mass;
   int    *type = atom->type;
   double l_vr[3];
   double l_mr;
   double l_sfr[3];
   double l_sfvr;
   double l_K;
   double mi;
   double l_buf[9];
   double buf[9];

   // calculate partial sums

   l_vr[0]  = l_vr[1]  = l_vr[2] = 0;
   l_sfr[0] = l_sfr[1] = l_sfr[2] = 0;
   l_sfvr   = 0;
   l_mr     = 0;
   l_K      = 0;

   for (int i = 0; i < nlocal; i++) {
     if (scalingmask[i]) {

        // check if reduced mass is used

        mi    = (rmass) ? rmass[i] : mass[type[i]];

        // accumulate total mass

        l_mr += mi;

        // accumulate kinetic energy

        l_K  += mi/2. * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);

        // sum_j f_j * v_j

        l_sfvr  += f[i][0]*v[i][0] + f[i][1]*v[i][1] + f[i][2]*v[i][2];

        // accumulate com velocity and sum of forces

        for (int k=0; k<3; k++) {
          l_vr[k] += mi * v[i][k];
          l_sfr[k]+= f[i][k];
        }
     }
   }

   // reduce sums

   l_buf[0] = l_vr[0];
   l_buf[1] = l_vr[1];
   l_buf[2] = l_vr[2];
   l_buf[3] = l_K;
   l_buf[4] = l_mr;
   l_buf[5] = l_sfr[0];
   l_buf[6] = l_sfr[1];
   l_buf[7] = l_sfr[2];
   l_buf[8] = l_sfvr;

   MPI_Allreduce(l_buf, buf, 9, MPI_DOUBLE, MPI_SUM, world);

   // total mass of region

   *mr = buf[4];

   if (*mr < 1.e-14) {
      error->all(FLERR, "Fix ehex error mass of region is close to zero");
   }

   // total kinetic energy of region

   *K  = buf[3];

   // centre of mass velocity of region

   vr[0] = buf[0]/(*mr);
   vr[1] = buf[1]/(*mr);
   vr[2] = buf[2]/(*mr);

   // sum of forces

   sfr[0] = buf[5];
   sfr[1] = buf[6];
   sfr[2] = buf[7];

   // calculate non-translational kinetic energy

   *Kr = *K - 0.5* (*mr) * (vr[0]*vr[0]+vr[1]*vr[1]+vr[2]*vr[2]);

   // calculate sum_j f_j * (v_j - v_r) = sum_j f_j * v_j  - v_r * sum_j f_j

   *sfvr =  buf[8] - (vr[0]*sfr[0] + vr[1]*sfr[1] + vr[2]*sfr[2]);
}

