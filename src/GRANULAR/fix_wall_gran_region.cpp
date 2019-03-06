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
   Contributing authors: Dan Bolintineanu (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_wall_gran_region.h"
#include "region.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// same as FixWallGran

enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,BONDED_HISTORY};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallGranRegion::FixWallGranRegion(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg), region(NULL), region_style(NULL), ncontact(NULL),
  walls(NULL), shearmany(NULL), c2r(NULL)
{
  restart_global = 1;
  motion_resetflag = 0;

  int iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/gran/region does not exist");
  region = domain->regions[iregion];
  region_style = new char[strlen(region->style)+1];
  strcpy(region_style,region->style);
  nregion = region->nregion;

  tmax = domain->regions[iregion]->tmax;
  c2r = new int[tmax];

  // re-allocate atom-based arrays with nshear
  // do not register with Atom class, since parent class did that

  memory->destroy(shearone);
  shearone = NULL;

  ncontact = NULL;
  walls = NULL;
  shearmany = NULL;
  grow_arrays(atom->nmax);

  // initialize shear history as if particle is not touching region

  if (history) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      ncontact[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

FixWallGranRegion::~FixWallGranRegion()
{
  delete [] c2r;
  delete [] region_style;

  memory->destroy(ncontact);
  memory->destroy(walls);
  memory->destroy(shearmany);
}

/* ---------------------------------------------------------------------- */

void FixWallGranRegion::init()
{
  FixWallGran::init();

  int iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix wall/gran/region does not exist");
  region = domain->regions[iregion];

  // check if region properties changed between runs
  // reset if restart info was inconsistent

  if (strcmp(idregion,region->id) != 0 ||
      strcmp(region_style,region->style) != 0 ||
      nregion != region->nregion) {
    char str[256];
    snprintf(str,256,"Region properties for region %s changed between runs, "
             "resetting its motion",idregion);
    error->warning(FLERR,str);
    region->reset_vel();
  }

  if (motion_resetflag){
    char str[256];
    snprintf(str,256,"Region properties for region %s are inconsistent "
             "with restart file, resetting its motion",idregion);
    error->warning(FLERR,str);
    region->reset_vel();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranRegion::post_force(int /*vflag*/)
{
  int i,m,nc,iwall;
  double dx,dy,dz,rsq,meff;
  double vwall[3];

  // do not update shear history during setup

  shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"wall/gran:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    }
  }

  int regiondynamic = region->dynamic_check();
  if (!regiondynamic) vwall[0] = vwall[1] = vwall[2] = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // set current motion attributes of region
  // set_velocity() also updates prev to current step

  if (regiondynamic) {
    region->prematch();
    region->set_velocity();
  }

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (!region->match(x[i][0],x[i][1],x[i][2])) continue;

      nc = region->surface(x[i][0],x[i][1],x[i][2],radius[i]);
      if (nc > tmax)
        error->one(FLERR,"Too many wall/gran/region contacts for one particle");

      // shear history maintenance
      // update ncontact,walls,shear2many for particle I
      //   to reflect new and persistent shear history values
      // also set c2r[] = indices into region->contact[] for each of N contacts
      // process zero or one contact here, otherwise invoke update_contacts()

      if (history) {
        if (nc == 0) {
          ncontact[i] = 0;
          continue;
        }
        if (nc == 1) {
          c2r[0] = 0;
          iwall = region->contact[0].iwall;
          if (ncontact[i] == 0) {
            ncontact[i] = 1;
            walls[i][0] = iwall;
            for (m = 0; m < sheardim; m++)
              shearmany[i][0][m] = 0.0;
          } else if (ncontact[i] > 1 || iwall != walls[i][0])
            update_contacts(i,nc);
        } else update_contacts(i,nc);
      }

      // process current contacts

      for (int ic = 0; ic < nc; ic++) {

        // rsq = squared contact distance
        // xc = contact point

        rsq = region->contact[ic].r*region->contact[ic].r;

        dx = region->contact[ic].delx;
        dy = region->contact[ic].dely;
        dz = region->contact[ic].delz;

        if (regiondynamic) region->velocity_contact(vwall, x[i], ic);


        // meff = effective mass of sphere
        // if I is part of rigid body, use body mass

        meff = rmass[i];
        if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

        // invoke sphere/wall interaction

        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,vwall,v[i],f[i],
                omega[i],torque[i],radius[i],meff);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,vwall,v[i],f[i],
                        omega[i],torque[i],radius[i],meff,
                        shearmany[i][c2r[ic]]);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,vwall,region->contact[ic].radius,
                        v[i],f[i],omega[i],torque[i],
                        radius[i],meff,shearmany[i][c2r[ic]]);
        else if (pairstyle == BONDED_HISTORY)
          bonded_history(rsq,dx,dy,dz,vwall,region->contact[ic].radius,
                         v[i],f[i],omega[i],torque[i],
                         radius[i],meff,shearmany[i][c2r[ic]]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   update contact info in ncontact, walls, shear2many for particle I
   based on ncontacts[i] old contacts and N new contacts
     matched via their associated walls
   delete/zero shear history for broken/new contacts
   also set c2r[i] = index of Ith contact in region list of contacts
------------------------------------------------------------------------- */

void FixWallGranRegion::update_contacts(int i, int nc)
{
  int j,m,iold,nold,ilast,inew,iadd,iwall;

  // loop over old contacts
  // if not in new contact list:
  //   delete old contact by copying last contact over it

  iold = 0;
  while (iold < ncontact[i]) {
    for (m = 0; m < nc; m++)
      if (region->contact[m].iwall == walls[i][iold]) break;
    if (m >= nc) {
      ilast = ncontact[i]-1;
      for (j = 0; j < sheardim; j++)
        shearmany[i][iold][j] = shearmany[i][ilast][j];
      walls[i][iold] = walls[i][ilast];
      ncontact[i]--;
    } else iold++;
  }

  // loop over new contacts
  // if not in newly compressed contact list of length nold:
  //   add it with zeroed shear history
  // set all values in c2r

  nold = ncontact[i];

  for (inew = 0; inew < nc; inew++) {
    iwall = region->contact[inew].iwall;
    for (m = 0; m < nold; m++)
      if (walls[i][m] == iwall) break;
    if (m < nold) c2r[m] = inew;
    else {
      iadd = ncontact[i];

      c2r[iadd] = inew;
      for (j = 0; j < sheardim; j++)
        shearmany[i][iadd][j] = 0.0;
      walls[i][iadd] = iwall;
      ncontact[i]++;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGranRegion::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  if (history) {                                   // shear history
    bytes += nmax * sizeof(int);                   // ncontact
    bytes += nmax*tmax * sizeof(int);              // walls
    bytes += nmax*tmax*sheardim * sizeof(double);  // shearmany
  }
  if (fix_rigid) bytes += nmax * sizeof(int);      // mass_rigid
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranRegion::grow_arrays(int nmax)
{
  if (history) {
    memory->grow(ncontact,nmax,"fix_wall_gran:ncontact");
    memory->grow(walls,nmax,tmax,"fix_wall_gran:walls");
    memory->grow(shearmany,nmax,tmax,sheardim,"fix_wall_gran:shearmany");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranRegion::copy_arrays(int i, int j, int /*delflag*/)
{
  int m,n,iwall;

  if (!history) return;

  n = ncontact[i];

  for (iwall = 0; iwall < n; iwall++) {
    walls[j][iwall] = walls[i][iwall];
    for (m = 0; m < sheardim; m++)
      shearmany[j][iwall][m] = shearmany[i][iwall][m];
  }
  ncontact[j] = ncontact[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGranRegion::set_arrays(int i)
{
  if (!history) return;
  ncontact[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGranRegion::pack_exchange(int i, double *buf)
{
  int m;

  if (!history) return 0;

  int n = 0;
  int count = ncontact[i];

  buf[n++] = ubuf(count).d;
  for (int iwall = 0; iwall < count; iwall++) {
    buf[n++] = ubuf(walls[i][iwall]).d;
    for (m = 0; m < sheardim; m++)
      buf[n++] = shearmany[i][iwall][m];
  }

  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGranRegion::unpack_exchange(int nlocal, double *buf)
{
  int m;

  if (!history) return 0;

  int n = 0;
  int count = ncontact[nlocal] = (int) ubuf(buf[n++]).i;

  for (int iwall = 0; iwall < count; iwall++) {
    walls[nlocal][iwall] = (int) ubuf(buf[n++]).i;
    for (m = 0; m < sheardim; m++)
      shearmany[nlocal][iwall][m] = buf[n++];
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGranRegion::pack_restart(int i, double *buf)
{
  int m;

  if (!history) return 0;

  int n = 1;
  int count = ncontact[i];

  buf[n++] = ubuf(count).d;
  for (int iwall = 0; iwall < count; iwall++) {
    buf[n++] = ubuf(walls[i][iwall]).d;
    for (m = 0; m < sheardim; m++)
      buf[n++] = shearmany[i][iwall][m];
  }
  buf[0] = n;
  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGranRegion::unpack_restart(int nlocal, int nth)
{
  int k;

  if (!history) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  int count = ncontact[nlocal] = (int) ubuf(extra[nlocal][m++]).i;
  for (int iwall = 0; iwall < count; iwall++) {
    walls[nlocal][iwall] = (int) ubuf(extra[nlocal][m++]).i;
    for (k = 0; k < sheardim; k++)
      shearmany[nlocal][iwall][k] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGranRegion::maxsize_restart()
{
  if (!history) return 0;
  return 2 + tmax*(sheardim+1);
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGranRegion::size_restart(int nlocal)
{
  if (!history) return 0;
  return 2 + ncontact[nlocal]*(sheardim+1);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixWallGranRegion::write_restart(FILE *fp)
{
  if (comm->me) return;
  int len = 0;
  region->length_restart_string(len);
  fwrite(&len, sizeof(int),1,fp);
  region->write_restart(fp);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixWallGranRegion::restart(char *buf)
{
  int n = 0;
  if (!region->restart(buf,n)) motion_resetflag = 1;
}
