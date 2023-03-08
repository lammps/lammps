/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Dan Bolintineanu (SNL)
------------------------------------------------------------------------- */

#include "fix_wall_gran_region.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "neighbor.h"
#include "region.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallGranRegion::FixWallGranRegion(LAMMPS *lmp, int narg, char **arg) :
    FixWallGran(lmp, narg, arg), region(nullptr), ncontact(nullptr), walls(nullptr),
    history_many(nullptr), c2r(nullptr)
{
  restart_global = 1;
  motion_resetflag = 0;

  region = domain->get_region_by_id(idregion);
  if (!region) error->all(FLERR, "Region {} for fix wall/gran/region does not exist", idregion);
  nregion = region->nregion;
  tmax = region->tmax;
  c2r = new int[tmax];

  // re-allocate atom-based arrays with nshear
  // do not register with Atom class, since parent class did that

  memory->destroy(history_one);
  history_one = nullptr;

  ncontact = nullptr;
  walls = nullptr;
  history_many = nullptr;
  FixWallGranRegion::grow_arrays(atom->nmax);

  // initialize shear history as if particle is not touching region

  if (use_history) {
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) ncontact[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

FixWallGranRegion::~FixWallGranRegion()
{
  delete[] c2r;

  memory->destroy(ncontact);
  memory->destroy(walls);
  memory->destroy(history_many);
}

/* ---------------------------------------------------------------------- */

void FixWallGranRegion::init()
{
  FixWallGran::init();

  auto newregion = domain->get_region_by_id(idregion);
  if (!newregion) error->all(FLERR, "Region {} for fix wall/gran/region does not exist", idregion);

  // check if region properties changed between runs
  // reset if restart info was inconsistent

  if (newregion != region) {
    region = newregion;
    if (comm->me == 0)
      error->warning(FLERR,
                     "Region properties for region {} changed between runs, resetting its motion",
                     idregion);
    nregion = region->nregion;
    tmax = region->tmax;
    delete[] c2r;
    c2r = new int[tmax];
    region = newregion;
    region->reset_vel();
  }

  if (motion_resetflag) {
    if (comm->me == 0)
      error->warning(FLERR,
                     "Region properties for region {} are inconsistent with restart file, "
                     "resetting its motion",
                     idregion);
    region->reset_vel();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranRegion::post_force(int /*vflag*/)
{
  int i, m, nc, iwall;
  double dx, dy, dz, rsq, meff;
  double vwall[3];

  // do not update shear history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body", tmp);
    auto mass_body = (double *) fix_rigid->extract("masstotal", tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid, nmax, "wall/gran:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0)
        mass_rigid[i] = mass_body[body[i]];
      else
        mass_rigid[i] = 0.0;
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

  if (peratom_flag) { clear_stored_contacts(); }

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (!region->match(x[i][0], x[i][1], x[i][2])) continue;

      if (pairstyle == FixWallGran::GRANULAR && normal_model == FixWallGran::JKR) {
        nc = region->surface(x[i][0], x[i][1], x[i][2], radius[i] + pulloff_distance(radius[i]));
      } else {
        nc = region->surface(x[i][0], x[i][1], x[i][2], radius[i]);
      }
      if (nc > tmax) error->one(FLERR, "Too many wall/gran/region contacts for one particle");

      // shear history maintenance
      // update ncontact,walls,shear2many for particle I
      //   to reflect new and persistent shear history values
      // also set c2r[] = indices into region->contact[] for each of N contacts
      // process zero or one contact here, otherwise invoke update_contacts()

      if (use_history) {
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
            for (m = 0; m < size_history; m++) history_many[i][0][m] = 0.0;
          } else if (ncontact[i] > 1 || iwall != walls[i][0])
            update_contacts(i, nc);
        } else
          update_contacts(i, nc);
      }

      // process current contacts
      for (int ic = 0; ic < nc; ic++) {

        // rsq = squared contact distance
        // xc = contact point

        rsq = region->contact[ic].r * region->contact[ic].r;

        if (pairstyle == FixWallGran::GRANULAR && normal_model == FixWallGran::JKR) {
          if (history_many[i][c2r[ic]][0] == 0.0 && rsq > radius[i] * radius[i]) {
            for (m = 0; m < size_history; m++) history_many[i][0][m] = 0.0;
            continue;
          }
        }

        dx = region->contact[ic].delx;
        dy = region->contact[ic].dely;
        dz = region->contact[ic].delz;

        if (regiondynamic) region->velocity_contact(vwall, x[i], ic);

        // meff = effective mass of sphere
        // if I is part of rigid body, use body mass

        meff = rmass[i];
        if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

        // store contact info
        if (peratom_flag) {
          array_atom[i][0] = 1.0;
          array_atom[i][4] = x[i][0] - dx;
          array_atom[i][5] = x[i][1] - dy;
          array_atom[i][6] = x[i][2] - dz;
          array_atom[i][7] = radius[i];
        }

        // invoke sphere/wall interaction
        double *contact;
        if (peratom_flag)
          contact = array_atom[i];
        else
          contact = nullptr;

        if (pairstyle == FixWallGran::HOOKE)
          hooke(rsq, dx, dy, dz, vwall, v[i], f[i], omega[i], torque[i], radius[i], meff, contact);
        else if (pairstyle == FixWallGran::HOOKE_HISTORY)
          hooke_history(rsq, dx, dy, dz, vwall, v[i], f[i], omega[i], torque[i], radius[i], meff,
                        history_many[i][c2r[ic]], contact);
        else if (pairstyle == FixWallGran::HERTZ_HISTORY)
          hertz_history(rsq, dx, dy, dz, vwall, region->contact[ic].radius, v[i], f[i], omega[i],
                        torque[i], radius[i], meff, history_many[i][c2r[ic]], contact);
        else if (pairstyle == FixWallGran::GRANULAR)
          granular(rsq, dx, dy, dz, vwall, region->contact[ic].radius, v[i], f[i], omega[i],
                   torque[i], radius[i], meff, history_many[i][c2r[ic]], contact);
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
  int j, m, iold, nold, ilast, inew, iadd, iwall;

  // loop over old contacts
  // if not in new contact list:
  //   delete old contact by copying last contact over it

  iold = 0;
  while (iold < ncontact[i]) {
    for (m = 0; m < nc; m++)
      if (region->contact[m].iwall == walls[i][iold]) break;
    if (m >= nc) {
      ilast = ncontact[i] - 1;
      for (j = 0; j < size_history; j++) history_many[i][iold][j] = history_many[i][ilast][j];
      walls[i][iold] = walls[i][ilast];
      ncontact[i]--;
    } else
      iold++;
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
    if (m < nold)
      c2r[m] = inew;
    else {
      iadd = ncontact[i];

      c2r[iadd] = inew;
      for (j = 0; j < size_history; j++) history_many[i][iadd][j] = 0.0;
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
  if (use_history) {                                                  // shear history
    bytes += (double) nmax * sizeof(int);                             // ncontact
    bytes += (double) nmax * tmax * sizeof(int);                      // walls
    bytes += (double) nmax * tmax * size_history * sizeof(double);    // history_many
  }
  if (fix_rigid) bytes += (double) nmax * sizeof(int);    // mass_rigid
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranRegion::grow_arrays(int nmax)
{
  if (use_history) {
    memory->grow(ncontact, nmax, "fix_wall_gran:ncontact");
    memory->grow(walls, nmax, tmax, "fix_wall_gran:walls");
    memory->grow(history_many, nmax, tmax, size_history, "fix_wall_gran:history_many");
  }
  if (peratom_flag) memory->grow(array_atom, nmax, size_peratom_cols, "fix_wall_gran:array_atom");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranRegion::copy_arrays(int i, int j, int /*delflag*/)
{
  int m, n, iwall;

  if (use_history) {
    n = ncontact[i];
    for (iwall = 0; iwall < n; iwall++) {
      walls[j][iwall] = walls[i][iwall];
      for (m = 0; m < size_history; m++) history_many[j][iwall][m] = history_many[i][iwall][m];
    }
    ncontact[j] = ncontact[i];
  }

  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) array_atom[j][m] = array_atom[i][m];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGranRegion::set_arrays(int i)
{
  if (use_history) ncontact[i] = 0;
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) array_atom[i][m] = 0;
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGranRegion::pack_exchange(int i, double *buf)
{
  int m;

  int n = 0;
  if (use_history) {
    int count = ncontact[i];
    buf[n++] = ubuf(count).d;
    for (int iwall = 0; iwall < count; iwall++) {
      buf[n++] = ubuf(walls[i][iwall]).d;
      for (m = 0; m < size_history; m++) buf[n++] = history_many[i][iwall][m];
    }
  }
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) buf[n++] = array_atom[i][m];
  }

  return n;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGranRegion::unpack_exchange(int nlocal, double *buf)
{
  int m;

  int n = 0;
  if (use_history) {
    int count = ncontact[nlocal] = (int) ubuf(buf[n++]).i;
    for (int iwall = 0; iwall < count; iwall++) {
      walls[nlocal][iwall] = (int) ubuf(buf[n++]).i;
      for (m = 0; m < size_history; m++) history_many[nlocal][iwall][m] = buf[n++];
    }
  }
  if (peratom_flag) {
    for (int m = 0; m < size_peratom_cols; m++) array_atom[nlocal][m] = buf[n++];
  }

  return n;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGranRegion::pack_restart(int i, double *buf)
{
  int m;

  if (!use_history) return 0;

  int n = 1;
  int count = ncontact[i];

  buf[n++] = ubuf(count).d;
  for (int iwall = 0; iwall < count; iwall++) {
    buf[n++] = ubuf(walls[i][iwall]).d;
    for (m = 0; m < size_history; m++) buf[n++] = history_many[i][iwall][m];
  }
  // pack buf[0] this way because other fixes unpack it
  buf[0] = n;
  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGranRegion::unpack_restart(int nlocal, int nth)
{
  int k;

  if (!use_history) return;

  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  int count = ncontact[nlocal] = (int) ubuf(extra[nlocal][m++]).i;
  for (int iwall = 0; iwall < count; iwall++) {
    walls[nlocal][iwall] = (int) ubuf(extra[nlocal][m++]).i;
    for (k = 0; k < size_history; k++) history_many[nlocal][iwall][k] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGranRegion::maxsize_restart()
{
  if (!use_history) return 0;
  return 2 + tmax * (size_history + 1);
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGranRegion::size_restart(int nlocal)
{
  if (!use_history) return 0;
  return 2 + ncontact[nlocal] * (size_history + 1);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixWallGranRegion::write_restart(FILE *fp)
{
  if (comm->me) return;
  int len = 0;
  region->length_restart_string(len);
  fwrite(&len, sizeof(int), 1, fp);
  region->write_restart(fp);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixWallGranRegion::restart(char *buf)
{
  int n = 0;
  if (!region->restart(buf, n)) motion_resetflag = 1;
}
