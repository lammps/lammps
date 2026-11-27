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

#include "region_intersect.h"

#include "domain.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegIntersect::RegIntersect(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), idsub(nullptr), contact_indx(nullptr)
{
  nregion = 0;

  if (narg < 5) utils::missing_cmd_args(FLERR, "region intersect", error);
  int n = utils::inumeric(FLERR, arg[2], false, lmp);
  if (n < 2) error->all(FLERR, "Illegal region intersect n: {}", n);
  options(narg - (n + 3), &arg[n + 3]);

  // build list of regions to intersect
  // store sub-region IDs in idsub

  idsub = new char *[n];
  reglist = new Region *[n];
  nregion = 0;

  for (int iarg = 0; iarg < n; iarg++) {
    idsub[nregion] = utils::strdup(arg[iarg + 3]);
    reglist[nregion] = domain->get_region_by_id(idsub[nregion]);
    if (!reglist[nregion])
      error->all(FLERR, "Region intersect region {} does not exist", idsub[nregion]);
    nregion++;
  }

  // this region is variable shape or dynamic if any of sub-regions are

  for (int ilist = 0; ilist < nregion; ilist++) {
    if (reglist[ilist]->varshape) varshape = 1;
    if (reglist[ilist]->dynamic) dynamic = 1;
    if (reglist[ilist]->moveflag) moveflag = 1;
    if (reglist[ilist]->rotateflag) rotateflag = 1;
  }

  // extent of intersection of regions
  // has bounding box if interior and any sub-region has bounding box

  bboxflag = 0;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (reglist[ilist]->bboxflag == 1) bboxflag = 1;
  if (!interior) bboxflag = 0;

  if (bboxflag) {
    int first = 1;
    for (int ilist = 0; ilist < nregion; ilist++) {
      if (reglist[ilist]->bboxflag == 0) continue;
      if (first) {
        extent_xlo = reglist[ilist]->extent_xlo;
        extent_ylo = reglist[ilist]->extent_ylo;
        extent_zlo = reglist[ilist]->extent_zlo;
        extent_xhi = reglist[ilist]->extent_xhi;
        extent_yhi = reglist[ilist]->extent_yhi;
        extent_zhi = reglist[ilist]->extent_zhi;
        first = 0;
      }

      extent_xlo = MAX(extent_xlo, reglist[ilist]->extent_xlo);
      extent_ylo = MAX(extent_ylo, reglist[ilist]->extent_ylo);
      extent_zlo = MAX(extent_zlo, reglist[ilist]->extent_zlo);
      extent_xhi = MIN(extent_xhi, reglist[ilist]->extent_xhi);
      extent_yhi = MIN(extent_yhi, reglist[ilist]->extent_yhi);
      extent_zhi = MIN(extent_zhi, reglist[ilist]->extent_zhi);
    }
  }

  // possible contacts = sum of possible contacts in all sub-regions
  // for near contacts and touching contacts

  cmax = 0;
  for (int ilist = 0; ilist < nregion; ilist++) cmax += reglist[ilist]->cmax;
  contact = new Contact[cmax];
  contact_indx = new ContactIndx[cmax];

  tmax = 0;
  for (int ilist = 0; ilist < nregion; ilist++) {
    if (interior)
      tmax += reglist[ilist]->tmax;
    else
      tmax++;
  }
}

/* ---------------------------------------------------------------------- */

RegIntersect::~RegIntersect()
{
  for (int ilist = 0; ilist < nregion; ilist++) delete[] idsub[ilist];
  delete[] idsub;
  delete[] reglist;
  delete[] contact;
  delete[] contact_indx;
}

/* ---------------------------------------------------------------------- */

void RegIntersect::init()
{
  Region::init();

  // re-build list of sub-regions in case regions were changed
  // error if a sub-region was deleted

  for (int ilist = 0; ilist < nregion; ilist++) {
    reglist[ilist] = domain->get_region_by_id(idsub[ilist]);
    if (!reglist[ilist])
      error->all(FLERR, "Region intersect region {} does not exist", idsub[ilist]);
  }

  // init the sub-regions

  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->init();
}

/* ----------------------------------------------------------------------
   overwrite parent to call for subregions
------------------------------------------------------------------------- */

void RegIntersect::prematch()
{
  int ilist;
  for (ilist = 0; ilist < nregion; ilist++) reglist[ilist]->prematch();
}

/* ----------------------------------------------------------------------
   overwrite parent to skip dynamic methods (may not be defined for union)
     will be performed by subregion match calls from inside()
------------------------------------------------------------------------- */

int RegIntersect::match(double x, double y, double z)
{
  return !(inside(x, y, z) ^ interior);
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with all sub-regions
   else inside = 0
------------------------------------------------------------------------- */

int RegIntersect::inside(double x, double y, double z)
{
  int ilist;
  for (ilist = 0; ilist < nregion; ilist++)
    if (!reglist[ilist]->match(x, y, z)) break;

  if (ilist == nregion) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   overwrite parent to skip dynamic methods (may not be defined for union)
     will be performed by surface interior/exterior calls separately
     for each subregion
------------------------------------------------------------------------- */

int RegIntersect::surface(double x, double y, double z, double cutoff)
{
  int ncontact;
  double xorig[3] = {x, y, z};

  if (!openflag) {
    if (interior)
      ncontact = surface_interior(xorig, cutoff);
    else
      ncontact = surface_exterior(xorig, cutoff);
  } else {
    // one of surface_int/ext() will return 0
    // so no need to worry about offset of contact indices
    ncontact = surface_exterior(xorig, cutoff) + surface_interior(xorig, cutoff);
  }

  return ncontact;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) compute contacts in each sub-region
   (2) only keep a contact if surface point is match() to all other regions
------------------------------------------------------------------------- */

int RegIntersect::surface_interior(double *xorig, double cutoff)
{
  int m, ilist, jlist, ncontacts;
  double xs, ys, zs;
  double xnear[3];

  int n = 0;
  int walloffset = 0;
  for (ilist = 0; ilist < nregion; ilist++) {
    auto *region = reglist[ilist];

    // copy coordinates b/c values will change
    xnear[0] = xorig[0];
    xnear[1] = xorig[1];
    xnear[2] = xorig[2];

    ncontacts = region->surface(xnear[0], xnear[1], xnear[2], cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = xnear[0] - region->contact[m].delx;
      ys = xnear[1] - region->contact[m].dely;
      zs = xnear[2] - region->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
        if (jlist == ilist) continue;
        if (!reglist[jlist]->match(xs, ys, zs)) break;
      }
      if (jlist == nregion) {
        contact[n].r = region->contact[m].r;
        contact[n].radius = region->contact[m].radius;
        contact[n].delx = region->contact[m].delx;
        contact[n].dely = region->contact[m].dely;
        contact[n].delz = region->contact[m].delz;
        contact[n].iwall = region->contact[m].iwall + walloffset;
        contact[n].varflag = region->contact[m].varflag;
        contact_indx[n].ilist = ilist;
        contact_indx[n].ic = m;
        n++;
      }
    }

    if (region->rotateflag && ncontacts) {
      for (int i = n - ncontacts; i < n; i++) {
        xs = xnear[0] - contact[i].delx;
        ys = xnear[1] - contact[i].dely;
        zs = xnear[2] - contact[i].delz;
        region->forward_transform(xs, ys, zs);
        contact[i].delx = xorig[0] - xs;
        contact[i].dely = xorig[1] - ys;
        contact[i].delz = xorig[2] - zs;
      }
    }

    // increment by cmax instead of tmax to ensure
    // possible wall IDs for sub-regions are non overlapping
    walloffset += region->cmax;
  }

  return n;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) flip interior/exterior flag of each sub-region
   (2) compute contacts in each sub-region
   (3) only keep a contact if surface point is not match() to all other regions
   (4) flip interior/exterior flags back to original settings
   this is effectively same algorithm as surface_interior() for RegUnion
------------------------------------------------------------------------- */

int RegIntersect::surface_exterior(double *xorig, double cutoff)
{
  int m, ilist, jlist, ncontacts;
  double xs, ys, zs;
  double xnear[3];

  int n = 0;
  for (ilist = 0; ilist < nregion; ilist++) reglist[ilist]->interior ^= 1;

  for (ilist = 0; ilist < nregion; ilist++) {
    auto *region = reglist[ilist];

    // copy coordinates b/c values will change
    xnear[0] = xorig[0];
    xnear[1] = xorig[1];
    xnear[2] = xorig[2];

    ncontacts = region->surface(xnear[0], xnear[1], xnear[2], cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = xnear[0] - region->contact[m].delx;
      ys = xnear[1] - region->contact[m].dely;
      zs = xnear[2] - region->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
        if (jlist == ilist) continue;
        if (reglist[jlist]->match(xs, ys, zs)) break;
      }
      if (jlist == nregion) {
        contact[n].r = region->contact[m].r;
        contact[n].radius = region->contact[m].radius;
        contact[n].delx = region->contact[m].delx;
        contact[n].dely = region->contact[m].dely;
        contact[n].delz = region->contact[m].delz;
        contact[n].iwall = ilist;
        contact[n].varflag = region->contact[m].varflag;
        contact_indx[n].ilist = ilist;
        contact_indx[n].ic = m;
        n++;
      }
    }

    if (region->rotateflag && ncontacts) {
      for (int i = n - ncontacts; i < n; i++) {
        xs = xnear[0] - contact[i].delx;
        ys = xnear[1] - contact[i].dely;
        zs = xnear[2] - contact[i].delz;
        region->forward_transform(xs, ys, zs);
        contact[i].delx = xorig[0] - xs;
        contact[i].dely = xorig[1] - ys;
        contact[i].delz = xorig[2] - zs;
      }
    }
  }

  for (ilist = 0; ilist < nregion; ilist++) reglist[ilist]->interior ^= 1;

  return n;
}

/* ----------------------------------------------------------------------
   change region shape of all sub-regions
------------------------------------------------------------------------- */

void RegIntersect::shape_update()
{
  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->shape_update();
}

/* ----------------------------------------------------------------------
   move/rotate all sub-regions
------------------------------------------------------------------------- */

void RegIntersect::pretransform()
{
  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->pretransform();
}

/* ----------------------------------------------------------------------
   get translational/angular velocities of all subregions
------------------------------------------------------------------------- */

void RegIntersect::set_velocity()
{
  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->set_velocity();
}

/* ----------------------------------------------------------------------
   call subregion method to compute velocity of wall for given contact
------------------------------------------------------------------------- */

void RegIntersect::velocity_contact(double *vwall, double *x, int ic)
{
  auto *region = reglist[contact_indx[ic].ilist];
  region->velocity_contact(vwall, x, contact_indx[ic].ic);
}

/* ----------------------------------------------------------------------
   increment length of restart buffer based on region info
   used by restart of fix/wall/gran/region
------------------------------------------------------------------------- */

void RegIntersect::length_restart_string(int &n)
{
  n += sizeof(int) + strlen(id) + 1 + sizeof(int) + strlen(style) + 1 + sizeof(int);
  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->length_restart_string(n);
}
/* ----------------------------------------------------------------------
   region writes its current position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

void RegIntersect::write_restart(FILE *fp)
{
  int sizeid = (strlen(id) + 1);
  int sizestyle = (strlen(style) + 1);
  fwrite(&sizeid, sizeof(int), 1, fp);
  fwrite(id, 1, sizeid, fp);
  fwrite(&sizestyle, sizeof(int), 1, fp);
  fwrite(style, 1, sizestyle, fp);
  fwrite(&nregion, sizeof(int), 1, fp);

  for (int ilist = 0; ilist < nregion; ilist++) { reglist[ilist]->write_restart(fp); }
}

/* ----------------------------------------------------------------------
   region reads its previous position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

int RegIntersect::restart(char *buf, int &n)
{
  int size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n], id) != 0)) return 0;
  n += size;

  size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n], style) != 0)) return 0;
  n += size;

  int restart_nreg = *((int *) (&buf[n]));
  n += sizeof(int);
  if (restart_nreg != nregion) return 0;

  for (int ilist = 0; ilist < nregion; ilist++)
    if (!reglist[ilist]->restart(buf, n)) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   set prev vector to zero
------------------------------------------------------------------------- */

void RegIntersect::reset_vel()
{
  for (int ilist = 0; ilist < nregion; ilist++) reglist[ilist]->reset_vel();
}
