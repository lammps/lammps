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

#include "region_union.h"

#include "domain.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegUnion::RegUnion(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg),
  idsub(nullptr)
{
  nregion = 0;
  if (narg < 5) error->all(FLERR,"Illegal region command");
  int n = utils::inumeric(FLERR,arg[2],false,lmp);
  if (n < 2) error->all(FLERR,"Illegal region command");
  options(narg-(n+3),&arg[n+3]);

  // build list of region indices to union
  // store sub-region IDs in idsub

  idsub = new char*[n];
  list = new int[n];
  nregion = 0;

  for (int iarg = 0; iarg < n; iarg++) {
    idsub[nregion] = utils::strdup(arg[iarg+3]);
    int iregion = domain->find_region(idsub[nregion]);
    if (iregion == -1)
      error->all(FLERR,"Region union region ID does not exist");
    list[nregion++] = iregion;
  }

  // this region is variable shape or dynamic if any of sub-regions are

  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++) {
    if (regions[list[ilist]]->varshape) varshape = 1;
    if (regions[list[ilist]]->dynamic) dynamic = 1;
    if (regions[list[ilist]]->moveflag) moveflag = 1;
    if (regions[list[ilist]]->rotateflag) rotateflag = 1;
  }

  // extent of union of regions
  // has bounding box if interior and all sub-regions have bounding box

  bboxflag = 1;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->bboxflag == 0) bboxflag = 0;
  if (!interior) bboxflag = 0;

  if (bboxflag) {
    extent_xlo = extent_ylo = extent_zlo = BIG;
    extent_xhi = extent_yhi = extent_zhi = -BIG;

    for (int ilist = 0; ilist < nregion; ilist++) {
      extent_xlo = MIN(extent_xlo,regions[list[ilist]]->extent_xlo);
      extent_ylo = MIN(extent_ylo,regions[list[ilist]]->extent_ylo);
      extent_zlo = MIN(extent_zlo,regions[list[ilist]]->extent_zlo);
      extent_xhi = MAX(extent_xhi,regions[list[ilist]]->extent_xhi);
      extent_yhi = MAX(extent_yhi,regions[list[ilist]]->extent_yhi);
      extent_zhi = MAX(extent_zhi,regions[list[ilist]]->extent_zhi);
    }
  }

  // possible contacts = sum of possible contacts in all sub-regions
  // for near contacts and touching contacts

  cmax = 0;
  for (int ilist = 0; ilist < nregion; ilist++)
    cmax += regions[list[ilist]]->cmax;
  contact = new Contact[cmax];

  tmax = 0;
  for (int ilist = 0; ilist < nregion; ilist++) {
    if (interior) tmax += regions[list[ilist]]->tmax;
    else tmax++;
  }
}

/* ---------------------------------------------------------------------- */

RegUnion::~RegUnion()
{
  for (int ilist = 0; ilist < nregion; ilist++) delete [] idsub[ilist];
  delete [] idsub;
  delete [] list;
  delete [] contact;
}

/* ---------------------------------------------------------------------- */

void RegUnion::init()
{
  Region::init();

  // re-build list of sub-regions in case other regions were deleted
  // error if a sub-region was deleted

  for (int ilist = 0; ilist < nregion; ilist++) {
    int iregion = domain->find_region(idsub[ilist]);
    if (iregion == -1)
      error->all(FLERR,"Region union region ID does not exist");
    list[ilist] = iregion;
  }

  // init the sub-regions

  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->init();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with any sub-region
   else inside = 0
------------------------------------------------------------------------- */

int RegUnion::inside(double x, double y, double z)
{
  int ilist;
  Region **regions = domain->regions;
  for (ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->match(x,y,z)) break;

  if (ilist == nregion) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of union of sub-regions
   (1) compute contacts in each sub-region
   (2) only keep a contact if surface point is not match() to all other regions
------------------------------------------------------------------------- */

int RegUnion::surface_interior(double *x, double cutoff)
{
  int m,ilist,jlist,iregion,jregion,ncontacts;
  double xs,ys,zs;

  Region **regions = domain->regions;
  int n = 0;

  int walloffset = 0 ;
  for (ilist = 0; ilist < nregion; ilist++) {
    iregion = list[ilist];
    ncontacts = regions[iregion]->surface(x[0],x[1],x[2],cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = x[0] - regions[iregion]->contact[m].delx;
      ys = x[1] - regions[iregion]->contact[m].dely;
      zs = x[2] - regions[iregion]->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
        if (jlist == ilist) continue;
        jregion = list[jlist];
        if (regions[jregion]->match(xs,ys,zs) &&
            !regions[jregion]->openflag) break;
      }
      if (jlist == nregion) {
        contact[n].r = regions[iregion]->contact[m].r;
        contact[n].radius = regions[iregion]->contact[m].radius;
        contact[n].delx = regions[iregion]->contact[m].delx;
        contact[n].dely = regions[iregion]->contact[m].dely;
        contact[n].delz = regions[iregion]->contact[m].delz;
        contact[n].iwall = regions[iregion]->contact[m].iwall + walloffset;
        contact[n].varflag = regions[iregion]->contact[m].varflag;
        n++;
      }
    }
    // increment by cmax instead of tmax to insure
    // possible wall IDs for sub-regions are non overlapping
    walloffset += regions[iregion]->cmax;
  }

  return n;
}

/* ----------------------------------------------------------------------
   compute contacts with exterior of union of sub-regions
   (1) flip interior/exterior flag of each sub-region
   (2) compute contacts in each sub-region
   (3) only keep a contact if surface point is match() to all other regions
   (4) flip interior/exterior flags back to original settings
   this is effectively same algorithm as surface_interior() for RegIntersect
------------------------------------------------------------------------- */

int RegUnion::surface_exterior(double *x, double cutoff)
{
  int m,ilist,jlist,iregion,jregion,ncontacts;
  double xs,ys,zs;

  Region **regions = domain->regions;
  int n = 0;

  for (ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->interior ^= 1;

  for (ilist = 0; ilist < nregion; ilist++) {
    iregion = list[ilist];
    ncontacts = regions[iregion]->surface(x[0],x[1],x[2],cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = x[0] - regions[iregion]->contact[m].delx;
      ys = x[1] - regions[iregion]->contact[m].dely;
      zs = x[2] - regions[iregion]->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
        if (jlist == ilist) continue;
        jregion = list[jlist];
        if (!regions[jregion]->match(xs,ys,zs)) break;
      }
      if (jlist == nregion) {
        contact[n].r = regions[iregion]->contact[m].r;
        contact[n].radius = regions[iregion]->contact[m].radius;
        contact[n].delx = regions[iregion]->contact[m].delx;
        contact[n].dely = regions[iregion]->contact[m].dely;
        contact[n].delz = regions[iregion]->contact[m].delz;
        contact[n].iwall = ilist;
        contact[n].r = regions[iregion]->contact[m].varflag;
        n++;
      }
    }
  }

  for (ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->interior ^= 1;

  return n;
}

/* ----------------------------------------------------------------------
   change region shape of all sub-regions
------------------------------------------------------------------------- */

void RegUnion::shape_update()
{
  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->shape_update();
}

/* ----------------------------------------------------------------------
   move/rotate all sub-regions
------------------------------------------------------------------------- */

void RegUnion::pretransform()
{
  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->pretransform();
}

/* ----------------------------------------------------------------------
   get translational/angular velocities of all subregions
------------------------------------------------------------------------- */

void RegUnion::set_velocity()
{
  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->set_velocity();
}

/* ----------------------------------------------------------------------
   increment length of restart buffer based on region info
   used by restart of fix/wall/gran/region
------------------------------------------------------------------------- */

void RegUnion::length_restart_string(int& n)
{
  n += sizeof(int) + strlen(id)+1 +
    sizeof(int) + strlen(style)+1 + sizeof(int);
  for (int ilist = 0; ilist < nregion; ilist++)
    domain->regions[list[ilist]]->length_restart_string(n);

}
/* ----------------------------------------------------------------------
   region writes its current position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

void RegUnion::write_restart(FILE *fp)
{
  int sizeid = (strlen(id)+1);
  int sizestyle = (strlen(style)+1);
  fwrite(&sizeid, sizeof(int), 1, fp);
  fwrite(id, 1, sizeid, fp);
  fwrite(&sizestyle, sizeof(int), 1, fp);
  fwrite(style, 1, sizestyle, fp);
  fwrite(&nregion,sizeof(int),1,fp);
  for (int ilist = 0; ilist < nregion; ilist++)
    domain->regions[list[ilist]]->write_restart(fp);
}

/* ----------------------------------------------------------------------
   region reads its previous position/angle
   needed by fix/wall/gran/region to compute velocity by differencing scheme
------------------------------------------------------------------------- */

int RegUnion::restart(char *buf, int &n)
{
  int size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n],id) != 0)) return 0;
  n += size;

  size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n],style) != 0)) return 0;
  n += size;

  int restart_nreg = *((int *) (&buf[n]));
  n += sizeof(int);
  if (restart_nreg != nregion) return 0;

  for (int ilist = 0; ilist < nregion; ilist++)
    if (!domain->regions[list[ilist]]->restart(buf,n)) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   set prev vector to zero
------------------------------------------------------------------------- */

void RegUnion::reset_vel()
{
  for (int ilist = 0; ilist < nregion; ilist++)
    domain->regions[list[ilist]]->reset_vel();
}
