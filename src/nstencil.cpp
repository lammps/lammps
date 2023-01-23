// clang-format off
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

#include "nstencil.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   NStencil classes
   each has method to create a stencil = list of bin offsets
     invoked each time simulation box size/shape changes
     since induces change in bins
   stencil = bins whose closest corner to central bin is within cutoff
   sx,sy,sz = bin bounds = furthest the stencil could possibly extend
     calculated below in create_setup()
   3d creates xyz stencil, 2d creates xy stencil
   for full list or half list with newton off
     use a full stencil
     stencil is all surrounding bins including self
     regardless of triclinic
   for half list with newton on
     use a half stencil
     stencil is bins to the "upper right" of central bin
     stencil does not include self
     no versions that allow ghost on (no callers need it?)
   for half list with newton on and triclinic:
     use a half stencil
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
     no versions that allow ghost on (no callers need it?)
   for multi/old:
     create one stencil for each atom type
     stencil follows same rules for half/full, newton on/off, triclinic
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
     no versions that allow ghost on (any need for it?)
   for multi:
     create one stencil for each icollection-jcollection pairing
     the full/half stencil label refers to the same-collection stencil
     a half list with newton on has a half same-collection stencil
     a full list or half list with newton off has a full same-collection stencil
     cross collection stencils are always full to allow small-to-large lookups
     for orthogonal boxes, a half stencil includes bins to the "upper right" of central bin
     for triclinic, a half stencil includes bins in the z (3D) or y (2D) plane of self and above
     cutoff is not cutneighmaxsq, but max cutoff for that atom collection
     no versions that allow ghost on (any need for it?)
------------------------------------------------------------------------- */

NStencil::NStencil(LAMMPS *lmp) : Pointers(lmp)
{
  last_stencil = -1;

  xyzflag = 0;
  maxstencil = maxstencil_multi_old = 0;
  stencil = nullptr;
  stencilxyz = nullptr;
  nstencil_multi_old = nullptr;
  stencil_multi_old = nullptr;
  distsq_multi_old = nullptr;

  nstencil_multi = nullptr;
  stencil_multi = nullptr;
  maxstencil_multi = nullptr;

  flag_half_multi = nullptr;
  flag_skip_multi = nullptr;
  bin_collection_multi = nullptr;

  maxcollections = 0;

  dimension = domain->dimension;
}

/* ---------------------------------------------------------------------- */

NStencil::~NStencil()
{
  memory->destroy(stencil);
  memory->destroy(stencilxyz);

  if (stencil_multi_old) {

    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      memory->destroy(stencil_multi_old[i]);
      memory->destroy(distsq_multi_old[i]);
    }
    delete [] nstencil_multi_old;
    delete [] stencil_multi_old;
    delete [] distsq_multi_old;
  }

  if (maxstencil_multi) {

    memory->destroy(nstencil_multi);
    for (int i = 0; i < maxcollections; i++) {
      for (int j = 0; j < maxcollections; j++)
          memory->destroy(stencil_multi[i][j]);
      delete [] stencil_multi[i];
    }
    delete [] stencil_multi;
    memory->destroy(maxstencil_multi);
    memory->destroy(flag_half_multi);
    memory->destroy(flag_skip_multi);
    memory->destroy(bin_collection_multi);

    memory->destroy(stencil_sx_multi);
    memory->destroy(stencil_sy_multi);
    memory->destroy(stencil_sz_multi);

    memory->destroy(stencil_mbinx_multi);
    memory->destroy(stencil_mbiny_multi);
    memory->destroy(stencil_mbinz_multi);

    memory->destroy(stencil_binsizex_multi);
    memory->destroy(stencil_binsizey_multi);
    memory->destroy(stencil_binsizez_multi);
  }
}

/* ---------------------------------------------------------------------- */

void NStencil::post_constructor(NeighRequest *nrq)
{
  cutoff_custom = 0.0;
  if (nrq->cut) cutoff_custom = nrq->cutoff;
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_neighbor_info()
{
  neighstyle = neighbor->style;
  cutneighmax = neighbor->cutneighmax;
  cutneighmaxsq = neighbor->cutneighmaxsq;
  cuttypesq = neighbor->cuttypesq;
  cutneighsq = neighbor->cutneighsq;

  ncollections = neighbor->ncollections;
  collection = neighbor->collection;
  cutcollectionsq = neighbor->cutcollectionsq;

  // overwrite Neighbor cutoff with custom value set by requestor
  // only works for style = BIN (checked by Neighbor class)

  if (cutoff_custom > 0.0) {
    cutneighmax = cutoff_custom;
    cutneighmaxsq = cutneighmax * cutneighmax;
  }
}

/* ----------------------------------------------------------------------
   copy needed info from NBin class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_bin_info()
{
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  binsizex = nb->binsizex;
  binsizey = nb->binsizey;
  binsizez = nb->binsizez;
  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;
}

/* ----------------------------------------------------------------------
   copy needed info for multi from NBin class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_bin_info_multi()
{
  mbinx_multi = nb->mbinx_multi;
  mbiny_multi = nb->mbiny_multi;
  mbinz_multi = nb->mbinz_multi;
  binsizex_multi = nb->binsizex_multi;
  binsizey_multi = nb->binsizey_multi;
  binsizez_multi = nb->binsizez_multi;
  bininvx_multi = nb->bininvx_multi;
  bininvy_multi = nb->bininvy_multi;
  bininvz_multi = nb->bininvz_multi;
}

/* ----------------------------------------------------------------------
   ensure NBin data is current
   ensure stencils are allocated large enough
------------------------------------------------------------------------- */

void NStencil::create_setup()
{

  if (neighstyle != Neighbor::MULTI) {
    if (nb) copy_bin_info();
    last_stencil = update->ntimestep;

    // sx,sy,sz = max range of stencil in each dim
    // smax = max possible size of entire 3d stencil
    // stencil will be empty if cutneighmax = 0.0

    sx = static_cast<int>(cutneighmax*bininvx);
    if (sx*binsizex < cutneighmax) sx++;
    sy = static_cast<int>(cutneighmax*bininvy);
    if (sy*binsizey < cutneighmax) sy++;
    sz = static_cast<int>(cutneighmax*bininvz);
    if (sz*binsizez < cutneighmax) sz++;
    if (dimension == 2) sz = 0;

    int smax = (2*sx+1) * (2*sy+1) * (2*sz+1);

    // reallocate stencil structs if necessary
    // for BIN and MULTI_OLD styles

    if (neighstyle == Neighbor::BIN) {
      if (smax > maxstencil) {
        maxstencil = smax;
        memory->destroy(stencil);
        memory->create(stencil,maxstencil,"neighstencil:stencil");
        if (xyzflag) {
          memory->destroy(stencilxyz);
          memory->create(stencilxyz,maxstencil,3,"neighstencil:stencilxyz");
        }
      }

    } else {
      int i;
      int n = atom->ntypes;
      if (maxstencil_multi_old == 0) {
        nstencil_multi_old = new int[n+1];
        stencil_multi_old = new int*[n+1];
        distsq_multi_old = new double*[n+1];
        for (i = 1; i <= n; i++) {
          nstencil_multi_old[i] = 0;
          stencil_multi_old[i] = nullptr;
          distsq_multi_old[i] = nullptr;
        }
      }
      if (smax > maxstencil_multi_old) {
        maxstencil_multi_old = smax;
        for (i = 1; i <= n; i++) {
          memory->destroy(stencil_multi_old[i]);
          memory->destroy(distsq_multi_old[i]);
          memory->create(stencil_multi_old[i],maxstencil_multi_old,
                         "neighstencil:stencil_multi_old");
          memory->create(distsq_multi_old[i],maxstencil_multi_old,
                         "neighstencil:distsq_multi_old");
        }
      }
    }
  } else {
    int i, j, bin_collection, smax;
    double stencil_range;
    int n = ncollections;

    if (nb) copy_bin_info_multi();

    // Deallocate arrays if previously allocated
    if((n > maxcollections) && stencil_multi){
      memory->destroy(nstencil_multi);
      for (i = 0; i < maxcollections; i++) {
        for (j = 0; j < maxcollections; j++)
          memory->destroy(stencil_multi[i][j]);
        delete [] stencil_multi[i];
      }
      delete [] stencil_multi;
      memory->destroy(maxstencil_multi);
      memory->destroy(flag_half_multi);
      memory->destroy(flag_skip_multi);
      memory->destroy(bin_collection_multi);
      memory->destroy(stencil_sx_multi);
      memory->destroy(stencil_sy_multi);
      memory->destroy(stencil_sz_multi);
      memory->destroy(stencil_binsizex_multi);
      memory->destroy(stencil_binsizey_multi);
      memory->destroy(stencil_binsizez_multi);
      memory->destroy(stencil_mbinx_multi);
      memory->destroy(stencil_mbiny_multi);
      memory->destroy(stencil_mbinz_multi);
    }

    // Allocate arrays
    if (!maxstencil_multi) {
      memory->create(flag_half_multi, n, n,
                     "neighstencil:flag_half_multi");
      memory->create(flag_skip_multi, n, n,
                     "neighstencil:flag_skip_multi");
      memory->create(bin_collection_multi, n, n,
                     "neighstencil:bin_collection_multi");

      memory->create(stencil_sx_multi, n, n,
                     "neighstencil:stencil_sx_multi");
      memory->create(stencil_sy_multi, n, n,
                     "neighstencil:stencil_sy_multi");
      memory->create(stencil_sz_multi, n, n,
                     "neighstencil:stencil_sz_multi");

      memory->create(stencil_binsizex_multi, n, n,
                     "neighstencil:stencil_binsizex_multi");
      memory->create(stencil_binsizey_multi, n, n,
                     "neighstencil:stencil_binsizey_multi");
      memory->create(stencil_binsizez_multi, n, n,
                     "neighstencil:stencil_binsizez_multi");

      memory->create(stencil_mbinx_multi, n, n,
                     "neighstencil:stencil_mbinx_multi");
      memory->create(stencil_mbiny_multi, n, n,
                     "neighstencil:stencil_mbiny_multi");
      memory->create(stencil_mbinz_multi, n, n,
                     "neighstencil:stencil_mbinz_multi");

      memory->create(maxstencil_multi, n, n, "neighstencil::maxstencil_multi");
      memory->create(nstencil_multi, n, n, "neighstencil::nstencil_multi");
      stencil_multi = new int**[n]();
      for (i = 0; i < n; ++i) {
        stencil_multi[i] = new int*[n]();
        for (j = 0; j < n; ++j) {
              maxstencil_multi[i][j] = 0;
          nstencil_multi[i][j] = 0;
          stencil_multi[i][j] = nullptr;
        }
      }
      maxcollections = n;
    }

    // Skip all stencils by default, initialize smax
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        flag_skip_multi[i][j] = true;
      }
    }

    // Determine which stencils need to be built
    set_stencil_properties();

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        // Skip creation of unused stencils
        if (flag_skip_multi[i][j]) continue;

        // Copy bin info for this pair of atom collections
        bin_collection = bin_collection_multi[i][j];

        stencil_binsizex_multi[i][j] = binsizex_multi[bin_collection];
        stencil_binsizey_multi[i][j] = binsizey_multi[bin_collection];
        stencil_binsizez_multi[i][j] = binsizez_multi[bin_collection];

        stencil_mbinx_multi[i][j] = mbinx_multi[bin_collection];
        stencil_mbiny_multi[i][j] = mbiny_multi[bin_collection];
        stencil_mbinz_multi[i][j] = mbinz_multi[bin_collection];

        stencil_range = sqrt(cutcollectionsq[i][j]);

        sx = static_cast<int> (stencil_range*bininvx_multi[bin_collection]);
        if (sx*binsizex_multi[bin_collection] < stencil_range) sx++;
        sy = static_cast<int> (stencil_range*bininvy_multi[bin_collection]);
        if (sy*binsizey_multi[bin_collection] < stencil_range) sy++;
        sz = static_cast<int> (stencil_range*bininvz_multi[bin_collection]);
        if (sz*binsizez_multi[bin_collection] < stencil_range) sz++;
        if (dimension == 2) sz = 0;

        stencil_sx_multi[i][j] = sx;
        stencil_sy_multi[i][j] = sy;
        stencil_sz_multi[i][j] = sz;

        smax = ((2*sx+1) * (2*sy+1) * (2*sz+1));

        if (smax > maxstencil_multi[i][j]) {
          maxstencil_multi[i][j] = smax;
          if(stencil_multi[i][j])
            memory->destroy(stencil_multi[i][j]);
          memory->create(stencil_multi[i][j], smax,
                                 "neighstencil::stencil_multi");
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double NStencil::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;

  return (delx*delx + dely*dely + delz*delz);
}


/* ----------------------------------------------------------------------
   compute closest distance for a given atom collection
------------------------------------------------------------------------- */

double NStencil::bin_distance_multi(int i, int j, int k, int ic)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex_multi[ic];
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex_multi[ic];

  if (j > 0) dely = (j-1)*binsizey_multi[ic];
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey_multi[ic];

  if (k > 0) delz = (k-1)*binsizez_multi[ic];
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez_multi[ic];

  return (delx*delx + dely*dely + delz*delz);
}

/* ---------------------------------------------------------------------- */

double NStencil::memory_usage()
{
  double bytes = 0;
  if (neighstyle == Neighbor::BIN) {
    bytes += memory->usage(stencil,maxstencil);
    bytes += memory->usage(stencilxyz,maxstencil,3);
  } else if (neighstyle == Neighbor::MULTI_OLD) {
    bytes += (double)atom->ntypes*maxstencil_multi_old * sizeof(int);
    bytes += (double)atom->ntypes*maxstencil_multi_old * sizeof(double);
  } else if (neighstyle == Neighbor::MULTI) {
    int n = ncollections;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        bytes += (double)maxstencil_multi[i][j] * sizeof(int);
      }
    }
  }
  return bytes;
}
