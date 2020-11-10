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

#include "nstencil.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "nbin.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"

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
   for half list with newton off:
     stencil is all surrounding bins including self
     regardless of triclinic
   for half list with newton on:
     stencil is bins to the "upper right" of central bin
     stencil does not include self
     no versions that allow ghost on (no callers need it?)
   for half list with newton on and triclinic:
     stencil is all bins in z-plane of self and above, but not below
     in 2d is all bins in y-plane of self and above, but not below
     stencil includes self
     no versions that allow ghost on (no callers need it?)
   for full list:
     stencil is all surrounding bins including self
     regardless of newton on/off or triclinic
   for multi:
     create one stencil for each atom type
     stencil follows same rules for half/full, newton on/off, triclinic
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
     no versions that allow ghost on (any need for it?)
   for multi/tiered:
     create one stencil for each itype-jtype pairing
     stencils do not generally follow the same rules for half/full or newton on/off
     whole stencils including all surrounding bins are always used except
     for same-type stencils with newton on which uses a split stencil
     for orthogonal boxes, a split stencil includes bins to the "upper right" of central bin 
     for triclinic, a split stencil includes bins in the z (3D) or y (2D) plane of self and above
     cutoff is not cutneighmaxsq, but max cutoff for that atom type
     no versions that allow ghost on (any need for it?)
------------------------------------------------------------------------- */

NStencil::NStencil(LAMMPS *lmp) : Pointers(lmp)
{
  last_stencil = -1;

  xyzflag = 0;
  maxstencil = maxstencil_multi = 0;
  stencil = nullptr;
  stencilxyz = nullptr;
  nstencil_multi = nullptr;
  stencil_multi = nullptr;
  distsq_multi = nullptr;
  
  stencil_split = nullptr;
  stencil_skip = nullptr;
  stencil_bin_type = nullptr;
  stencil_cut = nullptr;

  dimension = domain->dimension;
}

/* ---------------------------------------------------------------------- */

NStencil::~NStencil()
{
  memory->destroy(stencil);
  memory->destroy(stencilxyz);

  if (stencil_multi) {

    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      memory->destroy(stencil_multi[i]);
      memory->destroy(distsq_multi[i]);
    }
    delete [] nstencil_multi;
    delete [] stencil_multi;
    delete [] distsq_multi;
  }
  
  if (stencil_multi_tiered) {
      
    int n = atom->ntypes;      
    memory->destroy(nstencil_multi_tiered);
    for (int i = 1; i <= n; i++) {
      for (int j = 0; j <= n; j++) {
        if (! stencil_skip[i][j]) 
            memory->destroy(stencil_multi_tiered[i][j]);
      }
      delete [] stencil_multi_tiered[i];
    }
    delete [] stencil_multi_tiered;
    memory->destroy(maxstencil_multi_tiered);
    memory->destroy(stencil_split);
    memory->destroy(stencil_skip);
    memory->destroy(stencil_bin_type);
    memory->destroy(stencil_cut);
    
    memory->destroy(sx_multi_tiered);
    memory->destroy(sy_multi_tiered);
    memory->destroy(sz_multi_tiered);
    
    memory->destroy(binsizex_multi_tiered);
    memory->destroy(binsizey_multi_tiered);
    memory->destroy(binsizez_multi_tiered);
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
   copy needed info for a given type from NBin class to this stencil class
------------------------------------------------------------------------- */

void NStencil::copy_bin_info_multi_tiered(int type)
{
  mbinx = nb->mbinx_tiered[type];
  mbiny = nb->mbiny_tiered[type];
  mbinz = nb->mbinz_tiered[type];
  binsizex = nb->binsizex_tiered[type];
  binsizey = nb->binsizey_tiered[type];
  binsizez = nb->binsizez_tiered[type];
  bininvx = nb->bininvx_tiered[type];
  bininvy = nb->bininvy_tiered[type];
  bininvz = nb->bininvz_tiered[type];
}

/* ----------------------------------------------------------------------
   insure NBin data is current
   insure stencils are allocated large enough
------------------------------------------------------------------------- */

void NStencil::create_setup()
{
    
  if (neighstyle != Neighbor::MULTI_TIERED){
    if (nb) copy_bin_info();
    last_stencil = update->ntimestep;
    
    // sx,sy,sz = max range of stencil in each dim
    // smax = max possible size of entire 3d stencil
    // stencil will be empty if cutneighmax = 0.0
    
    sx = static_cast<int> (cutneighmax*bininvx);
    if (sx*binsizex < cutneighmax) sx++;
    sy = static_cast<int> (cutneighmax*bininvy);
    if (sy*binsizey < cutneighmax) sy++;
    sz = static_cast<int> (cutneighmax*bininvz);
    if (sz*binsizez < cutneighmax) sz++;
    if (dimension == 2) sz = 0;
    
    int smax = (2*sx+1) * (2*sy+1) * (2*sz+1);
    
    // reallocate stencil structs if necessary
    // for BIN and MULTI styles
    
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
      if (maxstencil_multi == 0) {
        nstencil_multi = new int[n+1];
        stencil_multi = new int*[n+1];
        distsq_multi = new double*[n+1];
        for (i = 1; i <= n; i++) {
          nstencil_multi[i] = 0;
          stencil_multi[i] = nullptr;
          distsq_multi[i] = nullptr;
        }
      }
      if (smax > maxstencil_multi) {
        maxstencil_multi = smax;
        for (i = 1; i <= n; i++) {
          memory->destroy(stencil_multi[i]);
          memory->destroy(distsq_multi[i]);
          memory->create(stencil_multi[i],maxstencil_multi,
                         "neighstencil:stencil_multi");
          memory->create(distsq_multi[i],maxstencil_multi,
                         "neighstencil:distsq_multi");
        }
      }
    }
  } else {
    int i, j, bin_type, smax;
    double stencil_range;
    int n = atom->ntypes;

    // Allocate arrays to store stencil information
    memory->create(stencil_split, n, n, 
                   "neighstencil:stencil_split");"
    memory->create(stencil_skip, n, n, 
                   "neighstencil:stencil_skip");"
    memory->create(stencil_bin_type, n, n, 
                   "neighstencil:stencil_bin_type");"
    memory->create(stencil_cut, n, n, 
                   "neighstencil:stencil_cut");"

    memory->create(sx_multi_tiered, n, n, 
                   "neighstencil:sx_multi_tiered");"               
    memory->create(sy_multi_tiered, n, n, 
                   "neighstencil:sy_multi_tiered");"
    memory->create(sz_multi_tiered, n, n, 
                   "neighstencil:sz_multi_tiered");"

    memory->create(binsizex_multi_tiered, n, n, 
                   "neighstencil:binsizex_multi_tiered");"
    memory->create(binsizey_multi_tiered, n, n, 
                   "neighstencil:binsizey_multi_tiered");"
    memory->create(binsizez_multi_tiered, n, n, 
                   "neighstencil:binsizez_multi_tiered");"
                   
    // Determine which stencils need to be built
    set_stencil_properties(); 
    
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        
        // Skip creation of unused stencils 
        if (stencil_skip[i][j]) continue;
           
        // Copy bin info for this particular pair of types
        bin_type = stencil_bin_type[i][j];
        copy_bin_info_bytype(bin_type);
        
        binsizex_multi_tiered[i][j] = binsizex;
        binsizey_multi_tiered[i][j] = binsizey;
        binsizez_multi_tiered[i][j] = binsizez;
        
        stencil_range = stencil_cut[i][j];
        
        sx = static_cast<int> (stencil_range*bininvx);
        if (sx*binsizex < stencil_range) sx++;
        sy = static_cast<int> (stencil_range*bininvy);
        if (sy*binsizey < stencil_range) sy++;
        sz = static_cast<int> (stencil_range*bininvz);
        if (sz*binsizez < stencil_range) sz++;
        
        sx_multi_tiered[i][j] = sx;
        sy_multi_tiered[i][j] = sy;
        sz_multi_tiered[i][j] = sz;

        smax = ((2*sx+1) * (2*sy+1) * (2*sz+1));        
        
        if (smax > maxstencil_multi_tiered[i][j]) {
          maxstencil_multi_tiered[i][j] = smax;
          memory->destroy(stencil_multi_tiered[i][j]);
          memory->create(stencil_multi_tiered[i][j], smax,
	  	                 "neighstencil::stencil_multi_tiered");
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

/* ---------------------------------------------------------------------- */

double NStencil::memory_usage()
{
  double bytes = 0;
  if (neighstyle == Neighbor::BIN) {
    bytes += memory->usage(stencil,maxstencil);
    bytes += memory->usage(stencilxyz,maxstencil,3);
  } else if (neighstyle == Neighbor::MULTI) {
    bytes += atom->ntypes*maxstencil_multi * sizeof(int);
    bytes += atom->ntypes*maxstencil_multi * sizeof(double);
  } else if (neighstyle == Neighbor::MULTI_TIERED) {
    bytes += atom->ntypes*maxstencil_multi * sizeof(int);
    bytes += atom->ntypes*maxstencil_multi * sizeof(int);
    bytes += atom->ntypes*maxstencil_multi * sizeof(double);
  }
  return bytes;
}
