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

#include "nstencil_CAC.h"
#include "neighbor.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;

enum{NSQ,BIN,MULTI};     // also in Neighbor

/* ---------------------------------------------------------------------- */

NStencilCAC::NStencilCAC(LAMMPS *lmp) : NStencil(lmp) {
  post_create_flag=1;
}

/*empty create setup */
void NStencilCAC::create_setup() {}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilCAC::create()
{
  int i,j,k;
  comm->stencil_pointer=this;

}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff just before npair is called
------------------------------------------------------------------------- */

void NStencilCAC::post_create()
{
  int i,j,k;

  nstencil = 0;

  for (k = -sz; k <= sz; k++)
    for (j = -sy; j <= sy; j++)
      for (i = -sx; i <= sx; i++)
        if (bin_distance(i,j,k) < cutneighmaxsq)
          stencil[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
}

/* ----------------------------------------------------------------------
   insure NBin data is current
   insure stencils are allocated large enough
------------------------------------------------------------------------- */

void NStencilCAC::post_create_setup()
{
  if (nb) copy_bin_info();
  last_stencil = update->ntimestep;
  double cutneighmax= neighbor->cutneighmax;

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

  if (neighstyle == BIN) {
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
        stencil_multi[i] = NULL;
        distsq_multi[i] = NULL;
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
}