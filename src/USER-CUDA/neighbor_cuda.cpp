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
#include "neighbor_cuda.h"
#include "cuda.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;




enum {NSQ, BIN, MULTI};  // also in neigh_list.cpp

/* ---------------------------------------------------------------------- */

NeighborCuda::NeighborCuda(LAMMPS* lmp) : Neighbor(lmp)
{
  cuda = lmp->cuda;

  if(cuda == NULL)
    error->all(FLERR, "You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");
}

/* ---------------------------------------------------------------------- */

void NeighborCuda::init()
{
  cuda->set_neighinit(dist_check, 0.25 * skin * skin);
  cudable = 1;

  Neighbor::init();
}

/* ----------------------------------------------------------------------
   overwrite either full_nsq or full_bin with CUDA-equivalent methods
   any other neighbor build method is unchanged
------------------------------------------------------------------------- */

void NeighborCuda::choose_build(int index, NeighRequest* rq)
{
  Neighbor::choose_build(index, rq);

  if(rq->full && style == NSQ && rq->cudable)
    pair_build[index] = (Neighbor::PairPtr) &NeighborCuda::full_nsq_cuda;
  else if(rq->full && style == BIN && rq->cudable)
    pair_build[index] = (Neighbor::PairPtr) &NeighborCuda::full_bin_cuda;
}

/* ---------------------------------------------------------------------- */

int NeighborCuda::check_distance()
{
  double delx, dely, delz, rsq;
  double delta, deltasq, delta1, delta2;

  if(boxcheck) {
    if(triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx * delx + dely * dely + delz * delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx * delx + dely * dely + delz * delz);
      delta = 0.5 * (skin - (delta1 + delta2));
      deltasq = delta * delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;

      for(int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx * delx + dely * dely + delz * delz);

        if(delta > delta1) delta1 = delta;
        else if(delta > delta2) delta2 = delta;
      }

      delta = 0.5 * (skin - (delta1 + delta2));
      deltasq = delta * delta;
    }
  } else deltasq = triggersq;

  double** x = atom->x;
  int nlocal = atom->nlocal;

  if(includegroup) nlocal = atom->nfirst;

  int flag = 0;

  if(not cuda->neighbor_decide_by_integrator) {
    cuda->cu_x_download();

    for(int i = 0; i < nlocal; i++) {
      delx = x[i][0] - xhold[i][0];
      dely = x[i][1] - xhold[i][1];
      delz = x[i][2] - xhold[i][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if(rsq > deltasq) flag = 1;
    }
  } else flag = cuda->shared_data.atom.reneigh_flag;

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);

  if(flagall && ago == MAX(every, delay)) ndanger++;

  return flagall;
}

/* ---------------------------------------------------------------------- */

void NeighborCuda::build(int topoflag)
{
  int i;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;
  // store current atom positions and box size if needed

  if(dist_check) {
    if(cuda->decide_by_integrator())
      cuda->update_xhold(maxhold, &xhold[0][0]);
    else {
      if(cuda->finished_setup) cuda->cu_x_download();

      double** x = atom->x;
      int nlocal = atom->nlocal;

      if(includegroup) nlocal = atom->nfirst;

      if(nlocal > maxhold) {
        maxhold = atom->nmax;
        memory->destroy(xhold);
        memory->create(xhold, maxhold, 3, "neigh:xhold");
      }

      for(i = 0; i < nlocal; i++) {
        xhold[i][0] = x[i][0];
        xhold[i][1] = x[i][1];
        xhold[i][2] = x[i][2];
      }

      if(boxcheck) {
        if(triclinic == 0) {
          boxlo_hold[0] = bboxlo[0];
          boxlo_hold[1] = bboxlo[1];
          boxlo_hold[2] = bboxlo[2];
          boxhi_hold[0] = bboxhi[0];
          boxhi_hold[1] = bboxhi[1];
          boxhi_hold[2] = bboxhi[2];
        } else {
          domain->box_corners();
          corners = domain->corners;

          for(i = 0; i < 8; i++) {
            corners_hold[i][0] = corners[i][0];
            corners_hold[i][1] = corners[i][1];
            corners_hold[i][2] = corners[i][2];
          }
        }
      }
    }
  }

  if(not cudable && cuda->finished_setup && atom->avec->cudable)
    cuda->downloadAll();

  if(cudable && (not cuda->finished_setup)) {
    cuda->checkResize();
    cuda->uploadAll();
  }

  // if any lists store neighbors of ghosts:
  // invoke grow() if nlocal+nghost exceeds previous list size
  // else only invoke grow() if nlocal exceeds previous list size
  // only done for lists with growflag set and which are perpetual

  if(anyghostlist && atom->nlocal + atom->nghost > maxatom) {
    maxatom = atom->nmax;

    for(i = 0; i < nglist; i++) lists[glist[i]]->grow(maxatom);
  } else if(atom->nlocal > maxatom) {
    maxatom = atom->nmax;

    for(i = 0; i < nglist; i++) lists[glist[i]]->grow(maxatom);
  }

  // extend atom bin list if necessary

  if(style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->destroy(bins);
    memory->create(bins, maxbin, "bins");
  }

  // check that neighbor list with special bond flags will not overflow

  if(atom->nlocal + atom->nghost > NEIGHMASK)
    error->one(FLERR, "Too many local+ghost atoms for neighbor list");

  // invoke building of pair and molecular neighbor lists
  // only for pairwise lists with buildflag set

  for(i = 0; i < nblist; i++)
    (this->*pair_build[blist[i]])(lists[blist[i]]);

  if(atom->molecular && topoflag) {
    if(force->bond)(this->*bond_build)();
    if(force->angle)(this->*angle_build)();
    if(force->dihedral)(this->*dihedral_build)();
    if(force->improper)(this->*improper_build)();
  }
}
