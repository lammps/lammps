// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_orient_fcc_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include <cmath>
#include <cstdlib>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int BIG = 1000000000;

/* ---------------------------------------------------------------------- */

FixOrientFCCOMP::FixOrientFCCOMP(LAMMPS *lmp, int narg, char **arg) :
  FixOrientFCC(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded version of FixOrientFCC::post_force().

   Phase 1 (build per-atom Nbr structure + order parameter) is
   embarrassingly parallel over owned atoms: each atom writes only its own
   nbr[i]/order[i] slots.  The shared sort scratch of the serial version is
   replaced by a per-thread buffer, and added_energy / neighbor-count
   statistics are reduced across threads.

   Phase 2 (apply grain-boundary force) writes ONLY f[i] for the owned
   atom i -- it reads the neighbor data nbr[m] (acquired for ghosts via the
   forward_comm between phases) but never writes f[m]/f[j].  Hence it is
   race-free over a per-atom partition and needs no per-thread force buffer.
------------------------------------------------------------------------- */

void FixOrientFCCOMP::post_force(int /*vflag*/)
{
  // set local ptrs

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;
  const int nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // ensure nbr and order data structures are adequate size

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(nbr);
    memory->destroy(order);
    nbr = (Nbr *) memory->smalloc(nmax*sizeof(Nbr),"orient/fcc:nbr");
    memory->create(order,nmax,2,"orient/fcc:order");
    array_atom = order;
  }

  added_energy = 0.0;
  int count = 0;
  int mincount = BIG;
  int maxcount = 0;

  // loop over owned atoms and build Nbr data structure of neighbors
  // use full neighbor list

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    Sort *sort = nullptr;
    int sortmax = 0;
    double e_thr = 0.0;
    int count_thr = 0, minc_thr = BIG, maxc_thr = 0;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      if (jnum < minc_thr) minc_thr = jnum;
      if (jnum > maxc_thr) maxc_thr = jnum;
      if (jnum > sortmax) {
        delete [] sort;
        sortmax = jnum;
        sort = new Sort[sortmax];
      }

      // loop over all neighbors of atom i
      // for those within cutsq, build sort data structure

      int nsort = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        count_thr++;

        const double dx = x[i][0] - x[j][0];
        const double dy = x[i][1] - x[j][1];
        const double dz = x[i][2] - x[j][2];
        const double rsq = dx*dx + dy*dy + dz*dz;

        if (rsq < cutsq) {
          sort[nsort].id = j;
          sort[nsort].rsq = rsq;
          sort[nsort].delta[0] = dx;
          sort[nsort].delta[1] = dy;
          sort[nsort].delta[2] = dz;
          if (use_xismooth) {
            const double xismooth = (xicutoffsq - 2.0*rsq/(a*a)) / (xicutoffsq - 1.0);
            sort[nsort].xismooth = 1.0 - fabs(1.0-xismooth);
          }
          nsort++;
        }
      }

      // sort neighbors by rsq distance
      // no need to sort if nsort <= 12

      if (nsort > 12) qsort(sort,nsort,sizeof(Sort),compare);

      // copy up to 12 nearest neighbors into nbr data structure
      // operate on delta vector via find_best_ref() to compute dxi

      const int n = MIN(12,nsort);
      nbr[i].n = n;
      if (n == 0) continue;

      double xi_total = 0.0;
      double xi_sq, dxi[3];
      for (int j = 0; j < n; j++) {
        find_best_ref(sort[j].delta,0,xi_sq,dxi);
        xi_total += sqrt(xi_sq);
        nbr[i].id[j] = sort[j].id;
        nbr[i].dxi[j][0] = dxi[0]/n;
        nbr[i].dxi[j][1] = dxi[1]/n;
        nbr[i].dxi[j][2] = dxi[2]/n;
        if (use_xismooth) nbr[i].xismooth[j] = sort[j].xismooth;
      }
      xi_total /= n;
      order[i][0] = xi_total;

      // compute potential derivative to xi

      double edelta;
      if (xi_total < xi0) {
        nbr[i].duxi = 0.0;
        edelta = 0.0;
        order[i][1] = 0.0;
      } else if (xi_total > xi1) {
        nbr[i].duxi = 0.0;
        edelta = Vxi;
        order[i][1] = 1.0;
      } else {
        const double omega = MY_PI2*(xi_total-xi0) / (xi1-xi0);
        nbr[i].duxi = MY_PI*Vxi*sin(2.0*omega) / (2.0*(xi1-xi0));
        edelta = Vxi*(1 - cos(2.0*omega)) / 2.0;
        order[i][1] = omega / MY_PI2;
      }
      e_thr += edelta;
    }

    delete [] sort;

#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      added_energy += e_thr;
      count += count_thr;
      if (minc_thr < mincount) mincount = minc_thr;
      if (maxc_thr > maxcount) maxcount = maxc_thr;
    }
  }

  // communicate to acquire nbr data for ghost atoms

  comm->forward_comm(this);

  // compute grain boundary force on each owned atom
  // skip atoms not in group
  // only f[i] is written, so the per-atom partition is race-free

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic) default(shared)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    const int n = nbr[i].n;
    const double duxi = nbr[i].duxi;

    for (int j = 0; j < n; j++) {
      double *dxiptr = &nbr[i].dxi[j][0];
      if (use_xismooth) {
        const double xismooth = nbr[i].xismooth[j];
        f[i][0] += duxi * dxiptr[0] * xismooth;
        f[i][1] += duxi * dxiptr[1] * xismooth;
        f[i][2] += duxi * dxiptr[2] * xismooth;
      } else {
        f[i][0] += duxi * dxiptr[0];
        f[i][1] += duxi * dxiptr[1];
        f[i][2] += duxi * dxiptr[2];
      }

      // m = local index of neighbor
      // id_self = ID for atom I in atom M's neighbor list
      // if M is local atom, id_self will be local ID of atom I
      // if M is ghost atom, id_self will be global ID of atom I

      const int m = nbr[i].id[j];
      tagint id_self;
      if (m < nlocal) id_self = i;
      else id_self = tag[i];
      bool found_myself = false;
      const int nn = nbr[m].n;

      for (int k = 0; k < nn; k++) {
        if (id_self == nbr[m].id[k]) {
          if (found_myself) error->one(FLERR,"Fix orient/fcc found self twice");
          found_myself = true;
          const double duxi_other = nbr[m].duxi;
          dxiptr = &nbr[m].dxi[k][0];
          if (use_xismooth) {
            const double xismooth = nbr[m].xismooth[k];
            f[i][0] -= duxi_other * dxiptr[0] * xismooth;
            f[i][1] -= duxi_other * dxiptr[1] * xismooth;
            f[i][2] -= duxi_other * dxiptr[2] * xismooth;
          } else {
            f[i][0] -= duxi_other * dxiptr[0];
            f[i][1] -= duxi_other * dxiptr[1];
            f[i][2] -= duxi_other * dxiptr[2];
          }
        }
      }
    }
  }

  // print statistics every nstats timesteps

  if (nstats && update->ntimestep % nstats == 0) {
    int total;
    MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,world);
    double ave = static_cast<double>(total)/atom->natoms;

    int min,max;
    MPI_Allreduce(&mincount,&min,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&maxcount,&max,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      std::string mesg = fmt::format("orient step {}: {} atoms have {} "
                                     "neighbors\n", update->ntimestep,
                                     atom->natoms,total);
      mesg += fmt::format("  neighs: min = {}, max ={}, ave = {}\n",
                          min,max,ave);
      utils::logmesg(lmp,mesg);
    }
  }
}
