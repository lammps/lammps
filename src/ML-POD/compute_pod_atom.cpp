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

#include "compute_pod_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

#include "eapod.h"

using namespace LAMMPS_NS;

enum { SCALAR, VECTOR, ARRAY };

ComputePODAtom::ComputePODAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), podptr(nullptr), pod(nullptr), tmpmem(nullptr),
    rij(nullptr), elements(nullptr), map(nullptr)
{
  int nargmin = 6;

  if (narg < nargmin) error->all(FLERR, "Illegal compute {} command", style);
  if (comm->nprocs > 1) error->all(FLERR, "compute command does not support multi processors");

  std::string pod_file = std::string(arg[3]);      // pod input file
  std::string coeff_file = std::string(arg[4]);    // coefficient input file
  podptr = new EAPOD(lmp, pod_file, coeff_file);

  int ntypes = atom->ntypes;
  memory->create(map, ntypes + 1, "compute_pod_global:map");
  map_element2type(narg - 5, arg + 5, podptr->nelements);

  cutmax = podptr->rcut;

  nmax = 0;
  nijmax = 0;
  pod = nullptr;
  elements = nullptr;

  size_peratom_cols = podptr->Mdesc * podptr->nClusters;
  peratom_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputePODAtom::~ComputePODAtom()
{
  memory->destroy(map);
  memory->destroy(pod);
  delete podptr;
}

/* ---------------------------------------------------------------------- */

void ComputePODAtom::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute pod requires a pair style be defined");

  if (cutmax > force->pair->cutforce)
    error->all(FLERR,"Compute pod cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("pod").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute pod");
}


/* ---------------------------------------------------------------------- */

void ComputePODAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePODAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow pod array if necessary

  if (atom->natoms > nmax) {
    memory->destroy(pod);
    nmax = atom->natoms;
    int numdesc = podptr->Mdesc * podptr->nClusters;
    memory->create(pod, nmax, numdesc,"sna/atom:sna");
    array_atom = pod;
  }

  for (int i = 0; i < atom->natoms; i++)
    for (int icoeff = 0; icoeff < size_peratom_cols; icoeff++) {
      pod[i][icoeff] = 0.0;
    }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  double **x = atom->x;
  int **firstneigh = list->firstneigh;
  int *numneigh = list->numneigh;
  int *type = atom->type;
  int *ilist = list->ilist;
  int inum = list->inum;
  int nClusters = podptr->nClusters;
  int Mdesc = podptr->Mdesc;
  double rcutsq = podptr->rcut*podptr->rcut;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int jnum = numneigh[i];

    // allocate temporary memory
    if (nijmax < jnum) {
      nijmax = MAX(nijmax, jnum);
      podptr->free_temp_memory();
      podptr->allocate_temp_memory(nijmax);
    }

    rij = &podptr->tmpmem[0];
    tmpmem = &podptr->tmpmem[3*nijmax];
    ai = &podptr->tmpint[0];
    aj = &podptr->tmpint[nijmax];
    ti = &podptr->tmpint[2*nijmax];
    tj = &podptr->tmpint[3*nijmax];

    // get neighbor list for atom i
    lammpsNeighborList(x, firstneigh, atom->tag, type, numneigh, rcutsq, i);

    if (nij > 0) {
      // peratom base descriptors
      double *bd = &podptr->bd[0];
      double *bdd = &podptr->bdd[0];
      podptr->peratombase_descriptors(bd, bdd, rij, tmpmem, tj, nij);

      if (nClusters>1) {
        // peratom env descriptors
        double *pd = &podptr->pd[0];
        double *pdd = &podptr->pdd[0];
        podptr->peratomenvironment_descriptors(pd, pdd, bd, bdd, tmpmem, ti[0] - 1,  nij);
        for (int k = 0; k < nClusters; k++)
          for (int m = 0; m < Mdesc; m++) {
            int mk = m + Mdesc*k;
            pod[i][mk] = pd[k]*bd[m];
          //   for (int n=0; n<nij; n++) {
          //     int ain = 3*ai[n];
          //     int ajn = 3*aj[n];
          //     int nm = 3*n + 3*nij*m;
          //     int nk = 3*n + 3*nij*k;
          //     pod[1 + ain][imk] += bdd[0 + nm]*pd[k] + bd[m]*pdd[0+nk];
          //     pod[2 + ain][imk] += bdd[1 + nm]*pd[k] + bd[m]*pdd[1+nk];
          //     pod[3 + ain][imk] += bdd[2 + nm]*pd[k] + bd[m]*pdd[2+nk];
          //     pod[1 + ajn][imk] -= bdd[0 + nm]*pd[k] + bd[m]*pdd[0+nk];
          //     pod[2 + ajn][imk] -= bdd[1 + nm]*pd[k] + bd[m]*pdd[1+nk];
          //     pod[3 + ajn][imk] -= bdd[2 + nm]*pd[k] + bd[m]*pdd[2+nk];
          //   }
          }
      }
      else {
        for (int m = 0; m < Mdesc; m++) {
         pod[i][m] = bd[m];
        //  for (int n=0; n<nij; n++) {
        //     int ain = 3*ai[n];
        //     int ajn = 3*aj[n];
        //     int nm = 3*n + 3*nij*m;
        //     pod[1 + ain][im] += bdd[0 + nm];
        //     pod[2 + ain][im] += bdd[1 + nm];
        //     pod[3 + ain][im] += bdd[2 + nm];
        //     pod[1 + ajn][im] -= bdd[0 + nm];
        //     pod[2 + ajn][im] -= bdd[1 + nm];
        //     pod[3 + ajn][im] -= bdd[2 + nm];
        //   }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputePODAtom::memory_usage()
{
  double bytes = 0.0;

  return bytes;
}

void ComputePODAtom::lammpsNeighborList(double **x, int **firstneigh, tagint *atomid, int *atomtypes,
                               int *numneigh, double rcutsq, int gi)
{
  nij = 0;
  int itype = map[atomtypes[gi]] + 1;
  int m = numneigh[gi];
  ti[nij] = itype;
  for (int l = 0; l < m; l++) {           // loop over each atom around atom i
    int gj = firstneigh[gi][l];           // atom j
    double delx = x[gj][0] - x[gi][0];    // xj - xi
    double dely = x[gj][1] - x[gi][1];    // xj - xi
    double delz = x[gj][2] - x[gi][2];    // xj - xi
    double rsq = delx * delx + dely * dely + delz * delz;
    if (rsq < rcutsq && rsq > 1e-20) {
      rij[nij * 3 + 0] = delx;
      rij[nij * 3 + 1] = dely;
      rij[nij * 3 + 2] = delz;
      ai[nij] = atomid[gi]-1;
      aj[nij] = atomid[gj]-1;
      ti[nij] = itype;
      tj[nij] = map[atomtypes[gj]] + 1;
      nij++;
    }
  }
}

void ComputePODAtom::map_element2type(int narg, char **arg, int nelements)
{
  int i,j;
  const int ntypes = atom->ntypes;

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"
  // nelements = # of unique elements
  // elements = list of element names

  if (narg != ntypes)
    error->all(FLERR, "Number of element to type mappings does not match number of atom types");

  if (elements) {
    for (i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;
  }
  elements = new char*[ntypes];
  for (i = 0; i < ntypes; i++) elements[i] = nullptr;

  nelements = 0;
  map[0] = -1;
  for (i = 1; i <= narg; i++) {
    std::string entry = arg[i-1];
    if (entry == "NULL") {
      map[i] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (entry == elements[j]) break;
    map[i] = j;
    if (j == nelements) {
      elements[j] = utils::strdup(entry);
      nelements++;
    }
  }
}
