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
   Contributing author: Daniel Schwen
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_voronoi_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

#include <vector>
#include "voro++.hh"

using namespace LAMMPS_NS;
using namespace voro;

/* ---------------------------------------------------------------------- */

ComputeVoronoi::ComputeVoronoi(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute voronoi/atom command");

  peratom_flag = 1;
  size_peratom_cols = 2;

  nmax = 0;
  voro = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeVoronoi::~ComputeVoronoi()
{
  memory->destroy(voro);
}

/* ---------------------------------------------------------------------- */

void ComputeVoronoi::init()
{
  if (domain->triclinic != 0)
    error->all(FLERR,"Compute voronoi/atom not allowed for triclinic boxes");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"voronoi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute voronoi/atom command");
}

/* ----------------------------------------------------------------------
   gather compute vector data from other nodes
------------------------------------------------------------------------- */

void ComputeVoronoi::compute_peratom()
{
  int i;

  invoked_peratom = update->ntimestep;

  // grow per atom array if necessary

  int nlocal = atom->nlocal;
  if (nlocal > nmax) {
    memory->destroy(voro);
    nmax = atom->nmax;
    memory->create(voro,nmax,size_peratom_cols,"voronoi/atom:voro");
    array_atom = voro;
  }

  // n = # of voro++ spatial hash cells
  // TODO: make square

  int nall = nlocal + atom->nghost;
  int n = int(floor( pow( double(nall)/8.0, 0.333333 ) ));
  n = (n==0) ? 1 : n;

  // initialize voro++ container
  // preallocates 8 atoms per cell
  // voro++ allocates more memory if needed

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  double *cut = comm->cutghost;

  container con(sublo[0]-cut[0],subhi[0]+cut[0],
                sublo[1]-cut[1],subhi[1]+cut[1],
                sublo[2]-cut[2],subhi[2]+cut[2],
                n,n,n,false,false,false,8); 

  // pass coordinates for local and ghost atoms to voro++

  double **x = atom->x;
  for (i = 0; i < nall; i++)
    con.put(i,x[i][0],x[i][1],x[i][2]);

  // invoke voro++ and fetch results for owned atoms in group

  int *mask = atom->mask;
  std::vector<int> neigh;

  voronoicell_neighbor c;
  c_loop_all cl(con);
  if (cl.start()) do if (con.compute_cell(c,cl)) {
    i = cl.pid();
    if (i < nlocal && (mask[i] & groupbit)) {
      voro[i][0] = c.volume();
      c.neighbors(neigh);
      voro[i][1] = neigh.size();
    } else if (i < nlocal) voro[i][0] = voro[i][1] = 0.0;
  } while (cl.inc());
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeVoronoi::memory_usage()
{
  double bytes = size_peratom_cols * nmax * sizeof(double);
  return bytes;
}
