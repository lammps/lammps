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
   ILVES constraint solver: undirected graph utilities.
   Ported from GROMACS 2021 ILVES (LGPL-2.1), src/gromacs/mdlib/graph.cpp.
   The reference's k-way partitioning paths (for multi-threaded solves) are
   omitted: this port is serial, so only the CSR adjacency structure remains.
   See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#include "ilves_graph.h"

#include <vector>

namespace LAMMPS_NS::ILVES {

Graph::Graph() : nnodes(0) {}

Graph::Graph(const int nodes, const int edges) : nnodes(nodes), xadj(nodes + 1), adj(edges)
{
  xadj[0] = 0;
}

int Graph::num_nodes() const
{
  return nnodes;
}

int Graph::num_edges() const
{
  return adj.size();
}

double Graph::memory_usage() const
{
  return (double) (xadj.size() + adj.size()) * sizeof(int);
}

}    // namespace LAMMPS_NS::ILVES
