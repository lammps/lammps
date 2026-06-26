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
   omitted: this port is serial, so only the CSR adjacency structure and the
   vertex renumbering used for fill reduction remain.  See ilves_graph.h for
   full attribution.
------------------------------------------------------------------------- */

#include "ilves_graph.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace LAMMPS_NS {
namespace ILVES {

Graph::Graph() : nnodes(0) {}

Graph::Graph(const int nodes, const int edges) : nnodes(nodes), xadj(nodes + 1), adj(edges)
{
  xadj[0] = 0;
}

int Graph::num_nodes() const { return nnodes; }

int Graph::num_edges() const { return adj.size(); }

void Graph::renumber_vertices(const std::vector<int> &perm, std::vector<int> &iperm)
{
  // Compute the inverse permutation if it is not provided.
  if (iperm.empty()) {
    iperm.resize(nnodes);
    for (int node = 0; node < nnodes; ++node) {
      iperm[perm[node]] = node;
    }
  }

  Graph perm_graph(nnodes, adj.size());

  for (int i = 0, edge = 0; i < nnodes; ++i) {
    const int node = perm[i];    // The ith node is the old perm[i] node.

    const int old_edge = edge;
    for (int k = xadj[node]; k < xadj[node + 1]; ++k, ++edge) {
      // The ith node is connected to iperm[adj[k]].
      perm_graph.adj[edge] = iperm[adj[k]];
    }
    // The adjacency list of node may not be in order.
    std::sort(perm_graph.adj.begin() + old_edge, perm_graph.adj.begin() + edge);

    perm_graph.xadj[i + 1] = edge;
  }

  *this = std::move(perm_graph);
}

double Graph::memory_usage() const
{
  return (double) (xadj.size() + adj.size()) * sizeof(int);
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
