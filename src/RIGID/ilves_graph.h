/* -*- c++ -*- ----------------------------------------------------------
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
   ILVES constraint solver: undirected graph utilities (CSR storage and
   vertex renumbering for fill reduction).

   Ported from the GROMACS 2021 reference implementation of ILVES
   (https://github.com/LorienLV/_PAPER_ILVES, LGPL-2.1),
   src/gromacs/mdlib/graph.{cpp,h}.  The reference GROMACS formatting and
   logic are preserved intentionally to keep the port diffable against
   upstream; only the namespace and the (disabled) METIS code paths differ.

   Reference:
     L. Lopez-Villellas, C. C. K. Mikkelsen, J. J. Galano-Frutos,
     S. Marco-Sola, J. Alastruey-Bende, P. Ibanez, P. Echenique,
     M. Moreto, M. C. De Rosa, and P. Garcia-Risueno, "ILVES: Accurate and
     Efficient Bond Length and Angle Constraints in Molecular Dynamics",
     J. Chem. Theory Comput. 21, 8711-8719 (2025),
     doi:10.1021/acs.jctc.5c01376

   LAMMPS port: Axel Kohlmeyer (Temple U) with Claude Code (Opus 4.8)
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_GRAPH_H
#define LMP_ILVES_GRAPH_H

#include <cstdint>
#include <vector>

/*
 * This is a module for manipulating general undirected graphs (vertices and
 * edges). Note: In this representation, a node is connected by other node with
 * exactly one edge.
 *
 */

namespace LAMMPS_NS::ILVES {

// Data structure for graphs
class Graph {
 public:
  /*
   * Here m is the number of vertices, and xadj[i] is the index in adj of the
   * first element in the ith adjency list. The adjacency lists are stored
   * contiguously in the array adj. This is a very compact way of storing
   * graphs.
   *
   * Example:
   *
   * 0 - 1
   * |
   * 2
   *
   * m = 3
   * xadj = [0, 3, 5, 7]
   * adj = [0, 1, 2, 0, 1, 0, 2]
   *
   */
  int nnodes;
  std::vector<int> xadj;
  std::vector<int> adj;

  /**
   * Constructs an empty graph.
   */
  Graph();

  /**
   * Constructs a graph and reserves space for the nodes and edges.
   *
   * @param nodes The number of nodes of the graph.
   * @param edges The number of edges of the graph.
   */
  Graph(int nodes, int edges);

  /**
   * Get the number of nodes in the graph.
   */
  [[nodiscard]] int num_nodes() const;

  /**
   * Get the number of edges in the graph.
   */
  [[nodiscard]] int num_edges() const;

  /**
   * Estimate the memory used by the graph (CSR index + adjacency arrays).
   *
   * @return The size of the graph storage in bytes.
   */
  [[nodiscard]] double memory_usage() const;
};

}    // namespace LAMMPS_NS::ILVES

#endif
