/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#pragma once
#include "neigh_list.h"
#include "atom.h"
#include <vector>

using namespace LAMMPS_NS;

//since LAMMPS is compiled with C++ 2003, define a substitution for std::array
template<typename T, int N>
class array2003{
public:
  T& operator[] (int idx){ return data[idx];};
  const T& operator[] (int idx) const{ return data[idx];};
private:
  T data[N];
};


class MESONTList {
public:
  MESONTList(const Atom* atom, const NeighList* nblist, double rc2);
  ~MESONTList() {};
  //list of segments
  const std::vector<array2003<int,2> >& get_segments() const;
  //list of triplets
  const std::vector<array2003<int,3> >& get_triplets() const;
  //list of neigbor chains [start,end] for segments
  //(use idx() to get real indexes)
  const std::vector<std::vector<array2003<int,2> > >& get_nbs() const;
  //convert idx from sorted representation to real idx
  int get_idx(int idx) const;
  //return list of indexes for conversion from sorted representation
  const std::vector<int>& get_idx_list() const;
  //convert idx from real idx to sorted representation
  int get_idxb(int idx) const;
  //return list of indexes for conversion to sorted representation
  const std::vector<int>& get_idxb_list() const;
  //check if the node is the end of the tube
  bool is_end(int idx) const;

  array2003<int, 2> get_segment(int idx) const;
  array2003<int, 3> get_triplet(int idx) const;

  static const int cnt_end = -1;
  static const int domain_end = -2;
  static const int not_cnt = -3;
private:
  std::vector<array2003<int, 2> > chain_list, segments;
  std::vector<array2003<int, 3> > triplets;
  std::vector<std::vector<array2003<int, 2> > > nb_chains;
  std::vector<int> index_list, index_list_b;
};

//=============================================================================

inline const std::vector<std::vector<array2003<int, 2> > > &
 MESONTList::get_nbs() const {
  return nb_chains;
}

inline int MESONTList::get_idx(int idx) const {
  return index_list[idx];
}

inline const std::vector<int>& MESONTList::get_idx_list() const {
  return index_list;
};


inline int MESONTList::get_idxb(int idx) const {
  return index_list_b[idx];
}

inline const std::vector<int>& MESONTList::get_idxb_list() const {
  return index_list_b;
};

inline const std::vector<array2003<int, 2> > & MESONTList::get_segments()
 const {
  return segments;
}

inline const std::vector<array2003<int, 3> > & MESONTList::get_triplets()
 const {
  return triplets;
}

inline array2003<int, 2> MESONTList::get_segment(int idx) const {
  array2003<int, 2> result;
  result[0] = chain_list[idx][0];
  result[1] = idx;
  return result;
}

inline array2003<int, 3> MESONTList::get_triplet(int idx) const {
  array2003<int, 3> result;
  result[0] = chain_list[idx][0];
  result[1] = idx;
  result[2] = chain_list[idx][1];
  return result;
}

inline bool MESONTList::is_end(int idx) const {
  return chain_list[idx][0] == cnt_end || chain_list[idx][1] == cnt_end;
};
