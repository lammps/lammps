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

#include "mesont_list.h"
#include <algorithm>
#include <cmath>

template<typename T>
void vector_union(std::vector<T>& v1, std::vector<T>& v2,
 std::vector<T>& merged) {
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  merged.reserve(v1.size() + v2.size());
  typename std::vector<T>::iterator it1 = v1.begin();
  typename std::vector<T>::iterator it2 = v2.begin();

  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 < *it2) {
      if (merged.empty() || merged.back() < *it1) merged.push_back(*it1);
        ++it1;
    }
    else {
      if (merged.empty() || merged.back() < *it2) merged.push_back(*it2);
      ++it2;
    }
  }
  while (it1 != v1.end()) {
    if (merged.empty() || merged.back() < *it1) merged.push_back(*it1);
    ++it1;
  }

  while (it2 != v2.end()) {
  if (merged.empty() || merged.back() < *it2) merged.push_back(*it2);
    ++it2;
  }
}

#ifndef NULL
#define NULL 0
#endif

MESONTList::MESONTList(const Atom* atom, const NeighList* nblist, double rc2){
  if (atom == NULL || nblist == NULL) return;
  //number of local atoms at the node
  int nlocal = atom->nlocal;
  //total number of atoms in the node and ghost shell
  int nall = nblist->inum + nblist->gnum;
  int ntot = atom->nlocal + atom->nghost;
  tagint* const g_id = atom->tag;
  tagint** const bonds = atom->bond_cnt;
  tagint* const chain_id = atom->molecule;
  int* ilist = nblist->ilist;

  //convert bonds to local id representation
  array2003<int, 2> tmp_arr;
  tmp_arr[0] = not_cnt; tmp_arr[1] = not_cnt;
  chain_list.resize(ntot, tmp_arr);
  for (int ii = 0; ii < nall; ii++) {
    int i = ilist[ii];
    chain_list[i][0] = domain_end;
    chain_list[i][1] = domain_end;
  }
  for (int ii = 0; ii < nall; ii++) {
    int i = ilist[ii];
    int nnb = nblist->numneigh[i];
    for (int m = 0; m < 2; m++)
      if (bonds[i][m] == cnt_end) chain_list[i][m] = cnt_end;
    for (int j = 0; j < nnb; j++) {
      int nb = nblist->firstneigh[i][j];
      if (bonds[i][0] == g_id[nb]){
        chain_list[i][0] = nb;
        chain_list[nb][1] = i;
        break;
      }
    }
  }

  //reorder chains: index list
  //list of indexes for conversion FROM reordered representation
  index_list.reserve(nall);
  index_list_b.resize(ntot, -1); // convert index TO reordered representation
  for (int i = 0; i < ntot; i++) {
    if (chain_list[i][0] == cnt_end || chain_list[i][0] == domain_end) {
      index_list.push_back(i);
      index_list_b[i] = index_list.size() - 1;
      int idx = i;
      while (1) {
        idx = chain_list[idx][1];
        if (idx == cnt_end || idx == domain_end) break;
        else index_list.push_back(idx);
        index_list_b[idx] = index_list.size() - 1;
      }
    }
  }

  //segment list
  for (int i = 0; i < nlocal; i++) {
    if (chain_list[i][0] == not_cnt) continue;
    if (chain_list[i][0] != cnt_end && chain_list[i][0] != domain_end &&
     g_id[i] < g_id[chain_list[i][0]]){
      array2003<int, 2> tmp_c;
      tmp_c[0] = i; tmp_c[1] = chain_list[i][0];
      segments.push_back(tmp_c);
    }
    if (chain_list[i][1] != cnt_end && chain_list[i][1] != domain_end &&
     g_id[i] < g_id[chain_list[i][1]]){
      array2003<int, 2> tmp_c;
       tmp_c[0] = i; tmp_c[1] = chain_list[i][1];
       segments.push_back(tmp_c);
    }
  }
  int nbonds = segments.size();

  //triplets
  for (int i = 0; i < nlocal; i++){
    if (chain_list[i][0] == not_cnt) continue;
    if (chain_list[i][0] != cnt_end && chain_list[i][0] != domain_end &&
     chain_list[i][1] != cnt_end && chain_list[i][1] != domain_end)
      triplets.push_back(get_triplet(i));
  }

  //segment neighbor list
  nb_chains.resize(nbonds);
  std::vector<int> nb_list_i[2], nb_list;
  for (int i = 0; i < nbonds; i++) {
    //union of nb lists
    for (int m = 0; m < 2; m++) {
      nb_list_i[m].resize(0);
      int idx = segments[i][m];
      if (idx >= nlocal) continue;
      int nnb = nblist->numneigh[idx];
      for (int j = 0; j < nnb; j++) {
        int jdx = nblist->firstneigh[idx][j];
        //no selfinteractions for nbs within the same tube
        if (chain_id[jdx] == chain_id[idx] &&
         std::abs(index_list_b[idx] - index_list_b[jdx]) <= 5) continue;
        nb_list_i[m].push_back(index_list_b[jdx]);
      }
    }
    vector_union(nb_list_i[0], nb_list_i[1], nb_list);

    int nnb = nb_list.size();
    if (nnb > 0) {
      int idx_s = nb_list[0];
      for (int j = 0; j < nnb; j++) {
        //if nodes are not continous in the sorted representation
        //or represent chain ends, create a new neighbor chain
        int idx_next = chain_list[index_list[nb_list[j]]][1];
        if ((j == nnb - 1) || (nb_list[j] + 1 != nb_list[j+1]) ||
         (idx_next == cnt_end) || (idx_next == domain_end)) {
          int idx_f = nb_list[j];
          array2003<int, 2> chain;
          chain[0] = idx_s;
          chain[1] = nb_list[j];
          //make sure that segments having at least one node
          //in the neighbor list are included
          int idx0 = index_list[chain[0]]; // real id of the ends
          int idx1 = index_list[chain[1]];
          if (chain_list[idx0][0] != cnt_end &&
           chain_list[idx0][0] != domain_end) chain[0] -= 1;
          if (chain_list[idx1][1] != cnt_end &&
           chain_list[idx1][1] != domain_end) chain[1] += 1;
          if(chain[0] != chain[1]) nb_chains[i].push_back(chain);
          idx_s = (j == nnb - 1) ? -1 : nb_list[j + 1];
        }
      }
    }
    nb_list.resize(0);
  }
}
