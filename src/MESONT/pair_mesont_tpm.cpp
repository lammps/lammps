// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "pair_mesont_tpm.h"
#include "export_mesont.h"


#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include <cstring>
#include <cmath>
#include <array>

#include <fstream>
#include <algorithm>

using namespace LAMMPS_NS;

class MESONTList {
public:
  MESONTList(const Atom* atom, const NeighList* nblist);
  ~MESONTList() {};
  //list of segments
  const std::vector<std::array<int,2>>& get_segments() const;
  //list of triplets
  const std::vector<std::array<int,3>>& get_triplets() const;
  //list of neighbor chains [start,end] for segments
  //(use idx() to get real indexes)
  const std::vector<std::vector<std::array<int,2>>>& get_nbs() const;
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

  std::array<int,2> get_segment(int idx) const;
  std::array<int,3> get_triplet(int idx) const;

  static const int cnt_end = -1;
  static const int domain_end = -2;
  static const int not_cnt = -3;
private:
  std::vector<std::array<int,2>> chain_list, segments;
  std::vector<std::array<int,3>> triplets;
  std::vector<std::vector<std::array<int,2>>> nb_chains;
  std::vector<int> index_list, index_list_b;
};

//=============================================================================

inline const std::vector<std::vector<std::array<int,2>>> &
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

inline const std::vector<std::array<int,2>> & MESONTList::get_segments()
 const {
  return segments;
}

inline const std::vector<std::array<int,3>> & MESONTList::get_triplets()
 const {
  return triplets;
}

inline std::array<int,2> MESONTList::get_segment(int idx) const {
  std::array<int,2> result;
  result[0] = chain_list[idx][0];
  result[1] = idx;
  return result;
}

inline std::array<int,3> MESONTList::get_triplet(int idx) const {
  std::array<int,3> result;
  result[0] = chain_list[idx][0];
  result[1] = idx;
  result[2] = chain_list[idx][1];
  return result;
}

inline bool MESONTList::is_end(int idx) const {
  return chain_list[idx][0] == cnt_end || chain_list[idx][1] == cnt_end;
};

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

MESONTList::MESONTList(const Atom* atom, const NeighList* nblist) {
  if (atom == nullptr || nblist == nullptr) return;
  //number of local atoms at the node
  int nlocal = atom->nlocal;
  //total number of atoms in the node and ghost shell treated as NTs
  int nall = nblist->inum + nblist->gnum;
  //total number of atoms in the node and ghost shell
  int ntot = atom->nlocal + atom->nghost;
  tagint* const g_id = atom->tag;
  tagint** const bonds = atom->bond_nt;
  tagint* const chain_id = atom->molecule;
  int* ilist = nblist->ilist;

  //convert bonds to local id representation
  chain_list.resize(ntot, {not_cnt,not_cnt});
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
      if (bonds[i][0] == g_id[nb]) {
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
     g_id[i] < g_id[chain_list[i][0]])
      segments.push_back({i,chain_list[i][0]});
    if (chain_list[i][1] != cnt_end && chain_list[i][1] != domain_end &&
     g_id[i] < g_id[chain_list[i][1]])
      segments.push_back({i,chain_list[i][1]});
  }
  int nbonds = segments.size();

  //triplets
  for (int i = 0; i < nlocal; i++) {
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
        //no self interactions for nbs within the same tube
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
        //if nodes are not continuous in the sorted representation
        //or represent chain ends, create a new neighbor chain
        int idx_next = chain_list[index_list[nb_list[j]]][1];
        if ((j == nnb - 1) || (nb_list[j] + 1 != nb_list[j+1]) ||
         (idx_next == cnt_end) || (idx_next == domain_end)) {
          std::array<int,2> chain;
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
          if (chain[0] != chain[1]) nb_chains[i].push_back(chain);
          idx_s = (j == nnb - 1) ? -1 : nb_list[j + 1];
        }
      }
    }
    nb_list.resize(0);
  }
}

/* ---------------------------------------------------------------------- */

// the cutoff distance between walls of tubes
static const double TPBRcutoff  = 3.0*3.4;
int PairMESONTTPM::instance_count = 0;
/* ---------------------------------------------------------------------- */

PairMESONTTPM::PairMESONTTPM(LAMMPS *lmp) : Pair(lmp) {
  writedata=1;
  BendingMode = 0;  // Harmonic bending model
  TPMType = 0;      // Inter-tube segment-segment interaction
  tab_path = nullptr;
  tab_path_length = 0;

  eatom_s = nullptr;
  eatom_b = nullptr;
  eatom_t = nullptr;
  nmax = 0;
  instance_count++;
  if (instance_count > 1) error->all(FLERR,
   "only a single instance of mesont/tpm pair style can be created");
}

/* ---------------------------------------------------------------------- */

PairMESONTTPM::~PairMESONTTPM()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);

    memory->destroy(eatom_s);
    memory->destroy(eatom_b);
    memory->destroy(eatom_t);
  }
  instance_count--;
  if (tab_path != nullptr) memory->destroy(tab_path);
}

/* ---------------------------------------------------------------------- */

void PairMESONTTPM::compute(int eflag, int vflag) {
  // set per atom values and accumulators
  // reallocate per-atom arrays if necessary
  ev_init(eflag,vflag);
  if (atom->nmax > nmax && eflag_atom) {
    memory->destroy(eatom_s);
    memory->create(eatom_s,comm->nthreads*maxeatom,"pair:eatom_s");
    memory->destroy(eatom_b);
    memory->create(eatom_b,comm->nthreads*maxeatom,"pair:eatom_b");
    memory->destroy(eatom_t);
    memory->create(eatom_t,comm->nthreads*maxeatom,"pair:eatom_t");
    nmax = atom->nmax;
  }
  //total number of atoms in the node and ghost shell treated as NTs
  int nall = list->inum + list->gnum;
  //total number of atoms in the node and ghost shell
  int ntot = atom->nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  if (!newton_pair)
   error->all(FLERR,"Pair style mesont/tpm requires newton pair on");

  double **x = atom->x;
  double **f = atom->f;
  double *r = atom->radius;
  double *l = atom->length;
  int *buckling = atom->buckling;
  tagint *g_id = atom->tag;

  //check if cutoff is chosen correctly
  double RT = mesont_lib_get_R();
  double Lmax = 0.0;
  for (int ii = 0; ii < list->inum; ii++) {
    int i = list->ilist[ii];
    if (Lmax < l[i]) Lmax = l[i];
  }
  double Rcut_min = std::max(2.0*Lmax, std::sqrt(0.5*Lmax*Lmax + std::pow((2.0*RT + TPBRcutoff),2)));
  if (cut_global < Rcut_min) {
    error->all(FLERR, "The selected cutoff is too small for the current system : "
               "L_max = {:.8}, R_max = {:.8}, Rc = {:.8}, Rcut_min = {:.8}",
               Lmax, RT, cut_global, Rcut_min);
  }

  //generate bonds and chain nblist
  MESONTList ntlist(atom, list);

  //reorder data to make it contiguous within tubes
  //and compatible with Fortran functions
  std::vector<double> x_sort(3*nall), f_sort(3*nall), s_sort(9*nall);
  std::vector<double> u_ts_sort(nall), u_tb_sort(nall), u_tt_sort(nall);
  std::vector<int> b_sort(nall);
  for (int i = 0; i < nall; i++) {
    int idx = ntlist.get_idx(i);
    for (int j = 0; j < 3; j++) x_sort[3*i+j] = x[idx][j];
    b_sort[i] = buckling[idx];
  }

  //bending potential
  int n_triplets = ntlist.get_triplets().size();
  for (int i = 0; i < n_triplets; i++) {
    const std::array<int,3>& t = ntlist.get_triplets()[i];
    //idx of nodes of a triplet in sorted representation
    int idx_s0 = ntlist.get_idxb(t[0]);
    int idx_s1 = ntlist.get_idxb(t[1]);
    int idx_s2 = ntlist.get_idxb(t[2]);

    double* X1 = &(x_sort[3*idx_s0]);
    double* X2 = &(x_sort[3*idx_s1]);
    double* X3 = &(x_sort[3*idx_s2]);
    double& U1b = u_tb_sort[idx_s0];
    double& U2b = u_tb_sort[idx_s1];
    double& U3b = u_tb_sort[idx_s2];
    double* F1 = &(f_sort[3*idx_s0]);
    double* F2 = &(f_sort[3*idx_s1]);
    double* F3 = &(f_sort[3*idx_s2]);
    double* S1 = &(s_sort[9*idx_s0]);
    double* S2 = &(s_sort[9*idx_s1]);
    double* S3 = &(s_sort[9*idx_s2]);
    double& R123 = r[t[1]];
    double& L123 = l[t[1]];
    int& BBF2 = b_sort[idx_s1];

    mesont_lib_TubeBendingForceField(U1b, U2b, U3b, F1, F2, F3, S1, S2, S3,
     X1, X2, X3, R123, L123, BBF2);
  }

  //share new values of buckling
  if (BendingMode == 1) {
    for (int i = 0; i < nall; i++) {
      int idx = ntlist.get_idx(i);
      buckling[idx] = b_sort[i];
    }
    comm->forward_comm_pair(this);
    for (int i = 0; i < nall; i++) {
      int idx = ntlist.get_idx(i);
      b_sort[i] = buckling[idx];
    }
  }

  //segment-segment and segment-tube interactions
  int n_segments = ntlist.get_segments().size();
  double Rmax = 0.0;
  Lmax = 0.0;
  for (int i = 0; i < n_segments; i++) {
    const std::array<int,2>& s = ntlist.get_segments()[i];
    //idx of a segment end 1 in sorted representation
    int idx_s0 = ntlist.get_idxb(s[0]);
    //idx of a segment end 2 in sorted representation
    int idx_s1 = ntlist.get_idxb(s[1]);
    double* X1 = &(x_sort[3*idx_s0]);
    double* X2 = &(x_sort[3*idx_s1]);
    double length = std::sqrt(std::pow(X1[0]-X2[0],2) +
     std::pow(X1[1]-X2[1],2) + std::pow(X1[2]-X2[2],2));
    if (length > Lmax) Lmax = length;
    double& U1t = u_tt_sort[idx_s0];
    double& U2t = u_tt_sort[idx_s1];
    double& U1s = u_ts_sort[idx_s0];
    double& U2s = u_ts_sort[idx_s1];
    double* F1 = &(f_sort[3*idx_s0]);
    double* F2 = &(f_sort[3*idx_s1]);
    double* S1 = &(s_sort[9*idx_s0]);
    double* S2 = &(s_sort[9*idx_s1]);
    double R12 = r[s[0]]; if (R12 > Rmax) Rmax = R12;
    if (std::abs(R12 - RT) > 1e-3)
        error->all(FLERR,"Inconsistent input and potential table");
    //assume that the length of the segment is defined by the node with
    //smallest global id
    double L12 = (g_id[s[0]] > g_id[s[1]]) ? l[s[1]] : l[s[0]];
    mesont_lib_TubeStretchingForceField(U1s, U2s, F1, F2, S1, S2, X1, X2,
     R12, L12);

    for (int nc = 0; nc < (int)ntlist.get_nbs()[i].size(); nc++) {
      //id of the beginning and end of the chain in the sorted representation
      const std::array<int,2>& chain = ntlist.get_nbs()[i][nc];
      int N = chain[1] - chain[0] + 1;  //number of elements in the chain
      int end1 = ntlist.get_idx(chain[0]);  //chain ends (real representation)
      int end2 = ntlist.get_idx(chain[1]);
      double* X = &(x_sort[3*chain[0]]);
      double* Ut = &(u_tt_sort[chain[0]]);
      double* F = &(f_sort[3*chain[0]]);
      double* S = &(s_sort[9*chain[0]]);
      double R = r[end1];
      int* BBF = &(b_sort[chain[0]]);
      int E1 = ntlist.is_end(end1);
      int E2 = ntlist.is_end(end2);

      int Ee = 0;
      double* Xe = X; double* Fe = F; double* Se = S;
      if (!E1 && ntlist.get_triplet(end1)[0] != MESONTList::domain_end &&
       ntlist.get_triplet(ntlist.get_triplet(end1)[0])[0] ==
       MESONTList::cnt_end) {
        Ee = 1;
        int idx = ntlist.get_idxb(ntlist.get_triplet(end1)[0]);
        Xe = &(x_sort[3*idx]);
        Fe = &(f_sort[3*idx]);
        Se = &(s_sort[9*idx]);
      }
      else if (!E2 && ntlist.get_triplet(end2)[2] != MESONTList::domain_end &&
       ntlist.get_triplet(ntlist.get_triplet(end2)[2])[2] ==
       MESONTList::cnt_end) {
        Ee = 2;
        int idx = ntlist.get_idxb(ntlist.get_triplet(end2)[2]);
        Xe = &(x_sort[3*idx]);
        Fe = &(f_sort[3*idx]);
        Se = &(s_sort[9*idx]);
      }

      mesont_lib_SegmentTubeForceField(U1t, U2t, Ut, F1, F2, F, Fe, S1, S2, S,
       Se, X1, X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType);
    }
  }

  //check if cutoff is chosen correctly
  Rcut_min = std::max(2.0*Lmax, std::sqrt(0.5*Lmax*Lmax + std::pow((2.0*Rmax + TPBRcutoff),2)));
  if (cut_global < Rcut_min) {
    error->all(FLERR, "The selected cutoff is too small for the current system : "
               "L_max = {:.8}, R_max = {:.8}, Rc = {:.8}, Rcut_min = {:.8}",
               Lmax, RT, cut_global, Rcut_min);
  }

  //convert from sorted representation
  for (int i = 0; i < nall; i++) {
      int idx = ntlist.get_idx(i);
      for (int j = 0; j < 3; j++) f[idx][j] += f_sort[3*i+j];
      buckling[idx] = b_sort[i];
  }
  if (eflag_global) {
    energy_s = energy_b = energy_t = 0.0;
    for (int i = 0; i < nall; i++) {
      energy_s += u_ts_sort[i];
      energy_b += u_tb_sort[i];
      energy_t += u_tt_sort[i];
    }
    eng_vdwl += energy_s + energy_b + energy_t;
  }
  if (eflag_atom) {
    for (int i = 0; i < ntot; i++)
      eatom_s[i] = eatom_b[i] = eatom_t[i] = 0.0;

    for (int i = 0; i < nall; i++) {
      int idx = ntlist.get_idx(i);
      eatom_s[idx] += u_ts_sort[i];
      eatom_b[idx] += u_tb_sort[i];
      eatom_t[idx] += u_tt_sort[i];
      eatom[idx] += u_ts_sort[i] + u_tb_sort[i] + u_tt_sort[i];
    }
  }
  if (vflag_global) {
    for (int i = 0; i < nall; i++) {
      virial[0] += s_sort[9*i+0]; //xx
      virial[1] += s_sort[9*i+4]; //yy
      virial[2] += s_sort[9*i+8]; //zz
      virial[3] += s_sort[9*i+1]; //xy
      virial[4] += s_sort[9*i+2]; //xz
      virial[5] += s_sort[9*i+5]; //yz
    }
  }
  if (vflag_atom) {
    for (int i = 0; i < nall; i++) {
      int idx = ntlist.get_idx(i);
      vatom[idx][0] += s_sort[9*i+0]; //xx
      vatom[idx][1] += s_sort[9*i+4]; //yy
      vatom[idx][2] += s_sort[9*i+8]; //zz
      vatom[idx][3] += s_sort[9*i+1]; //xy
      vatom[idx][4] += s_sort[9*i+2]; //xz
      vatom[idx][5] += s_sort[9*i+5]; //yz
    }
  }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMESONTTPM::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMESONTTPM::settings(int narg, char **arg) {
  if ((narg == 0) || (narg > 4))
    error->all(FLERR,"Illegal pair_style command");
  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        cut[i][j] = cut_global;
  }
  std::string TPMAFile = (narg > 1) ? arg[1] : "MESONT-TABTP.xrs";
  tab_path_length = TPMAFile.length();
  if (tab_path != nullptr) memory->destroy(tab_path);
  //c_str returns '\0' terminated string
  memory->create(tab_path,tab_path_length+1,"pair:path");
  std::memcpy(tab_path, TPMAFile.c_str(), tab_path_length+1);
  mesont_lib_SetTablePath(tab_path, tab_path_length);

  if (narg > 2) {
    BendingMode = utils::numeric(FLERR,arg[2],false,lmp);
    if ((BendingMode < 0) || (BendingMode > 1))
      error->all(FLERR,"Incorrect BendingMode");
  }
  if (narg > 3) {
    TPMType = utils::numeric(FLERR,arg[3],false,lmp);
    if ((TPMType < 0) || (TPMType > 1))
      error->all(FLERR,"Incorrect TPMType");
  }

  mesont_lib_TPBInit();
  int M, N;
  std::ifstream in(TPMAFile);
  if (!in.is_open()) error->all(FLERR,"Incorrect table path");
  std::string tmp;
  std::getline(in,tmp);
  std::getline(in,tmp);
  std::getline(in,tmp);
  in >> M >> N;
  in.close();
  mesont_lib_TPMInit(M, N);
  mesont_lib_InitCNTPotModule(1, 3, 0, BendingMode, mesont_lib_get_R());
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMESONTTPM::coeff(int narg, char **arg) {
  if ((narg < 2) || (narg > 3))
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double cut_one = cut_global;
  if (narg == 3) cut_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMESONTTPM::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMESONTTPM::write_restart(FILE *fp) {
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMESONTTPM::read_restart(FILE *fp) {
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMESONTTPM::write_restart_settings(FILE *fp) {
  fwrite(&BendingMode,sizeof(int),1,fp);
  fwrite(&TPMType,sizeof(int),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&tab_path_length,sizeof(int),1,fp);
  fwrite(tab_path,tab_path_length+1,1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMESONTTPM::read_restart_settings(FILE *fp) {
  int me = comm->me;
  if (me == 0) {
    fread(&BendingMode,sizeof(int),1,fp);
    fread(&TPMType,sizeof(int),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&tab_path_length,sizeof(int),1,fp);
  }
  MPI_Bcast(&BendingMode,1,MPI_INT,0,world);
  MPI_Bcast(&TPMType,1,MPI_INT,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tab_path_length,1,MPI_INT,0,world);

  if (tab_path != nullptr) memory->destroy(tab_path);
  memory->create(tab_path,tab_path_length+1,"pair:path");
  if (me == 0) fread(tab_path,tab_path_length+1,1,fp);
  MPI_Bcast(tab_path,tab_path_length+1,MPI_CHAR,0,world);
  mesont_lib_SetTablePath(tab_path,tab_path_length);
  mesont_lib_TPBInit();
  int M, N;
  std::ifstream in(tab_path);
  if (!in.is_open()) error->all(FLERR,"Incorrect table path");
  std::string tmp;
  std::getline(in,tmp);
  std::getline(in,tmp);
  std::getline(in,tmp);
  in >> M >> N;
  in.close();
  mesont_lib_TPMInit(M, N);
  mesont_lib_InitCNTPotModule(1, 3, 0, BendingMode, mesont_lib_get_R());
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMESONTTPM::write_data(FILE *fp) {
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\n",i);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMESONTTPM::write_data_all(FILE *fp) {
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}

/* ---------------------------------------------------------------------- */

void PairMESONTTPM::init_style() {
  //make sure that a full list is created (including ghost nodes)
  int r = neighbor->request(this,instance_me);
  neighbor->requests[r]->half = false;
  neighbor->requests[r]->full = true;
  neighbor->requests[r]->ghost = true;
}

void* PairMESONTTPM::extract(const char *str, int &) {
  if (strcmp(str,"mesonttpm_Es_tot") == 0) return &energy_s;
  else if (strcmp(str,"mesonttpm_Eb_tot") == 0) return &energy_b;
  else if (strcmp(str,"mesonttpm_Et_tot") == 0) return &energy_t;
  else if (strcmp(str,"mesonttpm_Es") == 0) return eatom_s;
  else if (strcmp(str,"mesonttpm_Eb") == 0) return eatom_b;
  else if (strcmp(str,"mesonttpm_Et") == 0) return eatom_t;
  else return nullptr;
};
