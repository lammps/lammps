/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NEIGHBOR_H
#define LMP_NEIGHBOR_H

#include "pointers.h"

namespace LAMMPS_NS {

// forward declarations
class NeighRequest;
class NeighList;

class Neighbor : protected Pointers {
 public:
  enum { NSQ, BIN, MULTI_OLD, MULTI };
  int style;           // 0,1,2,3 = nsq, bin, multi/old, multi
  int every;           // build every this many steps
  int delay;           // delay build for this many steps
  int dist_check;      // 0 = always build, 1 = only if 1/2 dist
  int ago;             // how many steps ago neighboring occurred
  int pgsize;          // size of neighbor page
  int oneatom;         // max # of neighbors for one atom
  int includegroup;    // only build pairwise lists for this group
  int build_once;      // 1 if only build lists once per run

  double skin;                    // skin distance
  double cutneighmin;             // min neighbor cutoff for all type pairs
  double cutneighmax;             // max neighbor cutoff for all type pairs
  double cutneighmaxsq;           // cutneighmax squared
  double **cutneighsq;            // neighbor cutneigh sq for each type pair
  double **cutneighghostsq;       // cutneigh sq for each ghost type pair
  double *cuttype;                // for each type, max neigh cut w/ others
  double *cuttypesq;              // cuttype squared
  double cut_inner_sq;            // outer cutoff for inner neighbor list
  double cut_middle_sq;           // outer cutoff for middle neighbor list
  double cut_middle_inside_sq;    // inner cutoff for middle neighbor list

  int binsizeflag;        // user-chosen bin size
  double binsize_user;    // set externally by some accelerator pkgs

  bigint ncalls;      // # of times build has been called
  bigint ndanger;     // # of dangerous builds
  bigint lastcall;    // timestep of last neighbor::build() call

  // geometry and static info, used by other Neigh classes

  double *bboxlo, *bboxhi;    // ptrs to full domain bounding box
                              // different for orthog vs triclinic

  // exclusion info, used by NeighPair

  int exclude;    // 0 if no type/group exclusions, 1 if yes

  int nex_type;                // # of entries in type exclusion list
  int *ex1_type, *ex2_type;    // pairs of types to exclude
  int **ex_type;               // 2d array of excluded type pairs

  int nex_group;                 // # of entries in group exclusion list
  int *ex1_group, *ex2_group;    // pairs of group #'s to exclude
  int *ex1_bit, *ex2_bit;        // pairs of group bits to exclude

  int nex_mol;          // # of entries in molecule exclusion list
  int *ex_mol_group;    // molecule group #'s to exclude
  int *ex_mol_bit;      // molecule group bits to exclude
  int *ex_mol_intra;    // 0 = exclude if in 2 molecules (inter)
                        // 1 = exclude if in same molecule (intra)

  // special info, used by NeighPair

  int special_flag[4];    // flags for 1-2, 1-3, 1-4 neighbors

  // cluster setting, used by NeighTopo

  int cluster_check;    // 1 if check bond/angle/etc satisfies minimg

  // pairwise neighbor lists and corresponding requests

  int nlist;           // # of pairwise neighbor lists
  int nrequest;        // # of requests, same as nlist
  int old_nrequest;    // # of requests for previous run

  NeighList **lists;
  NeighRequest **requests;        // from Pair,Fix,Compute,Command classes
  NeighRequest **old_requests;    // copy of requests to compare to
  int *j_sorted;                  // index of requests sorted by cutoff distance

  // data from topology neighbor lists

  int nbondlist;    // list of bonds to compute
  int **bondlist;
  int nanglelist;    // list of angles to compute
  int **anglelist;
  int ndihedrallist;    // list of dihedrals to compute
  int **dihedrallist;
  int nimproperlist;    // list of impropers to compute
  int **improperlist;

  // optional type grouping for multi

  int custom_collection_flag;      // 1 if custom collections are defined for multi
  int interval_collection_flag;    // 1 if custom collections use intervals
  int finite_cut_flag;             // 1 if multi considers finite atom size
  int ncollections;                // # of custom collections
  int nmax_collection;             // maximum atoms stored in collection array
  int *type2collection;            // ntype array mapping types to custom collections
  double *collection2cut;          // ncollection array with upper bounds on cutoff intervals
  double **cutcollectionsq;        // cutoffs for each combination of collections
  int *collection;                 // local per-atom array to store collection id

  // public methods

  Neighbor(class LAMMPS *);
  ~Neighbor() override;
  virtual void init();

  // old API for creating neighbor list requests
  int request(void *, int instance = 0);

  // new API for creating neighbor list requests
  NeighRequest *add_request(class Pair *, int flags = 0);
  NeighRequest *add_request(class Fix *, int flags = 0);
  NeighRequest *add_request(class Compute *, int flags = 0);
  NeighRequest *add_request(class Command *, const char *, int flags = 0);

  // set neighbor list request OpenMP flag
  void set_omp_neighbor(int);

  // report if we have INTEL package neighbor lists
  bool has_intel_request() const;

  int decide();                     // decide whether to build or not
  virtual int check_distance();     // check max distance moved since last build
  void setup_bins();                // setup bins based on box and cutoff
  virtual void build(int);          // build all perpetual neighbor lists
  virtual void build_topology();    // pairwise topology neighbor lists
  // create a one-time pairwise neigh list
  void build_one(class NeighList *list, int preflag = 0);
  void set(int, char **);                     // set neighbor style and skin distance
  void reset_timestep(bigint);                // reset of timestep counter
  void modify_params(int, char **);           // modify params that control builds
  void modify_params(const std::string &);    // convenience overload

  void exclusion_group_group_delete(int, int);    // rm a group-group exclusion
  int exclude_setting();                          // return exclude value to accelerator pkg

  // find a neighbor list based on requestor
  NeighList *find_list(void *, const int id = 0) const;
  // find a neighbor request based on requestor
  NeighRequest *find_request(void *, const int id = 0) const;

  const std::vector<NeighRequest *> get_pair_requests() const;
  int any_full();                // Check if any old requests had full neighbor lists
  void build_collection(int);    // build peratom collection array starting at the given index

  bigint get_nneigh_full();    // return number of neighbors in a regular full neighbor list
  bigint get_nneigh_half();    // return number of neighbors in a regular half neighbor list
  double memory_usage();

  bigint last_setup_bins;    // step of last neighbor::setup_bins() call

 protected:
  int me, nprocs;
  int firsttime;    // flag for calling init_styles() only once

  int dimension;      // 2/3 for 2d/3d
  int triclinic;      // 0 if domain is orthog, 1 if triclinic
  int newton_pair;    // 0 if newton off for pairwise, 1 if on

  int must_check;       // 1 if must check other classes to reneigh
  int restart_check;    // 1 if restart enabled, 0 if no
  int fix_check;        // # of fixes that induce reneigh
  int *fixchecklist;    // which fixes to check

  double triggersq;    // trigger = build when atom moves this dist

  double **xhold;    // atom coords at last neighbor build
  int maxhold;       // size of xhold array

  int boxcheck;                           // 1 if need to store box size
  double boxlo_hold[3], boxhi_hold[3];    // box size at last neighbor build
  double corners_hold[8][3];              // box corners at last neighbor build
  double (*corners)[3];                   // ptr to 8 corners of triclinic box

  double inner[2], middle[2];    // rRESPA cutoffs for extra lists

  int old_style, old_triclinic;    // previous run info
  int old_pgsize, old_oneatom;     // used to avoid re-creating neigh lists

  int nstencil_perpetual;    // # of perpetual NeighStencil classes
  int npair_perpetual;       // #x of perpetual NeighPair classes
  int *slist;                // indices of them in neigh_stencil
  int *plist;                // indices of them in neigh_pair

  int maxex_type;     // max # in exclusion type list
  int maxex_group;    // max # in exclusion group list
  int maxex_mol;      // max # in exclusion molecule list

  int maxatom;       // max size of atom-based NeighList arrays
  int maxrequest;    // max size of NeighRequest list

  // info for other Neigh classes: NBin,NStencil,NPair,NTopo

  int nbin, nstencil;
  int nbclass, nsclass, npclass;
  int bondwhich, anglewhich, dihedralwhich, improperwhich;

  typedef class NBin *(*BinCreator)(class LAMMPS *);
  BinCreator *binclass;
  char **binnames;
  int *binmasks;
  class NBin **neigh_bin;

  typedef class NStencil *(*StencilCreator)(class LAMMPS *);
  StencilCreator *stencilclass;
  char **stencilnames;
  int *stencilmasks;
  class NStencil **neigh_stencil;

  typedef class NPair *(*PairCreator)(class LAMMPS *);
  PairCreator *pairclass;
  char **pairnames;
  int *pairmasks;
  class NPair **neigh_pair;

  class NTopo *neigh_bond;
  class NTopo *neigh_angle;
  class NTopo *neigh_dihedral;
  class NTopo *neigh_improper;

  // internal methods
  // including creator methods for Nbin,Nstencil,Npair instances

  void init_styles();
  int init_pair();
  virtual void init_topology();

  void sort_requests();

  void morph_unique();
  void morph_skip();
  void morph_granular();
  void morph_halffull();
  void morph_copy_trim();

  void print_pairwise_info();
  void requests_new2old();

  int choose_bin(class NeighRequest *);
  int choose_stencil(class NeighRequest *);
  int choose_pair(class NeighRequest *);

  // dummy functions provided by NeighborKokkos, called in init()
  // otherwise NeighborKokkos would have to overwrite init()

  int copymode;

  virtual void init_cutneighsq_kokkos(int) {}
  virtual void create_kokkos_list(int) {}
  virtual void init_ex_type_kokkos(int) {}
  virtual void init_ex_bit_kokkos() {}
  virtual void init_ex_mol_bit_kokkos() {}
  virtual void grow_ex_mol_intra_kokkos() {}
  virtual void set_binsize_kokkos() {}
};

namespace NeighConst {

  enum {
    NB_INTEL = 1 << 0,
    NB_KOKKOS_DEVICE = 1 << 1,
    NB_KOKKOS_HOST = 1 << 2,
    NB_SSA = 1 << 3,
    NB_STANDARD = 1 << 4,
    NB_MULTI = 1 << 5
  };

  enum {
    NS_BIN = 1 << 0,
    NS_MULTI = 1 << 1,
    NS_HALF = 1 << 2,
    NS_FULL = 1 << 3,
    NS_2D = 1 << 4,
    NS_3D = 1 << 5,
    NS_ORTHO = 1 << 6,
    NS_TRI = 1 << 7,
    NS_GHOST = 1 << 8,
    NS_SSA = 1 << 9,
    NS_MULTI_OLD = 1 << 10
  };

  enum {
    NP_NSQ = 1 << 0,
    NP_BIN = 1 << 1,
    NP_MULTI = 1 << 2,
    NP_HALF = 1 << 3,
    NP_FULL = 1 << 4,
    NP_ORTHO = 1 << 5,
    NP_TRI = 1 << 6,
    NP_ATOMONLY = 1 << 7,
    NP_MOLONLY = 1 << 8,
    NP_NEWTON = 1 << 9,
    NP_NEWTOFF = 1 << 10,
    NP_GHOST = 1 << 11,
    NP_SIZE = 1 << 12,
    NP_ONESIDE = 1 << 13,
    NP_RESPA = 1 << 14,
    NP_BOND = 1 << 15,
    NP_OMP = 1 << 16,
    NP_INTEL = 1 << 17,
    NP_KOKKOS_DEVICE = 1 << 18,
    NP_KOKKOS_HOST = 1 << 19,
    NP_SSA = 1 << 20,
    NP_COPY = 1 << 21,
    NP_SKIP = 1 << 22,
    NP_HALF_FULL = 1 << 23,
    NP_OFF2ON = 1 << 24,
    NP_MULTI_OLD = 1 << 25,
    NP_TRIM = 1 << 26
  };

  enum {
    REQ_DEFAULT = 0,
    REQ_FULL = 1 << 0,
    REQ_GHOST = 1 << 1,
    REQ_SIZE = 1 << 2,
    REQ_HISTORY = 1 << 3,
    REQ_OCCASIONAL = 1 << 4,
    REQ_RESPA_INOUT = 1 << 5,
    REQ_RESPA_ALL = 1 << 6,
    REQ_NEWTON_ON = 1 << 8,
    REQ_NEWTON_OFF = 1 << 9,
    REQ_SSA = 1 << 10,
  };
}    // namespace NeighConst

}    // namespace LAMMPS_NS

#endif
