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

#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "pointers.h"

namespace LAMMPS_NS {

class Neighbor : protected Pointers {
 public:
  int style;                       // 0,1,2 = nsq, bin, multi
  int every;                       // build every this many steps
  int delay;                       // delay build for this many steps
  int dist_check;                  // 0 = always build, 1 = only if 1/2 dist
  int ago;                         // how many steps ago neighboring occurred
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  int includegroup;                // only build pairwise lists for this group

  double skin;                     // skin distance
  double cutneighmin;              // min neighbor cutoff for all type pairs
  double cutneighmax;              // max neighbor cutoff for all type pairs
  double *cuttype;                 // for each type, max neigh cut w/ others

  int ncalls;                      // # of times build has been called
  int ndanger;                     // # of dangerous builds

  int nrequest;                    // requests for pairwise neighbor lists
  class NeighRequest **requests;   // from Pair, Fix, Compute, Command classes
  int maxrequest;

  int old_style;                   // previous run info to avoid
  int old_nrequest;                // re-creation of pairwise neighbor lists
  class NeighRequest **old_requests;
  
  int nlist;                       // pairwise neighbor lists
  class NeighList **lists;

  int nbondlist;                   // list of bonds to compute
  int **bondlist;
  int nanglelist;                  // list of angles to compute
  int **anglelist;
  int ndihedrallist;               // list of dihedrals to compute
  int **dihedrallist;
  int nimproperlist;               // list of impropers to compute
  int **improperlist;

  Neighbor(class LAMMPS *);
  ~Neighbor();
  void init();
  int request(void *);         // another class requests a neighbor list
  void print_lists_of_lists(); // debug print out
  int decide();                // decide whether to build or not
  int check_distance();        // check max distance moved since last build
  void setup_bins();           // setup bins based on box and cutoff
  void build();                // create all neighbor lists (pair,bond)
  void build_one(int);         // create a single neighbor list
  void set(int, char **);      // set neighbor style and skin distance
  void modify_params(int, char**);  // modify parameters that control builds
  double memory_usage();
  
 private:
  int me,nprocs;

  int maxlocal;                    // size of atom-based NeighList arrays
  int maxbond,maxangle,maxdihedral,maximproper;   // size of bond lists

  int build_once;                  // 1 if only build lists once per run
  int must_check;                  // 1 if must check other classes to reneigh
  int restart_check;               // 1 if restart enabled, 0 if no
  int fix_check;                   // # of fixes that induce reneigh
  int *fixchecklist;               // which fixes to check

  double **cutneighsq;             // neighbor cutneigh sq for each type pair
  double cutneighmaxsq;            // cutneighmax squared
  double *cuttypesq;               // cuttype squared

  double triggersq;                // trigger = build when atom moves this dist

  double **xhold;                  // atom coords at last neighbor build
  int maxhold;                     // size of xhold array

  int nbinx,nbiny,nbinz;           // # of global bins
  int *bins;                       // ptr to next atom in each bin
  int maxbin;                      // size of bins array

  int *binhead;                    // ptr to 1st atom in each bin
  int maxhead;                     // size of binhead array

  int mbins;                       // # of local bins and offset
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  int binsizeflag;                 // user-chosen bin size
  double binsize_user;

  double binsizex,binsizey,binsizez;  // actual bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int sx,sy,sz,smax;               // bin stencil extents

  int dimension;                   // 2/3 for 2d/3d
  int triclinic;                   // 0 if domain is orthog, 1 if triclinic
  int newton_pair;                 // 0 if newton off, 1 if on for pairwise

  double *bboxlo,*bboxhi;          // copy of full domain bounding box

  double inner[2],middle[2];       // rRESPA cutoffs for extra lists
  double cut_inner_sq;		   // outer cutoff for inner neighbor list
  double cut_middle_sq;            // outer cutoff for middle neighbor list
  double cut_middle_inside_sq;     // inner cutoff for middle neighbor list

  int special_flag[4];             // flags for 1-2, 1-3, 1-4 neighbors

  int exclude;                     // 0 if no type/group exclusions, 1 if yes

  int nex_type;                    // # of entries in type exclusion list
  int maxex_type;                  // max # in type list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int maxex_group;                 // max # in group list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  int nex_mol;                     // # of entries in molecule exclusion list
  int maxex_mol;                   // max # in molecule list
  int *ex_mol_group;               // molecule group #'s to exclude
  int *ex_mol_bit;                 // molecule group bits to exclude

  int nblist,nglist,nslist;    // # of pairwise neigh lists of various kinds
  int *blist;                  // lists to build every reneighboring
  int *glist;                  // lists to grow atom arrays every reneigh
  int *slist;                  // lists to grow stencil arrays every reneigh

  void bin_atoms();                     // bin all atoms
  double bin_distance(int, int, int);   // distance between binx
  int coord2bin(double *);              // mapping atom coord to a bin
  int find_special(int, int);           // look for special neighbor
  int exclusion(int, int, int, int, int *, int *);  // test for pair exclusion
  void choose_build(int, class NeighRequest *);
  void choose_stencil(int, class NeighRequest *);

  // pairwise build functions

  typedef void (Neighbor::*PairPtr)(class NeighList *);
  PairPtr *pair_build;

  void half_nsq_no_newton(class NeighList *);
  void half_nsq_newton(class NeighList *);

  void half_bin_no_newton(class NeighList *);
  void half_bin_newton(class NeighList *);
  void half_bin_newton_tri(class NeighList *);

  void half_multi_no_newton(class NeighList *);
  void half_multi_newton(class NeighList *);
  void half_multi_newton_tri(class NeighList *);

  void full_nsq(class NeighList *);
  void full_bin(class NeighList *);
  void full_multi(class NeighList *);

  void half_from_full_no_newton(class NeighList *);
  void half_from_full_newton(class NeighList *);
  void skip_from(class NeighList *);
  void skip_from_granular(class NeighList *);
  void skip_from_respa(class NeighList *);
  void copy_from(class NeighList *);

  void granular_nsq_no_newton(class NeighList *);
  void granular_nsq_newton(class NeighList *);
  void granular_bin_no_newton(class NeighList *);
  void granular_bin_newton(class NeighList *);
  void granular_bin_newton_tri(class NeighList *);

  void respa_nsq_no_newton(class NeighList *);
  void respa_nsq_newton(class NeighList *);
  void respa_bin_no_newton(class NeighList *);
  void respa_bin_newton(class NeighList *);
  void respa_bin_newton_tri(class NeighList *);

  // pairwise stencil creation functions

  typedef void (Neighbor::*StencilPtr)(class NeighList *, int, int, int);
  StencilPtr *stencil_create;

  void stencil_half_bin_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_bin_3d_no_newton(class NeighList *, int, int, int);
  void stencil_half_bin_2d_newton(class NeighList *, int, int, int);
  void stencil_half_bin_3d_newton(class NeighList *, int, int, int);
  void stencil_half_bin_2d_newton_tri(class NeighList *, int, int, int);
  void stencil_half_bin_3d_newton_tri(class NeighList *, int, int, int);

  void stencil_half_multi_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_multi_3d_no_newton(class NeighList *, int, int, int);
  void stencil_half_multi_2d_newton(class NeighList *, int, int, int);
  void stencil_half_multi_3d_newton(class NeighList *, int, int, int);
  void stencil_half_multi_2d_newton_tri(class NeighList *, int, int, int);
  void stencil_half_multi_3d_newton_tri(class NeighList *, int, int, int);

  void stencil_full_bin_2d(class NeighList *, int, int, int);
  void stencil_full_bin_3d(class NeighList *, int, int, int);
  void stencil_full_multi_2d(class NeighList *, int, int, int);
  void stencil_full_multi_3d(class NeighList *, int, int, int);

  // topology build functions

  typedef void (Neighbor::*BondPtr)();   // ptrs to topology build functions

  BondPtr bond_build;                 // ptr to bond list functions
  void bond_all();                    // bond list with all bonds
  void bond_partial();                // exclude certain bonds

  BondPtr angle_build;                // ptr to angle list functions
  void angle_all();                   // angle list with all angles
  void angle_partial();               // exclude certain angles

  BondPtr dihedral_build;             // ptr to dihedral list functions
  void dihedral_all();                // dihedral list with all dihedrals
  void dihedral_partial();            // exclude certain dihedrals

  BondPtr improper_build;             // ptr to improper list functions
  void improper_all();                // improper list with all impropers
  void improper_partial();            // exclude certain impropers

};

}

#endif
