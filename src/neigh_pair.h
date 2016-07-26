/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NEIGH_PAIR_H
#define LMP_NEIGH_PAIR_H

#include "pointers.h"

namespace LAMMPS_NS {

class NeighPair : protected Pointers {
 public:
  int style;                       // ID for NeighPair method this is
  class NeighBin *nb;              // NeighBin instance I depend on
  class NeighStencil *ns;          // NeighStencil instance I depend on

  bigint last_build;            // timestep for last operations performed
  bigint last_copy_bin_setup;   // last timestep I invoked copy_bin_setup_info()
  bigint last_copy_bin;         // last step I invoked copy_bin_info()
  bigint last_copy_stencil;     // last step I invoked copy_bin_stencil_info()

  NeighPair(class LAMMPS *);
  virtual ~NeighPair() {}
  void copy_neighbor_info();
  void build_setup();
  virtual void build(class NeighList *) = 0;

 protected:

  // data from Neighbor class

  int includegroup;
  int exclude;
  double skin;
  double **cutneighsq;
  double **cutneighghostsq;
  double cut_inner_sq;
  double cut_middle_sq;
  double cut_middle_inside_sq;
  double *zeroes;
  double *bboxlo,*bboxhi;

  // exclusion data from Neighbor class

  int nex_type;                    // # of entries in type exclusion list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  int nex_mol;                     // # of entries in molecule exclusion list
  int *ex_mol_group;               // molecule group #'s to exclude
  int *ex_mol_bit;                 // molecule group bits to exclude

  // special data from Neighbor class

  int *special_flag;

  // data from NeighBin class

  int nbinx,nbiny,nbinz;
  int mbins;
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  double bininvx,bininvy,bininvz;
  int *bins;
  int *binhead;
  
  // data from NeighStencil class

  int nstencil;
  int *stencil;
  int **stencilxyz;
  int *nstencil_multi;
  int **stencil_multi;
  double **distsq_multi;

  // data common to all NeighPair variants

  int molecular;

  // methods for all NeighPair variants

  void copy_bin_setup_info();
  void copy_bin_info();
  void copy_stencil_info();

  int exclusion(int, int, int,
                int, int *, tagint *) const;   // test for pair exclusion
  int coord2bin(double *);                     // mapping atom coord to a bin
  int coord2bin(double *, int &, int &, int&); // ditto

  // find_special: determine if atom j is in special list of atom i
  // if it is not, return 0
  // if it is and special flag is 0 (both coeffs are 0.0), return -1
  // if it is and special flag is 1 (both coeffs are 1.0), return 0
  // if it is and special flag is 2 (otherwise), return 1,2,3
  //   for which level of neighbor it is (and which coeff it maps to)

  inline int find_special(const tagint *list, const int *nspecial,
                          const tagint tag) const {
    const int n1 = nspecial[0];
    const int n2 = nspecial[1];
    const int n3 = nspecial[2];

    for (int i = 0; i < n3; i++) {
      if (list[i] == tag) {
        if (i < n1) {
          if (special_flag[1] == 0) return -1;
          else if (special_flag[1] == 1) return 0;
          else return 1;
        } else if (i < n2) {
          if (special_flag[2] == 0) return -1;
          else if (special_flag[2] == 1) return 0;
          else return 2;
        } else {
          if (special_flag[3] == 0) return -1;
          else if (special_flag[3] == 1) return 0;
          else return 3;
        }
      }
    }
    return 0;
  };
};

}

#endif

/* ERROR/WARNING messages:

*/
