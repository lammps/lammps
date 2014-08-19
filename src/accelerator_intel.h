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

// NOTE: this file is *supposed* to be included multiple times

#ifdef LMP_USER_INTEL

// true interface to USER-INTEL

// this part is used inside the neighbor.h header file to
// add functions to the Neighbor class definition

#ifdef LMP_INSIDE_NEIGHBOR_H
  #ifdef LMP_INTEL_OFFLOAD
  #ifdef __INTEL_OFFLOAD
  template <class flt_t, class acc_t> friend class IntelBuffers;
  inline int * special_flag_alloc() { return special_flag; }
  #endif
  #endif

  friend class FixIntel;
  void *fix_intel;

  template <class flt_t, class acc_t>
  void bin_atoms(void *);

template <class flt_t, class acc_t, int>
  void hbni(const int, NeighList *, void *, const int, const int, void *,
	    const int offload_end = 0);
  template <class flt_t, class acc_t>
  void hbnni(const int, NeighList *, void *, const int, const int, void *);
  template <class flt_t, class acc_t, int>
  void hbnti(const int, NeighList *, void *, const int, const int, void *,
	     const int offload_end = 0);

  void half_bin_no_newton_intel(class NeighList *);
  void half_bin_newton_intel(class NeighList *);
  void half_bin_newton_tri_intel(class NeighList *);

#endif /* !LMP_INSIDE_NEIGHBOR_H */

#else /* !LMP_USER_INTEL */

// needed for compiling Neighbor class when USER-Intel is not installed

#ifdef LMP_INSIDE_NEIGHBOR_H

  void half_bin_no_newton_intel(class NeighList *) {}
  void half_bin_newton_intel(class NeighList *) {}
  void half_bin_newton_tri_intel(class NeighList *) {}

#endif

#endif /* !LMP_USER_INTEL */
