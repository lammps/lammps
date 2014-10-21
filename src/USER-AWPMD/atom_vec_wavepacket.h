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

/* ----------------------------------------------------------------------
   Contributing author: Ilya Valuev (JIHT RAS)
------------------------------------------------------------------------- */


#ifdef ATOM_CLASS

AtomStyle(wavepacket,AtomVecWavepacket)

#else

#ifndef LMP_ATOM_VEC_WAVEPACKET_H
#define LMP_ATOM_VEC_WAVEPACKET_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecWavepacket : public AtomVec {
public:
  AtomVecWavepacket(class LAMMPS *);
  ~AtomVecWavepacket() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  void force_clear(int, size_t);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  int pack_comm_hybrid(int, int *, double *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int unpack_comm_hybrid(int, int, double *);
  int pack_reverse(int, int, double *);
  int pack_reverse_hybrid(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int unpack_reverse_hybrid(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void data_vel(int, char **);
  int data_vel_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  void pack_vel(double **);
  int pack_vel_hybrid(int, double *);
  void write_vel(FILE *, int, double **);
  int write_vel_hybrid(FILE *, double *);
  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);
  bigint memory_usage();

private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;

  ///\en spin: -1 or 1 for electron, 0 for ion (compatible with eff)
  int *spin;
  ///\en charge: must be specified in the corresponding units (-1 for electron in real units, eff compatible)
  double *q;
  ///\en width of the wavepacket (compatible with eff)
  double *eradius;
  ///\en width velocity for the wavepacket (compatible with eff)
  double *ervel;
  ///\en (generalized) force on width  (compatible with eff)
  double *erforce;

  // AWPMD- specific:
  ///\en electron tag: must be the same for the WPs belonging to the same electron
  int *etag;
  ///\en wavepacket split coeffcients: cre, cim, size is 2*N
  double *cs;
  ///\en force on wavepacket split coeffcients: re, im, size is 2*N
  double *csforce;
  ///\en (generalized) force on velocity, size is 3*N
  double *vforce;
   ///\en (generalized) force on radius velocity, size is N
  double *ervelforce;
};

}

#endif
#endif
