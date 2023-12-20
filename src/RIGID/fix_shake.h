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

#ifdef FIX_CLASS
// clang-format off
FixStyle(shake,FixShake);
// clang-format on
#else

#ifndef LMP_FIX_SHAKE_H
#define LMP_FIX_SHAKE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixShake : public Fix {

  friend class FixEHEX;

 public:
  FixShake(class LAMMPS *, int, char **);
  ~FixShake() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void setup_pre_reverse(int, int) override;
  void min_setup(int) override;
  void pre_neighbor() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_pre_reverse(int, int) override;
  void min_post_force(int) override;
  void post_run() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void update_arrays(int, int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  virtual void shake_end_of_step(int vflag);
  virtual void correct_coordinates(int vflag);
  virtual void correct_velocities();

  int dof(int) override;
  void reset_dt() override;
  void *extract(const char *, int &) override;
  double compute_scalar() override;

 protected:
  int vflag_post_force;     // store the vflag of last post_force call
  int eflag_pre_reverse;    // store the eflag of last pre_reverse call
  int respa;                // 0 = vel. Verlet, 1 = respa
  int rattle;               // 0 = SHAKE, 1 = RATTLE
  double tolerance;         // SHAKE tolerance
  int max_iter;             // max # of SHAKE iterations
  int output_every;         // SHAKE stat output every so often
  bigint next_output;       // timestep for next output

  // settings from input command
  int *bond_flag, *angle_flag;    // bond/angle types to constrain
  int *type_flag;                 // constrain bonds to these types
  double *mass_list;              // constrain bonds to these masses
  int nmass;                      // # of masses in mass_list

  int molecular;                             // copy of atom->molecular
  double *bond_distance, *angle_distance;    // constraint distances
  double kbond;                              // force constant for restraint
  double ebond;                              // energy of bond restraints

  class FixRespa *fix_respa;    // rRESPA fix needed by SHAKE
  int nlevels_respa;            // copies of needed rRESPA variables
  int *loop_respa;
  double *step_respa;

  double **x, **v, **f;     // local ptrs to atom class quantities
  double **ftmp, **vtmp;    // pointers to temporary arrays for f,v

  double *mass, *rmass;
  int *type;
  int nlocal;
  // atom-based arrays
  int *shake_flag;        // 0 if atom not in SHAKE cluster
                          // 1 = size 3 angle cluster
                          // 2,3,4 = size of bond-only cluster
  tagint **shake_atom;    // global IDs of atoms in cluster
                          // central atom is 1st
                          // lowest global ID is 1st for size 2
  int **shake_type;       // bondtype of each bond in cluster
                          // for angle cluster, 3rd value
                          //   is angletype
  double **xshake;        // unconstrained atom coords
  int *nshake;            // count

  double dtv, dtfsq;                  // timesteps for trial move
  double dtf_inner, dtf_innerhalf;    // timesteps for rRESPA trial move

  int *list;             // list of clusters to SHAKE
  int **closest_list;    // list of closest atom indices in SHAKE clusters
  int nlist, maxlist;    // size and max-size of list

  // stat quantities
  int *b_count, *b_count_all;                   // counts for each bond type, atoms in bond cluster
  double *b_ave, *b_max, *b_min;                // ave/max/min dist for each bond type
  double *b_ave_all, *b_max_all, *b_min_all;    // MPI summing arrays
  int *a_count, *a_count_all;                   // ditto for angle types
  double *a_ave, *a_max, *a_min;
  double *a_ave_all, *a_max_all, *a_min_all;

  class Molecule **atommols;    // atom style template pointer
  class Molecule **onemols;     // molecule added on-the-fly
  int nmol;

  void find_clusters();
  void atom_owners();
  void partner_info(int *, tagint **, int **, int **, int **, int **);
  void nshake_info(int *, tagint **, int **);
  void shake_info(int *, tagint **, int **);

  int masscheck(double);
  virtual void unconstrained_update();
  void unconstrained_update_respa(int);
  void shake(int);
  void shake3(int);
  void shake4(int);
  void shake3angle(int);
  void bond_force(int, int, double);
  virtual void stats();
  int bondtype_findset(int, tagint, tagint, int);
  int angletype_findset(int, tagint, tagint, int);

  // data used by rendezvous callback methods

  int nrvous;
  tagint *atomIDs;
  int *procowner;

  struct IDRvous {
    int me;
    tagint atomID;
  };

  struct PartnerInfo {
    tagint atomID, partnerID;
    int mask, type, massflag, bondtype;
  };

  struct NShakeInfo {
    tagint atomID, partnerID;
    int nshake;
  };

  struct ShakeInfo {
    tagint atomID;
    tagint shake_atom[4];
    int shake_flag;
    int shake_type[3];
  };

  // callback functions for rendezvous communication

  static int rendezvous_ids(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_partners_info(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_nshake(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_shake(int, char *, int &, int *&, char *&, void *);
};

}    // namespace LAMMPS_NS

#endif
#endif
