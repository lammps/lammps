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

#ifdef FIX_CLASS

FixStyle(shake,FixShake)

#else

#ifndef LMP_FIX_SHAKE_H
#define LMP_FIX_SHAKE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixShake : public Fix {

 friend class FixEHEX;

 public:
  FixShake(class LAMMPS *, int, char **);
  virtual ~FixShake();
  virtual int setmask();
  virtual void init();
  void setup(int);
  void pre_neighbor();
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);

  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  void set_arrays(int);
  virtual void update_arrays(int, int);
  void set_molecule(int, tagint, int, double *, double *, double *);

  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(int, double *);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

  virtual void shake_end_of_step(int vflag);
  virtual void correct_coordinates(int vflag);
  virtual void correct_velocities();

  int dof(int);
  virtual void reset_dt();
  void *extract(const char *, int &);

 protected:
  int vflag_post_force;                  // store the vflag of last post_force call
  int respa;                             // 0 = vel. Verlet, 1 = respa
  int me,nprocs;
  int rattle;                            // 0 = SHAKE, 1 = RATTLE
  double tolerance;                      // SHAKE tolerance
  int max_iter;                          // max # of SHAKE iterations
  int output_every;                      // SHAKE stat output every so often
  bigint next_output;                    // timestep for next output

                                         // settings from input command
  int *bond_flag,*angle_flag;            // bond/angle types to constrain
  int *type_flag;                        // constrain bonds to these types
  double *mass_list;                     // constrain bonds to these masses
  int nmass;                             // # of masses in mass_list

  int molecular;                         // copy of atom->molecular
  double *bond_distance,*angle_distance; // constraint distances

  int ifix_respa;                        // rRESPA fix needed by SHAKE
  int nlevels_respa;                     // copies of needed rRESPA variables
  int *loop_respa;
  double *step_respa;

  double **x,**v,**f;                    // local ptrs to atom class quantities
  double **ftmp,**vtmp;                  // pointers to temporary arrays for f,v

  double *mass,*rmass;
  int *type;
  int nlocal;
                                         // atom-based arrays
  int *shake_flag;                       // 0 if atom not in SHAKE cluster
                                         // 1 = size 3 angle cluster
                                         // 2,3,4 = size of bond-only cluster
  tagint **shake_atom;                   // global IDs of atoms in cluster
                                         // central atom is 1st
                                         // lowest global ID is 1st for size 2
  int **shake_type;                      // bondtype of each bond in cluster
                                         // for angle cluster, 3rd value
                                         //   is angletype
  double **xshake;                       // unconstrained atom coords
  int *nshake;                           // count

  double dtv,dtfsq;                     // timesteps for trial move
  double dtf_inner,dtf_innerhalf;       // timesteps for rRESPA trial move

  int *list;                            // list of clusters to SHAKE
  int nlist,maxlist;                    // size and max-size of list

                                        // stat quantities
  int *b_count,*b_count_all;            // counts for each bond type
  double *b_ave,*b_max,*b_min;          // ave/max/min dist for each bond type
  double *b_ave_all,*b_max_all,*b_min_all;   // MPI summing arrays
  int *a_count,*a_count_all;            // ditto for angle types
  double *a_ave,*a_max,*a_min;
  double *a_ave_all,*a_max_all,*a_min_all;

  class Molecule **atommols;            // atom style template pointer
  class Molecule **onemols;             // molecule added on-the-fly
  int nmol;

  void find_clusters();
  void atom_owners();
  void partner_info(int *, tagint **, int **, int **, int **, int **);
  void nshake_info(int *, tagint **, int **);
  void shake_info(int *, tagint **, int **);
  
  int masscheck(double);
  void unconstrained_update();
  void unconstrained_update_respa(int);
  void shake(int);
  void shake3(int);
  void shake4(int);
  void shake3angle(int);
  void stats();
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
    tagint atomID,partnerID;
    int mask,type,massflag,bondtype;
  };

  struct NShakeInfo {
    tagint atomID,partnerID;
    int nshake;
  };

  struct ShakeInfo {
    tagint atomID;
    int shake_flag;
    int shake_atom[4];
    int shake_type[3];
  };

  // callback functions for rendezvous communication

  static int rendezvous_ids(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_partners_info(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_nshake(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_shake(int, char *, int &, int *&, char *&, void *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use fix shake with non-molecular system

Your choice of atom style does not have bonds.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid bond type index for fix shake

Self-explanatory.  Check the fix shake command in the input script.

E: Invalid angle type index for fix shake

Self-explanatory.

E: Invalid atom type index for fix shake

Atom types must range from 1 to Ntypes inclusive.

E: Invalid atom mass for fix shake

Mass specified in fix shake command must be > 0.0.

E: Too many masses for fix shake

The fix shake command cannot list more masses than there are atom
types.

E: Molecule template ID for fix shake does not exist

Self-explanatory.

W: Molecule template for fix shake has multiple molecules

The fix shake command will only recognize molecules of a single
type, i.e. the first molecule in the template.

E: Fix shake molecule template must have shake info

The defined molecule does not specify SHAKE information.

E: More than one fix shake

Only one fix shake can be defined.

E: Fix shake cannot be used with minimization

Cannot use fix shake while doing an energy minimization since
it turns off bonds that should contribute to the energy.

E: Shake fix must come before NPT/NPH fix

NPT fix must be defined in input script after SHAKE fix, else the
SHAKE fix contribution to the pressure virial is incorrect.

E: Bond potential must be defined for SHAKE

Cannot use fix shake unless bond potential is defined.

E: Angle potential must be defined for SHAKE

When shaking angles, an angle_style potential must be used.

E: Shake angles have different bond types

All 3-atom angle-constrained SHAKE clusters specified by the fix shake
command that are the same angle type, must also have the same bond
types for the 2 bonds in the angle.

E: Shake atoms %d %d missing on proc %d at step %ld

The 2 atoms in a single shake cluster specified by the fix shake
command are not all accessible to a processor.  This probably means
an atom has moved too far.

E: Shake atoms %d %d %d missing on proc %d at step %ld

The 3 atoms in a single shake cluster specified by the fix shake
command are not all accessible to a processor.  This probably means
an atom has moved too far.

E: Shake atoms %d %d %d %d missing on proc %d at step %ld

The 4 atoms in a single shake cluster specified by the fix shake
command are not all accessible to a processor.  This probably means
an atom has moved too far.

E: Did not find fix shake partner info

Could not find bond partners implied by fix shake command.  This error
can be triggered if the delete_bonds command was used before fix
shake, and it removed bonds without resetting the 1-2, 1-3, 1-4
weighting list via the special keyword.

E: Shake cluster of more than 4 atoms

A single cluster specified by the fix shake command can have no more
than 4 atoms.

E: Shake clusters are connected

A single cluster specified by the fix shake command must have a single
central atom with up to 3 other atoms bonded to it.

W: Shake determinant < 0.0

The determinant of the quadratic equation being solved for a single
cluster specified by the fix shake command is numerically suspect.  LAMMPS
will set it to 0.0 and continue.

E: Shake determinant = 0.0

The determinant of the matrix being solved for a single cluster
specified by the fix shake command is numerically invalid.

*/
