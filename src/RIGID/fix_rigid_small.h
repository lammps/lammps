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
FixStyle(rigid/small,FixRigidSmall);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_SMALL_H
#define LMP_FIX_RIGID_SMALL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRigidSmall : public Fix {
  friend class ComputeRigidLocal;

 public:
  FixRigidSmall(class LAMMPS *, int, char **);
  ~FixRigidSmall() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void write_restart_file(const char *) override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void setup_pre_neighbor() override;
  void pre_neighbor() override;
  int dof(int) override;
  void deform(int) override;
  void enforce2d() override;
  void reset_dt() override;
  void zero_momentum() override;
  void zero_rotation() override;
  int modify_param(int, char **) override;
  void *extract(const char *, int &) override;
  double extract_ke();
  double extract_erotational();
  double compute_scalar() override;
  double memory_usage() override;

 protected:
  int me, nprocs;
  double dtv, dtf, dtq;
  double *step_respa;
  int triclinic;

  char *inpfile;       // file to read rigid body attributes from
  int setupflag;       // 1 if body properties are setup, else 0
  int earlyflag;       // 1 if forces/torques are computed at post_force()
  int commflag;        // various modes of forward/reverse comm
  int customflag;      // 1 if custom property/variable define bodies
  int nbody;           // total # of rigid bodies
  int nlinear;         // total # of linear rigid bodies
  tagint maxmol;       // max mol-ID
  double maxextent;    // furthest distance from body owner to body atom

  struct Body {
    int natoms;            // total number of atoms in body
    int ilocal;            // index of owning atom
    double mass;           // total mass of body
    double xcm[3];         // COM position
    double xgc[3];         // geometric center position
    double vcm[3];         // COM velocity
    double fcm[3];         // force on COM
    double torque[3];      // torque around COM
    double quat[4];        // quaternion for orientation of body
    double inertia[3];     // 3 principal components of inertia
    double ex_space[3];    // principal axes in space coords
    double ey_space[3];
    double ez_space[3];
    double xgc_body[3];    // geometric center relative to xcm in body coords
    double angmom[3];      // space-frame angular momentum of body
    double omega[3];       // space-frame omega of body
    double conjqm[4];      // conjugate quaternion momentum
    int remapflag[4];      // PBC remap flags
    imageint image;        // image flags of xcm
    imageint dummy;        // dummy entry for better alignment
  };

  Body *body;         // list of rigid bodies, owned and ghost
  int nlocal_body;    // # of owned rigid bodies
  int nghost_body;    // # of ghost rigid bodies
  int nmax_body;      // max # of bodies that body can hold
  int bodysize;       // sizeof(Body) in doubles

  // per-atom quantities
  // only defined for owned atoms, except bodyown for own+ghost

  int *bodyown;          // index of body if atom owns a body, -1 if not
  tagint *bodytag;       // ID of body this atom is in, 0 if none
                         // ID = tag of atom that owns body
  int *atom2body;        // index of owned/ghost body this atom is in, -1 if not
                         // can point to original or any image of the body
  imageint *xcmimage;    // internal image flags for atoms in rigid bodies
                         // set relative to in-box xcm of each body
  double **displace;     // displacement of each atom in body coords
  int *eflags;           // flags for extended particles
  double **orient;       // orientation vector of particle wrt rigid body
  double **dorient;      // orientation of dipole mu wrt rigid body

  int extended;       // 1 if any particles have extended attributes
  int orientflag;     // 1 if particles store spatial orientation
  int dorientflag;    // 1 if particles store dipole orientation
  int reinitflag;     // 1 if re-initialize rigid bodies between runs

  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;

  // temporary per-body storage

  int **counts;        // counts of atom types in bodies
  double **itensor;    // 6 space-frame components of inertia tensor

  // mass per body, accessed by granular pair styles

  double *mass_body;
  int nmax_mass;

  // Langevin thermostatting

  int langflag;                        // 0/1 = no/yes Langevin thermostat
  double t_start, t_stop, t_period;    // thermostat params
  double **langextra;                  // Langevin thermostat forces and torques
  int maxlang;                         // max size of langextra
  class RanMars *random;               // RNG

  int tstat_flag, pstat_flag;    // 0/1 = no/yes thermostat/barostat

  int t_chain, t_iter, t_order;

  double p_start[3], p_stop[3];
  double p_period[3], p_freq[3];
  int p_flag[3];
  int pcouple, pstyle;
  int p_chain;

  int allremap;            // remap all atoms
  int dilate_group_bit;    // mask for dilation group
  char *id_dilate;         // group name to dilate

  char *id_gravity;    // ID of fix gravity command to add gravity forces
  double *gvec;        // ptr to gravity vector inside the fix

  double p_current[3], p_target[3];

  // molecules added on-the-fly as rigid bodies

  class Molecule **onemols;
  int nmol;

  // class data used by ring communication callbacks

  double rsqfar;

  struct InRvous {
    int me, ilocal;
    tagint atomID, bodyID;
    double x[3];
  };

  struct OutRvous {
    int ilocal;
    tagint atomID;
  };

  // local methods

  void image_shift();
  void set_xv();
  void set_v();
  void create_bodies(tagint *);
  void setup_bodies_static();
  void setup_bodies_dynamic();
  void apply_langevin_thermostat();
  void compute_forces_and_torques();
  void readfile(int, double **, int *);
  void grow_body();
  void reset_atom2body();

  // callback function for rendezvous communication

  static int rendezvous_body(int, char *, int &, int *&, char *&, void *);

  // debug

  //void check(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
