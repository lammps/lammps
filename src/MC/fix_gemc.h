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
FixStyle(gemc,FixGEMC);
// clang-format on
#else

#ifndef LMP_FIX_GEMC_H
#define LMP_FIX_GEMC_H

#include "fix.h"

#include <array>
#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

class FixGEMC : public Fix {
 public:
  FixGEMC(class LAMMPS *, int, char **);
  ~FixGEMC() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_exchange() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double compute_vector(int) override;

 private:
  // user provided inputs

  int ntranslate;             // number of translations attempted every nevery steps
  int nexchange;              // number of particle exchanges attempted every nevery steps
  int nvolume;                // number of volume exchanges attempted every nevery steps
  double box_temp;            // temperature of both boxes
  double displace;            // maximum displacement for translations
  double max_dlogvolratio;    // maximum change in log(V1/V2)
  double maxangle;            // maximum rotation angle of molecules (radians)
  int seed;                   // RNG seed
  int tune_every;             // adjust step sizes every this many invocations (0 = never)
  double tune_trans;          // target acceptance ratio of translations and rotations
  double tune_vol;            // target acceptance ratio of volume changes
  double tune_last[6];        // counters at the last adjustment
  int ntune_calls;            // number of invocations since the start of the run
  int full_flag;              // 1 if user requested full energy for all moves
  int local_flag;      // 1 if single-atom energies may be used for translations and exchanges
  int ghosts_stale;    // 1 if ghost atoms may be out of date
  int local_warned;    // 1 if the fallback to full energy at run time was reported

  // molecule exchange

  int molflag;                // 1 if whole molecules are moved and exchanged
  char *idmol;                // ID of molecule template
  class Molecule *onemol;     // molecule template
  int natoms_per_molecule;    // number of atoms in each molecule
  int group_charged;          // 1 if atoms in the fix group have charges (atoms only)
  int charge_warned;          // 1 if warning about charged exchanges was printed

  // for evaluating probability

  double beta;             // 1 / kT
  double energy_stored;    // current potential energy
  class Compute *c_pe;     // compute to get full potential energy

  // for determining which move to make

  int nmoves;            // total MC moves (translate + exchange + volume)
  double pc_exchange;    // cumulative probability MC move is an exchange
  double pc_volume;      // cumulative probability MC move is exchange or volume change

  // for tracking how many attempts/successes

  double ntranslation_attempts;
  double ntranslation_successes;
  double nrotation_attempts;
  double nrotation_successes;
  double nexchange_attempts;
  double nexchange_successes;
  double nvolume_attempts;
  double nvolume_successes;
  double nlast[6];       // counters at the time of the last progress message
  double logvolratio;    // log(V1/V2), identical in both boxes

  // particle - related props

  bigint natom_lower;                  // number of group atoms on lower ranks of my box
  int natom_local;                     // number of group atoms on this rank
  bigint natom_total;                  // number of group atoms in my box
  std::vector<int> gas_list;           // local indices of group atoms
  std::vector<int> gas_pos;            // position of a local index in gas_list or -1
  std::vector<double> exchange_buf;    // storage for atoms removed during a trial exchange

  // single-atom energies: owned and ghost atoms are sorted into a grid of bins at least
  // as large as the pair cutoff. accepted translations and exchanges update the grid and
  // all images of the atom on all ranks, including new images within the ghost cutoff,
  // instead of rebuilding the ghost atoms. atoms inserted this way and new images are
  // stored after the ghost atoms, and atoms removed this way stay in place until
  // flush_pending() applies the insertions and removals to the owned atoms.

  int grid_valid;                             // 1 if the grid matches the atom arrays
  int nstored;                                // number of owned, ghost, and pending atoms
  int nbin[3];                                // number of bins in each dimension
  double binlo[3];                            // lower corner of the grid
  double bininv[3];                           // inverse bin size in each dimension
  std::vector<std::vector<int>> bins;         // local indices of the atoms in each bin
  std::vector<int> atombin;                   // bin of each local index or -1
  std::vector<int> atombinpos;                // position of each local index in its bin
  std::unordered_map<tagint, int> taghead;    // first stored local index with an atom ID
  std::vector<int> tagnext;                   // next local index with the same atom ID
  int nbase;                      // number of owned and ghost atoms when the grid was built
  int use_map;                    // 1 if the atom map provides the images of owned and ghost atoms
  std::vector<int> image_list;    // scratch list of local indices
  std::vector<std::array<int, 3>> image_shifts;    // scratch list of periodic shifts
  std::vector<double> dacc;    // displacement of owned atoms since ghosts were built
  double ghost_skin;           // ghost cutoff minus pair cutoff
  double move_limit;           // largest displacement allowed by the subdomains
  tagint maxtag_box;           // largest atom ID in my box
  std::vector<int> pending;    // local indices of inserted atoms owned by this rank,
                               // not yet created, or -1 if removed again
  std::vector<int> removed;    // local indices of removed owned atoms, not yet deleted
  int pending_changes;         // 1 if atoms were inserted or removed since last flush

  // domain - related props

  int triclinic;                 // 1 if the box is triclinic
  double xlo, ylo, zlo;          // lower domain bounds
  double xhi, yhi, zhi;          // upper domain bounds
  double boxxy, boxxz, boxyz;    // tilt factors
  double voltot;                 // V1+V2, conserved

  // for communication

  int me;         // rank in my box
  int myworld;    // index of my box (0 or 1)

  MPI_Comm comm_replica;    // for communication between rank 0 of both boxes

  class RanPark *random_universe;    // RNG synchronized across all ranks of both boxes
  class RanPark *random_world;       // RNG synchronized across all ranks of one box
  class RanPark *random_proc;        // RNG unique to each rank

  int progress;    // last percentage of the run reported

  // rendezvous communication

  struct ComRvous {
    tagint mol;
    int proc;
    double sum[4];
  };
  struct CheckRvous {
    tagint mol, tag;
    int type;
    int ntopo[4];    // number of bonds, angles, dihedrals, impropers stored with the atom
    double q;
  };
  bigint rvous_count;    // number of molecules on this rank in the rendezvous decomposition
  int rvous_bad;         // 1 if a molecule on this rank does not match the template
  int rvous_charged;     // 1 if a molecule on this rank has a net charge
  static int rendezvous_centers(int, char *, int &, int *&, char *&, void *);
  static int rendezvous_check(int, char *, int &, int *&, char *&, void *);
  int group_topology();

  void attempt_atomic_translation_full();
  void attempt_volume_change_full();
  void attempt_atomic_exchange_full();
  void attempt_molecule_translation_full();
  void attempt_molecule_rotation_full();
  void attempt_molecule_exchange_full();

  int any_box(int);                // 1 if flag is set in either box
  int accept_both(double, int);    // joint acceptance decision of both boxes
  void reset_comm();               // re-distribute atoms and rebuild ghosts and neighbor lists
  void refresh_ghosts();           // re-distribute atoms and rebuild ghost atoms
  double energy_full();            // computes full potential energy
  double energy_local(int, int, tagint, double *);    // pair energy of one atom
  int use_local();        // 1 if the next move in my box may use energy_local()
  double local_skin();    // distance an atom may move before ghost atoms must be rebuilt
  double init_skin();     // estimate of local_skin() before the run is set up
  tagint insert_atom(int, int, double, double *, double *, int);    // insert atom into my box
  void update_gas_atoms_list();     // updates list of local group atoms
  void gas_add(int);                // add local index to gas_list
  void gas_remove(int);             // remove local index from gas_list
  int pick_random_gas_atom();       // picks random group atom
  tagint pick_random_molecule();    // picks random molecule of the group
  tagint first_atom(tagint);        // first atom ID of a molecule with consecutive IDs
  tagint gather_molecule(tagint, std::vector<double> &);    // gather atoms of one molecule

  void build_grid();                  // sort owned and ghost atoms into bins
  void bin_index(double *, int *);    // bin indices of a point
  int coord2bin(double *);            // bin of a point
  void bin_add(int);                  // add local index to its bin
  void bin_remove(int);               // remove local index from its bin
  void grow_stored(int);              // make room for local index
  int stored_atom(tagint, int, int, double, double *, double *);    // store pending atom
  void find_images(tagint, std::vector<int> &);                     // local indices of all images
  void remove_images(tagint);            // remove all images of an atom ID
  void move_images(tagint, double *);    // displace all images of an atom ID
  void insert_pending(tagint, int, int, double, double *, double *, int);
  void add_images(tagint, int, int, double, double *, double *);    // store missing images
  double subdomain_excess(double *);    // distance outside of my subdomain
  int need_refresh(int, tagint &);      // rebuild ghosts if an atom would move too far
  void flush_pending();                 // apply pending insertions and removals
  bigint molecule_centers(std::unordered_map<tagint, std::array<double, 3>> &);
  bigint scale_positions(double, int);    // scale molecule centers and atoms relative to box origin
  void set_box(double, double, double, double, double, double);    // change box size
  double box_volume();
  double min_box_width();
  void box_widths(double *);    // distances between opposite box faces
  double max_move();            // largest move of an atom allowed by the subdomain size
  double max_translation();     // displacement limited by the subdomain size
  int owns(double *);           // 1 if a point in the box is inside my subdomain
  int local_index(tagint);      // local index of owned atom with this ID or -1
  void random_point(double *);
  void changed_atoms();    // update after atoms or charges changed
  void check_molecules();
  void print_progress();
  void tune_steps();    // adjust maximum displacement and volume change
};

}    // namespace LAMMPS_NS

#endif
#endif
