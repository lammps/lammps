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

FixStyle(hyper/local,FixHyperLocal)

#else

#ifndef LMP_FIX_HYPER_LOCAL_H
#define LMP_FIX_HYPER_LOCAL_H

#include "fix_hyper.h"

namespace LAMMPS_NS {

class FixHyperLocal : public FixHyper {
 public:
  FixHyperLocal(class LAMMPS *, int, char **);
  ~FixHyperLocal();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup_pre_neighbor();
  void setup_pre_reverse(int, int);
  void pre_neighbor();
  void pre_reverse(int, int);
  void min_pre_neighbor();
  double compute_scalar();
  double compute_vector(int);
  double query(int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double memory_usage();

  // extra methods visible to callers

  void init_hyper();
  void build_bond_list(int);

 private:
  int me;
  double cutbond,qfactor,vmax,tequil,dcut;
  double alpha_user;         // timescale to apply boostostat (time units)
  double alpha;              // unitless dt/alpha_user
  double boosttarget;        // target value of boost
  int histoflag;
  int lostbond,lostbond_partner;
  double lostbond_coeff;
  int checkbias,checkbias_every,checkbias_flag,checkbias_count;
  int checkcoeff,checkcoeff_every,checkcoeff_flag,checkcoeff_count;

  int setupflag;             // 1 during setup, 0 during run
  int firstflag;             // set for first time bond_build takes place
  int nostrainyet;           // 1 until maxstrain is first computed

  int nboost_running,nobias_running;
  int nbondbuild;
  double time_bondbuild;
  bigint starttime;
  double sumboostcoeff;  // sum of aveboost at every timestep
  int allbonds;          // sum of bond count on this step
  double allboost;       // sum of boostcoeff on all bonds on this step

  int nnewbond;              // running tally of number of new bonds created
  int maxbondperatom;        // max # of bonds any atom ever has
  int commflag;              // flag for communication mode
  int nevent;                // # of events that trigger bond rebuild
  int nevent_atom;           // # of atoms that experienced an event
  double cutbondsq,dcutsq;
  double beta,invqfactorsq;
  double mybias;
  double maxbondlen;         // cummulative max length of any bond
  double maxdriftsq;         // max distance any atom drifts from original pos
  double maxboostcoeff;      // cummulative max boost coeff for any bond
  double minboostcoeff;      // cummulative min boost coeff for any bond
  double rmaxever,rmaxeverbig;
  int ghost_toofar;

  // extra timers

  //double timefirst,timesecond,timethird,timefourth;
  //double timefifth,timesixth,timeseventh,timetotal;

  // data structs for per-atom and per-bond info
  // all of these are for current owned and ghost atoms
  // except list and old2now are for atom indices at time of last bond build

  class NeighList *list;       // full neigh list up to Dcut distance
                               // created only when bonds are rebuilt

  int *old2now;                // o2n[i] = current local index of old atom i
                               //   stored for old owned and ghost atoms
                               //   I = old index when bonds were last created
                               //   old indices are stored in old neighbor list

  double **xold;               // coords of owned+ghost atoms when bonds created
  tagint *tagold;              // global IDs of owned+ghost atoms when b created

  int maxold;                  // allocated size of old2now
  int maxbond;                 // allocated size of bonds
  int old_nall;                // nlocal+nghost when old2now was last setup

  struct OneBond {             // single IJ bond, atom I is owner
    double r0;                 // original relaxed bond length
    double boostcoeff;         // boost coefficient
    tagint jtag;               // global index of J atom in bond IJ
    int j;                     // local index of J atom in bond IJ
  };

  struct OneBond **bonds;      // 2d array of bonds for owned atoms
  int *numbond;                // number of bonds for each owned atom

  double *maxstrain;           // max-strain of any bond atom I is part of
                               //   for owned and ghost atoms
  double *maxstrain_region;    // max-strain of any neighbor atom J of atom I
                               //   for owned and ghost atoms
  int *maxstrain_bondindex;    // index of max-strain bond of each atom I
                               //   just for owned atoms
  tagint *biasflag;            // atoms in biased bonds marked with bond partner
                               //   for owned and ghost atoms

  // list of boosted bonds that this proc will bias

  int maxboost;                // allocated size of boost list
  int nboost;                  // # of boosted bonds I own
  int *boost;                  // index of atom I in each boosted bond

  // histogramming of bond boost cooeficients

  int histo_every,histo_count,histo_print,histo_steps;
  double histo_delta,invhisto_delta,histo_lo;
  bigint *histo,*allhisto;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
