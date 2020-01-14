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
#include "my_page.h"

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
  int pack_reverse_comm_size(int, int);
  void unpack_reverse_comm(int, int *, double *);

  double memory_usage();

  // extra methods visible to callers

  void init_hyper();
  void build_bond_list(int);

 private:
  int me;

  // inputs and derived quantities

  double cutbond,qfactor,vmax,tequil,dcut;
  double alpha_user;         // timescale to apply boostostat (time units)
  double alpha;              // unitless dt/alpha_user
  double boost_target;       // target value of boost
  int checkghost,checkbias;  // flags for optional stats
  int checkbias_every;
  int checkbias_flag;
  int boundflag,resetfreq;   // bias coeff bounding and reset settings
  double boundfrac;

  bigint groupatoms;         // # of atoms in fix group
  double cutbondsq,dcutsq;
  double beta,invvmax,invqfactorsq;

  // DEBUG - 2 lines
  int overcount;
  double myboost;

  // flags

  int setupflag;             // 1 during setup, 0 during run
  int firstflag;             // set for first time bond_build takes place
  int nostrainyet;           // 1 until maxstrain is first compute
  bigint starttime;          // timestep when this fix was invoked
  int commflag;              // flag for communication mode

  // bias coeff bounds and reset

  double bound_lower,bound_upper;
  bigint lastreset;

  // stats

  int nbondbuild;            // # of rebuilds of bond list
  double time_bondbuild;     // CPU time for bond builds

  bigint allbonds;           // current total # of bonds
  int nnewbond;              // running tally of # of new bonds created
  int maxbondperatom;        // max # of bonds any atom ever has
  int nevent;                // # of events that trigger bond rebuild
  int nevent_atom;           // # of atoms that experienced an event

  int nbias_running;         // running count of biased bonds
  int nobias_running;        // ditto for bonds with bias = 0, b/c too long
  int negstrain_running;     // ditto for bonds with negative strain

  double mybias;             // sum of bias potentials for biased bonds
  double maxbondlen;         // cummulative max length of any bond
  double maxdriftsq;         // max distance any bond atom drifts from quenched x

  double sumboost;              // sum of all bond boosts at each timestep
  double aveboost_running;      // cummulative sumboost/allbonds across steps
  double aveboost_running_output;      // most recent output of ab_running
  double sumbiascoeff;          // sum of all bond bias coeffs at each timestep
  double avebiascoeff_running;  // cummulative sumbiascoeff/allbonds across steps
  double avebiascoeff_running_output;  // most recent output of abc_running
  double minbiascoeff;          // min bias coeff on this step for my bonds
  double maxbiascoeff;          // max bias coeff on this step for my bonds
  double minbiascoeff_running;  // cummulative min bias coeff for any bond
  double maxbiascoeff_running;  // cummulative max bias coeff for any bond

  double rmaxever,rmaxeverbig;
  int ghost_toofar;          // # of ghost atoms not found in Dcut neigh list
  int checkbias_count;       // count of too-close biased bonds

  // 2 neighbor lists

  class NeighList *listfull;   // full neigh list up to Dcut distance
  class NeighList *listhalf;   // half neigh list up to pair distance
                               // both created only when bonds are rebuilt

  // list of my owned bonds and bias coeffs
  // persists on a proc from one event until the next

  struct OneBond {             // single IJ bond, atom I is owner
    int i,j;                   // current local indices of 2 bond atoms
    int iold,jold;             // local indices when bonds were formed
    double r0;                 // relaxed bond length
  };

  OneBond *blist;              // list of owned bonds
  double *biascoeff;           // bias coefficient Cij for each bond
  int nblocal;                 // # of owned bonds
  int maxbond;                 // allocated size of blist

  // old data from last timestep bonds were formed
  // persists on a proc from one event until the next
  // first set of vectors are maxlocal in length
  // second set of vectors are maxall in length

  int nlocal_old;               // nlocal for old atoms
  int nall_old;                 // nlocal+nghost for old atoms
  int maxlocal;                 // allocated size of old local atom vecs
  int maxall;                   // allocated size of old all atom vecs

  int *numbond;                 // # of bonds owned by old owned atoms
  int *maxhalf;                 // bond index for maxstrain bond of old atoms
  int *eligible;                // 0/1 flag for bias on one of old atom's bonds
  double *maxhalfstrain;        // strain value for maxstrain bond of old atoms

  int *old2now;                 // o2n[i] = current local index of old atom I
                                // may be -1 if ghost atom has drifted
  tagint *tagold;               // IDs of atoms when bonds were formed
                                // 0 if a ghost atom is not in Dcut neigh list
  double **xold;                // coords of atoms when bonds were formed

  // vectors used to find maxstrain bonds within a local domain

  int maxatom;                 // size of these vectors, nlocal + nghost

  double *maxstrain;           // max-strain of any bond atom I is part of
                               //   for owned and ghost atoms
  double *maxstrain_domain;    // max-strain of any neighbor atom J of atom I
                               //   for owned and ghost atoms
  tagint *biasflag;            // atoms in biased bonds marked with bond partner
                               //   for owned and ghost atoms

  // list of biased bonds this proc owns

  int maxbias;                 // allocated size of bias list
  int nbias;                   // # of biased bonds I own
  int *bias;                   // index of biased bonds in my bond list

  // data structs for persisting bias coeffs when bond list is reformed

  struct OneCoeff {
    double biascoeff;
    tagint tag;
  };

  MyPage<OneCoeff> *cpage;     // pages of OneCoeff datums for clist
  OneCoeff **clist;            // ptrs to vectors of bias coeffs for each atom
  int *numcoeff;               // # of bias coeffs per atom (one per bond)
  int maxcoeff;                // allocate sized of clist and numcoeff

  // extra timers

  //double timefirst,timesecond,timethird,timefourth;
  //double timefifth,timesixth,timeseventh,timetotal;

  // private methods

  void grow_bond();
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
