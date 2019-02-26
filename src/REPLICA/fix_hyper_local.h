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

  double memory_usage();

  // extra methods visible to callers

  void init_hyper();
  void build_bond_list(int);

 private:
  int me;
  double cutbond,qfactor,vmax,tequil,dcut;
  double alpha_user;         // timescale to apply boostostat (time units)
  double alpha;              // unitless dt/alpha_user
  double boost_target;       // target value of boost
  int checkbias,checkbias_every,checkbias_flag,checkbias_count;

  int setupflag;             // 1 during setup, 0 during run
  int firstflag;             // set for first time bond_build takes place
  int nostrainyet;           // 1 until maxstrain is first computed

  int nbias_running,nobias_running;
  int nbondbuild;
  double time_bondbuild;
  bigint starttime;
  double sumbiascoeff;   // sum of aveboost at every timestep
  bigint allbonds;       // sum of bond count on this step
  double allbias;        // sum of biascoeff on all bonds on this step

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
  double maxbiascoeff;       // cummulative max bias coeff for any bond
  double minbiascoeff;       // cummulative min bias coeff for any bond
  double rmaxever,rmaxeverbig;
  int ghost_toofar;

  class NeighList *listfull;   // full neigh list up to Dcut distance
  class NeighList *listhalf;   // half neigh list up to pair distance
                               // both created only when bonds are rebuilt

  // list of my owned bonds
  // persists on a proc from one event until the next

  struct OneBond {             // single IJ bond, atom I is owner
    int i,j;                   // current local indices of 2 bond atoms
    int iold,jold;             // local indices when bonds were formed
    double r0;                 // relaxed bond length
    double biascoeff;          // biasing coefficient = prefactor Cij
  };

  struct OneBond *blist;       // list of owned bonds
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

  // data struct used to persist biascoeffs when bond list is re-created

  struct OneCoeff {
    double biascoeff;
    tagint jtag;
  };

  struct OneCoeff **clist;     // list of bond coeffs for each atom's bonds
  int *numcoeff;               // # of coeffs per atom
  int maxcoeff;                // allocate size of clist
  int maxcoeffperatom;         // allocated # of columns in clist

  // list of biased bonds this proc owns

  int maxbias;                 // allocated size of bias list
  int nbias;                   // # of biased bonds I own
  int *bias;                   // index of biased bonds in my bond list

  // extra timers

  //double timefirst,timesecond,timethird,timefourth;
  //double timefifth,timesixth,timeseventh,timetotal;

  // private methods

  void grow_bond();
  void grow_coeff();
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
