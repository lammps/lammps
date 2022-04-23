/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tine Curk (tcurk5@gmail.com) and Jiaxing Yuan (yuanjiaxing123@hotmail.com)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(charge/regulation,FixChargeRegulation);
// clang-format on
#else

#ifndef LMP_FIX_CHARGE_REGULATION_H
#define LMP_FIX_CHARGE_REGULATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChargeRegulation : public Fix {
 public:
  FixChargeRegulation(class LAMMPS *, int, char **);
  ~FixChargeRegulation() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  void forward_acid();
  void backward_acid();
  void forward_base();
  void backward_base();
  void forward_ions();
  void forward_ions_multival();
  void backward_ions();
  void backward_ions_multival();
  int get_random_particle(int, double, double, double *);
  int insert_particle(int, double, double, double *);
  double energy_full();
  int particle_number(int, double);
  int particle_number_xrd(int, double, double, double *);
  double compute_vector(int n) override;
  void assign_tags();
  void options(int, char **);
  void setThermoTemperaturePointer();
  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  int exclusion_group, exclusion_group_bit;
  int nevery, seed;    // begin MC cycle every nevery MD timesteps, random seed
  int nmc;             // MC move attempts per cycle
  double
      llength_unit_in_nm;    // LAMMPS unit of length in nm, needed since chemical potentials are in units of mol/l
  double pH, pKa, pKb, pKs, pI_plus,
      pI_minus;    // chemical potentials and equilibrium constant in log10 base
  double c10pH, c10pKa, c10pKb, c10pOH, c10pI_plus,
      c10pI_minus;    // 10 raised to chemical potential value, in units of concentration [mol/liter]
  double pmcmoves[3];    // mc move attempt probability: acid, base, ion pair exchange
  double pmcc;           // mc move cumulative attempt probability
  int npart_xrd;         // # of particles (ions) within xrd
  int npart_xrd2;        // # of particles (ions) within xrd
  double vlocal_xrd;     // # local volume within xrd
  bool
      only_salt_flag;    // true if performing only salt insertion/deletion, no acid/base dissociation.
  bool add_tags_flag;     // true if each inserted atom gets its unique atom tag
  int groupbitall;        // group bitmask for inserted atoms
  int ngroups;            // number of group-ids for inserted atoms
  char **groupstrings;    // list of group-ids for inserted atoms

  // counters
  unsigned long int nacid_attempts, nacid_successes, nbase_attempts, nbase_successes,
      nsalt_attempts, nsalt_successes;
  int nacid_neutral, nacid_charged, nbase_neutral, nbase_charged, ncation,
      nanion;     // particle type counts
  int cr_nmax;    //  max number of local particles
  double reservoir_temperature;
  double beta, sigma, volume,
      volume_rx;            // inverse temperature, speed, total volume, reacting volume
  int salt_charge[2];       // charge of salt ions: [0] - cation, [1] - anion
  int salt_charge_ratio;    // charge ratio when using multivalent ion exchange
  double xlo, xhi, ylo, yhi, zlo, zhi;    // box size
  double energy_stored;                   // full energy of old/current configuration
  int triclinic;                          // 0 = orthog box, 1 = triclinic
  double *sublo, *subhi;                  // triclinic size
  int *ptype_ID;                          // particle ID array
  double overlap_cutoffsq;                // square distance cutoff for overlap
  int overlap_flag;
  int acid_type, cation_type, base_type, anion_type;    // reacting atom types
  int reaction_distance_flag;                           // radial reaction restriction flag
  double reaction_distance;    // max radial distance from acid/base for ion insertion
  int pHvar, pHstyle;          // variable pH style
  char *pHstr;                 // variable pH input parsing

  class Compute *c_pe;               // energy compute pointer
  class RanPark *random_equal;       // random number generator
  class RanPark *random_unequal;     // random number generator
  char *idftemp;                     // pointer to the temperature fix
  double *target_temperature_tcp;    // current temperature of the thermostat
};
}    // namespace LAMMPS_NS

#endif
#endif
