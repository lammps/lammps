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
   Contributing authors: Ryan S. Elliott,
                         Valeriu Smirichinski,
                         Ellad Tadmor (U Minn)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-v1.6.0 (and newer) package
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(kim,PairKIM)

#else

#ifndef LMP_PAIR_KIM_H
#define LMP_PAIR_KIM_H

// includes from KIM & LAMMPS
class KIM_API_model;
#include "pair.h"
#include "KIM_SimulatorHeaders.hpp"


namespace LAMMPS_NS {

   class PairKIM : public Pair {
   public:
      PairKIM(class LAMMPS*);
      ~PairKIM();

      // LAMMPS Pair class virtual function prototypes
      virtual void compute(int, int);
      virtual void settings(int, char**);
      virtual void coeff(int, char**);
      virtual void init_style();
      virtual void init_list(int id, NeighList *ptr);
      virtual double init_one(int, int);
      virtual int pack_reverse_comm(int, int, double*);
      virtual void unpack_reverse_comm(int, int*, double*);
      virtual double memory_usage();

   private:
      // (nearly) all bool flags are not initialized in constructor, but set
      // explicitly in the indicated function.  All other data members are
      // initialized in constructor
      int settings_call_count;
      int init_style_call_count;

      // values set in settings()
      char* kim_modelname;

      // values set in coeff()

      // values set in allocate(), called by coeff()
      void allocate();
      int* lmps_map_species_to_unique;

      // values set in coeff(), after calling allocate()
      char** lmps_unique_elements;  // names of unique elements given
                                    // in pair_coeff command
      int lmps_num_unique_elements;

      // values set in set_lmps_flags(), called from init_style()
      bool lmps_using_newton;
      bool lmps_using_molecular;
      enum unit_sys {REAL, METAL, SI, CGS, ELECTRON};
      unit_sys lmps_units;
      KIM::LengthUnit lengthUnit;
      KIM::EnergyUnit energyUnit;
      KIM::ChargeUnit chargeUnit;
      KIM::TemperatureUnit temperatureUnit;
      KIM::TimeUnit timeUnit;


      KIM::Model * pkim;
      KIM::ComputeArguments * pargs;

      // values set in set_kim_model_has_flags(), called by kim_init()
      bool kim_model_has_energy;
      bool kim_model_has_forces;
      bool kim_model_has_virial;
      bool kim_model_has_particleEnergy;
      bool kim_model_has_particleVirial;

      // values set in kim_init()
      bool kim_init_ok;
      int lmps_local_tot_num_atoms;
      double kim_global_influence_distance;  // KIM Model cutoff value
      int kim_number_of_cutoffs;
      double const * kim_cutoff_values;
      class NeighList ** neighborLists;

      // values set in init_style()
      bool kim_model_init_ok;
      bool kim_particle_codes_ok;
      int *kim_particle_codes;

      // values set in compute()
      int lmps_maxalloc;              // max allocated memory value
      int* kim_particleSpecies;       // array of KIM particle species
      int* kim_particleContributing;  // array of KIM particle contributing
      int* lmps_stripped_neigh_list;  // neighbors of one atom, used when LAMMPS
                                      // is in molecular mode
      int** lmps_stripped_neigh_ptr;  // pointer into lists

      // KIM specific helper functions
      void kim_init();
      void kim_free();
      void set_argument_pointers();
      void set_lmps_flags();
      void set_kim_model_has_flags();
      // static methods used as callbacks from KIM
     static int get_neigh(
         void const * const dataObject,
         int const numberOfCutoffs, double const * const cutoffs,
         int const neighborListIndex, int const particleNumber,
         int * const numberOfNeighbors,
         int const ** const neighborsOfParticle);
   };
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unrecognized virial argument in pair_style command

Only two options are supported: LAMMPSvirial and KIMvirial

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Invalid args for non-hybrid pair coefficients

"NULL" is only supported in pair_coeff calls when using pair hybrid

E: PairKIM only works with 3D problems

This is a current limitation.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: KIM neighbor iterator exceeded range

This should not happen.  It likely indicates a bug
in the KIM implementation of the interatomic potential
where it is requesting neighbors incorrectly.

E: LAMMPS unit_style lj not supported by KIM models

Self-explanatory. Check the input script or data file.

E: Unknown unit_style

Self-explanatory. Check the input script or data file.

W: KIM Model does not provide `energy'; Potential energy will be zero

Self-explanatory.

W: KIM Model does not provide `forces'; Forces will be zero

Self-explanatory.

W: KIM Model does not provide `particleEnergy'; energy per atom will be zero

Self-explanatory.

W: KIM Model does not provide `particleVirial'; virial per atom will be zero



E: Test_descriptor_string already allocated

This is an internal error.  Contact the developers.

U: KIM Model does not provide 'energy'; Potential energy will be zero

Self-explanatory.

U: KIM Model does not provide 'forces'; Forces will be zero

Self-explanatory.

U: KIM Model does not provide 'particleEnergy'; energy per atom will be zero

Self-explanatory.

U: KIM Model does not provide 'particleVirial'; virial per atom will be zero

Self-explanatory.

*/
