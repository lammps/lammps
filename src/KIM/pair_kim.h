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
   Contributing authors: Ryan S. Elliott (UMinn)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-v2.0.0-beta.1 (and newer) package
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
#include <sstream>


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

   protected:
      // (nearly) all bool flags are not initialized in constructor, but set
      // explicitly in the indicated function.  All other data members are
      // initialized in constructor
      int settings_call_count;
      int init_style_call_count;

      // values set in settings()
      char* kim_modelname;

      // values set in coeff()

      // values set in allocate(), called by coeff()
      virtual void allocate();
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
      KIM::SupportStatus kim_model_support_for_energy;
      KIM::SupportStatus kim_model_support_for_forces;
      KIM::SupportStatus kim_model_support_for_particleEnergy;
      KIM::SupportStatus kim_model_support_for_particleVirial;

      // values set in kim_init()
      bool kim_init_ok;
      int lmps_local_tot_num_atoms;
      double kim_global_influence_distance;  // KIM Model cutoff value
      int kim_number_of_neighbor_lists;
      double const * kim_cutoff_values;
      int const * padding_neighbor_hints;
      int const * half_list_hints;
      class NeighList ** neighborLists;

      // values set in init_style()
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
      virtual void set_contributing();
      virtual void kim_init();
      virtual void kim_free();
      virtual void set_argument_pointers();
      virtual void set_lmps_flags();
      virtual void set_kim_model_has_flags();
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

E: Unable to set KIM particle species codes and/or contributing

A low-level kim-api error has occurred.

E: KIM Compute returned error

The KIM model was unable, for some reason, to complete the computation.

E: Illegal pair_style command

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.

E: create_kim_particle_codes: symbol not found: XX

The KIM model specified does not support the atomic species symbol

E: PairKIM only works with 3D problems

Self-explanatory.

E: All pair coeffs are not set

Self-explanatory.

E: Unable to destroy Compute Arguments Object

A low-level kim-api error has occurred.

E: KIM ModelCreate failed

The kim-api was not able to create a model object for the specified model.

E: KIM Model did not accept the requested unit system

The KIM Model does not support the specified LAMMPS unit system

E: KIM ComputeArgumentsCreate failed

A low-level kim-api error has occurred.

E: Unable to register KIM pointers

A low-level kim-api error has occurred.

E: Unable to set KIM argument pointers

A low-level kim-api error has occurred.

E: pair_kim does not support hybrid

Self-explanatory.

E: LAMMPS unit_style lj not suppored by KIM models

Self-explanatory.

E: KIM Model requires unsupported compute argument: XXX

A low-level kim-api error has occurred.

W: KIM Model does not provide `partialEnergy'; Potential energy will be zero

Self-explanatory.

W: KIM Model does not provide `partialForce'; Forces will be zero

Self-explanatory.

W: KIM Model does not provide `partialParticleEnergy'; energy per atom will be zero

Self-explanatory.

W: KIM Model does not provide `partialParticleVirial'; virial per atom will be zero

Self-explanatory.

E: KIM Model requires unsupported compute callback

A low-level kim-api error has occurred.

*/
