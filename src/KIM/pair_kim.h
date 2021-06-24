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
   Contributing authors: Ryan S. Elliott (UMinn), Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.0.2 (and newer) package
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(kim,PairKIM);
// clang-format on
#else

#ifndef LMP_PAIR_KIM_H
#define LMP_PAIR_KIM_H

// includes from KIM & LAMMPS
class KIM_API_model;
#include "pair.h"

extern "C" {
#include "KIM_SimulatorHeaders.h"
}

namespace LAMMPS_NS {

class PairKIM : public Pair {
 public:
  PairKIM(class LAMMPS *);
  ~PairKIM();

  // LAMMPS Pair class virtual function prototypes
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual void init_list(int id, NeighList *ptr);
  virtual double init_one(int, int);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);
  virtual double memory_usage();

  // Get the KIM_Model object
  KIM_Model *get_kim_model();
  // Get the atom type list
  std::string get_atom_type_list();

 protected:
  // (nearly) all bool flags are not initialized in constructor, but set
  // explicitly in the indicated function.  All other data members are
  // initialized in constructor
  int settings_call_count;
  int init_style_call_count;

  // values set in settings()
  char *kim_modelname;

  // list of args that map atom species to KIM elements
  std::string atom_type_list;

  // values set in coeff()

  // values set in allocate(), called by coeff()
  virtual void allocate();
  int *lmps_map_species_to_unique;

  // values set in coeff(), after calling allocate()
  char **lmps_unique_elements;    // names of unique elements given
                                  // in pair_coeff command
  int lmps_num_unique_elements;

  // values set in set_lmps_flags(), called from init_style()
  bool lmps_using_newton;
  bool lmps_using_molecular;
  enum unit_sys { REAL, METAL, SI, CGS, ELECTRON };
  unit_sys lmps_units;
  KIM_LengthUnit lengthUnit;
  KIM_EnergyUnit energyUnit;
  KIM_ChargeUnit chargeUnit;
  KIM_TemperatureUnit temperatureUnit;
  KIM_TimeUnit timeUnit;

  KIM_Model *pkim;
  KIM_ComputeArguments *pargs;

  // values set in set_kim_model_has_flags(), called by kim_init()
  KIM_SupportStatus kim_model_support_for_energy;
  KIM_SupportStatus kim_model_support_for_forces;
  KIM_SupportStatus kim_model_support_for_particleEnergy;
  KIM_SupportStatus kim_model_support_for_particleVirial;

  // values set in kim_init()
  bool kim_init_ok;
  int lmps_local_tot_num_atoms;
  double kim_global_influence_distance;    // KIM Model cutoff value
  int kim_number_of_neighbor_lists;
  double const *kim_cutoff_values;
  int const *modelWillNotRequestNeighborsOfNoncontributingParticles;
  class NeighList **neighborLists;

  // values set in init_style()
  bool kim_particle_codes_ok;
  int *kim_particle_codes;

  // values set in compute()
  int lmps_maxalloc;                // max allocated memory value
  int *kim_particleSpecies;         // array of KIM particle species
  int *kim_particleContributing;    // array of KIM particle contributing
  int *lmps_stripped_neigh_list;    // neighbors of one atom, used when LAMMPS
                                    // is in molecular mode
  int **lmps_stripped_neigh_ptr;    // pointer into lists

  // KIM specific helper functions
  virtual void set_contributing();
  virtual void kim_init();
  virtual void kim_free();
  virtual void set_argument_pointers();
  virtual void set_lmps_flags();
  virtual void set_kim_model_has_flags();
  virtual int check_for_routine_compatibility();
  // static methods used as callbacks from KIM
  static int get_neigh(void const *const dataObject, int const numberOfCutoffs,
                       double const *const cutoffs, int const neighborListIndex,
                       int const particleNumber, int *const numberOfNeighbors,
                       int const **const neighborsOfParticle);
};
}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Unable to set KIM particle species codes and/or contributing

A low-level kim-api error has occurred.

E: KIM Compute returned error

The KIM model was unable, for some reason, to complete the computation.

E: 'KIMvirial' or 'LAMMPSvirial' not supported with kim-api.

"KIMvirial or "LAMMPSvirial" found on the pair_style line.  These keys
are not supported kim-api.  (The virial computation is always performed
by LAMMPS.) Please remove these keys, make sure the KIM model you are
using supports kim-api, and rerun.

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
