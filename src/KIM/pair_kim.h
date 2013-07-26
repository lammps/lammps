/* ----------------------------------------------------------------------
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
   Designed for use with the openkim-api-v1.2.0 (and newer) package
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(kim,PairKIM)

#else

#ifndef LMP_PAIR_KIM_H
#define LMP_PAIR_KIM_H

// includes from KIM & LAMMPS
class KIM_API_model;
#include "pair.h"


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
      virtual double init_one(int, int);
      virtual int pack_reverse_comm(int, int, double*);
      virtual void unpack_reverse_comm(int, int*, double*);
      virtual double memory_usage();

   private:
      // (nearly) all bool flags are not initialized in constructor, but set
      // explicitly in the indicated function.  All other data members are
      // initialized in constructor

      // values set in settings()
      char* kim_modelname;

      // values set in coeff()

      // values set in allocate(), called by coeff()
      void allocate();
      int* lmps_map_types_to_unique;

      // values set in coeff(), after calling allocate()
      char** lmps_unique_elements;  // names of unique elements given
                                    // in pair_coeff command
      int lmps_num_unique_elements;

      // values set in set_lmps_flags(), called from init_style()
      bool lmps_using_newton;
      bool lmps_using_molecular;
      bool lmps_hybrid;             // true if running with pair hybrid
      bool lmps_support_cluster;    // true if running in mode compat.
                                    // with CLUSTER
      enum unit_sys {REAL, METAL, SI, CGS, ELECTRON};
      unit_sys lmps_units;

      // values set in set_kim_model_has_flags(), called by kim_init()
      KIM_API_model* pkim;
      bool kim_model_has_energy;
      bool kim_model_has_forces;
      bool kim_model_has_particleEnergy;
      bool kim_model_has_particleVirial;

      // values set in kim_init(), after call to string_init(_)
      bool kim_init_ok;
      bool kim_model_using_half;
      bool kim_model_using_cluster;
      bool kim_model_using_Rij;
      int kim_ind_coordinates;
      int kim_ind_numberOfParticles;
      int kim_ind_numberContributingParticles;
      int kim_ind_numberParticleTypes;
      int kim_ind_particleTypes;
      int kim_ind_get_neigh;
      int kim_ind_neighObject;
      int kim_ind_cutoff;
      int kim_ind_energy;
      int kim_ind_particleEnergy;
      int kim_ind_forces;
      int kim_ind_virial;
      int kim_ind_particleVirial;

      // values set in init_style(), after calling pkim->model_init()
      bool kim_model_init_ok;
      bool kim_particle_codes_ok;
      int *kim_particle_codes;

      // values set in set_statics(), called at end of kim_init(),
      //   then again in set_volatiles(), called in compute()
      int lmps_local_tot_num_atoms;
      double kim_global_cutoff;     // KIM Model cutoff value

      // values set in compute()
      int lmps_maxalloc;            // max allocated memory value
      int* kim_particleTypes;       // array of KIM particle types
      double** lmps_force_tmp;      // temp storage for f, when running in
                                    // hybrid mode needed to avoid reseting
                                    // f to zero in each object
      int* lmps_stripped_neigh_list;// neighbors of one atom, used when LAMMPS
                                    // is in molecular mode

      // values used in get_neigh()
      int kim_iterator_position;    //get_neigh iterator current position
      double *Rij;

      // KIM specific helper functions
      void kim_error(int, const char *, int);
      void kim_init();
      void kim_free();
      void set_statics();
      void set_volatiles();
      void set_lmps_flags();
      void set_kim_model_has_flags();
      void write_descriptor(char** test_descriptor_string);
      // static methods used as callbacks from KIM
      static int get_neigh(void** kimmdl, int* mode, int* request,
                           int* atom, int* numnei, int** nei1atom,
                           double** pRij);
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

The KIM API does not explicitly support anything other than 3D problems

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Internal KIM error

Self-explanatory. Check the output and kim.log file for more details.

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

Self-explanatory.

E: test_descriptor_string already allocated

This should not happen. It likely indicates a bug in the pair_kim
implementation.

*/
