/* -*- c++ -*- ----------------------------------------------------------
 * LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 * https://www.lammps.org/, Sandia National Laboratories
 * LAMMPS development team: developers@lammps.org
 *
 * Copyright (2003) Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.  This software is distributed under
 * the GNU General Public License.
 *
 * See the README file in the top-level LAMMPS directory.
 * ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Vishnu Raghuraman (UIUC)
   Email: vishnura@illinois.edu
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rmc/partial, FixRMCPartial)
// clang-format on
#else

#ifndef LMP_FIX_RMC_PARTIAL_H
#define LMP_FIX_RMC_PARTIAL_H

#include "fix.h"
#include <random>


namespace LAMMPS_NS {

class FixRMCPartial: public Fix {
   public:
   FixRMCPartial(class LAMMPS *, int, char **);
   ~FixRMCPartial();
   int setmask() override;
   int determine_molecule(int);
   int determine_dopant_or_semiconductor(int);
   double energy_full();
   double change_dihedral_parameters(int, int);
   double change_angle_parameters(int, int);
   double change_bond_parameters(int, int);
   void initial_integrate(int) override;
   void make_move();
   void post_mortem();

  struct Mol {
    Mol(int num_atoms_max);
    ~Mol();

    double **pos;
    double *charge;
    double *new_charge;
    double *mass;
    int *type;
    imageint *image;
    int *global_tag;
    int *local_tag;
    int *local_index;
    int local_atoms;
    int max_atoms;
    int charge_state;
  };
   int determine_charge_state(struct Mol*, double);
   void calculateMoleculeCOM(double*, struct Mol*);
   double calculateColoumbSelf(double **, int);
   void calculateSpecialBondMatrix(double **, struct Mol*, int);
   double calDistance(double*, double*);
   void calculateColoumbCross(double **, int, double **, int, double*);
   void bringMoleculeTogether(struct Mol*, double **, int);
   double getSpecialBondCoefficient(int, int);
   void modify_charge(struct Mol*, double*);
   void restore_charge(struct Mol*);

   private:
   int perform_step;
   int periodicity;
   int nmoves;
   int dopant_size;
   int semiconductor_size;
   int size_limit;
   int n_dopant;
   int n_semiconductor;
   int n_molecules;
   int type_threshold;
   int num_charge_states;
   int num_dihedrals;
   int num_angles;
   int num_bonds;
   int restart;
   int do_dihedral;
   int do_angle;
   int do_bond;
   std::string sysname;
   
   int acceptances;
   int rejections;
   int nmcsteps;
   int *num_dopant_charge;
   int *num_semiconductor_charge;
   double *molecule_charge_states;
   double *molecule_type;
  Mol *get_molecule(int, int);

   double acceptance_rate;
   double temperature;
   double beta;
   double delta_g; 
   double deltaQ;

   // Initialize random number generators for atoms
   std::random_device device_atom, device_type_source, device_type_destination, device_acc;
   std::mt19937 rng_atom, rng_type_source, rng_type_destination, rng_acc;
   std::uniform_int_distribution<> atom_dist;
   std::uniform_int_distribution<> type_dist;
   std::uniform_real_distribution<> acc_dist;



   int *dihedral_types;
   int **dihedral_list;
   int *angle_types;
   int **angle_list;
   int *bond_types;
   int **bond_list;
   double *delta_g_list;
   double *charges;
   double *dde;
   double *doping_efficiency;
   double **semiconductor_charges;
   double **dopant_charges;
   double **osc_mol_c, **dopant_mol_c;
   class Compute *c_pe;
};
};
#endif
#endif
