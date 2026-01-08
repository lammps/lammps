/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"

#include <string>

#include <torch/torch.h>
#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>


#ifndef LMP_METATOMIC_TYPES_H
#define LMP_METATOMIC_TYPES_H

namespace LAMMPS_NS {

struct CommonMetatomicData {
   CommonMetatomicData(std::string length_unit);
   void load_model(LAMMPS* lmp, const char* path, const char* extensions_directory);

   // the metatomic model
   std::unique_ptr<metatensor_torch::Module> model;
   // the path used to load the model
   std::string model_path;
   // device to use for the calculations
   torch::Device device;
   // model capabilities, declared by the model
   metatomic_torch::ModelCapabilities capabilities;
   // run-time evaluation options, decided by this class
   metatomic_torch::ModelEvaluationOptions evaluation_options;

   // should metatomic check the data LAMMPS send to the model
   // and the data the model returns?
   bool check_consistency = false;
   // how far away the model needs to know about neighbors
   double max_cutoff = -1;

   // allocation cache for the selected atoms
   torch::Tensor selected_atoms_values;
};

struct PairMetatomicData: public CommonMetatomicData {
   PairMetatomicData(std::string length_unit): CommonMetatomicData(std::move(length_unit)) {}

   // the energy output we'll request from a model
   metatomic_torch::ModelOutput energy_output;
   // wether the model capabilities say that it can do per-atom energies
   bool is_energy_output_per_atom = false;

   // energy uncertainty output we'll request from a model, or nullptr if the
   // model does not have such output
   metatomic_torch::ModelOutput uncertainty_output;
   // threshold for energy uncertainty warnings
   double uncertainty_threshold = 0.0;

   // non-conservative forces/stress outputs we'll request from a model, or
   // nullptr if the model does not have such outputs
   metatomic_torch::ModelOutput nc_forces_output;
   metatomic_torch::ModelOutput nc_stress_output;

   // whether non-conservative forces and stresses should be used
   bool non_conservative = false;

   // energy key for the model
   std::string energy_key;
   // energy uncertainty key for the model
   std::string energy_uq_key;
   // non-conservative forces key for the model
   std::string nc_forces_key;
   // non-conservative stress key for the model
   std::string nc_stress_key;
};

struct FixMetatomicData: public CommonMetatomicData {
   FixMetatomicData(std::string length_unit): CommonMetatomicData(std::move(length_unit)) {}
};

}    // namespace LAMMPS_NS

#endif
