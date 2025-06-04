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

struct PairMetatomicData {
   PairMetatomicData(std::string length_unit, std::string energy_unit);

   void load_model(LAMMPS* lmp, const char* path, const char* extensions_directory);

   // torch model in metatensor format
   std::unique_ptr<torch::jit::Module> model;
   // device to use for the calculations
   torch::Device device;
   // model capabilities, declared by the model
   metatomic_torch::ModelCapabilities capabilities;
   // run-time evaluation options, decided by this class
   metatomic_torch::ModelEvaluationOptions evaluation_options;
   // should metatomic check the data LAMMPS send to the model
   // and the data the model returns?
   bool check_consistency;
   // whether pairs should be remapped, removing pairs between ghosts if there
   // is an equivalent pair involving at least one local atom.
   bool remap_pairs;
   // whether non-conservative forces and stresses should be used
   bool non_conservative;
   // how far away the model needs to know about neighbors
   double max_cutoff;

   // allocation cache for the selected atoms
   torch::Tensor selected_atoms_values;
};

}    // namespace LAMMPS_NS

#endif
