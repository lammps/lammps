/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Guillaume Fraux <guillaume.fraux@epfl.ch>
------------------------------------------------------------------------- */
#include "metatomic_types.h"

#include "citeme.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

PairMetatomicData::PairMetatomicData(std::string length_unit, std::string energy_unit):
    device(torch::kCPU),
    check_consistency(false),
    remap_pairs(true),
    non_conservative(false),
    max_cutoff(-1)
{
    auto options = torch::TensorOptions().dtype(torch::kInt32);
    this->selected_atoms_values = torch::zeros({0, 2}, options);

    // Initialize evaluation_options
    this->evaluation_options = torch::make_intrusive<metatomic_torch::ModelEvaluationOptionsHolder>();
    this->evaluation_options->set_length_unit(std::move(length_unit));

    auto output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
    output->explicit_gradients = {};
    output->set_quantity("energy");
    output->set_unit(std::move(energy_unit));
    output->per_atom = false;

    this->evaluation_options->outputs.insert("energy", output);
}

void PairMetatomicData::load_model(
   LAMMPS* lmp,
   const char* path,
   const char* extensions_directory
) {
   // TODO: seach for the model & extensions inside `$LAMMPS_POTENTIALS`?

   if (this->model != nullptr) {
       lmp->error->all(FLERR, "torch model is already loaded");
   }

   torch::optional<std::string> extensions = torch::nullopt;
   if (extensions_directory != nullptr) {
       extensions = std::string(extensions_directory);
   }

   try {
       this->model = std::make_unique<metatensor_torch::Module>(
           metatomic_torch::load_atomistic_model(path, extensions)
       );
   } catch (const c10::Error& e) {
       lmp->error->all(FLERR, "failed to load metatomic model at '{}': {}", path, e.what());
   }

   auto capabilities_ivalue = this->model->run_method("capabilities");
   this->capabilities = capabilities_ivalue.toCustomClass<metatomic_torch::ModelCapabilitiesHolder>();

   if (!this->capabilities->outputs().contains("energy")) {
       lmp->error->all(FLERR, "the model at '{}' does not have an \"energy\" output, we can not use it in pair_style metatomic", path);
   }

   if (lmp->comm->me == 0) {
       auto metadata_ivalue = this->model->run_method("metadata");
       auto metadata = metadata_ivalue.toCustomClass<metatomic_torch::ModelMetadataHolder>();
       auto to_print = metadata->print();

       if (lmp->screen) {
           fprintf(lmp->screen, "\n%s\n", to_print.c_str());
       }
       if (lmp->logfile) {
           fprintf(lmp->logfile,"\n%s\n", to_print.c_str());
       }

       // add the model references to LAMMPS citation handling mechanism
       if (lmp->citeme) {
          for (const auto& it: metadata->references) {
             for (const auto& ref: it.value()) {
                lmp->citeme->add(ref + "\n");
             }
          }
       }
   }
}
