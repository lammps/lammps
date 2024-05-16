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
#ifdef PAIR_CLASS
// clang-format off
PairStyle(metatensor, PairMetatensor);
// clang-format on
#else

#ifndef LMP_PAIR_METATENSOR_H
#define LMP_PAIR_METATENSOR_H

#include <memory>

#include "pair.h"

#include <metatensor/torch/atomistic.hpp>

namespace LAMMPS_NS {
class MetatensorSystemAdaptor;

class PairMetatensor : public Pair {
public:
    PairMetatensor(class LAMMPS *);
    ~PairMetatensor();

    void compute(int, int) override;
    void settings(int, char **) override;
    void coeff(int, char **) override;
    void init_style() override;
    double init_one(int, int) override;
    void init_list(int id, NeighList *ptr) override;

    void allocate();

private:
    void load_torch_model(const char* path, const char* extensions_directory);

    // torch model in metatensor format
    std::unique_ptr<torch::jit::Module> torch_model;
    // device to use for the calculations
    torch::Device device;
    // model capabilities, declared by the model
    metatensor_torch::ModelCapabilities capabilities;
    // run-time evaluation options, decided by this class
    metatensor_torch::ModelEvaluationOptions evaluation_options;
    // should metatensor check the data LAMMPS send to the model
    // and the data the model returns?
    bool check_consistency;
    // how far away the model needs to know about neighbors
    double interaction_range = -1;

    // allocation cache for the selected atoms
    torch::Tensor selected_atoms_values;
    // adaptor from LAMMPS system to metatensor's
    std::unique_ptr<MetatensorSystemAdaptor> system_adaptor;
    // mapping from LAMMPS types to metatensor types
    int32_t* type_mapping = nullptr;
};

}    // namespace LAMMPS_NS

#endif
#endif
