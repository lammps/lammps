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
#ifdef PAIR_CLASS
// clang-format off
PairStyle(metatensor, PairMetatensor);
// clang-format on
#else

#ifndef LMP_PAIR_METATENSOR_H
#define LMP_PAIR_METATENSOR_H

#include <memory>

#include "pair.h"

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>

namespace LAMMPS_NS {

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
    void load_torch_model(const char* path);
    metatensor_torch::System system_from_lmp(bool do_virial);

    // cached allocations for the neighbors list TensorBlock
    struct NeighborsData {
        double cutoff = -1;
        std::vector<int32_t> samples;
        std::vector<double> distances;
    };

    std::vector<metatensor_torch::NeighborsListOptions> neigh_options;
    std::vector<NeighList*> neigh_lists;
    std::vector<NeighborsData> neigh_data_cache;

    // TODO: make this user-configurable
    bool check_consistency = true;

    std::unique_ptr<torch::jit::Module> torch_model;
    metatensor_torch::ModelCapabilities capabilities;
    metatensor_torch::ModelEvaluationOptions evaluation_options;
    // various allocation caches
    torch::Tensor selected_atoms_values;
    torch::Tensor atomic_types;

    // explicit strain for virial calculations
    torch::Tensor strain;

    // mapping from LAMMPS types to metatensor types
    int32_t* lammps_to_metatensor_types = nullptr;
};

}    // namespace LAMMPS_NS

#endif
#endif
