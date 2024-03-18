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

#include <vector>
#include <array>
#include <unordered_set>

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
    void pairs_with_mapping(metatensor_torch::System& system);

    // == data for the model

    // the model itself
    std::unique_ptr<torch::jit::Module> torch_model;
    // model capabilities
    metatensor_torch::ModelCapabilities capabilities;
    // interaction range of the model, in LAMMPS units
    double interaction_range;
    // run-time evaluation options, set by us
    metatensor_torch::ModelEvaluationOptions evaluation_options;
    // should metatensor check the data LAMMPS send to the model
    // and the data the model returns?
    bool check_consistency;

    // == data for neighbors lists
    struct NeighborsData {
        // single neighbors sample containing [i, j, S_a, S_b, S_c]
        using sample_t = std::array<int32_t, 5>;

        struct SampleHasher {
            static void hash_combine(std::size_t& seed, const int32_t& v) {
                seed ^= std::hash<int32_t>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }

            size_t operator()(const sample_t& s) const {
                size_t hash = 0;
                hash_combine(hash, s[0]);
                hash_combine(hash, s[1]);
                hash_combine(hash, s[2]);
                hash_combine(hash, s[3]);
                hash_combine(hash, s[4]);
                return hash;
            }
        };

        double cutoff;
        metatensor_torch::NeighborsListOptions options;
        // we keep the set of samples twice: once in `known_samples` to remove
        // duplicated pairs, and once in `samples` in a format that can be
        // used to create a torch::Tensor.
        std::unordered_set<sample_t, SampleHasher> known_samples;
        std::vector<sample_t> samples;
        // pairs distances vectors
        std::vector<std::array<double, 3>> distances;
    };

    // cached allocations for the LAMMPS -> metatensor NL translation
    // TODO: report memory usage for these?
    std::vector<NeighborsData> neigh_cache;

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
