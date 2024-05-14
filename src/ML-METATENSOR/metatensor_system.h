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

#ifndef LMP_METATENSOR_SYSTEM_H
#define LMP_METATENSOR_SYSTEM_H

#include <vector>
#include <array>
#include <unordered_set>

#include "pointers.h"
#include "pair.h"
#include "neigh_list.h"

#include <metatensor/torch/atomistic.hpp>


namespace LAMMPS_NS {

struct MetatensorSystemOptions {
    // Mapping from LAMMPS types to metatensor types
    const int32_t* types_mapping;
    // interaction range of the model, in LAMMPS units
    double interaction_range;
    // should we run extra checks on the neighbor lists?
    bool check_consistency;
};

// data for metatensor neighbors lists
struct MetatensorNeighborsData {
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

    // cutoff for this NL in LAMMPS units
    double cutoff;
    // options of the NL as requested by the model
    metatensor_torch::NeighborListOptions options;

    // Below are cached allocations for the LAMMPS -> metatensor NL translation
    // TODO: report memory usage for these?

    // we keep the set of samples twice: once in `known_samples` to remove
    // duplicated pairs, and once in `samples` in a format that can be
    // used to create a torch::Tensor.
    std::unordered_set<sample_t, SampleHasher> known_samples;
    std::vector<sample_t> samples;
    // pairs distances vectors
    std::vector<std::array<double, 3>> distances_f64;
    std::vector<std::array<float, 3>> distances_f32;
};

class MetatensorSystemAdaptor : public Pointers {
public:
    MetatensorSystemAdaptor(LAMMPS* lmp, Pair* requestor, MetatensorSystemOptions options);
    MetatensorSystemAdaptor(LAMMPS* lmp, Compute* requestor, MetatensorSystemOptions options);

    ~MetatensorSystemAdaptor();

    void init_list(int id, NeighList* ptr);

    void add_nl_request(double cutoff, metatensor_torch::NeighborListOptions request);

    // Create a metatensor system matching the LAMMPS system data
    metatensor_torch::System system_from_lmp(bool do_virial, torch::ScalarType dtype, torch::Device device);

    // Explicit strain for virial calculations. This uses the same dtype/device
    // as LAMMPS data (positions, …)
    torch::Tensor strain;
    // keep the positions as coming from LAMMPS (before any dtype/device
    // conversion) to access its gradient
    torch::Tensor positions;

private:
    // setup the metatensor neighbors list from the internal LAMMPS one
    void setup_neighbors(metatensor_torch::System& system);

    // options for this system adaptor
    MetatensorSystemOptions options_;

    // LAMMPS NL
    NeighList* list_;
    // allocations caches for all the NL requested by
    // the model
    std::vector<MetatensorNeighborsData> caches_;
    // allocation cache for the atomic types in the system
    torch::Tensor atomic_types_;
    // allocation cache holding the "original atom" id for all atoms in the
    // system. This is the same as the atom id for all local atoms. For ghost
    // atoms, this is either the id of the corresponding local atom if the ghost
    // is a periodic image of a local atom, the id of the first ghost we found
    // with a given atom tag if the ghost is a periodic image of another ghost;
    // or the id of the ghost in all other cases.
    std::vector<int> original_atom_id_;
    // allocation cache holding the map from atom tag to atom id for local
    // atoms.
    std::unordered_map<tagint, int> local_atoms_tags_;
    // allocation cache holding the map from atom tag to atom id for ghost
    // atoms. When there are multiple periodic images of the same atom, only one
    // will be included here.
    std::unordered_map<tagint, int> ghost_atoms_tags_;

    // TODO: should we use LAMMPS allocations/deallocation facilities for the
    // allocation caches? If we don't, should we report memory usage from the
    // allocations caches to LAMMPS one way or another?
};

}    // namespace LAMMPS_NS

#endif
