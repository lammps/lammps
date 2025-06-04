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
PairStyle(metatomic, PairMetatomic);

// additional versions of pair_style metatomic, to be used with pair_style
// hybrid when combining multiple metatomic potentials together
PairStyle(metatomic_1, PairMetatomic);
PairStyle(metatomic_2, PairMetatomic);
PairStyle(metatomic_3, PairMetatomic);
PairStyle(metatomic_4, PairMetatomic);
PairStyle(metatomic_5, PairMetatomic);
PairStyle(metatomic_6, PairMetatomic);
PairStyle(metatomic_7, PairMetatomic);
PairStyle(metatomic_8, PairMetatomic);
PairStyle(metatomic_9, PairMetatomic);
// clang-format on
#else

#ifndef LMP_PAIR_METATOMIC_H
#define LMP_PAIR_METATOMIC_H

#include "pair.h"

#include <vector>

// this is the actual namespace where `torch::Device` is defined
namespace c10 {
    class Device;

    enum class DeviceType: int8_t;
}

namespace LAMMPS_NS {
class MetatomicSystemAdaptor;
struct PairMetatomicData;

class PairMetatomic : public Pair {
public:
    PairMetatomic(class LAMMPS *);
    ~PairMetatomic();

    void compute(int, int) override;
    void settings(int, char **) override;
    void coeff(int, char **) override;
    void init_style() override;
    double init_one(int, int) override;
    void init_list(int id, NeighList *ptr) override;

    void allocate();

protected:
    // get the set of devices both available on the current machine and supported
    // by the model
    std::vector<c10::DeviceType> available_devices();

    // pick the correct device to use from the user request (or nullptr) in
    // `pair_style metatomic`
    virtual void pick_device(c10::Device* device, const char* requested);

    PairMetatomicData* mta_data;
    NeighList *mta_list;

    // mapping from LAMMPS types to metatomic types
    int32_t *type_mapping;
    // adaptor from LAMMPS system to metatomic's
    std::unique_ptr<MetatomicSystemAdaptor> system_adaptor;

    double scale;
};

}    // namespace LAMMPS_NS

#endif
#endif
