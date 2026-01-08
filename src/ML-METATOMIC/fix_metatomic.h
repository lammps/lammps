/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(metatomic, FixMetatomic);
// clang-format on
#else

#ifndef LMP_FIX_METATOMIC_H
#define LMP_FIX_METATOMIC_H

#include "fix.h"

#include <string>
#include <optional>
#include <memory>

namespace c10 {
    class Device;
    enum class DeviceType: int8_t;
}

namespace at {
    class Tensor;
}

namespace LAMMPS_NS {
class MetatomicSystemAdaptor;
class FixMetatomicData;

class FixMetatomic : public Fix {
public:
    FixMetatomic(class LAMMPS *, int, char **);
    ~FixMetatomic();

    int setmask() override;
    void init() override;

    // Integration methods for ML-driven dynamics
    void initial_integrate(int) override;  // ML prediction of positions/momenta
    void post_force(int) override;         // Snapshot forces for Langevin compatibility
    void final_integrate() override;       // Apply force corrections
    void init_list(int id, NeighList *ptr) override;

 protected:
    virtual void pick_device(c10::Device& device, const char* requested);

    double momentum_conversion_factor;    // Conversion factor for momenta
    double dt;                            // Timestep
    std::string model_path;               // Path to ML model file
    std::optional<std::string> extensions_directory;     // Directory for model extensions
    std::optional<std::string> requested_device;         // Device to run model on (cpu/cuda/mps)

    // Metatomic model data and configuration
    FixMetatomicData* mta_data;
    NeighList *mta_list;
    int mta_list_reqid;

    // Mapping from LAMMPS atom types to metatomic model types
    int32_t *type_mapping;

    // Helper class to convert between LAMMPS and metatomic representations
    std::unique_ptr<MetatomicSystemAdaptor> system_adaptor;
};

}    // namespace LAMMPS_NS

#endif
#endif
