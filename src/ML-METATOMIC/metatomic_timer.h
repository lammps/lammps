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

#ifndef LMP_METATOMIC_TIMER_H
#define LMP_METATOMIC_TIMER_H

#include <chrono>
#include <string>

namespace LAMMPS_NS {

/// Simple timer for profiling the LAMMPS/Metatensor integration. This starts
/// the timer when created, and print the elapsed time to stderr when going out
/// of scope.
class MetatomicTimer {
public:
    MetatomicTimer(std::string name);
    ~MetatomicTimer();

    // enable/disable profiling
    static void enable(bool toggle);

private:
    bool enabled_;
    std::string name_;
    size_t starting_counter_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

}    // namespace LAMMPS_NS

#endif
