/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Norbert Podhorszki (Oak Ridge National Laboratory)
------------------------------------------------------------------------- */

#ifndef LMP_ADIOS_COMMON_H
#define LMP_ADIOS_COMMON_H

#include <string_view>

// Default runtime configuration for all ADIOS2 dump classes.
//
// Engine upgrade: BP4 → BP5 (default since ADIOS2 2.9; lower memory
// footprint, faster metadata, better multi-step append performance).
//
// Parameter rename: "substreams" (BP4) → "NumAggregators" (BP5).
// NumAggregators defaults to one process per compute node when unset;
// the value "1" here is a conservative baseline that users can override
// by supplying their own adios2_config.xml before the run starts.
//
// Using inline constexpr std::string_view avoids the multiple-definition
// problem of a plain const char[] in a header included by several TUs.

inline constexpr std::string_view default_config =
    R"xml(<?xml version="1.0"?>
<adios-config>
    <io name="atom">
        <engine type="BP5">
            <parameter key="NumAggregators" value="1"/>
        </engine>
    </io>
    <io name="custom">
        <engine type="BP5">
            <parameter key="NumAggregators" value="1"/>
        </engine>
    </io>
    <io name="local">
        <engine type="BP5">
            <parameter key="NumAggregators" value="1"/>
        </engine>
    </io>
    <io name="read_dump">
        <engine type="BP5">
        </engine>
    </io>
</adios-config>
)xml";

#endif
