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

// common definitions for all ADIOS package classes

static const char default_config[] = "<?xml version=\"1.0\"?>\n"
                                     "<adios-config>\n"
                                     "    <io name=\"atom\">\n"
                                     "        <engine type=\"BP4\">\n"
                                     "            <parameter key=\"substreams\" value=\"1\"/>\n"
                                     "        </engine>\n"
                                     "    </io>\n"
                                     "    <io name=\"custom\">\n"
                                     "        <engine type=\"BP4\">\n"
                                     "            <parameter key=\"substreams\" value=\"1\"/>\n"
                                     "        </engine>\n"
                                     "    </io>\n"
                                     "    <io name=\"read_dump\">\n"
                                     "        <engine type=\"BP4\">\n"
                                     "        </engine>\n"
                                     "    </io>\n"
                                     "</adios-config>\n";

#endif
