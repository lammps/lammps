/*
- there are 7 basic functions
  - delete_file
  - cleanup_lammps
  - init_lammps
  - run_lammps
  - restart_lammps
  - data_lammps
  - generate_yaml_file
- move delete_file amd cleanup_lammps to a single file
- I don't understand utility lambda, but they are being reused multiple times
- add as many comments as possible, to show my understanding and document the code
- code for matching forces, energy and stress are repeated 3 times
- run_lammps looks to be same across all tests - it isn't, there's subtle difference
*/

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// unit tests for improper styles intended for molecular systems

#include "error_stats.h"
#include "test_config.h"
#include "test_config_reader.h"
#include "test_main.h"
#include "yaml_reader.h"
#include "yaml_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "atom.h"
#include "improper.h"
#include "compute.h"
#include "fmt/format.h"
#include "force.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "universe.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mpi.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

using ::testing::HasSubstr;
using ::testing::StartsWith;

using namespace LAMMPS_NS;

static void delete_file(const std::string &filename)
{
    remove(filename.c_str());
};

void cleanup_lammps(LAMMPS *lmp, const TestConfig &cfg)
{
    delete_file(cfg.basename + ".restart");
    delete_file(cfg.basename + ".data");
    delete_file(cfg.basename + "-coeffs.in");
    delete lmp;
}

void run_lammps(LAMMPS *lmp)
{
    // utility lambda to improve readability
    auto command = [&](const std::string &line) {
        lmp->input->one(line.c_str());
    };

    command("fix 1 all nve");
    command("compute pe all pe/atom improper");
    command("compute sum all reduce sum c_pe");
    command("thermo_style custom step temp pe press c_sum");
    command("thermo 2");
    command("run 4 post no");
}

LAMMPS *init_lammps(int argc, char **argv, const TestConfig &cfg, const bool newton = true){}

void restart_lammps(LAMMPS *lmp, const TestConfig &cfg){}

void data_lammps(LAMMPS *lmp, const TestConfig &cfg){}

void generate_yaml_file(const char *outfile, const TestConfig &config){}

TEST(ImproperStyle, plain){}
TEST(ImproperStyle, omp){}
TEST(ImproperStyle, single){}
