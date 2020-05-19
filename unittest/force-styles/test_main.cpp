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

#include "test_main.h"
#include "test_config.h"
#include "gtest/gtest.h"

#include <mpi.h>
#include <cstring>
#include <iostream>

// need to be defined in unit test body
extern bool read_yaml_file(const char *, TestConfig &);
extern void generate_yaml_file(const char *, const TestConfig &);

void usage(std::ostream &out, const char *name)
{
    out << "usage: " << name << " <testfile.yaml> [OPTIONS]\n\n"
        << "Available options:\n"
        << "  -g <newfile.yaml>   regenerate yaml file under a new name\n"
        << "  -u                  update the original yaml file\n"
        << "  -v                  run tests with verbose output\n"
        << "  -s                  run tests with error statistics output\n"
        << "  -h                  print this message\n"
        << std::endl;
}
    
// test configuration settings read from yaml file
TestConfig test_config;

// whether to print error statistics
bool print_stats = false;

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);

    if (argc < 2) {
        usage(std::cerr, argv[0]);
        return 1;
    }

    if (!read_yaml_file(argv[1], test_config)) {
        std::cerr << "Error parsing yaml file: " << argv[1] << std::endl;
        return 2;
    }

    int iarg=2;
    while (iarg < argc) {

        if (strcmp(argv[iarg],"-g") == 0) {
            if (iarg+1 < argc) {
                generate_yaml_file(argv[iarg+1], test_config);
                return 0;
            } else {
                usage(std::cerr,argv[0]);
                return 1;
            }
        } else if (strcmp(argv[iarg],"-u") == 0) {
            generate_yaml_file(argv[1], test_config);
            return 0;
        } else if (strcmp(argv[iarg],"-s") == 0) {
            print_stats = true;
            ++iarg;
        } else if (strcmp(argv[iarg],"-v") == 0) {
            verbose = true;
            ++iarg;
        } else {
            std::cerr << "unknown option: " << argv[iarg] << "\n\n";
            usage(std::cerr, argv[0]);
            return 1;
        }
    }
    return RUN_ALL_TESTS();
}
