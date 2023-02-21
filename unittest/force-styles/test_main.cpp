/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "test_main.h"
#include "atom.h"
#include "error_stats.h"
#include "lammps.h"
#include "pointers.h"
#include "test_config.h"
#include "test_config_reader.h"
#include "utils.h"
#include "yaml_writer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <mpi.h>
#include <vector>

using LAMMPS_NS::Atom;
using LAMMPS_NS::LAMMPS;
using LAMMPS_NS::tagint;
using LAMMPS_NS::utils::split_words;
using LAMMPS_NS::utils::trim;

void EXPECT_STRESS(const std::string &name, double *stress, const stress_t &expected_stress,
                   double epsilon)
{
    SCOPED_TRACE("EXPECT_STRESS: " + name);
    ErrorStats stats;
    EXPECT_FP_LE_WITH_EPS(stress[0], expected_stress.xx, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[1], expected_stress.yy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[2], expected_stress.zz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[3], expected_stress.xy, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[4], expected_stress.xz, epsilon);
    EXPECT_FP_LE_WITH_EPS(stress[5], expected_stress.yz, epsilon);
    if (print_stats) std::cerr << name << " stats" << stats << std::endl;
}

void EXPECT_FORCES(const std::string &name, Atom *atom, const std::vector<coord_t> &f_ref,
                   double epsilon)
{
    SCOPED_TRACE("EXPECT_FORCES: " + name);
    double **f       = atom->f;
    tagint *tag      = atom->tag;
    const int nlocal = atom->nlocal;
    ASSERT_EQ(nlocal + 1, f_ref.size());
    ErrorStats stats;
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(f[i][0], f_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][1], f_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(f[i][2], f_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << name << " stats" << stats << std::endl;
}

void EXPECT_POSITIONS(const std::string &name, Atom *atom, const std::vector<coord_t> &x_ref,
                      double epsilon)
{
    SCOPED_TRACE("EXPECT_POSITIONS: " + name);
    double **x       = atom->x;
    tagint *tag      = atom->tag;
    const int nlocal = atom->nlocal;
    ASSERT_EQ(nlocal + 1, x_ref.size());
    ErrorStats stats;
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(x[i][0], x_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][1], x_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(x[i][2], x_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << name << " stats" << stats << std::endl;
}

void EXPECT_VELOCITIES(const std::string &name, Atom *atom, const std::vector<coord_t> &v_ref,
                       double epsilon)
{
    SCOPED_TRACE("EXPECT_VELOCITIES: " + name);
    double **v       = atom->v;
    tagint *tag      = atom->tag;
    const int nlocal = atom->nlocal;
    ASSERT_EQ(nlocal + 1, v_ref.size());
    ErrorStats stats;
    for (int i = 0; i < nlocal; ++i) {
        EXPECT_FP_LE_WITH_EPS(v[i][0], v_ref[tag[i]].x, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][1], v_ref[tag[i]].y, epsilon);
        EXPECT_FP_LE_WITH_EPS(v[i][2], v_ref[tag[i]].z, epsilon);
    }
    if (print_stats) std::cerr << name << " stats" << stats << std::endl;
}

// common read_yaml_file function
bool read_yaml_file(const char *infile, TestConfig &config)
{
    auto reader = TestConfigReader(config);
    if (reader.parse_file(infile)) return false;

    config.basename = reader.get_basename();
    return true;
}

// write out common header items for yaml files
void write_yaml_header(YamlWriter *writer, TestConfig *cfg, const char *version)
{
    // lammps_version
    writer->emit("lammps_version", version);

    // tags
    writer->emit("tags", cfg->tags_line());

    // date_generated
    std::time_t now   = time(nullptr);
    std::string block = trim(ctime(&now));
    writer->emit("date_generated", block);

    // epsilon
    writer->emit("epsilon", cfg->epsilon);

    // skip tests
    block.clear();
    for (auto &skip : cfg->skip_tests) {
        if (block.empty())
            block = skip;
        else
            block += " " + skip;
    }
    writer->emit("skip_tests", block);

    // prerequisites
    block.clear();
    for (auto &prerequisite : cfg->prerequisites) {
        block += prerequisite.first + " " + prerequisite.second + "\n";
    }
    writer->emit_block("prerequisites", block);

    // pre_commands
    block.clear();
    for (auto &command : cfg->pre_commands) {
        block += command + "\n";
    }
    writer->emit_block("pre_commands", block);

    // post_commands
    block.clear();
    for (auto &command : cfg->post_commands) {
        block += command + "\n";
    }
    writer->emit_block("post_commands", block);

    // input_file
    writer->emit("input_file", cfg->input_file);
}

// need to be defined in unit test body
extern void generate_yaml_file(const char *, const TestConfig &);

void usage(std::ostream &out, const char *name)
{
    out << "usage: " << name << " <testfile.yaml> [OPTIONS]\n\n"
        << "Available options:\n"
        << "  -g <newfile.yaml>   regenerate yaml file under a new name\n"
        << "  -d <folder>         set folder where to find input files\n"
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

// location for 'in.*' and 'data.*' files required by tests
#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val
std::string INPUT_FOLDER = STRINGIFY(TEST_INPUT_FOLDER);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 2) {
        usage(std::cerr, argv[0]);
        return 1;
    }

    if (!read_yaml_file(argv[1], test_config)) {
        std::cerr << "Error parsing yaml file: " << argv[1] << std::endl;
        return 2;
    }

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-u") {
                generate_yaml_file(argv[1], test_config);
                return 0;
            } else if (arg == "-s") {
                print_stats = true;
            } else if (arg == "-v") {
                verbose = true;
            }
        }
    }

    int iarg = 2;
    while (iarg < argc) {

        if (strcmp(argv[iarg], "-g") == 0) {
            if (iarg + 1 < argc) {
                generate_yaml_file(argv[iarg + 1], test_config);
                MPI_Finalize();
                return 0;
            } else {
                usage(std::cerr, argv[0]);
                MPI_Finalize();
                return 1;
            }
        } else if (strcmp(argv[iarg], "-u") == 0) {
            generate_yaml_file(argv[1], test_config);
            MPI_Finalize();
            return 0;
        } else if (strcmp(argv[iarg], "-d") == 0) {
            if (iarg + 1 < argc) {
                INPUT_FOLDER = argv[iarg + 1];
                iarg += 2;
            } else {
                usage(std::cerr, argv[0]);
                MPI_Finalize();
                return 1;
            }
        } else if (strcmp(argv[iarg], "-s") == 0) {
            print_stats = true;
            ++iarg;
        } else if (strcmp(argv[iarg], "-v") == 0) {
            verbose = true;
            ++iarg;
        } else {
            std::cerr << "unknown option: " << argv[iarg] << "\n\n";
            usage(std::cerr, argv[0]);
            MPI_Finalize();
            return 1;
        }
    }

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
