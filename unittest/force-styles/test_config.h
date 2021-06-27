/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

#include <set>
#include <string>
#include <utility>
#include <vector>

struct coord_t {
    double x, y, z;
};

struct stress_t {
    double xx, yy, zz, xy, xz, yz;
};

class TestConfig {
public:
    std::string lammps_version;
    std::string date_generated;
    std::string basename;
    double epsilon;
    std::set<std::string> skip_tests;
    std::vector<std::pair<std::string, std::string>> prerequisites;
    std::vector<std::string> pre_commands;
    std::vector<std::string> post_commands;
    std::string input_file;
    std::string pair_style;
    std::string bond_style;
    std::string angle_style;
    std::string dihedral_style;
    std::string improper_style;
    std::string kspace_style;
    std::vector<std::string> pair_coeff;
    std::vector<std::string> bond_coeff;
    std::vector<std::string> angle_coeff;
    std::vector<std::string> dihedral_coeff;
    std::vector<std::string> improper_coeff;
    std::vector<double> equilibrium;
    std::vector<std::pair<std::string, int>> extract;
    int natoms;
    double init_energy;
    double run_energy;
    double init_vdwl;
    double run_vdwl;
    double init_coul;
    double run_coul;
    stress_t init_stress;
    stress_t run_stress;
    double global_scalar;
    std::vector<double> global_vector;
    std::vector<coord_t> init_forces;
    std::vector<coord_t> run_forces;
    std::vector<coord_t> run_pos;
    std::vector<coord_t> restart_pos;
    std::vector<coord_t> run_vel;
    std::vector<coord_t> restart_vel;

    TestConfig() :
        lammps_version(""), date_generated(""), basename(""), epsilon(1.0e-14), input_file(""),
        pair_style("zero"), bond_style("zero"), angle_style("zero"), dihedral_style("zero"),
        improper_style("zero"), kspace_style("none"), natoms(0), init_energy(0), run_energy(0),
        init_vdwl(0), run_vdwl(0), init_coul(0), run_coul(0), init_stress({0, 0, 0, 0, 0, 0}),
        run_stress({0, 0, 0, 0, 0, 0}), global_scalar(0)
    {
        skip_tests.clear();
        prerequisites.clear();
        pre_commands.clear();
        post_commands.clear();
        pair_coeff.clear();
        bond_coeff.clear();
        angle_coeff.clear();
        dihedral_coeff.clear();
        improper_coeff.clear();
        extract.clear();
        init_forces.clear();
        run_forces.clear();
        run_pos.clear();
        restart_pos.clear();
        run_vel.clear();
        restart_vel.clear();
        global_vector.clear();
    }
    virtual ~TestConfig(){};

private:
    TestConfig(const TestConfig &){};
};

#endif
