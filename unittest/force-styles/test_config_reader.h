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

#ifndef TEST_CONFIG_READER_H
#define TEST_CONFIG_READER_H

#include "test_config.h"
#include "yaml_reader.h"

class TestConfigReader : public YamlReader<TestConfigReader> {
    TestConfig &config;

public:
    TestConfigReader(TestConfig &config);
    TestConfigReader() = delete;
    const TestConfigReader & operator=(TestConfig &) = delete;

    void skip_tests(const yaml_event_t &event);
    void prerequisites(const yaml_event_t &event);
    void pre_commands(const yaml_event_t &event);
    void post_commands(const yaml_event_t &event);
    void lammps_version(const yaml_event_t &event);
    void date_generated(const yaml_event_t &event);
    void epsilon(const yaml_event_t &event);
    void input_file(const yaml_event_t &event);
    void extract(const yaml_event_t &event);
    void natoms(const yaml_event_t &event);
    void init_stress(const yaml_event_t &event);
    void run_stress(const yaml_event_t &event);
    void init_forces(const yaml_event_t &event);
    void run_forces(const yaml_event_t &event);
    void run_pos(const yaml_event_t &event);
    void run_vel(const yaml_event_t &event);
    void pair_style(const yaml_event_t &event);
    void pair_coeff(const yaml_event_t &event);
    void bond_style(const yaml_event_t &event);
    void bond_coeff(const yaml_event_t &event);
    void angle_style(const yaml_event_t &event);
    void angle_coeff(const yaml_event_t &event);
    void dihedral_style(const yaml_event_t &event);
    void dihedral_coeff(const yaml_event_t &event);
    void improper_style(const yaml_event_t &event);
    void improper_coeff(const yaml_event_t &event);
    void equilibrium(const yaml_event_t &event);
    void init_vdwl(const yaml_event_t &event);
    void init_coul(const yaml_event_t &event);
    void run_vdwl(const yaml_event_t &event);
    void run_coul(const yaml_event_t &event);
    void init_energy(const yaml_event_t &event);
    void run_energy(const yaml_event_t &event);
    void global_scalar(const yaml_event_t &event);
    void global_vector(const yaml_event_t &event);
    void tags(const yaml_event_t &event);
};

#endif
