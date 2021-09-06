/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "test_config_reader.h"
#include "test_config.h"
#include "utils.h"
#include "yaml.h"
#include "yaml_reader.h"

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using LAMMPS_NS::utils::split_words;

TestConfigReader::TestConfigReader(TestConfig &config) : YamlReader(), config(config)
{
    consumers["lammps_version"] = &TestConfigReader::lammps_version;
    consumers["date_generated"] = &TestConfigReader::date_generated;
    consumers["epsilon"]        = &TestConfigReader::epsilon;
    consumers["skip_tests"]     = &TestConfigReader::skip_tests;
    consumers["prerequisites"]  = &TestConfigReader::prerequisites;
    consumers["pre_commands"]   = &TestConfigReader::pre_commands;
    consumers["post_commands"]  = &TestConfigReader::post_commands;
    consumers["input_file"]     = &TestConfigReader::input_file;
    consumers["extract"]        = &TestConfigReader::extract;
    consumers["natoms"]         = &TestConfigReader::natoms;
    consumers["init_stress"]    = &TestConfigReader::init_stress;
    consumers["run_stress"]     = &TestConfigReader::run_stress;
    consumers["init_forces"]    = &TestConfigReader::init_forces;
    consumers["run_forces"]     = &TestConfigReader::run_forces;
    consumers["run_pos"]        = &TestConfigReader::run_pos;
    consumers["run_vel"]        = &TestConfigReader::run_vel;

    consumers["pair_style"] = &TestConfigReader::pair_style;
    consumers["pair_coeff"] = &TestConfigReader::pair_coeff;
    consumers["init_vdwl"]  = &TestConfigReader::init_vdwl;
    consumers["init_coul"]  = &TestConfigReader::init_coul;
    consumers["run_vdwl"]   = &TestConfigReader::run_vdwl;
    consumers["run_coul"]   = &TestConfigReader::run_coul;

    consumers["global_scalar"] = &TestConfigReader::global_scalar;
    consumers["global_vector"] = &TestConfigReader::global_vector;

    consumers["bond_style"]     = &TestConfigReader::bond_style;
    consumers["bond_coeff"]     = &TestConfigReader::bond_coeff;
    consumers["angle_style"]    = &TestConfigReader::angle_style;
    consumers["angle_coeff"]    = &TestConfigReader::angle_coeff;
    consumers["dihedral_style"] = &TestConfigReader::dihedral_style;
    consumers["dihedral_coeff"] = &TestConfigReader::dihedral_coeff;
    consumers["improper_style"] = &TestConfigReader::improper_style;
    consumers["improper_coeff"] = &TestConfigReader::improper_coeff;
    consumers["init_energy"]    = &TestConfigReader::init_energy;
    consumers["run_energy"]     = &TestConfigReader::run_energy;
    consumers["equilibrium"]    = &TestConfigReader::equilibrium;
}

void TestConfigReader::skip_tests(const yaml_event_t &event)
{
    config.skip_tests.clear();
    for (auto &word : split_words((char *)event.data.scalar.value))
        config.skip_tests.insert(word);
}

void TestConfigReader::prerequisites(const yaml_event_t &event)
{
    config.prerequisites.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string key, value;

    while (1) {
        data >> key >> value;
        if (data.eof()) break;
        config.prerequisites.push_back(std::make_pair(key, value));
    }
}

void TestConfigReader::pre_commands(const yaml_event_t &event)
{
    config.pre_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.pre_commands.push_back(line);
    }
}

void TestConfigReader::post_commands(const yaml_event_t &event)
{
    config.post_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.post_commands.push_back(line);
    }
}

void TestConfigReader::lammps_version(const yaml_event_t &event)
{
    config.lammps_version = (char *)event.data.scalar.value;
}

void TestConfigReader::date_generated(const yaml_event_t &event)
{
    config.date_generated = (char *)event.data.scalar.value;
}

void TestConfigReader::epsilon(const yaml_event_t &event)
{
    config.epsilon = atof((char *)event.data.scalar.value);
}

void TestConfigReader::input_file(const yaml_event_t &event)
{
    config.input_file = (char *)event.data.scalar.value;
}

void TestConfigReader::extract(const yaml_event_t &event)
{
    config.extract.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string name;
    int value;
    while (1) {
        data >> name >> value;
        if (data.eof()) break;
        config.extract.push_back(make_pair(name, value));
    }
}

void TestConfigReader::natoms(const yaml_event_t &event)
{
    config.natoms = atoi((char *)event.data.scalar.value);
}

void TestConfigReader::init_stress(const yaml_event_t &event)
{
    stress_t stress;
    sscanf((char *)event.data.scalar.value, "%lg %lg %lg %lg %lg %lg", &stress.xx, &stress.yy,
           &stress.zz, &stress.xy, &stress.xz, &stress.yz);
    config.init_stress = stress;
}

void TestConfigReader::run_stress(const yaml_event_t &event)
{
    stress_t stress;
    sscanf((char *)event.data.scalar.value, "%lg %lg %lg %lg %lg %lg", &stress.xx, &stress.yy,
           &stress.zz, &stress.xy, &stress.xz, &stress.yz);
    config.run_stress = stress;
}

void TestConfigReader::init_forces(const yaml_event_t &event)
{
    config.init_forces.clear();
    config.init_forces.resize(config.natoms + 1);
    std::stringstream data((const char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        int tag = 0;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.init_forces[tag] = xyz;
    }
}

void TestConfigReader::run_forces(const yaml_event_t &event)
{
    config.run_forces.clear();
    config.run_forces.resize(config.natoms + 1);
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        int tag;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.run_forces[tag] = xyz;
    }
}

void TestConfigReader::run_pos(const yaml_event_t &event)
{
    config.run_pos.clear();
    config.run_pos.resize(config.natoms + 1);
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        int tag;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.run_pos[tag] = xyz;
    }
}

void TestConfigReader::run_vel(const yaml_event_t &event)
{
    config.run_vel.clear();
    config.run_vel.resize(config.natoms + 1);
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        int tag;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.run_vel[tag] = xyz;
    }
}

void TestConfigReader::pair_style(const yaml_event_t &event)
{
    config.pair_style = (char *)event.data.scalar.value;
}

void TestConfigReader::pair_coeff(const yaml_event_t &event)
{
    config.pair_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.pair_coeff.push_back(line);
    }
}

void TestConfigReader::bond_style(const yaml_event_t &event)
{
    config.bond_style = (char *)event.data.scalar.value;
}

void TestConfigReader::bond_coeff(const yaml_event_t &event)
{
    config.bond_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.bond_coeff.push_back(line);
    }
}

void TestConfigReader::angle_style(const yaml_event_t &event)
{
    config.angle_style = (char *)event.data.scalar.value;
}

void TestConfigReader::angle_coeff(const yaml_event_t &event)
{
    config.angle_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.angle_coeff.push_back(line);
    }
}

void TestConfigReader::dihedral_style(const yaml_event_t &event)
{
    config.dihedral_style = (char *)event.data.scalar.value;
}

void TestConfigReader::dihedral_coeff(const yaml_event_t &event)
{
    config.dihedral_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.dihedral_coeff.push_back(line);
    }
}

void TestConfigReader::improper_style(const yaml_event_t &event)
{
    config.improper_style = (char *)event.data.scalar.value;
}

void TestConfigReader::improper_coeff(const yaml_event_t &event)
{
    config.improper_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.improper_coeff.push_back(line);
    }
}

void TestConfigReader::equilibrium(const yaml_event_t &event)
{
    std::stringstream data((char *)event.data.scalar.value);
    config.equilibrium.clear();
    double value;
    std::size_t num;
    data >> num;
    for (std::size_t i = 0; i < num; ++i) {
        data >> value;
        if (data.eof()) break;
        config.equilibrium.push_back(value);
    }
}

void TestConfigReader::init_vdwl(const yaml_event_t &event)
{
    config.init_vdwl = atof((char *)event.data.scalar.value);
}

void TestConfigReader::init_coul(const yaml_event_t &event)
{
    config.init_coul = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_vdwl(const yaml_event_t &event)
{
    config.run_vdwl = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_coul(const yaml_event_t &event)
{
    config.run_coul = atof((char *)event.data.scalar.value);
}

void TestConfigReader::init_energy(const yaml_event_t &event)
{
    config.init_energy = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_energy(const yaml_event_t &event)
{
    config.run_energy = atof((char *)event.data.scalar.value);
}

void TestConfigReader::global_scalar(const yaml_event_t &event)
{
    config.global_scalar = atof((char *)event.data.scalar.value);
}

void TestConfigReader::global_vector(const yaml_event_t &event)
{
    std::stringstream data((char *)event.data.scalar.value);
    config.global_vector.clear();
    double value;
    std::size_t num;
    data >> num;
    for (std::size_t i = 0; i < num; ++i) {
        data >> value;
        config.global_vector.push_back(value);
    }
}
