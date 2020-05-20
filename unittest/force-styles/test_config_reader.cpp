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

#include "test_config.h"
#include "test_config_reader.h"
#include "yaml_reader.h"
#include "yaml.h"

#include <cstring>
#include <cstdlib>

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

TestConfigReader::TestConfigReader(TestConfig & config)
        : YamlReader(), config(config) {
    consumers["lammps_version"] = &TestConfigReader::lammps_version;
    consumers["date_generated"] = &TestConfigReader::date_generated;
    consumers["epsilon"]        = &TestConfigReader::epsilon;
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

    consumers["bond_style"]     = &TestConfigReader::bond_style;
    consumers["bond_coeff"]     = &TestConfigReader::bond_coeff;
    consumers["init_energy"]    = &TestConfigReader::init_energy;
    consumers["run_energy"]     = &TestConfigReader::run_energy;

    consumers["angle_style"]    = &TestConfigReader::angle_style;
    consumers["angle_coeff"]    = &TestConfigReader::angle_coeff;
    consumers["init_energy"]    = &TestConfigReader::init_energy;
    consumers["run_energy"]     = &TestConfigReader::run_energy;
}

void TestConfigReader::prerequisites(const yaml_event_t & event) {
    config.prerequisites.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while(std::getline(data, line, '\n')) {
        std::size_t found = line.find_first_of(" \t");
        std::string key = line.substr(0,found);
        found = line.find_first_not_of(" \t",found);
        // skip invalid data
        if (found == std::string::npos) {
            std::cerr << "Skipping invalid prerequisite line:\n"
                      << line << std::endl;
            continue;
        }
        std::string value = line.substr(found,line.find_first_of(" \t",found));
        config.prerequisites.push_back(std::pair<std::string,std::string>(key,value));
    }
}
void TestConfigReader::pre_commands(const yaml_event_t & event) {
    config.pre_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while(std::getline(data, line, '\n')) {
        config.pre_commands.push_back(line);
    }
}

void TestConfigReader::post_commands(const yaml_event_t & event) {
    config.post_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.post_commands.push_back(line);
    }
}

void TestConfigReader::lammps_version(const yaml_event_t & event) {
    config.lammps_version = (char *)event.data.scalar.value;
}

void TestConfigReader::date_generated(const yaml_event_t & event) {
    config.date_generated = (char *)event.data.scalar.value;
}

void TestConfigReader::epsilon(const yaml_event_t & event) {
    config.epsilon = atof((char *)event.data.scalar.value);
}

void TestConfigReader::input_file(const yaml_event_t & event) {
    config.input_file = (char *)event.data.scalar.value;
}

void TestConfigReader::extract(const yaml_event_t & event) {
    config.extract.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        std::size_t found = line.find_first_of(" \t");
        std::pair<std::string,int> data;
        data.first = line.substr(0,found);
        data.second = atoi(line.substr(found).c_str());
        config.extract.push_back(data);
    }
}

void TestConfigReader::natoms(const yaml_event_t & event) {
    config.natoms = atoi((char *)event.data.scalar.value);
}

void TestConfigReader::init_stress(const yaml_event_t & event) {
    stress_t stress;
    sscanf((char *)event.data.scalar.value,
           "%lg %lg %lg %lg %lg %lg",
           &stress.xx, &stress.yy, &stress.zz,
           &stress.xy, &stress.xz, &stress.yz);
    config.init_stress = stress;
}

void TestConfigReader::run_stress(const yaml_event_t & event) {
    stress_t stress;
    sscanf((char *)event.data.scalar.value,
           "%lg %lg %lg %lg %lg %lg",
           &stress.xx, &stress.yy, &stress.zz,
           &stress.xy, &stress.xz, &stress.yz);
    config.run_stress = stress;
}

void TestConfigReader::init_forces(const yaml_event_t & event) {
    config.init_forces.clear();
    config.init_forces.resize(config.natoms+1);
    std::stringstream data((const char*)event.data.scalar.value);
    std::string line;

    while(std::getline(data, line, '\n')) {
        int tag = 0;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.init_forces[tag] = xyz;
    }
}

void TestConfigReader::run_forces(const yaml_event_t & event) {
    config.run_forces.clear();
    config.run_forces.resize(config.natoms+1);
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while(std::getline(data, line, '\n')) {
        int tag;
        coord_t xyz;
        sscanf(line.c_str(), "%d %lg %lg %lg", &tag, &xyz.x, &xyz.y, &xyz.z);
        config.run_forces[tag] = xyz;
    }
}

void TestConfigReader::pair_style(const yaml_event_t & event) {
        config.pair_style = (char *)event.data.scalar.value;
}

void TestConfigReader::pair_coeff(const yaml_event_t & event) {
    config.pair_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.pair_coeff.push_back(line);
    }
}

void TestConfigReader::bond_style(const yaml_event_t & event) {
        config.bond_style = (char *)event.data.scalar.value;
}

void TestConfigReader::bond_coeff(const yaml_event_t & event) {
    config.bond_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.bond_coeff.push_back(line);
    }
}

void TestConfigReader::angle_style(const yaml_event_t & event) {
        config.angle_style = (char *)event.data.scalar.value;
}

void TestConfigReader::angle_coeff(const yaml_event_t & event) {
    config.angle_coeff.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        config.angle_coeff.push_back(line);
    }
}

void TestConfigReader::init_vdwl(const yaml_event_t & event) {
    config.init_vdwl = atof((char *)event.data.scalar.value);
}

void TestConfigReader::init_coul(const yaml_event_t & event) {
    config.init_coul = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_vdwl(const yaml_event_t & event) {
    config.run_vdwl = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_coul(const yaml_event_t & event) {
    config.run_coul = atof((char *)event.data.scalar.value);
}

void TestConfigReader::init_energy(const yaml_event_t & event) {
    config.init_energy = atof((char *)event.data.scalar.value);
}

void TestConfigReader::run_energy(const yaml_event_t & event) {
    config.run_energy = atof((char *)event.data.scalar.value);
}

