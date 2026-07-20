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

#include "test_config_reader.h"
#include "test_config.h"
#include "tokenizer.h"
#include "utils.h"
#include "yaml.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using LAMMPS_NS::Tokenizer;
using LAMMPS_NS::TokenizerException;
using LAMMPS_NS::ValueTokenizer;
using LAMMPS_NS::utils::split_words;

// convert the scalar value of a yaml event to a string
static std::string event_string(const yaml_event_t &event)
{
    return {(const char *)event.data.scalar.value};
}

// abort with a clear message on any malformed field. the yaml files are
// machine generated, so a parse error indicates a corrupted file and
// continuing would compare against garbage reference data.
[[noreturn]] static void parse_error(const std::string &what, const std::string &text)
{
    std::cerr << "Error parsing yaml data: " << what << "\nOffending text: " << text << std::endl;
    exit(2);
}

// parse a block of per-atom data with lines of "tag x y z" into a
// vector of coord_t indexed by the atom tag

static void parse_coord_block(const yaml_event_t &event, std::vector<coord_t> &block, int natoms)
{
    block.clear();
    block.resize(natoms + 1);
    for (const auto &line : Tokenizer(event_string(event), "\n").as_vector()) {
        try {
            ValueTokenizer values(line);
            int tag = values.next_int();
            if ((tag < 1) || (tag > natoms)) parse_error("atom tag out of range", line);
            block[tag].x = values.next_double();
            block[tag].y = values.next_double();
            block[tag].z = values.next_double();
        } catch (TokenizerException &e) {
            parse_error(e.what(), line);
        }
    }
}

// parse a block of per-atom data with lines of "tag x y z w" into a
// vector of coord4_t indexed by the atom tag

static void parse_coord4_block(const yaml_event_t &event, std::vector<coord4_t> &block, int natoms)
{
    block.clear();
    block.resize(natoms + 1);
    for (const auto &line : Tokenizer(event_string(event), "\n").as_vector()) {
        try {
            ValueTokenizer values(line);
            int tag = values.next_int();
            if ((tag < 1) || (tag > natoms)) parse_error("atom tag out of range", line);
            block[tag].x = values.next_double();
            block[tag].y = values.next_double();
            block[tag].z = values.next_double();
            block[tag].w = values.next_double();
        } catch (TokenizerException &e) {
            parse_error(e.what(), line);
        }
    }
}

// parse a block scalar with 6 stress tensor components

static void parse_stress(const yaml_event_t &event, stress_t &stress)
{
    try {
        ValueTokenizer values(event_string(event));
        stress.xx = values.next_double();
        stress.yy = values.next_double();
        stress.zz = values.next_double();
        stress.xy = values.next_double();
        stress.xz = values.next_double();
        stress.yz = values.next_double();
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}

// parse a scalar double value

static double parse_double(const yaml_event_t &event)
{
    try {
        return ValueTokenizer(event_string(event)).next_double();
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}

// split a block scalar into lines

static std::vector<std::string> parse_lines(const yaml_event_t &event)
{
    return Tokenizer(event_string(event), "\n").as_vector();
}

TestConfigReader::TestConfigReader(TestConfig &config) : config(config)
{
    consumers["lammps_version"] = &TestConfigReader::lammps_version;
    consumers["tags"]           = &TestConfigReader::tags;
    consumers["date_generated"] = &TestConfigReader::date_generated;
    consumers["epsilon"]        = &TestConfigReader::epsilon;
    consumers["skip_tests"]     = &TestConfigReader::skip_tests;
    consumers["prerequisites"]  = &TestConfigReader::prerequisites;
    consumers["pre_commands"]   = &TestConfigReader::pre_commands;
    consumers["post_commands"]  = &TestConfigReader::post_commands;
    consumers["input_file"]     = &TestConfigReader::input_file;
    consumers["input_coeffs"]   = &TestConfigReader::input_coeffs;
    consumers["extract"]        = &TestConfigReader::extract;
    consumers["natoms"]         = &TestConfigReader::natoms;
    consumers["init_stress"]    = &TestConfigReader::init_stress;
    consumers["run_stress"]     = &TestConfigReader::run_stress;
    consumers["init_forces"]    = &TestConfigReader::init_forces;
    consumers["run_forces"]     = &TestConfigReader::run_forces;
    consumers["run_pos"]        = &TestConfigReader::run_pos;
    consumers["run_vel"]        = &TestConfigReader::run_vel;
    consumers["run_torque"]     = &TestConfigReader::run_torque;
    consumers["init_mag_forces"] = &TestConfigReader::init_mag_forces;
    consumers["run_mag_forces"] = &TestConfigReader::run_mag_forces;
    consumers["run_spin"]       = &TestConfigReader::run_spin;
    consumers["timestep"]       = &TestConfigReader::timestep;

    consumers["pair_style"] = &TestConfigReader::pair_style;
    consumers["pair_coeff"] = &TestConfigReader::pair_coeff;
    consumers["init_vdwl"]  = &TestConfigReader::init_vdwl;
    consumers["init_coul"]  = &TestConfigReader::init_coul;
    consumers["run_vdwl"]   = &TestConfigReader::run_vdwl;
    consumers["run_coul"]   = &TestConfigReader::run_coul;

    consumers["global_scalar"] = &TestConfigReader::global_scalar;
    consumers["global_vector"] = &TestConfigReader::global_vector;
    consumers["global_array"]  = &TestConfigReader::global_array;
    consumers["peratom_data"]  = &TestConfigReader::peratom_data;
    consumers["local_data"]    = &TestConfigReader::local_data;

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
    for (const auto &word : split_words(event_string(event)))
        config.skip_tests.insert(word);
}

void TestConfigReader::prerequisites(const yaml_event_t &event)
{
    config.prerequisites.clear();
    ValueTokenizer data(event_string(event));

    while (data.has_next()) {
        std::string key = data.next_string();
        if (!data.has_next()) break;
        std::string value = data.next_string();
        config.prerequisites.emplace_back(key, value);
    }
}

void TestConfigReader::pre_commands(const yaml_event_t &event)
{
    config.pre_commands = parse_lines(event);
}

void TestConfigReader::post_commands(const yaml_event_t &event)
{
    config.post_commands = parse_lines(event);
}

void TestConfigReader::lammps_version(const yaml_event_t &event)
{
    config.lammps_version = event_string(event);
}

void TestConfigReader::date_generated(const yaml_event_t &event)
{
    config.date_generated = event_string(event);
}

void TestConfigReader::epsilon(const yaml_event_t &event)
{
    config.epsilon = parse_double(event);
}

void TestConfigReader::input_file(const yaml_event_t &event)
{
    config.input_file = event_string(event);
}

void TestConfigReader::input_coeffs(const yaml_event_t &event)
{
    config.input_coeffs = event_string(event);
}

void TestConfigReader::extract(const yaml_event_t &event)
{
    config.extract.clear();
    ValueTokenizer data(event_string(event));

    try {
        while (data.has_next()) {
            std::string name = data.next_string();
            if (!data.has_next()) break;
            int value = data.next_int();
            config.extract.emplace_back(name, value);
        }
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}

void TestConfigReader::natoms(const yaml_event_t &event)
{
    try {
        config.natoms = ValueTokenizer(event_string(event)).next_int();
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}

void TestConfigReader::init_stress(const yaml_event_t &event)
{
    parse_stress(event, config.init_stress);
}

void TestConfigReader::run_stress(const yaml_event_t &event)
{
    parse_stress(event, config.run_stress);
}

void TestConfigReader::init_forces(const yaml_event_t &event)
{
    parse_coord_block(event, config.init_forces, config.natoms);
}

void TestConfigReader::run_forces(const yaml_event_t &event)
{
    parse_coord_block(event, config.run_forces, config.natoms);
}

void TestConfigReader::run_pos(const yaml_event_t &event)
{
    parse_coord_block(event, config.run_pos, config.natoms);
}

void TestConfigReader::run_vel(const yaml_event_t &event)
{
    parse_coord_block(event, config.run_vel, config.natoms);
}

void TestConfigReader::run_torque(const yaml_event_t &event)
{
    parse_coord_block(event, config.run_torque, config.natoms);
}

void TestConfigReader::init_mag_forces(const yaml_event_t &event)
{
    parse_coord_block(event, config.init_mag_forces, config.natoms);
}

void TestConfigReader::run_mag_forces(const yaml_event_t &event)
{
    parse_coord_block(event, config.run_mag_forces, config.natoms);
}

void TestConfigReader::run_spin(const yaml_event_t &event)
{
    parse_coord4_block(event, config.run_spin, config.natoms);
}

void TestConfigReader::timestep(const yaml_event_t &event)
{
    config.timestep = parse_double(event);
}

void TestConfigReader::pair_style(const yaml_event_t &event)
{
    config.pair_style = event_string(event);
}

void TestConfigReader::pair_coeff(const yaml_event_t &event)
{
    config.pair_coeff = parse_lines(event);
}

void TestConfigReader::bond_style(const yaml_event_t &event)
{
    config.bond_style = event_string(event);
}

void TestConfigReader::bond_coeff(const yaml_event_t &event)
{
    config.bond_coeff = parse_lines(event);
}

void TestConfigReader::angle_style(const yaml_event_t &event)
{
    config.angle_style = event_string(event);
}

void TestConfigReader::angle_coeff(const yaml_event_t &event)
{
    config.angle_coeff = parse_lines(event);
}

void TestConfigReader::dihedral_style(const yaml_event_t &event)
{
    config.dihedral_style = event_string(event);
}

void TestConfigReader::dihedral_coeff(const yaml_event_t &event)
{
    config.dihedral_coeff = parse_lines(event);
}

void TestConfigReader::improper_style(const yaml_event_t &event)
{
    config.improper_style = event_string(event);
}

void TestConfigReader::improper_coeff(const yaml_event_t &event)
{
    config.improper_coeff = parse_lines(event);
}

void TestConfigReader::equilibrium(const yaml_event_t &event)
{
    config.equilibrium.clear();
    try {
        ValueTokenizer data(event_string(event));
        std::size_t num = data.next_int();
        for (std::size_t i = 0; i < num; ++i) {
            if (!data.has_next()) break;
            config.equilibrium.push_back(data.next_double());
        }
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}

void TestConfigReader::init_vdwl(const yaml_event_t &event)
{
    config.init_vdwl = parse_double(event);
}

void TestConfigReader::init_coul(const yaml_event_t &event)
{
    config.init_coul = parse_double(event);
}

void TestConfigReader::run_vdwl(const yaml_event_t &event)
{
    config.run_vdwl = parse_double(event);
}

void TestConfigReader::run_coul(const yaml_event_t &event)
{
    config.run_coul = parse_double(event);
}

void TestConfigReader::init_energy(const yaml_event_t &event)
{
    config.init_energy = parse_double(event);
}

void TestConfigReader::run_energy(const yaml_event_t &event)
{
    config.run_energy = parse_double(event);
}

void TestConfigReader::global_scalar(const yaml_event_t &event)
{
    config.global_scalar = parse_double(event);
}

void TestConfigReader::global_vector(const yaml_event_t &event)
{
    config.global_vector.clear();
    try {
        ValueTokenizer data(event_string(event));
        std::size_t num = data.next_int();
        for (std::size_t i = 0; i < num; ++i) {
            if (!data.has_next()) break;
            config.global_vector.push_back(data.next_double());
        }
    } catch (TokenizerException &e) {
        parse_error(e.what(), event_string(event));
    }
}


// parse a block of rows of white-space separated numbers into a
// vector of vectors, so rows may have different numbers of columns

static void parse_rows(const yaml_event_t &event, std::vector<std::vector<double>> &rows)
{
    rows.clear();
    for (const auto &line : Tokenizer(event_string(event), "\n").as_vector()) {
        try {
            ValueTokenizer values(line);
            std::vector<double> row;
            while (values.has_next())
                row.push_back(values.next_double());
            if (!row.empty()) rows.push_back(row);
        } catch (TokenizerException &e) {
            parse_error(e.what(), line);
        }
    }
}

void TestConfigReader::global_array(const yaml_event_t &event)
{
    parse_rows(event, config.global_array);
}

void TestConfigReader::peratom_data(const yaml_event_t &event)
{
    parse_rows(event, config.peratom_data);
}

void TestConfigReader::local_data(const yaml_event_t &event)
{
    parse_rows(event, config.local_data);
}

void TestConfigReader::tags(const yaml_event_t &event)
{
    config.tags = split_words(event_string(event));
}
