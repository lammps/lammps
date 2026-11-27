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
#include "utils.h"
#include "yaml.h"

#include <sstream>
#include <string>

using LAMMPS_NS::utils::split_words;

TestConfigReader::TestConfigReader(TestConfig &config) : config(config)
{
    consumers["lammps_version"]     = &TestConfigReader::lammps_version;
    consumers["date_generated"]     = &TestConfigReader::date_generated;
    consumers["tags"]               = &TestConfigReader::tags;
    consumers["prerequisites"]      = &TestConfigReader::prerequisites;
    consumers["setup_commands"]     = &TestConfigReader::setup_commands;
    consumers["run_commands"]       = &TestConfigReader::run_commands;
    consumers["nprocs"]             = &TestConfigReader::nprocs;
    consumers["image_size"]         = &TestConfigReader::image_size;
    consumers["sampling"]           = &TestConfigReader::sampling;
    consumers["epsilon_projection"] = &TestConfigReader::epsilon_projection;
    consumers["epsilon_blocks"]     = &TestConfigReader::epsilon_blocks;
    consumers["row_means"]          = &TestConfigReader::row_means;
    consumers["col_means"]          = &TestConfigReader::col_means;
    consumers["sample_blocks"]      = &TestConfigReader::sample_blocks;
}

void TestConfigReader::lammps_version(const yaml_event_t &event)
{
    config.lammps_version = (char *)event.data.scalar.value;
}

void TestConfigReader::date_generated(const yaml_event_t &event)
{
    config.date_generated = (char *)event.data.scalar.value;
}

void TestConfigReader::tags(const yaml_event_t &event)
{
    config.tags = split_words((char *)event.data.scalar.value);
}

void TestConfigReader::prerequisites(const yaml_event_t &event)
{
    config.prerequisites.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string key, value;

    while (true) {
        data >> key >> value;
        if (data.eof()) break;
        config.prerequisites.emplace_back(key, value);
    }
}

void TestConfigReader::setup_commands(const yaml_event_t &event)
{
    config.setup_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n'))
        config.setup_commands.push_back(line);
}

void TestConfigReader::run_commands(const yaml_event_t &event)
{
    config.run_commands.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n'))
        config.run_commands.push_back(line);
}

void TestConfigReader::nprocs(const yaml_event_t &event)
{
    config.nprocs = atoi((char *)event.data.scalar.value);
}

void TestConfigReader::image_size(const yaml_event_t &event)
{
    std::stringstream data((char *)event.data.scalar.value);
    data >> config.image_width >> config.image_height;
}

void TestConfigReader::sampling(const yaml_event_t &event)
{
    std::stringstream data((char *)event.data.scalar.value);
    data >> config.block_size >> config.block_stride;
}

void TestConfigReader::epsilon_projection(const yaml_event_t &event)
{
    config.epsilon_projection = atof((char *)event.data.scalar.value);
}

void TestConfigReader::epsilon_blocks(const yaml_event_t &event)
{
    config.epsilon_blocks = atof((char *)event.data.scalar.value);
}

void TestConfigReader::parse_rgb_lines(const yaml_event_t &event, std::vector<rgb_t> &values)
{
    values.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        std::stringstream item(line);
        rgb_t rgb;
        item >> rgb.r >> rgb.g >> rgb.b;
        if (!item.fail()) values.push_back(rgb);
    }
}

void TestConfigReader::row_means(const yaml_event_t &event)
{
    parse_rgb_lines(event, config.row_means);
}

void TestConfigReader::col_means(const yaml_event_t &event)
{
    parse_rgb_lines(event, config.col_means);
}

void TestConfigReader::sample_blocks(const yaml_event_t &event)
{
    config.sample_blocks.clear();
    std::stringstream data((char *)event.data.scalar.value);
    std::string line;

    while (std::getline(data, line, '\n')) {
        std::stringstream item(line);
        block_t block;
        item >> block.y >> block.x >> block.mean.r >> block.mean.g >> block.mean.b;
        if (!item.fail()) config.sample_blocks.push_back(block);
    }
}
