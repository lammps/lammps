/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

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

protected:
    void lammps_version(const yaml_event_t &event);
    void date_generated(const yaml_event_t &event);
    void tags(const yaml_event_t &event);
    void prerequisites(const yaml_event_t &event);
    void setup_commands(const yaml_event_t &event);
    void run_commands(const yaml_event_t &event);
    void nprocs(const yaml_event_t &event);
    void image_size(const yaml_event_t &event);
    void sampling(const yaml_event_t &event);
    void epsilon_projection(const yaml_event_t &event);
    void epsilon_blocks(const yaml_event_t &event);
    void row_means(const yaml_event_t &event);
    void col_means(const yaml_event_t &event);
    void sample_blocks(const yaml_event_t &event);

private:
    void parse_rgb_lines(const yaml_event_t &event, std::vector<rgb_t> &values);
};

#endif
