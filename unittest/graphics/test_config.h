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

#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct rgb_t {
    double r, g, b;
};

struct block_t {
    int y, x;
    rgb_t mean;
};

class TestConfig {
public:
    std::string lammps_version;
    std::string date_generated;
    std::string basename;
    std::vector<std::string> tags;
    std::vector<std::pair<std::string, std::string>> prerequisites;
    // scene setup: box/atoms/styles plus the dump image and dump_modify
    // commands under test.  the image file name must be given as
    // ${imagefile}; the test driver defines that variable.
    std::vector<std::string> setup_commands;
    // commands that trigger the render, usually just "run 0"
    std::vector<std::string> run_commands;
    // required MPI rank count: rendered images legitimately differ with the
    // number of ranks (depth-tie compositing, SSAO), so each reference is
    // only valid for one rank count
    int nprocs;
    // reference image dimensions
    int image_width, image_height;
    // uniform sampling block edge length and grid stride (pixels)
    int block_size, block_stride;
    // max allowed deviation of per-row/per-column mean RGB values (0-255 scale)
    double epsilon_projection;
    // max allowed deviation of sampled block mean RGB values (0-255 scale)
    double epsilon_blocks;
    // reference data: per-row and per-column mean RGB plus uniform grid blocks
    std::vector<rgb_t> row_means;
    std::vector<rgb_t> col_means;
    std::vector<block_t> sample_blocks;

    TestConfig() :
        lammps_version(""), date_generated(""), basename(""), nprocs(1), image_width(0),
        image_height(0), block_size(3), block_stride(16), epsilon_projection(0.5),
        epsilon_blocks(2.0)
    {
        tags.clear();
        prerequisites.clear();
        setup_commands.clear();
        run_commands.clear();
        row_means.clear();
        col_means.clear();
        sample_blocks.clear();
    }
    TestConfig(const TestConfig &)            = delete;
    TestConfig &operator=(const TestConfig &) = delete;

    [[nodiscard]] std::string tags_line() const
    {
        if (tags.size() > 0) {
            std::stringstream line;
            line << tags[0];
            for (size_t i = 1; i < tags.size(); i++)
                line << ", " << tags[i];
            return line.str();
        }
        return "generated";
    }
};

#endif
