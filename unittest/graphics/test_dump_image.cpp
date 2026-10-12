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

// Test driver for dump image and GRAPHICS package rendering.
//
// A test YAML file describes a scene (setup commands including the dump
// image and dump_modify commands under test) and stores sampled reference
// pixel data from the rendered image directly in the YAML file:
// per-row and per-column mean RGB values (full-image coverage) plus mean
// RGB values of small blocks on a uniform grid (failure localization).
// See README.md in this folder for the YAML format and how the sampling
// scheme was validated.

#include "test_config.h"
#include "test_config_reader.h"
#include "yaml_writer.h"

#include "fmt/format.h"
#include "info.h"
#include "input.h"
#include "lammps.h"
#include "platform.h"
#include "utils.h"
#include "version.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <mpi.h>

#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using LAMMPS_NS::Info;
using LAMMPS_NS::LAMMPS;
using LAMMPS_NS::utils::split_words;
using LAMMPS_NS::utils::trim;
namespace platform = LAMMPS_NS::platform;

// common read-only data available to all tests
TestConfig test_config;
bool verbose = false;
#if defined(TEST_INPUT_FOLDER)
#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val
std::string INPUT_FOLDER = STRINGIFY(TEST_INPUT_FOLDER);
#undef STRINGIFY
#undef XSTR
#else
std::string INPUT_FOLDER = ".";
#endif

/* ---------------------------------------------------------------------- */
// minimal 8-bit RGB image with PPM (P6) reader

struct Image {
    int width  = 0;
    int height = 0;
    std::vector<unsigned char> data;

    double val(int y, int x, int c) const { return (double)data[(y * width + x) * 3 + c]; }
};

static bool read_ppm(const std::string &filename, Image &img)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) return false;
    std::string raw((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    // parse header: magic, width, height, maxval; '#' starts a comment line
    std::size_t pos = 0;
    std::vector<std::string> tokens;
    while ((pos < raw.size()) && (tokens.size() < 4)) {
        if (isspace(raw[pos])) {
            ++pos;
        } else if (raw[pos] == '#') {
            while ((pos < raw.size()) && (raw[pos] != '\n'))
                ++pos;
        } else {
            std::size_t begin = pos;
            while ((pos < raw.size()) && !isspace(raw[pos]))
                ++pos;
            tokens.push_back(raw.substr(begin, pos - begin));
        }
    }
    if ((tokens.size() < 4) || (tokens[0] != "P6") || (tokens[3] != "255")) return false;
    ++pos; // single whitespace character after maxval
    img.width  = std::stoi(tokens[1]);
    img.height = std::stoi(tokens[2]);
    const std::size_t npixels = (std::size_t)img.width * img.height * 3;
    if (raw.size() - pos < npixels) return false;
    img.data.assign(raw.begin() + pos, raw.begin() + pos + npixels);
    return true;
}

/* ---------------------------------------------------------------------- */
// sampling: row/column mean RGB projections and uniform grid of blocks

static std::vector<rgb_t> project_rows(const Image &img)
{
    std::vector<rgb_t> means(img.height);
    for (int y = 0; y < img.height; ++y) {
        double sum[3] = {0.0, 0.0, 0.0};
        for (int x = 0; x < img.width; ++x)
            for (int c = 0; c < 3; ++c)
                sum[c] += img.val(y, x, c);
        means[y] = {sum[0] / img.width, sum[1] / img.width, sum[2] / img.width};
    }
    return means;
}

static std::vector<rgb_t> project_cols(const Image &img)
{
    std::vector<rgb_t> means(img.width);
    for (int x = 0; x < img.width; ++x) {
        double sum[3] = {0.0, 0.0, 0.0};
        for (int y = 0; y < img.height; ++y)
            for (int c = 0; c < 3; ++c)
                sum[c] += img.val(y, x, c);
        means[x] = {sum[0] / img.height, sum[1] / img.height, sum[2] / img.height};
    }
    return means;
}

static rgb_t block_mean(const Image &img, int cy, int cx, int size)
{
    const int half = size / 2;
    double sum[3]  = {0.0, 0.0, 0.0};
    for (int y = cy - half; y <= cy + half; ++y)
        for (int x = cx - half; x <= cx + half; ++x)
            for (int c = 0; c < 3; ++c)
                sum[c] += img.val(y, x, c);
    const double norm = (double)size * size;
    return {sum[0] / norm, sum[1] / norm, sum[2] / norm};
}

static std::vector<std::pair<int, int>> block_centers(int width, int height, int size, int stride)
{
    std::vector<std::pair<int, int>> centers;
    const int offset = size / 2 + 1;
    for (int y = offset; y < height - offset; y += stride)
        for (int x = offset; x < width - offset; x += stride)
            centers.emplace_back(y, x);
    return centers;
}

static double max_rgb_diff(const rgb_t &one, const rgb_t &two)
{
    double diff = fabs(one.r - two.r);
    diff        = std::max(diff, fabs(one.g - two.g));
    diff        = std::max(diff, fabs(one.b - two.b));
    return diff;
}

/* ---------------------------------------------------------------------- */
// rendering: set up a LAMMPS instance, run the scene, locate the image

// image files are written as <basename>-step<N>.ppm in the working directory

static void purge_images(const TestConfig &cfg)
{
    const std::string prefix = cfg.basename + "-step";
    for (const auto &file : platform::list_directory(".")) {
        if ((file.rfind(prefix, 0) == 0) && LAMMPS_NS::utils::strmatch(file, "\\.ppm$"))
            platform::unlink(file);
    }
}

// find the rendered image with the highest timestep number

static std::string find_image(const TestConfig &cfg)
{
    const std::string prefix = cfg.basename + "-step";
    std::string found;
    long beststep = -1;
    for (const auto &file : platform::list_directory(".")) {
        if ((file.rfind(prefix, 0) != 0) || !LAMMPS_NS::utils::strmatch(file, "\\.ppm$")) continue;
        const std::string num = file.substr(prefix.size(), file.size() - prefix.size() - 4);
        if (!LAMMPS_NS::utils::is_integer(num)) continue;
        const long step = std::stol(num);
        if (step > beststep) {
            beststep = step;
            found    = file;
        }
    }
    return found;
}

static LAMMPS *init_lammps(const TestConfig &cfg)
{
    LAMMPS::argv args = {"DumpImageTest", "-log", "none", "-nocite"};
    if (verbose) {
        args.emplace_back("-echo");
        args.emplace_back("screen");
    } else {
        args.emplace_back("-screen");
        args.emplace_back("none");
    }
    auto *lmp = new LAMMPS(args, MPI_COMM_WORLD);

    // check if the GRAPHICS package and all prerequisite styles are available
    Info info(lmp);
    int nfail = 0;
    if (!info.has_style("dump", "image")) ++nfail;
    for (const auto &prerequisite : cfg.prerequisites) {
        if (!info.has_style(prerequisite.first, prerequisite.second)) ++nfail;
    }
    if (nfail > 0) {
        delete lmp;
        return nullptr;
    }

    lmp->input->one(fmt::format("variable input_dir string \"{}\"", INPUT_FOLDER));
    lmp->input->one(fmt::format("variable imagefile string {}-step*.ppm", cfg.basename));
    for (const auto &line : cfg.setup_commands)
        lmp->input->one(line);
    for (const auto &line : cfg.run_commands)
        lmp->input->one(line);
    return lmp;
}

// render the scene of the current test config and read back the image.
// returns false when prerequisite styles are missing.  all file handling
// happens on MPI rank 0; other ranks return an empty image.

static bool render_scene(const TestConfig &cfg, Image &img, std::string &imagefile)
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0) purge_images(cfg);
    MPI_Barrier(MPI_COMM_WORLD);

    LAMMPS *lmp = init_lammps(cfg);
    int okay    = (lmp != nullptr) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    // delete the LAMMPS instance so all files are flushed and closed
    delete lmp;
    if (!okay) return false;

    if (me == 0) {
        imagefile = find_image(cfg);
        if (!imagefile.empty() && !read_ppm(imagefile, img)) img = Image();
    }
    return true;
}

/* ---------------------------------------------------------------------- */

TEST(DumpImage, compare)
{
    int me, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // rendered images legitimately depend on the rank count: pixel-level
    // depth ties composite differently and SSAO slices its RNG stream by rank
    if (nprocs != test_config.nprocs)
        GTEST_SKIP() << "test reference requires exactly " << test_config.nprocs
                     << " MPI ranks, running on " << nprocs;

    Image img;
    std::string imagefile;
    if (!render_scene(test_config, img, imagefile)) {
        GTEST_SKIP() << "One or more prerequisite styles are not available in this LAMMPS build";
    }

    // image checks are done on rank 0 only, which writes the dump image file
    if (me != 0) return;

    ASSERT_FALSE(imagefile.empty()) << "no rendered image file found";
    ASSERT_GT(img.width, 0) << "could not read image file " << imagefile;
    ASSERT_EQ(img.width, test_config.image_width);
    ASSERT_EQ(img.height, test_config.image_height);
    ASSERT_EQ(test_config.row_means.size(), (std::size_t)img.height);
    ASSERT_EQ(test_config.col_means.size(), (std::size_t)img.width);

    // full-image coverage: per-row and per-column mean RGB projections

    double maxdiff = 0.0;
    int worst      = -1;
    auto rows      = project_rows(img);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const double diff = max_rgb_diff(rows[i], test_config.row_means[i]);
        if (diff > maxdiff) {
            maxdiff = diff;
            worst   = (int)i;
        }
    }
    EXPECT_LE(maxdiff, test_config.epsilon_projection)
        << "row mean RGB projections differ most at row " << worst;

    maxdiff   = 0.0;
    worst     = -1;
    auto cols = project_cols(img);
    for (std::size_t i = 0; i < cols.size(); ++i) {
        const double diff = max_rgb_diff(cols[i], test_config.col_means[i]);
        if (diff > maxdiff) {
            maxdiff = diff;
            worst   = (int)i;
        }
    }
    EXPECT_LE(maxdiff, test_config.epsilon_projection)
        << "column mean RGB projections differ most at column " << worst;

    // failure localization: mean RGB of small blocks on a uniform grid

    maxdiff      = 0.0;
    int worsty   = -1;
    int worstx   = -1;
    for (const auto &block : test_config.sample_blocks) {
        const auto mean   = block_mean(img, block.y, block.x, test_config.block_size);
        const double diff = max_rgb_diff(mean, block.mean);
        if (diff > maxdiff) {
            maxdiff = diff;
            worsty  = block.y;
            worstx  = block.x;
        }
    }
    EXPECT_LE(maxdiff, test_config.epsilon_blocks)
        << "sample blocks differ most at pixel y=" << worsty << " x=" << worstx;

    // keep the rendered image for inspection when the comparison failed

    if (::testing::Test::HasFailure()) {
        const std::string keep = test_config.basename + ".failed.ppm";
        std::rename(imagefile.c_str(), keep.c_str());
        std::cerr << "Rendered image kept as: " << keep << std::endl;
    }
    purge_images(test_config);
}

/* ---------------------------------------------------------------------- */
// generate reference sample data from the rendered scene into a new YAML file

static void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    int me, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Image img;
    std::string imagefile;
    if (!render_scene(config, img, imagefile)) {
        if (me == 0) {
            std::cerr << "One or more prerequisite styles are not available "
                         "in this LAMMPS configuration:\n";
            for (const auto &prerequisite : config.prerequisites)
                std::cerr << "  " << prerequisite.first << " style " << prerequisite.second
                          << "\n";
        }
        return;
    }
    if (me != 0) return;
    if (imagefile.empty() || (img.width == 0)) {
        std::cerr << "Failed to render or read an image for this scene\n";
        return;
    }

    YamlWriter writer(outfile);

    // header
    writer.emit("lammps_version", LAMMPS_VERSION);
    writer.emit("tags", config.tags_line());
    time_t now = time(nullptr);
    writer.emit("date_generated", trim(ctime(&now)));
    writer.emit("epsilon_projection", config.epsilon_projection);
    writer.emit("epsilon_blocks", config.epsilon_blocks);
    // the reference data is only valid for the rank count it was created with
    writer.emit("nprocs", nprocs);

    std::string block;
    for (const auto &prerequisite : config.prerequisites)
        block += fmt::format("{} {}\n", prerequisite.first, prerequisite.second);
    writer.emit_block("prerequisites", block);

    block.clear();
    for (const auto &line : config.setup_commands)
        block += line + "\n";
    writer.emit_block("setup_commands", block);

    block.clear();
    for (const auto &line : config.run_commands)
        block += line + "\n";
    writer.emit_block("run_commands", block);

    writer.emit("image_size", fmt::format("{} {}", img.width, img.height));
    writer.emit("sampling", fmt::format("{} {}", config.block_size, config.block_stride));

    block.clear();
    for (const auto &mean : project_rows(img))
        block += fmt::format("{:.2f} {:.2f} {:.2f}\n", mean.r, mean.g, mean.b);
    writer.emit_block("row_means", block);

    block.clear();
    for (const auto &mean : project_cols(img))
        block += fmt::format("{:.2f} {:.2f} {:.2f}\n", mean.r, mean.g, mean.b);
    writer.emit_block("col_means", block);

    block.clear();
    for (const auto &center :
         block_centers(img.width, img.height, config.block_size, config.block_stride)) {
        const auto mean = block_mean(img, center.first, center.second, config.block_size);
        block += fmt::format("{} {} {:.2f} {:.2f} {:.2f}\n", center.first, center.second, mean.r,
                             mean.g, mean.b);
    }
    writer.emit_block("sample_blocks", block);

    purge_images(config);
}

/* ---------------------------------------------------------------------- */

static bool read_yaml_file(const char *infile, TestConfig &config)
{
    auto reader = TestConfigReader(config);
    if (reader.parse_file(infile)) return false;
    config.basename = reader.get_basename();
    return true;
}

static void usage(std::ostream &out, const char *name)
{
    out << "usage: " << name << " <testfile.yaml> [OPTIONS]\n\n"
        << "Available options:\n"
        << "  -g <newfile.yaml>  regenerate sample data in new YAML file\n"
        << "  -u                 update sample data in original YAML file\n"
        << "  -d <folder>        set folder for scene input files\n"
        << "  -v                 run tests with verbose output\n"
        << std::endl;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (argc < 2) {
        usage(std::cerr, argv[0]);
        MPI_Finalize();
        return 1;
    }

    if (!read_yaml_file(argv[1], test_config)) {
        std::cerr << "Error parsing yaml file: " << argv[1] << std::endl;
        MPI_Finalize();
        return 2;
    }

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (const auto &arg : env) {
            if (arg == "-u") {
                generate_yaml_file(argv[1], test_config);
                MPI_Finalize();
                return 0;
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
