/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stl_reader.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using ::testing::DoubleEq;

namespace {

// write the given text to a file and return its name

std::string write_text(const std::string &name, const std::string &content)
{
    FILE *fp = fopen(name.c_str(), "wb");
    fwrite(content.data(), 1, content.size(), fp);
    fclose(fp);
    return name;
}

// write a minimal binary STL file with the given triangles (12 floats each:
// normal + 3 vertices) and an 80-byte title header

std::string write_binary(const std::string &name, const std::string &title,
                         const std::vector<std::array<float, 12>> &tris)
{
    FILE *fp = fopen(name.c_str(), "wb");
    char header[80] = {0};
    for (std::size_t i = 0; (i < title.size()) && (i < 79); ++i) header[i] = title[i];
    fwrite(header, 1, 80, fp);
    auto ntri = (uint32_t) tris.size();
    fwrite(&ntri, sizeof(ntri), 1, fp);
    for (const auto &t : tris) {
        fwrite(t.data(), sizeof(float), 12, fp);
        uint16_t attr = 0;
        fwrite(&attr, sizeof(attr), 1, fp);
    }
    fclose(fp);
    return name;
}

// a single reference triangle: normal (0,0,1), unit right triangle in z=0 plane

const char *ascii_one =
    "solid mycube\n"
    "  facet normal 0.0 0.0 1.0\n"
    "    outer loop\n"
    "      vertex 0.0 0.0 0.0\n"
    "      vertex 1.0 0.0 0.0\n"
    "      vertex 0.0 1.0 0.0\n"
    "    endloop\n"
    "  endfacet\n"
    "endsolid mycube\n";

}    // namespace

TEST(STLReader, ascii_single_triangle)
{
    auto file = write_text("test_stl_ascii.stl", ascii_one);
    std::string title;
    auto tris = STLReader::parse(file, &title);
    remove(file.c_str());

    ASSERT_EQ(tris.size(), 1u);
    EXPECT_EQ(title, "mycube");
    EXPECT_THAT(tris[0].normal[2], DoubleEq(1.0));
    EXPECT_THAT(tris[0].vert[1][0], DoubleEq(1.0));
    EXPECT_THAT(tris[0].vert[2][1], DoubleEq(1.0));
    EXPECT_THAT(tris[0].vert[0][0], DoubleEq(0.0));
}

TEST(STLReader, ascii_two_triangles)
{
    std::string two = std::string(ascii_one);
    // insert a second facet before endsolid
    two = "solid two\n"
          "facet normal 0 0 1\n outer loop\n"
          " vertex 0 0 0\n vertex 1 0 0\n vertex 0 1 0\n endloop\n endfacet\n"
          "facet normal 0 0 -1\n outer loop\n"
          " vertex 0 0 1\n vertex 1 0 1\n vertex 0 1 1\n endloop\n endfacet\n"
          "endsolid two\n";
    auto file = write_text("test_stl_two.stl", two);
    auto tris = STLReader::parse(file);
    remove(file.c_str());

    ASSERT_EQ(tris.size(), 2u);
    EXPECT_THAT(tris[1].normal[2], DoubleEq(-1.0));
    EXPECT_THAT(tris[1].vert[0][2], DoubleEq(1.0));
}

TEST(STLReader, ascii_empty_solid_name)
{
    auto file = write_text("test_stl_noname.stl",
                           "solid\n"
                           "facet normal 0 0 1\n outer loop\n"
                           " vertex 0 0 0\n vertex 1 0 0\n vertex 0 1 0\n endloop\n endfacet\n"
                           "endsolid\n");
    std::string title;
    auto tris = STLReader::parse(file, &title);
    remove(file.c_str());

    ASSERT_EQ(tris.size(), 1u);
    EXPECT_EQ(title, "");
}

TEST(STLReader, binary_single_triangle)
{
    std::vector<std::array<float, 12>> data = {
        {{0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f}}};
    auto file = write_binary("test_stl_bin.stl", "binary header", data);
    std::string title;
    auto tris = STLReader::parse(file, &title);
    remove(file.c_str());

    ASSERT_EQ(tris.size(), 1u);
    EXPECT_EQ(title, "binary header");
    EXPECT_THAT(tris[0].normal[2], DoubleEq(1.0));
    EXPECT_THAT(tris[0].vert[1][0], DoubleEq(1.0));
    EXPECT_THAT(tris[0].vert[2][1], DoubleEq(1.0));
}

TEST(STLReader, missing_file)
{
    EXPECT_THROW(STLReader::parse("this_file_does_not_exist.stl"), STLReaderException);
}

TEST(STLReader, not_an_stl_file)
{
    auto file = write_text("test_stl_bad.stl", "this is not an STL file\n");
    EXPECT_THROW(STLReader::parse(file), STLReaderException);
    remove(file.c_str());
}

TEST(STLReader, malformed_missing_endloop)
{
    auto file = write_text("test_stl_malformed.stl",
                           "solid bad\n"
                           "facet normal 0 0 1\n outer loop\n"
                           " vertex 0 0 0\n vertex 1 0 0\n vertex 0 1 0\n endfacet\n"
                           "endsolid bad\n");
    EXPECT_THROW(STLReader::parse(file), STLReaderException);
    remove(file.c_str());
}

TEST(STLReader, malformed_bad_vertex_keyword)
{
    auto file = write_text("test_stl_badvert.stl",
                           "solid bad\n"
                           "facet normal 0 0 1\n outer loop\n"
                           " vertex 0 0 0\n point 1 0 0\n vertex 0 1 0\n endloop\n endfacet\n"
                           "endsolid bad\n");
    EXPECT_THROW(STLReader::parse(file), STLReaderException);
    remove(file.c_str());
}
