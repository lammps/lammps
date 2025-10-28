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

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

#include "json.h"

using namespace LAMMPS_NS;
using ::testing::Eq;

#define STRINGIFY(val) XSTR(val)
#define XSTR(val) #val

// this tests a subset of the JSON class that is most relevant to LAMMPS

TEST(JSON, namespace)
{
    std::string expected = "nlohmann_lmp::json_abi";
    expected += "_v" STRINGIFY(LMP_NLOHMANN_JSON_VERSION_MAJOR);
    expected += "_" STRINGIFY(LMP_NLOHMANN_JSON_VERSION_MINOR);
    expected += "_" STRINGIFY(LMP_NLOHMANN_JSON_VERSION_PATCH);

    const std::string ns{STRINGIFY(LMP_NLOHMANN_JSON_NAMESPACE)};
    ASSERT_THAT(expected, Eq(ns));
}

TEST(JSON, serialize_deserialize)
{
    json j1;
    j1["pi"]      = 3.141;
    j1["happy"]   = true;
    j1["name"]    = "Niels";
    j1["nothing"] = nullptr;

    std::string expected = "{\"pi\":3.141,\"happy\":true,\"name\":\"Niels\",\"nothing\":null}";
    std::string dumped   = j1.dump(-1);
    ASSERT_THAT(expected, Eq(dumped));

    expected = "{\n  \"pi\": 3.141,\n  \"happy\": true,\n  \"name\": \"Niels\",\n  \"nothing\": null\n}";
    dumped   = j1.dump(2, ' ');
    ASSERT_THAT(expected, Eq(dumped));

    json j2 = json::parse(expected);
    ASSERT_TRUE(j1 == j2);
}

TEST(JSON, init_vs_incremental)
{
    json j1;
    j1["pi"]                   = 3.141;
    j1["happy"]                = true;
    j1["name"]                 = "Niels";
    j1["nothing"]              = nullptr;
    j1["answer"]["everything"] = 42;
    j1["list"]                 = {1, 0, 2};
    j1["object"]               = {{"currency", "USD"}, {"value", 42.99}};

    json j2 = {{"pi", 3.141},
               {"happy", true},
               {"name", "Niels"},
               {"nothing", nullptr},
               {"answer", {{"everything", 42}}},
               {"list", {1, 0, 2}},
               {"object", {{"currency", "USD"}, {"value", 42.99}}}};

    ASSERT_TRUE(j1 == j2);
}
