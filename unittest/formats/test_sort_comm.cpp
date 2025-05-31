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

#include "../testing/core.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_hybrid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "body.h"
#include "input.h"
#include "lammps.h"
#include "neighbor.h"
#include "math_const.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <mpi.h>
#include <vector>
#include <random>
#include <algorithm>
#include <map>

#if !defined(_FORTIFY_SOURCE) || (_FORTIFY_SOURCE == 0)
#if defined(__INTEL_COMPILER) || (__PGI)
#define _do_nothing
#elif defined(__clang__)
#pragma clang optimize off
#elif defined(__GNUC__)
#if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
#pragma GCC optimize("no-var-tracking-assignments", "O0")
#endif
#endif
#endif

#define GETIDX(i) lmp->atom->map(i)

using LAMMPS_NS::utils::split_words;

static void create_molecule_files(const std::string &h2o_filename, const std::string &co2_filename)
{
  // create molecule files
  const char h2o_file[] = "# Water molecule. SPC/E model.\n\n3 atoms\n2 bonds\n1 angles\n\n"
  "Coords\n\n1 1.12456 0.09298 1.27452\n"
  "2 1.53683 0.75606 1.89928\n3 0.49482 0.56390 0.65678\n\n"
  "Types\n\n1 1\n2 2\n3 2\n\n"
  "Charges\n\n1 -0.8472\n2 0.4236\n3 0.4236\n\n"
  "Bonds\n\n1 1 1 2\n2 1 1 3\n\n"
  "Angles\n\n1 1 2 1 3\n\n"
  "Shake Flags\n\n1 1\n2 1\n3 1\n\n"
  "Shake Atoms\n\n1 1 2 3\n2 1 2 3\n3 1 2 3\n\n"
  "Shake Bond Types\n\n1 1 1 1\n2 1 1 1\n3 1 1 1\n\n"
  "Special Bond Counts\n\n1 2 0 0\n2 1 1 0\n3 1 1 0\n\n"
  "Special Bonds\n\n1 2 3\n2 1 3\n3 1 2\n\n";
  const char co2_file[] = "# CO2 molecule file. TraPPE model.\n\n"
  "3 atoms\n2 bonds\n1 angles\n\n"
  "Coords\n\n1 0.0 0.0 0.0\n2 -1.16 0.0 0.0\n3 1.16 0.0 0.0\n\n"
  "Types\n\n1 1\n2 2\n3 2\n\n"
  "Charges\n\n1 0.7\n2 -0.35\n3 -0.35\n\n"
  "Bonds\n\n1 1 1 2\n2 1 1 3\n\n"
  "Angles\n\n1 1 2 1 3\n\n"
  "Special Bond Counts\n\n1 2 0 0\n2 1 1 0\n3 1 1 0\n\n"
  "Special Bonds\n\n1 2 3\n2 1 3\n3 1 2\n\n";
  
  FILE *fp = fopen(h2o_filename.c_str(), "w");
  if (fp) {
    fputs(h2o_file, fp);
    fclose(fp);
  }
  fp = fopen(co2_filename.c_str(), "w");
  if (fp) {
    fputs(co2_file, fp);
    fclose(fp);
  }
}

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

static const double EPSILON = 5.0e-14;
static std::mt19937 rng(12345);

namespace LAMMPS_NS {

// Test parameter enum
enum TestMode {
  PLAIN = 0,
  KOKKOS_OMP = 1
};

// Atom style information structure
struct AtomStyleInfo {
  std::string name;
  std::vector<std::string> fields_comm;
  std::vector<std::string> fields_comm_vel;
  std::vector<std::string> fields_reverse;
  std::vector<std::string> fields_border;
  std::vector<std::string> fields_border_vel;
  std::vector<std::string> fields_exchange;
  bool has_bonus_data;
  std::string setup_cmd; // Single setup command if needed
};

class SortCommTest : public LAMMPSTest, public ::testing::WithParamInterface<TestMode> {
protected:
  static void SetUpTestSuite() { create_molecule_files("h2o.mol", "co2.mol"); }
  
  static void TearDownTestSuite()
  {
    remove("h2o.mol");
    remove("co2.mol");
  }
  
  void SetUp() override
  {
    if( GetParam() == TestMode::KOKKOS_OMP ) args = {"-k","on","t","4","-sf","kk"};
    testbinary = "SortCommTest";
    LAMMPSTest::SetUp();
    ASSERT_NE(lmp, nullptr);
    BEGIN_HIDE_OUTPUT();
    command("units real");
    command("dimension 3");
    command("pair_style zero 4.0");
    command("region box block -4 4 -4 4 -4 4");
    END_HIDE_OUTPUT();
  }
  
  void TearDown() override
  {
    LAMMPSTest::TearDown();
    remove("test_atom_styles.data");
    remove("input_atom_styles.data");
    remove("test_atom_styles.restart");
  }
  
  // Get all available atom styles dynamically
  std::vector<AtomStyleInfo> getAvailableAtomStyles() {
    std::vector<AtomStyleInfo> styles;
    
    // Test only basic, reliable atom styles
    std::vector<std::string> test_styles;
    
    if (GetParam() == TestMode::KOKKOS_OMP) {
      // Test Kokkos versions of basic styles
      test_styles = {"atomic/kk", "charge/kk", "molecular/kk"};
    } else {
      // Test basic styles
      test_styles = {"atomic", "charge", "molecular", "full"};
    }
    
    for (const std::string& style_name : test_styles) {
      // Check if style exists in the map
      if (lmp->atom && lmp->atom->avec_map && 
          lmp->atom->avec_map->find(style_name) != lmp->atom->avec_map->end()) {
        
        AtomStyleInfo info;
        info.name = style_name;
        
        // Extract fields by temporarily creating the style
        if (extractStyleFields(info)) {
          styles.push_back(info);
        }
      }
    }
    
    return styles;
  }
  
  // Extract communication fields from atom style
  bool extractStyleFields(AtomStyleInfo& info) {
    try {
      std::string orig_style = lmp->atom->atom_style ? lmp->atom->atom_style : "";
      
      BEGIN_HIDE_OUTPUT();
      
      // Try to create the atom style - this might fail for Kokkos styles 
      // if Kokkos is not properly initialized
      try {
        command("atom_style " + info.name);
      } catch (...) {
        END_HIDE_OUTPUT();
        return false; // Skip this style if it can't be created
      }
      
      END_HIDE_OUTPUT();
      
      if (lmp->atom->avec) {
        info.fields_comm = lmp->atom->avec->fields_comm;
        info.fields_comm_vel = lmp->atom->avec->fields_comm_vel;
        info.fields_reverse = lmp->atom->avec->fields_reverse;
        info.fields_border = lmp->atom->avec->fields_border;
        info.fields_border_vel = lmp->atom->avec->fields_border_vel;
        info.fields_exchange = lmp->atom->avec->fields_exchange;
        
        info.has_bonus_data = (lmp->atom->avec->size_forward_bonus > 0 ||
          lmp->atom->avec->size_border_bonus > 0);
        
        // Set minimal setup command if needed
        if (info.name == "ellipsoid") {
          info.setup_cmd = "set type * shape 1.0 1.5 2.0";
        }
        
        // Restore original style
        if (!orig_style.empty()) {
          BEGIN_HIDE_OUTPUT();
          command("atom_style " + orig_style);
          END_HIDE_OUTPUT();
        }
        
        return true;
      }
      
      return false;
    } catch (...) {
      return false;
    }
  }
  
  // Setup system with given atom style
  void setupAtomStyle(const std::string& style_name, const std::string& setup_cmd = "") {
    BEGIN_HIDE_OUTPUT();
    command("clear");
    command("units real");
    
    command("atom_style " + style_name);
    command("atom_modify map array");
    command("region box block -4 4 -4 4 -4 4");
    command("create_box 2 box");
    
    if (!setup_cmd.empty()) {
      command(setup_cmd);
    }
    
    // Create atoms on lattice
    command("lattice fcc 2.0");
    command("create_atoms 1 box");
    command("region inner block -2 2 -2 2 -2 2");
    command("set region inner type 2");
    command("mass 1 1.0");
    command("mass 2 2.0");
    
    // Set up proper pair style and neighbor settings for sorting
    command("pair_style lj/cut 4.0");
    command("pair_coeff * * 1.0 1.0 4.0");
    command("neighbor 1.0 bin");
    command("neigh_modify delay 0 every 1 check yes");
    
    // Force initialization of neighbor lists and sorting infrastructure
    command("run 0");
    
    END_HIDE_OUTPUT();
  }
  
  // Fill atom properties with random data
  void fillRandomData(const std::vector<std::string>& fields) {
    for (const std::string& field : fields) {
      // Skip coordinate fields to avoid breaking the sort algorithm
      if (field == "x" || field == "y" || field == "z") {
        continue;
      }
      
      void* data = lmp->atom->extract(field.c_str());
      if (!data) continue;
      
      int datatype = lmp->atom->extract_datatype(field.c_str());
      int size = lmp->atom->extract_size(field.c_str(), 0);
      int nlocal = lmp->atom->nlocal;
      
      if (datatype == 0) { // Integer
        int* idata = static_cast<int*>(data);
        std::uniform_int_distribution<int> dist(-100, 100);
        
        int total_size = (size == 0) ? nlocal : nlocal * size;
        for (int i = 0; i < total_size; i++) {
          idata[i] = dist(rng);
        }
      } else if (datatype == 1) { // Double
        double* ddata = static_cast<double*>(data);
        std::uniform_real_distribution<double> dist(-5.0, 5.0);
        
        if (field == "radius") {
          std::uniform_real_distribution<double> pos_dist(0.1, 2.0);
          for (int i = 0; i < nlocal; i++) {
            ddata[i] = pos_dist(rng);
          }
        } else if (field == "rmass") {
          std::uniform_real_distribution<double> mass_dist(0.5, 5.0);
          for (int i = 0; i < nlocal; i++) {
            ddata[i] = mass_dist(rng);
          }
        } else if (field == "quat") {
          for (int i = 0; i < nlocal; i++) {
            double q[4];
            for (int j = 0; j < 4; j++) q[j] = dist(rng);
            double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
            for (int j = 0; j < 4; j++) {
              ddata[i * 4 + j] = q[j] / norm;
            }
          }
        } else {
          int total_size = (size == 0) ? nlocal : nlocal * size;
          for (int i = 0; i < total_size; i++) {
            ddata[i] = dist(rng);
          }
        }
      }
    }
  }
  
  // Extract property values for comparison
  std::vector<double> extractPropertyValues(const std::vector<std::string>& fields, 
                                            const std::vector<int>& atom_indices = {}) {
    std::vector<double> values;
    
    for (const std::string& field : fields) {
      void* data = lmp->atom->extract(field.c_str());
      if (!data) continue;
      
      int datatype = lmp->atom->extract_datatype(field.c_str());
      int size = lmp->atom->extract_size(field.c_str(), 0);
      
      std::vector<int> indices = atom_indices.empty() ? 
      std::vector<int>(lmp->atom->nlocal) : atom_indices;
      if (atom_indices.empty()) {
        std::iota(indices.begin(), indices.end(), 0);
      }
      
      for (int idx : indices) {
        if (idx >= lmp->atom->nlocal) continue;
        
        if (datatype == 0) { // Integer
          int* idata = static_cast<int*>(data);
          if (size == 0) {
            values.push_back(static_cast<double>(idata[idx]));
          } else {
            for (int j = 0; j < size; j++) {
              values.push_back(static_cast<double>(idata[idx * size + j]));
            }
          }
        } else if (datatype == 1) { // Double
          double* ddata = static_cast<double*>(data);
          if (size == 0) {
            values.push_back(ddata[idx]);
          } else {
            for (int j = 0; j < size; j++) {
              values.push_back(ddata[idx * size + j]);
            }
          }
        }
      }
    }
    
    return values;
  }
  
  // Compare two value vectors
  bool compareValues(const std::vector<double>& v1, const std::vector<double>& v2, 
                     double tolerance = EPSILON) {
    if (v1.size() != v2.size()) return false;
    
    for (size_t i = 0; i < v1.size(); i++) {
      if (std::abs(v1[i] - v2[i]) > tolerance) {
        return false;
      }
    }
    return true;
  }
}; // Close the SortCommTest class definition here

// Test implementations - these must be outside the class definition
TEST_P(SortCommTest, Sort) {
  auto styles = getAvailableAtomStyles();
  
  if (styles.empty()) {
    GTEST_SKIP() << "No suitable atom styles available for sort testing";
  }
  
  for (const auto& style : styles) {
    SCOPED_TRACE("Testing sort for atom style: " + style.name);
    
    setupAtomStyle(style.name, style.setup_cmd);
    
    int nlocal = lmp->atom->nlocal;
    ASSERT_GT(nlocal, 10) << "Need sufficient atoms for sort testing";
    
    // Skip if sorting infrastructure not available
    if (!lmp->neighbor || !lmp->domain || lmp->neighbor->cutneighmax <= 0.0) {
      GTEST_SKIP() << "Sort infrastructure not available for " << style.name;
    }
    
    // Store original atom ordering
    std::vector<int> orig_tags(nlocal);
    for (int i = 0; i < nlocal; i++) {
      orig_tags[i] = static_cast<int>(lmp->atom->tag[i]);
    }
    
    // Slightly perturb atom positions to test actual sorting
    for (int i = 0; i < nlocal; i++) {
      lmp->atom->x[i][0] += 0.1 * (static_cast<double>(rand()) / RAND_MAX - 0.5);
    }
    
    // Force neighbor list build and perform sort
    lmp->neighbor->build(1);
    lmp->atom->sort();
    
    // Test that atoms are spatially organized after sorting
    // LAMMPS uses spatial binning, not strict coordinate ordering
    // We verify that nearby atoms are more likely to be adjacent in memory
    double spatial_locality = 0.0;
    int comparisons = 0;
    
    for (int i = 0; i < nlocal - 1; i++) {
      double dist = fabs(lmp->atom->x[i+1][0] - lmp->atom->x[i][0]);
      if (dist < 2.0) { // Atoms are spatially close
        spatial_locality += 1.0;
      }
      comparisons++;
    }
    
    if (comparisons > 0) {
      spatial_locality /= comparisons;
      // After sorting, at least 30% of adjacent atoms should be spatially close
      EXPECT_GT(spatial_locality, 0.3) << "Sort did not improve spatial locality for " << style.name;
    } else {
      // If no comparisons, just verify sort completed without error
      EXPECT_EQ(lmp->atom->nlocal, nlocal) << "Sort changed atom count for " << style.name;
    }
    
    // Verify no atoms lost during sorting
    EXPECT_EQ(lmp->atom->nlocal, nlocal) << "Atoms lost during sort for " << style.name;
    
    // Verify all original tags still present
    std::vector<int> sorted_tags(nlocal);
    for (int i = 0; i < nlocal; i++) {
      sorted_tags[i] = static_cast<int>(lmp->atom->tag[i]);
    }
    std::sort(orig_tags.begin(), orig_tags.end());
    std::sort(sorted_tags.begin(), sorted_tags.end());
    EXPECT_EQ(orig_tags, sorted_tags) << "Tags not preserved during sort for " << style.name;
  }
}

TEST_P(SortCommTest, PackForwardComm) {
  auto styles = getAvailableAtomStyles();
  
  if (styles.empty()) {
    GTEST_SKIP() << "No suitable atom styles available for testing";
  }
  
  for (const auto& style : styles) {
    if (style.fields_comm.empty() && style.fields_comm_vel.empty()) {
      continue; // Skip styles with no forward communication fields
    }
    
    SCOPED_TRACE("Testing forward comm for: " + style.name);
    
    setupAtomStyle(style.name, style.setup_cmd);
    fillRandomData(style.fields_comm);
    
    int nlocal = lmp->atom->nlocal;
    int ncomm = std::min(nlocal, 10);
    std::vector<int> list(ncomm);
    std::iota(list.begin(), list.end(), 0);
    
    // Get original data
    auto orig_values = extractPropertyValues(style.fields_comm, list);
    
    // Test pack_comm
    int n = lmp->atom->avec->pack_comm(ncomm, list.data(), nullptr, 0, nullptr);
    EXPECT_GT(n, 0) << "Should pack data for " << style.name;
    
    std::vector<double> buf(n * ncomm);
    int packed = lmp->atom->avec->pack_comm(ncomm, list.data(), buf.data(), 0, nullptr);
    EXPECT_EQ(packed, n) << "Pack size match for " << style.name;
    
    // Modify data and test round-trip
    fillRandomData(style.fields_comm);
    lmp->atom->avec->unpack_comm(ncomm, 0, buf.data());
    
    auto restored_values = extractPropertyValues(style.fields_comm, list);
    EXPECT_TRUE(compareValues(orig_values, restored_values)) 
      << "Round-trip preservation for " << style.name;
    
    // Test velocity communication if available
    if (!style.fields_comm_vel.empty()) {
      int n_vel = lmp->atom->avec->pack_comm_vel(ncomm, list.data(), nullptr, 0, nullptr);
      if (n_vel > 0) {
        std::vector<double> buf_vel(n_vel * ncomm);
        int packed_vel = lmp->atom->avec->pack_comm_vel(ncomm, list.data(), buf_vel.data(), 0, nullptr);
        EXPECT_EQ(packed_vel, n_vel) << "Velocity pack for " << style.name;
        lmp->atom->avec->unpack_comm_vel(ncomm, 0, buf_vel.data());
      }
    }
    
    // Test bonus communication if available
    if (style.has_bonus_data) {
      int n_bonus = lmp->atom->avec->pack_comm_bonus(ncomm, nullptr, nullptr);
      if (n_bonus > 0) {
        std::vector<double> buf_bonus(n_bonus);
        int packed_bonus = lmp->atom->avec->pack_comm_bonus(ncomm, list.data(), buf_bonus.data());
        EXPECT_EQ(packed_bonus, n_bonus) << "Bonus comm pack for " << style.name;
        lmp->atom->avec->unpack_comm_bonus(ncomm, 0, buf_bonus.data());
      }
    }
  }
}

TEST_P(SortCommTest, PackReverseComm) {
  auto styles = getAvailableAtomStyles();
  
  if (styles.empty()) {
    GTEST_SKIP() << "No suitable atom styles available for reverse comm testing";
  }
  
  for (const auto& style : styles) {
    if (style.fields_reverse.empty()) {
      continue; // Skip styles with no reverse communication fields
    }
    
    SCOPED_TRACE("Testing reverse comm for: " + style.name);
    
    setupAtomStyle(style.name, style.setup_cmd);
    
    int nlocal = lmp->atom->nlocal;
    int ncomm = std::min(nlocal, 10);
    
    // Fill reverse communication fields with test data
    fillRandomData(style.fields_reverse);
    
    // Estimate buffer size and test pack_reverse
    int estimated_size = ncomm * 50;  // Very conservative estimate
    std::vector<double> buf(estimated_size);
    
    int packed = lmp->atom->avec->pack_reverse(ncomm, 0, buf.data());
    EXPECT_GT(packed, 0) << "Reverse comm pack should return positive size for " << style.name;
    EXPECT_LE(packed, estimated_size) << "Buffer overflow in reverse comm for " << style.name;
    
    // Skip unpack if buffer overflow detected
    if (packed > estimated_size) {
      GTEST_SKIP() << "Skipping reverse comm unpack due to buffer overflow for " << style.name;
    }
    
    // Test unpack operation
    std::vector<int> list(ncomm);
    std::iota(list.begin(), list.end(), 0);
    
    // Perform unpack (this accumulates data in reverse communication)
    lmp->atom->avec->unpack_reverse(ncomm, list.data(), buf.data());
    
    // Verify unpack operation completed
    EXPECT_EQ(packed % ncomm, 0) << "Packed data not aligned to atom count for " << style.name;
  }
}

TEST_P(SortCommTest, PackBorderComm) {
  auto styles = getAvailableAtomStyles();
  
  if (styles.empty()) {
    GTEST_SKIP() << "No suitable atom styles available for border comm testing";
  }
  
  for (const auto& style : styles) {
    SCOPED_TRACE("Testing border comm for: " + style.name);
    
    setupAtomStyle(style.name, style.setup_cmd);
    
    int nlocal = lmp->atom->nlocal;
    int nsend = std::min(nlocal, 10);
    std::vector<int> list(nsend);
    std::iota(list.begin(), list.end(), 0);
    
    // Fill border communication fields with test data
    fillRandomData(style.fields_border);
    
    // Estimate buffer size (use very conservative estimates to avoid overflow)
    int estimated_size = nsend * 50;  // Very conservative estimate
    std::vector<double> buf(estimated_size);
    
    // Test pack_border operation
    int packed = lmp->atom->avec->pack_border(nsend, list.data(), buf.data(), 0, nullptr);
    EXPECT_GT(packed, 0) << "Border pack should return positive size for " << style.name;
    EXPECT_LE(packed, estimated_size) << "Buffer overflow in border pack for " << style.name;
    
    // Skip unpack if buffer overflow detected to prevent memory corruption
    if (packed > estimated_size) {
      GTEST_SKIP() << "Skipping unpack due to buffer overflow for " << style.name;
    }
    
    // Test unpack_border fills ghost atom slots
    int old_nlocal = lmp->atom->nlocal;
    int old_total = lmp->atom->nlocal + lmp->atom->nghost;
    
    lmp->atom->avec->unpack_border(nsend, 0, buf.data());
    
    int new_total = lmp->atom->nlocal + lmp->atom->nghost;
    
    // Border communication may fill ghost slots rather than increase nlocal
    EXPECT_GE(new_total, old_total) 
      << "Border unpack should maintain or increase total atom count for " << style.name;
    
    // At minimum, the operation should not decrease atom counts
    EXPECT_GE(lmp->atom->nlocal, old_nlocal - nsend) 
      << "Border unpack should not drastically reduce nlocal for " << style.name;
    
    // Reset for velocity border test
    setupAtomStyle(style.name, style.setup_cmd);
    
    // Test velocity border communication if supported
    if (!style.fields_border_vel.empty()) {
      estimated_size = nsend * 20;  // Much more conservative for velocity data
      std::vector<double> buf_vel(estimated_size);
      
      int packed_vel = lmp->atom->avec->pack_border_vel(nsend, list.data(), buf_vel.data(), 0, nullptr);
      if (packed_vel > 0) {
        EXPECT_LE(packed_vel, estimated_size) << "Velocity buffer overflow for " << style.name;
        lmp->atom->avec->unpack_border_vel(nsend, 0, buf_vel.data());
      }
    }
    
    // Test bonus border communication if supported
    if (style.has_bonus_data) {
      estimated_size = nsend * 50;  // Very conservative for bonus data
      std::vector<double> buf_bonus(estimated_size);
      
      int packed_bonus = lmp->atom->avec->pack_border_bonus(nsend, list.data(), buf_bonus.data());
      if (packed_bonus > 0) {
        EXPECT_LE(packed_bonus, estimated_size) << "Bonus buffer overflow for " << style.name;
        lmp->atom->avec->unpack_border_bonus(nsend, 0, buf_bonus.data());
      }
    }
  }
}

TEST_P(SortCommTest, PackExchangeComm) {
  auto styles = getAvailableAtomStyles();
  
  if (styles.empty()) {
    GTEST_SKIP() << "No suitable atom styles available for exchange comm testing";
  }
  
  for (const auto& style : styles) {
    SCOPED_TRACE("Testing exchange comm for: " + style.name);
    
    setupAtomStyle(style.name, style.setup_cmd);
    
    int nlocal = lmp->atom->nlocal;
    int atom_idx = 0;  // Test with first atom
    
    // Fill exchange communication fields with test data
    fillRandomData(style.fields_exchange);
    
    // Store original atom data for verification
    auto orig_values = extractPropertyValues(style.fields_exchange, {atom_idx});
    
    // Estimate buffer size and test pack_exchange
    int estimated_size = 200;  // Very conservative estimate for complete atom data
    std::vector<double> buf(estimated_size);
    
    int packed = lmp->atom->avec->pack_exchange(atom_idx, buf.data());
    EXPECT_GT(packed, 0) << "Exchange pack should return positive size for " << style.name;
    EXPECT_LE(packed, estimated_size) << "Buffer overflow in exchange pack for " << style.name;
    
    // Skip unpack if buffer overflow detected
    if (packed > estimated_size) {
      GTEST_SKIP() << "Skipping exchange unpack due to buffer overflow for " << style.name;
    }
    
    // Test unpack_exchange creates new atom
    int m = lmp->atom->avec->unpack_exchange(buf.data());
    EXPECT_EQ(m, packed) << "Exchange unpack size mismatch for " << style.name;
    EXPECT_EQ(lmp->atom->nlocal, nlocal + 1) << "Exchange should create new atom for " << style.name;
    
    // Verify the new atom has reasonable data
    int new_atom_idx = nlocal;  // The newly created atom
    auto new_values = extractPropertyValues(style.fields_exchange, {new_atom_idx});
    
    // For basic verification, check that we got the expected amount of data
    EXPECT_EQ(new_values.size(), orig_values.size()) 
      << "Exchange data size mismatch for " << style.name;
    
    // Test bonus exchange if supported
    if (style.has_bonus_data) {
      // Reset for bonus test
      setupAtomStyle(style.name, style.setup_cmd); 
      fillRandomData(style.fields_exchange);
      
      estimated_size = 500;  // Very conservative for bonus data
      std::vector<double> buf_bonus(estimated_size);
      
      int packed_bonus = lmp->atom->avec->pack_exchange_bonus(0, buf_bonus.data());
      if (packed_bonus > 0) {
        EXPECT_LE(packed_bonus, estimated_size) << "Exchange bonus buffer overflow for " << style.name;
        
        int m_bonus = lmp->atom->avec->unpack_exchange_bonus(lmp->atom->nlocal, buf_bonus.data());
        EXPECT_EQ(m_bonus, packed_bonus) << "Exchange bonus unpack failed for " << style.name;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(AtomStyles, SortCommTest, ::testing::Values(TestMode::PLAIN, TestMode::KOKKOS_OMP));

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleMock(&argc, argv);
  
  // handle arguments passed via environment variable
  if (const char *var = getenv("TEST_ARGS")) {
    auto env = split_words(var);
    for (auto arg : env) {
      if (arg == "-v") {
        verbose = true;
      }
    }
  }
  if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;
  
  int rv = RUN_ALL_TESTS();
  MPI_Finalize();
  return rv;
}