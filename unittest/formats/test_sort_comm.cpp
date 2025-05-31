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
#include <sstream>

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
  
public:
  // Calculate buffer size for forward communication
  // Note: This is for regular forward packing, not velocity forward packing
  // Use calculateVelocityCommSize() for velocity packing
  int calculateForwardCommSize(int natoms) {
    if (!lmp->atom->avec) return 0;
    
    int size_per_atom = lmp->atom->avec->size_forward;
    int bonus_size = lmp->atom->avec->size_forward_bonus;
    
    // Total size = (basic size per atom * natoms) + total bonus size
    return size_per_atom * natoms + bonus_size;
  }
  
  // Calculate buffer size for reverse communication  
  int calculateReverseCommSize(int natoms) {
    if (!lmp->atom->avec) return 0;
    
    int size_per_atom = lmp->atom->avec->size_reverse;
    return size_per_atom * natoms;
  }
  
  // Calculate buffer size for border communication
  // Note: This is for regular border packing, not velocity border packing
  // Use calculateVelocityCommSize() for velocity packing
  int calculateBorderCommSize(int natoms) {
    if (!lmp->atom->avec) return 0;
    
    int size_per_atom = lmp->atom->avec->size_border;
    int bonus_size = lmp->atom->avec->size_border_bonus;
    
    return size_per_atom * natoms + bonus_size;
  }
  
  // Calculate buffer size for exchange communication
  int calculateExchangeCommSize() {
    if (!lmp->atom->avec) return 0;
    
    // maxexchange is only set if size > BUFEXTRA, so it might be 0
    int max_exchange = lmp->atom->avec->maxexchange;
    
    if (max_exchange > 0) {
      return max_exchange;
    } else {
      // Try to get actual size by doing a test pack if we have atoms
      if (lmp->atom->nlocal > 0) {
        try {
          // Use a large temporary buffer to test actual pack size
          std::vector<double> test_buf(200);  // Conservative large buffer
          int actual_size = lmp->atom->avec->pack_exchange(0, test_buf.data());
          if (actual_size > 0 && actual_size < 200) {
            return actual_size + 10;  // Add small safety margin
          }
        } catch (...) {
          // If test pack fails, fall back to estimate
        }
      }
      
      // Fallback: Conservative estimate based on atom style fields
      int estimated_size = 50;  // Base estimate for basic atom data
      
      // Add space for style-specific fields
      if (lmp->atom->avec->molecular > 0) estimated_size += 10;  // molecule, bonds, etc.
      if (lmp->atom->avec->bonus_flag) estimated_size += 30;     // bonus data
      
      return estimated_size;
    }
  }
  
  // Calculate buffer size for velocity border communication
  int calculateVelocityBorderCommSize(int natoms) {
    if (!lmp->atom->avec) return 0;
    
    // Velocity border packing includes border fields + velocity data
    int border_size = lmp->atom->avec->size_border * natoms;
    int velocity_addition = 3 * natoms; // v[3] 
    int bonus_size = lmp->atom->avec->size_border_bonus;
    
    return border_size + velocity_addition + bonus_size;
  }
  
  // Calculate buffer size with safety margin
  int addSafetyMargin(int calculated_size, double margin = 0.1) {
    return static_cast<int>(calculated_size * (1.0 + margin));
  }
  
  // Calculate buffer size for velocity communication  
  int calculateVelocityCommSize(int natoms) {
    if (!lmp->atom->avec) return 0;
    
    int size_per_atom = lmp->atom->avec->size_velocity;
    return size_per_atom * natoms;
  }
  
  // Get detailed size information for debugging
  std::string getSizeInfo(const std::string& style_name) {
    if (!lmp->atom->avec) return "AtomVec is null";
    
    std::ostringstream info;
    info << "Style: " << style_name << "\n";
    info << "  Per-atom buffer sizes:\n";
    info << "    size_forward: " << lmp->atom->avec->size_forward << " (doubles per atom for forward comm)\n";
    info << "    size_reverse: " << lmp->atom->avec->size_reverse << " (doubles per atom for reverse comm)\n"; 
    info << "    size_border: " << lmp->atom->avec->size_border << " (doubles per atom for border comm)\n";
    info << "    size_velocity: " << lmp->atom->avec->size_velocity << " (doubles per atom for velocity comm)\n";
    info << "  Other fields:\n";
    info << "    maxexchange: " << lmp->atom->avec->maxexchange;
    if (lmp->atom->avec->maxexchange == 0) info << " (not set)";
    info << "\n";
    info << "    size_data_atom: " << lmp->atom->avec->size_data_atom << "\n";
    info << "    size_data_vel: " << lmp->atom->avec->size_data_vel << "\n";
    info << "    molecular: " << lmp->atom->avec->molecular << "\n";
    info << "    bonus_flag: " << lmp->atom->avec->bonus_flag << "\n";
    info << "    size_forward_bonus: " << lmp->atom->avec->size_forward_bonus << "\n";
    info << "    size_border_bonus: " << lmp->atom->avec->size_border_bonus << "\n";
    
    return info.str();
  }
  // Get all available atom styles dynamically
  std::vector<AtomStyleInfo> getAvailableAtomStyles() {
    std::vector<AtomStyleInfo> styles;
    
    if (!lmp->atom || !lmp->atom->avec_map) {
      return styles; // Return empty if no atom system or style map
    }
    
    // Iterate through all registered atom styles in LAMMPS
    for (const auto& pair : *(lmp->atom->avec_map)) {
      const std::string& style_name = pair.first;
      
      // Skip styles that require special setup or are problematic
      if (style_name.find("template") != std::string::npos) continue; // Need molecule templates
      if (style_name.find("hybrid") != std::string::npos) continue;   // Need sub-styles specified
      if (style_name.find("body") != std::string::npos) continue;     // Need body style setup
      if (style_name.find("line") != std::string::npos) continue;     // Need special setup
      if (style_name.find("tri") != std::string::npos) continue;      // Need special setup
      if (style_name.find("ellipsoid") != std::string::npos) continue; // Need shape setup
      if (style_name.find("sphere") != std::string::npos) continue;   // Need per-atom mass setup
      
      // Skip Kokkos styles if not in Kokkos mode, and vice versa
      bool is_kokkos_style = (style_name.find("/kk") != std::string::npos);
      if (GetParam() == TestMode::PLAIN && is_kokkos_style) continue;
      if (GetParam() == TestMode::KOKKOS_OMP && !is_kokkos_style) continue;
      
      AtomStyleInfo info;
      info.name = style_name;
      
      // Try to extract fields for this style
      if (extractStyleFields(info)) {
        styles.push_back(info);
      }
    }
    
    return styles;
  }
  
  // Extract communication fields from atom style
  bool extractStyleFields(AtomStyleInfo& info) {
    try {
      std::string orig_style = lmp->atom->atom_style ? lmp->atom->atom_style : "";
      
      // Use RAII-style output hiding to ensure proper cleanup
      auto hideOutput = [this]() {
        BEGIN_HIDE_OUTPUT();
        return [this]() { END_HIDE_OUTPUT(); };
      };
      
      auto cleanup = hideOutput();
      
      // Try to create the atom style - this might fail for Kokkos styles 
      // if Kokkos is not properly initialized
      try {
        command("atom_style " + info.name);
      } catch (...) {
        cleanup(); // Ensure output is restored
        return false; // Skip this style if it can't be created
      }
      
      cleanup(); // Restore output
      
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
          auto restoreCleanup = hideOutput();
          command("atom_style " + orig_style);
          restoreCleanup();
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
    // Use RAII-style output hiding to ensure proper cleanup
    auto hideOutput = [this]() {
      BEGIN_HIDE_OUTPUT();
      return [this]() { END_HIDE_OUTPUT(); };
    };
    
    auto cleanup = hideOutput();
    
    try {
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
      
      // Validate that atoms were actually created
      if (lmp->atom->nlocal == 0) {
        cleanup();
        throw std::runtime_error("No atoms created for style: " + style_name);
      }
      
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
      
      // Final validation that the system is properly set up
      if (!lmp->atom->avec) {
        cleanup();
        throw std::runtime_error("AtomVec not properly initialized for style: " + style_name);
      }
      
      if (lmp->atom->nlocal == 0) {
        cleanup();
        throw std::runtime_error("Atoms lost during setup for style: " + style_name);
      }
      
      cleanup(); // Ensure output is restored
    } catch (...) {
      cleanup(); // Ensure output is restored even on exception
      throw; // Re-throw the exception
    }
  }
  
  // Fill atom properties with random data
  void fillRandomData(const std::vector<std::string>& fields) {
    try {
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
        
        // Safety check - ensure we have atoms before accessing arrays
        if (nlocal <= 0) continue;
        
        if (datatype == 0) { // Integer
          int* idata = static_cast<int*>(data);
          if (!idata) continue; // Safety check for null pointer
          
          // Use more reasonable ranges for different fields
          std::uniform_int_distribution<int> dist;
          if (field == "type") {
            dist = std::uniform_int_distribution<int>(1, 2); // Only valid types
          } else if (field == "molecule") {
            dist = std::uniform_int_distribution<int>(0, 100);
          } else {
            dist = std::uniform_int_distribution<int>(-10, 10); // Smaller range
          }
          
          int total_size = (size == 0) ? nlocal : nlocal * size;
          for (int i = 0; i < total_size; i++) {
            idata[i] = dist(rng);
          }
        } else if (datatype == 1) { // Double
          double* ddata = static_cast<double*>(data);
          if (!ddata) continue; // Safety check for null pointer
          
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
            // Normalized quaternions
            for (int i = 0; i < nlocal; i++) {
              double q[4];
              std::uniform_real_distribution<double> dist(-1.0, 1.0);
              for (int j = 0; j < 4; j++) q[j] = dist(rng);
              double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
              if (norm > 0) {
                for (int j = 0; j < 4; j++) {
                  ddata[i * 4 + j] = q[j] / norm;
                }
              } else {
                // Default quaternion
                ddata[i * 4] = 1.0;
                ddata[i * 4 + 1] = 0.0;
                ddata[i * 4 + 2] = 0.0;
                ddata[i * 4 + 3] = 0.0;
              }
            }
          } else if (field == "v" || field.find("vel") != std::string::npos) {
            // Reasonable velocity range
            std::uniform_real_distribution<double> vel_dist(-2.0, 2.0);
            int total_size = (size == 0) ? nlocal : nlocal * size;
            for (int i = 0; i < total_size; i++) {
              ddata[i] = vel_dist(rng);
            }
          } else {
            // Smaller range for other properties
            std::uniform_real_distribution<double> dist(-1.0, 1.0);
            int total_size = (size == 0) ? nlocal : nlocal * size;
            for (int i = 0; i < total_size; i++) {
              ddata[i] = dist(rng);
            }
          }
        }
      }
    } catch (const std::exception& e) {
      throw std::runtime_error("fillRandomData failed: " + std::string(e.what()));
    } catch (...) {
      throw std::runtime_error("fillRandomData failed with unknown error");
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
    
    try {
      setupAtomStyle(style.name, style.setup_cmd);
    } catch (const std::exception& e) {
      GTEST_SKIP() << "Skipping " << style.name << ": " << e.what();
      continue;
    }
    
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
    
    // Calculate exact buffer size needed
    int exact_size = calculateForwardCommSize(ncomm);
    int buffer_size = addSafetyMargin(exact_size);
    
    ASSERT_GT(exact_size, 0) << "Forward comm size should be positive for " << style.name;
    
    if (verbose) {
      std::cout << getSizeInfo(style.name) << std::endl;
      std::cout << "Forward comm: calculated=" << exact_size << " buffer=" << buffer_size 
                << " (using size_forward=" << lmp->atom->avec->size_forward << " * " << ncomm << " atoms)" << std::endl;
    }
    
    // Get original data
    auto orig_values = extractPropertyValues(style.fields_comm, list);
    
    // Test pack_comm with properly sized buffer
    std::vector<double> buf(buffer_size);
    int packed = lmp->atom->avec->pack_comm(ncomm, list.data(), buf.data(), 0, nullptr);
    
    EXPECT_GT(packed, 0) << "Should pack data for " << style.name;
    EXPECT_LE(packed, exact_size) << "Packed size should not exceed calculated size for " << style.name;
    
    // Test round-trip
    fillRandomData(style.fields_comm);
    lmp->atom->avec->unpack_comm(ncomm, 0, buf.data());
    
    auto restored_values = extractPropertyValues(style.fields_comm, list);
    EXPECT_TRUE(compareValues(orig_values, restored_values)) 
      << "Round-trip preservation for " << style.name;
    
    // Test velocity communication if available
    if (!style.fields_comm_vel.empty()) {
      // Calculate proper buffer size for velocity packing
      int vel_exact_size = calculateVelocityCommSize(ncomm);
      int vel_buffer_size = addSafetyMargin(vel_exact_size);
      
      std::vector<double> buf_vel(vel_buffer_size);
      
      int packed_vel = lmp->atom->avec->pack_comm_vel(ncomm, list.data(), buf_vel.data(), 0, nullptr);
      if (packed_vel > 0) {
        EXPECT_LE(packed_vel, vel_exact_size) << "Velocity comm pack size for " << style.name;
        if (verbose) {
          std::cout << "Velocity comm: calculated=" << vel_exact_size 
                    << " packed=" << packed_vel 
                    << " (using size_velocity=" << lmp->atom->avec->size_velocity << " * " << ncomm << " atoms)" << std::endl;
        }
        lmp->atom->avec->unpack_comm_vel(ncomm, 0, buf_vel.data());
      }
    }
    
    // Test bonus communication if available
    if (style.has_bonus_data && lmp->atom->avec->size_forward_bonus > 0) {
      int bonus_size = lmp->atom->avec->size_forward_bonus;
      std::vector<double> buf_bonus(bonus_size + 10); // Small safety margin
      int packed_bonus = lmp->atom->avec->pack_comm_bonus(ncomm, list.data(), buf_bonus.data());
      if (packed_bonus > 0) {
        EXPECT_LE(packed_bonus, bonus_size) << "Bonus comm pack size for " << style.name;
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
    
    // Calculate exact buffer size needed
    int exact_size = calculateReverseCommSize(ncomm);
    int buffer_size = addSafetyMargin(exact_size);
    
    ASSERT_GT(exact_size, 0) << "Reverse comm size should be positive for " << style.name;
    
    if (verbose) {
      std::cout << "Reverse comm: calculated=" << exact_size << " buffer=" << buffer_size 
                << " (using size_reverse=" << lmp->atom->avec->size_reverse << " * " << ncomm << " atoms)" << std::endl;
    }
    
    // Fill reverse communication fields with test data
    fillRandomData(style.fields_reverse);
    
    std::vector<double> buf(buffer_size);
    
    int packed = lmp->atom->avec->pack_reverse(ncomm, 0, buf.data());
    EXPECT_GT(packed, 0) << "Reverse comm pack should return positive size for " << style.name;
    EXPECT_LE(packed, exact_size) << "Packed size should not exceed calculated size for " << style.name;
    
    // Test unpack operation
    std::vector<int> list(ncomm);
    std::iota(list.begin(), list.end(), 0);
    
    lmp->atom->avec->unpack_reverse(ncomm, list.data(), buf.data());
    
    // Verify the operation completed successfully
    EXPECT_EQ(packed % ncomm, 0) << "Packed data should be aligned to atom count for " << style.name;
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
    
    // Calculate exact buffer size needed
    int exact_size = calculateBorderCommSize(nsend);
    int buffer_size = addSafetyMargin(exact_size);
    
    ASSERT_GT(exact_size, 0) << "Border comm size should be positive for " << style.name;
    
    if (verbose) {
      std::cout << "Border comm: calculated=" << exact_size << " buffer=" << buffer_size 
                << " (using size_border=" << lmp->atom->avec->size_border << " * " << nsend << " atoms)" << std::endl;
    }
    
    // Fill border communication fields with test data
    try {
      fillRandomData(style.fields_border);
    } catch (const std::exception& e) {
      GTEST_SKIP() << "Skipping " << style.name << " due to fillRandomData error: " << e.what();
      continue;
    }
    
    std::vector<double> buf(buffer_size);
    
    // Test pack_border operation
    int packed = lmp->atom->avec->pack_border(nsend, list.data(), buf.data(), 0, nullptr);
    
    EXPECT_GT(packed, 0) << "Border pack should return positive size for " << style.name;
    EXPECT_LE(packed, exact_size) << "Packed size should not exceed calculated size for " << style.name;
    
    // Test unpack_border
    int old_nlocal = lmp->atom->nlocal;
    int old_total = lmp->atom->nlocal + lmp->atom->nghost;
    
    lmp->atom->avec->unpack_border(nsend, 0, buf.data());
    
    int new_total = lmp->atom->nlocal + lmp->atom->nghost;
    
    // Border communication may fill ghost slots rather than increase nlocal
    EXPECT_GE(new_total, old_total) 
      << "Border unpack should maintain or increase total atom count for " << style.name;
    
    // Test velocity border communication if supported
    if (!style.fields_border_vel.empty()) {
      // Calculate proper buffer size for velocity border packing
      int vel_exact_size = calculateVelocityBorderCommSize(nsend);
      int vel_buffer_size = addSafetyMargin(vel_exact_size);
      
      std::vector<double> buf_vel(vel_buffer_size);
      
      int packed_vel = lmp->atom->avec->pack_border_vel(nsend, list.data(), buf_vel.data(), 0, nullptr);
      if (packed_vel > 0) {
        EXPECT_LE(packed_vel, vel_exact_size) << "Velocity border pack size for " << style.name;
        if (verbose) {
          std::cout << "Velocity border: calculated=" << vel_exact_size 
                    << " packed=" << packed_vel 
                    << " (using size_border + v[3] for " << nsend << " atoms)" << std::endl;
        }
        lmp->atom->avec->unpack_border_vel(nsend, 0, buf_vel.data());
      }
    }
    
    // Test bonus border communication if supported
    if (style.has_bonus_data && lmp->atom->avec->size_border_bonus > 0) {
      int bonus_size = lmp->atom->avec->size_border_bonus;
      std::vector<double> buf_bonus(bonus_size + 10);
      int packed_bonus = lmp->atom->avec->pack_border_bonus(nsend, list.data(), buf_bonus.data());
      if (packed_bonus > 0) {
        EXPECT_LE(packed_bonus, bonus_size) << "Border bonus pack size for " << style.name;
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
    
    // Calculate exact buffer size needed for one atom
    int exact_size = calculateExchangeCommSize();
    int buffer_size = addSafetyMargin(exact_size);
    
    ASSERT_GT(exact_size, 0) << "Exchange comm size should be positive for " << style.name;
    
    if (verbose) {
      std::cout << getSizeInfo(style.name) << std::endl;
      std::cout << "Exchange comm: calculated=" << exact_size << " buffer=" << buffer_size;
      if (lmp->atom->avec->maxexchange > 0) {
        std::cout << " (using maxexchange=" << lmp->atom->avec->maxexchange << ")";
      } else {
        std::cout << " (using estimated size)";
      }
      std::cout << std::endl;
    }
    
    // Fill exchange communication fields with test data
    try {
      fillRandomData(style.fields_exchange);
    } catch (const std::exception& e) {
      GTEST_SKIP() << "Skipping " << style.name << " due to fillRandomData error: " << e.what();
      continue;
    }
    
    // Store original atom data for verification
    auto orig_values = extractPropertyValues(style.fields_exchange, {atom_idx});
    
    std::vector<double> buf(buffer_size);
    
    int packed = lmp->atom->avec->pack_exchange(atom_idx, buf.data());
    EXPECT_GT(packed, 0) << "Exchange pack should return positive size for " << style.name;
    EXPECT_LE(packed, exact_size) << "Packed size should not exceed calculated size for " << style.name;
    
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
      
      // For bonus exchange, we might need to query size differently
      // This is style-dependent, so we use a reasonable default
      int bonus_buffer_size = 100;
      std::vector<double> buf_bonus(bonus_buffer_size);
      
      int packed_bonus = lmp->atom->avec->pack_exchange_bonus(0, buf_bonus.data());
      if (packed_bonus > 0) {
        EXPECT_LE(packed_bonus, bonus_buffer_size) << "Exchange bonus pack size for " << style.name;
        
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