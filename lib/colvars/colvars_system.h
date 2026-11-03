// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARS_SYSTEM_H
#define COLVARS_SYSTEM_H

#include "colvartypes.h"

/// Class to store the system's boundary conditions
class colvarmodule::system_boundary_conditions {
public:

  /// Type of boundary conditions defined for the current computation
  enum class types {
    non_periodic,   /// All three dimensions are non-periodic
    mixed,          /// Some dimensions are periodic, but others are not
    pbc_orthogonal, /// All three dimensions are periodic, lattice vectors are orthogonal
    pbc_triclinic,  /// All three dimensions are periodic, lattice vectors are *not* orthogonal
    unsupported     /// Unsupported boundary conditions
  };

  /// Default constructor
  inline COLVARS_HOST_DEVICE system_boundary_conditions() { reset(); }

  /// Copy constructor
  system_boundary_conditions(system_boundary_conditions const &) = default;

  /// Type of boundary conditions in the current computation
  inline COLVARS_HOST_DEVICE types type() const { return type_; }

  /// Set the type of boundary explicitly
  inline COLVARS_HOST_DEVICE void set_type(types t) { type_ = t; }

  /// Compute the distance between two positions
  inline COLVARS_HOST_DEVICE cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                                           cvm::atom_pos const &pos2) const;

  /// Compute a shift vector that accounts for tilt factors up to 0.5
  inline COLVARS_HOST_DEVICE cvm::rvector get_triclinic_shift(cvm::rvector const &diff) const;

  /// Reset to defaults (non-periodic)
  inline COLVARS_HOST_DEVICE void reset() {
    periodic_x = periodic_y = periodic_z = false;
    type_ = types::non_periodic;
    unit_cell_x.reset();
    unit_cell_y.reset();
    unit_cell_z.reset();
    reciprocal_cell_x.reset();
    reciprocal_cell_y.reset();
    reciprocal_cell_z.reset();
  }

  /// Set from explicit boundary configuration
  inline COLVARS_HOST_DEVICE void set_boundaries(bool periodic_x_in, bool periodic_y_in,
                                                bool periodic_z_in, cvm::rvector const &A,
                                                cvm::rvector const &B, cvm::rvector const &C);
protected:

  /// Type of boundary conditions in the current computation
  types type_ = types::non_periodic;

  /// Bravais lattice vectors
  cvm::rvector unit_cell_x, unit_cell_y, unit_cell_z;

  /// Reciprocal lattice vectors
  cvm::rvector reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z;

  /// Periodic flags in each dimension
  bool periodic_x = false, periodic_y = false, periodic_z = false;
};


/// Set from explicit boundary configuration
inline COLVARS_HOST_DEVICE void
cvm::system_boundary_conditions::set_boundaries(bool periodic_x_in, bool periodic_y_in,
                                                bool periodic_z_in, cvm::rvector const &A,
                                                cvm::rvector const &B, cvm::rvector const &C)
{
  constexpr double diagonal_tol2 = 1.0e-10;

  periodic_x = periodic_x_in;
  periodic_y = periodic_y_in;
  periodic_z = periodic_z_in;

  // Avoid using stale reciprocal vectors when switching boundary types (especially for mixed periodicity)
  reciprocal_cell_x.reset();
  reciprocal_cell_y.reset();
  reciprocal_cell_z.reset();
  if ((periodic_x == periodic_y) && (periodic_x == periodic_z)) {
    if (periodic_x) {
      // Temporarily set as fully-periodic & orthogonal; will check below for triclinic
      set_type(types::pbc_orthogonal);
    } else {
      reset();
      return;
    }
  } else {
    set_type(types::mixed);
  }

  bool off_diagonal = false;

  if (periodic_x) {
    unit_cell_x = A;
    if ((A.y * A.y) > diagonal_tol2 || (A.z * A.z) > diagonal_tol2) {
      off_diagonal = true;
    } else {
      reciprocal_cell_x = unit_cell_x/unit_cell_x.norm2();
    }
  } else {
    unit_cell_x.reset();
    reciprocal_cell_x.reset();
  }

  if (periodic_y) {
    unit_cell_y = B;
    if ((B.x * B.x) > diagonal_tol2 || (B.z * B.z) > diagonal_tol2) {
      off_diagonal = true;
    } else {
      reciprocal_cell_y = unit_cell_y/unit_cell_y.norm2();
    }
  } else {
    unit_cell_y.reset();
    reciprocal_cell_y.reset();
  }

  if (periodic_z) {
    unit_cell_z = C;
    if ((C.x * C.x) > diagonal_tol2 || (C.y * C.y) > diagonal_tol2) {
      off_diagonal = true;
    } else {
      reciprocal_cell_z = unit_cell_z/unit_cell_z.norm2();
    }
  } else {
    unit_cell_z.reset();
    reciprocal_cell_z.reset();
  }

  if (type() == types::pbc_orthogonal && off_diagonal) {
    set_type(types::pbc_triclinic);
  }

  if (type() == types::pbc_triclinic) {
    cvm::rvector const v_yz = cvm::rvector::outer(unit_cell_y, unit_cell_z);
    reciprocal_cell_x = v_yz / (v_yz * unit_cell_x);
    cvm::rvector const v_zx = cvm::rvector::outer(unit_cell_z, unit_cell_x);
    reciprocal_cell_y = v_zx / (v_zx * unit_cell_y);
    cvm::rvector const v_xy = cvm::rvector::outer(unit_cell_x, unit_cell_y);
    reciprocal_cell_z = v_xy / (v_xy * unit_cell_z);
  }
}


inline COLVARS_HOST_DEVICE
cvm::rvector cvm::system_boundary_conditions::position_distance(cvm::atom_pos const &pos1,
                                                                cvm::atom_pos const &pos2) const
{
  cvm::rvector diff = (pos2 - pos1);

  if (type() == types::non_periodic) {
    return diff;
  }

  if (type() == types::unsupported) {
#if !(defined(__HIP_DEVICE_COMPILE__)) && !(defined(__CUDA_ARCH__))
    cvm::error_static("Error: unsupported boundary conditions.\n", COLVARS_INPUT_ERROR);
#endif
    return diff;
  }

  cvm::real const x_shift = ::floor(reciprocal_cell_x * diff + 0.5);
  cvm::real const y_shift = ::floor(reciprocal_cell_y * diff + 0.5);
  cvm::real const z_shift = ::floor(reciprocal_cell_z * diff + 0.5);

  diff.x -= x_shift * unit_cell_x.x + y_shift * unit_cell_y.x + z_shift * unit_cell_z.x;
  diff.y -= x_shift * unit_cell_x.y + y_shift * unit_cell_y.y + z_shift * unit_cell_z.y;
  diff.z -= x_shift * unit_cell_x.z + y_shift * unit_cell_y.z + z_shift * unit_cell_z.z;

  if (type() != types::pbc_orthogonal) {
    // Matches both "mixed" and "pbc_triclinic", because reciprocal cell vectors are not used
    diff += get_triclinic_shift(diff);
  }

  return diff;
}


inline COLVARS_HOST_DEVICE cvm::rvector
cvm::system_boundary_conditions::get_triclinic_shift(cvm::rvector const &diff) const
{
  cvm::real min_dist2 = diff.norm2();
  cvm::rvector result{0.0, 0.0, 0.0};

  int const nx = periodic_x ? 1 : 0;
  int const ny = periodic_y ? 1 : 0;
  int const nz = periodic_z ? 1 : 0;

  // Loop over neighboring cells to find a shorter distance
  for (int ix = -nx; ix <= nx; ix++) {
    for (int iy = -ny; iy <= ny; iy++) {
      for (int iz = -nz; iz <= nz; iz++) {
        cvm::rvector const shift = ix * unit_cell_x + iy * unit_cell_y + iz * unit_cell_z;
        cvm::real const this_dist2 = (diff + shift).norm2();
        if (this_dist2 < min_dist2) {
          result = shift;
          min_dist2 = this_dist2;
        }
      }
    }
  }

  return result;
}


#endif
