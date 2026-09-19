/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

//
// Contributing author, Richard Meng, Queen's University at Kingston, 22.11.24, contact@richardzjm.com
//

#include "mtp_radial_basis.h"

#include "error.h"
#include "memory.h"
#include "text_file_reader.h"

using namespace LAMMPS_NS;

RadialMTPBasis::RadialMTPBasis(TextFileReader &tfr, LAMMPS *lmp) :
    Pointers(lmp), size(0), min_cutoff(0.0), max_cutoff(0.0), scaling(1.0),
    radial_basis_vals(nullptr), radial_basis_ders(nullptr)
{ read_basis_properties(tfr); }

RadialMTPBasis::RadialMTPBasis(int size, LAMMPS *lmp) :
    Pointers(lmp), size(size), min_cutoff(0.0), max_cutoff(0.0), scaling(1.0),
    radial_basis_vals(nullptr), radial_basis_ders(nullptr)
{
  if (size < 2) error->one(FLERR, "MTP radial basis size must be at least 2.");
  memory->create(radial_basis_vals, size, "pair:mtp_radial_vals");
  memory->create(radial_basis_ders, size, "pair:mtp_radial_ders");
}

void RadialMTPBasis::read_basis_properties(TextFileReader &tfr)
{
  const std::string new_separators = "=, ";
  const std::string separators = TOKENIZER_DEFAULT_SEPARATORS + new_separators;

  ValueTokenizer line_tokens{std::string(tfr.next_line()), separators};
  std::string keyword = line_tokens.next_string();

  // First check if scaling is available
  if (keyword == "scaling") {
    scaling = line_tokens.next_double();
    line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
    keyword = line_tokens.next_string();
  }

  // Read the lower cutoff
  if (keyword != "min_val" && keyword != "min_dist")
    error->one(FLERR, "Error in reading MTP file. Cannot read lower cutoff.");
  min_cutoff = line_tokens.next_double();

  // Read the upper cutoff
  line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
  keyword = line_tokens.next_string();
  if (keyword != "max_val" && keyword != "max_dist")
    error->one(FLERR, "Error in reading MTP file. Cannot read upper cutoff.");
  max_cutoff = line_tokens.next_double();
  if (!(max_cutoff > min_cutoff))
    error->one(FLERR, "MTP maximum cutoff must exceed the minimum cutoff.");

  // Read the basis size set value
  line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
  keyword = line_tokens.next_string();
  if (keyword != "radial_basis_size")
    error->one(FLERR, "Error in reading MTP file. Cannot read radial basis set size.");
  size = line_tokens.next_int();
  if (size < 2) error->one(FLERR, "MTP radial basis size must be at least 2.");

  //Allocate the memory for the basis set values and derivatives.
  memory->create(radial_basis_vals, size, "pair:mtp_radial_vals");
  memory->create(radial_basis_ders, size, "pair:mtp_radial_ders");
}

RadialMTPBasis::~RadialMTPBasis()
{
  memory->destroy(radial_basis_vals);
  memory->destroy(radial_basis_ders);
}

void RBChebyshev::calc_radial_basis_ders(double dist)
{
  const double delta = dist - max_cutoff;
  const double span = max_cutoff - min_cutoff;
  const double mult = 2.0 / span;
  const double ksi = (2 * dist - (min_cutoff + max_cutoff)) / span;

  double val0 = scaling * delta * delta;
  double val1 = scaling * (ksi * delta * delta);
  double der0 = scaling * 2 * delta;
  double der1 = scaling * (mult * delta * delta + 2 * ksi * delta);
  radial_basis_vals[0] = val0;
  radial_basis_vals[1] = val1;
  radial_basis_ders[0] = der0;
  radial_basis_ders[1] = der1;

  for (int i = 2; i < size; i++) {
    const double val = 2 * ksi * val1 - val0;
    const double der = 2 * (mult * val1 + ksi * der1) - der0;
    radial_basis_vals[i] = val;
    radial_basis_ders[i] = der;
    val0 = val1;
    val1 = val;
    der0 = der1;
    der1 = der;
  }
}
