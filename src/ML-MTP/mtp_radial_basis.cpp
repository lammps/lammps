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
#include "utils.h"

using namespace LAMMPS_NS;

RadialMTPBasis::RadialMTPBasis(TextFileReader &tfr, LAMMPS *lmp) : Pointers(lmp)
{ read_basis_properties(tfr); }

RadialMTPBasis::RadialMTPBasis(int size, LAMMPS *lmp) : Pointers(lmp), size(size)
{
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
    error->all(FLERR, "Error in reading MTP file. Cannot read lower cutoff.");
  min_cutoff = line_tokens.next_double();

  // Read the upper cutoff
  line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
  keyword = line_tokens.next_string();
  if (keyword != "max_val" && keyword != "max_dist")
    error->all(FLERR, "Error in reading MTP file. Cannot read upper cutoff.");
  max_cutoff = line_tokens.next_double();

  // Read the basis size set value
  line_tokens = ValueTokenizer(std::string(tfr.next_line()), separators);
  keyword = line_tokens.next_string();
  if (keyword != "radial_basis_size")
    error->all(FLERR, "Error in reading MTP file. Cannot read radial basis set size.");
  size = line_tokens.next_int();

  //Allocate the memory for the basis set values and derivatives.
  memory->create(radial_basis_vals, size, "pair:mtp_radial_vals");
  memory->create(radial_basis_ders, size, "pair:mtp_radial_ders");
}

RadialMTPBasis::~RadialMTPBasis()
{
  memory->destroy(radial_basis_vals);
  memory->destroy(radial_basis_ders);
}

void RBChebyshev::calc_radial_basis(double dist)
{
  double ksi = (2 * dist - (min_cutoff + max_cutoff)) / (max_cutoff - min_cutoff);

  radial_basis_vals[0] = scaling * (dist - max_cutoff) * (dist - max_cutoff);
  radial_basis_vals[1] = scaling * (ksi * (dist - max_cutoff) * (dist - max_cutoff));
  for (int i = 2; i < size; i++) {
    radial_basis_vals[i] = 2 * ksi * radial_basis_vals[i - 1] - radial_basis_vals[i - 2];
  }
}

void RBChebyshev::calc_radial_basis_ders(double dist)
{
  RBChebyshev::calc_radial_basis(dist);

  double mult = 2.0 / (max_cutoff - min_cutoff);
  double ksi = (2 * dist - (min_cutoff + max_cutoff)) / (max_cutoff - min_cutoff);

  radial_basis_ders[0] = scaling * 2 * (dist - max_cutoff);
  radial_basis_ders[1] =
      scaling * (mult * (dist - max_cutoff) * (dist - max_cutoff) + 2 * ksi * (dist - max_cutoff));
  for (int i = 2; i < size; i++) {
    radial_basis_ders[i] = 2 * (mult * radial_basis_vals[i - 1] + ksi * radial_basis_ders[i - 1]) -
        radial_basis_ders[i - 2];
  }
}