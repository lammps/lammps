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

#include "uf3_pair_bspline.h"

#include "uf3_bspline_basis2.h"
#include "uf3_bspline_basis3.h"

#include "utils.h"
#include "error.h"
#include <vector>

using namespace LAMMPS_NS;

// Dummy constructor
uf3_pair_bspline::uf3_pair_bspline() {}

// Constructor
// Passing vectors by reference
uf3_pair_bspline::uf3_pair_bspline(LAMMPS *ulmp, const std::vector<double> &uknot_vect,
                                   const std::vector<double> &ucoeff_vect,
                                   const int &uknot_spacing_type)
{
  lmp = ulmp;
  knot_vect = uknot_vect;
  coeff_vect = ucoeff_vect;

  knot_spacing_type = uknot_spacing_type;
  if (knot_spacing_type==0){
    knot_spacing = knot_vect[4]-knot_vect[3];
    get_starting_index=&uf3_pair_bspline::get_starting_index_uniform;
  }
  else if (knot_spacing_type==1){
    knot_spacing = 0;
    get_starting_index=&uf3_pair_bspline::get_starting_index_nonuniform;
  }

  else
    lmp->error->all(FLERR, "UF3: Expected either '0'(uniform-knots) or \n\
            '1'(non-uniform knots)");

  knot_vect_size = uknot_vect.size();
  coeff_vect_size = ucoeff_vect.size();

  // Initialize B-Spline Basis Functions
  for (int i = 0; i < knot_vect.size() - 4; i++)
    bspline_bases.push_back(uf3_bspline_basis3(lmp, &knot_vect[i], coeff_vect[i]));

  // Initialize Coefficients and Knots for Derivatives
  // The last coefficient needs to be droped
  for (int i = 0; i < coeff_vect_size - 1; i++) {
    double dntemp4 = 3 / (knot_vect[i + 4] - knot_vect[i + 1]);
    dncoeff_vect.push_back((coeff_vect[i + 1] - coeff_vect[i]) * dntemp4);
  }
  //What we have is a clamped bspline -->i.e value of the bspline curve at the
  //knots with multiplicity equal to the degree of bspline is equal to the coefficient
  //
  //Therefore for the derivative bspline the very first and last knot needs to be droped
  //to change their multiplicity from 4 (necessary condition for clamped cubic bspline)
  //to 3 (necessary condition for clamped quadratic bspline)
  //
  //Also if the coeff vector size of decreases by 1 for the derivative bspline
  //knots size needs to go down by 2 as ==> knots = coefficient + degree + 1
  for (int i = 1; i < knot_vect_size - 1; i++) dnknot_vect.push_back(knot_vect[i]);

  // Initialize B-Spline Derivative Basis Functions
  for (int i = 0; i < dnknot_vect.size() - 3; i++)
    dnbspline_bases.push_back(uf3_bspline_basis2(lmp, &dnknot_vect[i], dncoeff_vect[i]));
}

// Constructor
// Passing arrays
uf3_pair_bspline::uf3_pair_bspline(LAMMPS *ulmp, const double* uknot_array,
                                   const int uknot_array_size,
                                   const double* ucoeff_array,
                                   const int ucoeff_array_size,
                                   const int uknot_spacing_type)
{
  lmp = ulmp;
  
  knot_vect = std::vector<double> (uknot_array, uknot_array + uknot_array_size);
  coeff_vect = std::vector<double> (ucoeff_array, ucoeff_array + ucoeff_array_size);

  knot_spacing_type = uknot_spacing_type;
  if (knot_spacing_type==0){
    knot_spacing = knot_vect[4]-knot_vect[3];
    get_starting_index=&uf3_pair_bspline::get_starting_index_uniform;
  }
  else if (knot_spacing_type==1){
    knot_spacing = 0;
    get_starting_index=&uf3_pair_bspline::get_starting_index_nonuniform;
  }

  else
    lmp->error->all(FLERR, "UF3: Expected either '0'(uniform-knots) or \n\
            '1'(non-uniform knots)");

  knot_vect_size = uknot_array_size;
  coeff_vect_size = ucoeff_array_size;

  // Initialize B-Spline Basis Functions
  for (int i = 0; i < knot_vect.size() - 4; i++)
    bspline_bases.push_back(uf3_bspline_basis3(lmp, &knot_vect[i], coeff_vect[i]));

  // Initialize Coefficients and Knots for Derivatives
  // The last coefficient needs to be droped
  for (int i = 0; i < coeff_vect_size - 1; i++) {
    double dntemp4 = 3 / (knot_vect[i + 4] - knot_vect[i + 1]);
    dncoeff_vect.push_back((coeff_vect[i + 1] - coeff_vect[i]) * dntemp4);
  }
  //What we have is a clamped bspline -->i.e value of the bspline curve at the
  //knots with multiplicity equal to the degree of bspline is equal to the coefficient
  //
  //Therefore for the derivative bspline the very first and last knot needs to be droped
  //to change their multiplicity from 4 (necessary condition for clamped cubic bspline)
  //to 3 (necessary condition for clamped quadratic bspline)
  //
  //Also if the coeff vector size of decreases by 1 for the derivative bspline
  //knots size needs to go down by 2 as ==> knots = coefficient + degree + 1
  for (int i = 1; i < knot_vect_size - 1; i++) dnknot_vect.push_back(knot_vect[i]);

  // Initialize B-Spline Derivative Basis Functions
  for (int i = 0; i < dnknot_vect.size() - 3; i++)
    dnbspline_bases.push_back(uf3_bspline_basis2(lmp, &dnknot_vect[i], dncoeff_vect[i]));
}

uf3_pair_bspline::~uf3_pair_bspline() {}

int uf3_pair_bspline::get_starting_index_uniform(double r)
{
  return 3+(int)((r-knot_vect[0])/knot_spacing);
}

int uf3_pair_bspline::get_starting_index_nonuniform(double r)
{
  if (knot_vect.front() <= r && r < knot_vect.back()) {
    //Determine the interval for value_rij
    for (int i = 3; i < knot_vect_size - 1; ++i) {
      if (knot_vect[i] <= r && r < knot_vect[i + 1]) {
        return i;
      }
    }
  }
}

double *uf3_pair_bspline::eval(double r)
{

  // Find knot starting position

  int start_index=(this->*get_starting_index)(r);
  /*if (knot_vect.front() <= r && r < knot_vect.back()) {
    //Determine the interval for value_rij
    for (int i = 3; i < knot_vect_size - 1; ++i) {
      if (knot_vect[i] <= r && r < knot_vect[i + 1]) {
        start_index = i;
        break;
      }
    }
  }*/

  int knot_affect_start = start_index - 3;

  double rsq = r * r;
  double rth = rsq * r;

  // Calculate energy

  ret_val[0] = bspline_bases[knot_affect_start + 3].eval0(rth, rsq, r);
  ret_val[0] += bspline_bases[knot_affect_start + 2].eval1(rth, rsq, r);
  ret_val[0] += bspline_bases[knot_affect_start + 1].eval2(rth, rsq, r);
  ret_val[0] += bspline_bases[knot_affect_start].eval3(rth, rsq, r);

  // Calculate force

  ret_val[1] = dnbspline_bases[knot_affect_start + 2].eval0(rsq, r);
  ret_val[1] += dnbspline_bases[knot_affect_start + 1].eval1(rsq, r);
  ret_val[1] += dnbspline_bases[knot_affect_start].eval2(rsq, r);

  return ret_val;
}

double uf3_pair_bspline::memory_usage()
{
  double bytes = 0;

  bytes += (double)2*sizeof(int);                           //knot_vect_size,
                                                            //coeff_vect_size
  bytes += (double)knot_vect.size()*sizeof(double);         //knot_vect
  bytes += (double)dnknot_vect.size()*sizeof(double);       //dnknot_vect
  bytes += (double)coeff_vect.size()*sizeof(double);        //coeff_vect
  bytes += (double)dncoeff_vect.size()*sizeof(double);      //dncoeff_vect

  for (int i = 0; i < knot_vect.size() - 4; i++)
    bytes += (double)bspline_bases[i].memory_usage();       //bspline_basis3

  for (int i = 0; i < dnknot_vect.size() - 3; i++)
    bytes += (double)dnbspline_bases[i].memory_usage();     //bspline_basis2

  bytes += (double)2*sizeof(double);    //ret_val

  return bytes;
}
