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

#include "uf3_triplet_bspline.h"
#include "error.h"

#include <vector>

using namespace LAMMPS_NS;

// Dummy constructor
uf3_triplet_bspline::uf3_triplet_bspline(){};

// Construct a new 3D B-Spline
uf3_triplet_bspline::uf3_triplet_bspline(
    LAMMPS *ulmp, const std::vector<std::vector<double>> &uknot_matrix,
    const std::vector<std::vector<std::vector<double>>> &ucoeff_matrix,
    const int &uknot_spacing_type)
{
  lmp = ulmp;
  knot_matrix = uknot_matrix;
  coeff_matrix = ucoeff_matrix;

  knot_spacing_type = uknot_spacing_type;
  if (knot_spacing_type==0){
    knot_spacing_ij = knot_matrix[2][4]-knot_matrix[2][3];
    knot_spacing_ik = knot_matrix[1][4]-knot_matrix[1][3];
    knot_spacing_jk = knot_matrix[0][4]-knot_matrix[0][3];
    get_starting_index=&uf3_triplet_bspline::get_starting_index_uniform;
  }
  else if (knot_spacing_type==1){
    knot_spacing_ij = 0;
    knot_spacing_ik = 0;
    knot_spacing_jk = 0;
    get_starting_index=&uf3_triplet_bspline::get_starting_index_nonuniform;
  }

  else
    lmp->error->all(FLERR, "UF3: Expected either '0'(uniform-knots) or \n\
            '1'(non-uniform knots)");

  knot_vect_size_ij = knot_matrix[2].size();
  knot_vect_size_ik = knot_matrix[1].size();
  knot_vect_size_jk = knot_matrix[0].size();

  int resolution_ij = knot_vect_size_ij - 4;
  int resolution_ik = knot_vect_size_ik - 4;
  int resolution_jk = knot_vect_size_jk - 4;

  // Cache Spline Basis Functions
  for (int l = 0; l < resolution_ij; l++) {
    bsplines_ij.push_back(uf3_bspline_basis3(lmp, &knot_matrix[2][l], 1));
  }

  for (int l = 0; l < resolution_ik; l++) {
    // Reuse jk Basis if Knots match
    if (knot_matrix[1][l] == knot_matrix[2][l] && knot_matrix[1][l + 1] == knot_matrix[2][l + 1] &&
        knot_matrix[1][l + 2] == knot_matrix[2][l + 2] &&
        knot_matrix[1][l + 3] == knot_matrix[2][l + 3])
      bsplines_ik.push_back(bsplines_ij[l]);
    else
      bsplines_ik.push_back(uf3_bspline_basis3(lmp, &knot_matrix[1][l], 1));
  }

  for (int l = 0; l < resolution_jk; l++) {
    bsplines_jk.push_back(uf3_bspline_basis3(lmp, &knot_matrix[0][l], 1));
  }

  // Initialize Coefficients for Derivatives
  for (int i = 0; i < coeff_matrix.size(); i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    for (int j = 0; j < coeff_matrix[0].size(); j++) {
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size() - 1; k++) {
        double dntemp4 = 3 / (knot_matrix[0][k + 4] - knot_matrix[0][k + 1]);
        dncoeff_vect.push_back((coeff_matrix[i][j][k + 1] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_jk.push_back(dncoeff_vect2);
  }

  for (int i = 0; i < coeff_matrix.size(); i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    for (int j = 0; j < coeff_matrix[0].size() - 1; j++) {
      double dntemp4 = 3 / (knot_matrix[1][j + 4] - knot_matrix[1][j + 1]);
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size(); k++) {
        dncoeff_vect.push_back((coeff_matrix[i][j + 1][k] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_ik.push_back(dncoeff_vect2);
  }

  for (int i = 0; i < coeff_matrix.size() - 1; i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    double dntemp4 = 3 / (knot_matrix[2][i + 4] - knot_matrix[2][i + 1]);
    for (int j = 0; j < coeff_matrix[0].size(); j++) {
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size(); k++) {
        dncoeff_vect.push_back((coeff_matrix[i + 1][j][k] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_ij.push_back(dncoeff_vect2);
  }

  std::vector<std::vector<double>> dnknot_matrix;
  for (int i = 0; i < knot_matrix.size(); i++) {
    std::vector<double> dnknot_vect;
    for (int j = 1; j < knot_matrix[0].size() - 1; j++) {
      dnknot_vect.push_back(knot_matrix[i][j]);
    }
    dnknot_matrix.push_back(dnknot_vect);
  }

  // Cache Derivative Spline Basis Functions
  for (int l = 0; l < resolution_ij - 1; l++) {
    dnbsplines_ij.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[2][l], 1));
  }

  for (int l = 0; l < resolution_ik - 1; l++) {
    // Reuse jk Basis if Knots match
    if (dnknot_matrix[1][l] == dnknot_matrix[2][l] &&
        dnknot_matrix[1][l + 1] == dnknot_matrix[2][l + 1] &&
        dnknot_matrix[1][l + 2] == dnknot_matrix[2][l + 2])
      dnbsplines_ik.push_back(dnbsplines_ij[l]);
    else
      dnbsplines_ik.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[1][l], 1));
  }

  for (int l = 0; l < resolution_jk - 1; l++) {
    dnbsplines_jk.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[0][l], 1));
  }
}

// Construct a new 3D B-Spline from arrays
uf3_triplet_bspline::uf3_triplet_bspline(
    LAMMPS *ulmp, double **uknot_array, const int *uknot_array_size,
    double ***ucoeff_array, const int *ucoeff_array_size,
    const int &uknot_spacing_type)
{
  lmp = ulmp;

  knot_matrix.resize(3);
  //utils::logmesg(lmp, "knot_matrix dim = {} {} {}\nknots = ",uknot_array_size[0],
  //               uknot_array_size[1], uknot_array_size[2]);
  for (int i = 0; i < 3; i++) {
    knot_matrix[i].resize(uknot_array_size[i]);
    //utils::logmesg(lmp, "{}= ",i);
    for (int j = 0; j < uknot_array_size[i]; j++) {
      //utils::logmesg(lmp, "{} ", uknot_array[i][j]);
      knot_matrix[i][j] = uknot_array[i][j];
    }
    //utils::logmesg(lmp,"\n");
  }

  coeff_matrix.resize(ucoeff_array_size[0]);
  for (int i = 0; i < ucoeff_array_size[0]; i++) {
    coeff_matrix[i].resize(ucoeff_array_size[1]);
    for (int j = 0; j < ucoeff_array_size[1]; j++) {
      coeff_matrix[i][j].resize(ucoeff_array_size[2]);
      for (int k = 0; k < ucoeff_array_size[2]; k++){
        coeff_matrix[i][j][k] = ucoeff_array[i][j][k];
      }
    }
  }

  knot_spacing_type = uknot_spacing_type;
  if (knot_spacing_type==0){
    knot_spacing_ij = knot_matrix[2][4]-knot_matrix[2][3];
    knot_spacing_ik = knot_matrix[1][4]-knot_matrix[1][3];
    knot_spacing_jk = knot_matrix[0][4]-knot_matrix[0][3];
    get_starting_index=&uf3_triplet_bspline::get_starting_index_uniform;
  }
  else if (knot_spacing_type==1){
    knot_spacing_ij = 0;
    knot_spacing_ik = 0;
    knot_spacing_jk = 0;
    get_starting_index=&uf3_triplet_bspline::get_starting_index_nonuniform;
  }

  else
    lmp->error->all(FLERR, "UF3: Expected either '0'(uniform-knots) or \n\
            '1'(non-uniform knots)");

  knot_vect_size_ij = knot_matrix[2].size();
  knot_vect_size_ik = knot_matrix[1].size();
  knot_vect_size_jk = knot_matrix[0].size();

  int resolution_ij = knot_vect_size_ij - 4;
  int resolution_ik = knot_vect_size_ik - 4;
  int resolution_jk = knot_vect_size_jk - 4;

  // Cache Spline Basis Functions
  for (int l = 0; l < resolution_ij; l++) {
    bsplines_ij.push_back(uf3_bspline_basis3(lmp, &knot_matrix[2][l], 1));
  }

  for (int l = 0; l < resolution_ik; l++) {
    // Reuse jk Basis if Knots match
    if (knot_matrix[1][l] == knot_matrix[2][l] && knot_matrix[1][l + 1] == knot_matrix[2][l + 1] &&
        knot_matrix[1][l + 2] == knot_matrix[2][l + 2] &&
        knot_matrix[1][l + 3] == knot_matrix[2][l + 3])
      bsplines_ik.push_back(bsplines_ij[l]);
    else
      bsplines_ik.push_back(uf3_bspline_basis3(lmp, &knot_matrix[1][l], 1));
  }

  for (int l = 0; l < resolution_jk; l++) {
    bsplines_jk.push_back(uf3_bspline_basis3(lmp, &knot_matrix[0][l], 1));
  }

  // Initialize Coefficients for Derivatives
  for (int i = 0; i < coeff_matrix.size(); i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    for (int j = 0; j < coeff_matrix[0].size(); j++) {
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size() - 1; k++) {
        double dntemp4 = 3 / (knot_matrix[0][k + 4] - knot_matrix[0][k + 1]);
        dncoeff_vect.push_back((coeff_matrix[i][j][k + 1] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_jk.push_back(dncoeff_vect2);
  }

  for (int i = 0; i < coeff_matrix.size(); i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    for (int j = 0; j < coeff_matrix[0].size() - 1; j++) {
      double dntemp4 = 3 / (knot_matrix[1][j + 4] - knot_matrix[1][j + 1]);
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size(); k++) {
        dncoeff_vect.push_back((coeff_matrix[i][j + 1][k] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_ik.push_back(dncoeff_vect2);
  }

  for (int i = 0; i < coeff_matrix.size() - 1; i++) {
    std::vector<std::vector<double>> dncoeff_vect2;
    double dntemp4 = 3 / (knot_matrix[2][i + 4] - knot_matrix[2][i + 1]);
    for (int j = 0; j < coeff_matrix[0].size(); j++) {
      std::vector<double> dncoeff_vect;
      for (int k = 0; k < coeff_matrix[0][0].size(); k++) {
        dncoeff_vect.push_back((coeff_matrix[i + 1][j][k] - coeff_matrix[i][j][k]) * dntemp4);
      }
      dncoeff_vect2.push_back(dncoeff_vect);
    }
    dncoeff_matrix_ij.push_back(dncoeff_vect2);
  }

  std::vector<std::vector<double>> dnknot_matrix;
  for (int i = 0; i < knot_matrix.size(); i++) {
    std::vector<double> dnknot_vect;
    for (int j = 1; j < knot_matrix[0].size() - 1; j++) {
      dnknot_vect.push_back(knot_matrix[i][j]);
    }
    dnknot_matrix.push_back(dnknot_vect);
  }

  // Cache Derivative Spline Basis Functions
  for (int l = 0; l < resolution_ij - 1; l++) {
    dnbsplines_ij.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[2][l], 1));
  }

  for (int l = 0; l < resolution_ik - 1; l++) {
    // Reuse jk Basis if Knots match
    if (dnknot_matrix[1][l] == dnknot_matrix[2][l] &&
        dnknot_matrix[1][l + 1] == dnknot_matrix[2][l + 1] &&
        dnknot_matrix[1][l + 2] == dnknot_matrix[2][l + 2])
      dnbsplines_ik.push_back(dnbsplines_ij[l]);
    else
      dnbsplines_ik.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[1][l], 1));
  }

  for (int l = 0; l < resolution_jk - 1; l++) {
    dnbsplines_jk.push_back(uf3_bspline_basis2(lmp, &dnknot_matrix[0][l], 1));
  }
}


// Destructor
uf3_triplet_bspline::~uf3_triplet_bspline() {}

// Evaluate 3D B-Spline value
double *uf3_triplet_bspline::eval(double value_rij, double value_rik, double value_rjk)
{

  // Find starting knots

  //int iknot_ij = starting_knot(knot_matrix[2], knot_vect_size_ij, value_rij) - 3;
  //int iknot_ik = starting_knot(knot_matrix[1], knot_vect_size_ik, value_rik) - 3;
  //int iknot_jk = starting_knot(knot_matrix[0], knot_vect_size_jk, value_rjk) - 3;
  int iknot_ij = (this->*get_starting_index)(knot_matrix[2], knot_vect_size_ij, value_rij,knot_spacing_ij) - 3;
  int iknot_ik = (this->*get_starting_index)(knot_matrix[1], knot_vect_size_ik, value_rik,knot_spacing_ik) - 3;
  int iknot_jk = (this->*get_starting_index)(knot_matrix[0], knot_vect_size_jk, value_rjk,knot_spacing_jk) - 3;

  double rsq_ij = value_rij * value_rij;
  double rsq_ik = value_rik * value_rik;
  double rsq_jk = value_rjk * value_rjk;
  double rth_ij = rsq_ij * value_rij;
  double rth_ik = rsq_ik * value_rik;
  double rth_jk = rsq_jk * value_rjk;

  // Calculate energies

  double basis_ij[4];
  basis_ij[0] = bsplines_ij[iknot_ij].eval3(rth_ij, rsq_ij, value_rij);
  basis_ij[1] = bsplines_ij[iknot_ij + 1].eval2(rth_ij, rsq_ij, value_rij);
  basis_ij[2] = bsplines_ij[iknot_ij + 2].eval1(rth_ij, rsq_ij, value_rij);
  basis_ij[3] = bsplines_ij[iknot_ij + 3].eval0(rth_ij, rsq_ij, value_rij);

  double basis_ik[4];
  basis_ik[0] = bsplines_ik[iknot_ik].eval3(rth_ik, rsq_ik, value_rik);
  basis_ik[1] = bsplines_ik[iknot_ik + 1].eval2(rth_ik, rsq_ik, value_rik);
  basis_ik[2] = bsplines_ik[iknot_ik + 2].eval1(rth_ik, rsq_ik, value_rik);
  basis_ik[3] = bsplines_ik[iknot_ik + 3].eval0(rth_ik, rsq_ik, value_rik);

  double basis_jk[4];
  basis_jk[0] = bsplines_jk[iknot_jk].eval3(rth_jk, rsq_jk, value_rjk);
  basis_jk[1] = bsplines_jk[iknot_jk + 1].eval2(rth_jk, rsq_jk, value_rjk);
  basis_jk[2] = bsplines_jk[iknot_jk + 2].eval1(rth_jk, rsq_jk, value_rjk);
  basis_jk[3] = bsplines_jk[iknot_jk + 3].eval0(rth_jk, rsq_jk, value_rjk);

  ret_val[0] = 0;
  ret_val[1] = 0;
  ret_val[2] = 0;
  ret_val[3] = 0;

  for (int i = 0; i < 4; i++) {
    const double basis_iji = basis_ij[i]; // prevent repeated access of same memory location
    for (int j = 0; j < 4; j++) {
      const double factor = basis_iji * basis_ik[j]; // prevent repeated access of same memory location
      const double* slice = &coeff_matrix[i + iknot_ij][j + iknot_ik][iknot_jk]; // declare a contigues 1D slice of memory
      double tmp[4]; // declare tmp array that holds the 4 tmp values so the can be computed simultaniously in 4 separate registeres.
      tmp[0] = slice[0] * basis_jk[0];
      tmp[1] = slice[1] * basis_jk[1];
      tmp[2] = slice[2] * basis_jk[2];
      tmp[3] = slice[3] * basis_jk[3];
      double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
      ret_val[0] += factor * sum; // use 1 fused multiply-add (FMA)
    }
  }

  // Calculate forces

  double dnbasis_ij[4];
  dnbasis_ij[0] = dnbsplines_ij[iknot_ij].eval2(rsq_ij, value_rij);
  dnbasis_ij[1] = dnbsplines_ij[iknot_ij + 1].eval1(rsq_ij, value_rij);
  dnbasis_ij[2] = dnbsplines_ij[iknot_ij + 2].eval0(rsq_ij, value_rij);
  dnbasis_ij[3] = 0;

  double dnbasis_ik[4];
  dnbasis_ik[0] = dnbsplines_ik[iknot_ik].eval2(rsq_ik, value_rik);
  dnbasis_ik[1] = dnbsplines_ik[iknot_ik + 1].eval1(rsq_ik, value_rik);
  dnbasis_ik[2] = dnbsplines_ik[iknot_ik + 2].eval0(rsq_ik, value_rik);
  dnbasis_ik[3] = 0;

  double dnbasis_jk[4];
  dnbasis_jk[0] = dnbsplines_jk[iknot_jk].eval2(rsq_jk, value_rjk);
  dnbasis_jk[1] = dnbsplines_jk[iknot_jk + 1].eval1(rsq_jk, value_rjk);
  dnbasis_jk[2] = dnbsplines_jk[iknot_jk + 2].eval0(rsq_jk, value_rjk);
  dnbasis_jk[3] = 0;

  for (int i = 0; i < 3; i++) {
    const double dnbasis_iji = dnbasis_ij[i];
    for (int j = 0; j < 4; j++) {
      const double factor = dnbasis_iji * basis_ik[j];
      const double* slice = &dncoeff_matrix_ij[iknot_ij + i][iknot_ik + j][iknot_jk];
      double tmp[4];
      tmp[0] = slice[0] * basis_jk[0];
      tmp[1] = slice[1] * basis_jk[1];
      tmp[2] = slice[2] * basis_jk[2];
      tmp[3] = slice[3] * basis_jk[3];
      double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
      ret_val[1] += factor * sum;
    }
  }

  for (int i = 0; i < 4; i++) {
    const double basis_iji = basis_ij[i];
    for (int j = 0; j < 3; j++) {
      const double factor = basis_iji * dnbasis_ik[j];
      const double* slice = &dncoeff_matrix_ik[iknot_ij + i][iknot_ik + j][iknot_jk];
      double tmp[4];
      tmp[0] = slice[0] * basis_jk[0];
      tmp[1] = slice[1] * basis_jk[1];
      tmp[2] = slice[2] * basis_jk[2];
      tmp[3] = slice[3] * basis_jk[3];
      double sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
      ret_val[2] += factor * sum;
    }
  }

  for (int i = 0; i < 4; i++) {
    const double basis_iji = basis_ij[i];
    for (int j = 0; j < 4; j++) {
      const double factor = basis_iji * basis_ik[j];
      const double* slice = &dncoeff_matrix_jk[iknot_ij + i][iknot_ik + j][iknot_jk];
      double tmp[3];
      tmp[0] = slice[0] * dnbasis_jk[0];
      tmp[1] = slice[1] * dnbasis_jk[1];
      tmp[2] = slice[2] * dnbasis_jk[2];
      double sum = tmp[0] + tmp[1] + tmp[2];
      ret_val[3] += factor * sum;
    }
  }

  return ret_val;
}

// Find starting knot for spline evaluation

int uf3_triplet_bspline::starting_knot(const std::vector<double> knot_vect, int knot_vect_size,
                                       double r)
{
  if (knot_vect.front() <= r && r < knot_vect.back()) {
    for (int i = 3; i < knot_vect_size - 1; i++) {
      if (knot_vect[i] <= r && r < knot_vect[i + 1]) return i;
    }
  }

  return 0;
}

int uf3_triplet_bspline::get_starting_index_uniform(const std::vector<double> knot_vect, int knot_vect_size,
                                                    double r, double knot_spacing)
{
  return 3+(int)((r-knot_vect[0])/knot_spacing);
}

int uf3_triplet_bspline::get_starting_index_nonuniform(const std::vector<double> knot_vect, int knot_vect_size,
                                                       double r, double knot_spacing)
{
  if (knot_vect.front() <= r && r < knot_vect.back()) {
    //Determine the interval for value_rij
    for (int i = 3; i < knot_vect_size - 1; ++i) {
      if (knot_vect[i] <= r && r < knot_vect[i + 1]) {
        return i;
      }
    }
  }
  return -1;
}

double uf3_triplet_bspline::memory_usage()
{
  double bytes = 0;

  bytes += (double) 3*sizeof(int);          //knot_vect_size_ij,
                                            //knot_vect_size_ik,
                                            //knot_vect_size_jk;

  for (int i=0; i<coeff_matrix.size(); i++)
    for (int j=0; j<coeff_matrix[i].size(); j++)
      bytes += (double)coeff_matrix[i][j].size()*sizeof(double);


  for (int i=0; i<dncoeff_matrix_ij.size(); i++)
    for (int j=0; j<dncoeff_matrix_ij[i].size(); j++)
      bytes += (double)dncoeff_matrix_ij[i][j].size()*sizeof(double);

  for (int i=0; i<dncoeff_matrix_ik.size(); i++)
    for (int j=0; j<dncoeff_matrix_ik[i].size(); j++)
      bytes += (double)dncoeff_matrix_ik[i][j].size()*sizeof(double);

  for (int i=0; i<dncoeff_matrix_jk.size(); i++)
    for (int j=0; j<dncoeff_matrix_jk[i].size(); j++)
      bytes += (double)dncoeff_matrix_jk[i][j].size()*sizeof(double);

  bytes += (double)knot_matrix[0].size()*sizeof(double);
  bytes += (double)knot_matrix[1].size()*sizeof(double);
  bytes += (double)knot_matrix[2].size()*sizeof(double);

  for (int i=0; i < bsplines_ij.size(); i++)
    bytes += (double)bsplines_ij[i].memory_usage();

  for (int i=0; i < bsplines_ik.size(); i++)
    bytes += (double)bsplines_ik[i].memory_usage();

  for (int i=0; i < bsplines_jk.size(); i++)
    bytes += (double)bsplines_jk[i].memory_usage();

  for (int i=0; i < dnbsplines_ij.size(); i++)
    bytes += (double)dnbsplines_ij[i].memory_usage();

  for (int i=0; i < dnbsplines_ik.size(); i++)
    bytes += (double)dnbsplines_ik[i].memory_usage();

  for (int i=0; i < dnbsplines_jk.size(); i++)
    bytes += (double)dnbsplines_jk[i].memory_usage();

  bytes += (double)4*sizeof(double);                        //ret_val
  return bytes;
}
