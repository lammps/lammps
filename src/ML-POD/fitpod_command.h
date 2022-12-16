/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(fitpod,FitPOD);
// clang-format on
#else

#ifndef LMP_FITPOD_COMMAND_H
#define LMP_FITPOD_COMMAND_H

#include "command.h"

namespace LAMMPS_NS {

class FitPOD : public Command {
 public:
  FitPOD(LAMMPS *);
  void command(int, char **) override;

 private:
  struct datastruct {
    std::string file_format = "extxyz";
    std::string file_extension = "xyz";
    std::string data_path;
    std::vector<std::string> data_files;
    std::vector<std::string> filenames;
    std::string filenametag = "pod";

    std::vector<int> num_atom;
    std::vector<int> num_atom_cumsum;
    std::vector<int> num_atom_each_file;
    std::vector<int> num_config;
    std::vector<int> num_config_cumsum;
    int num_atom_sum;
    int num_atom_min;
    int num_atom_max;
    int num_config_sum;

    double *lattice;
    double *energy;
    double *stress;
    double *position;
    double *force;
    int *atomtype;

    int training = 1;
    int normalizeenergy = 1;
    int training_analysis = 1;
    int test_analysis = 1;
    int training_calculation = 0;
    int test_calculation = 0;
    int randomize = 1;
    int precision = 8;
    double fraction = 1.0;

    double fitting_weights[12] = {100.0, 1.0, 0.0, 1, 1, 0, 0, 1, 1, 1, 1, 1e-10};

    void copydatainfo(datastruct &data) const
    {
      data.data_path = data_path;
      data.file_format = file_format;
      data.file_extension = file_extension;
      data.data_files = data_files;
      data.filenametag = filenametag;
      data.filenames = filenames;
      data.training_analysis = training_analysis;
      data.test_analysis = test_analysis;
      data.training_calculation = training_calculation;
      data.test_calculation = test_calculation;
      data.fraction = fraction;
      data.randomize = randomize;
      data.precision = precision;
      data.training = training;
      data.normalizeenergy = normalizeenergy;
      for (int i = 0; i < 12; i++) data.fitting_weights[i] = fitting_weights[i];
    }
  };

  struct neighborstruct {
    int *alist;
    int *pairnum;
    int *pairnum_cumsum;
    int *pairlist;
    double *y;

    int natom;
    int nalist;
    int natom_max;
    int sze;
    int sza;
    int szy;
    int szp;
  };

  struct descriptorstruct {
    double *gd;     // global descriptors
    double *gdd;    // derivatives of global descriptors and peratom descriptors
    double *A;      // least-square matrix for all descriptors
    double *b;      // least-square vector for all descriptors
    double *c;      // coefficents of descriptors
    int *tmpint;
    int szd;
    int szi;
  };

  datastruct traindata;
  datastruct testdata;
  descriptorstruct desc;
  neighborstruct nb;
  class MLPOD *podptr;

  // functions for collecting/collating arrays

  void print_matrix(const char *desc, int m, int n, int *a, int lda);
  void print_matrix(const char *desc, int m, int n, double *a, int lda);
  void print_matrix(const char *desc, int m, int n, double **a, int lda);
  void podCumsum(int *output, int *input, int length);
  double podArrayNorm(double *a, int n);
  double podArrayErrorNorm(double *a, double *b, int n);
  void podArraySetValue(double *y, double a, int n);
  void podArrayCopy(double *y, double *x, int n);
  void podArrayFill(int *output, int start, int length);
  double podArrayMin(double *a, int n);
  double podArrayMax(double *a, int n);
  double podArraySum(double *a, int n);
  int podArrayMin(int *a, int n);
  int podArrayMax(int *a, int n);
  void podKron(double *C, double *A, double *B, double alpha, int M1, int M2);
  void rotation_matrix(double *Rmat, double alpha, double beta, double gamma);
  void triclinic_lattice_conversion(double *a, double *b, double *c, double *A, double *B,
                                    double *C);
  void matrix33_multiplication(double *xrot, double *Rmat, double *x, int natom);
  void matrix33_inverse(double *invA, double *A1, double *A2, double *A3);

  // functions for reading input files and fitting

  int read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension,
                     std::string &test_path, std::string &training_path, std::string &filenametag,
                     const std::string &data_file);
  void get_exyz_files(std::vector<std::string> &, const std::string &, const std::string &);
  int get_number_atom_exyz(std::vector<int> &num_atom, int &num_atom_sum, std::string file);
  int get_number_atoms(std::vector<int> &num_atom, std::vector<int> &num_atom_sum,
                       std::vector<int> &num_config, std::vector<std::string> training_files);
  void read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces,
                      int *atomtype, std::string file, std::vector<std::string> species);
  void get_data(datastruct &data, const std::vector<std::string>& species);
  std::vector<int> linspace(int start_in, int end_in, int num_in);
  std::vector<int> shuffle(int start_in, int end_in, int num_in);
  std::vector<int> select(int n, double fraction, int randomize);
  void select_data(datastruct &newdata, const datastruct &data);
  void read_data_files(const std::string& data_file, const std::vector<std::string>& species);
  int latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3,
                    double rcut, int *pbc, int nx);
  int podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N,
                      int dim);
  int podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum,
                          double *x, double *a1, double *a2, double *a3, double rcut, int *pbc,
                          int nx);
  void allocate_memory(const datastruct &data);
  void linear_descriptors(const datastruct &data, int ci);
  void quadratic_descriptors(const datastruct &data, int ci);
  void cubic_descriptors(const datastruct &data, int ci);
  void least_squares_matrix(const datastruct &data, int ci);
  void least_squares_fit(const datastruct &data);
  void print_analysis(const datastruct &data, double *outarray, double *errors);
  void error_analysis(const datastruct &data, double *coeff);
  double energyforce_calculation(double *force, double *coeff, const datastruct &data, int ci);
  void energyforce_calculation(const datastruct &data, double *coeff);
};

}    // namespace LAMMPS_NS

#endif
#endif
