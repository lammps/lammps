/* ----------------------------------------------------------------------
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
CommandStyle(fitpod,CFITPOD);
// clang-format on
#else

#ifndef LMP_FITPOD_H
#define LMP_FITPOD_H

#include "command.h"

namespace LAMMPS_NS {

class CFITPOD : public Command {
private:

public:
  struct datastruct {
    std::string file_format;
    std::string file_extension;
    std::string data_path;
    std::vector<std::string> data_files;
    std::vector<std::string> filenames;

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
    double fraction = 1.0;

    double fitting_weights[12] = {0.0, 0.0, 0.0, 1, 1, 0, 0, 1, 1, 1, 1, 0};

    void copydatainfo(datastruct &data) {
      data.data_path = data_path;
      data.file_format = file_format;
      data.file_extension = file_extension;
      data.data_files = data_files;
      data.filenames = filenames;
      data.training_analysis = training_analysis;
      data.test_analysis = test_analysis;
      data.training_calculation = training_calculation;
      data.test_calculation = test_calculation;
      data.fraction = fraction;
      data.randomize = randomize;
      data.training = training;
      data.normalizeenergy = normalizeenergy;
      for (int i = 0; i < 12; i++)
        data.fitting_weights[i] = fitting_weights[i];
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
    double *gd;  // global descriptors
    double *gdd; // derivatives of global descriptors and peratom descriptors
    double *A;  // least-square matrix for all descriptors
    double *b;  // least-square vector for all descriptors
    double *c;  // coefficents of descriptors
    int *tmpint;
    int szd;
    int szi;
  };

  datastruct traindata;
  datastruct testdata;
  descriptorstruct desc;
  neighborstruct nb;
  class MLPOD *podptr;

  CFITPOD(LAMMPS *lmp) : Command(lmp) {}

  void command(int, char **) override;
  void read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension,
    std::string &test_path, std::string &training_path, const std::string &data_file);
  void get_exyz_files(std::vector<std::string> &, const std::string &, const std::string &);
  int get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file);
  int get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files);
  void read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces,
    int *atomtype, std::string file, std::vector<std::string> species);
  void get_data(datastruct &data, std::vector<std::string> species);
  std::vector<int> linspace(int start_in, int end_in, int num_in);
  std::vector<int> shuffle(int start_in, int end_in, int num_in);
  std::vector<int> select(int n, double fraction, int randomize);
  void select_data(datastruct &newdata, datastruct data);
  void read_data_files(std::string data_file, std::vector<std::string> species);
  int latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);
  int podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim);
  int podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum,
    double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);
  void allocate_memory(datastruct data);
  void linear_descriptors(datastruct data, int ci);
  void quadratic_descriptors(datastruct data, int ci);
  void cubic_descriptors(datastruct data, int ci);
  void least_squares_matrix(datastruct data, int ci);
  void least_squares_fit(datastruct data);
  void print_analysis(datastruct data, double *outarray, double *errors);
  void error_analysis(datastruct data, double *coeff);
  double energyforce_calculation(double *force, double *coeff, datastruct data, int ci);
  void energyforce_calculation(datastruct data, double *coeff);
};

}  // namespace LAMMPS_NS

#endif
#endif

