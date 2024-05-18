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
#include <unordered_map>

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
    std::vector<std::string> data_files; // sorted file names
    std::vector<std::string> group_names; // sorted group names
    std::vector<std::string> filenames;
    std::string filenametag = "pod";
    std::string group_weight_type = "global";

    std::vector<int> num_atom;
    std::vector<int> num_atom_cumsum;
    std::vector<int> num_atom_each_file;
    std::vector<int> num_config;
    std::vector<int> num_config_cumsum;
    int num_atom_sum;
    int num_atom_min;
    int num_atom_max;
    int num_config_sum;

    double *lattice=nullptr;
    double *energy=nullptr;
    double *stress=nullptr;
    double *position=nullptr;
    double *force=nullptr;
    int *atomtype=nullptr;
    // Group weights will have same size as energy.
    double *we=nullptr;
    double *wf=nullptr;

    int training = 1;
    int normalizeenergy = 1;
    int training_analysis = 1;
    int test_analysis = 1;
    int training_calculation = 0;
    int test_calculation = 0;
    int randomize = 1;
    int precision = 8;
    double fraction = 1.0;

    std::unordered_map<std::string, double> we_map;
    std::unordered_map<std::string, double> wf_map;

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
      for (int i = 0; i < 12; i++)
        data.fitting_weights[i] = fitting_weights[i];
      data.we_map = we_map;
      data.wf_map = wf_map;
    }
  };

  struct neighborstruct {
    int *alist=nullptr;
    int *pairnum=nullptr;
    int *pairnum_cumsum=nullptr;
    int *pairlist=nullptr;
    double *y=nullptr;

    //int natom;
    //int nalist;
    int natom_max = 0;
    int sze = 0;
    int sza = 0;
    int szy = 0;
    int szp = 0;
  };

  struct descriptorstruct {
    double *bd=nullptr;  // base descriptors
    double *pd=nullptr;  // multi-environment descriptors (probabilities)
    double *gd=nullptr;  // global descriptors
    double *gdd=nullptr; // derivatives of global descriptors and peratom descriptors
    double *A=nullptr;  // least-square matrix for all descriptors
    double *b=nullptr;  // least-square vector for all descriptors
    double *c=nullptr;  // coefficents of descriptors
    int szd = 0;
    int nCoeffAll = 0; // number of global descriptors
    int nClusters = 0; // number of environment clusters
  };

  int save_descriptors = 0;
  int compute_descriptors = 0;
  datastruct traindata;
  datastruct testdata;
  datastruct envdata;
  descriptorstruct desc;
  neighborstruct nb;
  class EAPOD *fastpodptr;

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

  double squareDistance(const double *a, const double *b, int DIMENSIONS);
  void assignPointsToClusters(double *points, double *centroids, int *assignments, int *clusterSizes, int NUM_POINTS, int NUM_CLUSTERS, int DIMENSION);
  void updateCentroids(double *points, double *centroids, int *assignments, int *clusterSizes, int NUM_POINTS, int NUM_CLUSTERS, int DIMENSIONS);
  void KmeansClustering(double *points, double *centroids, int *assignments, int *clusterSizes, int NUM_POINTS, int NUM_CLUSTERS, int DIMENSIONS, int MAX_ITER);

  void savedata2textfile(std::string filename, std::string text, double *A, int n, int m, int dim);
  void savematrix2binfile(std::string filename, double *A, int nrows, int ncols);
  void saveintmatrix2binfile(std::string filename, int *A, int nrows, int ncols);

  // functions for reading input files and fitting

  int read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension, std::string &env_path,
    std::string &test_path, std::string &training_path, std::string &filenametag, const std::string &data_file, std::string &group_weight_type,
    std::unordered_map<std::string, double> &we_map, std::unordered_map<std::string, double> &wf_map);
  void get_exyz_files(std::vector<std::string> &, std::vector<std::string> &, const std::string &, const std::string &);
  int get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file);
  int get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files);
  void read_exyz_file(double *lattice, double *stress, double *energy, double *we, double *wf, double *pos, double *forces,
    int *atomtype, std::string file, std::vector<std::string> species, double we_group, double wf_group);
  void get_data(datastruct &data, const std::vector<std::string> &species);
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
    double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);
  void estimate_memory_neighborstruct(const datastruct &data, int *pbc, double rcut, int nelements);
  void allocate_memory_neighborstruct();
  void allocate_memory_descriptorstruct(int nd);
  void estimate_memory_fastpod(const datastruct &data);
  void local_descriptors_fastpod(const datastruct &data, int ci);
  void base_descriptors_fastpod(const datastruct &data, int ci);
  void least_squares_matrix(const datastruct &data, int ci);
  void least_squares_fit(const datastruct &data);
  void descriptors_calculation(const datastruct &data);
  void environment_cluster_calculation(const datastruct &data);
  void print_analysis(const datastruct &data, double *outarray, double *errors);
  void error_analysis(const datastruct &data, double *coeff);
  double energyforce_calculation_fastpod(double *force, const datastruct &data, int ci);
  void energyforce_calculation(const datastruct &data, double *coeff);
};

}    // namespace LAMMPS_NS

#endif
#endif
