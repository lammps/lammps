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

/* ----------------------------------------------------------------------
   Contributing authors: Ngoc Cuong Nguyen (MIT) and Andrew Rohskopf (SNL)
------------------------------------------------------------------------- */

#include "fitpod_command.h"

#include "mlpod.h"

#include "comm.h"
#include "error.h"
#include "math_special.h"
#include "memory.h"
#include "tokenizer.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>

using namespace LAMMPS_NS;
using MathSpecial::powint;

static constexpr int MAXLINE = 1024;
static constexpr double SMALL = 1.0e-10;

FitPOD::FitPOD(LAMMPS *_lmp) : Command(_lmp), podptr(nullptr)
{
}

void FitPOD::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "fitpod", error);

  std::string pod_file = std::string(arg[0]);  // pod input file
  std::string data_file = std::string(arg[1]); // data input file
  std::string coeff_file; // coefficient input file

  if (narg > 2)
    coeff_file = std::string(arg[2]); // coefficient input file
  else
    coeff_file = "";

  podptr = new MLPOD(lmp, pod_file, coeff_file);
  read_data_files(data_file, podptr->pod.species);

  if ((int) traindata.data_path.size() > 1)
    allocate_memory(traindata);
  else if ((int) testdata.data_path.size() > 1)
    allocate_memory(testdata);

  // get POD coefficients from an input file

  if (coeff_file != "") podArrayCopy(desc.c, podptr->pod.coeff, podptr->pod.nd);

  // compute POD coefficients using least-squares method

  least_squares_fit(traindata);

  // calculate errors for the training data set

  if ((traindata.training_analysis) && ((int) traindata.data_path.size() > 1) )
    error_analysis(traindata, desc.c);

  // calculate errors for the test data set

  if ((testdata.test_analysis) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path))
    error_analysis(testdata, desc.c);

  // calculate energy and force for the training data set

  if ((traindata.training_calculation) && ((int) traindata.data_path.size() > 1) )
    energyforce_calculation(traindata, desc.c);

  // calculate energy and force for the test data set

  if ((testdata.test_calculation) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path) )
    energyforce_calculation(testdata, desc.c);

  // deallocate training data

  if ((int) traindata.data_path.size() > 1){
    memory->destroy(traindata.lattice);
    memory->destroy(traindata.energy);
    memory->destroy(traindata.stress);
    memory->destroy(traindata.position);
    memory->destroy(traindata.force);
    memory->destroy(traindata.atomtype);
  }

  // deallocate testing data

  if ((int) testdata.data_path.size() > 1 && (testdata.data_path != traindata.data_path)){
    memory->destroy(testdata.lattice);
    memory->destroy(testdata.energy);
    memory->destroy(testdata.stress);
    memory->destroy(testdata.position);
    memory->destroy(testdata.force);
    memory->destroy(testdata.atomtype);
  }

  // deallocate descriptors

  memory->destroy(desc.gd);
  memory->destroy(desc.gdd);
  memory->destroy(desc.A);
  memory->destroy(desc.b);
  memory->destroy(desc.c);
  memory->destroy(desc.tmpint);

  // deallocate neighbor data

  memory->destroy(nb.alist);
  memory->destroy(nb.pairnum);
  memory->destroy(nb.pairnum_cumsum);
  memory->destroy(nb.pairlist);
  memory->destroy(nb.y);
  delete podptr;
}

/* ---------------------------------------------------------------------- */

int FitPOD::read_data_file(double *fitting_weights, std::string &file_format,
                             std::string &file_extension, std::string &test_path,
                             std::string &training_path, std::string &filenametag,
                             const std::string &data_file)
{
  int precision = 8;

  std::string datafilename = data_file;
  FILE *fpdata;
  if (comm->me == 0) {

    fpdata = utils::open_potential(datafilename,lmp,nullptr);
    if (fpdata == nullptr)
      error->one(FLERR,"Cannot open training data file {}: ", datafilename, utils::getsyserror());
  }

  // loop through lines of training data file and parse keywords

  char line[MAXLINE] = {'\0'};
  char *ptr;
  int eof = 0;
  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpdata);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpdata);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // words = ptrs to all words in line
    // strip single and double quotes from words

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &) {
      // ignore
    }

    if (words.size() == 0) continue;

    const auto &keywd = words[0];

    if (words.size() != 2)
      error->one(FLERR,"Improper POD file.", utils::getsyserror());

    // settings for fitting weights

    if (keywd == "fitting_weight_energy") fitting_weights[0] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "fitting_weight_force") fitting_weights[1] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "fitting_weight_stress") fitting_weights[2] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "error_analysis_for_training_data_set") fitting_weights[3] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "error_analysis_for_test_data_set") fitting_weights[4] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "energy_force_calculation_for_training_data_set") fitting_weights[5] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "energy_force_calculation_for_test_data_set") fitting_weights[6] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "fraction_training_data_set") fitting_weights[7] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "fraction_test_data_set") fitting_weights[8] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "randomize_training_data_set") fitting_weights[9] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "randomize_test_data_set") fitting_weights[10] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "fitting_regularization_parameter") fitting_weights[11] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "precision_for_pod_coefficients") precision = utils::inumeric(FLERR,words[1],false,lmp);

    // other settings

    if (keywd == "file_format") file_format = words[1];
    if (keywd == "file_extension") file_extension = words[1];
    if (keywd == "path_to_training_data_set") training_path = words[1];
    if (keywd == "path_to_test_data_set") test_path = words[1];
    if (keywd == "basename_for_output_files") filenametag = words[1];
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "**************** Begin of Data File ****************\n");
    utils::logmesg(lmp, "file format: {}\n", file_format);
    utils::logmesg(lmp, "file extension: {}\n", file_extension);
    utils::logmesg(lmp, "path to training data set: {}\n", training_path);
    utils::logmesg(lmp, "path to test data set: {}\n", test_path);
    utils::logmesg(lmp, "training fraction: {}\n", fitting_weights[7]);
    utils::logmesg(lmp, "test fraction: {}\n", fitting_weights[8]);
    utils::logmesg(lmp, "randomize training data set: {}\n", fitting_weights[9]);
    utils::logmesg(lmp, "randomize test data set: {}\n", fitting_weights[10]);
    utils::logmesg(lmp, "error analysis for training data set: {}\n", fitting_weights[3]);
    utils::logmesg(lmp, "error analysis for test data set: {}\n", fitting_weights[4]);
    utils::logmesg(lmp, "energy/force calculation for training data set: {}\n", fitting_weights[5]);
    utils::logmesg(lmp, "energy/force calculation for test data set: {}\n", fitting_weights[6]);
    utils::logmesg(lmp, "fitting weight for energy: {}\n", fitting_weights[0]);
    utils::logmesg(lmp, "fitting weight for force: {}\n", fitting_weights[1]);
    utils::logmesg(lmp, "fitting weight for stress: {}\n", fitting_weights[2]);
    utils::logmesg(lmp, "fitting regularization parameter: {}\n", fitting_weights[11]);
    utils::logmesg(lmp, "**************** End of Data File ****************\n");
  }

  return precision;
}

void FitPOD::get_exyz_files(std::vector<std::string>& files, const std::string &datapath,
                             const std::string &extension)
{
  auto allfiles = platform::list_directory(datapath);
  std::sort(allfiles.begin(), allfiles.end());
  for (const auto &fname : allfiles) {
    if (utils::strmatch(fname, fmt::format(".*\\.{}$", extension)))
      files.push_back(datapath + platform::filepathsep + fname);
  }
}

int FitPOD::get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)
{
  std::string filename = std::move(file);
  FILE *fp;
  if (comm->me == 0) {
    fp = utils::open_potential(filename,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ", filename, utils::getsyserror());
  }

  char line[MAXLINE] = {'\0'};
  char *ptr;
  int eof = 0;
  int num_configs = 0;
  num_atom_sum = 0;

  // loop over all lines of this xyz file and extract number of atoms and number of configs

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fp);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // words = ptrs to all words in line
    // strip single and double quotes from words

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &) {
      // ignore
    }

    if (words.size() == 0) continue;

    int natom;
    if (words.size() == 1) {
      natom = utils::inumeric(FLERR,words[0],false,lmp);
      num_atom.push_back(natom);
      num_configs += 1;
      num_atom_sum += natom;
    }
  }
  return num_configs;
}

int FitPOD::get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)
{
  int nfiles = training_files.size(); // number of files
  int d, n;

  for (int i=0; i<nfiles; i++) {
    d = get_number_atom_exyz(num_atom, n, training_files[i]);
    num_config.push_back(d);
    num_atom_sum.push_back(n);
  }

  int num_atom_all = 0;
  for (int i=0; i< (int) num_atom.size(); i++)
    num_atom_all += num_atom[i];

  return num_atom_all;
}

void FitPOD::read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces,
    int *atomtype, std::string file, std::vector<std::string> species)
{

  std::string filename = std::move(file);
  FILE *fp;
  if (comm->me == 0) {
    fp = utils::open_potential(filename,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ", filename, utils::getsyserror());
  }

  char line[MAXLINE] = {'\0'};
  char *ptr;
  int eof = 0;
  int cfi = 0;
  int nat = 0;
  int ns = species.size();

  // loop over all lines of this xyz file and extract training data

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fp);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // words = ptrs to all words in line
    // strip single and double quotes from words

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &) {
      // ignore
    }

    if (words.size() == 0) continue;

    ValueTokenizer text(utils::trim_comment(line),"\"' \t\n\r\f");
    if (text.contains("attice")) {

      // find the word containing "lattice"

      auto it = std::find_if(words.begin(), words.end(), [](const std::string& str) { return str.find("attice") != std::string::npos; });

      // get index of element from iterator

      int index = std::distance(words.begin(), it);

      if (words[index].find("=") != std::string::npos) {

        // lattice numbers start at index + 1

        for (int k = 0; k < 9; k++) {
          lattice[k + 9*cfi] = utils::numeric(FLERR,words[index+1+k],false,lmp);
        }
      } else {

        // lattice numbers start at index + 2

        for (int k = 0; k < 9; k++) {
          lattice[k + 9*cfi] = utils::numeric(FLERR,words[index+2+k],false,lmp);
        }
      }

      // find the word containing "energy"

      it = std::find_if(words.begin(), words.end(), [](const std::string& str) { return str.find("nergy") != std::string::npos; });

      // get index of element from iterator

      index = std::distance(words.begin(), it);

      if (words[index].find("=") != std::string::npos) {

        // energy is after "=" inside this string

        std::size_t found = words[index].find("=");
        energy[cfi] = utils::numeric(FLERR,words[index].substr(found+1),false,lmp);
      } else {

        // energy is at index + 2

        energy[cfi] = utils::numeric(FLERR,words[index+2],false,lmp);

      }

      // find the word containing "stress"

      it = std::find_if(words.begin(), words.end(), [](const std::string& str) { return str.find("tress") != std::string::npos; });

      // get index of element from iterator

      index = std::distance(words.begin(), it);

      if (words[index].find("=") != std::string::npos) {

        // stress numbers start at index + 1

        for (int k = 0; k < 9; k++) {
          stress[k + 9*cfi] = utils::numeric(FLERR,words[index+1+k],false,lmp);
        }
      } else {

        // lattice numbers start at index + 2

        for (int k = 0; k < 9; k++) {
          stress[k + 9*cfi] = utils::numeric(FLERR,words[index+2+k],false,lmp);
        }
      }

      cfi += 1;

    }

    // loop over atoms

    else if (words.size() > 1) {

      for (int ii = 0; ii < ns; ii++)
        if (species[ii] == words[0])
          atomtype[nat] = ii+1;

      for (int k = 0; k < 6; k++) {
        if (k <= 2) pos[k + 3*nat] = utils::numeric(FLERR,words[1+k],false,lmp);
        if (k > 2 ) forces[k-3 + 3*nat] = utils::numeric(FLERR,words[1+k],false,lmp);
      }
      nat += 1;
    }
  }
}

void FitPOD::get_data(datastruct &data, const std::vector<std::string>& species)
{
  get_exyz_files(data.data_files, data.data_path, data.file_extension);
  data.num_atom_sum = get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);
  data.num_config_sum = data.num_atom.size();
  size_t maxname = 9;
  for (const auto &fname : data.data_files) maxname = MAX(maxname,fname.size());
  maxname -= data.data_path.size()+1;
  const std::string sepline(maxname+46, '-');
  if (comm->me == 0)
    utils::logmesg(lmp, "{}\n {:^{}} | number of configurations | number of atoms\n{}\n",
                   sepline, "data file", maxname, sepline);
  int i = 0;
  for (const auto &fname : data.data_files) {
    std::string filename = fname.substr(data.data_path.size()+1);
    data.filenames.push_back(filename);
    if (comm->me == 0)
      utils::logmesg(lmp, " {:<{}} |        {:>10}        |    {:>8}\n",
                     filename, maxname, data.num_config[i], data.num_atom_each_file[i]);
    ++i;
  }
  if (comm->me == 0) {
    utils::logmesg(lmp, "{}\n", sepline);
    utils::logmesg(lmp, "number of files: {}\n", data.data_files.size());
    utils::logmesg(lmp, "number of configurations in all files: {}\n", data.num_config_sum);
    utils::logmesg(lmp, "number of atoms in all files: {}\n", data.num_atom_sum);
  }

  if (data.data_files.size() < 1) error->all(FLERR, "Cannot fit potential without data files");

  int n = data.num_config_sum;
  memory->create(data.lattice, 9*n, "fitpod:lattice");
  memory->create(data.stress, 9*n, "fitpod:stress");
  memory->create(data.energy, n, "fitpod:energy");
  n = data.num_atom_sum;
  memory->create(data.position, 3*n, "fitpod:position");
  memory->create(data.force, 3*n, "fitpod:force");
  memory->create(data.atomtype, n, "fitpod:atomtype");

  int nfiles = data.data_files.size(); // number of files
  int nconfigs = 0;
  int natoms = 0;
  for (int i=0; i<nfiles; i++) {
    read_exyz_file(&data.lattice[9*nconfigs], &data.stress[9*nconfigs], &data.energy[nconfigs],
        &data.position[3*natoms], &data.force[3*natoms], &data.atomtype[natoms],
        data.data_files[i], species);
    nconfigs += data.num_config[i];
    natoms += data.num_atom_each_file[i];
  }

  int len = data.num_atom.size();
  data.num_atom_min = podArrayMin(&data.num_atom[0], len);
  data.num_atom_max = podArrayMax(&data.num_atom[0], len);
  data.num_atom_cumsum.resize(len+1);
  podCumsum(&data.num_atom_cumsum[0], &data.num_atom[0], len+1);

  data.num_config_cumsum.resize(nfiles+1);
  podCumsum(&data.num_config_cumsum[0], &data.num_config[0], nfiles+1);

  // convert all structures to triclinic system

  constexpr int DIM = 3;
  double Qmat[DIM*DIM];
  for (int ci=0; ci<len; ci++) {
    int natom = data.num_atom[ci];
    int natom_cumsum = data.num_atom_cumsum[ci];
    double *x = &data.position[DIM*natom_cumsum];
    double *f = &data.force[DIM*natom_cumsum];
    double *lattice = &data.lattice[9*ci];
    double *a1 = &lattice[0];
    double *a2 = &lattice[3];
    double *a3 = &lattice[6];

    matrix33_inverse(Qmat, a1, a2, a3);
    triclinic_lattice_conversion(a1, a2, a3, a1, a2, a3);
    matrix33_multiplication(Qmat, lattice, Qmat, DIM);
    matrix33_multiplication(x, Qmat, x, natom);
    matrix33_multiplication(f, Qmat, f, natom);
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "minimum number of atoms: {}\n", data.num_atom_min);
    utils::logmesg(lmp, "maximum number of atoms: {}\n", data.num_atom_max);
  }
}

std::vector<int> FitPOD::linspace(int start_in, int end_in, int num_in)
{

  std::vector<int> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  int elm;

  if (num == 0) { return linspaced; }
  if (num == 1)
  {
    elm = (int) std::round(start);
    linspaced.push_back(elm);
    return linspaced;
  }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
  {
    elm = (int) std::round(start + delta * i);
    linspaced.push_back(elm);
  }

  elm = (int) std::round(end);
  linspaced.push_back(elm);

  return linspaced;
}

std::vector<int> FitPOD::shuffle(int start_in, int end_in, int num_in)
{
  int sz = end_in - start_in + 1;
  std::vector<int> myvector(sz);

  for (int i = 0; i<sz; i++)
    myvector[i] = start_in + i;

  unsigned seed = (unsigned) platform::walltime()*1.0e9;
  std::shuffle (myvector.begin(), myvector.end(), std::default_random_engine(seed));

  std::vector<int> shuffle_vec(num_in);
  for (int i = 0; i<num_in; i++)
    shuffle_vec[i] = myvector[i];

  return shuffle_vec;
}

std::vector<int> FitPOD::select(int n, double fraction, int randomize)
{
  std::vector<int> selected;

  int m = (int) std::round(n*fraction);
  m = MAX(m, 1);

  selected = (randomize==1) ? shuffle(1, n, m) : linspace(1, n, m);

  return selected;
}

void FitPOD::select_data(datastruct &newdata, const datastruct &data)
{
  double fraction = data.fraction;
  int randomize = data.randomize;

  if (comm->me == 0) {
    if (randomize==1)
      utils::logmesg(lmp, "Select {} fraction of the data set at random using shuffle\n", data.fraction);
    else
      utils::logmesg(lmp, "Select {} fraction of the data set deterministically using linspace\n", data.fraction);
  }

  int nfiles = data.data_files.size();  // number of files
  std::vector<std::vector<int>> selected(nfiles);

  newdata.num_config.resize(nfiles);
  newdata.num_config_cumsum.resize(nfiles+1);
  newdata.num_atom_each_file.resize(nfiles);

  for (int file = 0; file < nfiles; file++) {
    int nconfigs = data.num_config[file];
    selected[file] = select(nconfigs, fraction, randomize);
    int ns = (int) selected[file].size(); // number of selected configurations

    newdata.num_config[file] = ns;
    int num_atom_sum = 0;
    for (int ii=0; ii < ns; ii++) { // loop over each selected configuration in a file
      int ci =  data.num_config_cumsum[file] + selected[file][ii] - 1;
      int natom = data.num_atom[ci];
      newdata.num_atom.push_back(natom);
      num_atom_sum += natom;
    }
    newdata.num_atom_each_file[file] = num_atom_sum;
  }
  int len = newdata.num_atom.size();
  newdata.num_atom_min = podArrayMin(&newdata.num_atom[0], len);
  newdata.num_atom_max = podArrayMax(&newdata.num_atom[0], len);
  newdata.num_atom_cumsum.resize(len+1);
  podCumsum(&newdata.num_atom_cumsum[0], &newdata.num_atom[0], len+1);
  newdata.num_atom_sum = newdata.num_atom_cumsum[len];
  podCumsum(&newdata.num_config_cumsum[0], &newdata.num_config[0], nfiles+1);
  newdata.num_config_sum = newdata.num_atom.size();

  int n = data.num_config_sum;
  memory->create(newdata.lattice, 9*n, "fitpod:newdata_lattice");
  memory->create(newdata.stress, 9*n, "fitpod:newdata_stress");
  memory->create(newdata.energy, n, "fitpod:newdata_energy");
  n = data.num_atom_sum;
  memory->create(newdata.position, 3*n, "fitpod:newdata_position");
  memory->create(newdata.force, 3*n, "fitpod:newdata_force");
  memory->create(newdata.atomtype, n, "fitpod:newdata_atomtype");

  int cn = 0;
  int dim = 3;
  for (int file = 0; file < nfiles; file++) {
    int ns = (int) selected[file].size(); // number of selected configurations
    for (int ii=0; ii < ns; ii++) { // loop over each selected configuration in a file
      int ci =  data.num_config_cumsum[file] + selected[file][ii] - 1;
      int natom = data.num_atom[ci];
      int natom_cumsum = data.num_atom_cumsum[ci];

      int natomnew = newdata.num_atom[cn];
      int natomnew_cumsum = newdata.num_atom_cumsum[cn];

      if (natom != natomnew)
        error->all(FLERR,"number of atoms in the new data set must be the same as that in the old data set.");

      int *atomtype = &data.atomtype[natom_cumsum];
      double *position = &data.position[dim*natom_cumsum];
      double *force = &data.force[dim*natom_cumsum];

      newdata.energy[cn] = data.energy[ci];
      for (int j=0; j<9; j++) {
        newdata.stress[j+9*cn] = data.stress[j+9*ci];
        newdata.lattice[j+9*cn] = data.lattice[j+9*ci];
      }

      for (int na=0; na<natom; na++) {
        newdata.atomtype[na+natomnew_cumsum] = atomtype[na];
        for (int j=0; j<dim; j++) {
          newdata.position[j + 3*na + dim*natomnew_cumsum] = position[j + 3*na];
          newdata.force[j + 3*na + dim*natomnew_cumsum] = force[j + 3*na];
        }
      }
      cn += 1;
    }
  }

  data.copydatainfo(newdata);
  size_t maxname = 9;
  for (const auto &fname : data.data_files) maxname = MAX(maxname,fname.size());
  maxname -= data.data_path.size()+1;

  if (comm->me == 0)
    utils::logmesg(lmp, "{:-<{}}\n {:^{}} | # configs (selected) | # atoms (selected) "
                   "| # configs (original) | # atoms (original)\n{:-<{}}\n",
                   "", maxname+90, "data_file", maxname, "", maxname+90);
  for (int i=0; i< (int) newdata.data_files.size(); i++) {
    std::string filename = newdata.data_files[i].substr(newdata.data_path.size()+1,newdata.data_files[i].size());
    newdata.filenames.emplace_back(filename.c_str());
    if (comm->me == 0)
      utils::logmesg(lmp, " {:<{}} |       {:>8}       |      {:>8}      |       {:>8}       |     {:>8}\n",
                     newdata.filenames[i], maxname, newdata.num_config[i], newdata.num_atom_each_file[i],
                     data.num_config[i], data.num_atom_each_file[i]);
  }
  if (comm->me == 0) {
    utils::logmesg(lmp, "{:-<{}}\nnumber of files: {}\n", "", maxname+90, newdata.data_files.size());
    utils::logmesg(lmp, "number of configurations in all files (selected and original): {} and {}\n", newdata.num_config_sum, data.num_config_sum);
    utils::logmesg(lmp, "number of atoms in all files (selected and original: {} and {}\n", newdata.num_atom_sum, data.num_atom_sum);
  }
}

void FitPOD::read_data_files(const std::string& data_file, const std::vector<std::string>& species)
{
  datastruct data;

  // read data input file to datastruct

  data.precision = read_data_file(data.fitting_weights, data.file_format, data.file_extension,
                      testdata.data_path, data.data_path, data.filenametag, data_file);

  data.training_analysis = (int) data.fitting_weights[3];
  data.test_analysis = (int) data.fitting_weights[4];
  data.training_calculation = (int) data.fitting_weights[5];
  data.test_calculation = (int) data.fitting_weights[6];
  data.fraction = data.fitting_weights[7];
  data.randomize = (int) data.fitting_weights[9];

  data.copydatainfo(traindata);

  if (data.fraction >= 1.0) {
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** Begin of Training Data Set ****************\n");
    if (traindata.data_path.size() > 1)
      get_data(traindata, species);
    else
      error->all(FLERR,"data set is not found");
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");
  } else {
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** Begin of Training Data Set ****************\n");
    if (data.data_path.size() > 1)
      get_data(data, species);
    else
      error->all(FLERR,"data set is not found");
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");

    if (comm->me == 0)
      utils::logmesg(lmp, "**************** Begin of Select Training Data Set ****************\n");
    select_data(traindata, data);
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** End of Select Training Data Set ****************\n");

    memory->destroy(data.lattice);
    memory->destroy(data.energy);
    memory->destroy(data.stress);
    memory->destroy(data.position);
    memory->destroy(data.force);
    memory->destroy(data.atomtype);
  }

  if (((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path)) {
    testdata.training = 0;
    testdata.file_format = traindata.file_format;
    testdata.file_extension = traindata.file_extension;
    testdata.training_analysis = traindata.training_analysis;
    testdata.test_analysis = traindata.test_analysis;
    testdata.training_calculation = traindata.training_calculation;
    testdata.test_calculation = traindata.test_calculation;
    testdata.fraction = traindata.fitting_weights[8];
    testdata.randomize = (int) traindata.fitting_weights[10];
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** Begin of Test Data Set ****************\n");
    get_data(testdata, species);
    if (comm->me == 0)
      utils::logmesg(lmp, "**************** End of Test Data Set ****************\n");
  }
  else {
    testdata.data_path = traindata.data_path;
  }
}

int FitPOD::latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  int m=0, n=0, p=0;
  if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);
  if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);
  if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);

  // index for the center lattice

  int ind = m + (2*m+1)*(n) + (2*m+1)*(2*n+1)*(p);

  // number of lattices

  int nl = (2*m+1)*(2*n+1)*(2*p+1);

  for (int j=0; j<3*nx; j++)
    y[j] = x[j];
  int q = nx;

  for (int i = 0; i < (2*p+1); i++)
    for (int j = 0; j < (2*n+1); j++)
      for (int k = 0; k < (2*m+1); k++) {
        int ii = k + (2*m+1)*j + (2*m+1)*(2*n+1)*i;
        if (ii != ind) {
          double x0 = a1[0]*(k - m) + a2[0]*(j - n) + a3[0]*(i - p);
          double x1 = a1[1]*(k - m) + a2[1]*(j - n) + a3[1]*(i - p);
          double x2 = a1[2]*(k - m) + a2[2]*(j - n) + a3[2]*(i - p);
          for (int jj=0; jj<nx; jj++) {
            y[0+3*q] = x0 + x[0+3*jj];
            y[1+3*q] = x1 + x[1+3*jj];
            y[2+3*q] = x2 + x[2+3*jj];
            q = q + 1;
          }
        }
      }

  for (int i=0; i <nl; i++)
    for (int j=0; j<nx; j++)
      alist[j + nx*i] = j;

  return nl;
}

int FitPOD::podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
{
  int k = 0;
  for (int i = 0; i<nx; i++) {
    double *ri = &r[i*dim];
    int inc = 0;
    for (int j=0; j<N; j++) {
      double *rj = &r[dim*j];
      double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
      if  ((rijsq > SMALL) && (rijsq <= rcutsq))  {
        inc += 1;
        neighlist[k] = j;
        k += 1;
      }
    }
    numneigh[i] = inc;
  }
  return k;
}

int FitPOD::podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum,
    double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  double rcutsq = rcut*rcut;
  int dim = 3, nl = 0, nn = 0;

  // number of lattices

  nl = latticecoords(y, alist, x, a1, a2, a3, rcut, pbc, nx);
  int N = nx*nl;

  // total number of neighbors

  nn = podneighborlist(neighlist, numneigh, y, rcutsq, nx, N, dim);

  podCumsum(numneighsum, numneigh, nx+1);

  return nn;
}

void FitPOD::allocate_memory(const datastruct &data)
{
  int nd = podptr->pod.nd;
  memory->create(desc.gd, nd, "fitpod:desc_gd");
  memory->create(desc.A, nd*nd, "fitpod:desc_A");
  memory->create(desc.b, nd, "fitpod:desc_b");
  memory->create(desc.c, nd, "fitpod:desc_c");
  podArraySetValue(desc.A, 0.0, nd*nd);
  podArraySetValue(desc.b, 0.0, nd);
  podArraySetValue(desc.c, 0.0, nd);

  int dim = 3;
  int natom_max = data.num_atom_max;
  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nelements = podptr->pod.nelements;
  int nbesselpars = podptr->pod.nbesselpars;
  int nrbf2 = podptr->pod.nbf2;
  int nabf3 = podptr->pod.nabf3;
  int nrbf3 = podptr->pod.nrbf3;
  int *pdegree2 = podptr->pod.twobody;
  int *pdegree3 = podptr->pod.threebody;
  int *pbc = podptr->pod.pbc;
  double rcut = podptr->pod.rcut;

  int Nj=0, Nij=0;
  int m=0, n=0, p=0, nl=0, ny=0, na=0, np=0;

  for (int ci=0; ci<(int) data.num_atom.size(); ci++)
  {
    int natom = data.num_atom[ci];
    double *lattice = &data.lattice[9*ci];
    double *a1 = &lattice[0];
    double *a2 = &lattice[3];
    double *a3 = &lattice[6];
    if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);
    if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);
    if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);

    // number of lattices

    nl = (2*m+1)*(2*n+1)*(2*p+1);
    ny = MAX(ny,dim*natom*nl);
    na = MAX(na, natom*nl);
    np = MAX(np, natom*natom*nl);
  }

  memory->create(nb.y, ny, "fitpod:nb_y");
  memory->create(nb.alist, na, "fitpod:nb_alist");
  memory->create(nb.pairnum, natom_max, "fitpod:nb_pairnum");
  memory->create(nb.pairnum_cumsum, natom_max+1, "fitpod:nb_pairnum_cumsum");
  memory->create(nb.pairlist, np, "fitpod:nb_pairlist");

  nb.natom_max = natom_max;
  nb.sze = nelements*nelements;
  nb.sza = na;
  nb.szy = ny;
  nb.szp = np;

  if (comm->me == 0)
    utils::logmesg(lmp,"**************** Begin of Memory Allocation ****************\n");

  int szd = 0, szi=0, szsnap=0;
  for (int ci=0; ci<(int) data.num_atom.size(); ci++)
  {
    int natom = data.num_atom[ci];
    int natom_cumsum = data.num_atom_cumsum[ci];
    double *x = &data.position[dim*natom_cumsum];
    double *lattice = &data.lattice[9*ci];
    double *a1 = &lattice[0];
    double *a2 = &lattice[3];
    double *a3 = &lattice[6];

    Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, x, a1, a2, a3, rcut, pbc, natom);

    int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
    int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

    int szd1 = 3*Nij+ (1+dim)*Nij*MAX(nrbf2+ns2,nrbf3+ns3) + (nabf3+1)*7;
    int szi1 = 6*Nij + 2*natom+1 + (Nj-1)*Nj;
    szd = MAX(szd, szd1);
    szi = MAX(szi, szi1);

    if (podptr->sna.twojmax>0) {
      szd1 = 0;
      szd1 += Nij*dim; // rij
      szd1 += MAX(2*podptr->sna.idxu_max*Nij, 2*podptr->sna.idxz_max*podptr->sna.ndoubles*natom); // (Ur, Ui) and (Zr, Zi)
      szd1 += 2*podptr->sna.idxu_max*dim*Nij; // dUr, dUi
      szd1 += MAX(podptr->sna.idxb_max*podptr->sna.ntriples*dim*Nij, 2*podptr->sna.idxu_max*podptr->sna.nelements*natom); // dblist and (Utotr, Utoti)
      szsnap = MAX(szsnap, szd1);
    }
  }

  szd = MAX(szsnap, szd);
  szd = MAX(natom_max*(nd1+nd2+nd3+nd4) + szd, dim*natom_max*(nd-nd1-nd2-nd3-nd4));
  szd = dim*natom_max*(nd1+nd2+nd3+nd4) + szd;

  // gdd includes linear descriptors derivatives, quadratic descriptors derivatives and temporary memory

  memory->create(desc.gdd, szd, "fitpod:desc_gdd");
  memory->create(desc.tmpint, szi, "fitpod:desc_tmpint");
  desc.szd = szd;
  desc.szi = szi;

  if (comm->me == 0) {
    utils::logmesg(lmp, "maximum number of atoms in periodic domain: {}\n", natom_max);
    utils::logmesg(lmp, "maximum number of atoms in extended domain: {}\n", nb.sza);
    utils::logmesg(lmp, "maximum number of neighbors in extended domain: {}\n", nb.szp);
    utils::logmesg(lmp, "size of double memory: {}\n", szd);
    utils::logmesg(lmp, "size of int memory: {}\n", szi);
    utils::logmesg(lmp, "size of descriptor matrix: {} x {}\n", nd, nd);
    utils::logmesg(lmp, "**************** End of Memory Allocation ****************\n");
  }
}

void FitPOD::linear_descriptors(const datastruct &data, int ci)
{
  int dim = 3;
  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nd1234 = nd1+nd2+nd3+nd4;
  int *pbc = podptr->pod.pbc;
  double rcut = podptr->pod.rcut;

  int natom = data.num_atom[ci];
  int natom_cumsum = data.num_atom_cumsum[ci];
  int *atomtype = &data.atomtype[natom_cumsum];
  double *position = &data.position[dim*natom_cumsum];
  double *lattice = &data.lattice[9*ci];
  double *a1 = &lattice[0];
  double *a2 = &lattice[3];
  double *a3 = &lattice[6];

  // neighbor list
  int Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum,
        position, a1, a2, a3, rcut, pbc, natom);

  int *tmpint = &desc.tmpint[0];
  double *tmpmem = &desc.gdd[dim*natom*nd1234+natom*nd1234];
  podptr->linear_descriptors(desc.gd, desc.gdd, nb.y, tmpmem, atomtype, nb.alist,
      nb.pairlist, nb.pairnum, nb.pairnum_cumsum, tmpint, natom, Nij);

}

void FitPOD::quadratic_descriptors(const datastruct &data, int ci)
{
  int dim = 3;
  int natom = data.num_atom[ci];
  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nd22 = podptr->pod.nd22;
  int nd23 = podptr->pod.nd23;
  int nd24 = podptr->pod.nd24;
  int nd33 = podptr->pod.nd33;
  int nd34 = podptr->pod.nd34;
  int nd44 = podptr->pod.nd44;
  int nd123 = nd1+nd2+nd3;
  int nd1234 = nd1+nd2+nd3+nd4;

  double *fatom2 = &desc.gdd[dim*natom*(nd1)];
  double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];
  double *fatom4 = &desc.gdd[dim*natom*(nd123)];

  // global descriptors for four-body quadratic22 potential

  if (nd22 > 0) {
    int nq2 = podptr->pod.quadratic22[0]*podptr->pod.nc2;
    podptr->quadratic_descriptors(&desc.gd[nd1234], &desc.gdd[dim*natom*nd1234],
        &desc.gd[nd1], fatom2, nq2, dim*natom);
  }

  // global descriptors for four-body quadratic23 potential

  if (nd23 > 0) {
    int nq2 = podptr->pod.quadratic23[0]*podptr->pod.nc2;
    int nq3 = podptr->pod.quadratic23[1]*podptr->pod.nc3;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22], &desc.gdd[dim*natom*(nd1234+nd22)],
        &desc.gd[nd1], &desc.gd[nd1+nd2], fatom2, fatom3, nq2, nq3, dim*natom);
  }

  // global descriptors for five-body quadratic24 potential

  if (nd24 > 0) {
    int nq2 = podptr->pod.quadratic24[0]*podptr->pod.nc2;
    int nq4 = podptr->pod.quadratic24[1]*podptr->pod.nc4;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23], &desc.gdd[dim*natom*(nd1234+nd22+nd23)],
        &desc.gd[nd1], &desc.gd[nd1+nd2+nd3], fatom2, fatom4, nq2, nq4, dim*natom);
  }

  // global descriptors for five-body quadratic33 potential

  if (nd33 > 0) {
    int nq3 = podptr->pod.quadratic33[0]*podptr->pod.nc3;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24)],
        &desc.gd[nd1+nd2], fatom3, nq3, dim*natom);
  }

  // global descriptors for six-body quadratic34 potential

  if (nd34 > 0) {
    int nq3 = podptr->pod.quadratic34[0]*podptr->pod.nc3;
    int nq4 = podptr->pod.quadratic34[1]*podptr->pod.nc4;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24+nd33], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24+nd33)],
        &desc.gd[nd1+nd2], &desc.gd[nd1+nd2+nd3], fatom3, fatom4, nq3, nq4, dim*natom);
  }

  // global descriptors for seven-body quadratic44 potential

  if (nd44 > 0) {
    int nq4 = podptr->pod.quadratic44[0]*podptr->pod.nc4;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24+nd33+nd34], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24+nd33+nd34)],
        &desc.gd[nd1+nd2+nd3], fatom4, nq4, dim*natom);
  }

  // normalize quadratic descriptors

  for (int i=0; i<(nd22+nd23+nd24+nd33+nd34+nd44); i++)
    desc.gd[nd1234+i] = desc.gd[nd1234+i]/(natom);

  for (int i=0; i<dim*natom*(nd22+nd23+nd24+nd33+nd34+nd44); i++)
    desc.gdd[dim*natom*nd1234+i] = desc.gdd[dim*natom*nd1234+i]/(natom);
}

void FitPOD::cubic_descriptors(const datastruct &data, int ci)
{
  int dim = 3;
  int natom = data.num_atom[ci];
  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nd22 = podptr->pod.nd22;
  int nd23 = podptr->pod.nd23;
  int nd24 = podptr->pod.nd24;
  int nd33 = podptr->pod.nd33;
  int nd34 = podptr->pod.nd34;
  int nd44 = podptr->pod.nd44;
  int nd234 = podptr->pod.nd234;
  int nd333 = podptr->pod.nd333;
  int nd444 = podptr->pod.nd444;
  int nd123 = nd1+nd2+nd3;
  int nd1234 = nd1+nd2+nd3+nd4;

  // global descriptors for seven-body cubic234 potential
  if (nd234 > 0) {
    int nq2 = podptr->pod.cubic234[0]*podptr->pod.nc2;
    int nq3 = podptr->pod.cubic234[1]*podptr->pod.nc3;
    int nq4 = podptr->pod.cubic234[2]*podptr->pod.nc4;
    int np3 = nd1234+nd22+nd23+nd24+nd33+nd34+nd44;
    double *eatom2 = &desc.gd[nd1];
    double *eatom3 = &desc.gd[nd1+nd2];
    double *eatom4 = &desc.gd[nd123];
    double *fatom2 = &desc.gdd[dim*natom*(nd1)];
    double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];
    double *fatom4 = &desc.gdd[dim*natom*(nd123)];
    podptr->cubic_descriptors(&desc.gd[np3], &desc.gdd[dim*natom*np3],
        eatom2, eatom3, eatom4, fatom2, fatom3, fatom4, nq2, nq3, nq4, dim*natom);
  }

  // global descriptors for seven-body cubic333 potential

  if (nd333 > 0) {
    int nq3 = podptr->pod.cubic333[0]*podptr->pod.nc3;
    int np3 = nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234;
    double *eatom3 = &desc.gd[nd1+nd2];
    double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];
    podptr->cubic_descriptors(&desc.gd[np3], &desc.gdd[dim*natom*np3],
        eatom3, fatom3, nq3, dim*natom);
  }

  // global descriptors for ten-body cubic444 potential

  if (nd444 > 0) {
    int nq4 = podptr->pod.cubic444[0]*podptr->pod.nc4;
    int np4 = nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234+nd333;
    double *eatom4 = &desc.gd[nd123];
    double *fatom4 = &desc.gdd[dim*natom*(nd123)];
    podptr->cubic_descriptors(&desc.gd[np4], &desc.gdd[dim*natom*(np4)],
        eatom4, fatom4, nq4, dim*natom);
  }

  // normalize cubic descriptors
  int nd = podptr->pod.nd;
  for (int i=(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); i<nd; i++)
    desc.gd[i] = desc.gd[i]/(natom*natom);

  for (int i=dim*natom*(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); i<dim*natom*nd; i++)
    desc.gdd[i] = desc.gdd[i]/(natom*natom);
}

void FitPOD::least_squares_matrix(const datastruct &data, int ci)
{
  int dim = 3;
  int natom = data.num_atom[ci];
  int natom_cumsum = data.num_atom_cumsum[ci];
  int nd = podptr->pod.nd;
  int nforce = dim*natom;

  // compute energy weight and force weight

  double normconst = 1.0;
  if (data.normalizeenergy==1) normconst = 1.0/natom;
  double we = data.fitting_weights[0];
  double wf = data.fitting_weights[1];
  double we2 = (we*we)*(normconst*normconst);
  double wf2 = (wf*wf);

  // get energy and force from the training data set

  double energy = data.energy[ci];
  double *force = &data.force[dim*natom_cumsum];

  // least-square matrix for all descriptors: A = A + (we*we)*(gd^T * gd)

  podKron(desc.A, desc.gd, desc.gd, we2, nd, nd);

  // least-square matrix for all descriptors derivatives: A =  A + (wf*wf) * (gdd^T * gdd)

  char cht = 'T';
  char chn = 'N';
  double one = 1.0;
  int inc1 = 1;
  DGEMM(&cht, &chn, &nd, &nd, &nforce, &wf2, desc.gdd, &nforce, desc.gdd, &nforce, &one, desc.A, &nd);

  // least-square vector for all descriptors: b = b + (we*we*energy)*gd

  double wee = we2*energy;
  for (int i = 0; i< nd; i++)
    desc.b[i] += wee*desc.gd[i];

  // least-square vector for all descriptors derivatives: b = b + (wf*wf) * (gdd^T * f)

  DGEMV(&cht, &nforce, &nd, &wf2, desc.gdd, &nforce, force, &inc1, &one, desc.b, &inc1);

}

void FitPOD::least_squares_fit(const datastruct &data)
{
  if (comm->me == 0)
    utils::logmesg(lmp, "**************** Begin of Least-Squares Fitting ****************\n");

  // loop over each configuration in the training data set

  for (int ci=0; ci < (int) data.num_atom.size(); ci++) {

    if ((ci % 100)==0) {
      if (comm->me == 0)
        utils::logmesg(lmp, "Configuration: # {}\n", ci+1);
    }

    if ((ci % comm->nprocs) == comm->me) {

      // compute linear POD descriptors

      linear_descriptors(data, ci);

      // compute quadratic POD descriptors

      quadratic_descriptors(data, ci);

      // compute cubic POD descriptors

      cubic_descriptors(data, ci);

      // assemble the least-squares linear system

      least_squares_matrix(data, ci);
    }
  }

  int nd = podptr->pod.nd;

  MPI_Allreduce(MPI_IN_PLACE, desc.b, nd, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, desc.A, nd*nd, MPI_DOUBLE, MPI_SUM, world);

  if (comm->me == 0) {

    // symmetrize A

    for (int i = 0; i<nd; i++)
      for (int j = i; j<nd; j++) {
        double a1 = desc.A[i + nd*j];
        double a2 = desc.A[j + nd*i];
        desc.A[i + nd*j] = 0.5*(a1+a2);
        desc.A[j + nd*i] = 0.5*(a1+a2);
      }

    // scale A and b

    double maxb = 0.0;
    for (int i = 0; i<nd; i++)
      maxb = (maxb > fabs(desc.b[i])) ? maxb : fabs(desc.b[i]);

    maxb = 1.0/maxb;
    for (int i = 0; i<nd; i++)
      desc.b[i] = desc.b[i]*maxb;

    for (int i = 0; i<nd*nd; i++)
      desc.A[i] = desc.A[i]*maxb;

    double regularizing_parameter = data.fitting_weights[11];

    for (int i = 0; i<nd; i++) {
      desc.c[i] = desc.b[i];
      desc.A[i + nd*i] = desc.A[i + nd*i]*(1.0 + regularizing_parameter);
      if (desc.A[i + nd*i] < regularizing_parameter) desc.A[i + nd*i] = regularizing_parameter;
    }

    // solving the linear system A * c = b

    int nrhs=1, info;
    char chu = 'U';
    DPOSV(&chu, &nd, &nrhs, desc.A, &nd, desc.c, &nd, &info);
  }

  MPI_Bcast(desc.c, nd, MPI_DOUBLE, 0, world);

  if (comm->me == 0) {     // save coefficients into a text file
    std::string filename = data.filenametag + "_coefficients"  + ".pod";
    FILE *fp = fopen(filename.c_str(), "w");

    fmt::print(fp, "POD_coefficients: {}\n", nd);
    for (int count = 0; count < nd; count++) {
      fmt::print(fp, "{:<10.{}f}\n",  desc.c[count], data.precision);
    }
    fclose(fp);
    utils::logmesg(lmp, "**************** End of Least-Squares Fitting ****************\n");
  }
}

double FitPOD::energyforce_calculation(double *force, double *coeff, const datastruct &data, int ci)
{
  int dim = 3;
  int *pbc = podptr->pod.pbc;
  double rcut = podptr->pod.rcut;
  int nd1234 = podptr->pod.nd1 + podptr->pod.nd2 + podptr->pod.nd3 + podptr->pod.nd4;

  int natom = data.num_atom[ci];
  int natom_cumsum2 = data.num_atom_cumsum[ci];
  int *atomtype = &data.atomtype[natom_cumsum2];
  double *position = &data.position[dim*natom_cumsum2];
  double *lattice = &data.lattice[9*ci];
  double *a1 = &lattice[0];
  double *a2 = &lattice[3];
  double *a3 = &lattice[6];

  // neighbor list

  int Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum,
        position, a1, a2, a3, rcut, pbc, natom);

  double *tmpmem = &desc.gdd[0];
  int *tmpint = &desc.tmpint[0];
  double *rij = &tmpmem[0]; // 3*Nij
  int *ai = &tmpint[0];   // Nij
  int *aj = &tmpint[Nij];   // Nij
  int *ti = &tmpint[2*Nij]; // Nij
  int *tj = &tmpint[3*Nij]; // Nij
  int *idxi = &tmpint[4*Nij]; // Nij
  podptr->podNeighPairs(rij, nb.y, idxi, ai, aj, ti, tj, nb.pairnum_cumsum, atomtype, nb.pairlist, nb.alist, natom);

  double *effectivecoeff = &tmpmem[3*Nij]; // 3*Nij
  podArraySetValue(effectivecoeff, 0.0, nd1234);

  double energy = podptr->energyforce_calculation(force, coeff, effectivecoeff, desc.gd, rij,
    &tmpmem[3*Nij+nd1234], nb.pairnum_cumsum, atomtype, idxi, ai, aj, ti, tj, natom, Nij);

  return energy;
}

void FitPOD::print_analysis(const datastruct &data, double *outarray, double *errors)
{
  int nfiles = data.data_files.size();  // number of files
  int lm = 10;
  for (int i = 0; i < nfiles; i++)
    lm = MAX(lm, (int) data.filenames[i].size());
  lm = lm + 2;

  std::string filename_errors = fmt::format("{}_{}_errors.pod", data.filenametag, data.training ? "training" : "test");
  std::string filename_analysis = fmt::format("{}_{}_analysis.pod", data.filenametag, data.training ? "training" : "test");

  FILE *fp_errors = nullptr;
  FILE *fp_analysis = nullptr;
  fp_errors = fopen(filename_errors.c_str(), "w");
  fp_analysis = fopen(filename_analysis.c_str(), "w");

  std::string mystr = fmt::format("**************** Begin of Error Analysis for the {} Data Set ****************\n",
                                  data.training ? "Training" : "Test");

  utils::logmesg(lmp, mystr);
  fmt::print(fp_errors, mystr);

  std::string sa(lm+80,'-');
  sa += '\n';
  std::string sb = fmt::format(" {:^{}} | # configs |  # atoms  | MAE energy  | RMSE energy | MAE force  | RMSE force\n",
                               "File", lm);
  utils::logmesg(lmp, sa + sb + sa);
  fmt::print(fp_errors, sa + sb + sa);

  int ci=0, m=8, nc=0, nf=0;
  for (int file = 0; file < nfiles; file++) {
    fmt::print(fp_analysis, "# {}\n", data.filenames[file]);
    fmt::print(fp_analysis, "  config   # atoms      energy        DFT energy     energy error   "
               "  force          DFT force       force error\n");

    int nforceall = 0;
    int nconfigs = data.num_config[file];
    nc += nconfigs;
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      fmt::print(fp_analysis, "{:6}   {:8}    ", outarray[m*ci], outarray[1 + m*ci]);
      for(int count = 2; count < m; count ++)
        fmt::print(fp_analysis, "{:<15.10} ", outarray[count + m*ci]);
      fmt::print(fp_analysis, "\n");

      nforceall += 3*data.num_atom[ci];
      ci += 1;
    }
    nf += nforceall;

    int q = file+1;
    auto s = fmt::format("{:<{}} {:>10} {:>11}     {:<10.6f}    {:<10.6f}    {:<10.6f}    {:<10.6f}\n",
                         data.filenames[file], lm, nconfigs, nforceall/3,
                         errors[0 + 4*q], errors[1 + 4*q], errors[2 + 4*q], errors[3 + 4*q]);
    utils::logmesg(lmp, s);
    fmt::print(fp_errors, s);
  }
  utils::logmesg(lmp, sa);
  fmt::print(fp_errors, sa);

  auto s = fmt::format("{:<{}} {:>10} {:>11}     {:<10.6f}    {:<10.6f}    {:<10.6f}    {:<10.6f}\n",
                       "All files", lm, nc, nf/3, errors[0], errors[1], errors[2], errors[3]);
  utils::logmesg(lmp, s + sa);
  fmt::print(fp_errors, "{}", s + sa);

  mystr = fmt::format("**************** End of Error Analysis for the {} Data Set ****************\n",
                      data.training ? "Training" : "Test");

  utils::logmesg(lmp, mystr);
  fmt::print(fp_errors, mystr);

  fclose(fp_errors);
  fclose(fp_analysis);
}

void FitPOD::error_analysis(const datastruct &data, double *coeff)
{
  int dim = 3;
  double energy;
  std::vector<double> force(dim*data.num_atom_max);

  int nfiles = data.data_files.size();  // number of files
  int num_configs = data.num_atom.size(); // number of configurations in all files

  int m = 8;
  std::vector<double> outarray(m*num_configs);
  for (int i=0; i<m*num_configs; i++)
    outarray[i] = 0.0;

  std::vector<double> ssrarray(num_configs);
  for (int i=0; i<num_configs; i++)
    ssrarray[i] = 0.0;

  std::vector<double> errors(4*(nfiles+1));
  for (int i=0; i<4*(nfiles+1); i++)
    errors[i] = 0.0;

  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nd22 = podptr->pod.nd22;
  int nd23 = podptr->pod.nd23;
  int nd24 = podptr->pod.nd24;
  int nd33 = podptr->pod.nd33;
  int nd34 = podptr->pod.nd34;
  int nd44 = podptr->pod.nd44;
  int nd1234 = nd1+nd2+nd3+nd4;
  int nd = podptr->pod.nd;

  std::vector<double> newcoeff(nd);
  for (int j=0; j<nd; j++)
    newcoeff[j] = coeff[j];

  if (comm->me == 0)
    utils::logmesg(lmp, "**************** Begin of Error Calculation ****************\n");

  int ci = 0; // configuration counter
  for (int file = 0; file < nfiles; file++) { // loop over each file in the training data set

    int nconfigs = data.num_config[file];
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file

      if ((ci % 100)==0) {
        if (comm->me == 0)
          utils::logmesg(lmp, "Configuration: # {}\n", ci+1);
      }

      if ((ci % comm->nprocs) == comm->me) {
        int natom = data.num_atom[ci];
        int nforce = dim*natom;

        for (int j=nd1234; j<(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j++)
          newcoeff[j] = coeff[j]/(natom);

        for (int j=(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j<nd; j++)
          newcoeff[j] = coeff[j]/(natom*natom);

        energy = energyforce_calculation(force.data(), newcoeff.data(), data, ci);

        double DFTenergy = data.energy[ci];
        int natom_cumsum = data.num_atom_cumsum[ci];
        double *DFTforce = &data.force[dim*natom_cumsum];

        outarray[0 + m*ci] = ci+1;
        outarray[1 + m*ci] = natom;
        outarray[2 + m*ci] = energy;
        outarray[3 + m*ci] = DFTenergy;
        outarray[4 + m*ci] = fabs(DFTenergy-energy)/natom;
        outarray[5 + m*ci] = podArrayNorm(force.data(), nforce);
        outarray[6 + m*ci] = podArrayNorm(DFTforce, nforce);

        double diff, sum = 0.0, ssr = 0.0;
        for (int j=0; j<dim*natom; j++) {
          diff = DFTforce[j] - force[j];
          sum += fabs(diff);
          ssr += diff*diff;
        }
        outarray[7 + m*ci] = sum/nforce;
        ssrarray[ci] = ssr;
      }

      ci += 1;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &outarray[0], m*num_configs, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &ssrarray[0], num_configs, MPI_DOUBLE, MPI_SUM, world);

  ci = 0; // configuration counter
  int nc = 0, nf = 0;
  for (int file = 0; file < nfiles; file++) { // loop over each file in the training data set

    double emae=0.0, essr=0.0, fmae=0.0, fssr=0.0;
    int nforceall = 0;

    int nconfigs = data.num_config[file];
    nc += nconfigs;
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file

      int natom = data.num_atom[ci];
      int nforce = dim*natom;

      emae += outarray[4 + m*ci];
      essr += outarray[4 + m*ci]*outarray[4 + m*ci];
      fmae += outarray[7 + m*ci]*nforce;
      fssr += ssrarray[ci];
      nforceall += nforce;
      ci += 1;
    }

    int q = file + 1;
    if (nconfigs == 0) nconfigs = 1;
    if (nforceall == 0) nforceall = 1;
    errors[0 + 4*q] = emae/nconfigs;
    errors[1 + 4*q] = sqrt(essr/nconfigs);
    errors[2 + 4*q] = fmae/nforceall;
    errors[3 + 4*q] = sqrt(fssr/nforceall);

    nf += nforceall;
    errors[0] += emae;
    errors[1] += essr;
    errors[2] += fmae;
    errors[3] += fssr;
  }

  if (nc == 0) nc = 1;
  if (nf == 0) nf = 1;
  errors[0] = errors[0]/nc;
  errors[1] = sqrt(errors[1]/nc);
  errors[2] = errors[2]/nf;
  errors[3] = sqrt(errors[3]/nf);

  if (comm->me == 0) {
    utils::logmesg(lmp, "**************** End of Error Calculation ****************\n");
    print_analysis(data, outarray.data(), errors.data());
  }
}

void FitPOD::energyforce_calculation(const datastruct &data, double *coeff)
{
  int dim = 3;
  double energy;
  std::vector<double> force(1+dim*data.num_atom_max);

  int nfiles = data.data_files.size();  // number of files

  if (comm->me == 0)
    utils::logmesg(lmp, "**************** Begin of Energy/Force Calculation ****************\n");

  int ci = 0; // configuration counter
  for (int file = 0; file < nfiles; file++) { // loop over each file in the data set

    int nconfigs = data.num_config[file];
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      if ((ci % 100)==0) {
        if (comm->me == 0) utils::logmesg(lmp, "Configuration: # {}\n", ci+1);
      }

      int natom = data.num_atom[ci];
      int nforce = dim*natom;

      if ((ci % comm->nprocs) == comm->me) {
        energy = energyforce_calculation(force.data()+1, coeff, data, ci);

        // save energy and force into a binary file

        force[0] = energy;
        std::string filename = "energyforce_config" + std::to_string(ci+1) + ".bin";

        FILE *fp = fopen(filename.c_str(), "wb");

        fwrite( reinterpret_cast<char*>( force.data() ), sizeof(double) * (1 + nforce), 1, fp);

        fclose(fp);
      }
      ci += 1;
    }
  }
  if (comm->me == 0)
    utils::logmesg(lmp, "**************** End of Energy/Force Calculation ****************\n");
}

void FitPOD::print_matrix(const char *desc, int m, int n, double **a, int /*lda*/ )
{
  int i, j;
  printf( "\n %s\n", desc );

  for( i = 0; i < m; i++ )
  {
    for( j = 0; j < n; j++ ) printf( " %6.12f", a[j][i] );
    printf( "\n" );
  }
}

void FitPOD::print_matrix(const char *desc, int m, int n, double *a, int lda )
{
  int i, j;
  printf( "\n %s\n", desc );

  for( i = 0; i < m; i++ )
  {
    for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
    printf( "\n" );
  }
}

void FitPOD::print_matrix(const char *desc, int m, int n, int *a, int lda)
{
  int i, j;
  printf( "\n %s\n", desc );

  for( i = 0; i < m; i++ )
  {
    for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
    printf( "\n" );
  }
}

void FitPOD::podArrayFill(int* output, int start, int length)
{
        for (int j = 0; j < length; ++j)
                output[j] = start + j;
}

double FitPOD::podArraySum(double *a, int n)
{
  double e = a[0];
  for (int i=1; i<n; i++)
    e += a[i];
  return e;
}

double FitPOD::podArrayMin(double *a, int n)
{
  double b = a[0];
  for (int i=1; i<n; i++)
    if (a[i]<b)
      b = a[i];
  return b;
}

double FitPOD::podArrayMax(double *a, int n)
{
  double b = a[0];
  for (int i=1; i<n; i++)
    if (a[i]>b)
      b = a[i];
  return b;
}

int FitPOD::podArrayMin(int *a, int n)
{
  int b = a[0];
  for (int i=1; i<n; i++)
    if (a[i]<b)
      b = a[i];
  return b;
}

int FitPOD::podArrayMax(int *a, int n)
{
  int b = a[0];
  for (int i=1; i<n; i++)
    if (a[i]>b)
      b = a[i];
  return b;
}

void FitPOD::podKron(double *C, double *A, double *B, double alpha, int M1, int M2)
{
  int M = M1*M2;
  for (int idx=0; idx<M; idx++)
  {
    int ib = idx%M2;
    int ia = (idx-ib)/M2;
    C[idx] += alpha*A[ia]*B[ib];
  }
}

void FitPOD::podCumsum(int* output, int* input, int length)
{
        output[0] = 0;
        for (int j = 1; j < length; ++j)
                output[j] = input[j - 1] + output[j - 1];
}

double FitPOD::podArrayNorm(double *a, int n)
{
  double e = a[0]*a[0];
  for (int i=1; i<n; i++)
    e += a[i]*a[i];
  return sqrt(e);
}

double FitPOD::podArrayErrorNorm(double *a, double *b, int n)
{
  double e = (a[0]-b[0])*(a[0]-b[0]);
  for (int i=1; i<n; i++)
    e += (a[i]-b[i])*(a[i]-b[i]);
  return sqrt(e);
}

void FitPOD::podArraySetValue(double *y, double a, int n)
{
  for (int i=0; i<n; i++)
    y[i] = a;
}

void FitPOD::podArrayCopy(double *y, double *x, int n)
{
  for (int i=0; i<n; i++)
    y[i] = x[i];
}

void FitPOD::rotation_matrix(double *Rmat, double alpha, double beta, double gamma)
{
  double ca = cos(alpha);
  double cb = cos(beta);
  double cg = cos(gamma);
  double sa = sin(alpha);
  double sb = sin(beta);
  double sg = sin(gamma);

  Rmat[0] = ca*cg*cb - sa*sg;
  Rmat[3] = -ca*cb*sg - sa*cg;
  Rmat[6] = ca*sb;

  Rmat[1] = sa*cg*cb + ca*sg;
  Rmat[4] = -sa*cb*sg + ca*cg;
  Rmat[7] = sa*sb;

  Rmat[2] = -sb*cg;
  Rmat[5] = sb*sg;
  Rmat[8] = cb;
}

void FitPOD::matrix33_multiplication(double *xrot, double *Rmat, double *x, int natom)
{
  double x1, x2, x3;
  for (int i=0; i < natom; i++) {
    x1 = x[0 + 3*i];
    x2 = x[1 + 3*i];
    x3 = x[2 + 3*i];
    xrot[0 + 3*i] = Rmat[0]*x1 + Rmat[3]*x2 + Rmat[6]*x3;
    xrot[1 + 3*i] = Rmat[1]*x1 + Rmat[4]*x2 + Rmat[7]*x3;
    xrot[2 + 3*i] = Rmat[2]*x1 + Rmat[5]*x2 + Rmat[8]*x3;
  }
}

void FitPOD::matrix33_inverse(double *invA, double *A1, double *A2, double *A3)
{
  double a11 = A1[0];
  double a21 = A1[1];
  double a31 = A1[2];
  double a12 = A2[0];
  double a22 = A2[1];
  double a32 = A2[2];
  double a13 = A3[0];
  double a23 = A3[1];
  double a33 = A3[2];
  double detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);

  invA[0] = (a22*a33 - a23*a32)/detA;
  invA[1] = (a23*a31 - a21*a33)/detA;
  invA[2] = (a21*a32 - a22*a31)/detA;
  invA[3] = (a13*a32 - a12*a33)/detA;
  invA[4] = (a11*a33 - a13*a31)/detA;
  invA[5] = (a12*a31 - a11*a32)/detA;
  invA[6] = (a12*a23 - a13*a22)/detA;
  invA[7] = (a13*a21 - a11*a23)/detA;
  invA[8] = (a11*a22 - a12*a21)/detA;
}

void FitPOD::triclinic_lattice_conversion(double *a, double *b, double *c, double *A, double *B, double *C)
{
  double Anorm = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  double Bnorm = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  double Cnorm = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);

  double Ahat[3];
  Ahat[0] = A[0]/Anorm; Ahat[1] = A[1]/Anorm; Ahat[2] = A[2]/Anorm;

  double ax = Anorm;
  double bx = B[0]*Ahat[0] + B[1]*Ahat[1] + B[2]*Ahat[2]; //dot(B,Ahat);
  double by = sqrt(Bnorm*Bnorm - bx*bx); //sqrt(Bnorm^2 - bx^2);// #norm(cross(Ahat,B));
  double cx = C[0]*Ahat[0] + C[1]*Ahat[1] + C[2]*Ahat[2]; // dot(C,Ahat);
  double cy = (B[0]*C[0] + B[1]*C[1] + B[2]*C[2] - bx*cx)/by; // (dot(B, C) - bx*cx)/by;
  double cz = sqrt(Cnorm*Cnorm - cx*cx - cy*cy); // sqrt(Cnorm^2 - cx^2 - cy^2);

  a[0] = ax; a[1] = 0.0; a[2] = 0.0;
  b[0] = bx; b[1] = by;  b[2] = 0.0;
  c[0] = cx; c[1] = cy;  c[2] = cz;
}
