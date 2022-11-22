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

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "pod.h"
#include "tokenizer.h"
#include "update.h"
#include "math_const.h"

#include <algorithm>
#include <cmath>
#include <glob.h>
#include <random>
#include <string>
#include <vector>

using namespace LAMMPS_NS;

#define MAXLINE 1024

std::vector<std::string> static globVector(const std::string& pattern, std::vector<std::string> & files)
{
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    std::string s = std::string(glob_result.gl_pathv[i]);
    files.push_back(s);
  }
  globfree(&glob_result);
  return files;
}

void CFITPOD::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "fitpod", error);

  std::string pod_file = std::string(arg[0]);  // pod input file
  std::string data_file = std::string(arg[1]); // data input file
  std::string coeff_file; // coefficient input file

  if (narg > 2)
    coeff_file = std::string(arg[2]); // coefficient input file
  else
    coeff_file = "";

  podptr = new CPOD(lmp, pod_file, coeff_file);
  read_data_files(data_file, podptr->pod.species);

  if ((int) traindata.data_path.size() > 1)
    allocate_memory(traindata);
  else if ((int) testdata.data_path.size() > 1)
    allocate_memory(testdata);

  // get POD coefficients from an input file

  if (coeff_file != "") podptr->podArrayCopy(desc.c, podptr->pod.coeff, podptr->pod.nd);

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
}

/* ---------------------------------------------------------------------- */

void CFITPOD::read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension,
    std::string &test_path, std::string &training_path, std::string data_file)
{
  std::string datafilename = data_file;
  FILE *fpdata;
  if (comm->me == 0){

    fpdata = utils::open_potential(datafilename,lmp,nullptr);
    if (fpdata == nullptr)
      error->one(FLERR,"Cannot open training data file {}: ",
                                   datafilename, utils::getsyserror());
  }

  // loop through lines of training data file and parse keywords

  char line[MAXLINE],*ptr;
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

    auto keywd = words[0];

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
    if (keywd == "percentage_training_data_set") fitting_weights[7] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "percentage_test_data_set") fitting_weights[8] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "randomize_training_data_set") fitting_weights[9] = utils::numeric(FLERR,words[1],false,lmp);
    if (keywd == "randomize_test_data_set") fitting_weights[10] = utils::numeric(FLERR,words[1],false,lmp);

    // other settings

    if (keywd == "file_format") file_format = words[1];
    if (keywd == "file_extension") file_extension = words[1];
    if (keywd == "path_to_training_data_set") training_path = words[1];
    if (keywd == "path_to_test_data_set") test_path = words[1];
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "**************** Begin of Data File ****************\n");
    utils::logmesg(lmp, "file format: {}\n", file_format);
    utils::logmesg(lmp, "file extension: {}\n", file_extension);
    utils::logmesg(lmp, "path to training data set: {}\n", training_path);
    utils::logmesg(lmp, "path to test data set: {}\n", test_path);
    utils::logmesg(lmp, "training percentage: {}\n", fitting_weights[7]);
    utils::logmesg(lmp, "test percentage: {}\n", fitting_weights[8]);
    utils::logmesg(lmp, "randomize training data set: {}\n", fitting_weights[9]);
    utils::logmesg(lmp, "randomize test data set: {}\n", fitting_weights[10]);
    utils::logmesg(lmp, "error analysis for training data set: {}\n", fitting_weights[3]);
    utils::logmesg(lmp, "error analysis for test data set: {}\n", fitting_weights[4]);
    utils::logmesg(lmp, "energy/force calculation for training data set: {}\n", fitting_weights[5]);
    utils::logmesg(lmp, "energy/force calculation for test data set: {}\n", fitting_weights[6]);
    utils::logmesg(lmp, "fitting weight for energy: {}\n", fitting_weights[0]);
    utils::logmesg(lmp, "fitting weight for force: {}\n", fitting_weights[1]);
    utils::logmesg(lmp, "fitting weight for stress: {}\n", fitting_weights[2]);
    utils::logmesg(lmp, "**************** End of Data File ****************\n");
  }
}

void CFITPOD::get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension)
{
  std::vector<std::string> res = globVector(datapath + "/*." + extension, files);
}

int CFITPOD::get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)
{
  std::string filename = file;
  FILE *fp;
  if (comm->me == 0){

    fp = utils::open_potential(filename,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ", filename, utils::getsyserror());
  }

  char line[MAXLINE],*ptr;
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
    if (words.size() == 1){
      natom = utils::inumeric(FLERR,words[0],false,lmp);
      num_atom.push_back(natom);
      num_configs += 1;
      num_atom_sum += natom;
    }
  }
  return num_configs;
}

int CFITPOD::get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)
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

void CFITPOD::read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces,
    int *atomtype, std::string file, std::vector<std::string> species)
{

  std::string filename = file;
  FILE *fp;
  if (comm->me == 0){

    fp = utils::open_potential(filename,lmp,nullptr);
    if (fp == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ",
                                   filename, utils::getsyserror());
  }

  char line[MAXLINE],*ptr;
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
    if (text.contains("attice")){

      // find the word containing "lattice"

      auto it = std::find_if(words.begin(), words.end(), [](const std::string& str) { return str.find("attice") != std::string::npos; });

      // get index of element from iterator

      int index = std::distance(words.begin(), it);

      if (words[index].find("=") != std::string::npos) {

        // lattice numbers start at index + 1

        for (int k = 0; k < 9; k++){
          lattice[k + 9*cfi] = utils::numeric(FLERR,words[index+1+k],false,lmp);
        }
      } else {

        // lattice numbers start at index + 2

        for (int k = 0; k < 9; k++){
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

        // lattice numbers start at index + 1

        for (int k = 0; k < 9; k++){
          stress[k + 9*cfi] = utils::numeric(FLERR,words[index+1+k],false,lmp);
        }
      } else {

        // lattice numbers start at index + 2

        for (int k = 0; k < 9; k++){
          stress[k + 9*cfi] = utils::numeric(FLERR,words[index+2+k],false,lmp);
        }
      }

      cfi += 1;

    }

    // loop over atoms

    else if (words.size() > 1){

      for (int ii = 0; ii < ns; ii++)
        if (species[ii] == words[0])
          atomtype[nat] = ii+1;

      for (int k = 0; k < 6; k++){
        if (k <= 2) pos[k + 3*nat] = utils::numeric(FLERR,words[1+k],false,lmp);
        if (k > 2 ) forces[k-3 + 3*nat] = utils::numeric(FLERR,words[1+k],false,lmp);
      }

      nat += 1;
    }
  }

}

void CFITPOD::get_data(datastruct &data, std::vector<std::string> species)
{
  get_exyz_files(data.data_files, data.data_path, data.file_extension);
  data.num_atom_sum = get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);
  data.num_config_sum = data.num_atom.size();
  size_t maxname = 9;
  for (auto fname : data.data_files) maxname = MAX(maxname,fname.size());
  maxname -= data.data_path.size()+1;
  utils::logmesg(lmp, " {:^{}} | number of configurations | number of atoms\n{:=<{}}\n", "data file",
                 maxname, "", maxname+46);
  int i = 0;
  for (auto fname : data.data_files) {
    std::string filename = fname.substr(data.data_path.size()+1);
    data.filenames.push_back(filename);
    utils::logmesg(lmp, " {:<{}} |        {:>10}        |    {:>8}\n",
                   filename, maxname, data.num_config[i], data.num_atom_each_file[i]);
    ++i;
  }
  utils::logmesg(lmp, "{:=<{}}\n", "", maxname+46);
  utils::logmesg(lmp, "number of files: {}\n", data.data_files.size());
  utils::logmesg(lmp, "number of configurations in all files: {}\n", data.num_config_sum);
  utils::logmesg(lmp, "number of atoms in all files: {}\n", data.num_atom_sum);

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
  data.num_atom_min = podptr->podArrayMin(&data.num_atom[0], len);
  data.num_atom_max = podptr->podArrayMax(&data.num_atom[0], len);
  data.num_atom_cumsum.resize(len+1);
  podptr->podCumsum(&data.num_atom_cumsum[0], &data.num_atom[0], len+1);

  data.num_config_cumsum.resize(nfiles+1);
  podptr->podCumsum(&data.num_config_cumsum[0], &data.num_config[0], nfiles+1);

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

    podptr->matrix33_inverse(Qmat, a1, a2, a3);
    podptr->triclinic_lattice_conversion(a1, a2, a3, a1, a2, a3);
    podptr->matrix33_multiplication(Qmat, lattice, Qmat, DIM);
    podptr->matrix33_multiplication(x, Qmat, x, natom);
    podptr->matrix33_multiplication(f, Qmat, f, natom);
  }

  utils::logmesg(lmp, "minimum number of atoms: {}\n", data.num_atom_min);
  utils::logmesg(lmp, "maximum number of atoms: {}\n", data.num_atom_max);
}

std::vector<int> CFITPOD::linspace(int start_in, int end_in, int num_in)
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

std::vector<int> CFITPOD::shuffle(int start_in, int end_in, int num_in)
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

std::vector<int> CFITPOD::select(int n, double percentage, int randomize)
{
  std::vector<int> selected;

  int m = (int) std::round(n*percentage);
  m = MAX(m, 1);

  selected = (randomize==1) ? shuffle(1, n, m) : linspace(1, n, m);

  return selected;
}

void CFITPOD::select_data(datastruct &newdata, datastruct data)
{
  double percentage = data.percentage;
  int randomize = data.randomize;

  if (randomize==1)
    utils::logmesg(lmp, "Select {} percent of the data set at random using shuffle\n", data.percentage*100);
  else
    utils::logmesg(lmp, "Select {} percent of the data set deterministically using linspace\n", data.percentage*100);

  int nfiles = data.data_files.size();  // number of files
  std::vector<std::vector<int>> selected(nfiles);

  newdata.num_config.resize(nfiles);
  newdata.num_config_cumsum.resize(nfiles+1);
  newdata.num_atom_each_file.resize(nfiles);

  for (int file = 0; file < nfiles; file++) {
    int nconfigs = data.num_config[file];
    selected[file] = select(nconfigs, percentage, randomize);
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
  newdata.num_atom_min = podptr->podArrayMin(&newdata.num_atom[0], len);
  newdata.num_atom_max = podptr->podArrayMax(&newdata.num_atom[0], len);
  newdata.num_atom_cumsum.resize(len+1);
  podptr->podCumsum(&newdata.num_atom_cumsum[0], &newdata.num_atom[0], len+1);
  newdata.num_atom_sum = newdata.num_atom_cumsum[len];
  podptr->podCumsum(&newdata.num_config_cumsum[0], &newdata.num_config[0], nfiles+1);
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

  utils::logmesg(lmp, "data file  | # configs (selected) | # atoms (selected) | # configs (original) | # atoms (original)\n");
  for (int i=0; i< (int) newdata.data_files.size(); i++) {
    std::string filename = newdata.data_files[i].substr(newdata.data_path.size()+1,newdata.data_files[i].size());
    newdata.filenames.push_back(filename.c_str());
    utils::logmesg(lmp, "{}   |   {}   |   {}   |   {}   |   {}\n", newdata.filenames[i], newdata.num_config[i], newdata.num_atom_each_file[i], data.num_config[i], data.num_atom_each_file[i]);
  }
  utils::logmesg(lmp, "number of files: {}\n", newdata.data_files.size());
  utils::logmesg(lmp, "number of configurations in all files (selected and original): {} and {}\n", newdata.num_config_sum, data.num_config_sum);
  utils::logmesg(lmp, "number of atoms in all files (selected and original: {} and {}\n", newdata.num_atom_sum, data.num_atom_sum);
}

void CFITPOD::read_data_files(std::string data_file, std::vector<std::string> species)
{
  datastruct data;

  // read data input file to datastruct

  read_data_file(data.fitting_weights, data.file_format, data.file_extension,
      testdata.data_path, data.data_path, data_file);

  data.training_analysis = (int) data.fitting_weights[3];
  data.test_analysis = (int) data.fitting_weights[4];
  data.training_calculation = (int) data.fitting_weights[5];
  data.test_calculation = (int) data.fitting_weights[6];
  data.percentage = data.fitting_weights[7];
  data.randomize = (int) data.fitting_weights[9];

  data.copydatainfo(traindata);

  if (data.percentage >= 1.0) {
    utils::logmesg(lmp, "**************** Begin of Training Data Set ****************\n");
    if ((int) traindata.data_path.size() > 1)
      get_data(traindata, species);
    else
      error->all(FLERR,"data set is not found");
    utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");
  }
  else {
    utils::logmesg(lmp, "**************** Begin of Training Data Set ****************\n");
    if ((int) data.data_path.size() > 1)
      get_data(data, species);
    else
      error->all(FLERR,"data set is not found");
    utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");

    utils::logmesg(lmp, "**************** Begin of Select Training Data Set ****************\n");
    select_data(traindata, data);
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
    testdata.percentage = traindata.fitting_weights[8];
    testdata.randomize = (int) traindata.fitting_weights[10];
    utils::logmesg(lmp, "**************** Begin of Test Data Set ****************\n");
    get_data(testdata, species);
    utils::logmesg(lmp, "**************** End of Test Data Set ****************\n");
  }
  else {
    testdata.data_path = traindata.data_path;
  }
}

int CFITPOD::latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
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

int CFITPOD::podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
{
  int k = 0;
  for (int i = 0; i<nx; i++) {
    double *ri = &r[i*dim];
    int inc = 0;
    for (int j=0; j<N; j++) {
      double *rj = &r[dim*j];
      double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
      if  ((rijsq > 1e-12) && (rijsq <= rcutsq))  {
        inc += 1;
        neighlist[k] = j;
        k += 1;
      }
    }
    numneigh[i] = inc;
  }
  return k;
}

int CFITPOD::podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum,
    double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  double rcutsq = rcut*rcut;
  int dim = 3, nl = 0, nn = 0;

  // number of lattices

  nl = latticecoords(y, alist, x, a1, a2, a3, rcut, pbc, nx);
  int N = nx*nl;

  // total number of neighbors

  nn = podneighborlist(neighlist, numneigh, y, rcutsq, nx, N, dim);

  podptr->podCumsum(numneighsum, numneigh, nx+1);

  return nn;
}

void CFITPOD::allocate_memory(datastruct data)
{
  int nd = podptr->pod.nd;
  memory->create(desc.gd, nd, "fitpod:desc_gd");
  memory->create(desc.A, nd*nd, "fitpod:desc_A");
  memory->create(desc.b, nd, "fitpod:desc_b");
  memory->create(desc.c, nd, "fitpod:desc_c");
  podptr->podArraySetValue(desc.A, 0.0, nd*nd);
  podptr->podArraySetValue(desc.b, 0.0, nd);
  podptr->podArraySetValue(desc.c, 0.0, nd);

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

  utils::logmesg(lmp, "maximum number of atoms in periodic domain: {}\n", natom_max);
  utils::logmesg(lmp, "maximum number of atoms in extended domain: {}\n", nb.sza);
  utils::logmesg(lmp, "maximum number of neighbors in extended domain: {}\n", nb.szp);
  utils::logmesg(lmp, "size of double memory: {}\n", szd);
  utils::logmesg(lmp, "size of int memory: {}\n", szi);
  utils::logmesg(lmp, "size of descriptor matrix: {} x {}\n", nd, nd);
  utils::logmesg(lmp, "**************** End of Memory Allocation ****************\n");
}

void CFITPOD::linear_descriptors(datastruct data, int ci)
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

void CFITPOD::quadratic_descriptors(datastruct data, int ci)
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

  if (nd22>0) {
    int nq2 = podptr->pod.quadratic22[0]*podptr->pod.nc2;
    podptr->quadratic_descriptors(&desc.gd[nd1234], &desc.gdd[dim*natom*nd1234],
        &desc.gd[nd1], fatom2, nq2, dim*natom);
  }

  // global descriptors for four-body quadratic23 potential

  if (nd23>0) {
    int nq2 = podptr->pod.quadratic23[0]*podptr->pod.nc2;
    int nq3 = podptr->pod.quadratic23[1]*podptr->pod.nc3;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22], &desc.gdd[dim*natom*(nd1234+nd22)],
        &desc.gd[nd1], &desc.gd[nd1+nd2], fatom2, fatom3, nq2, nq3, dim*natom);
  }

  // global descriptors for five-body quadratic24 potential

  if (nd24>0) {
    int nq2 = podptr->pod.quadratic24[0]*podptr->pod.nc2;
    int nq4 = podptr->pod.quadratic24[1]*podptr->pod.nc4;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23], &desc.gdd[dim*natom*(nd1234+nd22+nd23)],
        &desc.gd[nd1], &desc.gd[nd1+nd2+nd3], fatom2, fatom4, nq2, nq4, dim*natom);
  }

  // global descriptors for five-body quadratic33 potential

  if (nd33>0) {
    int nq3 = podptr->pod.quadratic33[0]*podptr->pod.nc3;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24)],
        &desc.gd[nd1+nd2], fatom3, nq3, dim*natom);
  }

  // global descriptors for six-body quadratic34 potential

  if (nd34>0) {
    int nq3 = podptr->pod.quadratic34[0]*podptr->pod.nc3;
    int nq4 = podptr->pod.quadratic34[1]*podptr->pod.nc4;
    podptr->quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24+nd33], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24+nd33)],
        &desc.gd[nd1+nd2], &desc.gd[nd1+nd2+nd3], fatom3, fatom4, nq3, nq4, dim*natom);
  }

  // global descriptors for seven-body quadratic44 potential

  if (nd44>0) {
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

void CFITPOD::cubic_descriptors(datastruct data, int ci)
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
  if (nd234>0) {
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

  if (nd333>0) {
    int nq3 = podptr->pod.cubic333[0]*podptr->pod.nc3;
    int np3 = nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234;
    double *eatom3 = &desc.gd[nd1+nd2];
    double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];
    podptr->cubic_descriptors(&desc.gd[np3], &desc.gdd[dim*natom*np3],
        eatom3, fatom3, nq3, dim*natom);
  }

  // global descriptors for ten-body cubic444 potential

  if (nd444>0) {
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

void CFITPOD::least_squares_matrix(datastruct data, int ci)
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

  podptr->podKron(desc.A, desc.gd, desc.gd, we2, nd, nd);

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

void CFITPOD::least_squares_fit(datastruct data)
{
  utils::logmesg(lmp, "**************** Begin of Least-Squares Fitting ****************\n");

  // loop over each configuration in the training data set

  for (int ci=0; ci < (int) data.num_atom.size(); ci++) {
    if ((ci % 100)==0) utils::logmesg(lmp, "Configuration: # {}\n", ci+1);

    // compute linear POD descriptors

    linear_descriptors(data, ci);

    // compute quadratic POD descriptors

    quadratic_descriptors(data, ci);

    // compute cubic POD descriptors

    cubic_descriptors(data, ci);

    // assemble the least-squares linear system

    least_squares_matrix(data, ci);
  }

  int nd = podptr->pod.nd;

  for (int i = 0; i<nd; i++) {
    desc.c[i] = desc.b[i];
    desc.A[i + nd*i] = desc.A[i + nd*i]*(1.0 + 1e-12);
  }

  // solving the linear system A * c = b

  int nrhs=1, info;
  char chu = 'U';
  DPOSV(&chu, &nd, &nrhs, desc.A, &nd, desc.c, &nd, &info);

  podptr->print_matrix( "Least-squares coefficient vector:", 1, nd, desc.c, 1);

  // save coefficients into a text file

  FILE *fp = fopen("coefficients.txt", "w");

  fmt::print(fp, "POD_coefficients: {}\n", nd);
  for (int count = 0; count < nd; count++){

    fmt::print(fp, "{:.20}\n", desc.c[count]);
  }

  fclose(fp);

  utils::logmesg(lmp, "**************** End of Least-Squares Fitting ****************\n");
}

double CFITPOD::energyforce_calculation(double *force, double *coeff, datastruct data, int ci)
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
  podptr->podArraySetValue(effectivecoeff, 0.0, nd1234);

  double energy = podptr->energyforce_calculation(force, coeff, effectivecoeff, desc.gd, rij,
    &tmpmem[3*Nij+nd1234], nb.pairnum_cumsum, atomtype, idxi, ai, aj, ti, tj, natom, Nij);

  return energy;
}

void CFITPOD::print_analysis(datastruct data, double *outarray, double *errors)
{
  std::string s = "All files";
  int nfiles = data.data_files.size();  // number of files
  int lm = s.size();
  for (int i = 0; i < nfiles; i++)
    lm = MAX(lm, (int) data.filenames[i].size());
  lm = lm + 2;

  std::string filename_errors = data.training ? "training_errors.txt" : "test_errors.txt";
  std::string filename_analysis = data.training ? "training_analysis.txt" : "test_analysis.txt";

  FILE *fp_errors = fopen(filename_errors.c_str(), "w");
  FILE *fp_analysis = fopen(filename_errors.c_str(), "w");

  std::string sa = "**************** Begin of Error Analysis for the Training Data Set ****************";
  std::string sb = "**************** Begin of Error Analysis for the Test Data Set ****************";
  std::string mystr = (data.training) ? sa : sb;

  utils::logmesg(lmp, "{}\n", mystr);
  fmt::print(fp_errors, mystr + "\n");

  sa = "----------------------------------------------------------------------------------------\n";
  sb = "  File    | # configs | # atoms  | MAE energy | RMSE energy | MAE force | RMSE force |\n";
  utils::logmesg(lmp, "{}", sa);
  utils::logmesg(lmp, "{}", sb);
  utils::logmesg(lmp, "{}", sa);
  fmt::print(fp_errors, sa);
  fmt::print(fp_errors, sb);
  fmt::print(fp_errors, sa);

  int ci=0, m=8, nc=0, nf=0;
  for (int file = 0; file < nfiles; file++) {
    fmt::print(fp_analysis, "# {}\n", data.filenames[file]);
    sb = "|  config  |  # atoms  |  energy  | DFT energy | energy error |  force  | DFT force | force error |\n";
    fmt::print(fp_analysis, sb);

    int nforceall = 0;
    int nconfigs = data.num_config[file];
    nc += nconfigs;
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      fmt::print(fp_analysis, "   ");
      for(int count = 0; count < m; count ++)
        fmt::print(fp_analysis, "{}   ", outarray[count + m*ci]);
      fmt::print(fp_analysis, "\n");

      nforceall += 3*data.num_atom[ci];
      ci += 1;
    }
    nf += nforceall;

    int q = file+1;
    std::string s = data.filenames[file];
    s = s + std::string(lm-s.size(), ' ');
    std::string s1 = std::to_string(nconfigs);
    s1 = s1 + std::string(MAX(6- (int) s1.size(),1), ' ');
    s = s + "   " + s1;
    s1 = std::to_string(nforceall/3);
    s1 = s1 + std::string(MAX(7 - (int) s1.size(),1), ' ');
    s = s + "   " + s1;
    s1 = std::to_string(errors[0 + 4*q]);
    s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
    s = s + "   " + s1;
    s1 = std::to_string(errors[1 + 4*q]);
    s1 = s1 + std::string(MAX(10 - (int)  s1.size(),1), ' ');
    s = s + "   " + s1;
    s1 = std::to_string(errors[2 + 4*q]);
    s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
    s = s + "   " + s1;
    s1 = std::to_string(errors[3 + 4*q]);
    s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
    s = s + "   " + s1 + "\n";
    utils::logmesg(lmp, "{}", s);
    fmt::print(fp_errors, "{}", s);
  }
  utils::logmesg(lmp, "{}", sa);
  fmt::print(fp_errors, "{}", sa);

  s = s + std::string(MAX(lm - (int) s.size(),1), ' ');
  std::string s1 = std::to_string(nc);
  s1 = s1 + std::string(MAX(6- (int) s1.size(),1), ' ');
  s = s + "   " + s1;
  s1 = std::to_string(nf/3);
  s1 = s1 + std::string(MAX(7 - (int) s1.size(),1), ' ');
  s = s + "   " + s1;
  s1 = std::to_string(errors[0]);
  s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
  s = s + "   " + s1;
  s1 = std::to_string(errors[1]);
  s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
  s = s + "   " + s1;
  s1 = std::to_string(errors[2]);
  s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
  s = s + "   " + s1;
  s1 = std::to_string(errors[3]);
  s1 = s1 + std::string(MAX(10 - (int) s1.size(),1), ' ');
  s = s + "   " + s1 + "\n";
  utils::logmesg(lmp, "{}", s);
  utils::logmesg(lmp, "{}", sa);
  fmt::print(fp_errors, "{}", s);
  fmt::print(fp_errors, "{}", sa);

  sa = "**************** End of Error Analysis for the Training Data Set ****************";
  sb = "**************** End of Error Analysis for the Test Data Set ****************";
  mystr = (data.training) ? sa : sb;
  utils::logmesg(lmp, "{}\n", mystr);
  fmt::print(fp_errors, "{}\n", mystr);

  fclose(fp_errors);
  fclose(fp_analysis);

}

void CFITPOD::error_analysis(datastruct data, double *coeff)
{
  int dim = 3;
  double energy;
  double force[dim*data.num_atom_max];

  int nfiles = data.data_files.size();  // number of files
  int num_configs = data.num_atom.size(); // number of configurations in all files
  //int nd12 = podptr->pod.nd1 + podptr->pod.nd2;
  //int nd123 = podptr->pod.nd1 + podptr->pod.nd2 + podptr->pod.nd3;
  //int nd1234 = podptr->pod.nd1 + podptr->pod.nd2 + podptr->pod.nd3 + podptr->pod.nd4;
  //double effectivecoeff[nd1234];

  int m = 8;
  double outarray[m*num_configs];
  double errors[4*(nfiles+1)];
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

  double newcoeff[nd];
  for (int j=0; j<nd; j++)
    newcoeff[j] = coeff[j];

  utils::logmesg(lmp, "**************** Begin of Error Calculation ****************\n");

  int ci = 0; // configuration counter
  int nc = 0, nf = 0;
  for (int file = 0; file < nfiles; file++) { // loop over each file in the training data set
    double emae=0.0, essr=0.0, fmae=0.0, fssr=0.0;
    int nforceall = 0;

    int nconfigs = data.num_config[file];
    nc += nconfigs;
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      if ((ci % 100)==0) utils::logmesg(lmp, "Configuration: # {}\n", ci+1);

      int natom = data.num_atom[ci];
      int nforce = dim*natom;

      for (int j=nd1234; j<(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j++)
        newcoeff[j] = coeff[j]/(natom);

      for (int j=(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j<nd; j++)
        newcoeff[j] = coeff[j]/(natom*natom);

      energy = energyforce_calculation(force, newcoeff, data, ci);

      double DFTenergy = data.energy[ci];
      int natom_cumsum = data.num_atom_cumsum[ci];
      double *DFTforce = &data.force[dim*natom_cumsum];

      outarray[0 + m*ci] = ci+1;
      outarray[1 + m*ci] = natom;
      outarray[2 + m*ci] = energy;
      outarray[3 + m*ci] = DFTenergy;
      outarray[4 + m*ci] = fabs(DFTenergy-energy)/natom;
      outarray[5 + m*ci] = podptr->podArrayNorm(force, nforce);
      outarray[6 + m*ci] = podptr->podArrayNorm(DFTforce, nforce);

      double diff, sum = 0.0, ssr = 0.0;
      for (int j=0; j<dim*natom; j++) {
        diff = DFTforce[j] - force[j];
        sum += fabs(diff);
        ssr += diff*diff;
      }
      outarray[7 + m*ci] = sum/nforce;

      emae += outarray[4 + m*ci];
      essr += outarray[4 + m*ci]*outarray[4 + m*ci];
      fmae += sum;
      fssr += ssr;
      nforceall += nforce;
      ci += 1;

    }
    int q = file + 1;
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
  errors[0] = errors[0]/nc;
  errors[1] = sqrt(errors[1]/nc);
  errors[2] = errors[2]/nf;
  errors[3] = sqrt(errors[3]/nf);

  utils::logmesg(lmp, "**************** End of Error Calculation ****************\n");

  print_analysis(data, outarray, errors);
}

void CFITPOD::energyforce_calculation(datastruct data, double *coeff)
{
  int dim = 3;
  double energy;
  double force[1+dim*data.num_atom_max];

  int nfiles = data.data_files.size();  // number of files


  utils::logmesg(lmp, "**************** Begin of Energy/Force Calculation ****************\n");

  int ci = 0; // configuration counter
  for (int file = 0; file < nfiles; file++) { // loop over each file in the data set

    int nconfigs = data.num_config[file];
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      if ((ci % 100)==0) utils::logmesg(lmp, "Configuration: # {}\n", ci+1);

      int natom = data.num_atom[ci];
      int nforce = dim*natom;

      energy = energyforce_calculation(&force[1], coeff, data, ci);

      ci += 1;

      // save energy and force into a binary file

      force[0] = energy;
      std::string filename = "energyforce_config" + std::to_string(ci) + ".bin";

      FILE *fp = fopen(filename.c_str(), "wb");

      fwrite( reinterpret_cast<char*>( &force[0] ), sizeof(double) * (1 + nforce), 1, fp);

      fclose(fp);
    }
  }
  utils::logmesg(lmp, "**************** End of Energy/Force Calculation ****************\n");
}
