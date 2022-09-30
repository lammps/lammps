/***************************************************************************                 
               CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __COMPUTEPODFIT
#define __COMPUTEPODFIT

#include "compute_podfit.h"

#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <glob.h>
#include <random>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>

using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using namespace LAMMPS_NS;

enum{SCALAR,VECTOR,ARRAY};

CPODFIT::CPODFIT(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),  podptr(nullptr)
{   
  if (narg < 5) error->all(FLERR,"Illegal compute podfit command");
    
  std::string pod_file = std::string(arg[3]);  // pod input file
  std::string data_file = std::string(arg[4]); // data input file       
  std::string coeff_file;             // coefficient input file       
    
  if (narg > 5)
    coeff_file = std::string(arg[5]); // coefficient input file       
  else
    coeff_file = "";
  
  podptr = new CPOD(lmp, pod_file, coeff_file);
  this->read_data_files(data_file, podptr->pod.species);     
  
  if ((int) traindata.data_path.size() > 1) 
    this->allocate_memory(traindata);  
  else if ((int) testdata.data_path.size() > 1)
    this->allocate_memory(testdata);  
        
  // get POD coefficients from an input file       
  if (coeff_file != "") podptr->podArrayCopy(desc.c, podptr->pod.coeff, podptr->pod.nd);    
}
  
// destructor    
CPODFIT::~CPODFIT()
{
  traindata.freememory(1);
  testdata.freememory(1);
  desc.freememory(1);
  nb.freememory(1);
}

void CPODFIT::init()
{
  //cout<<"POD init"<<endl;
  //cout<<modify->get_compute_by_style("podfit").size()<<endl;
  
  if (force->pair == nullptr)
  error->all(FLERR,"Compute podfit requires a pair style be defined");

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("podfit").size() > 1 && comm->me == 0)
  error->warning(FLERR,"More than one compute podfit");
  
//   memory->create(coeff, podptr->pod.nd, podptr->pod.nelements+1, "podfit:coeff");
//   for (int i=0; i<podptr->pod.nd; i++)
//     coeff[i][0] = 0.0;  
//   array = coeff;  
  
//   // find compute for reference energy
//   std::string id_pe = std::string("thermo_pe");
//   int ipe = modify->find_compute(id_pe);
//   if (ipe == -1)
//   error->all(FLERR,"compute thermo_pe does not exist.");
//   c_pe = modify->compute[ipe];
  
  // compute POD coefficients using least-squares method
  this->least_squares_fit(traindata);
  
  // calculate errors for the training data set
  if ((traindata.training_analysis) && ((int) traindata.data_path.size() > 1) )
    this->error_analsysis(traindata, desc.c);  
  
  // calculate errors for the test data set
  if ((testdata.test_analysis) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path)) 
    this->error_analsysis(testdata, desc.c);  
  
  // calculate energy and force for the training data set
  if ((traindata.training_calculation) && ((int) traindata.data_path.size() > 1) )
    this->energyforce_calculation(traindata, desc.c);   
  
  // calculate energy and force for the test data set
  if ((testdata.test_calculation) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path) )
    this->energyforce_calculation(testdata, desc.c);     
}

/* ---------------------------------------------------------------------- */

void CPODFIT::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

void CPODFIT::compute_array()
{
//   // compute POD coefficients using least-squares method
//   this->least_squares_fit(traindata);
//   
//   // calculate errors for the training data set
//   if ((traindata.training_analysis) && ((int) traindata.data_path.size() > 1) )
//     this->error_analsysis(traindata, desc.c);  
//   
//   // calculate errors for the test data set
//   if ((testdata.test_analysis) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path)) 
//     this->error_analsysis(testdata, desc.c);  
//   
//   // calculate energy and force for the training data set
//   if ((traindata.training_calculation) && ((int) traindata.data_path.size() > 1) )
//     this->energyforce_calculation(traindata, desc.c);   
//   
//   // calculate energy and force for the test data set
//   if ((testdata.test_calculation) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path) )
//     this->energyforce_calculation(testdata, desc.c);     
}

template <typename T> void CPODFIT::writearray2file(const char* filename, T *a, int N, int backend)
{
  if (N>0) {    
    // Open file to read
    ofstream out(filename, ios::out | ios::binary);

    if (!out) {
      error->all(FLERR,"Unable to open file ");            
    }

    if (backend==2) { //GPU
#ifdef  USE_CUDA            
      T *a_host;      
      a_host = (T*) malloc (sizeof (T)*N);      
      
      // transfer data from GPU to CPU to save in a file
      cudaMemcpy(&a_host[0], &a[0], N*sizeof(T), cudaMemcpyDeviceToHost);  
      
      out.write( reinterpret_cast<char*>( &a_host[0] ), sizeof(T) * N );
      
      free(a_host);
#endif      
    }
    else 
      out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );          
    
    out.close();
  }
}

bool CPODFIT::is_a_number(std::string line)
{  
  return isdigit(line.at(0));
}

void CPODFIT::read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension, 
    std::string &test_path, std::string &training_path, std::string data_file)
{
  std::ifstream file_in(data_file);
  if (!file_in) error->all(FLERR,"Error: Data input file is not found");       
  
  std::string line;    
  while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
  {                      
    if (line != "") {
      std::string s, s2;
      double d;
      
      std::istringstream ss_line(line);                  
      ss_line >> s;
                  
      if (s == "fitting_weight_energy") {
        ss_line >> d;       
        fitting_weights[0] = d;         
      }
      else if (s == "fitting_weight_force") {
        ss_line >> d;       
        fitting_weights[1] = d;         
      }
      else if (s == "fitting_weight_stress") {
        ss_line >> d;       
        fitting_weights[2] = d;         
      }
      else if (s == "error_analysis_for_training_data_set") {
        ss_line >> d;       
        fitting_weights[3] = d;         
      }
      else if (s == "error_analysis_for_test_data_set") {
        ss_line >> d;       
        fitting_weights[4] = d;         
      }
      else if (s == "energy_force_calculation_for_training_data_set") {
        ss_line >> d;       
        fitting_weights[5] = d;         
      }
      else if (s == "energy_force_calculation_for_test_data_set") {
        ss_line >> d;       
        fitting_weights[6] = d;         
      }
      else if (s == "percentage_training_data_set") {
        ss_line >> d;       
        fitting_weights[7] = d;         
      }
      else if (s == "percentage_test_data_set") {
        ss_line >> d;       
        fitting_weights[8] = d;         
      }
      else if (s == "randomize_training_data_set") {
        ss_line >> d;       
        fitting_weights[9] = d;         
      }
      else if (s == "randomize_test_data_set") {
        ss_line >> d;       
        fitting_weights[10] = d;         
      }
      else if (s != "#") {
        ss_line >> s2;        
        if (s == "file_format") file_format = s2;
        if (s == "file_extension") file_extension = s2;
        if (s == "path_to_training_data_set") {
          training_path = s2;          
          while(ss_line >> s2){  
            training_path = training_path + " " + s2;   
          }            
          //training_path.erase(std::remove(training_path.begin(), training_path.end(), '"'), training_path.end());                    
          training_path.erase(training_path.begin());          
          training_path.erase(training_path.end()-1);                    
        }        
        if (s == "path_to_test_data_set") {
          test_path = s2;  
          while (ss_line >> s2) {       
            test_path = test_path + " " + s2;
          }
          //test_path.erase(remove(test_path.begin(), test_path.end(), '"'), test_path.end());
          test_path.erase(test_path.begin());          
          test_path.erase(test_path.end()-1);                    
        }
      }
    }
  }    
  file_in.close();
  
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

std::vector<std::string> CPODFIT::globVector(const std::string& pattern, std::vector<std::string> & files)
{
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    std::string s = string(glob_result.gl_pathv[i]);
    //std::cout << s << "\n";
    files.push_back(s);
  }
  globfree(&glob_result);
  return files;
}

void CPODFIT::get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension)  
{
  std::vector<std::string> res = this->globVector(datapath + "/*." + extension, files);
} 

int CPODFIT::get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)  
{
  std::ifstream datafile(file);
  if (!datafile) {/*error*/}  
         
  std::string line; 
  int num_configs = 0;
  num_atom_sum = 0;
  while (std::getline(datafile, line)) // Read next line to `line`, stop if no more lines.
  {
    int d;
    if (this->is_a_number(line)) {      
      d = std::stoi(line);
      num_atom.push_back(d);                    
      num_configs += 1;
      num_atom_sum += d;
    }    
  }  
  datafile.close();
  
  return num_configs;
}

int CPODFIT::get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)  
{
  int nfiles = training_files.size(); // number of files
  int d, n;
  
  for (int i=0; i<nfiles; i++) {
    d = this->get_number_atom_exyz(num_atom, n, training_files[i]);  
    num_config.push_back(d);  
    num_atom_sum.push_back(n);  
  }
    
  int num_atom_all = 0;
  for (int i=0; i< (int) num_atom.size(); i++)
    num_atom_all += num_atom[i];
  
  return num_atom_all;
}

void CPODFIT::read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces, 
    int *atomtype, std::string file, std::vector<std::string> species)  
{
  std::ifstream datafile(file);
  if (!datafile) {/*error*/}  
         
  std::string substr1 = "nergy";
  //std::string substr2 = "ttice";
  std::string substr3 = "tress";
    
  int cfi = 0, nat=0, k = 0, ns = species.size();
  double d; 
  
  std::string line;       
  while (std::getline(datafile, line)) // Read next line to `line`, stop if no more lines.
  {
    if (line.substr(1,6) == "attice") {
      int index1 = line.find(substr1);
      //int index2 = line.find(substr2);
      int index3 = line.find(substr3);
                  
      if (line.substr(index1,6) == "nergy=") {    
        std::string s1 = line.substr(index1+6,-1);
        int ind = s1.find(" ");
        energy[cfi] = stod(line.substr(index1+6,ind));
      }
      
      //std::cout<<line<<std::endl;
      if (line.substr(1,6) == "attice") {    
        //std::string s1 = line.substr(0,-1);
        int ind1 = line.find("\"");
        int ind2 = line.find("\"",ind1+1);
        std::istringstream ss(line.substr(ind1+1,ind2-ind1));  
        //std::cout<<line.substr(ind1+1,ind2-ind1-1)<<std::endl;
        k = 0;
        while(ss >> d){      
          lattice[k + 9*cfi] = d;
          k += 1;
        }        
      }
      
      if (line.substr(index3,6) == "tress=") {    
        std::string s1 = line.substr(index3+7,-1);
        int ind = s1.find("\"");
        std::istringstream ss(line.substr(index3+7,ind));  
        k = 0;
        while(ss >> d){      
          stress[k + 9*cfi] = d;
          k += 1;
        }        
      }
            
      cfi += 1;      
    }
    else if (!this->is_a_number(line)) {      
      std::string s0;
      std::istringstream ss(line);            
      ss >> s0;
          
      for (int ii=0; ii<ns; ii++)
        if (species[ii] == s0)
          atomtype[nat] = ii+1;
        
      k = 0;
      while(ss >> d){      
        if (k <=2 ) pos[k + 3*nat] = d;
        if (k > 2 ) forces[k-3 + 3*nat] = d;
        k += 1;
      }      
      nat += 1;       
    }        
  }  
  datafile.close();  
}

void CPODFIT::get_data(datastruct &data, std::vector<std::string> species)
{  
  this->get_exyz_files(data.data_files, data.data_path, data.file_extension);
  data.num_atom_sum = this->get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);      
  data.num_config_sum = data.num_atom.size();
  
  utils::logmesg(lmp, "data file   |  number of configurations   |   number of atoms\n");
  for (int i=0; i< (int) data.data_files.size(); i++) {
    string filename = data.data_files[i].substr(data.data_path.size()+1,data.data_files[i].size());
    data.filenames.push_back(filename.c_str());
    utils::logmesg(lmp, "{}   |   {}   |   {}\n", data.filenames[i], data.num_config[i], data.num_atom_each_file[i]);
  }  
  utils::logmesg(lmp, "number of files: {}\n", data.data_files.size());
  utils::logmesg(lmp, "number of configurations in all files: {}\n", data.num_config_sum);
  utils::logmesg(lmp, "number of atoms in all files: {}\n", data.num_atom_sum);    
  
  int n = data.num_config_sum;
  data.lattice = (double *) malloc(9*n*sizeof(double));
  data.stress = (double *) malloc(9*n*sizeof(double));
  data.energy = (double *) malloc(n*sizeof(double));  
  n = data.num_atom_sum;
  data.position = (double *) malloc(3*n*sizeof(double));
  data.force = (double *) malloc(3*n*sizeof(double));
  data.atomtype = (int *) malloc(n*sizeof(int));  
  
  int nfiles = data.data_files.size(); // number of files  
  int nconfigs = 0;
  int natoms = 0;
  for (int i=0; i<nfiles; i++) {    
    this->read_exyz_file(&data.lattice[9*nconfigs], &data.stress[9*nconfigs], &data.energy[nconfigs], 
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
  int dim = 3;
  double Qmat[dim*dim];
  for (int ci=0; ci<len; ci++) {    
    int natom = data.num_atom[ci];
    int natom_cumsum = data.num_atom_cumsum[ci];  
    double *x = &data.position[dim*natom_cumsum];
    double *f = &data.force[dim*natom_cumsum]; 
    double *lattice = &data.lattice[9*ci];
    double *a1 = &lattice[0];
    double *a2 = &lattice[3];
    double *a3 = &lattice[6];
    
    podptr->matrix33_inverse(Qmat, a1, a2, a3);   
    podptr->triclinic_lattice_conversion(a1, a2, a3, a1, a2, a3);                 
    podptr->matrix33_multiplication(Qmat, lattice, Qmat, dim);         
    podptr->matrix33_multiplication(x, Qmat, x, natom);   
    podptr->matrix33_multiplication(f, Qmat, f, natom);   
  }
    
  utils::logmesg(lmp, "minimum number of atoms: {}\n", data.num_atom_min);
  utils::logmesg(lmp, "maximum number of atoms: {}\n", data.num_atom_max); 
}

std::vector<int> CPODFIT::linspace(int start_in, int end_in, int num_in)
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

std::vector<int> CPODFIT::shuffle(int start_in, int end_in, int num_in)
{
  int sz = end_in - start_in + 1;  
  std::vector<int> myvector(sz);

  for (int i = 0; i<sz; i++)
    myvector[i] = start_in + i;

//   auto rng = std::default_random_engine {};
//   std::shuffle (myvector.begin(), myvector.end(), rng);  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle (myvector.begin(), myvector.end(), std::default_random_engine(seed));
  
  std::vector<int> shuffle_vec(num_in);
  for (int i = 0; i<num_in; i++)
    shuffle_vec[i] = myvector[i];  
  
  return shuffle_vec;
}  
  
std::vector<int> CPODFIT::select(int n, double percentage, int randomize)
{
  std::vector<int> selected;
    
  int m = (int) std::round(n*percentage);  
  m = PODMAX(m, 1);
  
  selected = (randomize==1) ? this->shuffle(1, n, m) : this->linspace(1, n, m);
    
  return selected;
}

void CPODFIT::select_data(datastruct &newdata, datastruct data)
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
    //cout<<data.filenames[file]<<endl;
    int nconfigs = data.num_config[file];
    selected[file] = select(nconfigs, percentage, randomize);    
    int ns = (int) selected[file].size(); // number of selected configurations   
//     cout<<nconfigs<<"  "<<ns<<endl;
//     for (int ii=0; ii < ns; ii++) 
//       cout<<selected[file][ii]<<"  ";
//     cout<<endl;
    
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
  newdata.lattice = (double *) malloc(9*n*sizeof(double));
  newdata.stress = (double *) malloc(9*n*sizeof(double));
  newdata.energy = (double *) malloc(n*sizeof(double));  
  n = data.num_atom_sum;
  newdata.position = (double *) malloc(3*n*sizeof(double));
  newdata.force = (double *) malloc(3*n*sizeof(double));
  newdata.atomtype = (int *) malloc(n*sizeof(int));  
    
  int cn = 0, dim=3;
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
    string filename = newdata.data_files[i].substr(newdata.data_path.size()+1,newdata.data_files[i].size());
    newdata.filenames.push_back(filename.c_str());
    utils::logmesg(lmp, "{}   |   {}   |   {}   |   {}   |   {}\n", newdata.filenames[i], newdata.num_config[i], newdata.num_atom_each_file[i], data.num_config[i], data.num_atom_each_file[i]);              
  }
  utils::logmesg(lmp, "number of files: {}\n", newdata.data_files.size());
  utils::logmesg(lmp, "number of configurations in all files (selected and original): {} and {}\n", newdata.num_config_sum, data.num_config_sum);
  utils::logmesg(lmp, "number of atoms in all files (selected and original: {} and {}\n", newdata.num_atom_sum, data.num_atom_sum);    
}

void CPODFIT::read_data_files(std::string data_file, std::vector<std::string> species)
{  
  datastruct data;
      
  // read data input file to datastruct
  this->read_data_file(data.fitting_weights, data.file_format, data.file_extension, 
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
      this->get_data(traindata, species);
    else 
      error->all(FLERR,"data set is not found");     
    utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");               
  }
  else {
    utils::logmesg(lmp, "**************** Begin of Training Data Set ****************\n");    
    if ((int) data.data_path.size() > 1) 
      this->get_data(data, species);
    else 
      error->all(FLERR,"data set is not found");            
    utils::logmesg(lmp, "**************** End of Training Data Set ****************\n");          
    
    utils::logmesg(lmp, "**************** Begin of Select Training Data Set ****************\n");
    select_data(traindata, data);
    utils::logmesg(lmp, "**************** End of Select Training Data Set ****************\n");
    
    data.freememory(1);
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
    this->get_data(testdata, species);
    utils::logmesg(lmp, "**************** End of Test Data Set ****************\n"); 
  }
  else {
    testdata.data_path = traindata.data_path;
  }      
}

int CPODFIT::latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  int m=0, n=0, p=0;
  if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);  
  if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);  
  if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);        
    
  // index for the center lattice
  int ind = m + (2*m+1)*(n) + (2*m+1)*(2*n+1)*(p);
  
  // number of lattices
  int nl = (2*m+1)*(2*n+1)*(2*p+1);      
    
  //y = zeros(3, nx*nl)
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
  
  //alist = zeros(Int32,nx*nl)
  for (int i=0; i <nl; i++)
    for (int j=0; j<nx; j++) 
      alist[j + nx*i] = j;      
  
  return nl;
}

int CPODFIT::podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
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

int CPODFIT::podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum, 
    double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  double rcutsq = rcut*rcut;  
  int dim = 3, nl = 0, nn = 0;
  
  // number of lattices
  nl = this->latticecoords(y, alist, x, a1, a2, a3, rcut, pbc, nx);    
  int N = nx*nl;
      
  // total number of neighbors
   nn = this->podneighborlist(neighlist, numneigh, y, rcutsq, nx, N, dim);
  
   podptr->podCumsum(numneighsum, numneigh, nx+1); 
     
   return nn;
}   

void CPODFIT::allocate_memory(datastruct data)
{
  int nd = podptr->pod.nd;
  desc.gd = (double *) malloc(nd*sizeof(double));
  desc.A = (double *) malloc(nd*nd*sizeof(double));  
  desc.b = (double *) malloc(nd*sizeof(double));  
  desc.c = (double *) malloc(nd*sizeof(double));  
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
    ny = PODMAX(ny,dim*natom*nl);
    na = PODMAX(na, natom*nl);      
    np = PODMAX(np, natom*natom*nl);    
  }
    
  nb.y = (double*) malloc (sizeof (double)*(ny));
  nb.alist = (int*) malloc (sizeof (int)*(na));  
  nb.pairnum = (int*) malloc (sizeof (int)*(natom_max));
  nb.pairnum_cumsum = (int*) malloc (sizeof (int)*(natom_max+1));
  nb.pairlist = (int*) malloc (sizeof (int)*(np));       
  
//   nb.elemindex = (int*) malloc (sizeof (int)*(nelements*nelements));       
//   int k = 1;
//   for (int i=0; i < nelements; i++) 
//     for (int j=i; j < nelements; j++) {
//       nb.elemindex[i + nelements*j] = k;
//       nb.elemindex[j + nelements*i] = k;
//       k += 1;
//     }      
  
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
      
    // neighbor list
    Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, x, a1, a2, a3, rcut, pbc, natom);
  
    int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
    int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

//     Nj=0, Nijk=0;
//     for (int i=0; i < natom; i++) {
//       Nj = (Nj > nb.pairnum[i]) ? Nj : nb.pairnum[i];
//       Nijk +=  (nb.pairnum[i]-1)*nb.pairnum[i]/2;
//     }

    //int szd1 = 3*Nij+ (1+dim)*Nij*(nrbf3+ns3) + 2*(1+dim)*Nijk*nrbf3 + 4*Nj*nrbf3 + 2*dim*Nijk + 7*Nijk*nabf3;
    //int szi1 = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;        
    int szd1 = 3*Nij+ (1+dim)*Nij*PODMAX(nrbf2+ns2,nrbf3+ns3) + (nabf3+1)*7;
    int szi1 = 6*Nij + 2*natom+1 + (Nj-1)*Nj;        
    szd = PODMAX(szd, szd1);   
    szi = PODMAX(szi, szi1);   
    
    if (podptr->sna.twojmax>0) {
      szd1 = 0;
      szd1 += Nij*dim; // rij
      szd1 += PODMAX(2*podptr->sna.idxu_max*Nij, 2*podptr->sna.idxz_max*podptr->sna.ndoubles*natom); // (Ur, Ui) and (Zr, Zi)       
      szd1 += 2*podptr->sna.idxu_max*dim*Nij; // dUr, dUi      
      szd1 += PODMAX(podptr->sna.idxb_max*podptr->sna.ntriples*dim*Nij, 2*podptr->sna.idxu_max*podptr->sna.nelements*natom); // dblist and (Utotr, Utoti)                  
      szsnap = PODMAX(szsnap, szd1);   
    }    
    //std::cout<<"Ni, Nj, Nij, Nijk: "<<natom<<", "<<Nj<<", "<<Nij<<", "<<Nijk<<std::endl;
  }    
  
  szd = PODMAX(szsnap, szd);   
  szd = PODMAX(natom_max*(nd1+nd2+nd3+nd4) + szd, dim*natom_max*(nd-nd1-nd2-nd3-nd4));
  szd = dim*natom_max*(nd1+nd2+nd3+nd4) + szd;         
  
  // gdd includes linear descriptors derivatives, quadratic descriptors derivatives and temporary memory
  desc.gdd = (double*) malloc (sizeof (double)*(szd)); // [ldd qdd]
  desc.tmpint = (int*) malloc (sizeof (int)*(szi));  
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

void CPODFIT::linear_descriptors(datastruct data, int ci)
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
  int Nij = this->podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, 
        position, a1, a2, a3, rcut, pbc, natom);

  int *tmpint = &desc.tmpint[0];   
  double *tmpmem = &desc.gdd[dim*natom*nd1234+natom*nd1234];
  podptr->linear_descriptors(desc.gd, desc.gdd, nb.y, tmpmem, atomtype, nb.alist, 
      nb.pairlist, nb.pairnum, nb.pairnum_cumsum, tmpint, natom, Nij);   
  
}

void CPODFIT::quadratic_descriptors(datastruct data, int ci)
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

void CPODFIT::cubic_descriptors(datastruct data, int ci)
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

void CPODFIT::least_squares_matrix(datastruct data, int ci)
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

void CPODFIT::least_squares_fit(datastruct data)
{    
  utils::logmesg(lmp, "**************** Begin of Least-Squares Fitting ****************\n");
  
  // loop over each configuration in the training data set

  for (int ci=0; ci < (int) data.num_atom.size(); ci++) {
    if ((ci % 100)==0) utils::logmesg(lmp, "Configuration: # {}\n", ci+1);
    
    // compute linear POD descriptors

    this->linear_descriptors(data, ci);
    
    // compute quadratic POD descriptors

    this->quadratic_descriptors(data, ci);    
    
    // compute cubic POD descriptors

    this->cubic_descriptors(data, ci);  
    
    // assemble the least-squares linear system

    this->least_squares_matrix(data, ci);      
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

  std::ofstream myfile ("coefficients.txt");
  if (myfile.is_open())
  {
    myfile << "POD_coefficients: " + std::to_string(nd) + "\n";
    for(int count = 0; count < nd; count ++)
      myfile << std::setprecision(20)<< desc.c[count] << std::endl;    
    myfile.close();
  }
  else utils::logmesg(lmp, "Unable to open file\n");
  
  utils::logmesg(lmp, "**************** End of Least-Squares Fitting ****************\n");
}

double CPODFIT::energyforce_calculation(double *force, double *coeff, datastruct data, int ci)
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

  int Nij = this->podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, 
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

void CPODFIT::print_analysis(datastruct data, double *outarray, double *errors)
{        
  string s = "All files";
  int nfiles = data.data_files.size();  // number of files  
  int lm = s.size();  
  for (int i = 0; i < nfiles; i++) 
    lm = PODMAX(lm, (int) data.filenames[i].size());        
  lm = lm + 2;
  
  std::string filename = data.training ? "training_errors.txt" : "test_errors.txt";   
  std::ofstream myfile (filename);
  if (!myfile.is_open()) utils::logmesg(lmp, "Unable to open file");
  
  filename = data.training ? "training_analysis.txt" : "test_analysis.txt";
  std::ofstream mfile (filename);
  if (!mfile.is_open()) utils::logmesg(lmp, "Unable to open file");
    
  std::string sa = "**************** Begin of Error Analysis for the Training Data Set ****************";
  std::string sb = "**************** Begin of Error Analysis for the Test Data Set ****************";      
  std::string mystr = (data.training) ? sa : sb;   
  utils::logmesg(lmp, "{}\n", mystr);    
  myfile <<mystr+"\n";
  
  sa = "----------------------------------------------------------------------------------------\n";
  sb = "  File    | # configs | # atoms  | MAE energy | RMSE energy | MAE force | RMSE force |\n";
  utils::logmesg(lmp, "{}", sa); 
  utils::logmesg(lmp, "{}", sb);
  utils::logmesg(lmp, "{}", sa);
  myfile << sa;
  myfile << sb;
  myfile << sa;   
  
  int ci=0, m=8, nc=0, nf=0;  
  for (int file = 0; file < nfiles; file++) {    
    mfile<<"# " + data.filenames[file] + "\n";
    sb = "|  config  |  # atoms  |  energy  | DFT energy | energy error |  force  | DFT force | force error |\n";
    mfile <<sb;
    
    int nforceall = 0;
    int nconfigs = data.num_config[file];
    nc += nconfigs;
    for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
      
      mfile <<"   "; 
      for(int count = 0; count < m; count ++)
        mfile << outarray[count + m*ci] << "    ";  
      mfile<<std::endl;
      
      nforceall += 3*data.num_atom[ci];      
      ci += 1;
    }
    nf += nforceall;
    
    int q = file+1;
    string s = data.filenames[file];
    s = s + std::string(lm-s.size(), ' ');     
    string s1 = std::to_string(nconfigs);
    s1 = s1 + std::string(PODMAX(6- (int) s1.size(),1), ' ');  
    s = s + "   " + s1;
    s1 = std::to_string(nforceall/3);
    s1 = s1 + std::string(PODMAX(7 - (int) s1.size(),1), ' ');  
    s = s + "   " + s1;
    s1 = std::to_string(errors[0 + 4*q]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
    s = s + "   " + s1;
    s1 = std::to_string(errors[1 + 4*q]);
    s1 = s1 + std::string(PODMAX(10 - (int)  s1.size(),1), ' ');  
    s = s + "   " + s1;
    s1 = std::to_string(errors[2 + 4*q]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
    s = s + "   " + s1;
    s1 = std::to_string(errors[3 + 4*q]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
    s = s + "   " + s1 + "\n";    
    utils::logmesg(lmp, "{}", s);
    myfile <<s;    
  }
  utils::logmesg(lmp, "{}", sa);
  myfile << sa;
    
  s = s + std::string(PODMAX(lm - (int) s.size(),1), ' ');     
  string s1 = std::to_string(nc);
  s1 = s1 + std::string(PODMAX(6- (int) s1.size(),1), ' ');  
  s = s + "   " + s1;
  s1 = std::to_string(nf/3);
  s1 = s1 + std::string(PODMAX(7 - (int) s1.size(),1), ' ');  
  s = s + "   " + s1;
  s1 = std::to_string(errors[0]);
  s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
  s = s + "   " + s1;
  s1 = std::to_string(errors[1]);
  s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
  s = s + "   " + s1;
  s1 = std::to_string(errors[2]);
  s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
  s = s + "   " + s1;
  s1 = std::to_string(errors[3]);
  s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');  
  s = s + "   " + s1 + "\n";  
  utils::logmesg(lmp, "{}", s);
  utils::logmesg(lmp, "{}", sa);
  myfile << s;
  myfile << sa;
  
  sa = "**************** End of Error Analysis for the Training Data Set ****************";    
  sb = "**************** End of Error Analysis for the Test Data Set ****************";   
  mystr = (data.training) ? sa : sb;       
  utils::logmesg(lmp, "{}\n", mystr); 
  myfile <<mystr+"\n";  
  myfile.close();
  mfile.close();  
}

void CPODFIT::error_analsysis(datastruct data, double *coeff)
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
                  
      energy = this->energyforce_calculation(force, newcoeff, data, ci);
      
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

void CPODFIT::energyforce_calculation(datastruct data, double *coeff)
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
      
      energy = this->energyforce_calculation(&force[1], coeff, data, ci);
      
      ci += 1;       
      
      // save energy and force into a binary file

      force[0] = energy;
      string filename = "energyforce_config" + std::to_string(ci) + ".bin";
      writearray2file(filename.c_str(), force, 1 + nforce, 1);
    }
  }     
  utils::logmesg(lmp, "**************** End of Energy/Force Calculation ****************\n");
}

#endif    
