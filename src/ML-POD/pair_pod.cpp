
#include "pair_pod.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "tokenizer.h"

#include <glob.h>

#include <cmath>
#include <cstring>
#include <sys/time.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

CPairPOD::CPairPOD(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

CPairPOD::~CPairPOD()
{
  data.freememory(1);
  this->free_memory();
  TemplateFree(gd, backend);
  TemplateFree(energycoeff, backend);
  TemplateFree(forcecoeff, backend);
  TemplateFree(podcoeff, backend);
  TemplateFree(newpodcoeff, backend);

  delete podptr;

  if (allocated) {
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(scale);
  }
}

void CPairPOD::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);
  vflag_fdotr = 1;

  double **x = atom->x;
  double **f = atom->f;
  int **firstneigh = list->firstneigh;
  int *numneigh = list->numneigh;
  int *type = atom->type;
  int *ilist = list->ilist;
  int nlocal = atom->nlocal;
  int inum = list->inum;
  int nall = inum + atom->nghost;
  
  // initialize global descriptors to zero
  
  int nd1234 = podptr->pod.nd1234;
  podptr->podArraySetValue(gd, 0.0, nd1234);

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int jnum = numneigh[i];

    // allocate temporary memory
    
    if (nijmax < jnum) {
      nijmax = PODMAX(nijmax, jnum);
      nablockmax = 1;
      this->free_tempmemory();
      this->estimate_tempmemory();
      this->allocate_tempmemory();
    }

    // get neighbor pairs for atom i
    
    this->lammpsNeighPairs(x, firstneigh, type, numneigh, i);

    // compute global POD descriptors for atom i
    
    podptr->linear_descriptors_ij(gd, tmpmem, rij, &tmpmem[nd1234], numneighsum,
      typeai, idxi, ti, tj, 1, nij);
  }

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
  int nd = podptr->pod.nd;  
  bigint natom = atom->natoms;
  
  for (int j=nd1234; j<(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j++)
    newpodcoeff[j] = podcoeff[j]/(natom);

  for (int j=(nd1234+nd22+nd23+nd24+nd33+nd34+nd44); j<nd; j++)
    newpodcoeff[j] = podcoeff[j]/(natom*natom);
  
  // compute energy and effective coefficients
  
  eng_vdwl = podptr->calculate_energy(energycoeff, forcecoeff, gd, newpodcoeff);

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];

    // get neighbor pairs for atom i
    
    this->lammpsNeighPairs(x, firstneigh, type, numneigh, i);

    // compute atomic force for atom i
    
    podptr->calculate_force(f, forcecoeff, rij, tmpmem, numneighsum,
      typeai, idxi, ai, aj, ti, tj, 1, nij);
  }

  if (vflag_fdotr) virial_fdotr_compute();


//  printf("%d %d %d %d\n", eflag, vflag, vflag_fdotr, eflag_atom);
//
//   if (eflag_atom) {
//     eng_vdwl = this->lammpseatom(eatom, atom->x, list->firstneigh, atom->type, list->numneigh,
//             list->ilist, list->inum, list->inum + atom->nghost);
//
//     this->lammpsforce(atom->f, atom->x, list->firstneigh, atom->type, list->numneigh,
//             list->ilist, list->inum, list->inum + atom->nghost);
//   }
//   else {
//     eng_vdwl = this->lammpsenergyforce(atom->f, atom->x, list->firstneigh, atom->type, list->numneigh,
//           list->ilist, list->inum, list->inum + atom->nghost);
//   }
//
// //   podptr->print_matrix("atom->f", 3, list->inum + atom->nghost, atom->f, 3);
// //   podptr->print_matrix("atom->x", 3, list->inum + atom->nghost, atom->x, 3);
// //   cout<<"Energy: "<< eng_vdwl<<endl;
// //   podptr->print_matrix("virial", 1, 6, virial, 1);
//
//   if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void CPairPOD::settings(int narg, char ** /* arg */)
{
  if (narg > 0)
  error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void CPairPOD::coeff(int narg, char **arg)
{
  // set default scaling  
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(scale,n+1,n+1,"pair:scale");
  map = new int[n+1];
  for (int ii = 0; ii < n+1; ii++)
    for (int jj = 0; jj < n+1; jj++)
      scale[ii][jj] = 1.0;
  allocated = 1;

  if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");
//  if (narg != 5 + atom->ntypes) error->all(FLERR,"Incorrect args for pair coefficients");

  map_element2type(narg-4,&arg[4]);
//   cout<<map[0]<<endl;
//   cout<<map[1]<<endl;

  std::string pod_file = std::string(arg[2]);  // pod input file
  std::string coeff_file = std::string(arg[3]); // coefficient input file
  std::string data_file;             // data input file

//   if (narg>4) data_file = std::string(arg[4]); // data input file
//   else data_file = "";

  data_file = "";
  this->InitPairPOD(pod_file, coeff_file, data_file);

  for (int ii = 0; ii < n+1; ii++)
    for (int jj = 0; jj < n+1; jj++)
      cutsq[ii][jj] = podptr->pod.rcut*podptr->pod.rcut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void CPairPOD::init_style()
{
  if (force->newton_pair == 0)
  error->all(FLERR,"Pair style POD requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double CPairPOD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  scale[j][i] = scale[i][j];
  return podptr->pod.rcut;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double CPairPOD::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes;
}

void *CPairPOD::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}

void CPairPOD::InitPairPOD(std::string pod_file, std::string coeff_file, std::string data_file)
{
  podptr = new CPOD(lmp, pod_file, coeff_file);

  if (data_file != "")
    this->read_data_files(data_file, podptr->pod.species);
  //else
  //  error->all(FLERR,"\n Error: Data filename is empty");

  if ((int) data.data_path.size() > 1) {
    podpairlist = 1;
    lammpspairlist = 0;
  }
  else {
    podpairlist = 0;
    lammpspairlist = 1;
  }

  // get POD coefficients from an input file

  if (coeff_file != "") {
    TemplateMalloc(&podcoeff, podptr->pod.nd, backend);
    TemplateMalloc(&newpodcoeff, podptr->pod.nd, backend);
    TemplateMalloc(&energycoeff, podptr->pod.nd1234, backend);
    TemplateMalloc(&forcecoeff, podptr->pod.nd1234, backend);
    TemplateMalloc(&gd, podptr->pod.nd1234, backend);
    podptr->podArrayCopy(podcoeff, podptr->pod.coeff, podptr->pod.nd);
    podptr->podArrayCopy(newpodcoeff, podptr->pod.coeff, podptr->pod.nd);
  }
}

bool CPairPOD::is_a_number(std::string line)
{
  return isdigit(line.at(0));
}

void CPairPOD::read_data_file(double *inputs, std::string &file_format, std::string &file_extension,
  std::string &data_path, std::string data_file)
{
  std::ifstream file_in(data_file);
  if (!file_in) error->all(FLERR,"Error: Data input file is not found");

  std::string line;

  // Read next line to `line`, stop if no more lines

  while (std::getline(file_in, line))
  {
  if (line != "") {
    std::string s, s2;
    double d;

    std::istringstream ss_line(line);
    ss_line >> s;

    if (s == "error_analysis_for_data_set") {
    ss_line >> d;
    inputs[0] = d;
    }
    else if (s == "save_calculation_in_binary_files") {
    ss_line >> d;
    inputs[1] = d;
    }
    else if (s == "save_frequency") {
    ss_line >> d;
    inputs[2] = d;
    }
    else if (s == "run_molecular_dynamics_simulation") {
    ss_line >> d;
    inputs[3] = d;
    }
    else if (s == "number_of_atoms_per_computation_block") {
    ss_line >> d;
    inputs[4] = d;
    }
    else if (s == "random_rotation") {
    ss_line >> d;
    inputs[5] = d;
    }
    else if (s != "#") {
    ss_line >> s2;
    if (s == "file_format") file_format = s2;
    if (s == "file_extension") file_extension = s2;
    if (s == "path_to_data_set") {
      data_path = s2;
      while(ss_line >> s2){
      data_path = data_path + " " + s2;
      }
      //data_path.erase(std::remove(data_path.begin(), data_path.end(), '"'), data_path.end());
      data_path.erase(data_path.begin());
      data_path.erase(data_path.end()-1);
    }
    }
  }
  }
  file_in.close();

  if (inputs[4] == 0) inputs[4] = blocksize;

  std::cout<<"**************** Begin of Data File ****************"<<std::endl;
  std::cout<<"file format: "<<file_format<<std::endl;
  std::cout<<"file extension: "<<file_extension<<std::endl;
  std::cout<<"path to data set: "<<data_path<<std::endl;
  std::cout<<"error analysis for data set: "<<inputs[0]<<std::endl;
  std::cout<<"run MD simulation: "<<inputs[3]<<std::endl;
  std::cout<<"save calculation fin binary files: "<<inputs[1]<<std::endl;
  std::cout<<"save frequency: "<<inputs[2]<<std::endl;
  std::cout<<"number of atoms per computation block: "<<inputs[4]<<std::endl;
  std::cout<<"Apply random rotation to configurations: "<<inputs[5]<<std::endl;
  std::cout<<"**************** End of Data File ****************"<<std::endl<<std::endl;
}

std::vector<std::string> CPairPOD::globVector(const std::string& pattern, std::vector<std::string> & files)
{
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
  std::string s = string(glob_result.gl_pathv[i]);
  files.push_back(s);
  }
  globfree(&glob_result);
  return files;
}

void CPairPOD::get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension)
{
  std::vector<std::string> res = this->globVector(datapath + "/*." + extension, files);
}

int CPairPOD::get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)
{
  std::ifstream datafile(file);
  if (!datafile) {/*error*/}

  std::string line;
  int num_configs = 0;
  num_atom_sum = 0;

  // Read next line to `line`, stop if no more lines

  while (std::getline(datafile, line))
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

int CPairPOD::get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)
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

void CPairPOD::read_exyz_file(double *lattice, double *stress, double *energy, double *pos,
  double *vel, double *forces, int *atomtype, std::string file, std::vector<std::string> species)
{
  std::ifstream datafile(file);
  if (!datafile) {/*error*/}

  std::string substr1 = "nergy";
  std::string substr3 = "tress";

  int cfi = 0, nat=0, k = 0, ns = species.size();
  double d;

  std::string line;

  // Read next line to `line`, stop if no more lines

  while (std::getline(datafile, line))
  {
  if (line.substr(1,6) == "attice") {
    int index1 = line.find(substr1);
    int index3 = line.find(substr3);

    if (line.substr(index1,6) == "nergy=") {
    std::string s1 = line.substr(index1+6,-1);
    int ind = s1.find(" ");
    energy[cfi] = stod(line.substr(index1+6,ind));
    }

    if (line.substr(1,6) == "attice") {
    int ind1 = line.find("\"");
    int ind2 = line.find("\"",ind1+1);
    std::istringstream ss(line.substr(ind1+1,ind2-ind1));
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
    if (k >= 0 && k <=2) pos[k + 3*nat] = d;
    if (k >= 3 && k <=5) forces[k-3 + 3*nat] = d;
    if (k >= 6 && k <=8) vel[k-6 + 3*nat] = d;
    k += 1;
    }
    nat += 1;
  }
  }
  datafile.close();
}

void CPairPOD::get_data(std::vector<std::string> species)
{
  this->get_exyz_files(data.data_files, data.data_path, data.file_extension);
  data.num_atom_sum = this->get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);
  data.num_config_sum = data.num_atom.size();

  std::cout<<"data file   |  number of configurations   |   number of atoms "<<std::endl;
  for (int i=0; i< (int) data.data_files.size(); i++) {
  string filename = data.data_files[i].substr(data.data_path.size()+1,data.data_files[i].size());
  data.filenames.push_back(filename.c_str());
  std::cout<<data.filenames[i]<<"   |   "<<data.num_config[i]<<"   |   "<<data.num_atom_each_file[i]<< std::endl;
  //std::cout<<data.data_files[i].substr(data.data_path.size()+1,data.data_files[i].size())<<std::endl;
  }
  std::cout << "number of files: " <<data.data_files.size() << std::endl;
  std::cout << "number of configurations in all files: " <<data.num_config_sum << std::endl;
  std::cout << "number of atoms in all files: " <<data.num_atom_sum << std::endl;

  int n = data.num_config_sum;
  data.lattice = (double *) malloc(9*n*sizeof(double));
  data.stress = (double *) malloc(9*n*sizeof(double));
  data.energy = (double *) malloc(n*sizeof(double));
  n = data.num_atom_sum;
  data.position = (double *) malloc(3*n*sizeof(double));
  data.velocity = (double *) malloc(3*n*sizeof(double));
  data.force = (double *) malloc(3*n*sizeof(double));
  data.atomtype = (int *) malloc(n*sizeof(int));

  int nfiles = data.data_files.size(); // number of files
  int nconfigs = 0;
  int natoms = 0;
  for (int i=0; i<nfiles; i++) {
  this->read_exyz_file(&data.lattice[9*nconfigs], &data.stress[9*nconfigs], &data.energy[nconfigs],
    &data.position[3*natoms], &data.velocity[3*natoms], &data.force[3*natoms], &data.atomtype[natoms],
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

  std::cout << "minimum number of atoms: " <<data.num_atom_min << std::endl;
  std::cout << "maximum number of atoms: " <<data.num_atom_max << std::endl;
}

void CPairPOD::read_data_files(std::string data_file, std::vector<std::string> species)
{
  double inputs[100];
  for (int i=0; i<100; i++) inputs[i] = 0;

  // read data input file to datastruct

  this->read_data_file(inputs, data.file_format, data.file_extension, data.data_path, data_file);

  analysis = (int) inputs[0];
  savecalculation = (int) inputs[1];
  savefrequency = (int) inputs[2];
  runMD = (int) inputs[3];
  blocksize = (int) inputs[4];
  randomrotation = (int) inputs[5];

  if ((int) data.data_path.size() > 1)
  this->get_data(species);
  else
  error->all(FLERR,"data set is not found");
}

void CPairPOD::get_atomblocks(int natom)
{
  if (blocksize >= natom) {
  numblocks = 1;
  atomblocks[0] = 0;
  atomblocks[1] = natom;
  }
  else {
  numblocks = (int) ceil( ((double) natom)/((double) blocksize) );

  double delta = ((double) natom) / ((double) numblocks);
  for(int i=0; i < numblocks; ++i)
    atomblocks[i]= (int) delta * i;
  atomblocks[numblocks] = natom;
  }
  if (numblocks > 1000) error->all(FLERR,"number of computation blocks can not be more than 1000. This error can be fixed by increasing the number of atoms per computation block.");
}

int CPairPOD::latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
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

int CPairPOD::podneighborcount(double *r, double rcutsq, int nx, int N, int dim)
{
  int k = 0;
  for (int i = 0; i<nx; i++) {
  double *ri = &r[i*dim];
  for (int j=0; j<N; j++) {
    double *rj = &r[dim*j];
    double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
    if  ((rijsq > 1e-12) && (rijsq <= rcutsq)) k += 1;
  }
  }
  return k;
}

int CPairPOD::podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
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

int CPairPOD::podfullneighborlist(double *xy, int *alist, int *neighlist, int *numneigh, int *numneighsum,
  double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
  double rcutsq = rcut*rcut;
  int dim = 3, nl = 0, nn = 0;

  // number of lattices

  nl = this->latticecoords(xy, alist, x, a1, a2, a3, rcut, pbc, nx);
  int N = nx*nl;

  // total number of neighbors

   nn = this->podneighborlist(neighlist, numneigh, xy, rcutsq, nx, N, dim);

   podptr->podCumsum(numneighsum, numneigh, nx+1);

   return nn;
}

void CPairPOD::estimate_memory(datastruct data)
{

  int dim = 3;
  int natom_max = data.num_atom_max;
  int nd1 = podptr->pod.nd1;
  int nd2 = podptr->pod.nd2;
  int nd3 = podptr->pod.nd3;
  int nd4 = podptr->pod.nd4;
  int nbesselpars = podptr->pod.nbesselpars;
  int nrbf2 = podptr->pod.nbf2;
  int nabf3 = podptr->pod.nabf3;
  int nrbf3 = podptr->pod.nrbf3;
  int *pdegree2 = podptr->pod.twobody;
  int *pdegree3 = podptr->pod.threebody;
  int *pbc = podptr->pod.pbc;
  double rcut = podptr->pod.rcut;

  int Nij=0;
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

  double *y = (double*) malloc (sizeof (double)*(ny));
  int *atomID = (int*) malloc (sizeof (int)*(na));
  int *pairnum = (int*) malloc (sizeof (int)*(natom_max));
  int *pairnumsum = (int*) malloc (sizeof (int)*(natom_max+1));
  int *pairlist = (int*) malloc (sizeof (int)*(np));

  szd = 0, nijmax = 0;
  int szsnap=0;
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

  Nij = podfullneighborlist(y, atomID, pairlist, pairnum, pairnumsum, x, a1, a2, a3, rcut, pbc, natom);

  int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
  int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

  int szd1 = 3*Nij+ (1+dim)*Nij*PODMAX(nrbf2+ns2,nrbf3+ns3) + (nabf3+1)*7;
  szd = PODMAX(szd, szd1);
  nijmax = PODMAX(nijmax, Nij);

  if (podptr->sna.twojmax>0) {
    szd1 = 0;
    szd1 += Nij*dim; // rij
    szd1 += PODMAX(2*podptr->sna.idxu_max*Nij, 2*podptr->sna.idxz_max*podptr->sna.ndoubles*natom); // (Ur, Ui) and (Zr, Zi)
    szd1 += 2*podptr->sna.idxu_max*dim*Nij; // dUr, dUi
    szd1 += PODMAX(podptr->sna.idxb_max*podptr->sna.ntriples*dim*Nij, 2*podptr->sna.idxu_max*podptr->sna.nelements*natom); // dblist and (Utotr, Utoti)
    szsnap = PODMAX(szsnap, szd1);
  }
  }

  szd = PODMAX(szsnap, szd);
  szd = natom_max*(nd1+nd2+nd3+nd4) + szd;
  nlocalmax = natom_max;
  nmaxatom = ny;

  free(y); free(atomID); free(pairlist); free(pairnum); free(pairnumsum);
}

void CPairPOD::free_tempmemory()
{
  TemplateFree(rij, backend);
  TemplateFree(idxi, backend);
  TemplateFree(ai, backend);
  TemplateFree(aj, backend);
  TemplateFree(ti, backend);
  TemplateFree(tj, backend);
  TemplateFree(numneighsum, backend);
  TemplateFree(typeai, backend);
  TemplateFree(tmpmem, backend);
}

void CPairPOD::free_atommemory()
{
  TemplateFree(forces, backend);
  TemplateFree(stress, backend);
  if (atommemory) {
  TemplateFree(atomtype, backend);
  TemplateFree(pos, backend);
  TemplateFree(vel, backend);
  }
}

void CPairPOD::free_pairmemory()
{
  if (podpairlist) {
  TemplateFree(y, backend);
  TemplateFree(pairlist, backend);
  TemplateFree(pairnum, backend);
  TemplateFree(pairnumsum, backend);
  TemplateFree(atomID, backend);
  }
}

void CPairPOD::free_memory()
{
  this->free_tempmemory();
  this->free_atommemory();
  this->free_pairmemory();
}

void CPairPOD::allocate_tempmemory()
{
  TemplateMalloc(&rij, dim*nijmax, backend);
  TemplateMalloc(&idxi, nijmax, backend);
  TemplateMalloc(&ai, nijmax, backend);
  TemplateMalloc(&aj, nijmax, backend);
  TemplateMalloc(&ti, nijmax, backend);
  TemplateMalloc(&tj, nijmax, backend);
  TemplateMalloc(&numneighsum, nablockmax+1, backend);
  TemplateMalloc(&typeai, nablockmax, backend);
  TemplateMalloc(&tmpmem, szd, backend);
}

void CPairPOD::allocate_atommemory()
{
  TemplateMalloc(&forces, dim*nmaxatom, backend);
  TemplateMalloc(&stress, 9, backend);
  if (atommemory) {
  TemplateMalloc(&atomtype, nmaxatom, backend);
  TemplateMalloc(&pos, dim*nmaxatom, backend);
  TemplateMalloc(&vel, dim*nmaxatom, backend);
  }
}

void CPairPOD::allocate_pairmemory()
{
  if (podpairlist) {
  TemplateMalloc(&y, dim*nmaxatom, backend);
  TemplateMalloc(&atomID, nmaxatom, backend);
  TemplateMalloc(&pairnum, nlocalmax, backend);
  TemplateMalloc(&pairnumsum, nlocalmax+1, backend);
  TemplateMalloc(&pairlist, nmaxpairs, backend);
  }
}

void CPairPOD::allocate_memory()
{
  this->allocate_tempmemory();
  this->allocate_atommemory();
  this->allocate_pairmemory();
}

void CPairPOD::check_atommemory(int inum, int nall)
{
  if (nmaxatom < nall) {
  nmaxatom = nall;
  this->free_atommemory();
  this->allocate_atommemory();
  }
  nlocalatom = inum;
  nghostatom = nall - inum;
  ntotalatom = nall;
  nlocalmax = PODMAX(nlocalmax, nlocalatom);
}

void CPairPOD::check_pairmemory(double *x, double *a1, double *a2, double *a3, int natom)
{
  double rcut = podptr->pod.rcut;
  int m=0, n=0, p=0;
  if (podptr->pod.pbc[0] == 1) m = (int) ceil(rcut/a1[0]);
  if (podptr->pod.pbc[1] == 1) n = (int) ceil(rcut/a2[1]);
  if (podptr->pod.pbc[2] == 1) p = (int) ceil(rcut/a3[2]);

  // number of lattices

  int nl = (2*m+1)*(2*n+1)*(2*p+1);
  int nall = natom*nl;

  nlocalatom = natom;
  nghostatom = nall - natom;
  ntotalatom = nall;

  if (nlocalmax < natom) {
  nlocalmax = natom;
  printf("reallocate memory for pairnum and pairnumsum arrays\n");
  TemplateFree(pairnum, backend);
  TemplateFree(pairnumsum, backend);
  TemplateMalloc(&pairnum, nlocalmax, backend);
  TemplateMalloc(&pairnumsum, nlocalmax+1, backend);
  }

  if (nmaxatom < nall) {
  nmaxatom = nall;
  printf("reallocate memory for y and atomID arrays\n");
  TemplateFree(y, backend);
  TemplateFree(atomID, backend);
  TemplateMalloc(&y, dim*nmaxatom, backend);
  TemplateMalloc(&atomID, nmaxatom, backend);

  // allocate memory for atom arrays

  this->free_atommemory();
  this->allocate_atommemory();
  }

  double rcutsq = rcut*rcut;
  this->latticecoords(y, atomID, x, a1, a2, a3, rcut, podptr->pod.pbc, natom);

  natompairs = this->podneighborcount(y, rcutsq, natom, nall, dim);
  if (nmaxpairs < natompairs) {
  nmaxpairs = natompairs;
  printf("reallocate memory for pairlist arrays\n");
  TemplateFree(pairlist, backend);
  TemplateMalloc(&pairlist, nmaxpairs, backend);
  }

  // total number of neighbors

   natompairs = this->podneighborlist(pairlist, pairnum, y, rcutsq, natom, nall, dim);
   podptr->podCumsum(pairnumsum, pairnum, natom+1);
}

void CPairPOD::estimate_tempmemory()
{
  int nrbf2 = podptr->pod.nbf2;
  int nabf3 = podptr->pod.nabf3;
  int nrbf3 = podptr->pod.nrbf3;
  int ns2 = podptr->pod.ns2;
  int ns3 = podptr->pod.ns3;

  szd = dim*nijmax+ (1+dim)*nijmax*PODMAX(nrbf2+ns2,nrbf3+ns3) + (nabf3+1)*7;
  int szsnap = 0;
  if (podptr->sna.twojmax>0) {
  szsnap += nijmax*dim;
  szsnap += PODMAX(2*podptr->sna.idxu_max*nijmax, 2*podptr->sna.idxz_max*podptr->sna.ndoubles*nablockmax); // (Ur, Ui) and (Zr, Zi)
  szsnap += 2*podptr->sna.idxu_max*dim*nijmax; // dUr, dUi
  szsnap += PODMAX(podptr->sna.idxb_max*podptr->sna.ntriples*dim*nijmax, 2*podptr->sna.idxu_max*podptr->sna.nelements*nablockmax); // dblist and (Utotr, Utoti)
  }

  szd = PODMAX(szsnap, szd);
  szd = nablockmax*(podptr->pod.nd1234) + szd;
}

void CPairPOD::check_tempmemory(int start, int end)
{
  nablock = end - start;
  nij = 0;
  for (int ii=0; ii<nablock; ii++) {
  int gi = start + ii;
  nij += pairnumsum[gi+1] - pairnumsum[gi];
  }

  if ( (nij > nijmax) || (nablock > nablockmax) ) {
  nijmax = PODMAX(nijmax, nij);
  nablockmax = PODMAX(nablockmax, nablock);
  this->estimate_tempmemory();
  this->free_tempmemory();
  this->allocate_tempmemory();
  }
}

void CPairPOD::podNeighPairs(int *atomtypes, int start, int end)
{
  this->check_tempmemory(start, end);

  nablock = end - start;
  int k = 0;

  // loop over atoms ini a simulation block, used for GPU

  for (int ii=0; ii<nablock; ii++) {  
    int gi = start + ii; // atom i
    int itype = atomtypes[gi];
    int s = pairnumsum[gi];
    int m = pairnumsum[gi+1] - s;
    typeai[ii] = itype;
    numneighsum[ii+1] = m;
    for (int l=0; l<m ; l++) {
      int gj = pairlist[s + l]; // atom j
      idxi[k]  = ii;
      ai[k]  = atomID[gi];
      aj[k]  = atomID[gj];
      ti[k]  = itype;
      tj[k]  = atomtypes[aj[k]];
      rij[k*3+0]   = y[gj*3+0] -  y[gi*3+0];  // xj - xi
      rij[k*3+1]   = y[gj*3+1] -  y[gi*3+1];  // xj - xi
      rij[k*3+2]   = y[gj*3+2] -  y[gi*3+2];  // xj - xi
      k += 1;
    }
  }

  numneighsum[0] = 0;
  for (int ii=0; ii<nablock; ii++)
    numneighsum[ii+1] = numneighsum[ii+1] + numneighsum[ii];
}

double CPairPOD::podenergy(double *x, double *a1, double *a2, double *a3, int *atomtypes, int inum)
{
  // determine computation blocks

  this->get_atomblocks(inum);

  // check and allocate memory for atom/pair arrays, and create full neighbor list

  this->check_pairmemory(x, a1, a2, a3, inum);

  // initialize global descriptors to zero

  int nd1234 = podptr->pod.nd1234;
  podptr->podArraySetValue(gd, 0.0, nd1234);

  for (int i = 0; i< numblocks; i++) {

  // number of atoms in this computation block

  int nat = atomblocks[i+1] - atomblocks[i];

  // get POD neighbor pairs for this computation block

  podNeighPairs(atomtypes, atomblocks[i], atomblocks[i+1]);

  // compute global POD descriptors for this computation block

  podptr->linear_descriptors_ij(gd, tmpmem, rij, &tmpmem[nat*nd1234], numneighsum,
    typeai, idxi, ti, tj, nat, nij);

  }

  // compute energy and effective coefficients

  energy = podptr->calculate_energy(energycoeff, forcecoeff, gd, podcoeff);

  return energy;
}

double CPairPOD::podeatom(double *eatom, double *x, double *a1, double *a2, double *a3, int *atomtypes, int inum)
{
  int nd1234 = podptr->pod.nd1234;

  // compute energy and effective coefficients

  energy = this->podenergy(x, a1, a2, a3, atomtypes, inum);

  // initialize force to zero

  podptr->podArraySetValue(eatom, 0.0, inum);

  for (int i = 0; i< numblocks; i++) { // loop over each computation block

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get POD neighbor pairs for this computation block

    podNeighPairs(atomtypes, atomblocks[i], atomblocks[i+1]);

    // compute global POD descriptors for this computation block

    podptr->linear_descriptors_ij(gd, tmpmem, rij, &tmpmem[nat*nd1234], numneighsum,
      typeai, idxi, ti, tj, nat, nij);

    // calculate eatom = ld * energycoeff

    char chn = 'N';
    double one = 1.0, zero = 0.0;
    int inc1 = 1;
    DGEMV(&chn, &nat, &nd1234, &one, tmpmem, &nat, energycoeff, &inc1, &zero, &eatom[atomblocks[i]], &inc1);
  }

  return energy;
}

void CPairPOD::podforce(double *f, double *x, double *a1, double *a2, double *a3, int *atomtypes, int inum)
{
  // initialize force to zero

  podptr->podArraySetValue(f, 0.0, dim*inum);

  for (int i = 0; i< numblocks; i++) { // loop over each computation block

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get POD neighbor pairs for this computation block

    podNeighPairs(atomtypes, atomblocks[i], atomblocks[i+1]);

    // compute atomic force for this computation block

    podptr->calculate_force(f, forcecoeff, rij, tmpmem, numneighsum,
      typeai, idxi, ai, aj, ti, tj, nat, nij);
    }
}

double CPairPOD::podenergyforce(double *f, double *x, double *a1, double *a2, double *a3, int *atomtypes, int inum)
{
  // compute energy and effective coefficients

  energy = this->podenergy(x, a1, a2, a3, atomtypes, inum);

  // initialize force to zero

  podptr->podArraySetValue(f, 0.0, dim*inum);

  for (int i = 0; i< numblocks; i++) { // loop over each computation block

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get POD neighbor pairs for this computation block

    podNeighPairs(atomtypes, atomblocks[i], atomblocks[i+1]);

    // compute atomic force for this computation block

    podptr->calculate_force(f, forcecoeff, rij, tmpmem, numneighsum,
      typeai, idxi, ai, aj, ti, tj, nat, nij);
    }

  return energy;
}

void CPairPOD::lammpsNeighPairs(double **x, int **firstneigh, int *atomtypes, int *numneigh, int gi)
{
  double rcutsq = podptr->pod.rcut*podptr->pod.rcut;

  nij = 0;
  int itype = atomtypes[gi];
  int m = numneigh[gi];
  typeai[0] = itype;
  for (int l=0; l<m ; l++) {   // loop over each atom around atom i
    int gj = firstneigh[gi][l];  // atom j
    double delx   = x[gj][0] -  x[gi][0];  // xj - xi
    double dely   = x[gj][1] -  x[gi][1];  // xj - xi
    double delz   = x[gj][2] -  x[gi][2];  // xj - xi
    double rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rcutsq && rsq > 1e-20) {
      rij[nij*3 + 0] = delx;
      rij[nij*3 + 1] = dely;
      rij[nij*3 + 2] = delz;
      idxi[nij]    = 0;
      ai[nij]    = gi;
      aj[nij]    = gj;
      ti[nij]    = itype;
      tj[nij]    = atomtypes[gj];
      nij++;
    }
  }

  numneighsum[0] = 0;
  numneighsum[1] = nij;
}

void CPairPOD::check_tempmemory(double **x, int **firstneigh, int *numneigh, int *ilist, int start, int end)
{
  double rcutsq = podptr->pod.rcut*podptr->pod.rcut;
  nablock = end - start;
  nij = 0;
  for (int ii=0; ii<nablock; ii++) {  // for each atom i in the simulation box
    int gi = ilist[start+ii];   // atom i
    int m = numneigh[gi];
    for (int l=0; l<m ; l++) {   // loop over each atom around atom i
      int gj = firstneigh[gi][l];  // atom j
      double delx   = x[gj][0] -  x[gi][0];  // xj - xi
      double dely   = x[gj][1] -  x[gi][1];  // xj - xi
      double delz   = x[gj][2] -  x[gi][2];  // xj - xi
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rcutsq && rsq > 1e-20) nij++;
    }
  }

  if ( (nij > nijmax) || (nablock > nablockmax) ) {
    nijmax = PODMAX(nijmax, nij);
    nablockmax = PODMAX(nablockmax, nablock);
    this->estimate_tempmemory();
    this->free_tempmemory();
    this->allocate_tempmemory();
  }
}

void CPairPOD::lammpsNeighPairs(double **x, int **firstneigh, int *atomtypes, int *numneigh,
  int *ilist, int start, int end)
{
  this->check_tempmemory(x, firstneigh, numneigh, ilist, start, end);

  nablock = end - start;
  double rcutsq = podptr->pod.rcut*podptr->pod.rcut;

  nij = 0;
  for (int ii=0; ii<nablock; ii++) {  // for each atom i in the simulation box
    int gi = ilist[start+ii];   // atom i
    int itype = atomtypes[gi];
    int m = numneigh[gi];
    numneighsum[ii+1] = 0;
    typeai[ii] = itype;
    for (int l=0; l<m ; l++) {   // loop over each atom around atom i
      int gj = firstneigh[gi][l];  // atom j
      double delx   = x[gj][0] -  x[gi][0];  // xj - xi
      double dely   = x[gj][1] -  x[gi][1];  // xj - xi
      double delz   = x[gj][2] -  x[gi][2];  // xj - xi
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rcutsq && rsq > 1e-20) {
        rij[nij*3 + 0] = delx;
        rij[nij*3 + 1] = dely;
        rij[nij*3 + 2] = delz;
        idxi[nij]  = ii;
        ai[nij]  = gi;
        aj[nij]  = gj;
        ti[nij]  = itype;
        tj[nij]  = atomtypes[gj];
        nij++;
        numneighsum[ii+1] += 1;
      }
    }
  }

  numneighsum[0] = 0;
  for (int ii=0; ii<nablock; ii++)
    numneighsum[ii+1] = numneighsum[ii+1] + numneighsum[ii];
}

double CPairPOD::lammpsenergy(double **x, int **firstneigh, int *atomtypes, int *numneigh,
  int *ilist, int inum, int nall)
{
  // determine computation blocks

  this->get_atomblocks(inum);

  // check atom memory

  this->check_atommemory(inum, nall);

  int nd1234 = podptr->pod.nd1234;
  podptr->podArraySetValue(gd, 0.0, nd1234);

  // loop over computation blocks

  for (int i = 0; i< numblocks; i++) {

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get LAMMPS neighbor pairs for this computation block

    lammpsNeighPairs(x, firstneigh, atomtypes, numneigh, ilist, atomblocks[i], atomblocks[i+1]);

    // compute global POD descriptors for this computation block

    podptr->linear_descriptors_ij(gd, tmpmem, rij, &tmpmem[nat*nd1234], numneighsum,
      typeai, idxi, ti, tj, nat, nij);
  }

  // compute energy and effective coefficients

  energy = podptr->calculate_energy(energycoeff, forcecoeff, gd, podcoeff);

  return energy;
}

double CPairPOD::lammpseatom(double *eatom, double **x, int **firstneigh, int *atomtypes, int *numneigh,
  int *ilist, int inum, int nall)
{
  // compute energy and effective coefficients

  energy = this->lammpsenergy(x, firstneigh, atomtypes, numneigh, ilist, inum, nall);

  int nd1234 = podptr->pod.nd1234;
  podptr->podArraySetValue(gd, 0.0, nd1234);

  for (int i = 0; i< numblocks; i++) { // loop over each computation block

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    double *localdesc = tmpmem;
    double *ea = &tmpmem[nat*nd1234];

    // get LAMMPS neighbor pairs for this computation block

    lammpsNeighPairs(x, firstneigh, atomtypes, numneigh, ilist, atomblocks[i], atomblocks[i+1]);

    // compute global POD descriptors for this computation block

    podptr->linear_descriptors_ij(gd, localdesc, rij, ea, numneighsum,
      typeai, idxi, ti, tj, nat, nij);

    // calculate eatom = localdesc * energycoeff

    char chn = 'N';
    double one = 1.0, zero = 0.0;
    int inc1 = 1;
    DGEMV(&chn, &nat, &nd1234, &one, localdesc, &nat, energycoeff, &inc1, &zero, ea, &inc1);

    for (int j = 0; j<nat; j++)
      eatom[ilist[atomblocks[i] + j]] = ea[j];

  }

  return energy;
}

void CPairPOD::lammpsforce(double **f, double **x, int **firstneigh, int *atomtypes,
  int *numneigh, int *ilist, int inum, int nall)
{
  podptr->podArraySetValue(forces, 0.0, dim*nall);

  // loop over computation blocks

  for (int i = 0; i< numblocks; i++) {

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get LAMMPS neighbor pairs for this computation block

    lammpsNeighPairs(x, firstneigh, atomtypes, numneigh, ilist, atomblocks[i], atomblocks[i+1]);

    // compute atomic force for this computation block

    podptr->calculate_force(forces, forcecoeff, rij, tmpmem, numneighsum,
      typeai, idxi, ai, aj, ti, tj, nat, nij);
  }

  // copy force to lammps force array

  for (int i = 0; i<nall; i++) {
    f[i][0] = forces[0+3*i];
    f[i][1] = forces[1+3*i];
    f[i][2] = forces[2+3*i];
  }
}

double CPairPOD::lammpsenergyforce(double **f, double **x, int **firstneigh, int *atomtypes,
  int *numneigh, int *ilist, int inum, int nall)
{

  // compute energy and effective coefficients

  energy = this->lammpsenergy(x, firstneigh, atomtypes, numneigh, ilist, inum, nall);

  podptr->podArraySetValue(forces, 0.0, dim*nall);

  // loop over computation blocks

  for (int i = 0; i< numblocks; i++) {

    // # of atoms in this computation block

    int nat = atomblocks[i+1] - atomblocks[i];

    // get LAMMPS neighbor pairs for this computation block

    lammpsNeighPairs(x, firstneigh, atomtypes, numneigh, ilist, atomblocks[i], atomblocks[i+1]);

    // compute atomic force for this computation block

    podptr->calculate_force(forces, forcecoeff, rij, tmpmem, numneighsum,
      typeai, idxi, ai, aj, ti, tj, nat, nij);
  }

  // copy force to lammps force array

  for (int i = 0; i<nall; i++) {
    f[i][0] = forces[0+3*i];
    f[i][1] = forces[1+3*i];
    f[i][2] = forces[2+3*i];
  }

  return energy;
}

void CPairPOD::error_analsysis()
{
  double f[dim*data.num_atom_max];

  double Rmat[dim*dim]; // Rotation matrix
  double alpha = 38.0*M_PI/180.0;
  double beta = 62.0*M_PI/180.0;
  double gamma = 123.0*M_PI/180.0;
  if (randomrotation==1) {
    podptr->rotation_matrix(Rmat, alpha, beta, gamma);
    podptr->print_matrix("Rotation matrix", dim, dim, Rmat, dim);
  }

  double tricliniclattice[dim*dim]; // triclinic lattice vectors
  double *a = &tricliniclattice[0];
  double *b = &tricliniclattice[dim];
  double *c = &tricliniclattice[2*dim];
  double invmat[dim*dim];
  double Qmat[dim*dim];

  int nfiles = data.data_files.size();  // number of files
  int num_configs = data.num_atom.size(); // number of configurations in all files

  int m = 8;
  double outarray[m*num_configs];
  double errors[4*(nfiles+1)];
  for (int i=0; i<4*(nfiles+1); i++)
    errors[i] = 0.0;

  std::cout<<"**************** Begin of Error Calculation ****************"<<std::endl;

  int ci = 0; // configuration counter
  int nc = 0, nf = 0;

  // loop over files in training set

  for (int file = 0; file < nfiles; file++) {
    double emae=0.0, essr=0.0, fmae=0.0, fssr=0.0;
    int nforceall = 0;

    int nconfigs = data.num_config[file];
    nc += nconfigs;

    // loop over configurations in a file

    for (int ii=0; ii < nconfigs; ii++) {
      if ((ci % 100)==0) std::cout<<"Configuration: # "<<ci+1<<std::endl;

      int natom = data.num_atom[ci];
      int nforce = dim*natom;
      int natom_cumsum2 = data.num_atom_cumsum[ci];
      int *atomtypes = &data.atomtype[natom_cumsum2];
      double *x = &data.position[dim*natom_cumsum2];
      double *lattice = &data.lattice[9*ci];
      double *a1 = &lattice[0];
      double *a2 = &lattice[3];
      double *a3 = &lattice[6];

      if (randomrotation==1) {
        podptr->matrix33_multiplication(x, Rmat, x, natom);
        podptr->matrix33_multiplication(a1, Rmat, a1, 1);
        podptr->matrix33_multiplication(a2, Rmat, a2, 1);
        podptr->matrix33_multiplication(a3, Rmat, a3, 1);
      }

      // convert the atom structure to triclinic system

      podptr->triclinic_lattice_conversion(a, b, c, a1, a2, a3);
      podptr->matrix33_inverse(invmat, a1, a2, a3);
      podptr->matrix33_multiplication(Qmat, tricliniclattice, invmat, dim);
      podptr->matrix33_multiplication(x, Qmat, x, natom);

      // triclinic system is required

      energy = this->podenergyforce(f, x, a, b, c, atomtypes, natom);

      double DFTenergy = data.energy[ci];
      int natom_cumsum = data.num_atom_cumsum[ci];
      double *DFTforce = &data.force[dim*natom_cumsum];

      if (randomrotation==1)
        podptr->matrix33_multiplication(DFTforce, Rmat, DFTforce, natom);

      // convert the atom structure to triclinic system

      podptr->matrix33_multiplication(DFTforce, Qmat, DFTforce, natom);

      outarray[0 + m*ci] = ci+1;
      outarray[1 + m*ci] = natom;
      outarray[2 + m*ci] = energy;
      outarray[3 + m*ci] = DFTenergy;
      outarray[4 + m*ci] = fabs(DFTenergy-energy)/natom;
      outarray[5 + m*ci] = podptr->podArrayNorm(f, nforce);
      outarray[6 + m*ci] = podptr->podArrayNorm(DFTforce, nforce);

      double diff, sum = 0.0, ssr = 0.0;
      for (int j=0; j<dim*natom; j++) {
        diff = DFTforce[j] - f[j];
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

  std::cout<<"**************** End of Error Calculation ****************"<<std::endl<<std::endl;

  this->print_analysis(outarray, errors);
}

void CPairPOD::print_analysis(double *outarray, double *errors)
{
  string s = "All files";
  int nfiles = data.data_files.size();  // number of files
  int lm = s.size();
  for (int i = 0; i < nfiles; i++)
    lm = PODMAX(lm, (int) data.filenames[i].size());
  lm = lm + 2;

  std::string filename = "errors.txt";
  std::ofstream myfile (filename);
  if (!myfile.is_open()) std::cout << "Unable to open file";

  filename = "analysis.txt";
  std::ofstream mfile (filename);
  if (!mfile.is_open()) std::cout << "Unable to open file";

  std::string sa = "**************** Begin of Error Analysis for the Data Set ****************";
  std::string mystr = sa;
  std::cout<<mystr<<std::endl;
  myfile <<mystr+"\n";

  sa = "----------------------------------------------------------------------------------------\n";
  std::string sb = "  File  | # configs | # atoms  | MAE energy | RMSE energy | MAE force | RMSE force |\n";
  std::cout<<sa; myfile <<sa;
  std::cout<<sb; myfile <<sb;
  std::cout<<sa; myfile <<sa;

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
        mfile << outarray[count + m*ci] << "  ";
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
    std::cout<<s;
    myfile <<s;
  }
  std::cout<<sa; myfile <<sa;

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
  std::cout<<s; myfile <<s;
  std::cout<<sa; myfile <<sa;

  sa = "**************** End of Error Analysis for the Data Set ****************";
  mystr = sa;
  std::cout<<mystr<<std::endl;
  myfile <<mystr+"\n";
  myfile.close();
  mfile.close();
}
