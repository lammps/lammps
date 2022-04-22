#include <iostream>
#include <string.h>
#include <iomanip>
#include <limits>
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "output.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "fix.h"
#include "citeme.h"
#ifdef USE_TTM
#include "fix_ttm_mod.h"
#endif

#include "pair_deepmd.h"

using namespace LAMMPS_NS;
using namespace std;

static const char cite_user_deepmd_package[] =
	"USER-DEEPMD package:\n\n"
    "@article{Wang_ComputPhysCommun_2018_v228_p178,\n"
    "  author = {Wang, Han and Zhang, Linfeng and Han, Jiequn and E, Weinan},\n"
    "  doi = {10.1016/j.cpc.2018.03.016},\n"
    "  url = {https://doi.org/10.1016/j.cpc.2018.03.016},\n"
    "  year = 2018,\n"
    "  month = {jul},\n"
    "  publisher = {Elsevier {BV}},\n"
    "  volume = 228,\n"
    "  journal = {Comput. Phys. Commun.},\n"
    "  title = {{DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics}},\n"
    "  pages = {178--184}\n"
	"}\n\n";


static int stringCmp(const void *a, const void* b)
{
    char* m = (char*)a;
    char* n = (char*)b;
    int i, sum = 0;

    for(i = 0; i < MPI_MAX_PROCESSOR_NAME; i++)
        if (m[i] == n[i])
            continue;
        else
        {
            sum = m[i] - n[i];
            break;
        }
    return sum;
}

int PairDeepMD::get_node_rank() {
    char host_name[MPI_MAX_PROCESSOR_NAME];
    memset(host_name, '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    char (*host_names)[MPI_MAX_PROCESSOR_NAME];
    int n, namelen, color, rank, nprocs, myrank;
    size_t bytes;
    MPI_Comm nodeComm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Get_processor_name(host_name,&namelen);

    bytes = nprocs * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
    host_names = (char (*)[MPI_MAX_PROCESSOR_NAME]) malloc(bytes);
    for (int ii = 0; ii < nprocs; ii++) {
        memset(host_names[ii], '\0', sizeof(char) * MPI_MAX_PROCESSOR_NAME);
    }
    
    strcpy(host_names[rank], host_name);

    for (n=0; n<nprocs; n++)
        MPI_Bcast(&(host_names[n]),MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n, MPI_COMM_WORLD);
    qsort(host_names, nprocs,  sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);

    color = 0;
    for (n=0; n<nprocs-1; n++)
    {
        if(strcmp(host_name, host_names[n]) == 0)
        {
            break;
        }
        if(strcmp(host_names[n], host_names[n+1]))
        {
            color++;
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &nodeComm);
    MPI_Comm_rank(nodeComm, &myrank);

    MPI_Barrier(MPI_COMM_WORLD);
    int looprank=myrank;
    // printf (" Assigning device %d  to process on node %s rank %d, OK\n",looprank,  host_name, rank );
    free(host_names);
    return looprank;
}

std::string PairDeepMD::get_file_content(const std::string & model) {
  int myrank = 0, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int nchar = 0;
  std::string file_content;
  if (myrank == root) {
    deepmd::check_status(tensorflow::ReadFileToString(tensorflow::Env::Default(), model, &file_content));
    nchar = file_content.size();
  }
  MPI_Bcast(&nchar, 1, MPI_INT, root, MPI_COMM_WORLD);  
  char * buff = (char *)malloc(sizeof(char) * nchar);  
  if (myrank == root) {  
    memcpy(buff, file_content.c_str(), sizeof(char) * nchar);
  }
  MPI_Bcast(buff, nchar, MPI_CHAR, root, MPI_COMM_WORLD);
  file_content.resize(nchar);
  for (unsigned ii = 0; ii < nchar; ++ii) {
    file_content[ii] = buff[ii];
  }
  free(buff);
  return file_content;
}

std::vector<std::string> PairDeepMD::get_file_content(const std::vector<std::string> & models) {
  std::vector<std::string> file_contents(models.size());
  for (unsigned ii = 0; ii < models.size(); ++ii) {
    file_contents[ii] = get_file_content(models[ii]);
  }
  return file_contents;
}

static void 
ana_st (double & max, 
	double & min, 
	double & sum, 
	const vector<double> & vec, 
	const int & nloc) 
{
  if (nloc == 0) return;
  max = vec[0];
  min = vec[0];
  sum = vec[0];
  for (unsigned ii = 1; ii < nloc; ++ii){
    if (vec[ii] > max) max = vec[ii];
    if (vec[ii] < min) min = vec[ii];
    sum += vec[ii];
  }
}

static void 
make_uniform_aparam(
#ifdef HIGH_PREC    
    vector<double > & daparam,
    const vector<double > & aparam,
    const int & nlocal
#else
    vector<float > & daparam,
    const vector<float > & aparam,
    const int & nlocal
#endif
    )
{
  unsigned dim_aparam = aparam.size();
  daparam.resize(dim_aparam * nlocal);
  for (int ii = 0; ii < nlocal; ++ii){
    for (int jj = 0; jj < dim_aparam; ++jj){
      daparam[ii*dim_aparam+jj] = aparam[jj];
    }
  }
}

#ifdef USE_TTM
void PairDeepMD::make_ttm_aparam(
#ifdef HIGH_PREC
    vector<double > & daparam
#else
    vector<float > & daparam
#endif
    )
{
  assert(do_ttm);
  // get ttm_fix
  const FixTTMMod * ttm_fix = NULL;
  for (int ii = 0; ii < modify->nfix; ii++) {
    if (string(modify->fix[ii]->id) == ttm_fix_id){
      ttm_fix = dynamic_cast<FixTTMMod*>(modify->fix[ii]);
    }
  }
  assert(ttm_fix);
  // modify
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  vector<int> nnodes = ttm_fix->get_nodes();
  int nxnodes = nnodes[0];
  int nynodes = nnodes[1];
  int nznodes = nnodes[2];
  double *** const T_electron = ttm_fix->get_T_electron();
  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nynodes;
  // resize daparam
  daparam.resize(nlocal);
  // loop over atoms to assign aparam
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & ttm_fix->groupbit) {
      double xscale = (x[ii][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[ii][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[ii][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;
      daparam[ii] = T_electron[ixnode][iynode][iznode];
    }
  }
}
#endif

PairDeepMD::PairDeepMD(LAMMPS *lmp) 
    : Pair(lmp)
      
{
  if (lmp->citeme) lmp->citeme->add(cite_user_deepmd_package);
  if (strcmp(update->unit_style,"metal") != 0) {
    error->all(FLERR,"Pair deepmd requires metal unit, please set it by \"units metal\"");
  }
  restartinfo = 1;
#if LAMMPS_VERSION_NUMBER>=20201130
  centroidstressflag = CENTROID_AVAIL ; // set centroidstressflag = CENTROID_AVAIL to allow the use of the centroid/stress/atom. Added by Davide Tisi
#else 
  centroidstressflag = 2 ; // set centroidstressflag = 2 to allow the use of the centroid/stress/atom. Added by Davide Tisi
#endif
  pppmflag = 1;
  respa_enable = 0;
  writedata = 0;
  cutoff = 0.;
  numb_types = 0;
  numb_models = 0;
  out_freq = 0;
  out_each = 0;
  out_rel = 0;
  out_rel_v = 0;
  eps = 0.;
  eps_v = 0.;
  scale = NULL;
  do_ttm = false;
  single_model = false;
  multi_models_mod_devi = false;
  multi_models_no_mod_devi = false;
  is_restart = false;
  // set comm size needed by this Pair
  comm_reverse = 1;

  print_summary("  ");
}

void
PairDeepMD::print_summary(const string pre) const
{
  if (comm->me == 0){
    // capture cout to a string, then call LAMMPS's utils::logmesg
    // https://stackoverflow.com/a/4043813/9567349
    std::stringstream buffer;
    std::streambuf *sbuf = std::cout.rdbuf();
    std::cout.rdbuf(buffer.rdbuf());

    cout << "Summary of lammps deepmd module ..." << endl;
    cout << pre << ">>> Info of deepmd-kit:" << endl;
    deep_pot.print_summary(pre);
    cout << pre << ">>> Info of lammps module:" << endl;
    cout << pre << "use deepmd-kit at:  " << STR_DEEPMD_ROOT << endl;
    cout << pre << "source:             " << STR_GIT_SUMM << endl;
    cout << pre << "source branch:      " << STR_GIT_BRANCH << endl;
    cout << pre << "source commit:      " << STR_GIT_HASH << endl;
    cout << pre << "source commit at:   " << STR_GIT_DATE << endl;
    cout << pre << "build float prec:   " << STR_FLOAT_PREC << endl;
    cout << pre << "build with tf inc:  " << STR_TensorFlow_INCLUDE_DIRS << endl;
    cout << pre << "build with tf lib:  " << STR_TensorFlow_LIBRARY << endl;

    std::cout.rdbuf(sbuf);
    utils::logmesg(lmp, buffer.str());
  }
}


PairDeepMD::~PairDeepMD()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
}

void PairDeepMD::compute(int eflag, int vflag)
{
  if (numb_models == 0) return;
  if (eflag || vflag) ev_setup(eflag,vflag);
  bool do_ghost = true;
  
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = 0;
  if (do_ghost) {
    nghost = atom->nghost;
  }
  int nall = nlocal + nghost;
  int newton_pair = force->newton_pair;

  vector<int > dtype (nall);
  for (int ii = 0; ii < nall; ++ii){
    dtype[ii] = type[ii] - 1;
  }  

  double dener (0);
  vector<double > dforce (nall * 3);
  vector<double > dvirial (9, 0);
  vector<double > dcoord (nall * 3, 0.);
  vector<double > dbox (9, 0) ;
#ifdef HIGH_PREC
  vector<double > daparam;
#else 
  vector<float > daparam;
#endif

  // get box
  dbox[0] = domain->h[0];	// xx
  dbox[4] = domain->h[1];	// yy
  dbox[8] = domain->h[2];	// zz
  dbox[7] = domain->h[3];	// zy
  dbox[6] = domain->h[4];	// zx
  dbox[3] = domain->h[5];	// yx

  // get coord
  for (int ii = 0; ii < nall; ++ii){
    for (int dd = 0; dd < 3; ++dd){
      dcoord[ii*3+dd] = x[ii][dd] - domain->boxlo[dd];
    }
  }

  // uniform aparam
  if (aparam.size() > 0){
    make_uniform_aparam(daparam, aparam, nlocal);
  }
  else if (do_ttm) {
#ifdef USE_TTM
    make_ttm_aparam(daparam);
#endif
  }

  // int ago = numb_models > 1 ? 0 : neighbor->ago;
  int ago = neighbor->ago;
  if (numb_models > 1) {
      if (multi_models_no_mod_devi && (out_freq > 0 && update->ntimestep % out_freq == 0)) {
          ago = 0;
      }
      else if (multi_models_mod_devi && (out_freq == 0 || update->ntimestep % out_freq != 0)) {
        ago = 0;
      }
  }
  // compute
  single_model = (numb_models == 1);
  multi_models_no_mod_devi = (numb_models > 1 && (out_freq == 0 || update->ntimestep % out_freq != 0));
  multi_models_mod_devi = (numb_models > 1 && (out_freq > 0 && update->ntimestep % out_freq == 0));
  if (do_ghost) {
    deepmd::InputNlist lmp_list (list->inum, list->ilist, list->numneigh, list->firstneigh);
    if (single_model || multi_models_no_mod_devi) {
      //cvflag_atom is the right flag for the cvatom matrix 
      if ( ! (eflag_atom || cvflag_atom) ) {      
#ifdef HIGH_PREC
  try {
	deep_pot.compute (dener, dforce, dvirial, dcoord, dtype, dbox, nghost, lmp_list, ago, fparam, daparam);
  } catch(deepmd::deepmd_exception& e) {
    error->all(FLERR, e.what());
  }
#else
	vector<float> dcoord_(dcoord.size());
	vector<float> dbox_(dbox.size());
	for (unsigned dd = 0; dd < dcoord.size(); ++dd) dcoord_[dd] = dcoord[dd];
	for (unsigned dd = 0; dd < dbox.size(); ++dd) dbox_[dd] = dbox[dd];
	vector<float> dforce_(dforce.size(), 0);
	vector<float> dvirial_(dvirial.size(), 0);
	double dener_ = 0;
  try {
	deep_pot.compute (dener_, dforce_, dvirial_, dcoord_, dtype, dbox_, nghost, lmp_list, ago, fparam, daparam);
  } catch(deepmd::deepmd_exception& e) {
    error->all(FLERR, e.what());
  }
	for (unsigned dd = 0; dd < dforce.size(); ++dd) dforce[dd] = dforce_[dd];	
	for (unsigned dd = 0; dd < dvirial.size(); ++dd) dvirial[dd] = dvirial_[dd];	
	dener = dener_;
#endif
      }
      // do atomic energy and virial
      else {
	vector<double > deatom (nall * 1, 0);
	vector<double > dvatom (nall * 9, 0);
#ifdef HIGH_PREC
	deep_pot.compute (dener, dforce, dvirial, deatom, dvatom, dcoord, dtype, dbox, nghost, lmp_list, ago, fparam, daparam);
#else 
	vector<float> dcoord_(dcoord.size());
	vector<float> dbox_(dbox.size());
	for (unsigned dd = 0; dd < dcoord.size(); ++dd) dcoord_[dd] = dcoord[dd];
	for (unsigned dd = 0; dd < dbox.size(); ++dd) dbox_[dd] = dbox[dd];
	vector<float> dforce_(dforce.size(), 0);
	vector<float> dvirial_(dvirial.size(), 0);
	vector<float> deatom_(dforce.size(), 0);
	vector<float> dvatom_(dforce.size(), 0);
	double dener_ = 0;
  try {
	deep_pot.compute (dener_, dforce_, dvirial_, deatom_, dvatom_, dcoord_, dtype, dbox_, nghost, lmp_list, ago, fparam, daparam);
  } catch(deepmd::deepmd_exception& e) {
    error->all(FLERR, e.what());
  }
	for (unsigned dd = 0; dd < dforce.size(); ++dd) dforce[dd] = dforce_[dd];	
	for (unsigned dd = 0; dd < dvirial.size(); ++dd) dvirial[dd] = dvirial_[dd];	
	for (unsigned dd = 0; dd < deatom.size(); ++dd) deatom[dd] = deatom_[dd];	
	for (unsigned dd = 0; dd < dvatom.size(); ++dd) dvatom[dd] = dvatom_[dd];	
	dener = dener_;
#endif	
	if (eflag_atom) {
	  for (int ii = 0; ii < nlocal; ++ii) eatom[ii] += deatom[ii];
	}
	// Added by Davide Tisi 2020
	// interface the atomic virial computed by DeepMD 
	// with the one used in centroid atoms
	if (cvflag_atom) {
	  for (int ii = 0; ii < nall; ++ii){
	    //vatom[ii][0] += 1.0 * dvatom[9*ii+0];
	    //vatom[ii][1] += 1.0 * dvatom[9*ii+4];
	    //vatom[ii][2] += 1.0 * dvatom[9*ii+8];
	    //vatom[ii][3] += 1.0 * dvatom[9*ii+3];
	    //vatom[ii][4] += 1.0 * dvatom[9*ii+6];
	    //vatom[ii][5] += 1.0 * dvatom[9*ii+7];
            cvatom[ii][0] += -1.0 * dvatom[9*ii+0]; // xx
            cvatom[ii][1] += -1.0 * dvatom[9*ii+4]; // yy 
            cvatom[ii][2] += -1.0 * dvatom[9*ii+8]; // zz
            cvatom[ii][3] += -1.0 * dvatom[9*ii+3]; // xy
            cvatom[ii][4] += -1.0 * dvatom[9*ii+6]; // xz
            cvatom[ii][5] += -1.0 * dvatom[9*ii+7]; // yz
            cvatom[ii][6] += -1.0 * dvatom[9*ii+1]; // yx
            cvatom[ii][7] += -1.0 * dvatom[9*ii+2]; // zx
            cvatom[ii][8] += -1.0 * dvatom[9*ii+5]; // zy
	  }
	}
      }
    }
    else if (multi_models_mod_devi) {
      vector<double > deatom (nall * 1, 0);
      vector<double > dvatom (nall * 9, 0);
      vector<vector<double>> 	all_virial;	       
#ifdef HIGH_PREC
      vector<double> 		all_energy;
      vector<vector<double>> 	all_atom_energy;
      vector<vector<double>> 	all_atom_virial;
      try {
      deep_pot_model_devi.compute(all_energy, all_force, all_virial, all_atom_energy, all_atom_virial, dcoord, dtype, dbox, nghost, lmp_list, ago, fparam, daparam);
      } catch(deepmd::deepmd_exception& e) {
        error->all(FLERR, e.what());
      }
      // deep_pot_model_devi.compute_avg (dener, all_energy);
      // deep_pot_model_devi.compute_avg (dforce, all_force);
      // deep_pot_model_devi.compute_avg (dvirial, all_virial);
      // deep_pot_model_devi.compute_avg (deatom, all_atom_energy);
      // deep_pot_model_devi.compute_avg (dvatom, all_atom_virial);
      dener = all_energy[0];
      dforce = all_force[0];
      dvirial = all_virial[0];
      deatom = all_atom_energy[0];
      dvatom = all_atom_virial[0];
#else 
      vector<float> dcoord_(dcoord.size());
      vector<float> dbox_(dbox.size());
      for (unsigned dd = 0; dd < dcoord.size(); ++dd) dcoord_[dd] = dcoord[dd];
      for (unsigned dd = 0; dd < dbox.size(); ++dd) dbox_[dd] = dbox[dd];
      vector<float> dforce_(dforce.size(), 0);
      vector<float> dvirial_(dvirial.size(), 0);
      vector<float> deatom_(dforce.size(), 0);
      vector<float> dvatom_(dforce.size(), 0);
      double dener_ = 0;
      vector<double> 		all_energy_;
      vector<vector<float>>	all_force_;
      vector<vector<float>> 	all_virial_;	       
      vector<vector<float>> 	all_atom_energy_;
      vector<vector<float>> 	all_atom_virial_;
      try {
      deep_pot_model_devi.compute(all_energy_, all_force_, all_virial_, all_atom_energy_, all_atom_virial_, dcoord_, dtype, dbox_, nghost, lmp_list, ago, fparam, daparam);
      } catch(deepmd::deepmd_exception& e) {
        error->all(FLERR, e.what());
      }
      // deep_pot_model_devi.compute_avg (dener_, all_energy_);
      // deep_pot_model_devi.compute_avg (dforce_, all_force_);
      // deep_pot_model_devi.compute_avg (dvirial_, all_virial_);
      // deep_pot_model_devi.compute_avg (deatom_, all_atom_energy_);
      // deep_pot_model_devi.compute_avg (dvatom_, all_atom_virial_);
      dener_ = all_energy_[0];
      dforce_ = all_force_[0];
      dvirial_ = all_virial_[0];
      deatom_ = all_atom_energy_[0];
      dvatom_ = all_atom_virial_[0];
      dener = dener_;
      for (unsigned dd = 0; dd < dforce.size(); ++dd) dforce[dd] = dforce_[dd];	
      for (unsigned dd = 0; dd < dvirial.size(); ++dd) dvirial[dd] = dvirial_[dd];	
      for (unsigned dd = 0; dd < deatom.size(); ++dd) deatom[dd] = deatom_[dd];	
      for (unsigned dd = 0; dd < dvatom.size(); ++dd) dvatom[dd] = dvatom_[dd];	
      all_force.resize(all_force_.size());
      for (unsigned ii = 0; ii < all_force_.size(); ++ii){
	all_force[ii].resize(all_force_[ii].size());
	for (unsigned jj = 0; jj < all_force_[ii].size(); ++jj){
	  all_force[ii][jj] = all_force_[ii][jj];
	}
      }
      all_virial.resize(all_virial_.size());
      for (unsigned ii = 0; ii < all_virial_.size(); ++ii){
        all_virial[ii].resize(all_virial_[ii].size());
        for (unsigned jj = 0; jj < all_virial_[ii].size(); ++jj){
          all_virial[ii][jj] = all_virial_[ii][jj];
        }
      }
#endif
      if (eflag_atom) {
	for (int ii = 0; ii < nlocal; ++ii) eatom[ii] += deatom[ii];
      }
	// Added by Davide Tisi 2020
	// interface the atomic virial computed by DeepMD 
	// with the one used in centroid atoms
      if (cvflag_atom) {
	for (int ii = 0; ii < nall; ++ii){
	  //vatom[ii][0] += 1.0 * dvatom[9*ii+0];
	  //vatom[ii][1] += 1.0 * dvatom[9*ii+4];
	  //vatom[ii][2] += 1.0 * dvatom[9*ii+8];
	  //vatom[ii][3] += 1.0 * dvatom[9*ii+3];
	  //vatom[ii][4] += 1.0 * dvatom[9*ii+6];
	  //vatom[ii][5] += 1.0 * dvatom[9*ii+7];
            cvatom[ii][0] += -1.0 * dvatom[9*ii+0]; // xx
            cvatom[ii][1] += -1.0 * dvatom[9*ii+4]; // yy 
            cvatom[ii][2] += -1.0 * dvatom[9*ii+8]; // zz
            cvatom[ii][3] += -1.0 * dvatom[9*ii+3]; // xy
            cvatom[ii][4] += -1.0 * dvatom[9*ii+6]; // xz
            cvatom[ii][5] += -1.0 * dvatom[9*ii+7]; // yz
            cvatom[ii][6] += -1.0 * dvatom[9*ii+1]; // yx
            cvatom[ii][7] += -1.0 * dvatom[9*ii+2]; // zx
            cvatom[ii][8] += -1.0 * dvatom[9*ii+5]; // zy
	}
      }      
      if (out_freq > 0 && update->ntimestep % out_freq == 0) {
	int rank = comm->me;
	// std force 
	if (newton_pair) {
#if LAMMPS_VERSION_NUMBER>=20220324
	  comm->reverse_comm(this);
#else
    comm->reverse_comm_pair(this);
#endif
	}
	vector<double> std_f;
#ifdef HIGH_PREC
	vector<double> tmp_avg_f;
	deep_pot_model_devi.compute_avg (tmp_avg_f, all_force);  
	deep_pot_model_devi.compute_std_f (std_f, tmp_avg_f, all_force);
	if (out_rel == 1){
	    deep_pot_model_devi.compute_relative_std_f (std_f, tmp_avg_f, eps);
	}
#else 
	vector<float> tmp_avg_f_, std_f_;
	for (unsigned ii = 0; ii < all_force_.size(); ++ii){
	  for (unsigned jj = 0; jj < all_force_[ii].size(); ++jj){
	    all_force_[ii][jj] = all_force[ii][jj];
	  }
	}
	deep_pot_model_devi.compute_avg (tmp_avg_f_, all_force_);  
	deep_pot_model_devi.compute_std_f (std_f_, tmp_avg_f_, all_force_);
	std_f.resize(std_f_.size());
	if (out_rel == 1){
	    deep_pot_model_devi.compute_relative_std_f (std_f_, tmp_avg_f_, eps);
	}
	for (int dd = 0; dd < std_f_.size(); ++dd) std_f[dd] = std_f_[dd];
#endif
	double min = numeric_limits<double>::max(), max = 0, avg = 0;
	ana_st(max, min, avg, std_f, nlocal);
	int all_nlocal = 0;
	MPI_Reduce (&nlocal, &all_nlocal, 1, MPI_INT, MPI_SUM, 0, world);
	double all_f_min = 0, all_f_max = 0, all_f_avg = 0;
	MPI_Reduce (&min, &all_f_min, 1, MPI_DOUBLE, MPI_MIN, 0, world);
	MPI_Reduce (&max, &all_f_max, 1, MPI_DOUBLE, MPI_MAX, 0, world);
	MPI_Reduce (&avg, &all_f_avg, 1, MPI_DOUBLE, MPI_SUM, 0, world);
	all_f_avg /= double(all_nlocal);
	// std energy
	vector<double > std_e;
#ifdef HIGH_PREC
	vector<double > tmp_avg_e;
	deep_pot_model_devi.compute_avg (tmp_avg_e, all_atom_energy);
	deep_pot_model_devi.compute_std_e (std_e, tmp_avg_e, all_atom_energy);
#else 
	vector<float> tmp_avg_e_, std_e_;
	deep_pot_model_devi.compute_avg (tmp_avg_e_, all_atom_energy_);
	deep_pot_model_devi.compute_std_e (std_e_, tmp_avg_e_, all_atom_energy_);
	std_e.resize(std_e_.size());
	for (int dd = 0; dd < std_e_.size(); ++dd) std_e[dd] = std_e_[dd];
#endif	
	max = avg = 0;
	min = numeric_limits<double>::max();
	ana_st(max, min, avg, std_e, nlocal);
	double all_e_min = 0, all_e_max = 0, all_e_avg = 0;
	MPI_Reduce (&min, &all_e_min, 1, MPI_DOUBLE, MPI_MIN, 0, world);
	MPI_Reduce (&max, &all_e_max, 1, MPI_DOUBLE, MPI_MAX, 0, world);
	MPI_Reduce (&avg, &all_e_avg, 1, MPI_DOUBLE, MPI_SUM, 0, world);
	all_e_avg /= double(all_nlocal);
	// std v
	std::vector<double> send_v(9 * numb_models);
	std::vector<double> recv_v(9 * numb_models);
	for(int kk = 0; kk < numb_models; ++kk){
	  for(int ii = 0; ii < 9; ++ii){
	    send_v[kk*9+ii] = all_virial[kk][ii] / double(atom->natoms);
	  }
	}
	MPI_Reduce(&send_v[0], &recv_v[0], 9 * numb_models, MPI_DOUBLE, MPI_SUM, 0, world);
#ifdef HIGH_PREC
	std::vector<std::vector<double>> all_virial_1(numb_models);
	std::vector<double> avg_virial, std_virial;
#else
	std::vector<std::vector<float>> all_virial_1(numb_models);
	std::vector<float> avg_virial, std_virial;
#endif
	for(int kk = 0; kk < numb_models; ++kk){
	  all_virial_1[kk].resize(9);
	  for(int ii = 0; ii < 9; ++ii){
	    all_virial_1[kk][ii] = recv_v[kk*9+ii];
	  }
	}	
	double all_v_min = numeric_limits<double>::max(), all_v_max = 0, all_v_avg = 0;
	if (rank == 0){
	  deep_pot_model_devi.compute_avg(avg_virial, all_virial_1);
	  deep_pot_model_devi.compute_std(std_virial, avg_virial, all_virial_1, 1);
	  if (out_rel_v == 1){
	    deep_pot_model_devi.compute_relative_std(std_virial, avg_virial, eps_v, 1);
	  }
	  for(int ii = 0; ii < 9; ++ii){
	    if(std_virial[ii] > all_v_max){
	      all_v_max = std_virial[ii];
	    }
	    if(std_virial[ii] < all_v_min){
	      all_v_min = std_virial[ii];
	    }
	    all_v_avg += std_virial[ii] * std_virial[ii];
	  }
	  all_v_avg = sqrt(all_v_avg / 9);
	}
	// // total e
	// vector<double > sum_e(numb_models, 0.);
	// MPI_Reduce (&all_energy[0], &sum_e[0], numb_models, MPI_DOUBLE, MPI_SUM, 0, world);
	if (rank == 0) {
	  // double avg_e = 0;
	  // deep_pot_model_devi.compute_avg(avg_e, sum_e);
	  // double std_e_1 = 0;
	  // deep_pot_model_devi.compute_std(std_e_1, avg_e, sum_e);	
	  fp << setw(12) << update->ntimestep 
	     << " " << setw(18) << all_v_max
	     << " " << setw(18) << all_v_min
	     << " " << setw(18) << all_v_avg
	     << " " << setw(18) << all_f_max 
	     << " " << setw(18) << all_f_min
	     << " " << setw(18) << all_f_avg;
	     // << " " << setw(18) << avg_e
	     // << " " << setw(18) << std_e_1 / all_nlocal
	}
	if (out_each == 1){
	  vector<double> std_f_all(all_nlocal);
	  // Gather std_f and tags
	  tagint *tag = atom->tag;
	  int nprocs = comm->nprocs;
	  for (int ii = 0; ii < nlocal; ii++) {
	    tagsend[ii] = tag[ii];
	    stdfsend[ii] = std_f[ii];
	  }
	  MPI_Gather(&nlocal, 1, MPI_INT, counts, 1, MPI_INT, 0, world);
	  displacements[0] = 0;
	  for (int ii = 0; ii < nprocs-1; ii++) displacements[ii+1] = displacements[ii] + counts[ii];
	  MPI_Gatherv(tagsend, nlocal, MPI_LMP_TAGINT,
	              tagrecv, counts, displacements, MPI_LMP_TAGINT, 0, world);
	  MPI_Gatherv(stdfsend, nlocal, MPI_DOUBLE,
	              stdfrecv, counts, displacements, MPI_DOUBLE, 0, world);
	  if (rank == 0) {
	    for (int dd = 0; dd < all_nlocal; ++dd) {
	      std_f_all[tagrecv[dd]-1] = stdfrecv[dd];
	    }
	    for (int dd = 0; dd < all_nlocal; ++dd) {
	      fp << " " << setw(18) << std_f_all[dd];	
	    }
	  }
	}
	if (rank == 0) {
	  fp << endl;
	}
      }
    }
    else {
      error->all(FLERR,"unknown computational branch");
    }
  }
  else {
    if (numb_models == 1) {
#ifdef HIGH_PREC
      try {
      deep_pot.compute (dener, dforce, dvirial, dcoord, dtype, dbox);
      } catch(deepmd::deepmd_exception& e) {
        error->all(FLERR, e.what());
      }
#else
      vector<float> dcoord_(dcoord.size());
      vector<float> dbox_(dbox.size());
      for (unsigned dd = 0; dd < dcoord.size(); ++dd) dcoord_[dd] = dcoord[dd];
      for (unsigned dd = 0; dd < dbox.size(); ++dd) dbox_[dd] = dbox[dd];
      vector<float> dforce_(dforce.size(), 0);
      vector<float> dvirial_(dvirial.size(), 0);
      double dener_ = 0;
      try {
      deep_pot.compute (dener_, dforce_, dvirial_, dcoord_, dtype, dbox_);
      } catch(deepmd::deepmd_exception& e) {
        error->all(FLERR, e.what());
      }
      for (unsigned dd = 0; dd < dforce.size(); ++dd) dforce[dd] = dforce_[dd];	
      for (unsigned dd = 0; dd < dvirial.size(); ++dd) dvirial[dd] = dvirial_[dd];	
      dener = dener_;      
#endif
    }
    else {
      error->all(FLERR,"Serial version does not support model devi");
    }
  }

  // get force
  for (int ii = 0; ii < nall; ++ii){
    for (int dd = 0; dd < 3; ++dd){
      f[ii][dd] += scale[1][1] * dforce[3*ii+dd];
    }
  }
  
  // accumulate energy and virial
  if (eflag) eng_vdwl += scale[1][1] * dener;
  if (vflag) {
    virial[0] += 1.0 * dvirial[0] * scale[1][1];
    virial[1] += 1.0 * dvirial[4] * scale[1][1];
    virial[2] += 1.0 * dvirial[8] * scale[1][1];
    virial[3] += 1.0 * dvirial[3] * scale[1][1];
    virial[4] += 1.0 * dvirial[6] * scale[1][1];
    virial[5] += 1.0 * dvirial[7] * scale[1][1];
  }
}


void PairDeepMD::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(scale,n+1,n+1,"pair:scale");

  for (int i = 1; i <= n; i++){
    for (int j = i; j <= n; j++){
      setflag[i][j] = 0;
      scale[i][j] = 0;
    }
  }
  for (int i = 1; i <= numb_types; ++i) {
    if (i > n) continue;
    for (int j = i; j <= numb_types; ++j) {
      if (j > n) continue;
      setflag[i][j] = 1;
      scale[i][j] = 1;
    }
  }
}


static bool 
is_key (const string& input) 
{
  vector<string> keys ;
  keys.push_back("out_freq");
  keys.push_back("out_file");
  keys.push_back("fparam");
  keys.push_back("aparam");
  keys.push_back("ttm");
  keys.push_back("atomic");
  keys.push_back("relative");
  keys.push_back("relative_v");

  for (int ii = 0; ii < keys.size(); ++ii){
    if (input == keys[ii]) {
      return true;
    }
  }
  return false;
}


void PairDeepMD::settings(int narg, char **arg)
{
  if (narg <= 0) error->all(FLERR,"Illegal pair_style command");

  vector<string> models;
  int iarg = 0;
  while (iarg < narg){
    if (is_key(arg[iarg])) {
      break;
    }
    iarg ++;
  }
  for (int ii = 0; ii < iarg; ++ii){
    models.push_back(arg[ii]);
  }
  numb_models = models.size();
  if (numb_models == 1) {
    try {
    deep_pot.init (arg[0], get_node_rank(), get_file_content(arg[0]));
    } catch(deepmd::deepmd_exception& e) {
      error->all(FLERR, e.what());
    }
    cutoff = deep_pot.cutoff ();
    numb_types = deep_pot.numb_types();
    dim_fparam = deep_pot.dim_fparam();
    dim_aparam = deep_pot.dim_aparam();
  }
  else {
    try {
    deep_pot.init (arg[0], get_node_rank(), get_file_content(arg[0]));
    deep_pot_model_devi.init(models, get_node_rank(), get_file_content(models));
    } catch(deepmd::deepmd_exception& e) {
      error->all(FLERR, e.what());
    }
    cutoff = deep_pot_model_devi.cutoff();
    numb_types = deep_pot_model_devi.numb_types();
    dim_fparam = deep_pot_model_devi.dim_fparam();
    dim_aparam = deep_pot_model_devi.dim_aparam();
    assert(cutoff == deep_pot.cutoff());
    assert(numb_types == deep_pot.numb_types());
    assert(dim_fparam == deep_pot.dim_fparam());
    assert(dim_aparam == deep_pot.dim_aparam());
  }

  out_freq = 100;
  out_file = "model_devi.out";
  out_each = 0;
  out_rel = 0;
  eps = 0.;
  fparam.clear();
  aparam.clear();
  while (iarg < narg) {
    if (! is_key(arg[iarg])) {
      error->all(FLERR,"Illegal pair_style command\nwrong number of parameters\n");
    }
    if (string(arg[iarg]) == string("out_freq")) {
      if (iarg+1 >= narg) error->all(FLERR,"Illegal out_freq, not provided");
      out_freq = atoi(arg[iarg+1]);
      iarg += 2;
    }
    else if (string(arg[iarg]) == string("out_file")) {
      if (iarg+1 >= narg) error->all(FLERR,"Illegal out_file, not provided");
      out_file = string(arg[iarg+1]);	
      iarg += 2;
    }
    else if (string(arg[iarg]) == string("fparam")) {
      for (int ii = 0; ii < dim_fparam; ++ii){
	if (iarg+1+ii >= narg || is_key(arg[iarg+1+ii])) {
	  char tmp[1024];
	  sprintf(tmp, "Illegal fparam, the dimension should be %d", dim_fparam);		  
	  error->all(FLERR, tmp);
	}
	fparam.push_back(atof(arg[iarg+1+ii]));
      }
      iarg += 1 + dim_fparam ;
    }
    else if (string(arg[iarg]) == string("aparam")) {
      for (int ii = 0; ii < dim_aparam; ++ii){
	if (iarg+1+ii >= narg || is_key(arg[iarg+1+ii])) {
	  char tmp[1024];
	  sprintf(tmp, "Illegal aparam, the dimension should be %d", dim_aparam);		  
	  error->all(FLERR, tmp);
	}
	aparam.push_back(atof(arg[iarg+1+ii]));
      }      
      iarg += 1 + dim_aparam ;
    }
    else if (string(arg[iarg]) == string("ttm")) {
#ifdef USE_TTM
      for (int ii = 0; ii < 1; ++ii){
	if (iarg+1+ii >= narg || is_key(arg[iarg+1+ii])) {
	  error->all(FLERR, "invalid ttm key: should be ttm ttm_fix_id(str)");
	}
      }	
      do_ttm = true;
      ttm_fix_id = arg[iarg+1];
      iarg += 1 + 1;
#else
      error->all(FLERR, "The deepmd-kit was compiled without support for TTM, please rebuild it with -DUSE_TTM");
#endif      
    }
    else if (string(arg[iarg]) == string("atomic")) {
      out_each = 1;
      iarg += 1;
    }
    else if (string(arg[iarg]) == string("relative")) {
      out_rel = 1;
#ifdef HIGH_PREC
      eps = atof(arg[iarg+1]);
#else
      eps = strtof(arg[iarg+1], NULL);
#endif
      iarg += 2;
    }
    else if (string(arg[iarg]) == string("relative_v")) {
      out_rel_v = 1;
#ifdef HIGH_PREC
      eps_v = atof(arg[iarg+1]);
#else
      eps_v = strtof(arg[iarg+1], NULL);
#endif
      iarg += 2;
    }
  }
  if (out_freq < 0) error->all(FLERR,"Illegal out_freq, should be >= 0");
  if (do_ttm && aparam.size() > 0) {
    error->all(FLERR,"aparam and ttm should NOT be set simultaneously");
  }
  
  if (comm->me == 0){
    if (numb_models > 1 && out_freq > 0){
      if (!is_restart) {
      fp.open (out_file);
      fp << scientific;
      fp << "#"
	 << setw(12-1) << "step" 
	 << setw(18+1) << "max_devi_v"
	 << setw(18+1) << "min_devi_v"
	 << setw(18+1) << "avg_devi_v"
	 << setw(18+1) << "max_devi_f"
	 << setw(18+1) << "min_devi_f"
	 << setw(18+1) << "avg_devi_f"
	 << endl;
      } else {
        fp.open (out_file, std::ofstream::out | std::ofstream::app);
        fp << scientific;
      }
    }
    string pre = "  ";
    cout << pre << ">>> Info of model(s):" << endl
	 << pre << "using " << setw(3) << numb_models << " model(s): ";
    if (narg == 1) {
      cout << arg[0] << " ";
    }
    else {
      for (int ii = 0; ii < models.size(); ++ii){
      	cout << models[ii] << " ";
      }
    }
    cout << endl
	 << pre << "rcut in model:      " << cutoff << endl
	 << pre << "ntypes in model:    " << numb_types << endl;
    if (dim_fparam > 0) {
      cout << pre << "using fparam(s):    " ;
      for (int ii = 0; ii < dim_fparam; ++ii){
	cout << fparam[ii] << "  " ;
      }
      cout << endl;
    }
    if (aparam.size() > 0) {
      cout << pre << "using aparam(s):    " ;
      for (int ii = 0; ii < aparam.size(); ++ii){
	cout << aparam[ii] << "  " ;
      }
      cout << endl;
    }
  }
  
  comm_reverse = numb_models * 3;
  all_force.resize(numb_models);
}

void PairDeepMD::read_restart(FILE *)
{
  is_restart = true;
}

void PairDeepMD::write_restart(FILE *)
{
  // pass
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDeepMD::coeff(int narg, char **arg)
{
  if (!allocated) {
    allocate();
  }

  int n = atom->ntypes;
  int ilo,ihi,jlo,jhi;
  ilo = 0;
  jlo = 0;
  ihi = n;
  jhi = n;
  if (narg == 2) {
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
    if (ilo != 1 || jlo != 1 || ihi != n || jhi != n) {
      error->all(FLERR,"deepmd requires that the scale should be set to all atom types, i.e. pair_coeff * *.");
    }
  }  
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      scale[i][j] = 1.0;
      if (i > numb_types || j > numb_types) {
	char warning_msg[1024];
	sprintf(warning_msg, "Interaction between types %d and %d is set with deepmd, but will be ignored.\n Deepmd model has only %d types, it only computes the mulitbody interaction of types: 1-%d.", i, j, numb_types, numb_types);
	error->warning(FLERR, warning_msg);
      }
    }
  }
}


void PairDeepMD::init_style()
{
#if LAMMPS_VERSION_NUMBER>=20220324
  neighbor->add_request(this, NeighConst::REQ_FULL);
#else
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;  
  // neighbor->requests[irequest]->newton = 2;  
#endif
  if (out_each == 1){
    int ntotal = atom->natoms;
    int nprocs = comm->nprocs;
    memory->create(counts, nprocs, "deepmd:counts");
    memory->create(displacements, nprocs, "deepmd:displacements");
    memory->create(stdfsend,ntotal,"deepmd:stdfsendall");
    memory->create(stdfrecv,ntotal,"deepmd:stdfrecvall");
    memory->create(tagsend,ntotal,"deepmd:tagsendall");
    memory->create(tagrecv,ntotal,"deepmd:tagrecvall");
  }
}


double PairDeepMD::init_one(int i, int j)
{
  if (i > numb_types || j > numb_types) {
    char warning_msg[1024];
    sprintf(warning_msg, "Interaction between types %d and %d is set with deepmd, but will be ignored.\n Deepmd model has only %d types, it only computes the mulitbody interaction of types: 1-%d.", i, j, numb_types, numb_types);
    error->warning(FLERR, warning_msg);
  }

  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  return cutoff;
}


/* ---------------------------------------------------------------------- */

int PairDeepMD::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (int dd = 0; dd < numb_models; ++dd){
      buf[m++] = all_force[dd][3*i+0];
      buf[m++] = all_force[dd][3*i+1];
      buf[m++] = all_force[dd][3*i+2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairDeepMD::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (int dd = 0; dd < numb_models; ++dd){
      all_force[dd][3*j+0] += buf[m++];
      all_force[dd][3*j+1] += buf[m++];
      all_force[dd][3*j+2] += buf[m++];
    }
  }
}

void *PairDeepMD::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cutoff;
  }
  if (strcmp(str,"scale") == 0) {
    dim = 2;
    return (void *) scale;
  }
  return NULL;
}
