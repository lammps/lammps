/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "cuda.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "universe.h"
#include "input.h"
#include "atom_masks.h"
#include "error.h"

#include "cuda_neigh_list.h"
//#include "pre_binning_cu.h"
//#include "reverse_binning_cu.h"
#include <ctime>
#include <cmath>
#include "cuda_pair_cu.h"
#include "cuda_cu.h"

using namespace LAMMPS_NS;



Cuda::Cuda(LAMMPS* lmp) : Pointers(lmp)
{
  cuda_exists = true;
  lmp->cuda = this;

  if(universe->me == 0)
    printf("# Using LAMMPS_CUDA \n");

  shared_data.me = universe->me;
  device_set = false;

  Cuda_Cuda_GetCompileSettings(&shared_data);

  if(shared_data.compile_settings.prec_glob != sizeof(CUDA_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: Global Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_glob, sizeof(CUDA_FLOAT) / 4);

  if(shared_data.compile_settings.prec_x != sizeof(X_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: X Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_x, sizeof(X_FLOAT) / 4);

  if(shared_data.compile_settings.prec_v != sizeof(V_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: V Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_v, sizeof(V_FLOAT) / 4);

  if(shared_data.compile_settings.prec_f != sizeof(F_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: F Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_f, sizeof(F_FLOAT) / 4);

  if(shared_data.compile_settings.prec_pppm != sizeof(PPPM_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: PPPM Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_pppm, sizeof(PPPM_FLOAT) / 4);

  if(shared_data.compile_settings.prec_fft != sizeof(FFT_FLOAT) / 4) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: FFT Precision: cuda %i cpp %i\n\n", shared_data.compile_settings.prec_fft, sizeof(FFT_FLOAT) / 4);

#ifdef FFT_CUFFT

  if(shared_data.compile_settings.cufft != 1) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: cufft: cuda %i cpp %i\n\n", shared_data.compile_settings.cufft, 1);

#else

  if(shared_data.compile_settings.cufft != 0) printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: cufft: cuda %i cpp %i\n\n", shared_data.compile_settings.cufft, 0);

#endif

  if(shared_data.compile_settings.arch != CUDA_ARCH)  printf("\n\n # CUDA WARNING: Compile Settings of cuda and cpp code differ! \n # CUDA WARNING: arch: cuda %i cpp %i\n\n", shared_data.compile_settings.cufft, CUDA_ARCH);

  cu_x          = 0;
  cu_v          = 0;
  cu_f          = 0;
  cu_tag        = 0;
  cu_type       = 0;
  cu_mask       = 0;
  cu_image      = 0;
  cu_xhold      = 0;
  cu_q          = 0;
  cu_rmass      = 0;
  cu_mass       = 0;
  cu_virial     = 0;
  cu_eatom      = 0;
  cu_vatom      = 0;
  cu_radius          = 0;
  cu_density          = 0;
  cu_omega          = 0;
  cu_torque          = 0;

  cu_special           = 0;
  cu_nspecial   = 0;

  cu_molecule   = 0;

  cu_x_type           = 0;
  x_type                  = 0;
  cu_v_radius          = 0;
  v_radius          = 0;
  cu_omega_rmass          = 0;
  omega_rmass          = 0;

  binned_id = 0;
  cu_binned_id  = 0;
  binned_idnew = 0;
  cu_binned_idnew = 0;

  cu_map_array = 0;

  copy_buffer = 0;
  copy_buffersize = 0;

  neighbor_decide_by_integrator = 0;
  pinned = true;

  debugdata = 0;

  finished_setup = false;
  begin_setup = false;
  finished_run = false;

  setSharedDataZero();

  uploadtime = 0;
  downloadtime = 0;
  dotiming = false;

  dotestatom = false;
  testatom = 0;
  oncpu = true;

  self_comm = 0;
  MYDBG(printf("# CUDA: Cuda::Cuda Done...\n");)
  //cCudaData<double, float, yx >
}

Cuda::~Cuda()
{

  print_timings();

  if(universe->me == 0) printf("# CUDA: Free memory...\n");

  delete cu_q;
  delete cu_x;
  delete cu_v;
  delete cu_f;
  delete cu_tag;
  delete cu_type;
  delete cu_mask;
  delete cu_image;
  delete cu_xhold;
  delete cu_mass;
  delete cu_rmass;
  delete cu_virial;
  delete cu_eng_vdwl;
  delete cu_eng_coul;
  delete cu_extent;
  delete cu_eatom;
  delete cu_vatom;
  delete cu_radius;
  delete cu_density;
  delete cu_omega;
  delete cu_torque;
  delete cu_molecule;

  delete cu_x_type;
  delete [] x_type;
  delete cu_v_radius;
  delete [] v_radius;
  delete cu_omega_rmass;
  delete [] omega_rmass;

  delete cu_debugdata;
  delete[] debugdata;

  delete cu_map_array;

  std::map<NeighList*, CudaNeighList*>::iterator p = neigh_lists.begin();

  while(p != neigh_lists.end()) {
    delete p->second;
    ++p;
  }
}

void Cuda::accelerator(int narg, char** arg)
{
  if(device_set) return;

  if(universe->me == 0)
    printf("# CUDA: Activate GPU \n");

  int* devicelist = NULL;
  int pppn = 2;

  for(int i = 0; i < narg; i++) {
    if(strcmp(arg[i], "gpu/node") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting a number after 'gpu/node' option.");

      pppn = force->inumeric(FLERR,arg[i]);
    }

    if(strcmp(arg[i], "gpu/node/special") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting number of GPUs to be used per node after keyword 'gpu/node/special'.");

      pppn = force->inumeric(FLERR,arg[i]);

      if(pppn < 1) error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting number of GPUs to be used per node after keyword 'gpu/node special'.");

      if(i + pppn == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting list of device ids after keyword 'gpu/node special'.");

      devicelist = new int[pppn];

      for(int k = 0; k < pppn; k++) {
        i++;
        devicelist[k] = force->inumeric(FLERR,arg[i]);
      }
    }

    if(strcmp(arg[i], "pinned") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting a number after 'pinned' option.");

      pinned = force->inumeric(FLERR,arg[i]) == 0 ? false : true;

      if((pinned == false) && (universe->me == 0)) printf(" #CUDA: Pinned memory is not used for communication\n");
    }

    if(strcmp(arg[i], "timing") == 0) {
      dotiming = true;
    }

    if(strcmp(arg[i], "suffix") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting a string after 'suffix' option.");

      strcpy(lmp->suffix, arg[i]);
    }

    if(strcmp(arg[i], "overlap_comm") == 0) {
      shared_data.overlap_comm = 1;
    }

    if(strcmp(arg[i], "test") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting a number after 'test' option.");

      testatom = force->numeric(FLERR,arg[i]);
      dotestatom = true;
    }

    if(strcmp(arg[i], "override/bpa") == 0) {
      if(++i == narg)
        error->all(FLERR, "Invalid Options for 'accelerator' command. Expecting a number after 'override/bpa' option.");

      shared_data.pair.override_block_per_atom = force->inumeric(FLERR,arg[i]);
    }
  }

  CudaWrapper_Init(0, (char**)0, universe->me, pppn, devicelist);
  //if(shared_data.overlap_comm)
  CudaWrapper_AddStreams(3);
  cu_x          = 0;
  cu_v          = 0;
  cu_f          = 0;
  cu_tag        = 0;
  cu_type       = 0;
  cu_mask       = 0;
  cu_image      = 0;
  cu_xhold      = 0;
  cu_q          = 0;
  cu_rmass      = 0;
  cu_mass       = 0;
  cu_virial     = 0;
  cu_eatom      = 0;
  cu_vatom      = 0;
  cu_radius            = 0;
  cu_density          = 0;
  cu_omega            = 0;
  cu_torque            = 0;

  cu_special           = 0;
  cu_nspecial   = 0;

  cu_molecule   = 0;

  cu_x_type           = 0;
  cu_v_radius          = 0;
  cu_omega_rmass          = 0;

  cu_binned_id  = 0;
  cu_binned_idnew = 0;
  device_set = true;
  allocate();
  delete devicelist;
}

void Cuda::setSharedDataZero()
{
  MYDBG(printf("# CUDA: Cuda::setSharedDataZero ...\n");)
  shared_data.atom.nlocal = 0;
  shared_data.atom.nghost = 0;
  shared_data.atom.nall = 0;
  shared_data.atom.nmax = 0;
  shared_data.atom.ntypes = 0;
  shared_data.atom.q_flag = 0;
  shared_data.atom.need_eatom = 0;
  shared_data.atom.need_vatom = 0;
  shared_data.atom.update_nmax = 1;
  shared_data.atom.update_nlocal = 1;
  shared_data.atom.update_neigh = 1;

  shared_data.pair.cudable_force = 0;
  shared_data.pair.collect_forces_later = 0;
  shared_data.pair.use_block_per_atom = 0;
  shared_data.pair.override_block_per_atom = -1;
  shared_data.pair.cut = 0;
  shared_data.pair.cutsq = 0;
  shared_data.pair.cut_inner = 0;
  shared_data.pair.cut_coul = 0;
  shared_data.pair.special_lj = 0;
  shared_data.pair.special_coul = 0;

  shared_data.pair.neighall = false;

  shared_data.pppm.cudable_force = 0;

  shared_data.buffersize = 0;
  shared_data.buffer_new = 1;
  shared_data.buffer = NULL;

  shared_data.comm.comm_phase = 0;
  shared_data.overlap_comm = 0;

  shared_data.comm.buffer = NULL;
  shared_data.comm.buffer_size = 0;
  shared_data.comm.overlap_split_ratio = 0;
  // setTimingsZero();
}

void Cuda::allocate()
{
  accelerator(0, NULL);
  MYDBG(printf("# CUDA: Cuda::allocate ...\n");)

  if(not cu_virial) {
    cu_virial    = new cCudaData<double, ENERGY_FLOAT, x > (NULL, & shared_data.pair.virial , 6);
    cu_eng_vdwl  = new cCudaData<double, ENERGY_FLOAT, x > (NULL, & shared_data.pair.eng_vdwl , 1);
    cu_eng_coul  = new cCudaData<double, ENERGY_FLOAT, x > (NULL, & shared_data.pair.eng_coul , 1);
    cu_extent          = new cCudaData<double, double, x> (extent, 6);
    shared_data.flag = CudaWrapper_AllocCudaData(sizeof(int));
    int size = 2 * CUDA_MAX_DEBUG_SIZE;
    debugdata = new int[size];
    cu_debugdata    = new cCudaData<int, int, x > (debugdata , size);
    shared_data.debugdata = cu_debugdata->dev_data();
  }

  checkResize();
  setSystemParams();
  MYDBG(printf("# CUDA: Cuda::allocate done...\n");)
}

void Cuda::setSystemParams()
{
  MYDBG(printf("# CUDA: Cuda::setSystemParams ...\n");)
  shared_data.atom.nlocal = atom->nlocal;
  shared_data.atom.nghost = atom->nghost;
  shared_data.atom.nall = atom->nlocal + atom->nghost;
  shared_data.atom.ntypes = atom->ntypes;
  shared_data.atom.q_flag = atom->q_flag;
  shared_data.atom.rmass_flag = atom->rmass_flag;
  MYDBG(printf("# CUDA: Cuda::setSystemParams done ...\n");)
}

void Cuda::setDomainParams()
{
  MYDBG(printf("# CUDA: Cuda::setDomainParams ...\n");)
  cuda_shared_domain* cu_domain = &shared_data.domain;

  cu_domain->triclinic = domain->triclinic;

  for(short i = 0; i < 3; ++i) {
    cu_domain->periodicity[i] = domain->periodicity[i];
    cu_domain->sublo[i] = domain->sublo[i];
    cu_domain->subhi[i] = domain->subhi[i];
    cu_domain->boxlo[i] = domain->boxlo[i];
    cu_domain->boxhi[i] = domain->boxhi[i];
    cu_domain->prd[i] = domain->prd[i];
  }

  if(domain->triclinic) {
    for(short i = 0; i < 3; ++i) {
      cu_domain->boxlo_lamda[i] = domain->boxlo_lamda[i];
      cu_domain->boxhi_lamda[i] = domain->boxhi_lamda[i];
      cu_domain->prd_lamda[i] = domain->prd_lamda[i];
      cu_domain->sublo[i] = domain->sublo_lamda[i];
      cu_domain->subhi[i] = domain->subhi_lamda[i];
    }

    cu_domain->xy = domain->xy;
    cu_domain->xz = domain->xz;
    cu_domain->yz = domain->yz;
  }

  for(int i = 0; i < 6; i++) {
    cu_domain->h[i] = domain->h[i];
    cu_domain->h_inv[i] = domain->h_inv[i];
    cu_domain->h_rate[i] = domain->h_rate[i];
  }

  cu_domain->update = 2;
  MYDBG(printf("# CUDA: Cuda::setDomainParams done ...\n");)
}

void Cuda::checkResize()
{
  MYDBG(printf("# CUDA: Cuda::checkResize ...\n");)
  accelerator(0, NULL);
  cuda_shared_atom* cu_atom = & shared_data.atom;
  cuda_shared_pair* cu_pair = & shared_data.pair;
  cu_atom->q_flag      = atom->q_flag;
  cu_atom->rmass_flag  = atom->rmass ? 1 : 0;
  cu_atom->nall = atom->nlocal + atom->nghost;
  cu_atom->nlocal      = atom->nlocal;
  cu_atom->nghost      = atom->nghost;

  // do we have more atoms to upload than currently allocated memory on device? (also true if nothing yet allocated)
  if(atom->nmax > cu_atom->nmax || cu_tag == NULL) {
    delete cu_x;
    cu_x         = new cCudaData<double, X_FLOAT, yx> ((double*)atom->x , & cu_atom->x        , atom->nmax, 3, 0, true); //cu_x->set_buffer(&(shared_data.buffer),&(shared_data.buffersize),true);
    delete cu_v;
    cu_v         = new cCudaData<double, V_FLOAT, yx> ((double*)atom->v, & cu_atom->v         , atom->nmax, 3);
    delete cu_f;
    cu_f         = new cCudaData<double, F_FLOAT, yx> ((double*)atom->f, & cu_atom->f         , atom->nmax, 3, 0, true);
    delete cu_tag;
    cu_tag       = new cCudaData<int   , int    , x > (atom->tag       , & cu_atom->tag       , atom->nmax, 0, true);
    delete cu_type;
    cu_type      = new cCudaData<int   , int    , x > (atom->type      , & cu_atom->type      , atom->nmax, 0, true);
    delete cu_mask;
    cu_mask      = new cCudaData<int   , int    , x > (atom->mask      , & cu_atom->mask      , atom->nmax, 0, true);
    delete cu_image;
    cu_image     = new cCudaData<int   , int    , x > (atom->image     , & cu_atom->image     , atom->nmax, 0, true);

    if(atom->rmass) {
      delete cu_rmass;
      cu_rmass     = new cCudaData<double, V_FLOAT, x > (atom->rmass     , & cu_atom->rmass     , atom->nmax);
    }

    if(cu_atom->q_flag) {
      delete cu_q;
      cu_q         = new cCudaData<double, F_FLOAT, x > ((double*)atom->q, & cu_atom->q         , atom->nmax, 0 , true);
    }// cu_q->set_buffer(&(copy_buffer),&(copy_buffersize),true);}

    if(atom->radius) {
      delete cu_radius;
      cu_radius    = new cCudaData<double, X_FLOAT, x > (atom->radius    , & cu_atom->radius     , atom->nmax);
      delete cu_v_radius;
      cu_v_radius  = new cCudaData<V_FLOAT, V_FLOAT, x> (v_radius , & cu_atom->v_radius      , atom->nmax * 4);
      delete cu_omega_rmass;
      cu_omega_rmass  = new cCudaData<V_FLOAT, V_FLOAT, x> (omega_rmass , & cu_atom->omega_rmass      , atom->nmax * 4);
    }

    if(atom->omega) {
      delete cu_omega;
      cu_omega     = new cCudaData<double, V_FLOAT, yx > (((double*) atom->omega)    , & cu_atom->omega     , atom->nmax, 3);
    }

    if(atom->torque) {
      delete cu_torque;
      cu_torque    = new cCudaData<double, F_FLOAT, yx > (((double*) atom->torque)   , & cu_atom->torque     , atom->nmax, 3);
    }

    if(atom->special) {
      delete cu_special;
      cu_special    = new cCudaData<int, int, yx > (((int*) & (atom->special[0][0]))   , & cu_atom->special     , atom->nmax, atom->maxspecial, 0 , true);
      shared_data.atom.maxspecial = atom->maxspecial;
    }

    if(atom->nspecial) {
      delete cu_nspecial;
      cu_nspecial    = new cCudaData<int, int, yx > (((int*) atom->nspecial)  , & cu_atom->nspecial     , atom->nmax, 3, 0, true);
    }

    if(atom->molecule) {
      delete cu_molecule;
      cu_molecule    = new cCudaData<int, int, x > (((int*) atom->molecule)  , & cu_atom->molecule     , atom->nmax, 0 , true);
    }

    shared_data.atom.special_flag = neighbor->special_flag;
    shared_data.atom.molecular = atom->molecular;

    cu_atom->update_nmax = 2;
    cu_atom->nmax        = atom->nmax;

    delete cu_x_type;
    cu_x_type   = new cCudaData<X_FLOAT, X_FLOAT, x> (x_type , & cu_atom->x_type      , atom->nmax * 4);
  }

  if(((cu_xhold == NULL) || (cu_xhold->get_dim()[0] < neighbor->maxhold)) && neighbor->xhold) {
    delete cu_xhold;
    cu_xhold     = new cCudaData<double, X_FLOAT, yx> ((double*)neighbor->xhold, & cu_atom->xhold         , neighbor->maxhold, 3);
    shared_data.atom.maxhold = neighbor->maxhold;
  }

  if(atom->mass && !cu_mass) {
    cu_mass      = new cCudaData<double, V_FLOAT, x > (atom->mass      , & cu_atom->mass      , atom->ntypes + 1);
  }

  cu_atom->mass_host   = atom->mass;

  if(atom->map_style == 1) {
    if((cu_map_array == NULL)) {
      cu_map_array   = new cCudaData<int, int, x > (atom->get_map_array()   , & cu_atom->map_array     , atom->get_map_size());
    } else if(cu_map_array->dev_size() / sizeof(int) < atom->get_map_size()) {
      delete cu_map_array;
      cu_map_array   = new cCudaData<int, int, x > (atom->get_map_array()   , & cu_atom->map_array     , atom->get_map_size());
    }
  }


  // if any of the host pointers have changed (e.g. re-allocated somewhere else), set to correct pointer
  if(cu_x   ->get_host_data() != atom->x)    cu_x   ->set_host_data((double*)(atom->x));

  if(cu_v   ->get_host_data() != atom->v)    cu_v   ->set_host_data((double*)(atom->v));

  if(cu_f   ->get_host_data() != atom->f)    cu_f   ->set_host_data((double*)(atom->f));

  if(cu_tag ->get_host_data() != atom->tag)  cu_tag ->set_host_data(atom->tag);

  if(cu_type->get_host_data() != atom->type) cu_type->set_host_data(atom->type);

  if(cu_mask->get_host_data() != atom->mask) cu_mask->set_host_data(atom->mask);

  if(cu_image->get_host_data() != atom->image) cu_mask->set_host_data(atom->image);

  if(cu_xhold)
    if(cu_xhold->get_host_data() != neighbor->xhold) cu_xhold->set_host_data((double*)(neighbor->xhold));

  if(atom->rmass)
    if(cu_rmass->get_host_data() != atom->rmass) cu_rmass->set_host_data((double*)(atom->rmass));

  if(cu_atom->q_flag)
    if(cu_q->get_host_data() != atom->q) cu_q->set_host_data((double*)(atom->q));

  if(atom->radius)
    if(cu_radius->get_host_data() != atom->radius) cu_radius->set_host_data((double*)(atom->radius));

  if(atom->omega)
    if(cu_omega->get_host_data() != atom->omega) cu_omega->set_host_data((double*)(atom->omega));

  if(atom->torque)
    if(cu_torque->get_host_data() != atom->torque) cu_torque->set_host_data((double*)(atom->torque));

  if(atom->special)
    if(cu_special->get_host_data() != atom->special) {
      delete cu_special;
      cu_special    = new cCudaData<int, int, yx > (((int*) atom->special)   , & cu_atom->special     , atom->nmax, atom->maxspecial);
      shared_data.atom.maxspecial = atom->maxspecial;
    }

  if(atom->nspecial)
    if(cu_nspecial->get_host_data() != atom->nspecial) cu_nspecial->set_host_data((int*)(atom->nspecial));

  if(atom->molecule)
    if(cu_molecule->get_host_data() != atom->molecule) cu_molecule->set_host_data((int*)(atom->molecule));

  if(force)
    if(cu_virial   ->get_host_data() != force->pair->virial)    cu_virial   ->set_host_data(force->pair->virial);

  if(force)
    if(cu_eng_vdwl ->get_host_data() != &force->pair->eng_vdwl)    cu_eng_vdwl  ->set_host_data(&force->pair->eng_vdwl);

  if(force)
    if(cu_eng_coul ->get_host_data() != &force->pair->eng_coul)    cu_eng_coul   ->set_host_data(&force->pair->eng_coul);

  cu_atom->update_nlocal = 2;
  MYDBG(printf("# CUDA: Cuda::checkResize done...\n");)
}

void Cuda::evsetup_eatom_vatom(int eflag_atom, int vflag_atom)
{
  if(eflag_atom) {
    if(not cu_eatom)
      cu_eatom         = new cCudaData<double, ENERGY_FLOAT, x > (force->pair->eatom, & (shared_data.atom.eatom)         , atom->nmax);  // cu_eatom->set_buffer(&(copy_buffer),&(copy_buffersize),true);}

    if(cu_eatom->get_dim()[0] != atom->nmax) {
      //delete cu_eatom;
      //cu_eatom         = new cCudaData<double, ENERGY_FLOAT, x > (force->pair->eatom, & (shared_data.atom.eatom)         , atom->nmax  );// cu_eatom->set_buffer(&(copy_buffer),&(copy_buffersize),true);}
      shared_data.atom.update_nmax = 2;
    }

    cu_eatom->set_host_data(force->pair->eatom);
    cu_eatom->memset_device(0);
  }

  if(vflag_atom) {
    if(not cu_vatom)
      cu_vatom         = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)force->pair->vatom, & (shared_data.atom.vatom)         , atom->nmax , 6);// cu_vatom->set_buffer(&(copy_buffer),&(copy_buffersize),true);}

    if(cu_vatom->get_dim()[0] != atom->nmax) {
      //delete cu_vatom;
      //cu_vatom         = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)force->pair->vatom, & (shared_data.atom.vatom)         , atom->nmax ,6 );// cu_vatom->set_buffer(&(copy_buffer),&(copy_buffersize),true);}
      shared_data.atom.update_nmax = 2;
    }

    cu_vatom->set_host_data((double*)force->pair->vatom);
    cu_vatom->memset_device(0);
  }
}

void Cuda::uploadAll()
{
  MYDBG(printf("# CUDA: Cuda::uploadAll() ... start\n");)
  my_times starttime;
  my_times endtime;

  if(atom->nmax != shared_data.atom.nmax) checkResize();

  my_gettime(CLOCK_REALTIME, &starttime);
  cu_x   ->upload();
  cu_v   ->upload();
  cu_f   ->upload();
  cu_tag ->upload();
  cu_type->upload();
  cu_mask->upload();
  cu_image->upload();

  if(shared_data.atom.q_flag) cu_q    ->upload();

  if(atom->rmass)             cu_rmass->upload();

  if(atom->radius)            cu_radius->upload();

  if(atom->omega)             cu_omega->upload();

  if(atom->torque)            cu_torque->upload();

  if(atom->special)           cu_special->upload();

  if(atom->nspecial)          cu_nspecial->upload();

  if(atom->molecule)          cu_molecule->upload();

  if(cu_eatom) cu_eatom->upload();

  if(cu_vatom) cu_vatom->upload();

  my_gettime(CLOCK_REALTIME, &endtime);
  uploadtime += (endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
  CUDA_IF_BINNING(Cuda_PreBinning(& shared_data);)
  CUDA_IF_BINNING(Cuda_Binning(& shared_data);)

  shared_data.atom.triggerneighsq = neighbor->triggersq;
  MYDBG(printf("# CUDA: Cuda::uploadAll() ... end\n");)
}

void Cuda::downloadAll()
{
  MYDBG(printf("# CUDA: Cuda::downloadAll() ... start\n");)
  my_times starttime;
  my_times endtime;

  if(atom->nmax != shared_data.atom.nmax) checkResize();

  CUDA_IF_BINNING(Cuda_ReverseBinning(& shared_data);)
  my_gettime(CLOCK_REALTIME, &starttime);
  cu_x   ->download();
  cu_v   ->download();
  cu_f   ->download();
  cu_type->download();
  cu_tag ->download();
  cu_mask->download();
  cu_image->download();

  //if(shared_data.atom.need_eatom) cu_eatom->download();
  //if(shared_data.atom.need_vatom) cu_vatom->download();

  if(shared_data.atom.q_flag) cu_q    ->download();

  if(atom->rmass)             cu_rmass->download();

  if(atom->radius)            cu_radius->download();

  if(atom->omega)             cu_omega->download();

  if(atom->torque)            cu_torque->download();

  if(atom->special)           cu_special->download();

  if(atom->nspecial)          cu_nspecial->download();

  if(atom->molecule)          cu_molecule->download();

  if(cu_eatom) cu_eatom->download();

  if(cu_vatom) cu_vatom->download();

  my_gettime(CLOCK_REALTIME, &endtime);
  downloadtime += (endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
  MYDBG(printf("# CUDA: Cuda::downloadAll() ... end\n");)
}

void Cuda::upload(int datamask)
{
  MYDBG(printf("# CUDA: Cuda::upload() ... start\n");)
  my_times starttime;
  my_times endtime;

  if(atom->nmax != shared_data.atom.nmax) checkResize();

  my_gettime(CLOCK_REALTIME, &starttime);
  if(X_MASK & datamask) cu_x   ->upload();
  if(V_MASK & datamask) cu_v   ->upload();
  if(F_MASK & datamask) cu_f   ->upload();
  if(TYPE_MASK & datamask) cu_type->upload();
  if(TAG_MASK & datamask) cu_tag ->upload();
  if(MASK_MASK & datamask) cu_mask->upload();
  if(IMAGE_MASK & datamask) cu_image->upload();

  //if(shared_data.atom.need_eatom) cu_eatom->upload();
  //if(shared_data.atom.need_vatom) cu_vatom->upload();

  if(shared_data.atom.q_flag)
	  if(Q_MASK & datamask) cu_q    ->upload();

  if(atom->rmass)
	  if(RMASS_MASK & datamask) cu_rmass->upload();

  if(atom->radius)
	  if(RADIUS_MASK & datamask) cu_radius->upload();

  if(atom->omega)
	  if(OMEGA_MASK & datamask) cu_omega->upload();

  if(atom->torque)
	  if(TORQUE_MASK & datamask) cu_torque->upload();

  if(atom->special)
	  if(SPECIAL_MASK & datamask) cu_special->upload();

  if(atom->nspecial)
	  if(SPECIAL_MASK & datamask) cu_nspecial->upload();

  if(atom->molecule)
	  if(MOLECULE_MASK & datamask) cu_molecule->upload();

  if(cu_eatom) cu_eatom->upload();

  if(cu_vatom) cu_vatom->upload();

  my_gettime(CLOCK_REALTIME, &endtime);
  uploadtime += (endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
  MYDBG(printf("# CUDA: Cuda::upload() ... end\n");)
}

void Cuda::download(int datamask)
{
  MYDBG(printf("# CUDA: Cuda::download() ... start\n");)
  my_times starttime;
  my_times endtime;

  if(atom->nmax != shared_data.atom.nmax) checkResize();

  CUDA_IF_BINNING(Cuda_ReverseBinning(& shared_data);)
  my_gettime(CLOCK_REALTIME, &starttime);
  if(X_MASK & datamask) cu_x   ->download();
  if(V_MASK & datamask) cu_v   ->download();
  if(F_MASK & datamask) cu_f   ->download();
  if(TYPE_MASK & datamask) cu_type->download();
  if(TAG_MASK & datamask) cu_tag ->download();
  if(MASK_MASK & datamask) cu_mask->download();
  if(IMAGE_MASK & datamask) cu_image->download();

  //if(shared_data.atom.need_eatom) cu_eatom->download();
  //if(shared_data.atom.need_vatom) cu_vatom->download();

  if(shared_data.atom.q_flag)
	  if(Q_MASK & datamask) cu_q    ->download();

  if(atom->rmass)
	  if(RMASS_MASK & datamask) cu_rmass->download();

  if(atom->radius)
	  if(RADIUS_MASK & datamask) cu_radius->download();

  if(atom->omega)
	  if(OMEGA_MASK & datamask) cu_omega->download();

  if(atom->torque)
	  if(TORQUE_MASK & datamask) cu_torque->download();

  if(atom->special)
	  if(SPECIAL_MASK & datamask) cu_special->download();

  if(atom->nspecial)
	  if(SPECIAL_MASK & datamask) cu_nspecial->download();

  if(atom->molecule)
	  if(MOLECULE_MASK & datamask) cu_molecule->download();

  if(cu_eatom) cu_eatom->download();

  if(cu_vatom) cu_vatom->download();

  my_gettime(CLOCK_REALTIME, &endtime);
  downloadtime += (endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000);
  MYDBG(printf("# CUDA: Cuda::download() ... end\n");)
}

void Cuda::downloadX()
{
  Cuda_Pair_RevertXType(& this->shared_data);
  cu_x->download();
}

CudaNeighList* Cuda::registerNeighborList(class NeighList* neigh_list)
{
  MYDBG(printf("# CUDA: Cuda::registerNeighborList() ... start a\n");)
  std::map<NeighList*, CudaNeighList*>::iterator p = neigh_lists.find(neigh_list);

  if(p != neigh_lists.end()) return p->second;
  else {
    CudaNeighList* neigh_list_cuda = new CudaNeighList(lmp, neigh_list);
    neigh_lists.insert(std::pair<NeighList*, CudaNeighList*>(neigh_list, neigh_list_cuda));
    return neigh_list_cuda;
  }

  MYDBG(printf("# CUDA: Cuda::registerNeighborList() ... end b\n");)
}

void Cuda::uploadAllNeighborLists()
{
  MYDBG(printf("# CUDA: Cuda::uploadAllNeighborList() ... start\n");)
  std::map<NeighList*, CudaNeighList*>::iterator p = neigh_lists.begin();

  while(p != neigh_lists.end()) {
    p->second->nl_upload();

    if(not(p->second->neigh_list->cuda_list->build_cuda))
      for(int i = 0; i < atom->nlocal; i++)
        p->second->sneighlist.maxneighbors = MAX(p->second->neigh_list->numneigh[i], p->second->sneighlist.maxneighbors) ;

    ++p;
  }

  MYDBG(printf("# CUDA: Cuda::uploadAllNeighborList() ... done\n");)
}

void Cuda::downloadAllNeighborLists()
{
  MYDBG(printf("# CUDA: Cuda::downloadAllNeighborList() ... start\n");)
  std::map<NeighList*, CudaNeighList*>::iterator p = neigh_lists.begin();

  while(p != neigh_lists.end()) {
    p->second->nl_download();
    ++p;
  }
}

void Cuda::update_xhold(int &maxhold, double* xhold)
{
  if(this->shared_data.atom.maxhold < atom->nmax) {
    maxhold = atom->nmax;
    delete this->cu_xhold;
    this->cu_xhold     = new cCudaData<double, X_FLOAT, yx> ((double*)xhold, & this->shared_data.atom.xhold         , maxhold, 3);
  }

  this->shared_data.atom.maxhold = maxhold;
  CudaWrapper_CopyData(this->cu_xhold->dev_data(), this->cu_x->dev_data(), 3 * atom->nmax * sizeof(X_FLOAT));
}

void Cuda::setTimingsZero()
{
  shared_data.cuda_timings.test1 = 0;
  shared_data.cuda_timings.test2 = 0;

  //communication
  shared_data.cuda_timings.comm_forward_total = 0;
  shared_data.cuda_timings.comm_forward_mpi_upper = 0;
  shared_data.cuda_timings.comm_forward_mpi_lower = 0;
  shared_data.cuda_timings.comm_forward_kernel_pack = 0;
  shared_data.cuda_timings.comm_forward_kernel_unpack = 0;
  shared_data.cuda_timings.comm_forward_upload = 0;
  shared_data.cuda_timings.comm_forward_download = 0;

  shared_data.cuda_timings.comm_exchange_total = 0;
  shared_data.cuda_timings.comm_exchange_mpi = 0;
  shared_data.cuda_timings.comm_exchange_kernel_pack = 0;
  shared_data.cuda_timings.comm_exchange_kernel_unpack = 0;
  shared_data.cuda_timings.comm_exchange_kernel_fill = 0;
  shared_data.cuda_timings.comm_exchange_cpu_pack = 0;
  shared_data.cuda_timings.comm_exchange_upload = 0;
  shared_data.cuda_timings.comm_exchange_download = 0;

  shared_data.cuda_timings.comm_border_total = 0;
  shared_data.cuda_timings.comm_border_mpi = 0;
  shared_data.cuda_timings.comm_border_kernel_pack = 0;
  shared_data.cuda_timings.comm_border_kernel_unpack = 0;
  shared_data.cuda_timings.comm_border_kernel_buildlist = 0;
  shared_data.cuda_timings.comm_border_kernel_self = 0;
  shared_data.cuda_timings.comm_border_upload = 0;
  shared_data.cuda_timings.comm_border_download = 0;

  //pair forces
  shared_data.cuda_timings.pair_xtype_conversion = 0;
  shared_data.cuda_timings.pair_kernel = 0;
  shared_data.cuda_timings.pair_virial = 0;
  shared_data.cuda_timings.pair_force_collection = 0;

  //neighbor
  shared_data.cuda_timings.neigh_bin = 0;
  shared_data.cuda_timings.neigh_build = 0;
  shared_data.cuda_timings.neigh_special = 0;

  //PPPM
  shared_data.cuda_timings.pppm_particle_map = 0;
  shared_data.cuda_timings.pppm_make_rho = 0;
  shared_data.cuda_timings.pppm_brick2fft = 0;
  shared_data.cuda_timings.pppm_poisson = 0;
  shared_data.cuda_timings.pppm_fillbrick = 0;
  shared_data.cuda_timings.pppm_fieldforce = 0;
  shared_data.cuda_timings.pppm_compute = 0;

  CudaWrapper_CheckUploadTime(true);
  CudaWrapper_CheckDownloadTime(true);
  CudaWrapper_CheckCPUBufUploadTime(true);
  CudaWrapper_CheckCPUBufDownloadTime(true);
}

void Cuda::print_timings()
{
  if(universe->me != 0) return;

  if(not dotiming) return;

  printf("\n # CUDA: Special timings\n\n");
  printf("\n Transfer Times\n");
  printf(" PCIe Upload:  \t %lf s\n", CudaWrapper_CheckUploadTime());
  printf(" PCIe Download:\t %lf s\n", CudaWrapper_CheckDownloadTime());
  printf(" CPU Tempbbuf Upload:   \t %lf \n", CudaWrapper_CheckCPUBufUploadTime());
  printf(" CPU Tempbbuf Download: \t %lf \n", CudaWrapper_CheckCPUBufDownloadTime());

  printf("\n Communication \n");

  printf(" Forward Total           \t %lf \n", shared_data.cuda_timings.comm_forward_total);
  printf(" Forward MPI Upper Bound \t %lf \n", shared_data.cuda_timings.comm_forward_mpi_upper);
  printf(" Forward MPI Lower Bound \t %lf \n", shared_data.cuda_timings.comm_forward_mpi_lower);
  printf(" Forward Kernel Pack     \t %lf \n", shared_data.cuda_timings.comm_forward_kernel_pack);
  printf(" Forward Kernel Unpack   \t %lf \n", shared_data.cuda_timings.comm_forward_kernel_unpack);
  printf(" Forward Kernel Self     \t %lf \n", shared_data.cuda_timings.comm_forward_kernel_self);
  printf(" Forward Upload          \t %lf \n", shared_data.cuda_timings.comm_forward_upload);
  printf(" Forward Download        \t %lf \n", shared_data.cuda_timings.comm_forward_download);
  printf(" Forward Overlap Split Ratio\t %lf \n", shared_data.comm.overlap_split_ratio);
  printf("\n");

  printf(" Exchange Total          \t %lf \n", shared_data.cuda_timings.comm_exchange_total);
  printf(" Exchange MPI            \t %lf \n", shared_data.cuda_timings.comm_exchange_mpi);
  printf(" Exchange Kernel Pack    \t %lf \n", shared_data.cuda_timings.comm_exchange_kernel_pack);
  printf(" Exchange Kernel Unpack  \t %lf \n", shared_data.cuda_timings.comm_exchange_kernel_unpack);
  printf(" Exchange Kernel Fill    \t %lf \n", shared_data.cuda_timings.comm_exchange_kernel_fill);
  printf(" Exchange CPU Pack             \t %lf \n", shared_data.cuda_timings.comm_exchange_cpu_pack);
  printf(" Exchange Upload         \t %lf \n", shared_data.cuda_timings.comm_exchange_upload);
  printf(" Exchange Download       \t %lf \n", shared_data.cuda_timings.comm_exchange_download);
  printf("\n");

  printf(" Border Total            \t %lf \n", shared_data.cuda_timings.comm_border_total);
  printf(" Border MPI              \t %lf \n", shared_data.cuda_timings.comm_border_mpi);
  printf(" Border Kernel Pack      \t %lf \n", shared_data.cuda_timings.comm_border_kernel_pack);
  printf(" Border Kernel Unpack    \t %lf \n", shared_data.cuda_timings.comm_border_kernel_unpack);
  printf(" Border Kernel Self      \t %lf \n", shared_data.cuda_timings.comm_border_kernel_self);
  printf(" Border Kernel BuildList \t %lf \n", shared_data.cuda_timings.comm_border_kernel_buildlist);
  printf(" Border Upload           \t %lf \n", shared_data.cuda_timings.comm_border_upload);
  printf(" Border Download              \t %lf \n", shared_data.cuda_timings.comm_border_download);
  printf("\n");

  //pair forces
  printf(" Pair XType Conversion   \t %lf \n", shared_data.cuda_timings.pair_xtype_conversion);
  printf(" Pair Kernel             \t %lf \n", shared_data.cuda_timings.pair_kernel);
  printf(" Pair Virial             \t %lf \n", shared_data.cuda_timings.pair_virial);
  printf(" Pair Force Collection   \t %lf \n", shared_data.cuda_timings.pair_force_collection);
  printf("\n");

  //neighbor
  printf(" Neighbor Binning        \t %lf \n", shared_data.cuda_timings.neigh_bin);
  printf(" Neighbor Build          \t %lf \n", shared_data.cuda_timings.neigh_build);
  printf(" Neighbor Special        \t %lf \n", shared_data.cuda_timings.neigh_special);
  printf("\n");

  //pppm
  if(force->kspace) {
    printf(" PPPM Total              \t %lf \n", shared_data.cuda_timings.pppm_compute);
    printf(" PPPM Particle Map       \t %lf \n", shared_data.cuda_timings.pppm_particle_map);
    printf(" PPPM Make Rho           \t %lf \n", shared_data.cuda_timings.pppm_make_rho);
    printf(" PPPM Brick2fft          \t %lf \n", shared_data.cuda_timings.pppm_brick2fft);
    printf(" PPPM Poisson            \t %lf \n", shared_data.cuda_timings.pppm_poisson);
    printf(" PPPM Fillbrick          \t %lf \n", shared_data.cuda_timings.pppm_fillbrick);
    printf(" PPPM Fieldforce         \t %lf \n", shared_data.cuda_timings.pppm_fieldforce);
    printf("\n");
  }

  printf(" Debug Test 1            \t %lf \n", shared_data.cuda_timings.test1);
  printf(" Debug Test 2            \t %lf \n", shared_data.cuda_timings.test2);

  printf("\n");
}
