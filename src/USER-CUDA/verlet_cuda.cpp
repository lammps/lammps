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
#include "verlet_cuda.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify_cuda.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "cuda_wrapper_cu.h"
#include "thermo.h"
#include "cuda_pair_cu.h"
#include "cuda.h"
#include <ctime>
#include <cmath>
#ifdef _OPENMP
#include "omp.h"
#endif

using namespace LAMMPS_NS;

#define MAKETIMEING


VerletCuda::VerletCuda(LAMMPS* lmp, int narg, char** arg) : Verlet(lmp, narg, arg)
{
  cuda = lmp->cuda;

  if(cuda == NULL)
    error->all(FLERR, "You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  modify_cuda = (ModifyCuda*) modify;
  int ifix = modify->find_fix("package_omp");

  if(ifix >= 0) external_force_clear = 1;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void VerletCuda::setup()
{
  //debug related variables
  cuda->debugdata[0] = 0;
  cuda->cu_debugdata->upload();
  dotestatom = cuda->dotestatom;
  int testatom = cuda->testatom; //48267;

  if(atom->nlocal == 0)
    error->warning(FLERR, "# CUDA: There are currently no atoms on one of the MPI processes. This is known to cause errors with the USER-CUDA package. Please use the 'processors' keyword to enforce more balanced processor layout.");

  MYDBG(printf("# CUDA VerletCuda::setup start\n");)

  cuda->oncpu = true;
  cuda->begin_setup = true;
  cuda->finished_setup = false;
  cuda->finished_run = false;

  time_pair = 0;
  time_kspace = 0;
  time_comm = 0;
  time_modify = 0;
  time_fulliterate = 0;

  atom->setup();

  cuda_shared_atom*   cu_atom   = & cuda->shared_data.atom;
  cuda_shared_domain* cu_domain = & cuda->shared_data.domain;
  cuda_shared_pair*   cu_pair   = & cuda->shared_data.pair;
  cu_atom->update_nlocal = 1;
  cu_atom->update_nmax = 1;

  if(atom->molecular || (force->kspace && (not cuda->shared_data.pppm.cudable_force))) cuda->shared_data.pair.collect_forces_later = true;

  cuda->setDomainParams();


  if(cuda->shared_data.me == 0)
    printf("# CUDA: VerletCuda::setup: Allocate memory on device for maximum of %i atoms...\n", atom->nmax);

  if(cuda->shared_data.me == 0)
    printf("# CUDA: Using precision: Global: %u X: %u V: %u F: %u PPPM: %u \n", CUDA_PRECISION == 1 ? 4 : 8, sizeof(X_FLOAT), sizeof(V_FLOAT), sizeof(F_FLOAT), sizeof(PPPM_FLOAT));

  cuda->allocate();

  if(comm->me == 0 && screen) fprintf(screen, "Setting up run ...\n");

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists
  modify->setup_pre_exchange();

  if(triclinic) domain->x2lamda(atom->nlocal);

  domain->pbc();
  domain->reset_box();
  comm->setup();

  if(neighbor->style) neighbor->setup_bins();

  comm->exchange();

  if(atom->sortfreq > 0) atom->sort();

  comm->borders();

  if(triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

  cuda->setSystemParams();
  cuda->checkResize();

  if(cuda->shared_data.me == 0)
    printf("# CUDA: VerletCuda::setup: Upload data...\n");

  cuda->uploadAll();
  neighbor->build();
  neighbor->ncalls = 0;

  if(atom->mass)
    cuda->cu_mass->upload();

  if(cuda->cu_map_array)
    cuda->cu_map_array->upload();

  // compute all forces

  ev_set(update->ntimestep);

  if(elist_atom) cuda->shared_data.atom.need_eatom = 1;

  if(vlist_atom) cuda->shared_data.atom.need_vatom = 1;

  if(elist_atom || vlist_atom) cuda->checkResize();

  int test_BpA_vs_TpA = true;
  timespec starttime;
  timespec endtime;
#ifdef NO_PREC_TIMING
  double startsec, endsec;
#endif

  //if(atom->molecular||(force->kspace&&(not cuda->shared_data.pppm.cudable_force))) cuda->shared_data.pair.collect_forces_later = false;
  if(test_BpA_vs_TpA && cuda->shared_data.pair.cudable_force && force->pair && (cuda->shared_data.pair.override_block_per_atom < 0)) {
    int StyleLoops = 10;

    if(cuda->shared_data.me == 0)
      printf("Test TpA\n");

    cuda->shared_data.pair.use_block_per_atom = 0;
    neighbor->build();
    Cuda_Pair_GenerateXType(&cuda->shared_data);

    if(cuda->cu_v_radius)
      Cuda_Pair_GenerateVRadius(&cuda->shared_data);

    if(cuda->cu_omega_rmass)
      Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

    force->pair->compute(eflag, vflag);
    CudaWrapper_Sync();
#ifdef NO_PREC_TIMING
    startsec = 1.0 * clock() / CLOCKS_PER_SEC;
#endif
    clock_gettime(CLOCK_REALTIME, &starttime);

    for(int i = 0; i < StyleLoops; i++) {
      Cuda_Pair_GenerateXType(&cuda->shared_data);

      if(cuda->cu_v_radius)
        Cuda_Pair_GenerateVRadius(&cuda->shared_data);

      if(cuda->cu_omega_rmass)
        Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

      force->pair->compute(eflag, vflag);
      CudaWrapper_Sync();
    }

    clock_gettime(CLOCK_REALTIME, &endtime);

    double TpAtime = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
#ifdef NO_PREC_TIMING
    endsec = 1.0 * clock() / CLOCKS_PER_SEC;
    TpAtime = endsec - startsec;
#endif

    if(cuda->shared_data.me == 0)
      printf("Test BpA\n");

    cuda->shared_data.pair.use_block_per_atom = 1;
    neighbor->build();
    Cuda_Pair_GenerateXType(&cuda->shared_data);

    if(cuda->cu_v_radius)
      Cuda_Pair_GenerateVRadius(&cuda->shared_data);

    if(cuda->cu_omega_rmass)
      Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

    force->pair->compute(eflag, vflag);
    CudaWrapper_Sync();

    clock_gettime(CLOCK_REALTIME, &starttime);
#ifdef NO_PREC_TIMING
    startsec = 1.0 * clock() / CLOCKS_PER_SEC;
#endif

    for(int i = 0; i < StyleLoops; i++) {
      Cuda_Pair_GenerateXType(&cuda->shared_data);

      if(cuda->cu_v_radius)
        Cuda_Pair_GenerateVRadius(&cuda->shared_data);

      if(cuda->cu_omega_rmass)
        Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

      force->pair->compute(eflag, vflag);
      CudaWrapper_Sync();
    }

    clock_gettime(CLOCK_REALTIME, &endtime);
    double BpAtime = endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
#ifdef NO_PREC_TIMING
    endsec = 1.0 * clock() / CLOCKS_PER_SEC;
    BpAtime = endsec - startsec;
#endif

    if(cuda->shared_data.me == 0)
      printf("\n# CUDA: Timing of parallelisation layout with %i loops:\n", StyleLoops);

    if(cuda->shared_data.me == 0)
      printf("# CUDA: BpA TpA\n %lf %lf\n", BpAtime, TpAtime);

    if(BpAtime > TpAtime) cuda->shared_data.pair.use_block_per_atom = 0;
  } else
    cuda->shared_data.pair.use_block_per_atom = cuda->shared_data.pair.override_block_per_atom;

  //cuda->shared_data.pair.use_block_per_atom = 0;
  if(atom->molecular || (force->kspace && (not cuda->shared_data.pppm.cudable_force))) cuda->shared_data.pair.collect_forces_later = true;

  neighbor->build();
  neighbor->ncalls = 0;

  force_clear();

  modify->setup_pre_force(vflag);

  cuda->cu_f->download();

  if(cuda->cu_torque)
    cuda->cu_torque->download();

  //printf("# Verlet::setup: g f[0] = (%f, %f, %f)\n", atom->f[0][0], atom->f[0][1], atom->f[0][2]);

  MYDBG(printf("# CUDA: VerletCuda::setup: initial force compute\n");)

  //test_atom(testatom,"pre pair force");

  if(cuda->shared_data.pair.cudable_force) {
    cuda->uploadAll();
    Cuda_Pair_GenerateXType(&cuda->shared_data);

    if(cuda->cu_v_radius)
      Cuda_Pair_GenerateVRadius(&cuda->shared_data);

    if(cuda->cu_omega_rmass)
      Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);
  }

  if(force->pair) force->pair->compute(eflag, vflag);

  if(cuda->shared_data.pair.cudable_force) {
    if(cuda->shared_data.pair.collect_forces_later) {
      if(eflag) cuda->cu_eng_vdwl->upload();

      if(eflag) cuda->cu_eng_coul->upload();

      if(vflag) cuda->cu_virial->upload();

      Cuda_Pair_CollectForces(&cuda->shared_data, eflag, vflag);

      if(eflag) cuda->cu_eng_vdwl->download();

      if(eflag) cuda->cu_eng_coul->download();

      if(vflag) cuda->cu_virial->download();
    }

    cuda->downloadAll();
  }

  test_atom(testatom, "post pair force");

  MYDBG(printf("# CUDA: VerletCuda::setup: initial force compute done\n");)
  //printf("# Verlet::setup: h f[0] = (%f, %f, %f)\n", atom->f[0][0], atom->f[0][1], atom->f[0][2]);

  if(atom->molecular) {
    if(force->bond) force->bond->compute(eflag, vflag);

    if(force->angle) force->angle->compute(eflag, vflag);

    if(force->dihedral) force->dihedral->compute(eflag, vflag);

    if(force->improper) force->improper->compute(eflag, vflag);
  }


  if(cuda->shared_data.pppm.cudable_force) {
    cuda->cu_tag ->upload();
    cuda->cu_type->upload();
    cuda->cu_x   ->upload();
    cuda->cu_v   ->upload();
    cuda->cu_f   ->upload();

    if(cu_atom->q_flag) cuda->cu_q->upload();
  }

  if(force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag, vflag);
  }

  if(cuda->shared_data.pppm.cudable_force) {
    cuda->cu_f   ->download();
  }

  test_atom(testatom, "post kspace");

  cuda->uploadAll();

  if(force->newton) comm->reverse_comm();

  cuda->downloadAll();

  test_atom(testatom, "post reverse comm");

  if(cuda->shared_data.me == 0)
    printf("# CUDA: Total Device Memory useage post setup: %lf MB\n", 1.0 * CudaWrapper_CheckMemUseage() / 1024 / 1024);

  MYDBG(printf("# CUDA: VerletCuda::setup: call modify setup\n");)
  modify->setup(vflag);

  MYDBG(printf("# CUDA: VerletCuda::setup: call modify setup done\n");)
  output->setup(1);

  test_atom(testatom, "post setup");

  MYDBG(printf("# CUDA: VerletCuda::setup: done\n");)
  cuda->finished_setup = true;
  cuda->oncpu = false;
}


//this routine is in a messy state
void VerletCuda::setup_minimal(int flag)
{

  printf("SetupMinimal\n");
  dotestatom = 0;
  int testatom = 104;
  cuda->oncpu = true;
  cuda->begin_setup = true;
  cuda->finished_run = false;
  MYDBG(printf("# CUDA VerletCuda::setup start\n");)
  time_pair = 0;
  time_kspace = 0;
  time_comm = 0;
  time_modify = 0;
  time_fulliterate = 0;

  //cuda->allocate();

  cuda_shared_atom*   cu_atom   = & cuda->shared_data.atom;
  cuda_shared_domain* cu_domain = & cuda->shared_data.domain;
  cuda_shared_pair*   cu_pair   = & cuda->shared_data.pair;
  cu_atom->update_nlocal = 1;
  cu_atom->update_nmax = 1;

  if(atom->molecular) cuda->shared_data.pair.collect_forces_later = true;

  cuda->setDomainParams();



  if(cuda->shared_data.me == 0)
    printf("# CUDA: VerletCuda::setup: Allocate memory on device for maximum of %i atoms...\n", atom->nmax);

  cuda->allocate();




  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if(flag) {
    if(triclinic) domain->x2lamda(atom->nlocal);

    domain->pbc();
    domain->reset_box();
    comm->setup();

    if(neighbor->style) neighbor->setup_bins();

    comm->exchange();
    comm->borders();

    if(triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

    cuda->setSystemParams();
    cuda->checkResize();
    neighbor->build();
    neighbor->ncalls = 0;
  }

  if(cuda->shared_data.me == 0)
    printf("# CUDA: VerletCuda::setup: Upload data...\n");

  cuda->uploadAll();
  cuda->uploadAllNeighborLists();

  if(atom->mass)
    cuda->cu_mass->upload();

  if(cuda->cu_map_array)
    cuda->cu_map_array->upload();

  // compute all forces

  ev_set(update->ntimestep);

  if(elist_atom) cuda->shared_data.atom.need_eatom = 1;

  if(vlist_atom) cuda->shared_data.atom.need_vatom = 1;

  if(elist_atom || vlist_atom) cuda->checkResize();

  force_clear();
  cuda->cu_f->download();

  //printf("# Verlet::setup: g f[0] = (%f, %f, %f)\n", atom->f[0][0], atom->f[0][1], atom->f[0][2]);

  cuda->cu_mass->upload();
  MYDBG(printf("# CUDA: VerletCuda::setup: initial force compute\n");)

  test_atom(testatom, "pre pair force");

  if(cuda->shared_data.pair.cudable_force) {
    cuda->uploadAll();
    Cuda_Pair_GenerateXType(&cuda->shared_data);

    if(cuda->cu_v_radius)
      Cuda_Pair_GenerateVRadius(&cuda->shared_data);

    if(cuda->cu_omega_rmass)
      Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);
  }

  if(force->pair) force->pair->compute(eflag, vflag);

  if(cuda->shared_data.pair.cudable_force) {
    if(cuda->shared_data.pair.collect_forces_later) {
      if(eflag) cuda->cu_eng_vdwl->upload();

      if(eflag) cuda->cu_eng_coul->upload();

      if(vflag) cuda->cu_virial->upload();

      Cuda_Pair_CollectForces(&cuda->shared_data, eflag, vflag);

      if(eflag) cuda->cu_eng_vdwl->download();

      if(eflag) cuda->cu_eng_coul->download();

      if(vflag) cuda->cu_virial->download();
    }

    cuda->downloadAll();
  }

  test_atom(testatom, "post pair force");

  MYDBG(printf("# CUDA: VerletCuda::setup: initial force compute done\n");)
  //printf("# Verlet::setup: h f[0] = (%f, %f, %f)\n", atom->f[0][0], atom->f[0][1], atom->f[0][2]);

  if(atom->molecular) {
    if(force->bond) force->bond->compute(eflag, vflag);

    if(force->angle) force->angle->compute(eflag, vflag);

    if(force->dihedral) force->dihedral->compute(eflag, vflag);

    if(force->improper) force->improper->compute(eflag, vflag);
  }


  if(cuda->shared_data.pppm.cudable_force) {
    cuda->cu_tag ->upload();
    cuda->cu_type->upload();
    cuda->cu_x   ->upload();
    cuda->cu_v   ->upload();
    cuda->cu_f   ->upload();

    if(cu_atom->q_flag) cuda->cu_q->upload();
  }

  if(force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag, vflag);
  }

  if(cuda->shared_data.pppm.cudable_force) {
    cuda->cu_f   ->download();
  }

  test_atom(testatom, "post kspace");

  cuda->uploadAll();

  if(force->newton) comm->reverse_comm();

  cuda->downloadAll();

  test_atom(testatom, "post reverse comm");

  if(cuda->shared_data.me == 0)
    printf("# CUDA: Total Device Memory useage post setup: %lf MB\n", 1.0 * CudaWrapper_CheckMemUseage() / 1024 / 1024);

  MYDBG(printf("# CUDA: VerletCuda::setup: call modify setup\n");)
  modify->setup(vflag);

  MYDBG(printf("# CUDA: VerletCuda::setup: done\n");)
  cuda->finished_setup = true;
  cuda->oncpu = false;
}

//#define TESTATOM
/* ----------------------------------------------------------------------
   iterate for n steps
------------------------------------------------------------------------- */

void VerletCuda::run(int n)
{
  dotestatom = cuda->dotestatom;
  int testatom = cuda->testatom; //48267;


  timespec starttime;
  timespec endtime;
  timespec starttotal;
  timespec endtotal;

  cuda->setTimingsZero();

  static double testtime = 0.0;
  //                                clock_gettime(CLOCK_REALTIME,&starttime);
  //                                  clock_gettime(CLOCK_REALTIME,&endtime);
  //                                testtime+=endtime.tv_sec-starttime.tv_sec+1.0*(endtime.tv_nsec-starttime.tv_nsec)/1000000000;
  //                                 printf("Time: %lf\n",testtime);*/


  cuda_shared_domain* cu_domain = & cuda->shared_data.domain;

  int nflag, ntimestep, sortflag;

  int n_initial_integrate = modify_cuda->n_initial_integrate;
  int n_post_integrate = modify_cuda->n_post_integrate;
  int n_final_integrate = modify_cuda->n_final_integrate;
  int n_pre_exchange = modify_cuda->n_pre_exchange;
  int n_pre_neighbor = modify_cuda->n_pre_neighbor;
  int n_pre_force = modify_cuda->n_pre_force;
  int n_post_force = modify_cuda->n_post_force;
  int n_end_of_step = modify_cuda->n_end_of_step;
  MYDBG(printf("# CUDA: Fixes: i_int: %i p_int: %i f_int: %i pr_exc: %i pr_neigh: %i pr_f: %i p_f: %i eos: %i\n",
               n_initial_integrate, n_post_integrate, n_final_integrate, n_pre_exchange, n_pre_neighbor, n_pre_force, n_post_force, n_end_of_step);)

  if(atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;


  if(cuda->shared_data.me == 0) {
    if((not cuda->shared_data.pair.cudable_force) && (force->pair))
      error->warning(FLERR, "# CUDA: You asked for a Verlet integration using Cuda, "
                     "but selected a pair force which has not yet been ported to Cuda");

    if((not cuda->shared_data.pppm.cudable_force) && (force->kspace))
      error->warning(FLERR, "# CUDA: You asked for a Verlet integration using Cuda, "
                     "but selected a kspace force which has not yet been ported to Cuda");

    if(modify_cuda->n_post_integrate_host + modify_cuda->n_pre_exchange_host + modify_cuda->n_pre_neighbor_host + modify_cuda->n_pre_force_host + modify_cuda->n_post_force_host + modify_cuda->n_end_of_step_host + modify_cuda->n_initial_integrate_host + modify_cuda->n_final_integrate_host)
      error->warning(FLERR, "# CUDA: You asked for a Verlet integration using Cuda, "
                     "but several fixes have not yet been ported to Cuda.\n"
                     "This can cause a severe speed penalty due to frequent data synchronization between host and GPU.");

    if(atom->firstgroupname)
      error->warning(FLERR, "Warning: firstgroupname is used, this will cause additional data transfers.");
  }

  cuda->uploadAll();

  if(cuda->neighbor_decide_by_integrator && cuda->cu_xhold) {
    const int n = cuda->shared_data.atom.maxhold;
    CudaWrapper_CopyData(cuda->cu_xhold->dev_data(), cuda->cu_x->dev_data(), n * sizeof(X_FLOAT));
    CudaWrapper_CopyData((void*) & ((X_FLOAT*)cuda->cu_xhold->dev_data())[n], (void*) & ((X_FLOAT*)cuda->cu_x->dev_data())[atom->nmax], n * sizeof(X_FLOAT));
    CudaWrapper_CopyData((void*) & ((X_FLOAT*)cuda->cu_xhold->dev_data())[2 * n], (void*) & ((X_FLOAT*)cuda->cu_x->dev_data())[2 * atom->nmax], n * sizeof(X_FLOAT));
  }

  cuda->shared_data.atom.reneigh_flag = 0;
  cuda->shared_data.atom.update_nlocal = 1;
  cuda->shared_data.atom.update_nmax = 1;
  cuda->shared_data.atom.update_neigh = 1;
  cuda->shared_data.domain.update = 1;
  cuda->shared_data.buffer_new = 1;
  cuda->uploadtime = 0;
  cuda->downloadtime = 0;
  int firstreneigh = 1;

  for(int i = 0; i < n; i++) {
    if(atom->nlocal == 0)
      error->warning(FLERR, "# CUDA: There are currently no atoms on one of the MPI processes. This is currently prone to encountering errors with USER-CUDA package. Please use the 'processors' keyword to use a more balanced processor layout.");

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    test_atom(testatom, "Pre initial");

    MYDBG(printf("# CUDA VerletCuda::iterate: before initial_integrate\n");)

    modify->initial_integrate(vflag);

    MYDBG(printf("# CUDA VerletCuda::iterate: after initial_integrate\n");)

    if(n_post_integrate) modify->post_integrate();



    // regular communication vs neighbor list rebuild

    test_atom(testatom, "Pre Exchange");

    MYDBG(printf("# CUDA VerletCuda::iterate: before neighbor decide\n");)
    nflag = neighbor->decide();

    if(nflag == 0) {
      MYDBG(printf("# CUDA VerletCuda::iterate: communicate\n");)
      timer->stamp();

      if((not(eflag || vflag)) && (cuda->shared_data.overlap_comm)) {
        //overlap forward communication of ghost atom positions with inner force calculation (interactions between local atoms)
        //build communication buffers
        //      printf("Pre forward_comm(1)\n");
        clock_gettime(CLOCK_REALTIME, &starttotal);
        cuda->shared_data.atom.reneigh_flag = 0;
        clock_gettime(CLOCK_REALTIME, &starttime);
        timer->stamp();
        comm->forward_comm(1);
        timer->stamp(TIME_COMM);
        clock_gettime(CLOCK_REALTIME, &endtime);
        cuda->shared_data.cuda_timings.comm_forward_total +=
          endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

        //prepare force calculation
        //     printf("Pre force_clear\n");
        force_clear();
        //     printf("Pre Generate XType\n");
        Cuda_Pair_GenerateXType(&cuda->shared_data);

        if(cuda->cu_v_radius)
          Cuda_Pair_GenerateVRadius(&cuda->shared_data);

        if(cuda->cu_omega_rmass)
          Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

        //start force calculation asynchronus
        cuda->shared_data.comm.comm_phase = 1;
        force->pair->compute(eflag, vflag);
        timer->stamp(TIME_PAIR);
        //CudaWrapper_Sync();

        //download comm buffers from GPU, perform MPI communication and upload buffers again
        clock_gettime(CLOCK_REALTIME, &starttime);
        comm->forward_comm(2);
        clock_gettime(CLOCK_REALTIME, &endtime);
        cuda->shared_data.cuda_timings.comm_forward_total +=
          endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
        timer->stamp(TIME_COMM);

        //wait for force calculation
        CudaWrapper_Sync();
        timer->stamp(TIME_PAIR);

        //unpack communication buffers
        clock_gettime(CLOCK_REALTIME, &starttime);
        comm->forward_comm(3);
        clock_gettime(CLOCK_REALTIME, &endtime);
        cuda->shared_data.cuda_timings.comm_forward_total +=
          endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

        timer->stamp(TIME_COMM);
        MYDBG(printf("# CUDA VerletCuda::iterate: communicate done\n");)
        cuda->shared_data.cuda_timings.test1 +=
          endtotal.tv_sec - starttotal.tv_sec + 1.0 * (endtotal.tv_nsec - starttotal.tv_nsec) / 1000000000;
      } else {
        //perform standard forward communication
        clock_gettime(CLOCK_REALTIME, &starttime);
        comm->forward_comm();
        clock_gettime(CLOCK_REALTIME, &endtime);
        cuda->shared_data.cuda_timings.comm_forward_total +=
          endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
        timer->stamp(TIME_COMM);
        MYDBG(printf("# CUDA VerletCuda::iterate: communicate done\n");)
      }
    } else {
      int nlocalold = cuda->shared_data.atom.nlocal;

      if(firstreneigh) {
        cuda->shared_data.atom.update_nlocal = 1;
        cuda->shared_data.atom.update_nmax = 1;
        firstreneigh = 0;
      }

      cuda->shared_data.buffer_new = 1;
      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor\n");)
      cuda->setDomainParams();

      if(n_pre_exchange) modify->pre_exchange();

      if(atom->nlocal != cuda->shared_data.atom.nlocal) { //did someone add atoms during pre_exchange?
        cuda->checkResize();
        cuda->uploadAll();
      }

      //check domain changes
      if(domain->triclinic) domain->x2lamda(atom->nlocal);

      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor pbc\n");)
      domain->pbc();

      if(domain->box_change) {
        domain->reset_box();
        comm->setup();

        if(neighbor->style) neighbor->setup_bins();

      }

      timer->stamp();
      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor exchange\n");)

      //perform exchange of local atoms
      clock_gettime(CLOCK_REALTIME, &starttime);
      comm->exchange();
      clock_gettime(CLOCK_REALTIME, &endtime);

      //special and nspecial fields of the atom data are not currently transfered via the GPU buffer might be changed in the future
      if(comm->nprocs > 1) {
        clock_gettime(CLOCK_REALTIME, &starttime);

        if(atom->special)
          cuda->cu_special->upload();

        if(atom->nspecial)
          cuda->cu_nspecial->upload();

        clock_gettime(CLOCK_REALTIME, &endtime);
        cuda->shared_data.cuda_timings.test1 +=
          endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
      }

      cuda->shared_data.cuda_timings.comm_exchange_total +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

      if(nlocalold != cuda->shared_data.atom.nlocal) cuda->shared_data.atom.update_nlocal = 2;

      //sort atoms
      if(sortflag && ntimestep >= atom->nextsort) atom->sort();

      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor borders\n");)

      //generate ghost atom lists, and transfer ghost atom data
      clock_gettime(CLOCK_REALTIME, &starttime);
      comm->borders();
      clock_gettime(CLOCK_REALTIME, &endtime);
      cuda->shared_data.cuda_timings.comm_border_total +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

      clock_gettime(CLOCK_REALTIME, &starttime);
      //atom index maps are generated on CPU, and need to be transfered to GPU if they are used
      if(cuda->cu_map_array)
        cuda->cu_map_array->upload();


      if(domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

      if(n_pre_neighbor) modify->pre_neighbor();

      cuda->shared_data.buffer_new = 2;

      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor build\n");)
      timer->stamp(TIME_COMM);
      clock_gettime(CLOCK_REALTIME, &endtime);
      cuda->shared_data.cuda_timings.test2 +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;

      //rebuild neighbor list
      test_atom(testatom, "Pre Neighbor");
      neighbor->build(0);
      timer->stamp(TIME_NEIGHBOR);
      MYDBG(printf("# CUDA VerletCuda::iterate: neighbor done\n");)
      //if bonded interactions are used (in this case collect_forces_later is true), transfer data which only changes upon exchange/border routines from GPU to CPU
      if(cuda->shared_data.pair.collect_forces_later) {
        if(cuda->cu_molecule) cuda->cu_molecule->downloadAsync(2);

        cuda->cu_tag->downloadAsync(2);
        cuda->cu_type->downloadAsync(2);
        cuda->cu_mask->downloadAsync(2);

        if(cuda->cu_q) cuda->cu_q->downloadAsync(2);
      }
      cuda->shared_data.comm.comm_phase = 3;
    }

    test_atom(testatom, "Post Exchange");

    // force computations

    //only do force_clear if it has not been done during overlap of communication with local interactions
    if(not((not(eflag || vflag)) && (cuda->shared_data.overlap_comm) && (cuda->shared_data.comm.comm_phase < 3)))
      force_clear();

    if(n_pre_force) modify->pre_force(vflag);

    timer->stamp();

    //if overlap of bonded interactions with nonbonded interactions takes place, download forces and positions
    /*            if(cuda->shared_data.pair.collect_forces_later)
               {
                 cuda->cu_x->downloadAsync(2);
                 cuda->cu_f->downloadAsync(2);
               }*/

    if(force->pair) {
      if((not(eflag || vflag)) && (cuda->shared_data.overlap_comm) && (cuda->shared_data.comm.comm_phase < 3) && cuda->shared_data.pair.cudable_force) {
        //second part of force calculations in case of overlaping it with commuincation. Only interactions between local and ghost atoms are done now
        //regenerate data layout for force computations, its actually only needed for the ghost atoms
        cuda->shared_data.comm.comm_phase = 2;

        timespec atime1, atime2;
        clock_gettime(CLOCK_REALTIME, &atime1);

        Cuda_Pair_GenerateXType(&cuda->shared_data);

        if(cuda->cu_v_radius)
          Cuda_Pair_GenerateVRadius(&cuda->shared_data);

        if(cuda->cu_omega_rmass)
          Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

        clock_gettime(CLOCK_REALTIME, &atime2);
        cuda->shared_data.cuda_timings.pair_xtype_conversion +=
          atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;
        force->pair->compute(eflag, vflag);

      } else {
        //calculate complete pair interactions
        if(not cuda->shared_data.pair.cudable_force) cuda->downloadAll();
        else {
          //regenerate data layout for force computations, its actually only needed for the ghost atoms
          timespec atime1, atime2;
          clock_gettime(CLOCK_REALTIME, &atime1);

          Cuda_Pair_GenerateXType(&cuda->shared_data);

          if(cuda->cu_v_radius)
            Cuda_Pair_GenerateVRadius(&cuda->shared_data);

          if(cuda->cu_omega_rmass)
            Cuda_Pair_GenerateOmegaRmass(&cuda->shared_data);

          clock_gettime(CLOCK_REALTIME, &atime2);
          cuda->shared_data.cuda_timings.pair_xtype_conversion +=
            atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;
        }

        cuda->shared_data.comm.comm_phase = 0;
        force->pair->compute(eflag, vflag);
      }

      if(not cuda->shared_data.pair.cudable_force) cuda->uploadAll();

      //wait for force calculation in case of not using overlap with bonded interactions
      if(not cuda->shared_data.pair.collect_forces_later)
        CudaWrapper_Sync();

      timer->stamp(TIME_PAIR);
    }

    //calculate bonded interactions
    if(atom->molecular) {
      cuda->cu_x->downloadAsync(2);

      if(n_pre_force == 0) Verlet::force_clear();
      else  cuda->cu_f->downloadAsync(2);

      timer->stamp(TIME_PAIR);

      if(neighbor->lastcall == update->ntimestep) {
        neighbor->build_topology();
        timer->stamp(TIME_NEIGHBOR);
      }

      test_atom(testatom, "pre bond force");

      if(force->bond) force->bond->compute(eflag, vflag);

      if(force->angle) force->angle->compute(eflag, vflag);

      if(force->dihedral) force->dihedral->compute(eflag, vflag);

      if(force->improper) force->improper->compute(eflag, vflag);

      timer->stamp(TIME_BOND);
    }

    //collect forces in case pair force and bonded interactions were overlapped, and either no KSPACE or a GPU KSPACE style is used
    if(cuda->shared_data.pair.collect_forces_later && cuda->shared_data.pair.cudable_force && (not(force->kspace && (not cuda->shared_data.pppm.cudable_force)))) {
      clock_gettime(CLOCK_REALTIME, &starttime);
      cuda->cu_f->uploadAsync(2);

      test_atom(testatom, "post molecular force");


      if(eflag) cuda->cu_eng_vdwl->upload();

      if(eflag) cuda->cu_eng_coul->upload();

      if(vflag) cuda->cu_virial->upload();

      Cuda_Pair_CollectForces(&cuda->shared_data, eflag, vflag);

      if(eflag) cuda->cu_eng_vdwl->download();

      if(eflag) cuda->cu_eng_coul->download();

      if(vflag) cuda->cu_virial->download();

      timer->stamp(TIME_PAIR);

      clock_gettime(CLOCK_REALTIME, &endtime);
      cuda->shared_data.cuda_timings.pair_force_collection +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
    }

    //compute kspace force
    if(force->kspace) {
      if((not cuda->shared_data.pppm.cudable_force) && (not cuda->shared_data.pair.collect_forces_later))
        cuda->downloadAll();

      if((not cuda->shared_data.pppm.cudable_force) && (cuda->shared_data.pair.collect_forces_later) && (not atom->molecular)) {
        cuda->cu_x->downloadAsync(2);

        if(n_pre_force == 0) Verlet::force_clear();
        else  cuda->cu_f->downloadAsync(2);

        timer->stamp(TIME_PAIR);
      }

      force->kspace->compute(eflag, vflag);

      if((not cuda->shared_data.pppm.cudable_force) && (not cuda->shared_data.pair.collect_forces_later))
        cuda->uploadAll();

      timer->stamp(TIME_KSPACE);
    }

    //collect forces in case pair forces and kspace was overlaped
    if(cuda->shared_data.pair.collect_forces_later && cuda->shared_data.pair.cudable_force && ((force->kspace && (not cuda->shared_data.pppm.cudable_force)))) {
      cuda->cu_f->uploadAsync(2);

      clock_gettime(CLOCK_REALTIME, &starttime);

      if(eflag) cuda->cu_eng_vdwl->upload();

      if(eflag) cuda->cu_eng_coul->upload();

      if(vflag) cuda->cu_virial->upload();

      Cuda_Pair_CollectForces(&cuda->shared_data, eflag, vflag);

      if(eflag) cuda->cu_eng_vdwl->download();

      if(eflag) cuda->cu_eng_coul->download();

      if(vflag) cuda->cu_virial->download();

      timer->stamp(TIME_PAIR);

      clock_gettime(CLOCK_REALTIME, &endtime);
      cuda->shared_data.cuda_timings.pair_force_collection +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
    }

    //send forces on ghost atoms back to other GPU: THIS SHOULD NEVER HAPPEN
    if(force->newton) {
      comm->reverse_comm();
      timer->stamp(TIME_COMM);
    }

    test_atom(testatom, "post force");
    // force modifications, final time integration, diagnostics

    if(n_post_force) modify->post_force(vflag);

    test_atom(testatom, "pre final");

    modify->final_integrate();

    test_atom(testatom, "post final");

    if(n_end_of_step) modify->end_of_step();

    // all output

    test_atom(testatom, "pre output");

    if(ntimestep == output->next) {
      if(not output->thermo->cudable)
        cuda->downloadAll();

      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }


    test_atom(testatom, "post output");

    if(cuda->shared_data.atom.update_nlocal > 0)
      cuda->shared_data.atom.update_nlocal--;

    if(cuda->shared_data.atom.update_nmax > 0)
      cuda->shared_data.atom.update_nmax--;

    if(cuda->shared_data.atom.update_neigh > 0)
      cuda->shared_data.atom.update_neigh--;

    if(cuda->shared_data.domain.update > 0)
      cuda->shared_data.domain.update--;

    if(cuda->shared_data.buffer_new > 0)
      cuda->shared_data.buffer_new--;

    cuda->shared_data.atom.reneigh_flag = 0;
  }


  cuda->downloadAll();
  cuda->downloadAllNeighborLists();
  cuda->shared_data.atom.update_nlocal = 1;
  cuda->shared_data.atom.update_nmax = 1;
  cuda->shared_data.atom.update_neigh = 1;
  cuda->shared_data.buffer_new = 1;
  cuda->shared_data.domain.update = 1;
  cuda->oncpu = true;
  cuda->finished_run = true;
}


/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void VerletCuda::force_clear()
{
  cuda->cu_f->memset_device(0);

  if(cuda->cu_torque) cuda->cu_torque->memset_device(0);

  return;

  //The rest should not be necessary
  int i;

  for(i = 0; i < atom->nlocal; i++) {
    atom->f[i][0] = 0.0;
    atom->f[i][1] = 0.0;
    atom->f[i][2] = 0.0;
  }

  // clear force on all particles
  // if either newton flag is set, also include ghosts

  if(neighbor->includegroup == 0) {
    int nall;

    if(force->newton) nall = atom->nlocal + atom->nghost;
    else nall = atom->nlocal;

    if(torqueflag) {
      double** torque = atom->torque;

      for(i = 0; i < nall; i++) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
        torque[i][2] = 0.0;
      }
    }

    // neighbor includegroup flag is set
    // clear force only on initial nfirst particles
    // if either newton flag is set, also include ghosts

  } else {
    int nall = atom->nfirst;


    if(torqueflag) {
      double** torque = atom->torque;

      for(i = 0; i < nall; i++) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
        torque[i][2] = 0.0;
      }
    }

    if(force->newton) {
      nall = atom->nlocal + atom->nghost;

      if(torqueflag) {
        double** torque = atom->torque;

        for(i = atom->nlocal; i < nall; i++) {
          torque[i][0] = 0.0;
          torque[i][1] = 0.0;
          torque[i][2] = 0.0;
        }
      }
    }
  }
}

void VerletCuda::test_atom(int aatom, char* string)  //printing properties of one atom for test purposes
{
  if(not dotestatom) return;

  bool check = false;

  if(cuda->finished_setup) cuda->downloadAll();

  for(int i = 0; i < atom->nlocal + atom->nghost; i++) {
    if((atom->tag[i] == aatom) && (i < atom->nlocal)) {

      printf("%i # CUDA %s: %i %i %e %e %e %i ", comm->me, string, update->ntimestep, atom->tag[i], atom->x[i][0], atom->v[i][0], atom->f[i][0], i);

      if(atom->molecular && (i < atom->nlocal)) {
        printf(" // %i %i %i ", atom->num_bond[i], atom->num_angle[i], atom->num_dihedral[i]);

        for(int k = 0; k < atom->num_bond[i]; k++)
          printf("// %i %i ", atom->bond_type[i][k], atom->bond_atom[i][k]);
      }

      printf("\n");
    }

    if(i < atom->nlocal) {
      if((atom->v[i][0] < -100 || atom->v[i][0] > 100) ||
          (atom->v[i][1] < -100 || atom->v[i][1] > 100) ||
          (atom->v[i][2] < -100 || atom->v[i][2] > 100) ||
          (atom->v[i][0] != atom->v[i][0]) ||
          (atom->v[i][1] != atom->v[i][1]) ||
          (atom->v[i][2] != atom->v[i][2])) {
        printf("%i # CUDA %s velocity: %i %e %e %e %i\n", comm->me, string, atom->tag[i], atom->x[i][0], atom->v[i][0], atom->f[i][0], i);
        check = true;
      }

      if((atom->f[i][0] < -10000 || atom->f[i][0] > 10000) ||
          (atom->f[i][1] < -10000 || atom->f[i][1] > 10000) ||
          (atom->f[i][2] < -10000 || atom->f[i][2] > 10000) ||
          (atom->f[i][0] != atom->f[i][0]) ||
          (atom->f[i][1] != atom->f[i][1]) ||
          (atom->f[i][2] != atom->f[i][2])) {
        printf("%i # CUDA %s force: %i %e %e %e %i\n", comm->me, string, atom->tag[i], atom->x[i][0], atom->v[i][0], atom->f[i][0], i);
        check = true;
      }

      if(atom->tag[i] <= 0)
        printf("%i # CUDA %s tag: %i %e %e %e %i\n", comm->me, string, atom->tag[i], atom->x[i][0], atom->v[i][0], atom->f[i][0], i);
    }
  }

  if(check) exit(0);
}
