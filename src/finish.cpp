// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "finish.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"             // IWYU pragma: keep
#include "min.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"           // IWYU pragma: keep
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "timer.h"              // IWYU pragma: keep
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

#ifdef LMP_OPENMP
#include "modify.h"
#include "fix_omp.h"
#include "thr_data.h"
#endif

using namespace LAMMPS_NS;

// local function prototypes, code at end of file

static void mpi_timings(const char *label, Timer *t, enum Timer::ttype tt,
                        MPI_Comm world, const int nprocs, const int nthreads,
                        const int me, double time_loop, FILE *scr, FILE *log);

#ifdef LMP_OPENMP
static void omp_times(FixOMP *fix, const char *label, enum Timer::ttype which,
                      const int nthreads,FILE *scr, FILE *log);
#endif

/* ---------------------------------------------------------------------- */

Finish::Finish(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Finish::end(int flag)
{
  int i,nneigh,nneighfull;
  int histo[10];
  int minflag,prdflag,tadflag,hyperflag;
  int timeflag,fftflag,histoflag,neighflag;
  double time,tmp,ave,max,min;
  double time_loop,time_other,cpu_loop;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  const int nthreads = comm->nthreads;

  // recompute natoms in case atoms have been lost

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // choose flavors of statistical output
  // flag determines caller
  // flag = 0 = just loop summary
  // flag = 1 = dynamics or minimization
  // flag = 2 = PRD
  // flag = 3 = TAD
  // flag = 4 = HYPER
  // turn off neighflag for Kspace partition of verlet/split integrator

  minflag = prdflag = tadflag = hyperflag = 0;
  timeflag = fftflag = histoflag = neighflag = 0;
  time_loop = cpu_loop = time_other = 0.0;

  if (flag == 1) {
    if (update->whichflag == 2) minflag = 1;
    timeflag = histoflag = 1;
    neighflag = 1;
    if (update->whichflag == 1 &&
        strncmp(update->integrate_style,"verlet/split",12) == 0 &&
        universe->iworld == 1) neighflag = 0;
    if (force->kspace && force->kspace_match("^pppm",0)
        && force->kspace->fftbench) fftflag = 1;
  }
  if (flag == 2) prdflag = timeflag = histoflag = neighflag = 1;
  if (flag == 3) tadflag = histoflag = neighflag = 1;
  if (flag == 4) hyperflag = timeflag = histoflag = neighflag = 1;

  // loop stats

  if (timer->has_loop()) {

    // overall loop time

    time_loop = timer->get_wall(Timer::TOTAL);
    cpu_loop = timer->get_cpu(Timer::TOTAL);
    MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_loop = tmp/nprocs;
    MPI_Allreduce(&cpu_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    cpu_loop = tmp/nprocs;
    if (time_loop > 0.0) cpu_loop = cpu_loop/time_loop*100.0;

    if (me == 0) {
      output->thermo->footer();
      int ntasks = nprocs * nthreads;
      utils::logmesg(lmp,"Loop time of {:.6g} on {} procs for {} steps with {} atoms\n\n",
                     time_loop,ntasks,update->nsteps,atom->natoms);

      // Gromacs/NAMD-style performance metric for suitable unit settings

      if ( timeflag && !minflag && !prdflag && !tadflag &&
           (update->nsteps > 0) && (update->dt != 0.0) &&
           ((strcmp(update->unit_style,"lj") == 0) ||
            (strcmp(update->unit_style,"metal") == 0) ||
            (strcmp(update->unit_style,"micro") == 0) ||
            (strcmp(update->unit_style,"nano") == 0) ||
            (strcmp(update->unit_style,"electron") == 0) ||
            (strcmp(update->unit_style,"real") == 0))) {
        double one_fs = force->femtosecond;
        double t_step = ((double) time_loop) / ((double) update->nsteps);
        double step_t = 1.0/t_step;

        if (strcmp(update->unit_style,"lj") == 0) {
          double tau_day = 24.0*3600.0 / t_step * update->dt / one_fs;
          utils::logmesg(lmp,"Performance: {:.3f} tau/day, {:.3f} timesteps/s\n",tau_day,step_t);
        } else if (strcmp(update->unit_style,"electron") == 0) {
          double hrs_fs = t_step / update->dt * one_fs / 3600.0;
          double fs_day = 24.0*3600.0 / t_step * update->dt / one_fs;
          utils::logmesg(lmp,"Performance: {:.3f} fs/day, {:.3f} hours/fs, "
                         "{:.3f} timesteps/s\n",fs_day,hrs_fs,step_t);
        } else {
          double hrs_ns = t_step / update->dt * 1000000.0 * one_fs / 3600.0;
          double ns_day = 24.0*3600.0 / t_step * update->dt / one_fs/1000000.0;
          utils::logmesg(lmp,"Performance: {:.3f} ns/day, {:.3f} hours/ns, "
                         "{:.3f} timesteps/s\n",ns_day,hrs_ns,step_t);
        }
      }

      // CPU use on MPI tasks and OpenMP threads

      if (timeflag) {
        if (lmp->kokkos) {
          utils::logmesg(lmp,"{:.1f}% CPU use with {} MPI tasks x {} OpenMP threads\n",
                         cpu_loop,nprocs,lmp->kokkos->nthreads);
        } else {
#if defined(_OPENMP)
          utils::logmesg(lmp,"{:.1f}% CPU use with {} MPI tasks x {} OpenMP threads\n",
                         cpu_loop,nprocs,nthreads);
#else
          utils::logmesg(lmp,"{:.1f}% CPU use with {} MPI tasks x no OpenMP threads\n",
                         cpu_loop,nprocs);
#endif
        }
      }
    }
  }

  // avoid division by zero for very short runs

  if (time_loop == 0.0) time_loop = 1.0;
  if (cpu_loop == 0.0) cpu_loop = 100.0;

  // get "Other" wall time for later use

  if (timer->has_normal())
    time_other = timer->get_wall(Timer::TOTAL) - timer->get_wall(Timer::ALL);

  // minimization stats

  if (minflag) {
    if (me == 0) {
      std::string mesg = "\nMinimization stats:\n";

      mesg += fmt::format("  Stopping criterion = {}\n",update->minimize->stopstr);
      mesg += fmt::format("  Energy initial, next-to-last, final = \n"
                          "    {:18.15g} {:18.15g} {:18.15g}\n",
                          update->minimize->einitial,
                          update->minimize->eprevious,
                          update->minimize->efinal);
      mesg += fmt::format("  Force two-norm initial, final = {:.8} {:.8}\n",
                          update->minimize->fnorm2_init,update->minimize->fnorm2_final);
      mesg += fmt::format("  Force max component initial, final = {:.8} {:.8}\n",
                          update->minimize->fnorminf_init,
                          update->minimize->fnorminf_final);
      mesg += fmt::format("  Final line search alpha, max atom move = {:.8} {:.8}\n",
                          update->minimize->alpha_final,
                          update->minimize->alpha_final*
                          update->minimize->fnorminf_final);
      mesg += fmt::format("  Iterations, force evaluations = {} {}\n",
                          update->minimize->niter,update->minimize->neval);
      utils::logmesg(lmp,mesg);
    }
  }

  // pair_style timing stats if provided

  if (force->pair) force->pair->finish();

  // PRD stats

  if (prdflag) {
    if (me == 0) utils::logmesg(lmp,"\nPRD stats:\n");

    time = timer->get_wall(Timer::DEPHASE);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Dephase  time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Dynamics time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Quench   time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::REPCOMM);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Comm     time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::REPOUT);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Output   time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Other    time (%) = {} ({})\n",time,time/time_loop*100.0);
  }

  // TAD stats

  if (tadflag) {
    if (me == 0) utils::logmesg(lmp,"\nTAD stats:\n");

    time = timer->get_wall(Timer::NEB);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  NEB      time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Dynamics time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Quench   time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::REPCOMM);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Comm     time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::REPOUT);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Output   time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Other    time (%) = {} ({})\n",time,time/time_loop*100.0);
  }

  // HYPER stats

  if (hyperflag) {
    if (me == 0) utils::logmesg(lmp,"\nHyper stats:\n");

    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Dynamics time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Quench   time (%) = {} ({})\n",time,time/time_loop*100.0);
    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) utils::logmesg(lmp,"  Other    time (%) = {} ({})\n",time,time/time_loop*100.0);
  }

  // further timing breakdowns

  if (timeflag && timer->has_normal()) {

    if (me == 0) {
      if (timer->has_full())
        utils::logmesg(lmp,"\nMPI task timing breakdown:\nSection |  min time "
                       " |  avg time  |  max time  |%varavg|  %CPU | %total\n"
                       "-----------------------------------------------------"
                       "------------------\n");
      else
        utils::logmesg(lmp,"\nMPI task timing breakdown:\nSection |  min time "
                       " |  avg time  |  max time  |%varavg| %total\n---------"
                       "------------------------------------------------------\n");
    }

    mpi_timings("Pair",timer,Timer::PAIR, world,nprocs,nthreads,me,time_loop,screen,logfile);

    if (atom->molecular != Atom::ATOMIC)
      mpi_timings("Bond",timer,Timer::BOND,world,nprocs,nthreads,me,time_loop,screen,logfile);

    if (force->kspace)
      mpi_timings("Kspace",timer,Timer::KSPACE,world,nprocs,nthreads,me,time_loop,screen,logfile);

    mpi_timings("Neigh",timer,Timer::NEIGH,world,nprocs,nthreads,me,time_loop,screen,logfile);
    mpi_timings("Comm",timer,Timer::COMM,world,nprocs,nthreads,me,time_loop,screen,logfile);
    mpi_timings("Output",timer,Timer::OUTPUT,world,nprocs,nthreads,me,time_loop,screen,logfile);
    mpi_timings("Modify",timer,Timer::MODIFY,world,nprocs,nthreads,me,time_loop,screen,logfile);
    if (timer->has_sync())
      mpi_timings("Sync",timer,Timer::SYNC,world,nprocs,nthreads,me,time_loop,screen,logfile);

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;

    if (me == 0) {
      if (timer->has_full())
        utils::logmesg(lmp,"Other   |            | {:<10.4g} |            |  "
                       "     |       |{:6.2f}\n",time,time/time_loop*100.0);
      else
        utils::logmesg(lmp,"Other   |            | {:<10.4g} |            |  "
                       "     |{:6.2f}\n",time,time/time_loop*100.0);
    }
  }

#ifdef LMP_OPENMP
  FixOMP *fixomp = dynamic_cast<FixOMP *>(modify->get_fix_by_id("package_omp"));

  // print thread breakdown only with full timer detail

  if (fixomp && timer->has_full() && me == 0) {
    double thr_total = 0.0;
    ThrData *td;
    for (i=0; i < nthreads; ++i) {
      td = fixomp->get_thr(i);
      thr_total += td->get_time(Timer::ALL);
    }
    thr_total /= (double) nthreads;

    if (thr_total > 0.0) {
      const std::string thr_fmt =
        "\nThread timing breakdown (MPI rank {}):\nTotal threaded time {:.4g} / {:.1f}%\n"
        "Section |  min time  |  avg time  |  max time  |%varavg| %total\n"
        "---------------------------------------------------------------\n";
      utils::logmesg(lmp,thr_fmt,me,thr_total,thr_total/time_loop*100.0);

      omp_times(fixomp,"Pair",Timer::PAIR,nthreads,screen,logfile);
      if (atom->molecular != Atom::ATOMIC)
        omp_times(fixomp,"Bond",Timer::BOND,nthreads,screen,logfile);
      if (force->kspace)
        omp_times(fixomp,"Kspace",Timer::KSPACE,nthreads,screen,logfile);
      omp_times(fixomp,"Neigh",Timer::NEIGH,nthreads,screen,logfile);
      omp_times(fixomp,"Reduce",Timer::COMM,nthreads,screen,logfile);
    }
  }
#endif

  if ((comm->me == 0) && lmp->kokkos && (lmp->kokkos->ngpus > 0))
    if (const char* env_clb = getenv("CUDA_LAUNCH_BLOCKING"))
      if (!(strcmp(env_clb,"1") == 0)) {
        error->warning(FLERR,"Timing breakdown may not be accurate "
                       "since GPU/CPU overlap is enabled\n"
                       "Using 'export CUDA_LAUNCH_BLOCKING=1' will give an "
                       "accurate timing breakdown but will reduce performance");
      }

  // FFT timing statistics
  // time3d,time1d = total time during run for 3d and 1d FFTs
  // loop on timing() until nsample FFTs require at least 1.0 CPU sec
  // time_kspace may be 0.0 if another partition is doing Kspace

  if (fftflag) {
    if (me == 0) utils::logmesg(lmp,"\n");
    int nsteps = update->nsteps;

    double time3d;
    int nsample = 1;
    int nfft = force->kspace->timing_3d(nsample,time3d);
    while (time3d < 1.0) {
      nsample *= 2;
      nfft = force->kspace->timing_3d(nsample,time3d);
    }

    time3d = nsteps * time3d / nsample;
    MPI_Allreduce(&time3d,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time3d = tmp/nprocs;

    double time1d;
    nsample = 1;
    nfft = force->kspace->timing_1d(nsample,time1d);
    while (time1d < 1.0) {
      nsample *= 2;
      nfft = force->kspace->timing_1d(nsample,time1d);
    }

    time1d = nsteps * time1d / nsample;
    MPI_Allreduce(&time1d,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time1d = tmp/nprocs;

    double time_kspace = timer->get_wall(Timer::KSPACE);
    MPI_Allreduce(&time_kspace,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_kspace = tmp/nprocs;

    double ntotal = 1.0 * force->kspace->nx_pppm *
      force->kspace->ny_pppm * force->kspace->nz_pppm;
    double nflops = 5.0 * ntotal * log(ntotal);

    double fraction,flop3,flop1;
    if (nsteps) {
      if (time_kspace) fraction = time3d/time_kspace*100.0;
      else fraction = 0.0;
      flop3 = nfft*nflops/1.0e9/(time3d/nsteps);
      flop1 = nfft*nflops/1.0e9/(time1d/nsteps);
    } else fraction = flop3 = flop1 = 0.0;

    if (me == 0)
      utils::logmesg(lmp,"FFT time (% of Kspce) = {:.6} ({:.4})\n"
                     "FFT Gflps 3d (1d only) = {:.8} {:.8}\n",
                     time3d,fraction,flop3,flop1);
  }

  nneigh = nneighfull = 0;
  if (histoflag) {
    std::string mesg = "\n";
    tmp = atom->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      mesg += fmt::format("Nlocal:    {:11.6} ave {:11.6g} max {:11.6g} min\n",ave,max,min);
      mesg += "Histogram:";
      for (i = 0; i < 10; i++) mesg += fmt::format(" {}",histo[i]);
      mesg += "\n";
    }

    tmp = atom->nghost;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      mesg += fmt::format("Nghost:    {:11.6} ave {:11.6g} max {:11.6g} min\n",ave,max,min);
      mesg += "Histogram:";
      for (i = 0; i < 10; i++) mesg += fmt::format(" {}",histo[i]);
      mesg += "\n";
    }

    tmp = nneigh = neighbor->get_nneigh_half();
    if (tmp < 0.0) tmp = 0.0;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      mesg += fmt::format("Neighs:    {:11.6} ave {:11.6g} max {:11.6g} min\n",ave,max,min);
      mesg += "Histogram:";
      for (i = 0; i < 10; i++) mesg += fmt::format(" {}",histo[i]);
      mesg += "\n";
    }

    tmp = nneighfull = neighbor->get_nneigh_full();
    if (tmp >= 0.0) {
      stats(1,&tmp,&ave,&max,&min,10,histo);
      if (me == 0) {
        mesg += fmt::format("FullNghs:  {:11.6} ave {:11.6g} max {:11.6g} min\n",ave,max,min);
        mesg += "Histogram:";
        for (i = 0; i < 10; i++) mesg += fmt::format(" {}",histo[i]);
        mesg += "\n";
      }
    }
    if (me == 0) utils::logmesg(lmp,mesg);
  }

  if (neighflag) {
    if (me == 0) utils::logmesg(lmp,"\n");

    tmp = MAX(MAX(nneigh,nneighfull),0.0);
    double nall;
    MPI_Allreduce(&tmp,&nall,1,MPI_DOUBLE,MPI_SUM,world);

    int nspec;
    double nspec_all = 0;
    if (atom->molecular == Atom::MOLECULAR) {
      int **nspecial = atom->nspecial;
      int nlocal = atom->nlocal;
      nspec = 0;
      for (i = 0; i < nlocal; i++) nspec += nspecial[i][2];
      tmp = nspec;
      MPI_Allreduce(&tmp,&nspec_all,1,MPI_DOUBLE,MPI_SUM,world);
    } else if (atom->molecular == Atom::TEMPLATE) {
      Molecule **onemols = atom->avec->onemols;
      int *molindex = atom->molindex;
      int *molatom = atom->molatom;
      int nlocal = atom->nlocal;
      int imol,iatom;
      nspec = 0;
      for (i = 0; i < nlocal; i++) {
        if (molindex[i] < 0) continue;
        imol = molindex[i];
        iatom = molatom[i];
        nspec += onemols[imol]->nspecial[iatom][2];
      }
      tmp = nspec;
      MPI_Allreduce(&tmp,&nspec_all,1,MPI_DOUBLE,MPI_SUM,world);
    }

    if (me == 0) {
      std::string mesg;

      mesg += fmt::format("Total # of neighbors = {:.8g}\n",nall);
      if (atom->natoms > 0)
        mesg += fmt::format("Ave neighs/atom = {:.8}\n",nall/atom->natoms);
      if ((atom->molecular != Atom::ATOMIC) && (atom->natoms > 0))
        mesg += fmt::format("Ave special neighs/atom = {:.8}\n",nspec_all/atom->natoms);
      mesg += fmt::format("Neighbor list builds = {}\n",neighbor->ncalls);
      if (neighbor->dist_check)
        mesg += fmt::format("Dangerous builds = {}\n",neighbor->ndanger);
      else mesg += "Dangerous builds not checked\n";
      utils::logmesg(lmp,mesg);
    }
  }

  if (logfile) fflush(logfile);
}

/* ---------------------------------------------------------------------- */

void Finish::stats(int n, double *data, double *pave, double *pmax, double *pmin, int nhisto, int *histo)
{
  int i,m;
  int *histotmp;

  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;
  for (i = 0; i < n; i++) {
    ave += data[i];
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  int ntotal;
  MPI_Allreduce(&n,&ntotal,1,MPI_INT,MPI_SUM,world);
  double tmp;
  MPI_Allreduce(&ave,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  ave = tmp/ntotal;
  MPI_Allreduce(&min,&tmp,1,MPI_DOUBLE,MPI_MIN,world);
  min = tmp;
  MPI_Allreduce(&max,&tmp,1,MPI_DOUBLE,MPI_MAX,world);
  max = tmp;

  for (i = 0; i < nhisto; i++) histo[i] = 0;

  double del = max - min;
  for (i = 0; i < n; i++) {
    if (del == 0.0) m = 0;
    else m = static_cast<int> ((data[i]-min)/del * nhisto);
    if (m > nhisto-1) m = nhisto-1;
    histo[m]++;
  }

  memory->create(histotmp,nhisto,"finish:histotmp");
  MPI_Allreduce(histo,histotmp,nhisto,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
  memory->destroy(histotmp);

  *pave = ave;
  *pmax = max;
  *pmin = min;
}

/* ---------------------------------------------------------------------- */

void mpi_timings(const char *label, Timer *t, enum Timer::ttype tt,
                        MPI_Comm world, const int nprocs, const int nthreads,
                        const int me, double time_loop, FILE *scr, FILE *log)
{
  double tmp, time_max, time_min, time_sq;
  double time = t->get_wall(tt);

  double time_cpu = t->get_cpu(tt);
  if (time/time_loop < 0.001)  // insufficient timer resolution!
    time_cpu = 1.0;
  else
    time_cpu = time_cpu / time;
  if (time_cpu > nthreads) time_cpu = nthreads;

  MPI_Allreduce(&time,&time_min,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&time,&time_max,1,MPI_DOUBLE,MPI_MAX,world);
  time_sq = time*time;
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  MPI_Allreduce(&time_sq,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_sq = tmp/nprocs;
  MPI_Allreduce(&time_cpu,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_cpu = tmp/nprocs*100.0;

  // % variance from the average as measure of load imbalance
  if ((time > 0.001) && ((time_sq/time - time) > 1.0e-10))
    time_sq = sqrt(time_sq/time - time)*100.0;
  else
    time_sq = 0.0;

  if (me == 0) {
    tmp = time/time_loop*100.0;
    std::string mesg;
    if (t->has_full())
      mesg = fmt::format("{:<8s}| {:<10.5g} | {:<10.5g} | {:<10.5g} |{:6.1f} |"
                         "{:6.1f} |{:6.2f}\n",
                         label,time_min,time,time_max,time_sq,time_cpu,tmp);
    else
      mesg = fmt::format("{:<8s}| {:<10.5g} | {:<10.5g} | {:<10.5g} |{:6.1f} |"
                         "{:6.2f}\n",label,time_min,time,time_max,time_sq,tmp);
    if (scr) fputs(mesg.c_str(),scr);
    if (log) fputs(mesg.c_str(),log);
  }
}

/* ---------------------------------------------------------------------- */

#ifdef LMP_OPENMP
void omp_times(FixOMP *fix, const char *label, enum Timer::ttype which,
                      const int nthreads,FILE *scr, FILE *log)
{
  double time_min, time_max, time_avg, time_total, time_std;

  time_min =  1.0e100;
  time_max = -1.0e100;
  time_total = time_avg = time_std = 0.0;

  for (int i=0; i < nthreads; ++i) {
    ThrData *thr = fix->get_thr(i);
    double tmp=thr->get_time(which);
    time_min = MIN(time_min,tmp);
    time_max = MAX(time_max,tmp);
    time_avg += tmp;
    time_std += tmp*tmp;
    time_total += thr->get_time(Timer::ALL);
  }

  time_avg /= nthreads;
  time_std /= nthreads;
  time_total /= nthreads;

  if ((time_avg > 0.001) && ((time_std/time_avg -time_avg) > 1.0e-10))
    time_std = sqrt(time_std/time_avg - time_avg)*100.0;
  else
    time_std = 0.0;

  std::string mesg = fmt::format("{:<8s}| {:10.5g} | {:10.5g} | {:10.5g} |"
                                 "{:6.1f} |{:6.2f}\n",label,time_min,time_avg,
                                 time_max,time_std,time_avg/time_total*100.0);
  if (scr) fputs(mesg.c_str(),scr);
  if (log) fputs(mesg.c_str(),log);
}
#endif

