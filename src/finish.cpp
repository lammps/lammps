/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "finish.h"
#include "timer.h"
#include "universe.h"
#include "accelerator_kokkos.h"
#include "accelerator_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "min.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "output.h"
#include "memory.h"
#include "error.h"

#ifdef LMP_USER_OMP
#include "modify.h"
#include "fix_omp.h"
#include "thr_data.h"
#endif

using namespace LAMMPS_NS;

// local function prototypes, code at end of file

static void mpi_timings(const char *label, Timer *t, enum Timer::ttype tt,
                        MPI_Comm world, const int nprocs, const int nthreads,
                        const int me, double time_loop, FILE *scr, FILE *log);

#ifdef LMP_USER_OMP
static void omp_times(FixOMP *fix, const char *label, enum Timer::ttype which,
                      const int nthreads,FILE *scr, FILE *log);
#endif

/* ---------------------------------------------------------------------- */

Finish::Finish(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Finish::end(int flag)
{
  int i,m,nneigh,nneighfull;
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
    if (force->kspace && force->kspace_match("pppm",0)
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
      int ntasks = nprocs * nthreads;
      const char fmt1[] = "Loop time of %g on %d procs "
        "for %d steps with " BIGINT_FORMAT " atoms\n\n";
      if (screen) fprintf(screen,fmt1,time_loop,ntasks,update->nsteps,
                          atom->natoms);
      if (logfile) fprintf(logfile,fmt1,time_loop,ntasks,update->nsteps,
                           atom->natoms);

      // Gromacs/NAMD-style performance metric for suitable unit settings

      if ( timeflag && !minflag && !prdflag && !tadflag &&
           (update->nsteps > 0) && (update->dt != 0.0) &&
           ((strcmp(update->unit_style,"lj") == 0) ||
            (strcmp(update->unit_style,"metal") == 0) ||
            (strcmp(update->unit_style,"micro") == 0) ||
            (strcmp(update->unit_style,"nano") == 0) ||
            (strcmp(update->unit_style,"electron") == 0) ||
            (strcmp(update->unit_style,"real") == 0)) ) {
        double one_fs = force->femtosecond;
        double t_step = ((double) time_loop) / ((double) update->nsteps);
        double step_t = 1.0/t_step;

        if (strcmp(update->unit_style,"lj") == 0) {
          double tau_day = 24.0*3600.0 / t_step * update->dt / one_fs;
          const char perf[] = "Performance: %.3f tau/day, %.3f timesteps/s\n";
          if (screen) fprintf(screen,perf,tau_day,step_t);
          if (logfile) fprintf(logfile,perf,tau_day,step_t);
        } else if (strcmp(update->unit_style,"electron") == 0) {
          double hrs_fs = t_step / update->dt * one_fs / 3600.0;
          double fs_day = 24.0*3600.0 / t_step * update->dt / one_fs;
          const char perf[] =
            "Performance: %.3f fs/day, %.3f hours/fs, %.3f timesteps/s\n";
          if (screen) fprintf(screen,perf,fs_day,hrs_fs,step_t);
          if (logfile) fprintf(logfile,perf,fs_day,hrs_fs,step_t);

        } else {
          double hrs_ns = t_step / update->dt * 1000000.0 * one_fs / 3600.0;
          double ns_day = 24.0*3600.0 / t_step * update->dt / one_fs/1000000.0;
          const char perf[] =
            "Performance: %.3f ns/day, %.3f hours/ns, %.3f timesteps/s\n";
          if (screen) fprintf(screen,perf,ns_day,hrs_ns,step_t);
          if (logfile) fprintf(logfile,perf,ns_day,hrs_ns,step_t);
        }
      }

      // CPU use on MPI tasks and OpenMP threads

      if (timeflag) {
        if (lmp->kokkos) {
          const char fmt2[] =
            "%.1f%% CPU use with %d MPI tasks x %d OpenMP threads\n";
          if (screen) fprintf(screen,fmt2,cpu_loop,nprocs,
                              lmp->kokkos->num_threads);
          if (logfile) fprintf(logfile,fmt2,cpu_loop,nprocs,
                               lmp->kokkos->num_threads);
        } else {
#if defined(_OPENMP)
          const char fmt2[] =
            "%.1f%% CPU use with %d MPI tasks x %d OpenMP threads\n";
          if (screen) fprintf(screen,fmt2,cpu_loop,nprocs,nthreads);
          if (logfile) fprintf(logfile,fmt2,cpu_loop,nprocs,nthreads);
#else
          const char fmt2[] =
            "%.1f%% CPU use with %d MPI tasks x no OpenMP threads\n";
          if (screen) fprintf(screen,fmt2,cpu_loop,nprocs);
          if (logfile) fprintf(logfile,fmt2,cpu_loop,nprocs);
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
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    if (me == 0) {
      if (screen) {
        fprintf(screen,"Minimization stats:\n");
        fprintf(screen,"  Stopping criterion = %s\n",
                update->minimize->stopstr);
        fprintf(screen,"  Energy initial, next-to-last, final = \n"
                "    %18.12g %18.12g %18.12g\n",
                update->minimize->einitial,update->minimize->eprevious,
                update->minimize->efinal);
        fprintf(screen,"  Force two-norm initial, final = %g %g\n",
                update->minimize->fnorm2_init,update->minimize->fnorm2_final);
        fprintf(screen,"  Force max component initial, final = %g %g\n",
                update->minimize->fnorminf_init,
                update->minimize->fnorminf_final);
        fprintf(screen,"  Final line search alpha, max atom move = %g %g\n",
                update->minimize->alpha_final,
                update->minimize->alpha_final*
                update->minimize->fnorminf_final);
        fprintf(screen,"  Iterations, force evaluations = %d %d\n",
                update->minimize->niter,update->minimize->neval);
      }
      if (logfile) {
        fprintf(logfile,"Minimization stats:\n");
        fprintf(logfile,"  Stopping criterion = %s\n",
                update->minimize->stopstr);
        fprintf(logfile,"  Energy initial, next-to-last, final = \n"
                "    %18.12g %18.12g %18.12g\n",
                update->minimize->einitial,update->minimize->eprevious,
                update->minimize->efinal);
        fprintf(logfile,"  Force two-norm initial, final = %g %g\n",
                update->minimize->fnorm2_init,update->minimize->fnorm2_final);
        fprintf(logfile,"  Force max component initial, final = %g %g\n",
                update->minimize->fnorminf_init,
                update->minimize->fnorminf_final);
        fprintf(logfile,"  Final line search alpha, max atom move = %g %g\n",
                update->minimize->alpha_final,
                update->minimize->alpha_final*
                update->minimize->fnorminf_final);
        fprintf(logfile,"  Iterations, force evaluations = %d %d\n",
                update->minimize->niter,update->minimize->neval);
      }
    }
  }

  // PRD stats

  if (prdflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\nPRD stats:\n");
      if (logfile) fprintf(logfile,"\nPRD stats:\n");
    }

    time = timer->get_wall(Timer::DEPHASE);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Dephase  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Dephase  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

      time = timer->get_wall(Timer::REPCOMM);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Comm     time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Comm     time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }


    time = timer->get_wall(Timer::REPOUT);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Output   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Output   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) { // XXXX: replica comm, replica output
      if (screen)
        fprintf(screen,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }
  }

  // TAD stats

  if (tadflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    if (screen) fprintf(screen,"TAD stats:\n");
    if (logfile) fprintf(logfile,"TAD stats:\n");

    time = timer->get_wall(Timer::NEB);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  NEB      time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  NEB      time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }


    time = timer->get_wall(Timer::REPCOMM);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Comm     time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Comm     time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }


    time = timer->get_wall(Timer::REPOUT);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Output   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Output   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }
  }

  // HYPER stats

  if (hyperflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\nHyper stats:\n");
      if (logfile) fprintf(logfile,"\nHyper stats:\n");
    }

    time = timer->get_wall(Timer::DYNAMICS);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Dynamics time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->get_wall(Timer::QUENCH);
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Quench   time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"  Other    time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }
  }

  // further timing breakdowns

  if (timeflag && timer->has_normal()) {

    if (timer->has_full()) {
      const char hdr[] = "\nMPI task timing breakdown:\n"
        "Section |  min time  |  avg time  |  max time  "
        "|%varavg|  %CPU | %total\n"
        "-----------------------------------------------"
        "------------------------\n";
      if (me == 0) {
        if (screen)  fputs(hdr,screen);
        if (logfile) fputs(hdr,logfile);
      }
    } else {
      const char hdr[] = "\nMPI task timing breakdown:\n"
        "Section |  min time  |  avg time  |  max time  |%varavg| %total\n"
        "---------------------------------------------------------------\n";
      if (me == 0) {
        if (screen)  fputs(hdr,screen);
        if (logfile) fputs(hdr,logfile);
      }
    }

    mpi_timings("Pair",timer,Timer::PAIR, world,nprocs,
                nthreads,me,time_loop,screen,logfile);

    if (atom->molecular)
      mpi_timings("Bond",timer,Timer::BOND,world,nprocs,
                  nthreads,me,time_loop,screen,logfile);

    if (force->kspace)
      mpi_timings("Kspace",timer,Timer::KSPACE,world,nprocs,
                  nthreads,me,time_loop,screen,logfile);

    mpi_timings("Neigh",timer,Timer::NEIGH,world,nprocs,
                nthreads,me,time_loop,screen,logfile);
    mpi_timings("Comm",timer,Timer::COMM,world,nprocs,
                nthreads,me,time_loop,screen,logfile);
    mpi_timings("Output",timer,Timer::OUTPUT,world,nprocs,
                nthreads,me,time_loop,screen,logfile);
    mpi_timings("Modify",timer,Timer::MODIFY,world,nprocs,
                nthreads,me,time_loop,screen,logfile);
    if (timer->has_sync())
      mpi_timings("Sync",timer,Timer::SYNC,world,nprocs,
                  nthreads,me,time_loop,screen,logfile);

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;

    const char *fmt;
    if (timer->has_full())
      fmt = "Other   |            |%- 12.4g|            |       |       |%6.2f\n";
    else
      fmt = "Other   |            |%- 12.4g|            |       |%6.2f\n";

    if (me == 0) {
      if (screen) fprintf(screen,fmt,time,time/time_loop*100.0);
      if (logfile) fprintf(logfile,fmt,time,time/time_loop*100.0);
    }
  }

#ifdef LMP_USER_OMP
  const char thr_hdr_fmt[] =
    "\nThread timing breakdown (MPI rank %d):\nTotal threaded time %.4g / %.1f%%\n";
  const char thr_header[] =
    "Section |  min time  |  avg time  |  max time  |%varavg| %total\n"
    "---------------------------------------------------------------\n";

  int ifix = modify->find_fix("package_omp");

  // print thread breakdown only with full timer detail

  if ((ifix >= 0) && timer->has_full() && me == 0) {
    double thr_total = 0.0;
    ThrData *td;
    FixOMP *fixomp = static_cast<FixOMP *>(lmp->modify->fix[ifix]);
    for (i=0; i < nthreads; ++i) {
      td = fixomp->get_thr(i);
      thr_total += td->get_time(Timer::ALL);
    }
    thr_total /= (double) nthreads;

    if (thr_total > 0.0) {
      if (screen) {
        fprintf(screen,thr_hdr_fmt,me,thr_total,thr_total/time_loop*100.0);
        fputs(thr_header,screen);
      }
      if (logfile) {
        fprintf(logfile,thr_hdr_fmt,me,thr_total,thr_total/time_loop*100.0);
        fputs(thr_header,logfile);
      }

      omp_times(fixomp,"Pair",Timer::PAIR,nthreads,screen,logfile);

      if (atom->molecular)
        omp_times(fixomp,"Bond",Timer::BOND,nthreads,screen,logfile);

      if (force->kspace)
        omp_times(fixomp,"Kspace",Timer::KSPACE,nthreads,screen,logfile);

      omp_times(fixomp,"Neigh",Timer::NEIGH,nthreads,screen,logfile);
      omp_times(fixomp,"Reduce",Timer::COMM,nthreads,screen,logfile);
    }
  }
#endif

  if (lmp->kokkos && lmp->kokkos->ngpu > 0)
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
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

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

    if (me == 0) {
      if (screen) {
        fprintf(screen,"FFT time (%% of Kspce) = %g (%g)\n",time3d,fraction);
        fprintf(screen,"FFT Gflps 3d (1d only) = %g %g\n",flop3,flop1);
      }
      if (logfile) {
        fprintf(logfile,"FFT time (%% of Kspce) = %g (%g)\n",time3d,fraction);
        fprintf(logfile,"FFT Gflps 3d (1d only) = %g %g\n",flop3,flop1);
      }
    }
  }

  if (histoflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    tmp = atom->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    tmp = atom->nghost;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"Nghost:    %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"Nghost:    %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    // find a non-skip neighbor list containing half pairwise interactions
    // count neighbors in that list for stats purposes
    // allow it to be Kokkos neigh list as well

    for (m = 0; m < neighbor->old_nrequest; m++)
      if (neighbor->old_requests[m]->half &&
          neighbor->old_requests[m]->skip == 0 &&
          neighbor->lists[m] && neighbor->lists[m]->numneigh) break;

    nneigh = 0;
    if (m < neighbor->old_nrequest) {
      if (!neighbor->lists[m]->kokkos) {
        int inum = neighbor->lists[m]->inum;
        int *ilist = neighbor->lists[m]->ilist;
        int *numneigh = neighbor->lists[m]->numneigh;
        for (i = 0; i < inum; i++)
          nneigh += numneigh[ilist[i]];
      } else if (lmp->kokkos) nneigh = lmp->kokkos->neigh_count(m);
    }

    tmp = nneigh;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"Neighs:    %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"Neighs:    %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    // find a non-skip neighbor list containing full pairwise interactions
    // count neighbors in that list for stats purposes
    // allow it to be Kokkos neigh list as well

    for (m = 0; m < neighbor->old_nrequest; m++)
      if (neighbor->old_requests[m]->full &&
          neighbor->old_requests[m]->skip == 0) break;

    nneighfull = 0;
    if (m < neighbor->old_nrequest) {
      if (!neighbor->lists[m]->kokkos && neighbor->lists[m]->numneigh) {
        int inum = neighbor->lists[m]->inum;
        int *ilist = neighbor->lists[m]->ilist;
        int *numneigh = neighbor->lists[m]->numneigh;
        for (i = 0; i < inum; i++)
          nneighfull += numneigh[ilist[i]];
      } else if (lmp->kokkos)
          nneighfull = lmp->kokkos->neigh_count(m);

      tmp = nneighfull;
      stats(1,&tmp,&ave,&max,&min,10,histo);
      if (me == 0) {
        if (screen) {
          fprintf(screen,"FullNghs:  %g ave %g max %g min\n",ave,max,min);
          fprintf(screen,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
          fprintf(screen,"\n");
        }
        if (logfile) {
          fprintf(logfile,"FullNghs:  %g ave %g max %g min\n",ave,max,min);
          fprintf(logfile,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
          fprintf(logfile,"\n");
        }
      }
    }
  }

  if (neighflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    tmp = MAX(nneigh,nneighfull);
    double nall;
    MPI_Allreduce(&tmp,&nall,1,MPI_DOUBLE,MPI_SUM,world);

    int nspec;
    double nspec_all = 0;
    if (atom->molecular == 1) {
      int **nspecial = atom->nspecial;
      int nlocal = atom->nlocal;
      nspec = 0;
      for (i = 0; i < nlocal; i++) nspec += nspecial[i][2];
      tmp = nspec;
      MPI_Allreduce(&tmp,&nspec_all,1,MPI_DOUBLE,MPI_SUM,world);
    } else if (atom->molecular == 2) {
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
      if (screen) {
        if (nall < 2.0e9)
          fprintf(screen,
                  "Total # of neighbors = %d\n",static_cast<int> (nall));
        else fprintf(screen,"Total # of neighbors = %g\n",nall);
        if (atom->natoms > 0)
          fprintf(screen,"Ave neighs/atom = %g\n",nall/atom->natoms);
        if (atom->molecular && atom->natoms > 0)
          fprintf(screen,"Ave special neighs/atom = %g\n",
                  nspec_all/atom->natoms);
        fprintf(screen,"Neighbor list builds = " BIGINT_FORMAT "\n",
                neighbor->ncalls);
        if (neighbor->dist_check)
          fprintf(screen,"Dangerous builds = " BIGINT_FORMAT "\n",
                  neighbor->ndanger);
        else fprintf(screen,"Dangerous builds not checked\n");
      }
      if (logfile) {
        if (nall < 2.0e9)
          fprintf(logfile,
                  "Total # of neighbors = %d\n",static_cast<int> (nall));
        else fprintf(logfile,"Total # of neighbors = %g\n",nall);
        if (atom->natoms > 0)
          fprintf(logfile,"Ave neighs/atom = %g\n",nall/atom->natoms);
        if (atom->molecular && atom->natoms > 0)
          fprintf(logfile,"Ave special neighs/atom = %g\n",
                  nspec_all/atom->natoms);
        fprintf(logfile,"Neighbor list builds = " BIGINT_FORMAT "\n",
                neighbor->ncalls);
        if (neighbor->dist_check)
          fprintf(logfile,"Dangerous builds = " BIGINT_FORMAT "\n",
                  neighbor->ndanger);
        else fprintf(logfile,"Dangerous builds not checked\n");
      }
    }
  }

  if (logfile) fflush(logfile);
}

/* ---------------------------------------------------------------------- */

void Finish::stats(int n, double *data,
                   double *pave, double *pmax, double *pmin,
                   int nhisto, int *histo)
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
    if (t->has_full()) {
      const char fmt[] = "%-8s|%- 12.5g|%- 12.5g|%- 12.5g|%6.1f |%6.1f |%6.2f\n";
      if (scr)
        fprintf(scr,fmt,label,time_min,time,time_max,time_sq,time_cpu,tmp);
      if (log)
        fprintf(log,fmt,label,time_min,time,time_max,time_sq,time_cpu,tmp);
      time_loop = 100.0/time_loop;
    } else {
      const char fmt[] = "%-8s|%- 12.5g|%- 12.5g|%- 12.5g|%6.1f |%6.2f\n";
      if (scr)
        fprintf(scr,fmt,label,time_min,time,time_max,time_sq,tmp);
      if (log)
        fprintf(log,fmt,label,time_min,time,time_max,time_sq,tmp);
    }
  }
}

/* ---------------------------------------------------------------------- */

#ifdef LMP_USER_OMP
void omp_times(FixOMP *fix, const char *label, enum Timer::ttype which,
                      const int nthreads,FILE *scr, FILE *log)
{
  const char fmt[] = "%-8s|%- 12.5g|%- 12.5g|%- 12.5g|%6.1f |%6.2f\n";
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

  if (scr) fprintf(scr,fmt,label,time_min,time_avg,time_max,time_std,
                   time_avg/time_total*100.0);
  if (log) fprintf(log,fmt,label,time_min,time_avg,time_max,time_std,
                   time_avg/time_total*100.0);
}
#endif

