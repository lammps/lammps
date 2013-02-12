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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "finish.h"
#include "timer.h"
#include "universe.h"
#include "atom.h"
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Finish::Finish(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Finish::end(int flag)
{
  int i,m,nneigh,nneighfull;
  int histo[10];
  int loopflag,minflag,prdflag,tadflag,timeflag,fftflag,histoflag,neighflag;
  double time,tmp,ave,max,min;
  double time_loop,time_other;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // recompute natoms in case atoms have been lost

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // choose flavors of statistical output
  // flag determines caller
  // flag = 0 = just loop summary
  // flag = 1 = dynamics or minimization
  // flag = 2 = PRD
  // flag = 3 = TAD
  // turn off neighflag for Kspace partition of verlet/split integrator

  loopflag = 1;
  minflag = prdflag = tadflag = timeflag = fftflag = histoflag = neighflag = 0;

  if (flag == 1) {
    if (update->whichflag == 2) minflag = 1;
    timeflag = histoflag = 1;
    neighflag = 1;
    if (update->whichflag == 1 &&
        strcmp(update->integrate_style,"verlet/split") == 0 &&
        universe->iworld == 1) neighflag = 0;
    if (force->kspace && force->kspace_match("pppm",0)
        && force->kspace->fftbench) fftflag = 1;
  }
  if (flag == 2) prdflag = histoflag = neighflag = 1;
  if (flag == 3) tadflag = histoflag = neighflag = 1;

  // loop stats

  if (loopflag) {
    time_other = timer->array[TIME_LOOP] -
      (timer->array[TIME_PAIR] + timer->array[TIME_BOND] +
       timer->array[TIME_KSPACE] + timer->array[TIME_NEIGHBOR] +
       timer->array[TIME_COMM] + timer->array[TIME_OUTPUT]);

    time_loop = timer->array[TIME_LOOP];
    MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_loop = tmp/nprocs;

    // overall loop time

#if defined(_OPENMP)
    if (me == 0) {
      int ntasks = nprocs * comm->nthreads;
      if (screen) fprintf(screen,
                          "Loop time of %g on %d procs (%d MPI x %d OpenMP) "
                          "for %d steps with " BIGINT_FORMAT " atoms\n",
                          time_loop,ntasks,nprocs,comm->nthreads,
                          update->nsteps,atom->natoms);
      if (logfile) fprintf(logfile,
                          "Loop time of %g on %d procs (%d MPI x %d OpenMP) "
                          "for %d steps with " BIGINT_FORMAT " atoms\n",
                          time_loop,ntasks,nprocs,comm->nthreads,
                          update->nsteps,atom->natoms);
    }
#else
    if (me == 0) {
      if (screen) fprintf(screen,
                          "Loop time of %g on %d procs for %d steps with "
                          BIGINT_FORMAT " atoms\n",
                          time_loop,nprocs,update->nsteps,atom->natoms);
      if (logfile) fprintf(logfile,
                           "Loop time of %g on %d procs for %d steps with "
                           BIGINT_FORMAT " atoms\n",
                           time_loop,nprocs,update->nsteps,atom->natoms);
    }
#endif

    if (time_loop == 0.0) time_loop = 1.0;
  }

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

  // PRD stats using PAIR,BOND,KSPACE for dephase,dynamics,quench

  if (prdflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    if (screen) fprintf(screen,"PRD stats:\n");
    if (logfile) fprintf(logfile,"PRD stats:\n");

    time = timer->array[TIME_PAIR];
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

    time = timer->array[TIME_BOND];
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

    time = timer->array[TIME_KSPACE];
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

  // TAD stats using PAIR,BOND,KSPACE for neb,dynamics,quench

  if (tadflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    if (screen) fprintf(screen,"TAD stats:\n");
    if (logfile) fprintf(logfile,"TAD stats:\n");

    time = timer->array[TIME_PAIR];
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

    time = timer->array[TIME_BOND];
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

    time = timer->array[TIME_KSPACE];
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


    time = timer->array[TIME_COMM];
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


    time = timer->array[TIME_OUTPUT];
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

  // timing breakdowns

  if (timeflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    time = timer->array[TIME_PAIR];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"Pair  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"Pair  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    if (atom->molecular) {
      time = timer->array[TIME_BOND];
      MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
      time = tmp/nprocs;
      if (me == 0) {
        if (screen)
          fprintf(screen,"Bond  time (%%) = %g (%g)\n",
                  time,time/time_loop*100.0);
        if (logfile)
          fprintf(logfile,"Bond  time (%%) = %g (%g)\n",
                  time,time/time_loop*100.0);
      }
    }

    if (force->kspace) {
      time = timer->array[TIME_KSPACE];
      MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
      time = tmp/nprocs;
      if (me == 0) {
        if (screen)
          fprintf(screen,"Kspce time (%%) = %g (%g)\n",
                  time,time/time_loop*100.0);
        if (logfile)
          fprintf(logfile,"Kspce time (%%) = %g (%g)\n",
                  time,time/time_loop*100.0);
      }
    }

    time = timer->array[TIME_NEIGHBOR];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"Neigh time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"Neigh time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->array[TIME_COMM];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"Comm  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"Comm  time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = timer->array[TIME_OUTPUT];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"Outpt time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"Outpt time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen)
        fprintf(screen,"Other time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
      if (logfile)
        fprintf(logfile,"Other time (%%) = %g (%g)\n",
                time,time/time_loop*100.0);
    }
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

    double time_kspace = timer->array[TIME_KSPACE];
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

    // find a non-skip neighbor list containing half the pairwise interactions
    // count neighbors in that list for stats purposes

    for (m = 0; m < neighbor->old_nrequest; m++)
      if ((neighbor->old_requests[m]->half ||
           neighbor->old_requests[m]->gran ||
           neighbor->old_requests[m]->respaouter ||
           neighbor->old_requests[m]->half_from_full) &&
          neighbor->old_requests[m]->skip == 0 &&
          neighbor->lists[m]->numneigh) break;

    nneigh = 0;
    if (m < neighbor->old_nrequest) {
      int inum = neighbor->lists[m]->inum;
      int *ilist = neighbor->lists[m]->ilist;
      int *numneigh = neighbor->lists[m]->numneigh;
      for (i = 0; i < inum; i++)
        nneigh += numneigh[ilist[i]];
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

    for (m = 0; m < neighbor->old_nrequest; m++)
      if (neighbor->old_requests[m]->full &&
          neighbor->old_requests[m]->skip == 0) break;

    nneighfull = 0;
    if (m < neighbor->old_nrequest) {
      if (neighbor->lists[m]->numneigh > 0) {
        int inum = neighbor->lists[m]->inum;
        int *ilist = neighbor->lists[m]->ilist;
        int *numneigh = neighbor->lists[m]->numneigh;
        for (i = 0; i < inum; i++)
          nneighfull += numneigh[ilist[i]];
      }

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
    double nspec_all;
    if (atom->molecular) {
      nspec = 0;
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++) nspec += atom->nspecial[i][2];
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
        fprintf(screen,"Dangerous builds = " BIGINT_FORMAT "\n",
                neighbor->ndanger);
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
        fprintf(logfile,"Dangerous builds = " BIGINT_FORMAT "\n",
                neighbor->ndanger);
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
