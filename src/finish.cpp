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

static void print_timings(const char *label, Timer *t, enum Timer::ttype tt,
			  MPI_Comm world, const int nprocs, const int me,
			  double time_loop, FILE *screen, FILE *logfile)
{
  double tmp;
  double time = t->get_wall(tt);

#if defined(_OPENMP)
  double time_cpu = t->get_cpu(tt);
  if (time == 0.0)
    time_cpu = 0.0;
  else
    time_cpu = time_cpu / time * 100.0;
#endif

  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;

#if defined(_OPENMP)
  MPI_Allreduce(&time_cpu,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_cpu = tmp/nprocs;

  if (me == 0) {
    if (screen)
      fprintf(screen,"%5s time (%%) = %g (%g)  (%%CPU = %g)\n",
	      label,time,time/time_loop*100.0,time_cpu);
    if (logfile)
      fprintf(logfile,"%5s  time (%%) = %g (%g)  (%%CPU = %g)\n",
	      label,time,time/time_loop*100.0,time_cpu);
  }
#else
  if (me == 0) {
    if (screen)
      fprintf(screen,"%5s time (%%) = %g (%g)\n",
	      label,time,time/time_loop*100.0);
    if (logfile)
      fprintf(logfile,"%5s time (%%) = %g (%g)\n",
	      label,time,time/time_loop*100.0);
  }
#endif
}

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

  // choose flavors of statistical output
  // flag determines caller
  // flag = 0 = just loop summary
  // flag = 1 = dynamics or minimization
  // flag = 2 = PRD
  // flag = 3 = TAD
  
  loopflag = 1;
  minflag = prdflag = tadflag = timeflag = fftflag = histoflag = neighflag = 0;

  if (flag == 1) {
    if (update->whichflag == 2) minflag = 1;
    timeflag = histoflag = neighflag = 1;
    if (strstr(force->kspace_style,"pppm")) fftflag = 1;
  }
  if (flag == 2) prdflag = histoflag = neighflag = 1;
  if (flag == 3) tadflag = histoflag = neighflag = 1;

  // loop stats

  if (loopflag) {
    time_other = timer->get_wall(Timer::LOOP) -
      (timer->get_wall(Timer::PAIR) + timer->get_wall(Timer::BOND) + 
       timer->get_wall(Timer::KSPACE) + timer->get_wall(Timer::NEIGHBOR) +
       timer->get_wall(Timer::COMM) + timer->get_wall(Timer::OUTPUT) +
       timer->get_wall(Timer::MODIFY));
    
    time_loop = timer->get_wall(Timer::LOOP);
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

    time = timer->get_wall(Timer::PAIR);
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

    time = timer->get_wall(Timer::BOND);
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

    time = timer->get_wall(Timer::KSPACE);
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

    time = timer->get_wall(Timer::PAIR);
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

    time = timer->get_wall(Timer::BOND);
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

    time = timer->get_wall(Timer::KSPACE);
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


    time = timer->get_wall(Timer::COMM);
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


    time = timer->get_wall(Timer::OUTPUT);
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

    print_timings("Pair ",timer,Timer::PAIR,world,nprocs,me,time_loop,screen,logfile);

    if (atom->molecular)
      print_timings("Bond ",timer,Timer::BOND,world,nprocs,me,time_loop,screen,logfile);
    
    if (force->kspace)
      print_timings("Kspce",timer,Timer::KSPACE,world,nprocs,me,time_loop,screen,logfile);

    print_timings("Neigh",timer,Timer::NEIGHBOR,world,nprocs,me,time_loop,screen,logfile);
    print_timings("Comm ",timer,Timer::COMM,    world,nprocs,me,time_loop,screen,logfile);
    print_timings("Outpt",timer,Timer::OUTPUT,  world,nprocs,me,time_loop,screen,logfile);
    print_timings("Modfy",timer,Timer::MODIFY,  world,nprocs,me,time_loop,screen,logfile);

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

  if (fftflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    int nsteps = update->nsteps;

    int nsample = 5;
    double time3d,time1d;
    force->kspace->timing(nsample,time3d,time1d);
    
    time3d = nsteps * time3d / nsample;
    MPI_Allreduce(&time3d,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time3d = tmp/nprocs;
    
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
      fraction = time3d/time_kspace*100.0;
      flop3 = nflops/1.0e9/(time3d/4.0/nsteps);
      flop1 = nflops/1.0e9/(time1d/4.0/nsteps);
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
	  fprintf(logfile,"FullNghs: %g ave %g max %g min\n",ave,max,min);
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
      for (i = 0; i < atom->nlocal; i++) nspec += atom->nspecial[i][2];
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
	fprintf(screen,"Neighbor list builds = %d\n",neighbor->ncalls);
	fprintf(screen,"Dangerous builds = %d\n",neighbor->ndanger);
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
	fprintf(logfile,"Neighbor list builds = %d\n",neighbor->ncalls);
	fprintf(logfile,"Dangerous builds = %d\n",neighbor->ndanger);
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
