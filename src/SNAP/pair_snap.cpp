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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_snap.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "sna.h"
#include "memory.h"
#include "error.h"
#include "openmp_snap.h"
#include "domain.h"
#include <cmath>

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

PairSNAP::PairSNAP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  radelem = NULL;
  wjelem = NULL;
  coeffelem = NULL;

  nmax = 0;
  nthreads = 1;

  schedule_user = 0;
  schedule_time_guided = -1;
  schedule_time_dynamic = -1;
  ncalls_neigh =-1;

  ilistmask_max = 0;
  ilistmask = NULL;
  ghostinum = 0;
  ghostilist_max = 0;
  ghostilist = NULL;
  ghostnumneigh_max = 0;
  ghostnumneigh = NULL;
  ghostneighs = NULL;
  ghostfirstneigh = NULL;
  ghostneighs_total = 0;
  ghostneighs_max = 0;

  i_max = 0;
  i_neighmax = 0;
  i_numpairs = 0;
  i_rij = NULL;
  i_inside = NULL;
  i_wj = NULL;
  i_rcutij = NULL;
  i_ninside = NULL;
  i_pairs = NULL;
  i_uarraytot_r = NULL;
  i_uarraytot_i = NULL;
  i_zarray_r = NULL;
  i_zarray_i =NULL;

  use_shared_arrays = 0;
  timers[0] = 0;
  timers[1] = 0;
  timers[2] = 0;
  timers[3] = 0;

  // Need to set this because restart not handled by PairHybrid

  sna = NULL;

}

/* ---------------------------------------------------------------------- */

PairSNAP::~PairSNAP()
{
  if (nelements) {
    for (int i = 0; i < nelements; i++) 
      delete[] elements[i];
    delete[] elements;
    memory->destroy(radelem);
    memory->destroy(wjelem);
    memory->destroy(coeffelem);
  }

  // Need to set this because restart not handled by PairHybrid

  if (sna) {

    double time[5];
    double timeave[5];
    double timeave_mpi[5];
    double timemax_mpi[5];
    
    for (int i = 0; i < 5; i++) {
      time[i] = 0;
      timeave[i] = 0;
      for (int tid = 0; tid<nthreads; tid++) {
	if (sna[tid]->timers[i]>time[i])
	  time[i] = sna[tid]->timers[i];
	timeave[i] += sna[tid]->timers[i];
      }
      timeave[i] /= nthreads;
    }
    MPI_Reduce(timeave, timeave_mpi, 5, MPI_DOUBLE, MPI_SUM, 0, world);
    MPI_Reduce(time, timemax_mpi, 5, MPI_DOUBLE, MPI_MAX, 0, world);
    
    for (int tid = 0; tid<nthreads; tid++)
      delete sna[tid];
    delete [] sna;

  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
  }

}

void PairSNAP::compute(int eflag, int vflag)
{
  if (use_optimized)
    compute_optimized(eflag, vflag);
  else
    compute_regular(eflag, vflag);
}

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

void PairSNAP::compute_regular(int eflag, int vflag)
{
  int i,j,inum,jnum,ninside;
  double delx,dely,delz,evdwl,rsq;
  double fij[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  evdwl = 0.0;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  class SNA* snaptr = sna[0];

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    const int ielem = map[itype];
    const double radi = radelem[ielem];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // insure rij, inside, wj, and rcutij are of size jnum
      
    snaptr->grow_rij(jnum);

    // rij[][3] = displacements between atom I and those neighbors
    // inside = indices of neighbors of I within cutoff
    // wj = weights for neighbors of I within cutoff
    // rcutij = cutoffs for neighbors of I within cutoff
    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    ninside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type[j];
      int jelem = map[jtype];
      
      if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
	snaptr->rij[ninside][0] = delx;
	snaptr->rij[ninside][1] = dely;
	snaptr->rij[ninside][2] = delz;
	snaptr->inside[ninside] = j;
	snaptr->wj[ninside] = wjelem[jelem];
	snaptr->rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
	ninside++;
      }
    }

    // compute Ui, Zi, and Bi for atom I
    
    snaptr->compute_ui(ninside);
    snaptr->compute_zi();
    if (!gammaoneflag) {
      snaptr->compute_bi();
      snaptr->copy_bi2bvec();
    }

    // for neighbors of I within cutoff:
    // compute dUi/drj and dBi/drj
    // Fij = dEi/dRj = -dEi/dRi => add to Fi, subtract from Fj

    double* coeffi = coeffelem[ielem];

    for (int jj = 0; jj < ninside; jj++) {
      int j = snaptr->inside[jj];
      snaptr->compute_duidrj(snaptr->rij[jj],
			     snaptr->wj[jj],snaptr->rcutij[jj]);

      snaptr->compute_dbidrj();
      snaptr->copy_dbi2dbvec();

      fij[0] = 0.0;
      fij[1] = 0.0;
      fij[2] = 0.0;

      for (int k = 1; k <= ncoeff; k++) {
	double bgb;
	if (gammaoneflag) 
	  bgb = coeffi[k];
	else bgb = coeffi[k]*
	       gamma*pow(snaptr->bvec[k-1],gamma-1.0);
	fij[0] += bgb*snaptr->dbvec[k-1][0];
	fij[1] += bgb*snaptr->dbvec[k-1][1];
	fij[2] += bgb*snaptr->dbvec[k-1][2];
      }

      f[i][0] += fij[0];
      f[i][1] += fij[1];
      f[i][2] += fij[2];
      f[j][0] -= fij[0];
      f[j][1] -= fij[1];
      f[j][2] -= fij[2];

      if (evflag)
	ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
		     fij[0],fij[1],fij[2],
		     snaptr->rij[jj][0],snaptr->rij[jj][1],
		     snaptr->rij[jj][2]);
    }

    if (eflag) {
	
      // evdwl = energy of atom I, sum over coeffs_k * Bi_k

      evdwl = coeffi[0];
      if (gammaoneflag) {
	snaptr->compute_bi();
	snaptr->copy_bi2bvec();
	for (int k = 1; k <= ncoeff; k++)
	  evdwl += coeffi[k]*snaptr->bvec[k-1];
      } else
      	for (int k = 1; k <= ncoeff; k++)
      	  evdwl += coeffi[k]*pow(snaptr->bvec[k-1],gamma);
      ev_tally_full(i,2.0*evdwl,0.0,0.0,delx,dely,delz);
    }

  }
  
  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   This version is optimized for threading, micro-load balancing
   ---------------------------------------------------------------------- */

void PairSNAP::compute_optimized(int eflag, int vflag)
{
  // if reneighboring took place do load_balance if requested
  if (do_load_balance > 0 &&
      (neighbor->ncalls != ncalls_neigh)) {
    ghostinum = 0;
    // reset local ghost neighbor lists
    ncalls_neigh = neighbor->ncalls;
    if (ilistmask_max < list->inum) {
      memory->grow(ilistmask,list->inum,"PairSnap::ilistmask");
      ilistmask_max = list->inum;
    }
    for (int i = 0; i < list->inum; i++)
      ilistmask[i] = 1;

    //multiple passes for loadbalancing
    for (int i = 0; i < do_load_balance; i++)
      load_balance();
  }

  int numpairs = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    if ((do_load_balance <= 0) || ilistmask[ii]) {
      int i = list->ilist[ii];
      int jnum = list->numneigh[i];
      numpairs += jnum;
    }
  }

  if (do_load_balance)
    for (int ii = 0; ii < ghostinum; ii++) {
      int i = ghostilist[ii];
      int jnum = ghostnumneigh[i];
      numpairs += jnum;
    }

  // optimized schedule setting

  int time_dynamic = 0;
  int time_guided = 0;

  if (schedule_user == 0) schedule_user = 4;

  switch (schedule_user) {
  case 1:
    omp_set_schedule(omp_sched_static,1);
    break;
  case 2:
    omp_set_schedule(omp_sched_dynamic,1);
    break;
  case 3:
    omp_set_schedule(omp_sched_guided,2);
    break;
  case 4:
    omp_set_schedule(omp_sched_auto,0);
    break;
  case 5:
    if (numpairs < 8*nthreads) omp_set_schedule(omp_sched_dynamic,1);
    else if (schedule_time_guided < 0.0) {
      omp_set_schedule(omp_sched_guided,2);
      if (!eflag && !vflag) time_guided = 1;
    } else if (schedule_time_dynamic<0.0) {
      omp_set_schedule(omp_sched_dynamic,1);
      if (!eflag && !vflag) time_dynamic = 1;
    } else if (schedule_time_guided<schedule_time_dynamic)
      omp_set_schedule(omp_sched_guided,2);
    else
      omp_set_schedule(omp_sched_dynamic,1);
    break;
  }

  if (use_shared_arrays)
    build_per_atom_arrays();

#if defined(_OPENMP)
#pragma omp parallel shared(eflag,vflag,time_dynamic,time_guided) firstprivate(numpairs) default(none)
#endif
  {
    // begin of pragma omp parallel

    int tid = omp_get_thread_num();
    int** pairs_tid_unique = NULL;

    int** pairs;
    if (use_shared_arrays) pairs = i_pairs;
    else {
      memory->create(pairs_tid_unique,numpairs,4,"numpairs");
      pairs = pairs_tid_unique;
    }

    if (!use_shared_arrays) {
      numpairs = 0;
      for (int ii = 0; ii < list->inum; ii++) {
        if ((do_load_balance <= 0) || ilistmask[ii]) {
          int i = list->ilist[ii];
          int jnum = list->numneigh[i];
          for (int jj = 0; jj<jnum; jj++) {
            pairs[numpairs][0] = i;
            pairs[numpairs][1] = jj;
            pairs[numpairs][2] = -1;
            numpairs++;
          }
        }
      }

      for (int ii = 0; ii < ghostinum; ii++) {
        int i = ghostilist[ii];
        int jnum = ghostnumneigh[i];
        for (int jj = 0; jj<jnum; jj++) {
          pairs[numpairs][0] = i;
          pairs[numpairs][1] = jj;
          pairs[numpairs][2] = -1;
          numpairs++;
        }
      }
    }

    int ielem;
    int jj,k,inum,jnum,jtype,ninside;
    double delx,dely,delz,evdwl,rsq;
    double fij[3];
    int *ilist,*jlist,*numneigh,**firstneigh;
    evdwl = 0.0;

#if defined(_OPENMP)
#pragma omp master
#endif
    {
      if (eflag || vflag) ev_setup(eflag,vflag);
      else evflag = vflag_fdotr = 0;
    }

#if defined(_OPENMP)
#pragma omp barrier
    { ; }
#endif

    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int newton_pair = force->newton_pair;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // only update micro timers after setup
    static int count=0;
    if (count<2) {
      sna[tid]->timers[0] = 0;
      sna[tid]->timers[1] = 0;
      sna[tid]->timers[2] = 0;
      sna[tid]->timers[3] = 0;
      sna[tid]->timers[4] = 0;
    }
    count++;

    // did thread start working on interactions of new atom
    int iold = -1;

    double starttime, endtime;
    if (time_dynamic || time_guided)
      starttime = MPI_Wtime();

#if defined(_OPENMP)
#pragma omp for schedule(runtime)
#endif
    for (int iijj = 0; iijj < numpairs; iijj++) {
      int i = 0;
      if (use_shared_arrays) {
        i = i_pairs[iijj][0];
        if (iold != i) {
          set_sna_to_shared(tid,i_pairs[iijj][3]);
	  ielem = map[type[i]];
	}
        iold = i;
      } else {
        i = pairs[iijj][0];
        if (iold != i) {
          iold = i;
          const double xtmp = x[i][0];
          const double ytmp = x[i][1];
          const double ztmp = x[i][2];
          const int itype = type[i];
	  ielem = map[itype];
	  const double radi = radelem[ielem];

          if (i < nlocal) {
            jlist = firstneigh[i];
            jnum = numneigh[i];
          } else {
            jlist = ghostneighs+ghostfirstneigh[i];
            jnum = ghostnumneigh[i];
          }

          // insure rij, inside, wj, and rcutij are of size jnum

          sna[tid]->grow_rij(jnum);

          // rij[][3] = displacements between atom I and those neighbors
          // inside = indices of neighbors of I within cutoff
          // wj = weights of neighbors of I within cutoff
          // rcutij = cutoffs of neighbors of I within cutoff
          // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

          ninside = 0;
          for (jj = 0; jj < jnum; jj++) {
            int j = jlist[jj];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp; //unitialised
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;
            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];
	    int jelem = map[jtype];

            if (rsq < cutsq[itype][jtype]&&rsq>1e-20) { //unitialised
              sna[tid]->rij[ninside][0] = delx;
              sna[tid]->rij[ninside][1] = dely;
              sna[tid]->rij[ninside][2] = delz;
              sna[tid]->inside[ninside] = j;
              sna[tid]->wj[ninside] = wjelem[jelem];
              sna[tid]->rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
	      ninside++;

              // update index list with inside index
              pairs[iijj + (jj - pairs[iijj][1])][2] =
                ninside-1; //unitialised
            }
          }

          // compute Ui and Zi for atom I

          sna[tid]->compute_ui(ninside); //unitialised
          sna[tid]->compute_zi();
        }
      }

      // for neighbors of I within cutoff:
      // compute dUi/drj and dBi/drj
      // Fij = dEi/dRj = -dEi/dRi => add to Fi, subtract from Fj

      // entry into loop if inside index is set

      double* coeffi = coeffelem[ielem];

      if (pairs[iijj][2] >= 0) {
        jj = pairs[iijj][2];
        int j = sna[tid]->inside[jj];
        sna[tid]->compute_duidrj(sna[tid]->rij[jj],
				 sna[tid]->wj[jj],sna[tid]->rcutij[jj]);

        sna[tid]->compute_dbidrj();
        sna[tid]->copy_dbi2dbvec();
	if (!gammaoneflag) {
	  sna[tid]->compute_bi();
	  sna[tid]->copy_bi2bvec();
	}

        fij[0] = 0.0;
        fij[1] = 0.0;
        fij[2] = 0.0;

        for (k = 1; k <= ncoeff; k++) {
	  double bgb;
	  if (gammaoneflag) 
	    bgb = coeffi[k];
	  else bgb = coeffi[k]*
		 gamma*pow(sna[tid]->bvec[k-1],gamma-1.0);
	  fij[0] += bgb*sna[tid]->dbvec[k-1][0];
	  fij[1] += bgb*sna[tid]->dbvec[k-1][1];
	  fij[2] += bgb*sna[tid]->dbvec[k-1][2];
        }

#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          f[i][0] += fij[0];
          f[i][1] += fij[1];
          f[i][2] += fij[2];
          f[j][0] -= fij[0];
          f[j][1] -= fij[1];
          f[j][2] -= fij[2];
          if (evflag)
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                         fij[0],fij[1],fij[2],
                         sna[tid]->rij[jj][0],sna[tid]->rij[jj][1],
                         sna[tid]->rij[jj][2]);
        }
      }

      // evdwl = energy of atom I, sum over coeffs_k * Bi_k
      // only call this for first pair of each atom i
      // if atom has no pairs, eatom=0, which is wrong

      if (eflag&&pairs[iijj][1] == 0) {
	evdwl = coeffi[0];
	if (gammaoneflag) {
	  sna[tid]->compute_bi();
	  sna[tid]->copy_bi2bvec();
	  for (int k = 1; k <= ncoeff; k++)
	    evdwl += coeffi[k]*sna[tid]->bvec[k-1];
	} else
	  for (int k = 1; k <= ncoeff; k++)
	    evdwl += coeffi[k]*pow(sna[tid]->bvec[k-1],gamma);

#if defined(_OPENMP)
#pragma omp critical
#endif
        ev_tally_full(i,2.0*evdwl,0.0,0.0,delx,dely,delz);
      }

    }
    if (time_dynamic || time_guided)
      endtime = MPI_Wtime();
    if (time_dynamic) schedule_time_dynamic = endtime - starttime;
    if (time_guided) schedule_time_guided = endtime - starttime;
    if (!use_shared_arrays) memory->destroy(pairs);

  }// end of pragma omp parallel

  if (vflag_fdotr) virial_fdotr_compute();

}

inline int PairSNAP::equal(double* x,double* y)
{
  double dist2 =
    (x[0]-y[0])*(x[0]-y[0]) +
    (x[1]-y[1])*(x[1]-y[1]) +
    (x[2]-y[2])*(x[2]-y[2]);
  if (dist2 < 1e-20) return 1;
  return 0;
}

inline double PairSNAP::dist2(double* x,double* y)
{
  return
    (x[0]-y[0])*(x[0]-y[0]) +
    (x[1]-y[1])*(x[1]-y[1]) +
    (x[2]-y[2])*(x[2]-y[2]);
}

// return extra communication cutoff
// extra_cutoff = max(subdomain_length)

double PairSNAP::extra_cutoff()
{
  double sublo[3],subhi[3];

  if (domain->triclinic == 0) {
    for (int dim = 0 ; dim < 3 ; dim++) {
      sublo[dim] = domain->sublo[dim];
      subhi[dim] = domain->subhi[dim];
    }
  } else {
    domain->lamda2x(domain->sublo_lamda,sublo);
    domain->lamda2x(domain->subhi_lamda,subhi);
  }

  double sub_size[3];
  for (int dim = 0; dim < 3; dim++)
    sub_size[dim] = subhi[dim] - sublo[dim];

  double max_sub_size = 0;
  for (int dim = 0; dim < 3; dim++)
    max_sub_size = MAX(max_sub_size,sub_size[dim]);

  // note: for triclinic, probably need something different
  // see Comm::setup()

  return max_sub_size;
}

// micro load_balancer: each MPI process will
// check with each of its 26 neighbors,
// whether an imbalance exists in the number
// of atoms to calculate forces for.
// If it does it will set ilistmask of one of
// its local atoms to zero, and send its Tag
// to the neighbor process. The neighboring process
// will check its ghost list for the
// ghost atom with the same Tag which is closest
// to its domain center, and build a
// neighborlist for this ghost atom. For this to work,
// the communication cutoff has to be
// as large as the neighbor cutoff +
// maximum subdomain length.

// Note that at most one atom is exchanged per processor pair.

// Also note that the local atom assignement
// doesn't change. This load balancer will cause
// some ghost atoms to have full neighborlists
// which are unique to PairSNAP.
// They are not part of the generally accessible neighborlist.
// At the same time corresponding local atoms on
// other MPI processes will not be
// included in the force computation since
// their ilistmask is 0. This does not effect
// any other classes which might
// access the same general neighborlist.
// Reverse communication (newton on) of forces is required.

// Currently the load balancer does two passes,
// since its exchanging atoms upstream and downstream.

void PairSNAP::load_balance()
{
  double sublo[3],subhi[3];
  if (domain->triclinic == 0) {
    double* sublotmp = domain->sublo;
    double* subhitmp = domain->subhi;
    for (int dim = 0 ; dim<3 ; dim++) {
      sublo[dim]=sublotmp[dim];
      subhi[dim]=subhitmp[dim];
    }
  } else {
    double* sublotmp = domain->sublo_lamda;
    double* subhitmp = domain->subhi_lamda;
    domain->lamda2x(sublotmp,sublo);
    domain->lamda2x(subhitmp,subhi);
  }

  //if (list->inum==0) list->grow(atom->nmax);

  int nlocal = ghostinum;
  for (int i=0; i < list->inum; i++)
    if (ilistmask[i]) nlocal++;
  int ***grid2proc = comm->grid2proc;
  int* procgrid = comm->procgrid;

  int nlocal_up,nlocal_down;
  MPI_Status status;
  MPI_Request request;

  double sub_mid[3];
  for (int dim=0; dim<3; dim++)
    sub_mid[dim] = (subhi[dim] + sublo[dim])/2;

  if (comm->cutghostuser <
      neighbor->cutneighmax+extra_cutoff())
    error->all(FLERR,"Communication cutoff is too small "
               "for SNAP micro load balancing.\n"
               "Typically this can happen, if you change "
               "the neighbor skin after your pair_style "
               "command or if your box dimensions grow "
               "during the run.\n"
               "You need to set it via "
               "'communicate single cutoff NUMBER' "
               "to the needed length.");

  int nrecv = ghostinum;
  int totalsend = 0;
  int nsend = 0;
  int depth = 1;

  for (int dx = -depth; dx < depth+1; dx++)
    for (int dy = -depth; dy < depth+1; dy++)
      for (int dz = -depth; dz < depth+1; dz++) {

        if (dx == dy && dy == dz && dz == 0) continue;

        int sendloc[3] = {comm->myloc[0],
                          comm->myloc[1], comm->myloc[2]
                         };
        sendloc[0] += dx;
        sendloc[1] += dy;
        sendloc[2] += dz;
        for (int dim = 0; dim < 3; dim++)
          if (sendloc[dim] >= procgrid[dim])
            sendloc[dim] = sendloc[dim] - procgrid[dim];
        for (int dim = 0; dim < 3; dim++)
          if (sendloc[dim] < 0)
            sendloc[dim] = procgrid[dim] + sendloc[dim];
        int recvloc[3] = {comm->myloc[0],
                          comm->myloc[1], comm->myloc[2]
                         };
        recvloc[0] -= dx;
        recvloc[1] -= dy;
        recvloc[2] -= dz;
        for (int dim = 0; dim < 3; dim++)
          if (recvloc[dim] < 0)
            recvloc[dim] = procgrid[dim] + recvloc[dim];
        for (int dim = 0; dim < 3; dim++)
          if (recvloc[dim] >= procgrid[dim])
            recvloc[dim] = recvloc[dim] - procgrid[dim];

        int sendproc = grid2proc[sendloc[0]][sendloc[1]][sendloc[2]];
        int recvproc = grid2proc[recvloc[0]][recvloc[1]][recvloc[2]];

        // two stage process, first upstream movement, then downstream

        MPI_Sendrecv(&nlocal,1,MPI_INT,sendproc,0,
                     &nlocal_up,1,MPI_INT,recvproc,0,world,&status);
        MPI_Sendrecv(&nlocal,1,MPI_INT,recvproc,0,
                     &nlocal_down,1,MPI_INT,sendproc,0,world,&status);
        nsend = 0;

        // send upstream

        if (nlocal > nlocal_up+1) {

          int i = totalsend++;
          while(i < list->inum && ilistmask[i] == 0)
            i = totalsend++;

          if (i < list->inum)
            MPI_Isend(&atom->tag[i],1,MPI_INT,recvproc,0,world,&request);
          else {
            int j = -1;
            MPI_Isend(&j,1,MPI_INT,recvproc,0,world,&request);
          }

          if (i < list->inum) {
            for (int j = 0; j < list->inum; j++)
              if (list->ilist[j] == i)
                ilistmask[j] = 0;
            nsend = 1;
          }
        }

        // recv downstream

        if (nlocal < nlocal_down-1) {
          nlocal++;
          int get_tag = -1;
          MPI_Recv(&get_tag,1,MPI_INT,sendproc,0,world,&status);

          // if get_tag -1 the other process didnt have local atoms to send

          if (get_tag >= 0) {
            if (ghostinum >= ghostilist_max) {
              memory->grow(ghostilist,ghostinum+10,
                           "PairSnap::ghostilist");
              ghostilist_max = ghostinum+10;
            }
            if (atom->nlocal + atom->nghost >= ghostnumneigh_max) {
              ghostnumneigh_max = atom->nlocal+atom->nghost+100;
              memory->grow(ghostnumneigh,ghostnumneigh_max,
                           "PairSnap::ghostnumneigh");
              memory->grow(ghostfirstneigh,ghostnumneigh_max,
                           "PairSnap::ghostfirstneigh");
            }

            // find closest ghost image of the transfered particle

            double mindist = 1e200;
            int closestghost = -1;
            for (int j = 0; j < atom->nlocal + atom->nghost; j++)
              if (atom->tag[j] == get_tag)
                if (dist2(sub_mid, atom->x[j]) < mindist) {
                  closestghost = j;
                  mindist = dist2(sub_mid, atom->x[j]);
                }

            // build neighborlist for this particular
            // ghost atom, and add it to list->ilist

            if (ghostneighs_max - ghostneighs_total <
                neighbor->oneatom) {
              memory->grow(ghostneighs,
                           ghostneighs_total + neighbor->oneatom,
                           "PairSnap::ghostneighs");
              ghostneighs_max = ghostneighs_total + neighbor->oneatom;
            }

            int j = closestghost;

            ghostilist[ghostinum] = j;
            ghostnumneigh[j] = 0;
            ghostfirstneigh[j] = ghostneighs_total;

            ghostinum++;
            int* jlist = ghostneighs + ghostfirstneigh[j];

            // find all neighbors by looping
            // over all local and ghost atoms

            for (int k = 0; k < atom->nlocal + atom->nghost; k++)
              if (dist2(atom->x[j],atom->x[k]) <
                  neighbor->cutneighmax*neighbor->cutneighmax) {
                jlist[ghostnumneigh[j]] = k;
                ghostnumneigh[j]++;
                ghostneighs_total++;
              }
          }

          if (get_tag >= 0) nrecv++;
        }

        // decrease nlocal later, so that it is the
        // initial number both for receiving and sending

        if (nsend) nlocal--;

        // second pass through the grid

        MPI_Sendrecv(&nlocal,1,MPI_INT,sendproc,0,
                     &nlocal_up,1,MPI_INT,recvproc,0,world,&status);
        MPI_Sendrecv(&nlocal,1,MPI_INT,recvproc,0,
                     &nlocal_down,1,MPI_INT,sendproc,0,world,&status);

        // send downstream

        nsend=0;
        if (nlocal > nlocal_down+1) {
          int i = totalsend++;
          while(i < list->inum && ilistmask[i]==0) i = totalsend++;

          if (i < list->inum)
            MPI_Isend(&atom->tag[i],1,MPI_INT,sendproc,0,world,&request);
          else {
            int j =- 1;
            MPI_Isend(&j,1,MPI_INT,sendproc,0,world,&request);
          }

          if (i < list->inum) {
            for (int j=0; j<list->inum; j++)
              if (list->ilist[j] == i) ilistmask[j] = 0;
            nsend = 1;
          }
        }

        // receive upstream

        if (nlocal < nlocal_up-1) {
          nlocal++;
          int get_tag = -1;

          MPI_Recv(&get_tag,1,MPI_INT,recvproc,0,world,&status);

          if (get_tag >= 0) {
            if (ghostinum >= ghostilist_max) {
              memory->grow(ghostilist,ghostinum+10,
                           "PairSnap::ghostilist");
              ghostilist_max = ghostinum+10;
            }
            if (atom->nlocal + atom->nghost >= ghostnumneigh_max) {
              ghostnumneigh_max = atom->nlocal + atom->nghost + 100;
              memory->grow(ghostnumneigh,ghostnumneigh_max,
                           "PairSnap::ghostnumneigh");
              memory->grow(ghostfirstneigh,ghostnumneigh_max,
                           "PairSnap::ghostfirstneigh");
            }

            // find closest ghost image of the transfered particle

            double mindist = 1e200;
            int closestghost = -1;
            for (int j = 0; j < atom->nlocal + atom->nghost; j++)
              if (atom->tag[j] == get_tag)
                if (dist2(sub_mid,atom->x[j])<mindist) {
                  closestghost = j;
                  mindist = dist2(sub_mid,atom->x[j]);
                }

            // build neighborlist for this particular ghost atom

            if (ghostneighs_max-ghostneighs_total < neighbor->oneatom) {
              memory->grow(ghostneighs,ghostneighs_total + neighbor->oneatom,
                           "PairSnap::ghostneighs");
              ghostneighs_max = ghostneighs_total + neighbor->oneatom;
            }

            int j = closestghost;

            ghostilist[ghostinum] = j;
            ghostnumneigh[j] = 0;
            ghostfirstneigh[j] = ghostneighs_total;

            ghostinum++;
            int* jlist = ghostneighs + ghostfirstneigh[j];

            for (int k = 0; k < atom->nlocal + atom->nghost; k++)
              if (dist2(atom->x[j],atom->x[k]) <
                  neighbor->cutneighmax*neighbor->cutneighmax) {
                jlist[ghostnumneigh[j]] = k;
                ghostnumneigh[j]++;
                ghostneighs_total++;
              }
          }

          if (get_tag >= 0) nrecv++;
        }
        if (nsend) nlocal--;
      }
}

void PairSNAP::set_sna_to_shared(int snaid,int i)
{
  sna[snaid]->rij = i_rij[i];
  sna[snaid]->inside = i_inside[i];
  sna[snaid]->wj = i_wj[i];
  sna[snaid]->rcutij = i_rcutij[i];
  sna[snaid]->zarray_r = i_zarray_r[i];
  sna[snaid]->zarray_i = i_zarray_i[i];
  sna[snaid]->uarraytot_r = i_uarraytot_r[i];
  sna[snaid]->uarraytot_i = i_uarraytot_i[i];
}

void PairSNAP::build_per_atom_arrays()
{

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&starttime);
#endif

  int i_list[list->inum+ghostinum];
  int count = 0;
  int neighmax = 0;
  for (int ii = 0; ii < list->inum; ii++)
    if ((do_load_balance <= 0) || ilistmask[ii]) {
      neighmax=MAX(neighmax,list->numneigh[list->ilist[ii]]);
      i_list[count++]=list->ilist[ii];
    }
  for (int ii = 0; ii < ghostinum; ii++) {
    neighmax=MAX(neighmax,ghostnumneigh[ghostilist[ii]]);
    i_list[count++]=ghostilist[ii];
  }

  if (i_max < count || i_neighmax < neighmax) {
    int i_maxt = MAX(count,i_max);
    i_neighmax = MAX(neighmax,i_neighmax);
    memory->destroy(i_rij);
    memory->destroy(i_inside);
    memory->destroy(i_wj);
    memory->destroy(i_rcutij);
    memory->destroy(i_ninside);
    memory->destroy(i_pairs);
    memory->create(i_rij,i_maxt,i_neighmax,3,"PairSNAP::i_rij");
    memory->create(i_inside,i_maxt,i_neighmax,"PairSNAP::i_inside");
    memory->create(i_wj,i_maxt,i_neighmax,"PairSNAP::i_wj");
    memory->create(i_rcutij,i_maxt,i_neighmax,"PairSNAP::i_rcutij");
    memory->create(i_ninside,i_maxt,"PairSNAP::i_ninside");
    memory->create(i_pairs,i_maxt*i_neighmax,4,"PairSNAP::i_pairs");
  }

  if (i_max < count) {
    int jdim = sna[0]->twojmax+1;
    memory->destroy(i_uarraytot_r);
    memory->destroy(i_uarraytot_i);
    memory->create(i_uarraytot_r,count,jdim,jdim,jdim,
                   "PairSNAP::i_uarraytot_r");
    memory->create(i_uarraytot_i,count,jdim,jdim,jdim,
                   "PairSNAP::i_uarraytot_i");
    if (i_zarray_r != NULL)
      for (int i = 0; i < i_max; i++) {
        memory->destroy(i_zarray_r[i]);
        memory->destroy(i_zarray_i[i]);
      }

    delete [] i_zarray_r;
    delete [] i_zarray_i;
    i_zarray_r = new double*****[count];
    i_zarray_i = new double*****[count];
    for (int i = 0; i < count; i++) {
      memory->create(i_zarray_r[i],jdim,jdim,jdim,jdim,jdim,
                     "PairSNAP::i_zarray_r");
      memory->create(i_zarray_i[i],jdim,jdim,jdim,jdim,jdim,
                     "PairSNAP::i_zarray_i");
    }
  }

  if (i_max < count)
    i_max = count;

  count = 0;
  i_numpairs = 0;
  for (int ii = 0; ii < list->inum; ii++) {
    if ((do_load_balance <= 0) || ilistmask[ii]) {
      int i = list->ilist[ii];
      int jnum = list->numneigh[i];
      int* jlist = list->firstneigh[i];
      const double xtmp = atom->x[i][0];
      const double ytmp = atom->x[i][1];
      const double ztmp = atom->x[i][2];
      const int itype = atom->type[i];
      const int ielem = map[itype];
      const double radi = radelem[ielem];
      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        const double delx = atom->x[j][0] - xtmp;
        const double dely = atom->x[j][1] - ytmp;
        const double delz = atom->x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
        int jtype = atom->type[j];
	int jelem = map[jtype];

        i_pairs[i_numpairs][0] = i;
        i_pairs[i_numpairs][1] = jj;
        i_pairs[i_numpairs][2] = -1;
        i_pairs[i_numpairs][3] = count;
        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          i_rij[count][ninside][0] = delx;
          i_rij[count][ninside][1] = dely;
          i_rij[count][ninside][2] = delz;
          i_inside[count][ninside] = j;
          i_wj[count][ninside] = wjelem[jelem];
          i_rcutij[count][ninside] = (radi + radelem[jelem])*rcutfac;

          // update index list with inside index
          i_pairs[i_numpairs][2] = ninside++;
        }
        i_numpairs++;
      }
      i_ninside[count] = ninside;
      count++;
    }
  }

  for (int ii = 0; ii < ghostinum; ii++) {
    int i = ghostilist[ii];
    int jnum = ghostnumneigh[i];
    int* jlist = ghostneighs+ghostfirstneigh[i];
    const double xtmp = atom->x[i][0];
    const double ytmp = atom->x[i][1];
    const double ztmp = atom->x[i][2];
    const int itype = atom->type[i];
    const int ielem = map[itype];
    const double radi = radelem[ielem];
    int ninside = 0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const double delx = atom->x[j][0] - xtmp;
      const double dely = atom->x[j][1] - ytmp;
      const double delz = atom->x[j][2] - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = atom->type[j];
      int jelem = map[jtype];

      i_pairs[i_numpairs][0] = i;
      i_pairs[i_numpairs][1] = jj;
      i_pairs[i_numpairs][2] = -1;
      i_pairs[i_numpairs][3] = count;
      if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
        i_rij[count][ninside][0] = delx;
        i_rij[count][ninside][1] = dely;
        i_rij[count][ninside][2] = delz;
        i_inside[count][ninside] = j;
        i_wj[count][ninside] = wjelem[jelem];
        i_rcutij[count][ninside] = (radi + radelem[jelem])*rcutfac;
        // update index list with inside index
        i_pairs[i_numpairs][2] = ninside++;
      }
      i_numpairs++;
    }
    i_ninside[count] = ninside;
    count++;
  }
#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&endtime);
  timers[0]+=(endtime.tv_sec-starttime.tv_sec+1.0*
              (endtime.tv_nsec-starttime.tv_nsec)/1000000000);
#endif
#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&starttime);
#endif

#if defined(_OPENMP)
#pragma omp parallel for shared(count) default(none)
#endif
  for (int ii=0; ii < count; ii++) {
    int tid = omp_get_thread_num();
    set_sna_to_shared(tid,ii);
    //sna[tid]->compute_ui(i_ninside[ii]);
#ifdef TIMING_INFO
    clock_gettime(CLOCK_REALTIME,&starttime);
#endif
    sna[tid]->compute_ui_omp(i_ninside[ii],MAX(int(nthreads/count),1));
#ifdef TIMING_INFO
    clock_gettime(CLOCK_REALTIME,&endtime);
    sna[tid]->timers[0]+=(endtime.tv_sec-starttime.tv_sec+1.0*
                          (endtime.tv_nsec-starttime.tv_nsec)/1000000000);
#endif
  }

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&starttime);
#endif
  for (int ii=0; ii < count; ii++) {
    int tid = 0;//omp_get_thread_num();
    set_sna_to_shared(tid,ii);
    sna[tid]->compute_zi_omp(MAX(int(nthreads/count),1));
  }
#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&endtime);
  sna[0]->timers[1]+=(endtime.tv_sec-starttime.tv_sec+1.0*
                      (endtime.tv_nsec-starttime.tv_nsec)/1000000000);
#endif

#ifdef TIMING_INFO
  clock_gettime(CLOCK_REALTIME,&endtime);
  timers[1]+=(endtime.tv_sec-starttime.tv_sec+1.0*
              (endtime.tv_nsec-starttime.tv_nsec)/1000000000);
#endif
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSNAP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(map,n+1,"pair:map");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSNAP::settings(int narg, char **arg)
{

  // set default values for optional arguments

  nthreads = -1;
  use_shared_arrays=-1;
  do_load_balance = 0;
  use_optimized = 1;

  // optional arguments

  for (int i=0; i < narg; i++) {
    if (i+2>narg) error->all(FLERR,"Illegal pair_style command."
			     " Too few arguments.");
    if (strcmp(arg[i],"nthreads")==0) {
      nthreads=force->inumeric(FLERR,arg[++i]);
#if defined(LMP_USER_OMP)
      error->all(FLERR,"Please set number of threads via package omp command");
#else
      omp_set_num_threads(nthreads);
      comm->nthreads=nthreads;
#endif
      continue;
    }
    if (strcmp(arg[i],"optimized")==0) {
      use_optimized=force->inumeric(FLERR,arg[++i]);
      continue;
    }
    if (strcmp(arg[i],"shared")==0) {
      use_shared_arrays=force->inumeric(FLERR,arg[++i]);
      continue;
    }
    if (strcmp(arg[i],"loadbalance")==0) {
      do_load_balance = force->inumeric(FLERR,arg[++i]);
      if (do_load_balance) {
	double mincutoff = extra_cutoff() +
	  rcutmax + neighbor->skin;
	if (comm->cutghostuser < mincutoff) {
	  char buffer[255];

	  //apparently mincutoff is 0 after sprintf command ?????

	  double tmp = mincutoff + 0.1;
	  sprintf(buffer, "Communication cutoff is too small "
		  "for SNAP micro load balancing. "
		  "It will be increased to: %lf.",mincutoff+0.1);
	  if (comm->me==0)
	    error->warning(FLERR,buffer);

	  comm->cutghostuser = tmp;

	}
      }
      continue;
    }
    if (strcmp(arg[i],"schedule")==0) {
      i++;
      if (strcmp(arg[i],"static")==0)
	schedule_user = 1;
      if (strcmp(arg[i],"dynamic")==0)
	schedule_user = 2;
      if (strcmp(arg[i],"guided")==0)
	schedule_user = 3;
      if (strcmp(arg[i],"auto")==0)
	schedule_user = 4;
      if (strcmp(arg[i],"determine")==0)
	schedule_user = 5;
      if (schedule_user == 0)
	error->all(FLERR,"Illegal pair_style command."
		   " Illegal schedule argument.");
      continue;
    }
    char buffer[255];
    sprintf(buffer, "Illegal pair_style command."
	    " Unrecognized argument: %s.\n",arg[i]);
    error->all(FLERR,buffer);
  }

  if (nthreads < 0)
    nthreads = comm->nthreads;

  if (use_shared_arrays < 0) {
    if (nthreads > 1 && atom->nlocal <= 2*nthreads)
      use_shared_arrays = 1;
    else use_shared_arrays = 0;
  }

  // check if running non-optimized code with
  // optimization flags set

  if (!use_optimized)
    if (nthreads > 1 || 
	use_shared_arrays || 
	do_load_balance ||
	schedule_user)
      error->all(FLERR,"Illegal pair_style command."
                 "Advanced options require setting 'optimized 1'.");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSNAP::coeff(int narg, char **arg)
{
  // read SNAP element names between 2 filenames
  // nelements = # of SNAP elements
  // elements = list of unique element names

  if (narg < 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  if (nelements) {
    for (int i = 0; i < nelements; i++) 
      delete[] elements[i];
    delete[] elements;
    memory->destroy(radelem);
    memory->destroy(wjelem);
    memory->destroy(coeffelem);
  }

  nelements = narg - 4 - atom->ntypes;
  if (nelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");

  char* type1 = arg[0];
  char* type2 = arg[1];
  char* coefffilename = arg[2];
  char** elemlist = &arg[3];
  char* paramfilename = arg[3+nelements];
  char** elemtypes = &arg[4+nelements];

  // insure I,J args are * *

  if (strcmp(type1,"*") != 0 || strcmp(type2,"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  elements = new char*[nelements];

  for (int i = 0; i < nelements; i++) {
    char* elemname = elemlist[i];
    int n = strlen(elemname) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],elemname);
  }

  // read snapcoeff and snapparam files

  read_files(coefffilename,paramfilename);

  // read args that map atom types to SNAP elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  for (int i = 1; i <= atom->ntypes; i++) {
    char* elemname = elemtypes[i-1];
    int jelem;
    for (jelem = 0; jelem < nelements; jelem++)
      if (strcmp(elemname,elements[jelem]) == 0)
	break;

    if (jelem < nelements)
      map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) map[i] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  sna = new SNA*[nthreads];

  // allocate memory for per OpenMP thread data which
  // is wrapped into the sna class

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int tid = omp_get_thread_num();
    sna[tid] = new SNA(lmp,rfac0,twojmax,
                       diagonalstyle,use_shared_arrays,
		       rmin0,switchflag);
    if (!use_shared_arrays)
      sna[tid]->grow_rij(nmax);
  }

  if (ncoeff != sna[0]->ncoeff) {
    printf("ncoeff = %d snancoeff = %d \n",ncoeff,sna[0]->ncoeff);
    error->all(FLERR,"Incorrect SNAP parameter file");
  }

  // Calculate maximum cutoff for all elements

  rcutmax = 0.0;
  for (int ielem = 0; ielem < nelements; ielem++)
    rcutmax = MAX(2.0*radelem[ielem]*rcutfac,rcutmax);

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSNAP::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int tid = omp_get_thread_num();
    sna[tid]->init();
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSNAP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return (radelem[map[i]] + 
  	  radelem[map[j]])*rcutfac;
}

/* ---------------------------------------------------------------------- */

void PairSNAP::read_files(char *coefffilename, char *paramfilename)
{

  // open SNAP ceofficient file on proc 0

  FILE *fpcoeff;
  if (comm->me == 0) {
    fpcoeff = force->open_potential(coefffilename);
    if (fpcoeff == NULL) {
      char str[128];
      sprintf(str,"Cannot open SNAP coefficient file %s",coefffilename);
      error->one(FLERR,str);
    }
  }

  char line[MAXLINE],*ptr;
  int eof = 0;

  int n;
  int nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpcoeff);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
  }
  if (nwords != 2) 
    error->all(FLERR,"Incorrect format in SNAP coefficient file");

  // words = ptrs to all words in line
  // strip single and double quotes from words

  char* words[MAXWORD];
  int iword = 0;
  words[iword] = strtok(line,"' \t\n\r\f");
  iword = 1;
  words[iword] = strtok(NULL,"' \t\n\r\f");

  int nelemfile = atoi(words[0]);
  ncoeff = atoi(words[1])-1;

  // Set up element lists 

  memory->create(radelem,nelements,"pair:radelem");
  memory->create(wjelem,nelements,"pair:wjelem");
  memory->create(coeffelem,nelements,ncoeff+1,"pair:coeffelem");

  int *found = new int[nelements];
  for (int ielem = 0; ielem < nelements; ielem++) 
    found[ielem] = 0;

  // Loop over elements in the SNAP coefficient file

  for (int ielemfile = 0; ielemfile < nelemfile; ielemfile++) {
 
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == NULL) {
	eof = 1;
	fclose(fpcoeff);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) 
      error->all(FLERR,"Incorrect format in SNAP coefficient file");
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    nwords = atom->count_words(line);
    if (nwords != 3) 
      error->all(FLERR,"Incorrect format in SNAP coefficient file");

    iword = 0;
    words[iword] = strtok(line,"' \t\n\r\f");
    iword = 1;
    words[iword] = strtok(NULL,"' \t\n\r\f");
    iword = 2;
    words[iword] = strtok(NULL,"' \t\n\r\f");

    char* elemtmp = words[0];
    double radtmp = atof(words[1]);
    double wjtmp = atof(words[2]);

    // skip if element name isn't in element list

    int ielem;
    for (ielem = 0; ielem < nelements; ielem++)
      if (strcmp(elemtmp,elements[ielem]) == 0) break;
    if (ielem == nelements) {
      if (comm->me == 0)
	for (int icoeff = 0; icoeff <= ncoeff; icoeff++) 
	  ptr = fgets(line,MAXLINE,fpcoeff);
      continue;
    }

    // skip if element already appeared

    if (found[ielem]) {
      if (comm->me == 0)
	for (int icoeff = 0; icoeff <= ncoeff; icoeff++) 
	  ptr = fgets(line,MAXLINE,fpcoeff);
      continue;
    }

    found[ielem] = 1;
    radelem[ielem] = radtmp;
    wjelem[ielem] = wjtmp;


    if (comm->me == 0) {
      if (screen) fprintf(screen,"SNAP Element = %s, Radius %g, Weight %g \n", 
			  elements[ielem], radelem[ielem], wjelem[ielem]);
      if (logfile) fprintf(logfile,"SNAP Element = %s, Radius %g, Weight %g \n", 
			  elements[ielem], radelem[ielem], wjelem[ielem]);
    }

    for (int icoeff = 0; icoeff <= ncoeff; icoeff++) { 
      if (comm->me == 0) {
	ptr = fgets(line,MAXLINE,fpcoeff);
	if (ptr == NULL) {
	  eof = 1;
	  fclose(fpcoeff);
	} else n = strlen(line) + 1;
      }

      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) 
	error->all(FLERR,"Incorrect format in SNAP coefficient file");
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);

      nwords = atom->count_words(line);
      if (nwords != 1) 
	error->all(FLERR,"Incorrect format in SNAP coefficient file");

      iword = 0;
      words[iword] = strtok(line,"' \t\n\r\f");

      coeffelem[ielem][icoeff] = atof(words[0]);
      
    }
  }

  // set flags for required keywords

  rcutfacflag = 0;
  twojmaxflag = 0;

  // Set defaults for optional keywords

  gamma = 1.0;
  gammaoneflag = 1;
  rfac0 = 0.99363;
  rmin0 = 0.0;
  diagonalstyle = 3;
  switchflag = 1;
  // open SNAP parameter file on proc 0

  FILE *fpparam;
  if (comm->me == 0) {
    fpparam = force->open_potential(paramfilename);
    if (fpparam == NULL) {
      char str[128];
      sprintf(str,"Cannot open SNAP parameter file %s",paramfilename);
      error->one(FLERR,str);
    }
  }

  eof = 0;
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpparam);
      if (ptr == NULL) {
        eof = 1;
        fclose(fpparam);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    if (nwords != 2) 
      error->all(FLERR,"Incorrect format in SNAP parameter file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    char* keywd = strtok(line,"' \t\n\r\f");
    char* keyval = strtok(NULL,"' \t\n\r\f");

    if (comm->me == 0) {
      if (screen) fprintf(screen,"SNAP keyword %s %s \n",keywd,keyval);
      if (logfile) fprintf(logfile,"SNAP keyword %s %s \n",keywd,keyval);
    }

    if (strcmp(keywd,"rcutfac") == 0) {
      rcutfac = atof(keyval);
      rcutfacflag = 1;
    } else if (strcmp(keywd,"twojmax") == 0) {
      twojmax = atoi(keyval);
      twojmaxflag = 1;
    } else if (strcmp(keywd,"gamma") == 0)
      gamma = atof(keyval);
    else if (strcmp(keywd,"rfac0") == 0)
      rfac0 = atof(keyval);
    else if (strcmp(keywd,"rmin0") == 0)
      rmin0 = atof(keyval);
    else if (strcmp(keywd,"diagonalstyle") == 0)
      diagonalstyle = atoi(keyval);
    else if (strcmp(keywd,"switchflag") == 0)
      switchflag = atoi(keyval);
    else
      error->all(FLERR,"Incorrect SNAP parameter file");
  }

  if (rcutfacflag == 0 || twojmaxflag == 0)
    error->all(FLERR,"Incorrect SNAP parameter file");

  if (gamma == 1.0) gammaoneflag = 1;
  else gammaoneflag = 0;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairSNAP::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes+1;
  bytes += n*n*sizeof(int);
  bytes += n*n*sizeof(double);
  bytes += 3*nmax*sizeof(double);
  bytes += nmax*sizeof(int);
  bytes += (2*ncoeff+1)*sizeof(double);
  bytes += (ncoeff*3)*sizeof(double);
  bytes += sna[0]->memory_usage()*nthreads;
  return bytes;
}

