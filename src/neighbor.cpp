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

/* ----------------------------------------------------------------------
   Contributing author (triclinic and multi-neigh) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "update.h"
#include "respa.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define RQDELTA 1
#define EXDELTA 1

#define LB_FACTOR 1.5
#define SMALL 1.0e-6
#define BIG 1.0e20
#define CUT2BIN_RATIO 100

enum{NSQ,BIN,MULTI};     // also in neigh_list.cpp

//#define NEIGH_LIST_DEBUG 1

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  style = BIN;
  every = 1;
  delay = 10;
  dist_check = 1;
  pgsize = 100000;
  oneatom = 2000;
  binsizeflag = 0;
  build_once = 0;

  cutneighsq = NULL;
  cutneighghostsq = NULL;
  cuttype = NULL;
  cuttypesq = NULL;
  fixchecklist = NULL;

  // coords at last neighboring

  maxhold = 0;
  xhold = NULL;

  // binning

  maxhead = 0;
  binhead = NULL;
  maxbin = 0;
  bins = NULL;

  // pair exclusion list info

  includegroup = 0;

  nex_type = maxex_type = 0;
  ex1_type = ex2_type = NULL;
  ex_type = NULL;

  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = NULL;

  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = NULL;

  // pair lists

  maxatom = 0;
  nblist = nglist = nslist = 0;

  nlist = 0;
  lists = NULL;
  pair_build = NULL;
  stencil_create = NULL;
  blist = glist = slist = NULL;
  anyghostlist = 0;

  nrequest = maxrequest = 0;
  requests = NULL;

  old_style = BIN;
  old_triclinic = 0;
  old_nrequest = 0;
  old_requests = NULL;

  // bond lists

  maxbond = 0;
  bondlist = NULL;
  maxangle = 0;
  anglelist = NULL;
  maxdihedral = 0;
  dihedrallist = NULL;
  maximproper = 0;
  improperlist = NULL;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{
  memory->destroy(cutneighsq);
  memory->destroy(cutneighghostsq);
  delete [] cuttype;
  delete [] cuttypesq;
  delete [] fixchecklist;

  memory->destroy(xhold);

  memory->destroy(binhead);
  memory->destroy(bins);

  memory->destroy(ex1_type);
  memory->destroy(ex2_type);
  memory->destroy(ex_type);

  memory->destroy(ex1_group);
  memory->destroy(ex2_group);
  delete [] ex1_bit;
  delete [] ex2_bit;

  memory->destroy(ex_mol_group);
  delete [] ex_mol_bit;

  for (int i = 0; i < nlist; i++) delete lists[i];
  delete [] lists;
  delete [] pair_build;
  delete [] stencil_create;
  delete [] blist;
  delete [] glist;
  delete [] slist;

  for (int i = 0; i < nrequest; i++) delete requests[i];
  memory->sfree(requests);
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);

  memory->destroy(bondlist);
  memory->destroy(anglelist);
  memory->destroy(dihedrallist);
  memory->destroy(improperlist);
}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,m,n;

  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
  newton_pair = force->newton_pair;

  // error check

  if (delay > 0 && (delay % every) != 0)
    error->all(FLERR,"Neighbor delay must be 0 or multiple of every setting");

  if (pgsize < 10*oneatom)
    error->all(FLERR,"Neighbor page size must be >= 10x the one atom setting");

  // ------------------------------------------------------------------
  // settings

  // bbox lo/hi = bounding box of entire domain, stored by Domain

  if (triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  // set neighbor cutoffs (force cutoff + skin)
  // trigger determines when atoms migrate and neighbor lists are rebuilt
  //   needs to be non-zero for migration distance check
  //   even if pair = NULL and no neighbor lists are used
  // cutneigh = force cutoff + skin if cutforce > 0, else cutneigh = 0

  triggersq = 0.25*skin*skin;
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic || 
			     (dimension == 3 && domain->zperiodic)))
      boxcheck = 1;
      
  n = atom->ntypes;
  if (cutneighsq == NULL) {
    memory->create(cutneighsq,n+1,n+1,"neigh:cutneighsq");
    memory->create(cutneighghostsq,n+1,n+1,"neigh:cutneighghostsq");
    cuttype = new double[n+1];
    cuttypesq = new double[n+1];
  }

  double cutoff,delta,cut;
  cutneighmin = BIG;
  cutneighmax = 0.0;

  for (i = 1; i <= n; i++) {
    cuttype[i] = cuttypesq[i] = 0.0;
    for (j = 1; j <= n; j++) {
      if (force->pair) cutoff = sqrt(force->pair->cutsq[i][j]);
      else cutoff = 0.0;
      if (cutoff > 0.0) delta = skin;
      else delta = 0.0;
      cut = cutoff + delta;

      cutneighsq[i][j] = cut*cut;
      cuttype[i] = MAX(cuttype[i],cut);
      cuttypesq[i] = MAX(cuttypesq[i],cut*cut);
      cutneighmin = MIN(cutneighmin,cut);
      cutneighmax = MAX(cutneighmax,cut);

      if (force->pair && force->pair->ghostneigh) {
	cut = force->pair->cutghost[i][j] + skin;
	cutneighghostsq[i][j] = cut*cut;
      }
    }
  }
  cutneighmaxsq = cutneighmax * cutneighmax;

  // check other classes that can induce reneighboring in decide()
  // don't check if build_once is set

  restart_check = 0;
  if (output->restart_every) restart_check = 1;

  delete [] fixchecklist;
  fixchecklist = NULL;
  fixchecklist = new int[modify->nfix];

  fix_check = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;

  must_check = 0;
  if (restart_check || fix_check) must_check = 1;
  if (build_once) must_check = 0;

  // set special_flag for 1-2, 1-3, 1-4 neighbors
  // flag[0] is not used, flag[1] = 1-2, flag[2] = 1-3, flag[3] = 1-4
  // flag = 0 if both LJ/Coulomb special values are 0.0
  // flag = 1 if both LJ/Coulomb special values are 1.0
  // flag = 2 otherwise or if KSpace solver is enabled
  // pairwise portion of KSpace solver uses all 1-2,1-3,1-4 neighbors

  if (force->special_lj[1] == 0.0 && force->special_coul[1] == 0.0) 
    special_flag[1] = 0;
  else if (force->special_lj[1] == 1.0 && force->special_coul[1] == 1.0) 
    special_flag[1] = 1;
  else special_flag[1] = 2;

  if (force->special_lj[2] == 0.0 && force->special_coul[2] == 0.0) 
    special_flag[2] = 0;
  else if (force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0) 
    special_flag[2] = 1;
  else special_flag[2] = 2;

  if (force->special_lj[3] == 0.0 && force->special_coul[3] == 0.0) 
    special_flag[3] = 0;
  else if (force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) 
    special_flag[3] = 1;
  else special_flag[3] = 2;

  if (force->kspace) special_flag[1] = special_flag[2] = special_flag[3] = 2;

  // maxwt = max multiplicative factor on atom indices stored in neigh list

  maxwt = 0;
  if (special_flag[1] == 2) maxwt = 2;
  if (special_flag[2] == 2) maxwt = 3;
  if (special_flag[3] == 2) maxwt = 4;

  // rRESPA cutoffs

  int respa = 0;
  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  if (respa) {
    double *cut_respa = ((Respa *) update->integrate)->cutoff;
    cut_inner_sq = (cut_respa[1] + skin) * (cut_respa[1] + skin);
    cut_middle_sq = (cut_respa[3] + skin) * (cut_respa[3] + skin);
    cut_middle_inside_sq = (cut_respa[0] - skin) * (cut_respa[0] - skin);
  }

  // ------------------------------------------------------------------
  // xhold, bins, exclusion lists

  // free xhold and bins if not needed for this run

  if (dist_check == 0) {
    memory->destroy(xhold);
    maxhold = 0;
    xhold = NULL;
  }

  if (style == NSQ) {
    memory->destroy(bins);
    memory->destroy(binhead);
    maxbin = maxhead = 0;
    binhead = NULL;
    bins = NULL;
  }

  // 1st time allocation of xhold and bins

  if (dist_check) {
    if (maxhold == 0) {
      maxhold = atom->nmax;
      memory->create(xhold,maxhold,3,"neigh:xhold");
    }
  }

  if (style != NSQ) {
    if (maxbin == 0) {
      maxbin = atom->nmax;
      memory->create(bins,maxbin,"bins");
    }
  }
    
  // exclusion lists for type, group, molecule settings from neigh_modify

  n = atom->ntypes;

  if (nex_type == 0 && nex_group == 0 && nex_mol == 0) exclude = 0;
  else exclude = 1;

  if (nex_type) {
    memory->destroy(ex_type);
    memory->create(ex_type,n+1,n+1,"neigh:ex_type");

    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
	ex_type[i][j] = 0;

    for (i = 0; i < nex_type; i++) {
      if (ex1_type[i] <= 0 || ex1_type[i] > n || 
	  ex2_type[i] <= 0 || ex2_type[i] > n)
	error->all(FLERR,"Invalid atom type in neighbor exclusion list");
      ex_type[ex1_type[i]][ex2_type[i]] = 1;
      ex_type[ex2_type[i]][ex1_type[i]] = 1;
    }
  }

  if (nex_group) {
    delete [] ex1_bit;
    delete [] ex2_bit;
    ex1_bit = new int[nex_group];
    ex2_bit = new int[nex_group];

    for (i = 0; i < nex_group; i++) {
      ex1_bit[i] = group->bitmask[ex1_group[i]];
      ex2_bit[i] = group->bitmask[ex2_group[i]];
    }
  }

  if (nex_mol) {
    delete [] ex_mol_bit;
    ex_mol_bit = new int[nex_mol];

    for (i = 0; i < nex_mol; i++)
      ex_mol_bit[i] = group->bitmask[ex_mol_group[i]];
  }

  // ------------------------------------------------------------------
  // pairwise lists

  // test if pairwise lists need to be re-created
  // no need to re-create if:
  //   neigh style and triclinic has not changed and
  //   current requests = old requests

  int same = 1;
  if (style != old_style) same = 0;
  if (triclinic != old_triclinic) same = 0;
  if (nrequest != old_nrequest) same = 0;
  else
    for (i = 0; i < nrequest; i++)
      if (requests[i]->identical(old_requests[i]) == 0) same = 0;

#ifdef NEIGH_LIST_DEBUG
  if (comm->me == 0) printf("SAME flag %d\n",same);
#endif

  // if old and new are not the same, create new pairwise lists

  if (!same) {

    // delete old lists and create new ones

    for (i = 0; i < nlist; i++) delete lists[i];
    delete [] lists;
    delete [] pair_build;
    delete [] stencil_create;

    nlist = nrequest;
    lists = new NeighList*[nlist];
    pair_build = new PairPtr[nlist];
    stencil_create = new StencilPtr[nlist];

    // create individual lists, one per request
    // copy dnum setting from request to list
    // pass list ptr back to requestor (except for Command class)

    for (i = 0; i < nlist; i++) {
      lists[i] = new NeighList(lmp,pgsize);
      lists[i]->index = i;
      lists[i]->dnum = requests[i]->dnum;

      if (requests[i]->pair) {
	Pair *pair = (Pair *) requests[i]->requestor;
	pair->init_list(requests[i]->id,lists[i]);
      } else if (requests[i]->fix) {
	Fix *fix = (Fix *) requests[i]->requestor;
	fix->init_list(requests[i]->id,lists[i]);
      } else if (requests[i]->compute) {
	Compute *compute = (Compute *) requests[i]->requestor;
	compute->init_list(requests[i]->id,lists[i]);
      }
    }

    // detect lists that are connected to other lists
    // if-then-else sequence is important
    //   since don't want to re-process skip or copy lists further down
    // copy: point this list at request->otherlist, could be a skip list
    // skip: point this list at request->otherlist, copy skip info from request
    // half_from_full: point this list at preceeding full list
    // granhistory: set preceeding list's listgranhistory to this list
    //   also set preceeding list's ptr to FixShearHistory
    // respaouter: point this list at preceeding 1/2 inner/middle lists
    // pair and half: if there is a full non-occasional non-skip list
    //   change this list to half_from_full and point at the full list
    //   parent could be copy list or pair or fix
    // fix/compute requests:
    //   kind of request = half or full, occasional or not doesn't matter
    //   if request = half and non-skip pair half/respaouter exists,
    //     become copy of that list if cudable flag matches
    //   if request = full and non-skip pair full exists,
    //     become copy of that list if cudable flag matches
    //   if request = half and non-skip pair full exists,
    //     become half_from_full of that list if cudable flag matches
    //   if no matches, do nothing, fix/compute list will be built directly
    //   ok if parent is copy list

    for (i = 0; i < nlist; i++) {
      if (requests[i]->copy)
	lists[i]->listcopy = lists[requests[i]->otherlist];

      else if (requests[i]->skip) {
	lists[i]->listskip = lists[requests[i]->otherlist];
	lists[i]->copy_skip_info(requests[i]->iskip,requests[i]->ijskip);

      } else if (requests[i]->half_from_full)
	lists[i]->listfull = lists[i-1];

      else if (requests[i]->granhistory) {
	lists[i-1]->listgranhistory = lists[i];
	for (int ifix = 0; ifix < modify->nfix; ifix++)
	  if (strcmp(modify->fix[ifix]->style,"SHEAR_HISTORY") == 0) 
	    lists[i-1]->fix_history = (FixShearHistory *) modify->fix[ifix];
 
      } else if (requests[i]->respaouter) {
	if (requests[i-1]->respainner) {
	  lists[i]->respamiddle = 0;
	  lists[i]->listinner = lists[i-1];
	} else {
	  lists[i]->respamiddle = 1;
	  lists[i]->listmiddle = lists[i-1];
	  lists[i]->listinner = lists[i-2];
	}

      } else if (requests[i]->pair && requests[i]->half) {
	for (j = 0; j < nlist; j++)
	  if (requests[j]->full && requests[j]->occasional == 0 &&
	      requests[j]->skip == 0) break;
	if (j < nlist) {
	  requests[i]->half = 0;
	  requests[i]->half_from_full = 1;
	  lists[i]->listfull = lists[j];
	}

      } else if (requests[i]->fix || requests[i]->compute) {
	for (j = 0; j < nlist; j++) {
	  if (requests[i]->half && requests[j]->pair &&
	      requests[j]->skip == 0 && requests[j]->half) break;
	  if (requests[i]->full && requests[j]->pair &&
	      requests[j]->skip == 0 && requests[j]->full) break;
	  if (requests[i]->half && requests[j]->pair &&
	      requests[j]->skip == 0 && requests[j]->respaouter) break;
	}
	if (j < nlist && requests[j]->cudable != requests[i]->cudable)
	  j = nlist;
	if (j < nlist) {
	  requests[i]->copy = 1;
	  lists[i]->listcopy = lists[j];
	} else {
	  for (j = 0; j < nlist; j++) {
	    if (requests[i]->half && requests[j]->pair &&
		requests[j]->skip == 0 && requests[j]->full) break;
	  }
	  if (j < nlist && requests[j]->cudable != requests[i]->cudable)
	    j = nlist;
	  if (j < nlist) {
	    requests[i]->half = 0;
	    requests[i]->half_from_full = 1;
	    lists[i]->listfull = lists[j];
	  }
	}
      }
    }

    // set ptrs to pair_build and stencil_create functions for each list
    // ptrs set to NULL if not set explicitly
    // also set cudable to 0 if any neigh list request is not cudable

    for (i = 0; i < nlist; i++) {
      choose_build(i,requests[i]);
      if (style != NSQ) choose_stencil(i,requests[i]);
      else stencil_create[i] = NULL;
      if (!requests[i]->cudable) cudable = 0;
    }

    // set each list's build/grow/stencil/ghost flags based on neigh request
    // buildflag = 1 if its pair_build() invoked every reneighbor
    // growflag = 1 if it stores atom-based arrays and pages
    // stencilflag = 1 if it stores stencil arrays
    // ghostflag = 1 if it stores neighbors of ghosts
    // anyghostlist = 1 if any non-occasional list stores neighbors of ghosts

    anyghostlist = 0;
    for (i = 0; i < nlist; i++) {
      lists[i]->buildflag = 1;
      if (pair_build[i] == NULL) lists[i]->buildflag = 0;
      if (requests[i]->occasional) lists[i]->buildflag = 0;

      lists[i]->growflag = 1;
      if (requests[i]->copy) lists[i]->growflag = 0;

      lists[i]->stencilflag = 1;
      if (style == NSQ) lists[i]->stencilflag = 0;
      if (stencil_create[i] == NULL) lists[i]->stencilflag = 0;

      lists[i]->ghostflag = 0;
      if (requests[i]->ghost) lists[i]->ghostflag = 1;
      if (requests[i]->ghost && !requests[i]->occasional) anyghostlist = 1;
    }

#ifdef NEIGH_LIST_DEBUG
    for (i = 0; i < nlist; i++) lists[i]->print_attributes();
#endif

    // allocate atom arrays and 1st pages of lists that store them

    maxatom = atom->nmax;
    for (i = 0; i < nlist; i++)
      if (lists[i]->growflag) {
	lists[i]->grow(maxatom);
	lists[i]->add_pages();
      }

    // setup 3 vectors of pairwise neighbor lists
    // blist = lists whose pair_build() is invoked every reneighbor
    // glist = lists who store atom arrays which are used every reneighbor
    // slist = lists who store stencil arrays which are used every reneighbor
    // blist and glist vectors are used by neighbor::build()
    // slist vector is used by neighbor::setup_bins()

    nblist = nglist = nslist = 0;
    delete [] blist;
    delete [] glist;
    delete [] slist;
    blist = new int[nlist];
    glist = new int[nlist];
    slist = new int[nlist];

    for (i = 0; i < nlist; i++) {
      if (lists[i]->buildflag) blist[nblist++] = i;
      if (lists[i]->growflag && requests[i]->occasional == 0)
	glist[nglist++] = i;
      if (lists[i]->stencilflag && requests[i]->occasional == 0)
	slist[nslist++] = i;
    }

#ifdef NEIGH_LIST_DEBUG
    print_lists_of_lists();
#endif

    // reorder build vector if necessary
    // relevant for lists that copy/skip/half-full from parent
    // the derived list must appear in blist after the parent list
    // no occasional lists are in build vector
    // swap two lists within blist when dependency is mis-ordered
    // done when entire pass thru blist results in no swaps

    int done = 0;
    while (!done) {
      done = 1;
      for (i = 0; i < nblist; i++) {
	NeighList *ptr = NULL;
	if (lists[blist[i]]->listfull) ptr = lists[blist[i]]->listfull;
	if (lists[blist[i]]->listcopy) ptr = lists[blist[i]]->listcopy;
	if (lists[blist[i]]->listskip) ptr = lists[blist[i]]->listskip;
	if (ptr == NULL) continue;
	for (m = 0; m < nlist; m++)
	  if (ptr == lists[m]) break;
	for (j = 0; j < nblist; j++)
	  if (m == blist[j]) break;
	if (j < i) continue;
	int tmp = blist[i];
	blist[i] = blist[j];
	blist[j] = tmp;
	done = 0;
	break;
      }
    }

#ifdef NEIGH_LIST_DEBUG
    print_lists_of_lists();
#endif
  }

  // delete old requests
  // copy current requests and style to old for next run
  
  for (i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);
  old_nrequest = nrequest;
  old_requests = requests;
  nrequest = maxrequest = 0;
  requests = NULL;
  old_style = style;
  old_triclinic = triclinic;

  // ------------------------------------------------------------------
  // topology lists

  // 1st time allocation of topology lists

  if (atom->molecular && atom->nbonds && maxbond == 0) {
    if (nprocs == 1) maxbond = atom->nbonds;
    else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / nprocs);
    memory->create(bondlist,maxbond,3,"neigh:bondlist");
  }

  if (atom->molecular && atom->nangles && maxangle == 0) {
    if (nprocs == 1) maxangle = atom->nangles;
    else maxangle = static_cast<int> (LB_FACTOR * atom->nangles / nprocs);
    memory->create(anglelist,maxangle,4,"neigh:anglelist");
  }

  if (atom->molecular && atom->ndihedrals && maxdihedral == 0) {
    if (nprocs == 1) maxdihedral = atom->ndihedrals;
    else maxdihedral = static_cast<int> 
	   (LB_FACTOR * atom->ndihedrals / nprocs);
    memory->create(dihedrallist,maxdihedral,5,"neigh:dihedrallist");
  }

  if (atom->molecular && atom->nimpropers && maximproper == 0) {
    if (nprocs == 1) maximproper = atom->nimpropers;
    else maximproper = static_cast<int>
	   (LB_FACTOR * atom->nimpropers / nprocs);
    memory->create(improperlist,maximproper,5,"neigh:improperlist");
  }

  // set flags that determine which topology neighboring routines to use
  // SHAKE sets bonds and angles negative
  // bond_quartic sets bonds to 0
  // delete_bonds sets all interactions negative

  int bond_off = 0;
  int angle_off = 0;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"shake") == 0)
      bond_off = angle_off = 1;
  if (force->bond && force->bond_match("quartic")) bond_off = 1;

  if (atom->avec->bonds_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (bond_off) break;
      for (m = 0; m < atom->num_bond[i]; m++)
	if (atom->bond_type[i][m] <= 0) bond_off = 1;
    }
  }

  if (atom->avec->angles_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (angle_off) break;
      for (m = 0; m < atom->num_angle[i]; m++)
	if (atom->angle_type[i][m] <= 0) angle_off = 1;
    }
  }

  int dihedral_off = 0;
  if (atom->avec->dihedrals_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (dihedral_off) break;
      for (m = 0; m < atom->num_dihedral[i]; m++)
	if (atom->dihedral_type[i][m] <= 0) dihedral_off = 1;
    }
  }

  int improper_off = 0;
  if (atom->avec->impropers_allow) {
    for (i = 0; i < atom->nlocal; i++) {
      if (improper_off) break;
      for (m = 0; m < atom->num_improper[i]; m++)
	if (atom->improper_type[i][m] <= 0) improper_off = 1;
    }
  }

  // set ptrs to topology build functions

  if (bond_off) bond_build = &Neighbor::bond_partial;
  else bond_build = &Neighbor::bond_all;

  if (angle_off) angle_build = &Neighbor::angle_partial;
  else angle_build = &Neighbor::angle_all;

  if (dihedral_off) dihedral_build = &Neighbor::dihedral_partial;
  else dihedral_build = &Neighbor::dihedral_all;

  if (improper_off) improper_build = &Neighbor::improper_partial;
  else improper_build = &Neighbor::improper_all;

  // set topology neighbor list counts to 0
  // in case all are turned off but potential is still defined

  nbondlist = nanglelist = ndihedrallist = nimproperlist = 0;
}

/* ---------------------------------------------------------------------- */

int Neighbor::request(void *requestor)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **) 
      memory->srealloc(requests,maxrequest*sizeof(NeighRequest *),
		       "neighbor:requests");
  }

  requests[nrequest] = new NeighRequest(lmp);
  requests[nrequest]->requestor = requestor;
  nrequest++;
  return nrequest-1;
}

/* ----------------------------------------------------------------------
   determine which pair_build function each neigh list needs
   based on settings of neigh request
   copy -> copy_from function
   skip -> granular function if gran with granhistory,
           respa function if respaouter,
	   skip_from function for everything else
   half_from_full, half, full, gran, respaouter ->
     choose by newton and rq->newton and tri settings
     style NSQ options = newton off, newton on
     style BIN options = newton off, newton on and not tri, newton on and tri
     stlye MULTI options = same options as BIN
   if none of these, ptr = NULL since pair_build is not invoked for this list
   use "else if" b/c skip,copy can be set in addition to half,full,etc
------------------------------------------------------------------------- */

void Neighbor::choose_build(int index, NeighRequest *rq)
{
  PairPtr pb = NULL;

  if (rq->omp == 0) {

    if (rq->copy) pb = &Neighbor::copy_from;

    else if (rq->skip) {
      if (rq->gran && lists[index]->listgranhistory)
	pb = &Neighbor::skip_from_granular;
      else if (rq->respaouter) pb = &Neighbor::skip_from_respa;
      else pb = &Neighbor::skip_from;

    } else if (rq->half_from_full) {
      if (newton_pair == 0) pb = &Neighbor::half_from_full_no_newton;
      else if (newton_pair == 1) pb = &Neighbor::half_from_full_newton;

    } else if (rq->half) {
      if (style == NSQ) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_nsq_no_newton;
	  else if (newton_pair == 1) pb = &Neighbor::half_nsq_newton;
	} else if (rq->newton == 1) {
	  pb = &Neighbor::half_nsq_newton;
	} else if (rq->newton == 2) {
	  pb = &Neighbor::half_nsq_no_newton;
	}
      } else if (style == BIN) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_bin_no_newton;
	  else if (triclinic == 0) pb = &Neighbor::half_bin_newton;
	  else if (triclinic == 1) pb = &Neighbor::half_bin_newton_tri;
	} else if (rq->newton == 1) {
	  if (triclinic == 0) pb = &Neighbor::half_bin_newton;
	  else if (triclinic == 1) pb = &Neighbor::half_bin_newton_tri;
	} else if (rq->newton == 2) pb = &Neighbor::half_bin_no_newton;
      } else if (style == MULTI) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_multi_no_newton;
	  else if (triclinic == 0) pb = &Neighbor::half_multi_newton;
	  else if (triclinic == 1) pb = &Neighbor::half_multi_newton_tri;
	} else if (rq->newton == 1) {
	  if (triclinic == 0) pb = &Neighbor::half_multi_newton;
	  else if (triclinic == 1) pb = &Neighbor::half_multi_newton_tri;
	} else if (rq->newton == 2) pb = &Neighbor::half_multi_no_newton;
      }

    } else if (rq->full) {
      if (style == NSQ) {
	if (rq->ghost == 0) pb = &Neighbor::full_nsq;
	else if (includegroup) 
	  error->all(FLERR,
		     "Neighbor include group not allowed with ghost neighbors");
	else if (rq->ghost == 1) pb = &Neighbor::full_nsq_ghost;
      } else if (style == BIN) {
	if (rq->ghost == 0) pb = &Neighbor::full_bin;
	else if (includegroup) 
	  error->all(FLERR,
		     "Neighbor include group not allowed with ghost neighbors");
	else if (rq->ghost == 1) pb = &Neighbor::full_bin_ghost;
      } else if (style == MULTI) {
	if (rq->ghost == 0) pb = &Neighbor::full_multi;
	else error->all(FLERR,
			"Neighbor multi not yet enabled for ghost neighbors");
      }

    } else if (rq->gran) {
      if (style == NSQ) {
	if (newton_pair == 0) pb = &Neighbor::granular_nsq_no_newton;
	else if (newton_pair == 1) pb = &Neighbor::granular_nsq_newton;
      } else if (style == BIN) {
	if (newton_pair == 0) pb = &Neighbor::granular_bin_no_newton;
	else if (triclinic == 0) pb = &Neighbor::granular_bin_newton;
	else if (triclinic == 1) pb = &Neighbor::granular_bin_newton_tri;
      } else if (style == MULTI)
	error->all(FLERR,"Neighbor multi not yet enabled for granular");

    } else if (rq->respaouter) {
      if (style == NSQ) {
	if (newton_pair == 0) pb = &Neighbor::respa_nsq_no_newton;
	else if (newton_pair == 1) pb = &Neighbor::respa_nsq_newton;
      } else if (style == BIN) {
	if (newton_pair == 0) pb = &Neighbor::respa_bin_no_newton;
	else if (triclinic == 0) pb = &Neighbor::respa_bin_newton;
	else if (triclinic == 1) pb = &Neighbor::respa_bin_newton_tri;
      } else if (style == MULTI)
	error->all(FLERR,"Neighbor multi not yet enabled for rRESPA");
    }
  } else {

    if (rq->copy) pb = &Neighbor::copy_from;

    else if (rq->skip) {
      if (rq->gran && lists[index]->listgranhistory)
	pb = &Neighbor::skip_from_granular;
      else if (rq->respaouter) pb = &Neighbor::skip_from_respa;
      else pb = &Neighbor::skip_from;

    } else if (rq->half_from_full) {
      if (newton_pair == 0) pb = &Neighbor::half_from_full_no_newton_omp;
      else if (newton_pair == 1) pb = &Neighbor::half_from_full_newton_omp;

    } else if (rq->half) {
      if (style == NSQ) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_nsq_no_newton_omp;
	  else if (newton_pair == 1) pb = &Neighbor::half_nsq_newton_omp;
	} else if (rq->newton == 1) {
	  pb = &Neighbor::half_nsq_newton_omp;
	} else if (rq->newton == 2) {
	  pb = &Neighbor::half_nsq_no_newton_omp;
	}
      } else if (style == BIN) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_bin_no_newton_omp;
	  else if (triclinic == 0) pb = &Neighbor::half_bin_newton_omp;
	  else if (triclinic == 1) pb = &Neighbor::half_bin_newton_tri_omp;
	} else if (rq->newton == 1) {
	  if (triclinic == 0) pb = &Neighbor::half_bin_newton_omp;
	  else if (triclinic == 1) pb = &Neighbor::half_bin_newton_tri_omp;
	} else if (rq->newton == 2) pb = &Neighbor::half_bin_no_newton_omp;
      } else if (style == MULTI) {
	if (rq->newton == 0) {
	  if (newton_pair == 0) pb = &Neighbor::half_multi_no_newton_omp;
	  else if (triclinic == 0) pb = &Neighbor::half_multi_newton_omp;
	  else if (triclinic == 1) pb = &Neighbor::half_multi_newton_tri_omp;
	} else if (rq->newton == 1) {
	  if (triclinic == 0) pb = &Neighbor::half_multi_newton_omp;
	  else if (triclinic == 1) pb = &Neighbor::half_multi_newton_tri_omp;
	} else if (rq->newton == 2) pb = &Neighbor::half_multi_no_newton_omp;
      }

    } else if (rq->full) {
      if (style == NSQ) {
	if (rq->ghost == 0) pb = &Neighbor::full_nsq_omp;
	else if (includegroup) 
	  error->all(FLERR,
		     "Neighbor include group not allowed with ghost neighbors");
	else if (rq->ghost == 1) pb = &Neighbor::full_nsq_ghost_omp;
      } else if (style == BIN) {
	if (rq->ghost == 0) pb = &Neighbor::full_bin_omp;
	else if (includegroup) 
	  error->all(FLERR,
		     "Neighbor include group not allowed with ghost neighbors");
	else if (rq->ghost == 1) pb = &Neighbor::full_bin_ghost_omp;
      } else if (style == MULTI) {
	if (rq->ghost == 0) pb = &Neighbor::full_multi_omp;
	else error->all(FLERR,
			"Neighbor multi not yet enabled for ghost neighbors");
      }

    } else if (rq->gran) {
      if (style == NSQ) {
	if (newton_pair == 0) pb = &Neighbor::granular_nsq_no_newton_omp;
	else if (newton_pair == 1) pb = &Neighbor::granular_nsq_newton_omp;
      } else if (style == BIN) {
	if (newton_pair == 0) pb = &Neighbor::granular_bin_no_newton_omp;
	else if (triclinic == 0) pb = &Neighbor::granular_bin_newton_omp;
	else if (triclinic == 1) pb = &Neighbor::granular_bin_newton_tri_omp;
      } else if (style == MULTI)
	error->all(FLERR,"Neighbor multi not yet enabled for granular");

    } else if (rq->respaouter) {
      if (style == NSQ) {
	if (newton_pair == 0) pb = &Neighbor::respa_nsq_no_newton_omp;
	else if (newton_pair == 1) pb = &Neighbor::respa_nsq_newton_omp;
      } else if (style == BIN) {
	if (newton_pair == 0) pb = &Neighbor::respa_bin_no_newton_omp;
	else if (triclinic == 0) pb = &Neighbor::respa_bin_newton_omp;
	else if (triclinic == 1) pb = &Neighbor::respa_bin_newton_tri_omp;
      } else if (style == MULTI)
	error->all(FLERR,"Neighbor multi not yet enabled for rRESPA");
    }
  }

  // general error check

  if (rq->ghost && !rq->full)
    error->all(FLERR,
	       "Neighbors of ghost atoms only allowed for full neighbor lists");

  pair_build[index] = pb;
}

/* ----------------------------------------------------------------------
   determine which stencil_create function each neigh list needs
   based on settings of neigh request, only called if style != NSQ
   skip or copy or half_from_full -> no stencil
   half, gran, respaouter, full -> choose by newton and tri and dimension
   if none of these, ptr = NULL since this list needs no stencils
   use "else if" b/c skip,copy can be set in addition to half,full,etc
------------------------------------------------------------------------- */

void Neighbor::choose_stencil(int index, NeighRequest *rq)
{
  StencilPtr sc = NULL;

  if (rq->skip || rq->copy || rq->half_from_full) sc = NULL;

  else if (rq->half || rq->gran || rq->respaouter) {
    if (style == BIN) {
      if (rq->newton == 0) {
	if (newton_pair == 0) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_bin_2d_no_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_bin_3d_no_newton;
	} else if (triclinic == 0) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_bin_2d_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_bin_3d_newton;
	} else if (triclinic == 1) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_bin_2d_newton_tri;
	  else if (dimension == 3) 
	    sc = &Neighbor::stencil_half_bin_3d_newton_tri;
	}
      } else if (rq->newton == 1) {
	if (triclinic == 0) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_bin_2d_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_bin_3d_newton;
	} else if (triclinic == 1) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_bin_2d_newton_tri;
	  else if (dimension == 3) 
	    sc = &Neighbor::stencil_half_bin_3d_newton_tri;
	}
      } else if (rq->newton == 2) {
	if (dimension == 2)
	  sc = &Neighbor::stencil_half_bin_2d_no_newton;
	else if (dimension == 3)
	  sc = &Neighbor::stencil_half_bin_3d_no_newton;
      }

    } else if (style == MULTI) {
      if (rq->newton == 0) {
	if (newton_pair == 0) {
	  if (dimension == 2) 
	    sc = &Neighbor::stencil_half_multi_2d_no_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_multi_3d_no_newton;
	} else if (triclinic == 0) {
	  if (dimension == 2) 
	    sc = &Neighbor::stencil_half_multi_2d_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_multi_3d_newton;
	} else if (triclinic == 1) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_multi_2d_newton_tri;
	  else if (dimension == 3) 
	    sc = &Neighbor::stencil_half_multi_3d_newton_tri;
	}
      } else if (rq->newton == 1) {
	if (triclinic == 0) {
	  if (dimension == 2) 
	    sc = &Neighbor::stencil_half_multi_2d_newton;
	  else if (dimension == 3)
	    sc = &Neighbor::stencil_half_multi_3d_newton;
	} else if (triclinic == 1) {
	  if (dimension == 2)
	    sc = &Neighbor::stencil_half_multi_2d_newton_tri;
	  else if (dimension == 3) 
	    sc = &Neighbor::stencil_half_multi_3d_newton_tri;
	}
      } else if (rq->newton == 2) {
	if (dimension == 2) 
	  sc = &Neighbor::stencil_half_multi_2d_no_newton;
	else if (dimension == 3)
	  sc = &Neighbor::stencil_half_multi_3d_no_newton;
      }
    }

  } else if (rq->full) {
    if (style == BIN) {
      if (dimension == 2) {
	if (rq->ghost) sc = &Neighbor::stencil_full_ghost_bin_2d;
	else sc = &Neighbor::stencil_full_bin_2d;
      }
      else if (dimension == 3) {
	if (rq->ghost) sc = &Neighbor::stencil_full_ghost_bin_3d;
	else sc = &Neighbor::stencil_full_bin_3d;
      }
    } else if (style == MULTI) {
      if (dimension == 2) sc = &Neighbor::stencil_full_multi_2d;
      else if (dimension == 3) sc = &Neighbor::stencil_full_multi_3d;
    }
  }

  stencil_create[index] = sc;
}

/* ---------------------------------------------------------------------- */

void Neighbor::print_lists_of_lists()
{
  if (comm->me == 0) {
    printf("Build lists = %d: ",nblist);
    for (int i = 0; i < nblist; i++) printf("%d ",blist[i]);
    printf("\n");
    printf("Grow lists = %d: ",nglist);
    for (int i = 0; i < nglist; i++) printf("%d ",glist[i]);
    printf("\n");
    printf("Stencil lists = %d: ",nslist);
    for (int i = 0; i < nslist; i++) printf("%d ",slist[i]);
    printf("\n");
  }
}

/* ---------------------------------------------------------------------- */

int Neighbor::decide()
{
  if (must_check) {
    int n = update->ntimestep;
    if (restart_check && n == output->next_restart) return 1;
    for (int i = 0; i < fix_check; i++)
      if (n == modify->fix[fixchecklist[i]]->next_reneighbor) return 1;
  }

  ago++;
  if (ago >= delay && ago % every == 0) {
    if (build_once) return 0;
    if (dist_check == 0) return 1;
    return check_distance();
  } else return 0;
}

/* ----------------------------------------------------------------------
   if any atom moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
     compute distance each of 8 corners of box has moved since last reneighbor
     reduce skin distance by sum of 2 largest of the 8 values
     new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
------------------------------------------------------------------------- */

int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;
  double delta,deltasq,delta1,delta2;

  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx*delx + dely*dely + delz*delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx*delx + dely*dely + delz*delz);
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
	delx = corners[i][0] - corners_hold[i][0];
	dely = corners[i][1] - corners_hold[i][1];
	delz = corners[i][2] - corners_hold[i][2];
	delta = sqrt(delx*delx + dely*dely + delz*delz);
	if (delta > delta1) delta1 = delta;
	else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1+delta2));
      deltasq = delta*delta;
    }
  } else deltasq = triggersq;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - xhold[i][0];
    dely = x[i][1] - xhold[i][1];
    delz = x[i][2] - xhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > deltasq) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

/* ----------------------------------------------------------------------
   build all perpetual neighbor lists every few timesteps
   pairwise & topology lists are created as needed
------------------------------------------------------------------------- */

void Neighbor::build()
{
  int i;

  ago = 0;
  ncalls++;

  // store current atom positions and box size if needed

  if (dist_check) {
    double **x = atom->x;
    int nlocal = atom->nlocal;
    if (includegroup) nlocal = atom->nfirst;
    if (nlocal > maxhold) {
      maxhold = atom->nmax;
      memory->destroy(xhold);
      memory->create(xhold,maxhold,3,"neigh:xhold");
    }
    for (i = 0; i < nlocal; i++) {
      xhold[i][0] = x[i][0];
      xhold[i][1] = x[i][1];
      xhold[i][2] = x[i][2];
    }
    if (boxcheck) {
      if (triclinic == 0) {
	boxlo_hold[0] = bboxlo[0];
	boxlo_hold[1] = bboxlo[1];
	boxlo_hold[2] = bboxlo[2];
	boxhi_hold[0] = bboxhi[0];
	boxhi_hold[1] = bboxhi[1];
	boxhi_hold[2] = bboxhi[2];
      } else {
	domain->box_corners();
	corners = domain->corners;
	for (i = 0; i < 8; i++) {
	  corners_hold[i][0] = corners[i][0];
	  corners_hold[i][1] = corners[i][1];
	  corners_hold[i][2] = corners[i][2];
	}
      }
    }
  }

  // if any lists store neighbors of ghosts:
  // invoke grow() if nlocal+nghost exceeds previous list size
  // else only invoke grow() if nlocal exceeds previous list size
  // only done for lists with growflag set and which are perpetual

  if (anyghostlist && atom->nlocal+atom->nghost > maxatom) {
    maxatom = atom->nmax;
    for (i = 0; i < nglist; i++) lists[glist[i]]->grow(maxatom);
  } else if (atom->nlocal > maxatom) {
    maxatom = atom->nmax;
    for (i = 0; i < nglist; i++) lists[glist[i]]->grow(maxatom);
  }
  
  // extend atom bin list if necessary

  if (style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->destroy(bins);
    memory->create(bins,maxbin,"bins");
  }

  // check that neighbor list with special bond flags will not overflow

  if (atom->nlocal+atom->nghost > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // invoke building of pair and molecular neighbor lists
  // only for pairwise lists with buildflag set

  for (i = 0; i < nblist; i++)
    (this->*pair_build[blist[i]])(lists[blist[i]]);

  if (atom->molecular) {
    if (force->bond) (this->*bond_build)();
    if (force->angle) (this->*angle_build)();
    if (force->dihedral) (this->*dihedral_build)();
    if (force->improper) (this->*improper_build)();
  }
}

/* ----------------------------------------------------------------------
   build a single occasional pairwise neighbor list indexed by I
   called by other classes
------------------------------------------------------------------------- */

void Neighbor::build_one(int i)
{
  // update stencils and grow atom arrays and bins as needed
  // only for relevant settings of stencilflag and growflag
  // grow atom array for this list to current size of perpetual lists

  if (lists[i]->stencilflag) {
    lists[i]->stencil_allocate(smax,style);
    (this->*stencil_create[i])(lists[i],sx,sy,sz);
  }

  if (lists[i]->growflag) lists[i]->grow(maxatom);

  if (style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->destroy(bins);
    memory->create(bins,maxbin,"bins");
  }

  // check that neighbor list with special bond flags will not overflow

  if (atom->nlocal+atom->nghost > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // when occasional list built, LAMMPS can crash if atoms have moved too far
  // why is this?, give warning if this is the case
  // no easy workaround b/c all neighbor lists really need to be rebuilt
  // solution is for input script to check more often for rebuild
  // only check_distance if running a simulation, not between simulations

  int flag = 0;
  if (dist_check && update->whichflag) flag = check_distance();
  if (flag && me == 0)
    error->warning(FLERR,"Building an occasional neighobr list when "
		   "atoms may have moved too far");

  (this->*pair_build[i])(lists[i]);
}

/* ----------------------------------------------------------------------
   setup neighbor binning parameters
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to bbox, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*binsize to -binsize, etc
   code will work for any binsize
     since next(xyz) and stencil extend as far as necessary
     binsize = 1/2 of cutoff is roughly optimal
   for orthogonal boxes:
     a dim must be filled exactly by integer # of bins
     in periodic, procs on both sides of PBC must see same bin boundary
     in non-periodic, coord2bin() still assumes this by use of nbin xyz
   for triclinic boxes:
     tilted simulation box cannot contain integer # of bins
     stencil & neigh list built differently to account for this
   mbinlo = lowest global bin any of my ghost atoms could fall into
   mbinhi = highest global bin any of my ghost atoms could fall into
   mbin = number of bins I need in a dimension
------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by comm->cutghost
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3],bsubboxlo[3],bsubboxhi[3];
  double *cutghost = comm->cutghost;

  if (triclinic == 0) {
    bsubboxlo[0] = domain->sublo[0] - cutghost[0];
    bsubboxlo[1] = domain->sublo[1] - cutghost[1];
    bsubboxlo[2] = domain->sublo[2] - cutghost[2];
    bsubboxhi[0] = domain->subhi[0] + cutghost[0];
    bsubboxhi[1] = domain->subhi[1] + cutghost[1];
    bsubboxhi[2] = domain->subhi[2] + cutghost[2];
  } else {
    double lo[3],hi[3];
    lo[0] = domain->sublo_lamda[0] - cutghost[0];
    lo[1] = domain->sublo_lamda[1] - cutghost[1];
    lo[2] = domain->sublo_lamda[2] - cutghost[2];
    hi[0] = domain->subhi_lamda[0] + cutghost[0];
    hi[1] = domain->subhi_lamda[1] + cutghost[1];
    hi[2] = domain->subhi_lamda[2] + cutghost[2];
    domain->bbox(lo,hi,bsubboxlo,bsubboxhi);
  }

  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  // optimal bin size is roughly 1/2 the cutoff
  // for BIN style, binsize = 1/2 of max neighbor cutoff
  // for MULTI style, binsize = 1/2 of min neighbor cutoff
  // special case of all cutoffs = 0.0, binsize = box size

  double binsize_optimal;
  if (binsizeflag) binsize_optimal = binsize_user;
  else if (style == BIN) binsize_optimal = 0.5*cutneighmax;
  else binsize_optimal = 0.5*cutneighmin;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain

  if (bbox[0]*binsizeinv > MAXSMALLINT || bbox[1]*binsizeinv > MAXSMALLINT ||
      bbox[2]*binsizeinv > MAXSMALLINT)
    error->all(FLERR,"Domain too large for neighbor bins");

  // create actual bins
  // always have one bin even if cutoff > bbox
  // for 2d, nbinz = 1

  nbinx = static_cast<int> (bbox[0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1]*binsizeinv);
  if (dimension == 3) nbinz = static_cast<int> (bbox[2]*binsizeinv);
  else nbinz = 1;

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  // compute actual bin size for nbins to fit into box exactly
  // error if actual bin size << cutoff, since will create a zillion bins
  // this happens when nbin = 1 and box size << cutoff
  // typically due to non-periodic, flat system in a particular dim
  // in that extreme case, should use NSQ not BIN neighbor style

  binsizex = bbox[0]/nbinx;
  binsizey = bbox[1]/nbiny;
  binsizez = bbox[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;
  
  if (binsize_optimal*bininvx > CUT2BIN_RATIO || 
      binsize_optimal*bininvy > CUT2BIN_RATIO || 
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");
  
  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  int mbinxhi,mbinyhi,mbinzhi;
  double coord;

  coord = bsubboxlo[0] - SMALL*bbox[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*bbox[0];
  mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*bbox[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*bbox[1];
  mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

  if (dimension == 3) {
    coord = bsubboxlo[2] - SMALL*bbox[2];
    mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
    if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
    coord = bsubboxhi[2] + SMALL*bbox[2];
    mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);
  }

  // extend bins by 1 to insure stencil extent is included
  // if 2d, only 1 bin in z

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  if (dimension == 3) {
    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;
  } else mbinzlo = mbinzhi = 0;
  mbinz = mbinzhi - mbinzlo + 1;

  // memory for bin ptrs

  bigint bbin = mbinx*mbiny*mbinz;
  if (bbin > MAXSMALLINT) error->one(FLERR,"Too many neighbor bins");
  mbins = bbin;
  if (mbins > maxhead) {
    maxhead = mbins;
    memory->destroy(binhead);
    memory->create(binhead,maxhead,"neigh:binhead");
  }

  // create stencil of bins to search over in neighbor list construction
  // sx,sy,sz = max range of stencil in each dim
  // smax = max possible size of entire 3d stencil
  // stencil is empty if cutneighmax = 0.0

  sx = static_cast<int> (cutneighmax*bininvx);
  if (sx*binsizex < cutneighmax) sx++;
  sy = static_cast<int> (cutneighmax*bininvy);
  if (sy*binsizey < cutneighmax) sy++;
  sz = static_cast<int> (cutneighmax*bininvz);
  if (sz*binsizez < cutneighmax) sz++;
  if (dimension == 2) sz = 0;
  smax = (2*sx+1) * (2*sy+1) * (2*sz+1);

  // create stencils for pairwise neighbor lists
  // only done for lists with stencilflag and buildflag set

  for (int i = 0; i < nslist; i++) {
    lists[slist[i]]->stencil_allocate(smax,style);
    (this->*stencil_create[slist[i]])(lists[slist[i]],sx,sy,sz);
  }
}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double Neighbor::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;
 
  return (delx*delx + dely*dely + delz*delz);
}

/* ----------------------------------------------------------------------
   set neighbor style and skin distance
------------------------------------------------------------------------- */

void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command");

  skin = atof(arg[0]);
  if (skin < 0.0) error->all(FLERR,"Illegal neighbor command");

  if (strcmp(arg[1],"nsq") == 0) style = NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = BIN;
  else if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else error->all(FLERR,"Illegal neighbor command");
}

/* ----------------------------------------------------------------------
   modify parameters of the pair-wise neighbor build
------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      every = atoi(arg[iarg+1]);
      if (every <= 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      delay = atoi(arg[iarg+1]);
      if (delay < 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) dist_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) dist_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"once") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) build_once = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) build_once = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      pgsize = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      oneatom = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      binsize_user = atof(arg[iarg+1]);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"include") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      includegroup = group->find(arg[iarg+1]);
      if (includegroup < 0)
	error->all(FLERR,"Invalid group ID in neigh_modify command");
      if (includegroup && (atom->firstgroupname == NULL ||
			    strcmp(arg[iarg+1],atom->firstgroupname) != 0))
	error->all(FLERR,"Neigh_modify include group != atom_modify first group");
      iarg += 2;
    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");

      if (strcmp(arg[iarg+1],"type") == 0) {
	if (iarg+4 > narg) error->all(FLERR,"Illegal neigh_modify command");
	if (nex_type == maxex_type) {
	  maxex_type += EXDELTA;
	  memory->grow(ex1_type,maxex_type,"neigh:ex1_type");
	  memory->grow(ex2_type,maxex_type,"neigh:ex2_type");
	}
	ex1_type[nex_type] = atoi(arg[iarg+2]);
	ex2_type[nex_type] = atoi(arg[iarg+3]);
	nex_type++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"group") == 0) {
	if (iarg+4 > narg) error->all(FLERR,"Illegal neigh_modify command");
	if (nex_group == maxex_group) {
	  maxex_group += EXDELTA;
	  memory->grow(ex1_group,maxex_group,"neigh:ex1_group");
	  memory->grow(ex2_group,maxex_group,"neigh:ex2_group");
	}
	ex1_group[nex_group] = group->find(arg[iarg+2]);
	ex2_group[nex_group] = group->find(arg[iarg+3]);
	if (ex1_group[nex_group] == -1 || ex2_group[nex_group] == -1)
	  error->all(FLERR,"Invalid group ID in neigh_modify command");
	nex_group++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"molecule") == 0) {
	if (iarg+3 > narg) error->all(FLERR,"Illegal neigh_modify command");
	if (atom->molecule_flag == 0)
	  error->all(FLERR,"Neigh_modify exclude molecule "
		     "requires atom attribute molecule");
	if (nex_mol == maxex_mol) {
	  maxex_mol += EXDELTA;
	  memory->grow(ex_mol_group,maxex_mol,"neigh:ex_mol_group");
	}
	ex_mol_group[nex_mol] = group->find(arg[iarg+2]);
	if (ex_mol_group[nex_mol] == -1)
	  error->all(FLERR,"Invalid group ID in neigh_modify command");
	nex_mol++;
	iarg += 3;

      } else if (strcmp(arg[iarg+1],"none") == 0) {
	nex_type = nex_group = nex_mol = 0;
	iarg += 2;
      } else error->all(FLERR,"Illegal neigh_modify command");

    } else error->all(FLERR,"Illegal neigh_modify command");
  }
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void Neighbor::bin_atoms()
{
  int i,ibin;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
	ibin = coord2bin(x[i]);
	bins[i] = binhead[ibin];
	binhead[ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      ibin = coord2bin(x[i]);
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }

  } else {
    for (i = nall-1; i >= 0; i--) {
      ibin = coord2bin(x[i]);
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  }
}

/* ----------------------------------------------------------------------
   convert atom coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int Neighbor::coord2bin(double *x)
{
  int ix,iy,iz;

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;
  
  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;
  
  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/* ----------------------------------------------------------------------
   same as coord2bin, but also return ix,iy,iz offsets in each dim
------------------------------------------------------------------------- */

int Neighbor::coord2bin(double *x, int &ix, int &iy, int &iz)
{
  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = MIN(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;
  
  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = MIN(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;
  
  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = MIN(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  return iz*mbiny*mbinx + iy*mbinx + ix;
}

/* ----------------------------------------------------------------------
   test if atom pair i,j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
------------------------------------------------------------------------- */

int Neighbor::exclusion(int i, int j, int itype, int jtype,
			int *mask, int *molecule) const {
  int m;

  if (nex_type && ex_type[itype][jtype]) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask[i] & ex1_bit[m] && mask[j] & ex2_bit[m]) return 1;
      if (mask[i] & ex2_bit[m] && mask[j] & ex1_bit[m]) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)
      if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
	  molecule[i] == molecule[j]) return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Neighbor::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(xhold,maxhold,3);

  if (style != NSQ) {
    bytes += memory->usage(bins,maxbin);
    bytes += memory->usage(binhead,maxhead);
  }

  for (int i = 0; i < nlist; i++) bytes += lists[i]->memory_usage();

  bytes += memory->usage(bondlist,maxbond,3);
  bytes += memory->usage(anglelist,maxangle,4);
  bytes += memory->usage(dihedrallist,maxdihedral,5);
  bytes += memory->usage(improperlist,maximproper,5);

  return bytes;
}

/* ----------------------------------------------------------------------
   return the value of exclude - used to check compatibility with GPU
------------------------------------------------------------------------- */

int Neighbor::exclude_setting()
{
  return exclude;
}
