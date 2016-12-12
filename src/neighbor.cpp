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

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "style_nbin.h"
#include "style_nstencil.h"
#include "style_npair.h"
#include "style_ntopo.h"
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
#include "citeme.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace LAMMPS_NS;
using namespace NeighConst;

#define RQDELTA 1
#define EXDELTA 1

#define BIG 1.0e20

enum{NSQ,BIN,MULTI};     // also in NBin, NeighList, NStencil
enum{NONE,ALL,PARTIAL,TEMPLATE};

static const char cite_neigh_multi[] =
  "neighbor multi command:\n\n"
  "@Article{Intveld08,\n"
  " author =  {P.{\\,}J.~in{\\,}'t~Veld and S.{\\,}J.~Plimpton"
  " and G.{\\,}S.~Grest},\n"
  " title =   {Accurate and Efficient Methods for Modeling Colloidal\n"
  "            Mixtures in an Explicit Solvent using Molecular Dynamics},\n"
  " journal = {Comp.~Phys.~Comm.},\n"
  " year =    2008,\n"
  " volume =  179,\n"
  " pages =   {320--329}\n"
  "}\n\n";

//#define NEIGH_LIST_DEBUG 1

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  firsttime = 1;

  style = BIN;
  every = 1;
  delay = 10;
  dist_check = 1;
  pgsize = 100000;
  oneatom = 2000;
  binsizeflag = 0;
  build_once = 0;
  cluster_check = 0;
  ago = -1;

  cutneighmax = 0.0;
  cutneighsq = NULL;
  cutneighghostsq = NULL;
  cuttype = NULL;
  cuttypesq = NULL;
  fixchecklist = NULL;

  // pairwise neighbor lists and associated data structs

  nlist = 0;
  lists = NULL;

  nbin = 0;
  neigh_bin = NULL;

  nstencil = 0;
  neigh_stencil = NULL;

  neigh_pair = NULL;

  nstencil_perpetual = 0;
  slist = NULL;

  npair_perpetual = 0;
  plist = NULL;

  nrequest = maxrequest = 0;
  requests = NULL;

  old_nrequest = 0;
  old_requests = NULL;

  old_style = style;
  old_triclinic = 0;
  old_pgsize = pgsize;
  old_oneatom = oneatom;

  zeroes = NULL;

  binclass = NULL;
  binnames = NULL;
  binmasks = NULL;
  stencilclass = NULL;
  stencilnames = NULL;
  stencilmasks = NULL;

  // topology lists

  bondwhich = anglewhich = dihedralwhich = improperwhich = NONE;

  neigh_bond = NULL;
  neigh_angle = NULL;
  neigh_dihedral = NULL;
  neigh_improper = NULL;

  // coords at last neighboring

  maxhold = 0;
  xhold = NULL;
  lastcall = -1;
  last_setup_bins = -1;

  // pair exclusion list info

  includegroup = 0;

  nex_type = maxex_type = 0;
  ex1_type = ex2_type = NULL;
  ex_type = NULL;

  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = NULL;

  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = NULL;

  // Kokkos setting

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{
  if (copymode) return;

  memory->destroy(cutneighsq);
  memory->destroy(cutneighghostsq);
  delete [] cuttype;
  delete [] cuttypesq;
  delete [] fixchecklist;

  for (int i = 0; i < nlist; i++) delete lists[i];
  for (int i = 0; i < nbin; i++) delete neigh_bin[i];
  for (int i = 0; i < nstencil; i++) delete neigh_stencil[i];
  for (int i = 0; i < nlist; i++) delete neigh_pair[i];
  delete [] lists;
  delete [] neigh_bin;
  delete [] neigh_stencil;
  delete [] neigh_pair;

  delete [] slist;
  delete [] plist;

  for (int i = 0; i < nrequest; i++) delete requests[i];
  memory->sfree(requests);
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);

  delete [] zeroes;

  delete [] binclass;
  delete [] binnames;
  delete [] binmasks;
  delete [] stencilclass;
  delete [] stencilnames;
  delete [] stencilmasks;
  delete [] pairclass;
  delete [] pairnames;
  delete [] pairmasks;

  delete neigh_bond;
  delete neigh_angle;
  delete neigh_dihedral;
  delete neigh_improper;

  memory->destroy(xhold);

  memory->destroy(ex1_type);
  memory->destroy(ex2_type);
  memory->destroy(ex_type);

  memory->destroy(ex1_group);
  memory->destroy(ex2_group);
  delete [] ex1_bit;
  delete [] ex2_bit;

  memory->destroy(ex_mol_group);
  delete [] ex_mol_bit;
}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,n;

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

  // bbox lo/hi ptrs = bounding box of entire domain, stored by Domain

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
  // cutneighghost = pair cutghost if it requests it, else same as cutneigh

  triggersq = 0.25*skin*skin;
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
                             (dimension == 3 && domain->zperiodic)))
      boxcheck = 1;

  n = atom->ntypes;
  if (cutneighsq == NULL) {
    if (lmp->kokkos) init_cutneighsq_kokkos(n);
    else memory->create(cutneighsq,n+1,n+1,"neigh:cutneighsq");
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
      } else cutneighghostsq[i][j] = cut*cut;
    }
  }
  cutneighmaxsq = cutneighmax * cutneighmax;

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
    if (cut_respa[0]-skin < 0) cut_middle_inside_sq = 0.0;
  }

  // fixchecklist = other classes that can induce reneighboring in decide()

  restart_check = 0;
  if (output->restart_flag) restart_check = 1;

  delete [] fixchecklist;
  fixchecklist = NULL;
  fixchecklist = new int[modify->nfix];

  fix_check = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;

  must_check = 0;
  if (restart_check || fix_check) must_check = 1;

  // set special_flag for 1-2, 1-3, 1-4 neighbors
  // flag[0] is not used, flag[1] = 1-2, flag[2] = 1-3, flag[3] = 1-4
  // flag = 0 if both LJ/Coulomb special values are 0.0
  // flag = 1 if both LJ/Coulomb special values are 1.0
  // flag = 2 otherwise or if KSpace solver is enabled
  // pairwise portion of KSpace solver uses all 1-2,1-3,1-4 neighbors
  // or selected Coulomb-approixmation pair styles require it

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

  if (force->kspace || force->pair_match("coul/wolf",0) ||
      force->pair_match("coul/dsf",0) || force->pair_match("thole",0))
     special_flag[1] = special_flag[2] = special_flag[3] = 2;

  // maxwt = max multiplicative factor on atom indices stored in neigh list

  maxwt = 0;
  if (special_flag[1] == 2) maxwt = 2;
  if (special_flag[2] == 2) maxwt = 3;
  if (special_flag[3] == 2) maxwt = 4;

  // ------------------------------------------------------------------
  // xhold array

  // free if not needed for this run

  if (dist_check == 0) {
    memory->destroy(xhold);
    maxhold = 0;
    xhold = NULL;
  }

  // first time allocation

  if (dist_check) {
    if (maxhold == 0) {
      maxhold = atom->nmax;
      memory->create(xhold,maxhold,3,"neigh:xhold");
    }
  }

  // ------------------------------------------------------------------
  // exclusion lists

  // depend on type, group, molecule settings from neigh_modify
  // warn if exclusions used with KSpace solver

  n = atom->ntypes;

  if (nex_type == 0 && nex_group == 0 && nex_mol == 0) exclude = 0;
  else exclude = 1;

  if (nex_type) {
    if (lmp->kokkos)
      init_ex_type_kokkos(n);
    else {
      memory->destroy(ex_type);
      memory->create(ex_type,n+1,n+1,"neigh:ex_type");
    }

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
    if (lmp->kokkos)
      init_ex_bit_kokkos();
    else {
      delete [] ex1_bit;
      delete [] ex2_bit;
      ex1_bit = new int[nex_group];
      ex2_bit = new int[nex_group];
    }

    for (i = 0; i < nex_group; i++) {
      ex1_bit[i] = group->bitmask[ex1_group[i]];
      ex2_bit[i] = group->bitmask[ex2_group[i]];
    }
  }

  if (nex_mol) {
    if (lmp->kokkos)
      init_ex_mol_bit_kokkos();
    else {
      delete [] ex_mol_bit;
      ex_mol_bit = new int[nex_mol];
    }

    for (i = 0; i < nex_mol; i++)
      ex_mol_bit[i] = group->bitmask[ex_mol_group[i]];
  }

  if (exclude && force->kspace && me == 0)
    error->warning(FLERR,"Neighbor exclusions used with KSpace solver "
                   "may give inconsistent Coulombic energies");

  // ------------------------------------------------------------------
  // create pairwise lists
  // one-time call to init_styles() to scan style files and setup
  // init_pair() creates auxiliary classes: NBin, NStencil, NPair

  if (firsttime) init_styles();
  firsttime = 0;

  init_pair();

  // invoke copy_neighbor_info() in Bin,Stencil,Pair classes
  // copied once per run in case any cutoff, exclusion, special info changed

  for (i = 0; i < nbin; i++) neigh_bin[i]->copy_neighbor_info();
  for (i = 0; i < nstencil; i++) neigh_stencil[i]->copy_neighbor_info();
  for (i = 0; i < nlist; i++)
    if (neigh_pair[i]) neigh_pair[i]->copy_neighbor_info();

  if (!same && comm->me == 0) print_pairwise_info();
  requests_new2old();

  // ------------------------------------------------------------------
  // create topology lists
  // instantiated topo styles can change from run to run

  init_topology();
}

/* ----------------------------------------------------------------------
   create and initialize lists of Nbin, Nstencil, NPair classes
   lists have info on all classes in 3 style*.h files
   cannot do this in constructor, b/c too early to instantiate classes
------------------------------------------------------------------------- */

void Neighbor::init_styles()
{
  // extract info from NBin classes listed in style_nbin.h

  nbclass = 0;

#define NBIN_CLASS
#define NBinStyle(key,Class,bitmasks) nbclass++;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS

  binclass = new BinCreator[nbclass]; 
  binnames = new char*[nbclass];
  binmasks = new int[nbclass];
  nbclass = 0;

#define NBIN_CLASS
#define NBinStyle(key,Class,bitmasks) \
  binnames[nbclass] = (char *) #key; \
  binclass[nbclass] = &bin_creator<Class>; \
  binmasks[nbclass++] = bitmasks;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS

  // extract info from NStencil classes listed in style_nstencil.h

  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) nsclass++;
#include "style_nstencil.h"
#undef NStencilStyle
#undef NSTENCIL_CLASS

  stencilclass = new StencilCreator[nsclass]; 
  stencilnames = new char*[nsclass];
  stencilmasks = new int[nsclass];
  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) \
  stencilnames[nsclass] = (char *) #key; \
  stencilclass[nsclass] = &stencil_creator<Class>; \
  stencilmasks[nsclass++] = bitmasks;
#include "style_nstencil.h"
#undef NStencilStyle
#undef NSTENCIL_CLASS

  // extract info from NPair classes listed in style_npair.h

  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) npclass++;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS

  pairclass = new PairCreator[npclass]; 
  pairnames = new char*[npclass];
  pairmasks = new int[npclass];
  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) \
  pairnames[npclass] = (char *) #key; \
  pairclass[npclass] = &pair_creator<Class>; \
  pairmasks[npclass++] = bitmasks;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS
}

/* ----------------------------------------------------------------------
   create and initialize NPair classes
------------------------------------------------------------------------- */

void Neighbor::init_pair()
{
  int i,j,k,m;

  // test if pairwise lists need to be re-created
  // no need to re-create if:
  //   neigh style, triclinic, pgsize, oneatom have not changed
  //   current requests = old requests
  // first archive request params for current requests
  //   before possibly changing them below

  for (i = 0; i < nrequest; i++) requests[i]->archive();

  same = 1;
  if (style != old_style) same = 0;
  if (triclinic != old_triclinic) same = 0;
  if (pgsize != old_pgsize) same = 0;
  if (oneatom != old_oneatom) same = 0;
  if (nrequest != old_nrequest) same = 0;
  else
    for (i = 0; i < nrequest; i++)
      if (requests[i]->identical(old_requests[i]) == 0) same = 0;

#ifdef NEIGH_LIST_DEBUG
  if (comm->me == 0) printf("SAME flag %d\n",same);
#endif

  if (same) return;
  
  // delete old lists and create new ones

  for (i = 0; i < nlist; i++) delete lists[i];
  for (i = 0; i < nbin; i++) delete neigh_bin[i];
  for (i = 0; i < nstencil; i++) delete neigh_stencil[i];
  for (i = 0; i < nlist; i++) delete neigh_pair[i];
  delete [] lists;
  delete [] neigh_bin;
  delete [] neigh_stencil;
  delete [] neigh_pair;

  nlist = nrequest;
  
  lists = new NeighList*[nrequest];
  neigh_bin = new NBin*[nrequest];
  neigh_stencil = new NStencil*[nrequest];
  neigh_pair = new NPair*[nrequest];
  
  // create individual lists, one per request
  // pass list ptr back to requestor (except for Command class)
  // wait to allocate initial pages until copy lists are detected
  
  for (i = 0; i < nrequest; i++) {
    if (requests[i]->kokkos_host || requests[i]->kokkos_device)
      create_kokkos_list(i);
    else
      lists[i] = new NeighList(lmp);
    lists[i]->index = i;
    
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

  // morph requests via A,B,C rules
  // this is to avoid duplicate or inefficient builds
  // update both request and list when morph

  // (A) rule: 
  // invoke post_constructor() for all lists
  // processes copy,skip,half_from_full,granhistory,respaouter lists
  // error checks and resets internal ptrs to other lists that now exist

  for (i = 0; i < nrequest; i++)
    lists[i]->post_constructor(requests[i]);
  
  // (B) rule:
  // if request = pair, half, newton != 2
  //    and full perpetual non-skip/copy list exists,
  // then morph to half_from_full of matching parent list
  // NOTE: should be OK if parent is skip list?
  //       see build method comments
  // parent can be pair or fix, so long as perpetual fix
  // NOTE: could remove newton != 2 restriction if added 
  //       half_from_full_newtoff_ghost NPair class
  //       this would require full list having ghost info
  //       would be useful when reax/c used in hybrid mode, e.g. with airebo

  for (i = 0; i < nrequest; i++) {
    if (requests[i]->pair && requests[i]->half && requests[i]->newton != 2) {
      for (j = 0; j < nrequest; j++) {
        // Kokkos doesn't yet support half from full
        if (requests[i]->kokkos_device || requests[j]->kokkos_device) continue;
        if (requests[i]->kokkos_host || requests[j]->kokkos_host) continue;

        if (requests[j]->full && requests[j]->occasional == 0 &&
            !requests[j]->skip && !requests[j]->copy) break;
      }
      if (j < nrequest) {
        requests[i]->half = 0;
        requests[i]->half_from_full = 1;
        lists[i]->listfull = lists[j];
      }
    }
  }
      
  // (C) rule: 
  // for fix/compute requests, occasional or not does not matter
  // 1st check:
  // if request = half and non-skip/copy pair half/respaouter request exists,
  // or if request = full and non-skip/copy pair full request exists,
  // or if request = gran and non-skip/copy pair gran request exists,
  // then morph to copy of the matching parent list
  // 2nd check: only if no match to 1st check
  // if request = half and non-skip/copy pair full request exists,
  // then morph to half_from_full of the matching parent list
  // for 1st or 2nd check, parent can be copy list or pair or fix
  
  for (i = 0; i < nrequest; i++) {
    if (!requests[i]->fix && !requests[i]->compute) continue;
    for (j = 0; j < nrequest; j++) {
      // Kokkos flags must match
      if (requests[i]->kokkos_device != requests[j]->kokkos_device) continue;
      if (requests[i]->kokkos_host != requests[j]->kokkos_host) continue;

      if (requests[i]->ssa != requests[j]->ssa) continue;

      if (requests[i]->half && requests[j]->pair && 
          !requests[j]->skip && requests[j]->half && !requests[j]->copy)
        break;
      if (requests[i]->half && requests[j]->pair &&
          !requests[j]->skip && requests[j]->respaouter && !requests[j]->copy)
        break;
      if (requests[i]->full && requests[j]->pair &&
          !requests[j]->skip && requests[j]->full && !requests[j]->copy) 
        break;
      if (requests[i]->gran && requests[j]->pair &&
          !requests[j]->skip && requests[j]->gran && !requests[j]->copy)
        break;
    }
    if (j < nrequest) {
      requests[i]->copy = 1;
      requests[i]->otherlist = j;
      lists[i]->copy = 1;
      lists[i]->listcopy = lists[j];
      continue;
    }
    for (j = 0; j < nrequest; j++) {
      // Kokkos doesn't yet support half from full
      if (requests[i]->kokkos_device || requests[j]->kokkos_device) continue;
      if (requests[i]->kokkos_host || requests[j]->kokkos_host) continue;

      if (requests[i]->half && requests[j]->pair &&
          !requests[j]->skip && requests[j]->full && !requests[j]->copy)
        break;
    }
    if (j < nrequest) {
      requests[i]->half = 0;
      requests[i]->half_from_full = 1;
      lists[i]->listfull = lists[j];
    }
  }

  // assign Bin,Stencil,Pair style to each list
  
  int flag;
  for (i = 0; i < nrequest; i++) {
    flag = choose_bin(requests[i]);
    lists[i]->bin_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor bin option does not exist");

    flag = choose_stencil(requests[i]);
    lists[i]->stencil_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor stencil method does not exist");

    flag = choose_pair(requests[i]);
    lists[i]->pair_method = flag;
    if (flag < 0) 
      error->all(FLERR,"Requested neighbor pair method does not exist");
  }

  // instantiate unique Bin,Stencil classes in neigh_bin & neigh_stencil vecs
  // instantiate one Pair class per list in neigh_pair vec

  nbin = 0;
  for (i = 0; i < nrequest; i++) {
    flag = lists[i]->bin_method;
    if (flag == 0) continue;
    for (j = 0; j < nbin; j++)
      if (neigh_bin[j]->istyle == flag) break;
    if (j < nbin) continue;
    BinCreator bin_creator = binclass[flag-1];
    neigh_bin[nbin] = bin_creator(lmp);
    neigh_bin[nbin]->istyle = flag;
    nbin++;
  }

  nstencil = 0;
  for (i = 0; i < nrequest; i++) {
    flag = lists[i]->stencil_method;
    if (flag == 0) continue;
    for (j = 0; j < nstencil; j++)
      if (neigh_stencil[j]->istyle == flag) break;
    if (j < nstencil) continue;
    StencilCreator stencil_creator = stencilclass[flag-1];
    neigh_stencil[nstencil] = stencil_creator(lmp);
    neigh_stencil[nstencil]->istyle = flag;
    int bin_method = lists[i]->bin_method;
    for (k = 0; k < nbin; k++) {
      if (neigh_bin[k]->istyle == bin_method) {
        neigh_stencil[nstencil]->nb = neigh_bin[k];
        break;
      }
    }
    if (k == nbin) 
      error->all(FLERR,"Could not assign bin method to neighbor stencil");
    nstencil++;
  }

  for (i = 0; i < nrequest; i++) {
    flag = lists[i]->pair_method;
    if (flag == 0) {
      neigh_pair[i] = NULL;
      continue;
    }
    PairCreator pair_creator = pairclass[flag-1];
    neigh_pair[i] = pair_creator(lmp);
    neigh_pair[i]->istyle = flag;

    int bin_method = lists[i]->bin_method;
    if (bin_method == 0) neigh_pair[i]->nb = NULL;
    else {
      for (k = 0; k < nbin; k++) {
        if (neigh_bin[k]->istyle == bin_method) {
          neigh_pair[i]->nb = neigh_bin[k];
          break;
        }
      }
      if (k == nbin) 
        error->all(FLERR,"Could not assign bin method to neighbor pair");
    }
    int stencil_method = lists[i]->stencil_method;
    if (stencil_method == 0) neigh_pair[i]->ns = NULL;
    else {
      for (k = 0; k < nstencil; k++) {
        if (neigh_stencil[k]->istyle == stencil_method) {
          neigh_pair[i]->ns = neigh_stencil[k];
          break;
        }
      }
      if (k == nstencil) 
        error->all(FLERR,"Could not assign stencil method to neighbor pair");
    }
  }

  // allocate initial pages for each list, except if copy flag set
  // allocate dnum vector of zeroes if set
  
  int dnummax = 0;
  for (i = 0; i < nlist; i++) {
    if (lists[i]->copy) continue;
    lists[i]->setup_pages(pgsize,oneatom);
    dnummax = MAX(dnummax,lists[i]->dnum);
  }
  
  if (dnummax) {
    delete [] zeroes;
    zeroes = new double[dnummax];
    for (i = 0; i < dnummax; i++) zeroes[i] = 0.0;
  }

  // first-time allocation of per-atom data for lists that are built and store
  // lists that are not built: granhistory, respa inner/middle (no neigh_pair)
  // lists that do not store: copy 
  // use atom->nmax for both grow() args
  //   i.e. grow first time to expanded size to avoid future reallocs
  // also Kokkos list initialization
  
  int maxatom = atom->nmax;
  for (i = 0; i < nlist; i++)
    if (neigh_pair[i] && !lists[i]->copy) lists[i]->grow(maxatom,maxatom);

  // plist = indices of perpetual NPair classes
  //         perpetual = non-occasional, re-built at every reneighboring
  // slist = indices of perpetual NStencil classes
  //         perpetual = used by any perpetual NPair class
  
  delete [] slist;
  delete [] plist;
  nstencil_perpetual = npair_perpetual = 0;
  slist = new int[nstencil];
  plist = new int[nlist];

  for (i = 0; i < nlist; i++) {
    if (lists[i]->occasional == 0 && lists[i]->pair_method)
      plist[npair_perpetual++] = i;
  }
  
  for (i = 0; i < nstencil; i++) {
    flag = 0;
    for (j = 0; j < npair_perpetual; j++)
      if (lists[plist[j]]->stencil_method == neigh_stencil[i]->istyle) flag = 1;
    if (flag) slist[nstencil_perpetual++] = i;
  }

  // reorder plist vector if necessary
  // relevant for lists that copy/skip/half-full from parent
  // the child index must appear in plist after the parent index
  // swap two indices within plist when dependency is mis-ordered
  // done when entire pass thru plist results in no swaps
  
  NeighList *ptr;

  int done = 0;
  while (!done) {
    done = 1;
    for (i = 0; i < npair_perpetual; i++) {
      ptr = NULL;
      if (lists[plist[i]]->listcopy) ptr = lists[plist[i]]->listcopy;
      if (lists[plist[i]]->listskip) ptr = lists[plist[i]]->listskip;
      if (lists[plist[i]]->listfull) ptr = lists[plist[i]]->listfull;
      if (ptr == NULL) continue;
      for (m = 0; m < nrequest; m++)
        if (ptr == lists[m]) break;
      for (j = 0; j < npair_perpetual; j++)
        if (m == plist[j]) break;
      if (j < i) continue;
      int tmp = plist[i];     // swap I,J indices
      plist[i] = plist[j];
      plist[j] = tmp;
      done = 0;
      break;
    }
  }

  // debug output

#ifdef NEIGH_LIST_DEBUG
  for (i = 0; i < nrequest; i++) lists[i]->print_attributes();
#endif
}

/* ----------------------------------------------------------------------
   create and initialize NTopo classes
------------------------------------------------------------------------- */

void Neighbor::init_topology()
{
  int i,m;

  if (!atom->molecular) return;

  // set flags that determine which topology neighbor classes to use
  // these settings could change from run to run, depending on fixes defined
  // bonds,etc can only be broken for atom->molecular = 1, not 2
  // SHAKE sets bonds and angles negative
  // gcmc sets all bonds, angles, etc negative
  // bond_quartic sets bonds to 0
  // delete_bonds sets all interactions negative

  int bond_off = 0;
  int angle_off = 0;
  for (i = 0; i < modify->nfix; i++)
    if ((strcmp(modify->fix[i]->style,"shake") == 0)
        || (strcmp(modify->fix[i]->style,"rattle") == 0))
      bond_off = angle_off = 1;
  if (force->bond && force->bond_match("quartic")) bond_off = 1;

  if (atom->avec->bonds_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (bond_off) break;
      for (m = 0; m < atom->num_bond[i]; m++)
        if (atom->bond_type[i][m] <= 0) bond_off = 1;
    }
  }

  if (atom->avec->angles_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (angle_off) break;
      for (m = 0; m < atom->num_angle[i]; m++)
        if (atom->angle_type[i][m] <= 0) angle_off = 1;
    }
  }

  int dihedral_off = 0;
  if (atom->avec->dihedrals_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (dihedral_off) break;
      for (m = 0; m < atom->num_dihedral[i]; m++)
        if (atom->dihedral_type[i][m] <= 0) dihedral_off = 1;
    }
  }

  int improper_off = 0;
  if (atom->avec->impropers_allow && atom->molecular == 1) {
    for (i = 0; i < atom->nlocal; i++) {
      if (improper_off) break;
      for (m = 0; m < atom->num_improper[i]; m++)
        if (atom->improper_type[i][m] <= 0) improper_off = 1;
    }
  }

  for (i = 0; i < modify->nfix; i++)
    if ((strcmp(modify->fix[i]->style,"gcmc") == 0))
      bond_off = angle_off = dihedral_off = improper_off = 1;

  // sync on/off settings across all procs

  int onoff = bond_off;
  MPI_Allreduce(&onoff,&bond_off,1,MPI_INT,MPI_MAX,world);
  onoff = angle_off;
  MPI_Allreduce(&onoff,&angle_off,1,MPI_INT,MPI_MAX,world);
  onoff = dihedral_off;
  MPI_Allreduce(&onoff,&dihedral_off,1,MPI_INT,MPI_MAX,world);
  onoff = improper_off;
  MPI_Allreduce(&onoff,&improper_off,1,MPI_INT,MPI_MAX,world);

  // instantiate NTopo classes

  if (atom->avec->bonds_allow) {
    int old_bondwhich = bondwhich;
    if (atom->molecular == 2) bondwhich = TEMPLATE;
    else if (bond_off) bondwhich = PARTIAL;
    else bondwhich = ALL;
    if (!neigh_bond || bondwhich != old_bondwhich) {
      delete neigh_bond;
      if (bondwhich == ALL) 
        neigh_bond = new NTopoBondAll(lmp);
      else if (bondwhich == PARTIAL) 
        neigh_bond = new NTopoBondPartial(lmp);
      else if (bondwhich == TEMPLATE) 
        neigh_bond = new NTopoBondTemplate(lmp);
    }
  }

  if (atom->avec->angles_allow) {
    int old_anglewhich = anglewhich;
    if (atom->molecular == 2) anglewhich = TEMPLATE;
    else if (angle_off) anglewhich = PARTIAL;
    else anglewhich = ALL;
    if (!neigh_angle || anglewhich != old_anglewhich) {
      delete neigh_angle;
      if (anglewhich == ALL) 
        neigh_angle = new NTopoAngleAll(lmp);
      else if (anglewhich == PARTIAL) 
        neigh_angle = new NTopoAnglePartial(lmp);
      else if (anglewhich == TEMPLATE) 
        neigh_angle = new NTopoAngleTemplate(lmp);
    }
  }

  if (atom->avec->dihedrals_allow) {
    int old_dihedralwhich = dihedralwhich;
    if (atom->molecular == 2) dihedralwhich = TEMPLATE;
    else if (dihedral_off) dihedralwhich = PARTIAL;
    else dihedralwhich = ALL;
    if (!neigh_dihedral || dihedralwhich != old_dihedralwhich) {
      delete neigh_dihedral;
      if (dihedralwhich == ALL) 
        neigh_dihedral = new NTopoDihedralAll(lmp);
      else if (dihedralwhich == PARTIAL) 
        neigh_dihedral = new NTopoDihedralPartial(lmp);
      else if (dihedralwhich == TEMPLATE) 
        neigh_dihedral = new NTopoDihedralTemplate(lmp);
    }
  }

  if (atom->avec->impropers_allow) {
    int old_improperwhich = improperwhich;
    if (atom->molecular == 2) improperwhich = TEMPLATE;
    else if (improper_off) improperwhich = PARTIAL;
    else improperwhich = ALL;
    if (!neigh_improper || improperwhich != old_improperwhich) {
      delete neigh_improper;
      if (improperwhich == ALL) 
        neigh_improper = new NTopoImproperAll(lmp);
      else if (improperwhich == PARTIAL) 
        neigh_improper = new NTopoImproperPartial(lmp);
      else if (improperwhich == TEMPLATE) 
        neigh_improper = new NTopoImproperTemplate(lmp);
    }
  }
}

/* ----------------------------------------------------------------------
   output summary of pairwise neighbor list info
   only called by proc 0
------------------------------------------------------------------------- */

void Neighbor::print_pairwise_info()
{
  int i,j,m;
  char str[128];
  const char *kind;
  FILE *out;

  const double cutghost = MAX(cutneighmax,comm->cutghostuser);

  double binsize, bbox[3];
  bbox[0] =  bboxhi[0]-bboxlo[0];
  bbox[1] =  bboxhi[1]-bboxlo[1];
  bbox[2] =  bboxhi[2]-bboxlo[2];
  if (binsizeflag) binsize = binsize_user;
  else if (style == BIN) binsize = 0.5*cutneighmax;
  else binsize = 0.5*cutneighmin;
  if (binsize == 0.0) binsize = bbox[0];

  int nperpetual = 0;
  int noccasional = 0;
  int nextra = 0;
  for (i = 0; i < nlist; i++) {
    if (lists[i]->pair_method == 0) nextra++;
    else if (lists[i]->occasional) noccasional++;
    else nperpetual++;
  }

  for (m = 0; m < 2; m++) {
    if (m == 0) out = screen;
    else out = logfile;

    if (out) {
      fprintf(out,"Neighbor list info ...\n");
      fprintf(out,"  update every %d steps, delay %d steps, check %s\n",
              every,delay,dist_check ? "yes" : "no");
      fprintf(out,"  max neighbors/atom: %d, page size: %d\n",
              oneatom, pgsize);
      fprintf(out,"  master list distance cutoff = %g\n",cutneighmax);
      fprintf(out,"  ghost atom cutoff = %g\n",cutghost);
      if (style != NSQ)
        fprintf(out,"  binsize = %g, bins = %g %g %g\n",binsize,
                ceil(bbox[0]/binsize), ceil(bbox[1]/binsize),
                ceil(bbox[2]/binsize));
      
      fprintf(out,"  %d neighbor lists, "
              "perpetual/occasional/extra = %d %d %d\n",
              nlist,nperpetual,noccasional,nextra);
      
      for (i = 0; i < nlist; i++) {
        if (requests[i]->pair) {
          char *pname = force->pair_match_ptr((Pair *) requests[i]->requestor);
          sprintf(str,"  (%d) pair %s",i+1,pname);
        } else if (requests[i]->fix) {
          sprintf(str,"  (%d) fix %s",i+1,
                  ((Fix *) requests[i]->requestor)->style);
        } else if (requests[i]->compute) {
          sprintf(str,"  (%d) compute %s",i+1,
                  ((Compute *) requests[i]->requestor)->style);
        } else {
          sprintf(str,"  (%d) command %s",i+1,requests[i]->command_style);
        }
        fprintf(out,"%s\n",str);
        
        if (requests[i]->half) kind = "half";
        else if (requests[i]->full) kind = "full";
        else if (requests[i]->gran) kind = "size";
        else if (requests[i]->granhistory) kind = "size/history";
        else if (requests[i]->respainner) kind = "respa/inner";
        else if (requests[i]->respamiddle) kind = "respa/middle";
        else if (requests[i]->respaouter) kind = "respa/outer";
        else if (requests[i]->half_from_full) kind = "half/from/full";
        if (requests[i]->occasional) fprintf(out,", occasional");
        else fprintf(out,", perpetual");
        if (requests[i]->ghost) fprintf(out,", ghost");
        if (requests[i]->ssa) fprintf(out,", ssa");
        if (requests[i]->omp) fprintf(out,", omp");
        if (requests[i]->intel) fprintf(out,", intel");
        if (requests[i]->kokkos_device) fprintf(out,", kokkos_device");
        if (requests[i]->kokkos_host) fprintf(out,", kokkos_host");
        if (requests[i]->copy) 
          fprintf(out,", copy from (%d)",requests[i]->otherlist+1);
        if (requests[i]->skip)
          fprintf(out,", skip from (%d)",requests[i]->otherlist+1);
        if (requests[i]->off2on) fprintf(out,", off2on");
        fprintf(out,"\n");
        
        if (lists[i]->pair_method == 0) fprintf(out,"      pair build: none\n");
        else fprintf(out,"      pair build: %s\n",
                     pairnames[lists[i]->pair_method-1]);
        
        if (lists[i]->stencil_method == 0) fprintf(out,"      stencil: none\n");
        else fprintf(out,"      stencil: %s\n",
                     stencilnames[lists[i]->stencil_method-1]);
        
        if (lists[i]->bin_method == 0) fprintf(out,"      bin: none\n");
        else fprintf(out,"      bin: %s\n",binnames[lists[i]->bin_method-1]);
      }
      
      /*
      fprintf(out,"  %d stencil methods\n",nstencil);
      for (i = 0; i < nstencil; i++)
        fprintf(out,"    (%d) %s\n",
        i+1,stencilnames[neigh_stencil[i]->istyle-1]);

      fprintf(out,"  %d bin methods\n",nbin);
      for (i = 0; i < nbin; i++)
        fprintf(out,"    (%d) %s\n",i+1,binnames[neigh_bin[i]->istyle-1]);
      */
    }
  }
}
  
/* ----------------------------------------------------------------------
   delete old NeighRequests
   copy current requests and params to old for next run
------------------------------------------------------------------------- */

void Neighbor::requests_new2old()
{
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);

  old_nrequest = nrequest;
  old_requests = requests;
  nrequest = maxrequest = 0;
  requests = NULL;
  old_style = style;
  old_triclinic = triclinic;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
}

/* ----------------------------------------------------------------------
   assign NBin class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known Nbin classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
------------------------------------------------------------------------- */

int Neighbor::choose_bin(NeighRequest *rq)
{
  // no binning needed

  if (style == NSQ) return 0;
  if (rq->skip || rq->copy || rq->half_from_full) return 0;
  if (rq->granhistory) return 0;
  if (rq->respainner || rq->respamiddle) return 0;

  // flags for settings the request + system requires of NBin class
  //   ssaflag = no/yes ssa request
  //   intelflag = no/yes intel request
  //   kokkos_device_flag = no/yes kokkos device request
  //   kokkos_host_flag = no/yes kokkos host request

  int ssaflag,intelflag,kokkos_device_flag,kokkos_host_flag;

  ssaflag = intelflag = kokkos_device_flag = kokkos_host_flag = 0;

  if (rq->ssa) ssaflag = NB_SSA;
  if (rq->intel) intelflag = NB_INTEL;
  if (rq->kokkos_device) kokkos_device_flag = NB_KOKKOS_DEVICE;
  if (rq->kokkos_host) kokkos_host_flag = NB_KOKKOS_HOST;

  // use flags to match exactly one of NBin class masks, bit by bit

  int mask;

  for (int i = 0; i < nbclass; i++) {
    mask = binmasks[i];

    if (ssaflag != (mask & NB_SSA)) continue;
    if (intelflag != (mask & NB_INTEL)) continue;
    if (kokkos_device_flag != (mask & NB_KOKKOS_DEVICE)) continue;
    if (kokkos_host_flag != (mask & NB_KOKKOS_HOST)) continue;

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   assign NStencil class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known NStencil classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
------------------------------------------------------------------------- */

int Neighbor::choose_stencil(NeighRequest *rq)
{
  // no stencil creation needed

  if (style == NSQ) return 0;
  if (rq->skip || rq->copy || rq->half_from_full) return 0;
  if (rq->granhistory) return 0;
  if (rq->respainner || rq->respamiddle) return 0;

  // flags for settings the request + system requires of NStencil class
  //   halfflag = half request (gran and respa are also half lists)
  //   fullflag = full request
  //   ghostflag = no/yes ghost request
  //   ssaflag = no/yes ssa request
  //   dimension = 2d/3d
  //   newtflag = newton off/on request
  //   triclinic = orthgonal/triclinic box

  int halfflag,fullflag,ghostflag,ssaflag;

  halfflag = fullflag = ghostflag = ssaflag = 0;

  if (rq->half) halfflag = 1;
  if (rq->full) fullflag = 1;
  if (rq->gran) halfflag = 1;
  if (rq->respaouter) halfflag = 1;

  if (rq->ghost) ghostflag = NS_GHOST;
  if (rq->ssa) ssaflag = NS_SSA;

  int newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = 1;
  else if (rq->newton == 0 && !newton_pair) newtflag = 0;
  else if (rq->newton == 1) newtflag = 1;
  else if (rq->newton == 2) newtflag = 0;



  // use flags to match exactly one of NStencil class masks, bit by bit
  // exactly one of halfflag,fullflag is set and thus must match

  int mask;

  for (int i = 0; i < nsclass; i++) {
    mask = stencilmasks[i];

    if (halfflag) {
      if (!(mask & NS_HALF)) continue;
    } else if (fullflag) {
      if (!(mask & NS_FULL)) continue;
    }

    if (ghostflag != (mask & NS_GHOST)) continue;
    if (ssaflag != (mask & NS_SSA)) continue;

    if (style == BIN && !(mask & NS_BIN)) continue;
    if (style == MULTI && !(mask & NS_MULTI)) continue;


    if (dimension == 2 && !(mask & NS_2D)) continue;
    if (dimension == 3 && !(mask & NS_3D)) continue;

    if (newtflag && !(mask & NS_NEWTON)) continue;
    if (!newtflag && !(mask & NS_NEWTOFF)) continue;

    if (!triclinic && !(mask & NS_ORTHO)) continue;
    if (triclinic && !(mask & NS_TRI)) continue;

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   assign NPair class to a NeighList
   use neigh request settings to build mask
   match mask to list of masks of known NPair classes
   return index+1 of match in list of masks
   return 0 for no binning
   return -1 if no match
------------------------------------------------------------------------- */

int Neighbor::choose_pair(NeighRequest *rq)
{
  // no NPair build performed

  if (rq->granhistory) return 0;
  if (rq->respainner || rq->respamiddle) return 0;

  // error check for includegroup with ghost neighbor request

  if (includegroup && rq->ghost)
    error->all(FLERR,"Neighbor include group not allowed "
               "with ghost neighbors");

  // flags for settings the request + system requires of NPair class
  //   copyflag = no/yes copy request
  //   skipflag = no/yes skip request
  //   halfflag = half request (gran and respa are also half lists)
  //   fullflag = full request
  //   halffullflag = half_from_full request
  //   sizeflag = no/yes gran request for finite-size particles
  //   ghostflag = no/yes ghost request
  //   respaflag = no/yes respa request
  //   off2onflag = no/yes off2on request
  //   onesideflag = no/yes granonesided request
  //   ssaflag = no/yes request
  //   ompflag = no/yes omp request
  //   intelflag = no/yes intel request
  //   kokkos_device_flag = no/yes Kokkos device request
  //   kokkos_host_flag = no/yes Kokkos host request
  //   newtflag = newton off/on request
  //   style = NSQ/BIN/MULTI neighbor style
  //   triclinic = orthgonal/triclinic box

  int copyflag,skipflag,halfflag,fullflag,halffullflag,sizeflag,respaflag,
    ghostflag,off2onflag,onesideflag,ssaflag,ompflag,intelflag,kokkos_device_flag,kokkos_host_flag;

  copyflag = skipflag = halfflag = fullflag = halffullflag = sizeflag = 
    ghostflag = respaflag = off2onflag = onesideflag = ssaflag = 
    ompflag = intelflag = kokkos_device_flag = kokkos_host_flag = 0;

  if (rq->copy) copyflag = NP_COPY;
  if (rq->skip) skipflag = NP_SKIP;

  // NOTE: exactly one of these request flags is set (see neigh_request.h)
  //       this requires gran/respaouter also set halfflag
  //       can simplify this logic, if follow NOTE in neigh_request.h
  //       all why do size/off2on and size/off2on/oneside set NP_HALF
  //         either should set both half & full, or half should be in file name
  //         to be consistent with how other NP classes use "half"

  if (rq->half) halfflag = 1;
  if (rq->full) fullflag = 1;
  if (rq->half_from_full) halffullflag = 1;
  if (rq->gran) {
    sizeflag = NP_SIZE;
    halfflag = 1;
  }
  if (rq->respaouter) {
    respaflag = NP_RESPA;
    halfflag = 1;
  }

  if (rq->ghost) ghostflag = NP_GHOST;
  if (rq->off2on) off2onflag = NP_OFF2ON;
  if (rq->granonesided) onesideflag = NP_ONESIDE;
  if (rq->ssa) ssaflag = NP_SSA;
  if (rq->omp) ompflag = NP_OMP;
  if (rq->intel) intelflag = NP_INTEL;
  if (rq->kokkos_device) kokkos_device_flag = NP_KOKKOS_DEVICE;
  if (rq->kokkos_host) kokkos_host_flag = NP_KOKKOS_HOST;

  int newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = 1;
  else if (rq->newton == 0 && !newton_pair) newtflag = 0;
  else if (rq->newton == 1) newtflag = 1;
  else if (rq->newton == 2) newtflag = 0;

  // use flags to match exactly one of NPair class masks, bit by bit
  // copyflag match returns with no further checks
  // exactly one of halfflag,fullflag,halffullflag is set and thus must match

  int mask;

  //printf("FLAGS: %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
  //       copyflag,skipflag,halfflag,fullflag,halffullflag,
  //       sizeflag,respaflag,ghostflag,off2onflag,onesideflag,ssaflag,
  //       ompflag,intelflag,newtflag);

  for (int i = 0; i < npclass; i++) {
    mask = pairmasks[i];

    if (copyflag && (mask & NP_COPY)) {
      if (kokkos_device_flag != (mask & NP_KOKKOS_DEVICE)) continue;
      if (kokkos_host_flag != (mask & NP_KOKKOS_HOST)) continue;
      return i+1;
    }
    if (skipflag != (mask & NP_SKIP)) continue;

    if (halfflag) {
      if (!(mask & NP_HALF)) continue;
    } else if (fullflag) {
      if (!(mask & NP_FULL)) continue;
    } else if (halffullflag) {
      if (!(mask & NP_HALFFULL)) continue;
    }

    if (sizeflag != (mask & NP_SIZE)) continue;
    if (respaflag != (mask & NP_RESPA)) continue;
    if (ghostflag != (mask & NP_GHOST)) continue;
    if (off2onflag != (mask & NP_OFF2ON)) continue;
    if (onesideflag != (mask & NP_ONESIDE)) continue;
    if (ssaflag != (mask & NP_SSA)) continue;
    if (ompflag != (mask & NP_OMP)) continue;
    if (intelflag != (mask & NP_INTEL)) continue;
    if (kokkos_device_flag != (mask & NP_KOKKOS_DEVICE)) continue;
    if (kokkos_host_flag != (mask & NP_KOKKOS_HOST)) continue;

    if (style == NSQ && !(mask & NP_NSQ)) continue;
    if (style == BIN && !(mask & NP_BIN)) continue;
    if (style == MULTI && !(mask & NP_MULTI)) continue;

    if (newtflag && !(mask & NP_NEWTON)) continue;
    if (!newtflag && !(mask & NP_NEWTOFF)) continue;

    if (!triclinic && !(mask & NP_ORTHO)) continue;
    if (triclinic && !(mask & NP_TRI)) continue;

    return i+1;
  }

  //printf("NO MATCH\n");

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   called by other classes to request a pairwise neighbor list
------------------------------------------------------------------------- */

int Neighbor::request(void *requestor, int instance)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)
      memory->srealloc(requests,maxrequest*sizeof(NeighRequest *),
                       "neighbor:requests");
  }

  requests[nrequest] = new NeighRequest(lmp);
  requests[nrequest]->index = nrequest;
  requests[nrequest]->requestor = requestor;
  requests[nrequest]->requestor_instance = instance;
  nrequest++;
  return nrequest-1;
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_bin.h
------------------------------------------------------------------------- */

template <typename T>
NBin *Neighbor::bin_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_stencil.h
------------------------------------------------------------------------- */

template <typename T>
NStencil *Neighbor::stencil_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   one instance per entry in style_neigh_pair.h
------------------------------------------------------------------------- */

template <typename T>
NPair *Neighbor::pair_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   setup neighbor binning and neighbor stencils
   called before run and every reneighbor if box size/shape changes
   only operates on perpetual lists
   build_one() operates on occasional lists
------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  // invoke setup_bins() for all NBin
  // actual binning is performed in build()

  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->setup_bins(style);

  // invoke create_setup() and create() for all perpetual NStencil
  // same ops performed for occasional lists in build_one()

  for (int i = 0; i < nstencil_perpetual; i++) {
    neigh_stencil[slist[i]]->create_setup();
    neigh_stencil[slist[i]]->create();
  }

  last_setup_bins = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int Neighbor::decide()
{
  if (must_check) {
    bigint n = update->ntimestep;
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
   build perpetual neighbor lists
   called at setup and every few timesteps during run or minimization
   topology lists also built if topoflag = 1 (Kokkos calls with topoflag=0)
------------------------------------------------------------------------- */

void Neighbor::build(int topoflag)
{
  int i,m;

  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // check that using special bond flags will not overflow neigh lists

  if (nall > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");

  // store current atom positions and box size if needed

  if (dist_check) {
    double **x = atom->x;
    if (includegroup) nlocal = atom->nfirst;
    if (atom->nmax > maxhold) {
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

  // bin atoms for all NBin instances
  // not just NBin associated with perpetual lists
  // b/c cannot wait to bin occasional lists in build_one() call
  // if bin then, atoms may have moved outside of proc domain & bin extent,
  //   leading to errors or even a crash

  if (style != NSQ) {
    for (int i = 0; i < nbin; i++) {
      neigh_bin[i]->bin_atoms_setup(nall);
      neigh_bin[i]->bin_atoms();
    }
  }

  // build pairwise lists for all perpetual NPair/NeighList
  // grow() with nlocal/nall args so that only realloc if have to

  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    lists[m]->grow(nlocal,nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }

  // build topology lists for bonds/angles/etc

  if (atom->molecular && topoflag) build_topology();
}

/* ----------------------------------------------------------------------
   build topology neighbor lists: bond, angle, dihedral, improper
   copy their list info back to Neighbor for access by bond/angle/etc classes
------------------------------------------------------------------------- */

void Neighbor::build_topology()
{
  if (force->bond) {
    neigh_bond->build();
    nbondlist = neigh_bond->nbondlist;
    bondlist = neigh_bond->bondlist;
  }
  if (force->angle) {
    neigh_angle->build();
    nanglelist = neigh_angle->nanglelist;
    anglelist = neigh_angle->anglelist;
  }
  if (force->dihedral) {
    neigh_dihedral->build();
    ndihedrallist = neigh_dihedral->ndihedrallist;
    dihedrallist = neigh_dihedral->dihedrallist;
  }
  if (force->improper) {
    neigh_improper->build();
    nimproperlist = neigh_improper->nimproperlist;
    improperlist = neigh_improper->improperlist;
  }
}

/* ----------------------------------------------------------------------
   build a single occasional pairwise neighbor list indexed by I
   called by other classes
------------------------------------------------------------------------- */

void Neighbor::build_one(class NeighList *mylist, int preflag)
{
  // check if list structure is initialized

  if (mylist == NULL)
    error->all(FLERR,"Trying to build an occasional neighbor list "
               "before initialization completed");

  // build_one() should never be invoked on a perpetual list

  if (!mylist->occasional) 
    error->all(FLERR,"Neighbor build one invoked on perpetual list");

  // no need to build if already built since last re-neighbor
  // preflag is set by fix bond/create and fix bond/swap
  //   b/c they invoke build_one() on same step neigh list is re-built,
  //   but before re-build, so need to use ">" instead of ">="

  NPair *np = neigh_pair[mylist->index];

  if (preflag) {
    if (np->last_build > lastcall) return;
  } else {
    if (np->last_build >= lastcall) return;
  }

  // if this is copy list and parent is occasional list,
  // or this is half_from_full and parent is occasional list,
  // insure parent is current 

  if (mylist->listcopy && mylist->listcopy->occasional)
    build_one(mylist->listcopy,preflag);
  if (mylist->listfull && mylist->listfull->occasional)
    build_one(mylist->listfull,preflag);

  // create stencil if hasn't been created since last setup_bins() call

  NStencil *ns = np->ns;
  if (ns && ns->last_create < last_setup_bins) {
    ns->create_setup();
    ns->create();
  }



  // build the list

  np->build_setup();
  np->build(mylist);
}

/* ----------------------------------------------------------------------
   set neighbor style and skin distance
------------------------------------------------------------------------- */

void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command");

  skin = force->numeric(FLERR,arg[0]);
  if (skin < 0.0) error->all(FLERR,"Illegal neighbor command");

  if (strcmp(arg[1],"nsq") == 0) style = NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = BIN;
  else if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else error->all(FLERR,"Illegal neighbor command");

  if (style == MULTI && lmp->citeme) lmp->citeme->add(cite_neigh_multi);
}

/* ----------------------------------------------------------------------
   reset timestamps in all NeignBin, NStencil, NPair classes
   so that neighbor lists will rebuild properly with timestep change
------------------------------------------------------------------------- */

void Neighbor::reset_timestep(bigint ntimestep)
{
  for (int i = 0; i < nbin; i++) {
    neigh_bin[i]->last_setup = -1;
    neigh_bin[i]->last_bin = -1;
    neigh_bin[i]->last_bin_memory = -1;
  }

  for (int i = 0; i < nstencil; i++) {
    neigh_stencil[i]->last_create = -1;
    neigh_stencil[i]->last_stencil_memory = -1;
    neigh_stencil[i]->last_copy_bin = -1;
  }

  for (int i = 0; i < nlist; i++) {
    if (!neigh_pair[i]) continue;
    neigh_pair[i]->last_build = -1;
    neigh_pair[i]->last_copy_bin_setup = -1;
    neigh_pair[i]->last_copy_bin = -1;
    neigh_pair[i]->last_copy_stencil = -1;
  }
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
      every = force->inumeric(FLERR,arg[iarg+1]);
      if (every <= 0) error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      delay = force->inumeric(FLERR,arg[iarg+1]);
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
      old_pgsize = pgsize;
      pgsize = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      old_oneatom = oneatom;
      oneatom = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      binsize_user = force->numeric(FLERR,arg[iarg+1]);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cluster") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) cluster_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) cluster_check = 0;
      else error->all(FLERR,"Illegal neigh_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"include") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal neigh_modify command");
      includegroup = group->find(arg[iarg+1]);
      if (includegroup < 0)
        error->all(FLERR,"Invalid group ID in neigh_modify command");
      if (includegroup && (atom->firstgroupname == NULL ||
                            strcmp(arg[iarg+1],atom->firstgroupname) != 0))
        error->all(FLERR,
                   "Neigh_modify include group != atom_modify first group");
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
        ex1_type[nex_type] = force->inumeric(FLERR,arg[iarg+2]);
        ex2_type[nex_type] = force->inumeric(FLERR,arg[iarg+3]);
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
   remove the first group-group exclusion matching group1, group2
------------------------------------------------------------------------- */

void Neighbor::exclusion_group_group_delete(int group1, int group2)
{
  int m, mlast;
  for (m = 0; m < nex_group; m++)
    if (ex1_group[m] == group1 && ex2_group[m] == group2 )
      break;

  mlast = m;
  if (mlast == nex_group) 
    error->all(FLERR,"Unable to find group-group exclusion");
  
  for (m = mlast+1; m < nex_group; m++) {
    ex1_group[m-1] = ex1_group[m];
    ex2_group[m-1] = ex2_group[m];
    ex1_bit[m-1] = ex1_bit[m];
    ex2_bit[m-1] = ex2_bit[m];
  }
  nex_group--;
}


/* ----------------------------------------------------------------------
   return the value of exclude - used to check compatibility with GPU
------------------------------------------------------------------------- */

int Neighbor::exclude_setting()
{
  return exclude;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Neighbor::memory_usage()
{
  bigint bytes = 0;
  bytes += memory->usage(xhold,maxhold,3);

  for (int i = 0; i < nlist; i++)
    if (lists[i]) bytes += lists[i]->memory_usage();
  for (int i = 0; i < nstencil; i++)
    bytes += neigh_stencil[i]->memory_usage();
  for (int i = 0; i < nbin; i++)
    bytes += neigh_bin[i]->memory_usage();

  if (neigh_bond) bytes += neigh_bond->memory_usage();
  if (neigh_angle) bytes += neigh_angle->memory_usage();
  if (neigh_dihedral) bytes += neigh_dihedral->memory_usage();
  if (neigh_improper) bytes += neigh_improper->memory_usage();

  return bytes;
}
