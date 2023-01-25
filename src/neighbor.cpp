// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (triclinic and multi-neigh) : Pieter in 't Veld (SNL)
   Contributing author (improved multi-neigh) : Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "neighbor.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "nbin.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "npair.h"
#include "nstencil.h"
#include "ntopo.h"
#include "output.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "respa.h"
#include "style_nbin.h"  // IWYU pragma: keep
#include "style_npair.h"  // IWYU pragma: keep
#include "style_nstencil.h"  // IWYU pragma: keep
#include "style_ntopo.h"  // IWYU pragma: keep
#include "suffix.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace NeighConst;

#define RQDELTA 1
#define EXDELTA 1
#define DELTA_PERATOM 64

#define BIG 1.0e20

enum{NONE,ALL,PARTIAL,TEMPLATE};

static const char cite_neigh_multi_old[] =
  "neighbor multi/old command: doi:10.1016/j.cpc.2008.03.005\n\n"
  "@Article{Intveld08,\n"
  " author =  {in 't Veld, P. J. and S. J. Plimpton and G. S. Grest},\n"
  " title =   {Accurate and Efficient Methods for Modeling Colloidal\n"
  "            Mixtures in an Explicit Solvent using Molecular Dynamics},\n"
  " journal = {Comput.\\ Phys.\\ Commun.},\n"
  " year =    2008,\n"
  " volume =  179,\n"
  " number =  5,\n"
  " pages =   {320--329}\n"
  "}\n\n";

static const char cite_neigh_multi[] =
  "neighbor multi command: doi:10.1016/j.cpc.2008.03.005, doi:10.1007/s40571-020-00361-2\n\n"
  "@Article{Intveld08,\n"
  " author =  {in 't Veld, P. J. and S. J.~Plimpton and G. S. Grest},\n"
  " title =   {Accurate and Efficient Methods for Modeling Colloidal\n"
  "            Mixtures in an Explicit Solvent using Molecular Dynamics},\n"
  " journal = {Comput.\\ Phys.\\ Commut.},\n"
  " year =    2008,\n"
  " volume =  179,\n"
  " pages =   {320--329}\n"
  "}\n\n"
  "@article{Shire2020,\n"
  " author = {Shire, Tom and Hanley, Kevin J. and Stratford, Kevin},\n"
  " title = {{DEM} Simulations of Polydisperse Media: Efficient Contact\n"
  "          Detection Applied to Investigate the Quasi-Static Limit},\n"
  " journal = {Computational Particle Mechanics},\n"
  " year = {2020}\n"
  "}\n\n";

// template for factory functions:
// there will be one instance for each style keyword in the respective style_xxx.h files

template <typename S, typename T> static S *style_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

//#define NEIGH_LIST_DEBUG 1

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor(LAMMPS *lmp) : Pointers(lmp),
pairclass(nullptr), pairnames(nullptr), pairmasks(nullptr)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  firsttime = 1;

  style = Neighbor::BIN;
  every = 1;
  delay = 0;
  dist_check = 1;
  pgsize = 100000;
  oneatom = 2000;
  binsizeflag = 0;
  build_once = 0;
  cluster_check = 0;
  ago = -1;

  cutneighmax = 0.0;
  cutneighsq = nullptr;
  cutneighghostsq = nullptr;
  cuttype = nullptr;
  cuttypesq = nullptr;
  fixchecklist = nullptr;

  // pairwise neighbor lists and associated data structs

  nlist = 0;
  lists = nullptr;

  nbin = 0;
  neigh_bin = nullptr;

  nstencil = 0;
  neigh_stencil = nullptr;

  neigh_pair = nullptr;

  nstencil_perpetual = 0;
  slist = nullptr;

  npair_perpetual = 0;
  plist = nullptr;

  nrequest = maxrequest = 0;
  requests = nullptr;
  j_sorted = nullptr;

  old_nrequest = 0;
  old_requests = nullptr;

  old_style = style;
  old_triclinic = 0;
  old_pgsize = pgsize;
  old_oneatom = oneatom;

  binclass = nullptr;
  binnames = nullptr;
  binmasks = nullptr;
  stencilclass = nullptr;
  stencilnames = nullptr;
  stencilmasks = nullptr;

  // topology lists

  bondwhich = anglewhich = dihedralwhich = improperwhich = NONE;

  neigh_bond = nullptr;
  neigh_angle = nullptr;
  neigh_dihedral = nullptr;
  neigh_improper = nullptr;

  // coords at last neighboring

  maxhold = 0;
  xhold = nullptr;
  lastcall = -1;
  last_setup_bins = -1;

  // pair exclusion list info

  includegroup = 0;

  nex_type = maxex_type = 0;
  ex1_type = ex2_type = nullptr;
  ex_type = nullptr;

  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = nullptr;

  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = ex_mol_intra = nullptr;

  // Multi data

  type2collection = nullptr;
  collection2cut = nullptr;
  collection = nullptr;
  cutcollectionsq = nullptr;
  custom_collection_flag = 0;
  interval_collection_flag = 0;
  nmax_collection = 0;

  // Kokkos setting

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{
  if (copymode) return;

  memory->destroy(cutneighsq);
  memory->destroy(cutneighghostsq);
  delete[] cuttype;
  delete[] cuttypesq;
  delete[] fixchecklist;

  for (int i = 0; i < nlist; i++) delete lists[i];
  for (int i = 0; i < nbin; i++) delete neigh_bin[i];
  for (int i = 0; i < nstencil; i++) delete neigh_stencil[i];
  for (int i = 0; i < nlist; i++) delete neigh_pair[i];
  delete[] lists;
  delete[] neigh_bin;
  delete[] neigh_stencil;
  delete[] neigh_pair;

  delete[] slist;
  delete[] plist;

  for (int i = 0; i < nrequest; i++)
    if (requests[i]) delete requests[i];
  memory->sfree(requests);
  for (int i = 0; i < old_nrequest; i++)
    if (old_requests[i]) delete old_requests[i];
  memory->sfree(old_requests);

  delete[] j_sorted;

  delete[] binclass;
  delete[] binnames;
  delete[] binmasks;
  delete[] stencilclass;
  delete[] stencilnames;
  delete[] stencilmasks;
  delete[] pairclass;
  delete[] pairnames;
  delete[] pairmasks;

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
  delete[] ex1_bit;
  delete[] ex2_bit;

  memory->destroy(ex_mol_group);
  delete[] ex_mol_bit;
  memory->destroy(ex_mol_intra);

  memory->destroy(type2collection);
  memory->destroy(collection2cut);
  memory->destroy(collection);
  memory->destroy(cutcollectionsq);
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
  //   even if pair = nullptr and no neighbor lists are used
  // cutneigh = force cutoff + skin if cutforce > 0, else cutneigh = 0
  // cutneighghost = pair cutghost if it requests it, else same as cutneigh

  triggersq = 0.25*skin*skin;
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
                             (dimension == 3 && domain->zperiodic)))
      boxcheck = 1;

  n = atom->ntypes;
  if (cutneighsq == nullptr) {
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

  // Define cutoffs for multi
  if (style == Neighbor::MULTI) {
    int icollection, jcollection;

    // If collections not yet defined, create default map using types
    if (!custom_collection_flag) {
      ncollections = n;
      interval_collection_flag = 0;
      if (!type2collection)
        memory->create(type2collection,n+1,"neigh:type2collection");
      for (i = 1; i <= n; i++)
        type2collection[i] = i-1;
    }

    memory->grow(cutcollectionsq, ncollections, ncollections, "neigh:cutcollectionsq");

    // 3 possible ways of defining collections
    // 1) Types are used to define collections
    //    Each collection loops through its owned types, and uses cutneighsq to calculate its cutoff
    // 2) Collections are defined by intervals, point particles
    //    Types are first sorted into collections based on cutneighsq[i][i]
    //    Each collection loops through its owned types, and uses cutneighsq to calculate its cutoff
    // 3) Collections are defined by intervals, finite particles
    //

    // Define collection cutoffs
    for (i = 0; i < ncollections; i++)
      for (j = 0; j < ncollections; j++)
        cutcollectionsq[i][j] = 0.0;

    if (!interval_collection_flag) {
      finite_cut_flag = 0;
      for (i = 1; i <= n; i++){
        icollection = type2collection[i];
        for (j = 1; j <= n; j++){
          jcollection = type2collection[j];
          if (cutneighsq[i][j] > cutcollectionsq[icollection][jcollection]) {
            cutcollectionsq[icollection][jcollection] = cutneighsq[i][j];
            cutcollectionsq[jcollection][icollection] = cutneighsq[i][j];
          }
        }
      }
    } else {
      if (force->pair->finitecutflag) {
        finite_cut_flag = 1;
        // If cutoffs depend on finite atom sizes, use radii of intervals to find cutoffs
        double ri, rj, tmp;
        for (i = 0; i < ncollections; i++){
          ri = collection2cut[i]*0.5;
          for (j = 0; j < ncollections; j++){
            rj = collection2cut[j]*0.5;
            tmp = force->pair->radii2cut(ri, rj) + skin;
            cutcollectionsq[i][j] = tmp*tmp;
          }
        }
      } else {
        finite_cut_flag = 0;

        // Map types to collections
        if (!type2collection)
          memory->create(type2collection,n+1,"neigh:type2collection");

        for (i = 1; i <= n; i++)
          type2collection[i] = -1;

        double cuttmp;
        for (i = 1; i <= n; i++){
          // Remove skin added to cutneighsq
          cuttmp = sqrt(cutneighsq[i][i]) - skin;
          for (icollection = 0; icollection < ncollections; icollection ++){
            if (collection2cut[icollection] >= cuttmp) {
              type2collection[i] = icollection;
              break;
            }
          }

          if (type2collection[i] == -1)
            error->all(FLERR, "Pair cutoff exceeds interval cutoffs for multi");
        }

        // Define cutoffs
        for (i = 1; i <= n; i++){
          icollection = type2collection[i];
          for (j = 1; j <= n; j++){
            jcollection = type2collection[j];
            if (cutneighsq[i][j] > cutcollectionsq[icollection][jcollection]) {
              cutcollectionsq[icollection][jcollection] = cutneighsq[i][j];
              cutcollectionsq[jcollection][icollection] = cutneighsq[i][j];
            }
          }
        }
      }
    }
  }

  // rRESPA cutoffs

  int respa = 0;
  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    if ((dynamic_cast<Respa *>(update->integrate))->level_inner >= 0) respa = 1;
    if ((dynamic_cast<Respa *>(update->integrate))->level_middle >= 0) respa = 2;
  }

  if (respa) {
    double *cut_respa = (dynamic_cast<Respa *>(update->integrate))->cutoff;
    cut_inner_sq = (cut_respa[1] + skin) * (cut_respa[1] + skin);
    cut_middle_sq = (cut_respa[3] + skin) * (cut_respa[3] + skin);
    cut_middle_inside_sq = (cut_respa[0] - skin) * (cut_respa[0] - skin);
    if (cut_respa[0]-skin < 0) cut_middle_inside_sq = 0.0;
  }

  // fixchecklist = other classes that can induce reneighboring in decide()

  restart_check = 0;
  if (output->restart_flag) restart_check = 1;

  delete[] fixchecklist;
  fixchecklist = nullptr;
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
  // b/c pairwise portion of KSpace solver uses all 1-2,1-3,1-4 neighbors
  // some Coulomb-approximation pair styles also require it (below)

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

  // cannot remove special neighbors with kspace or kspace-like pair styles
  //   b/c exclusion needs to remove the full coulomb and not the damped interaction
  // special treatment required for hybrid pair styles since Force::pair_match()
  //   will only return a non-NULL pointer if there is only one substyle of the kind

  if (force->kspace) {
    special_flag[1] = special_flag[2] = special_flag[3] = 2;
  } else {
    PairHybrid *ph = reinterpret_cast<PairHybrid *>(force->pair_match("^hybrid",0));
    if (ph) {
      int flag=0;
      for (int isub=0; isub < ph->nstyles; ++isub) {
        if (force->pair_match("amoeba",0,isub)
            || force->pair_match("hippo",0,isub)
            || force->pair_match("coul/wolf",0,isub)
            || force->pair_match("coul/dsf",0,isub)
            || force->pair_match("coul/exclude",0)
            || force->pair_match("thole",0,isub))
          ++flag;
      }
      if (flag)
        special_flag[1] = special_flag[2] = special_flag[3] = 2;
    } else {
      if (force->pair_match("amoeba",0)
          || force->pair_match("hippo",0)
          || force->pair_match("coul/wolf",0)
          || force->pair_match("coul/dsf",0)
          || force->pair_match("coul/exclude",0)
          || force->pair_match("thole",0))
        special_flag[1] = special_flag[2] = special_flag[3] = 2;
    }
  }

  // ------------------------------------------------------------------
  // xhold array

  // free if not needed for this run

  if (dist_check == 0) {
    memory->destroy(xhold);
    maxhold = 0;
    xhold = nullptr;
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
      delete[] ex1_bit;
      delete[] ex2_bit;
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
      delete[] ex_mol_bit;
      ex_mol_bit = new int[nex_mol];
    }

    for (i = 0; i < nex_mol; i++)
      ex_mol_bit[i] = group->bitmask[ex_mol_group[i]];
  }

  if (exclude && force->kspace && me == 0)
    error->warning(FLERR,"Neighbor exclusions used with KSpace solver "
                   "may give inconsistent Coulombic energies");

  if (lmp->kokkos)
    set_binsize_kokkos();

  // ------------------------------------------------------------------
  // create pairwise lists
  // one-time call to init_styles() to scan style files and setup
  // init_pair() creates auxiliary classes: NBin, NStencil, NPair

  if (firsttime) init_styles();
  firsttime = 0;

  int same = init_pair();

  // invoke copy_neighbor_info() in Bin,Stencil,Pair classes
  // copied once per run in case any cutoff, exclusion, special info changed

  for (i = 0; i < nbin; i++) neigh_bin[i]->copy_neighbor_info();
  for (i = 0; i < nstencil; i++) neigh_stencil[i]->copy_neighbor_info();
  for (i = 0; i < nlist; i++)
    if (neigh_pair[i]) neigh_pair[i]->copy_neighbor_info();

  if (!same && (nrequest > 0) && (comm->me == 0)) print_pairwise_info();

  // can now delete requests so next run can make new ones
  // print_pairwise_info() made use of requests
  // set of NeighLists now stores all needed info

  for (i = 0; i < nrequest; i++) {
    delete requests[i];
    requests[i] = nullptr;
  }
  nrequest = 0;

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
#include "style_nbin.h"  // IWYU pragma: keep
#undef NBinStyle
#undef NBIN_CLASS

  binclass = new BinCreator[nbclass];
  binnames = new char*[nbclass];
  binmasks = new int[nbclass];
  nbclass = 0;

#define NBIN_CLASS
#define NBinStyle(key,Class,bitmasks) \
  binnames[nbclass] = (char *) #key; \
  binclass[nbclass] = &style_creator<NBin, Class>;   \
  binmasks[nbclass++] = bitmasks;
#include "style_nbin.h"  // IWYU pragma: keep
#undef NBinStyle
#undef NBIN_CLASS

  // extract info from NStencil classes listed in style_nstencil.h

  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) nsclass++;
#include "style_nstencil.h"  // IWYU pragma: keep
#undef NStencilStyle
#undef NSTENCIL_CLASS

  stencilclass = new StencilCreator[nsclass];
  stencilnames = new char*[nsclass];
  stencilmasks = new int[nsclass];
  nsclass = 0;

#define NSTENCIL_CLASS
#define NStencilStyle(key,Class,bitmasks) \
  stencilnames[nsclass] = (char *) #key; \
  stencilclass[nsclass] = &style_creator<NStencil, Class>;   \
  stencilmasks[nsclass++] = bitmasks;
#include "style_nstencil.h"  // IWYU pragma: keep
#undef NStencilStyle
#undef NSTENCIL_CLASS

  // extract info from NPair classes listed in style_npair.h

  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) npclass++;
#include "style_npair.h"  // IWYU pragma: keep
#undef NPairStyle
#undef NPAIR_CLASS

  pairclass = new PairCreator[npclass];
  pairnames = new char*[npclass];
  pairmasks = new int[npclass];
  npclass = 0;

#define NPAIR_CLASS
#define NPairStyle(key,Class,bitmasks) \
  pairnames[npclass] = (char *) #key; \
  pairclass[npclass] = &style_creator<NPair, Class>;  \
  pairmasks[npclass++] = bitmasks;
#include "style_npair.h"  // IWYU pragma: keep
#undef NPairStyle
#undef NPAIR_CLASS
}

/* ----------------------------------------------------------------------
   create and initialize NPair classes
------------------------------------------------------------------------- */

int Neighbor::init_pair()
{
  int i,j,k,m;

  // test if pairwise lists need to be re-created
  // no need to re-create if:
  //   neigh style, triclinic, pgsize, oneatom have not changed
  //   current requests = old requests
  // so just return:
  //   delete requests so next run can make new ones
  //   current set of NeighLists already stores all needed info
  // requests are compared via identical() before:
  //   any requests are morphed using logic below
  //   any requests are added below, e.g. as parents of pair hybrid skip lists
  // copy them via requests_new2old() BEFORE any changes made to requests
  //   necessary b/c morphs can change requestor settings (see comment below)

  int same = 1;
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

  if (same) return same;
  requests_new2old();

  // delete old lists since creating new ones

  for (i = 0; i < nlist; i++) delete lists[i];
  for (i = 0; i < nbin; i++) delete neigh_bin[i];
  for (i = 0; i < nstencil; i++) delete neigh_stencil[i];
  for (i = 0; i < nlist; i++) delete neigh_pair[i];
  delete[] lists;
  delete[] neigh_bin;
  delete[] neigh_stencil;
  delete[] neigh_pair;

  // error check on requests
  // do not allow occasional, ghost, bin list
  //   b/c it still uses variant of coord2bin() in NPair() method
  //     instead of atom2bin, this could cause error b/c stoms have
  //     moved out of proc domain by time occasional list is built
  //   solution would be to use a different NBin variant
  //     that used Npair::coord2bin(x,ix,iy,iz) (then delete it from NPair)
  //     and stored the ix,iy,iz values for all atoms (including ghosts)
  //     at time of binning when neighbor lists are rebuilt,
  //     similar to what vanilla Nbin::coord2atom() does now in atom2bin

  if (style == Neighbor::BIN) {
    for (i = 0; i < nrequest; i++)
      if (requests[i]->occasional && requests[i]->ghost)
        error->all(FLERR,"Cannot request an occasional binned neighbor list "
                   "with ghost info");
  }

  // morph requests in various ways
  // purpose is to avoid duplicate or inefficient builds
  // may add new requests if a needed request to derive from does not exist
  // methods:
  //   (1) unique = create unique lists if cutoff is explicitly set
  //   (2) skip = create any new non-skip lists needed by pair hybrid skip lists
  //   (3) granular = adjust parent and skip lists for granular onesided usage
  //   (4) h/f = pair up any matching half/full lists
  //   (5) copy = convert as many lists as possible to copy lists
  // order of morph methods matters:
  //   (3) after (2), b/c it adjusts lists created by (2)
  //   (4) after (2) and (3),
  //       b/c (2) may create new full lists, (3) may change them
  //   (5) last, after all lists are finalized, so all possible copies/trims found

  int nrequest_original = nrequest;

  morph_unique();
  morph_skip();
  morph_granular();     // this method can change flags set by requestor

  // sort requests by cutoff distance for trimming, used by
  //  morph_halffull and morph_copy_trim. Must come after
  //  morph_skip() which change the number of requests

  sort_requests();

  morph_halffull();
  morph_copy_trim();

  // create new lists, one per request including added requests
  // wait to allocate initial pages until copy lists are detected
  // NOTE: can I allocate now, instead of down below?

  nlist = nrequest;

  lists = new NeighList*[nrequest];
  neigh_bin = new NBin*[nrequest];
  neigh_stencil = new NStencil*[nrequest];
  neigh_pair = new NPair*[nrequest];

  // allocate new lists
  // pass list ptr back to requestor (except for Command class)
  // only for original requests, not ones added by Neighbor class

  for (i = 0; i < nrequest; i++) {
    if (requests[i]->kokkos_host || requests[i]->kokkos_device)
      create_kokkos_list(i);
    else lists[i] = new NeighList(lmp);
    lists[i]->index = i;
    lists[i]->requestor = requests[i]->requestor;

    if (requests[i]->pair) {
        lists[i]->requestor_type = NeighList::PAIR;
    } else if (requests[i]->fix) {
        lists[i]->requestor_type = NeighList::FIX;
    } else if (requests[i]->compute) {
        lists[i]->requestor_type = NeighList::COMPUTE;
    }

    if (requests[i]->pair && i < nrequest_original) {
      auto pair = (Pair *) requests[i]->requestor;
      pair->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->fix && i < nrequest_original) {
      Fix *fix = (Fix *) requests[i]->requestor;
      fix->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->compute && i < nrequest_original) {
      auto compute = (Compute *) requests[i]->requestor;
      compute->init_list(requests[i]->id,lists[i]);
    }
  }

  // invoke post_constructor() for all lists
  // copies info from requests to lists, sets ptrs to related lists

  for (i = 0; i < nrequest; i++)
    lists[i]->post_constructor(requests[i]);

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
  // unique = only one of its style, or request unique flag set (custom cutoff)

  nbin = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_bin = -1;
    flag = lists[i]->bin_method;
    if (flag == 0) continue;
    if (!requests[i]->unique) {
      for (j = 0; j < nbin; j++)
        if (neigh_bin[j]->istyle == flag &&
            neigh_bin[j]->cutoff_custom == 0.0) break;
      if (j < nbin) {
        requests[i]->index_bin = j;
        continue;
      }
    }

    BinCreator &bin_creator = binclass[flag-1];
    neigh_bin[nbin] = bin_creator(lmp);
    neigh_bin[nbin]->post_constructor(requests[i]);
    neigh_bin[nbin]->istyle = flag;

    requests[i]->index_bin = nbin;
    nbin++;
  }

  nstencil = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_stencil = -1;
    flag = lists[i]->stencil_method;
    if (flag == 0) continue;
    if (!requests[i]->unique) {
      for (j = 0; j < nstencil; j++)
        if (neigh_stencil[j]->istyle == flag &&
            neigh_stencil[j]->cutoff_custom == 0.0) break;
      if (j < nstencil) {
        requests[i]->index_stencil = j;
        continue;
      }
    }

    StencilCreator &stencil_creator = stencilclass[flag-1];
    neigh_stencil[nstencil] = stencil_creator(lmp);
    neigh_stencil[nstencil]->post_constructor(requests[i]);
    neigh_stencil[nstencil]->istyle = flag;

    if (lists[i]->bin_method > 0) {
      neigh_stencil[nstencil]->nb = neigh_bin[requests[i]->index_bin];
      if (neigh_stencil[nstencil]->nb == nullptr)
        error->all(FLERR,"Could not assign bin method to neighbor stencil");
    }

    requests[i]->index_stencil = nstencil;
    nstencil++;
  }

  // instantiate one Pair class per list in neigh_pair vec

  for (i = 0; i < nrequest; i++) {
    requests[i]->index_pair = -1;
    flag = lists[i]->pair_method;
    if (flag == 0) {
      neigh_pair[i] = nullptr;
      continue;
    }

    PairCreator &pair_creator = pairclass[flag-1];
    lists[i]->np = neigh_pair[i] = pair_creator(lmp);
    neigh_pair[i]->post_constructor(requests[i]);
    neigh_pair[i]->istyle = flag;

    if (lists[i]->bin_method > 0) {
      neigh_pair[i]->nb = neigh_bin[requests[i]->index_bin];
      if (neigh_pair[i]->nb == nullptr)
        error->all(FLERR,"Could not assign bin method to neighbor pair");
    }
    if (lists[i]->stencil_method > 0) {
      neigh_pair[i]->ns = neigh_stencil[requests[i]->index_stencil];
      if (neigh_pair[i]->ns == nullptr)
        error->all(FLERR,"Could not assign stencil method to neighbor pair");
    }

    requests[i]->index_pair = i;
  }

  // allocate initial pages for each list, except if copy flag set

  for (i = 0; i < nlist; i++) {
    if (lists[i]->copy && !lists[i]->trim && !lists[i]->kk2cpu)
      continue;
    lists[i]->setup_pages(pgsize,oneatom);
  }

  // first-time allocation of per-atom data for lists that are built and store
  // lists that do not store: copy
  // use atom->nmax for both grow() args
  //   i.e. grow first time to expanded size to avoid future reallocs
  // also Kokkos list initialization

  int maxatom = atom->nmax;
  for (i = 0; i < nlist; i++) {
    if (neigh_pair[i] && (!lists[i]->copy || lists[i]->trim || lists[i]->kk2cpu))
      lists[i]->grow(maxatom,maxatom);
  }

  // plist = indices of perpetual NPair classes
  //         perpetual = non-occasional, re-built at every reneighboring
  // slist = indices of perpetual NStencil classes
  //         perpetual = used by any perpetual NPair class

  delete[] slist;
  delete[] plist;
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
      if (lists[plist[j]]->stencil_method == neigh_stencil[i]->istyle)
        flag = 1;
    if (flag) slist[nstencil_perpetual++] = i;
  }

  // reorder plist vector if necessary
  // relevant for lists that are derived from a parent list:
  //   half-full,copy,skip
  // the child index must appear in plist after the parent index
  // swap two indices within plist when dependency is mis-ordered
  // start double loop check again whenever a swap is made
  // done when entire double loop test results in no swaps

  NeighList *ptr;

  int done = 0;
  while (!done) {
    done = 1;
    for (i = 0; i < npair_perpetual; i++) {
      for (k = 0; k < 3; k++) {
        ptr = nullptr;
        if (k == 0) ptr = lists[plist[i]]->listcopy;
        if (k == 1) ptr = lists[plist[i]]->listskip;
        if (k == 2) ptr = lists[plist[i]]->listfull;
        if (ptr == nullptr) continue;
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
      if (!done) break;
    }
  }

  // debug output

#ifdef NEIGH_LIST_DEBUG
  for (i = 0; i < nrequest; i++) lists[i]->print_attributes();
#endif

  return same;
}

/* ----------------------------------------------------------------------
   sort NeighRequests by cutoff distance
    to find smallest list for trimming
------------------------------------------------------------------------- */

void Neighbor::sort_requests()
{
  NeighRequest *jrq;
  int i,j,jmin;
  double jcut;

  delete[] j_sorted;
  j_sorted = new int[nrequest];

  for (i = 0; i < nrequest; i++)
    j_sorted[i] = i;

  for (i = 0; i < nrequest; i++) {
    double cutoff_min = cutneighmax;
    jmin = i;

    for (j = i; j < nrequest-1; j++) {
      jrq = requests[j_sorted[j]];
      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;

      if (jcut <= cutoff_min) {
        cutoff_min = jcut;
        jmin = j;
      }
    }
    int tmp = j_sorted[i];
    j_sorted[i] = j_sorted[jmin];
    j_sorted[jmin] = tmp;
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests to set additional flags:
   custom cutoff lists and accelerator lists
------------------------------------------------------------------------- */

void Neighbor::morph_unique()
{
  NeighRequest *irq;

  for (int i = 0; i < nrequest; i++) {
    irq = requests[i];

    // if cut flag set by requestor and cutoff is different than default,
    //  set unique flag, otherwise unset cut flag
    // this forces Pair,Stencil,Bin styles to be instantiated separately
    // also add skin to cutoff of perpetual lists

    if (irq->cut) {
      if (!irq->occasional)
        irq->cutoff += skin;

      if (irq->cutoff != cutneighmax) {
        irq->unique = 1;
      } else {
        irq->cut = 0;
        irq->cutoff = 0.0;
      }
    }

    // avoid flagging a neighbor list as both INTEL and OPENMP

    if (irq->intel) irq->omp = 0;

    // avoid flagging a neighbor list as both KOKKOS and INTEL or OPENMP

    if (irq->kokkos_host || irq->kokkos_device) irq->omp = irq->intel = 0;
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests to process all skip lists
   look for a matching non-skip list
   if one exists, point at it via skiplist
   else make new parent via copy_request() and point at it
------------------------------------------------------------------------- */

void Neighbor::morph_skip()
{
  int i,j,inewton,jnewton;
  NeighRequest *irq,*jrq,*nrq;

  for (i = 0; i < nrequest; i++) {
    irq = requests[i];

    // only processing skip lists

    if (!irq->skip) continue;

    // these lists are created other ways, no need for skipping
    // halffull list and its full parent may both skip,
    //   but are checked to ensure matching skip info

    if (irq->halffull) continue;
    if (irq->copy) continue;

    // check all other lists

    for (j = 0; j < nrequest; j++) {
      if (i == j) continue;
      jrq = requests[j];

      // can only skip from a perpetual non-skip list

      if (jrq->occasional) continue;
      if (jrq->skip) continue;

      // both lists must be half, or both full

      if (irq->half != jrq->half) continue;
      if (irq->full != jrq->full) continue;

      // both lists must be newton on, or both newton off
      // IJ newton = 1 for newton on, 2 for newton off

      inewton = irq->newton;
      if (inewton == 0) inewton = force->newton_pair ? 1 : 2;
      jnewton = jrq->newton;
      if (jnewton == 0) jnewton = force->newton_pair ? 1 : 2;
      if (inewton != jnewton) continue;

      // these flags must be same,
      //   else 2 lists do not store same pairs
      //   or their data structures are different
      // this includes custom cutoff set by requestor
      // NOTE: need check for 2 Kokkos flags?

      if (irq->ghost != jrq->ghost) continue;
      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->omp != jrq->omp) continue;
      if (irq->intel != jrq->intel) continue;
      if (irq->kokkos_host != jrq->kokkos_host) continue;
      if (irq->kokkos_device != jrq->kokkos_device) continue;
      if (irq->ssa != jrq->ssa) continue;
      if (irq->cut != jrq->cut) continue;
      if (irq->cutoff != jrq->cutoff) continue;

      // 2 lists are a match

      break;
    }

    // if matching list exists, point to it
    // else create a new identical list except non-skip
    // for new list, set neigh = 1, skip = 0, no skip vec/array,
    //   copy unique flag (since copy_request() will not do it)
    // note: parents of skip lists do not have associated history
    //   b/c child skip lists have the associated history

    if (j < nrequest) irq->skiplist = j;
    else {
      int newrequest = request(this,-1);
      irq->skiplist = newrequest;

      nrq = requests[newrequest];
      nrq->copy_request(irq,0);
      nrq->pair = nrq->fix = nrq->compute = nrq->command = 0;
      nrq->neigh = 1;
      nrq->skip = 0;
      if (irq->unique) nrq->unique = 1;
    }
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests just added by morph_skip for hybrid granular
   adjust newton/oneside parent settings if children require onesided skipping
   also set children off2on flag if parent becomes a newton off list
   this is needed because line/gran and tri/gran pair styles
     require onesided neigh lists and system newton on,
     but parent list must be newton off to enable the onesided skipping
------------------------------------------------------------------------- */

void Neighbor::morph_granular()
{
  int i,j;
  NeighRequest *irq,*jrq;

  for (i = 0; i < nrequest; i++) {
    irq = requests[i];

    // only examine NeighRequests added by morph_skip()
    // only those with size attribute for granular systems

    if (!irq->neigh) continue;
    if (!irq->size) continue;

    // check children of this list

    int onesided = -1;
    for (j = 0; j < nrequest; j++) {
      jrq = requests[j];

      // only consider JRQ pair, size lists that skip from Irq list

      if (!jrq->pair) continue;
      if (!jrq->size) continue;
      if (!jrq->skip || jrq->skiplist != i) continue;

      // onesided = -1 if no children
      // onesided = 0/1 = child granonesided value if same for all children
      // onesided = 2 if children have different granonesided values

      if (onesided < 0) onesided = jrq->granonesided;
      else if (onesided != jrq->granonesided) onesided = 2;
      if (onesided == 2) break;
    }

    // if onesided = 2, parent has children with both granonesided = 0/1
    // force parent newton off (newton = 2) to enable onesided skip by child
    // set parent granonesided = 0, so it stores all neighs in usual manner
    // set off2on = 1 for all children, since they expect newton on lists
    //   this is b/c granonesided only set by line/gran and tri/gran which
    //   both require system newton on

    if (onesided == 2) {
      irq->newton = 2;
      irq->granonesided = 0;

      for (j = 0; j < nrequest; j++) {
        jrq = requests[j];

        // only consider JRQ pair, size lists that skip from Irq list

        if (!jrq->pair) continue;
        if (!jrq->size) continue;
        if (!jrq->skip || jrq->skiplist != i) continue;

        jrq->off2on = 1;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests for possible half lists to derive from full lists
   if 2 requests match, set half list to derive from full list
------------------------------------------------------------------------- */

void Neighbor::morph_halffull()
{
  int i,j,jj;
  NeighRequest *irq,*jrq;
  double icut,jcut;

  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;

    // only processing half lists

    if (!irq->half) continue;

    // these lists are created other ways, no need for halffull
    // do want to process skip lists

    if (irq->copy) continue;

    // check all other lists

    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut) j = j_sorted[jj];
      else j = jj;

      jrq = requests[j];

      // can only derive from a perpetual full list
      // newton setting of derived list does not matter

      if (jrq->occasional) continue;
      if (!jrq->full) continue;

      // trim a list with longer cutoff

      if (irq->cut) icut = irq->cutoff;
      else icut = cutneighmax;

      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;

      if (icut > jcut) continue;
      else if (icut != jcut) trim_flag = 1;

      // these flags must be same,
      //   else 2 lists do not store same pairs
      //   or their data structures are different

      if (irq->ghost != jrq->ghost) continue;
      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->omp != jrq->omp) continue;
      if (irq->intel != jrq->intel) continue;
      if (irq->kokkos_host != jrq->kokkos_host) continue;
      if (irq->kokkos_device != jrq->kokkos_device) continue;
      if (irq->ssa != jrq->ssa) continue;

      // skip flag must be same
      // if both are skip lists, skip info must match

      if (irq->skip != jrq->skip) continue;
      if (irq->skip && irq->same_skip(jrq) == 0) continue;

      // 2 lists are a match

      break;
    }

    // if matching list exists, point to it

    if (jj < nrequest) {
      irq->halffull = 1;
      irq->halffulllist = j;
      irq->trim = trim_flag;
    }
  }
}

/* ----------------------------------------------------------------------
   scan NeighRequests for possible copies or trims
   if 2 requests match, turn one into a copy or trim of the other
------------------------------------------------------------------------- */

void Neighbor::morph_copy_trim()
{
  int i,j,jj,inewton,jnewton;
  NeighRequest *irq,*jrq;
  double icut,jcut;

  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;

    // this list is already a copy list due to another morph method

    if (irq->copy) continue;

    // check all other lists

    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut) j = j_sorted[jj];
      else j = jj;

      if (i == j) continue;
      jrq = requests[j];

      // other list is already copied from this one

      if (jrq->copy && jrq->copylist == i) continue;

      // trim a list with longer cutoff

      if (irq->cut) icut = irq->cutoff;
      else icut = cutneighmax;

      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;

      if (icut > jcut) continue;
      else if (icut != jcut) trim_flag = 1;

      // other list (jrq) to copy from must be perpetual
      // list that becomes a copy list (irq) can be perpetual or occasional
      // if both lists are perpetual, require j < i
      //   to prevent circular dependence with 3 or more copies of a list

      if (jrq->occasional) continue;
      if (!irq->occasional && !irq->cut && j > i) continue;

      // both lists must be half, or both full

      if (irq->half != jrq->half) continue;
      if (irq->full != jrq->full) continue;

      // both lists must be newton on, or both newton off
      // IJ newton = 1 for newton on, 2 for newton off

      inewton = irq->newton;
      if (inewton == 0) inewton = force->newton_pair ? 1 : 2;
      jnewton = jrq->newton;
      if (jnewton == 0) jnewton = force->newton_pair ? 1 : 2;
      if (inewton != jnewton) continue;

      // ok for non-ghost list to copy from ghost list, but not vice versa

      if (irq->ghost && !jrq->ghost) continue;

      // do not copy from a list with respa middle/inner
      // b/c its outer list will not be complete

      if (jrq->respamiddle) continue;
      if (jrq->respainner) continue;

      // these flags must be same,
      //   else 2 lists do not store same pairs
      //   or their data structures are different
      // no need to check omp b/c it stores same pairs
      // NOTE: need check for 2 Kokkos flags?

      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->intel != jrq->intel) continue;
      if (irq->kokkos_host && !jrq->kokkos_host) continue;
      if (irq->kokkos_device && !jrq->kokkos_device) continue;
      if (irq->ssa != jrq->ssa) continue;

      // skip flag must be same
      // if both are skip lists, skip info must match

      if (irq->skip != jrq->skip) continue;
      if (irq->skip && irq->same_skip(jrq) == 0) continue;

      // 2 lists are a match

      break;
    }

    // turn list I into a copy of list J
    // do not copy a list from another copy list, but from its parent list

    if (jj < nrequest) {
      irq->copy = 1;
      irq->trim = trim_flag;
      if (jrq->copy && irq->cutoff == requests[jrq->copylist]->cutoff)
        irq->copylist = jrq->copylist;
      else
        irq->copylist = j;
    }
  }
}

/* ----------------------------------------------------------------------
   create and initialize NTopo classes
------------------------------------------------------------------------- */

void Neighbor::init_topology()
{
  int i,m;

  if (atom->molecular == Atom::ATOMIC) return;

  // set flags that determine which topology neighbor classes to use
  // these settings could change from run to run, depending on fixes defined
  // bonds,etc can only be broken for atom->molecular = Atom::MOLECULAR, not Atom::TEMPLATE
  // SHAKE sets bonds and angles negative
  // gcmc sets all bonds, angles, etc negative
  // partial_flag sets bonds to 0
  // delete_bonds sets all interactions negative

  int bond_off = 0;
  int angle_off = 0;
  for (i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^shake")
        || utils::strmatch(modify->fix[i]->style,"^rattle"))
      bond_off = angle_off = 1;
  if (force->bond)
    if (force->bond->partial_flag)
      bond_off = 1;

  if (atom->avec->bonds_allow && atom->molecular == Atom::MOLECULAR) {
    for (i = 0; i < atom->nlocal; i++) {
      if (bond_off) break;
      for (m = 0; m < atom->num_bond[i]; m++)
        if (atom->bond_type[i][m] <= 0) bond_off = 1;
    }
  }

  if (atom->avec->angles_allow && atom->molecular == Atom::MOLECULAR) {
    for (i = 0; i < atom->nlocal; i++) {
      if (angle_off) break;
      for (m = 0; m < atom->num_angle[i]; m++)
        if (atom->angle_type[i][m] <= 0) angle_off = 1;
    }
  }

  int dihedral_off = 0;
  if (atom->avec->dihedrals_allow && atom->molecular == Atom::MOLECULAR) {
    for (i = 0; i < atom->nlocal; i++) {
      if (dihedral_off) break;
      for (m = 0; m < atom->num_dihedral[i]; m++)
        if (atom->dihedral_type[i][m] <= 0) dihedral_off = 1;
    }
  }

  int improper_off = 0;
  if (atom->avec->impropers_allow && atom->molecular == Atom::MOLECULAR) {
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
    if (atom->molecular == Atom::TEMPLATE) bondwhich = TEMPLATE;
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
    if (atom->molecular == Atom::TEMPLATE) anglewhich = TEMPLATE;
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
    if (atom->molecular == Atom::TEMPLATE) dihedralwhich = TEMPLATE;
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
    if (atom->molecular == Atom::TEMPLATE) improperwhich = TEMPLATE;
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
  int i;
  NeighRequest *rq;

  const double cutghost = MAX(cutneighmax,comm->cutghostuser);

  double binsize, bbox[3];
  bbox[0] =  bboxhi[0]-bboxlo[0];
  bbox[1] =  bboxhi[1]-bboxlo[1];
  bbox[2] =  bboxhi[2]-bboxlo[2];
  if (binsizeflag) binsize = binsize_user;
  else if (style == Neighbor::BIN) binsize = 0.5*cutneighmax;
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

  std::string out = "Neighbor list info ...\n";
  out += fmt::format("  update: every = {} steps, delay = {} steps, check = {}\n",
                     every,delay,dist_check ? "yes" : "no");
  out += fmt::format("  max neighbors/atom: {}, page size: {}\n",
                     oneatom, pgsize);
  out += fmt::format("  master list distance cutoff = {:.8g}\n",cutneighmax);
  out += fmt::format("  ghost atom cutoff = {:.8g}\n",cutghost);
  if (style != Neighbor::NSQ)
    out += fmt::format("  binsize = {:.8g}, bins = {:g} {:g} {:g}\n",binsize,
                       ceil(bbox[0]/binsize), ceil(bbox[1]/binsize),
                       ceil(bbox[2]/binsize));

  out += fmt::format("  {} neighbor lists, perpetual/occasional/extra = {} {} {}\n",
                     nlist,nperpetual,noccasional,nextra);

  for (i = 0; i < nlist; i++) {
    rq = requests[i];
    if (rq->pair) {
      char *pname = force->pair_match_ptr((Pair *) rq->requestor);
      if (pname) out += fmt::format("  ({}) pair {}",i+1,pname);
      else out += fmt::format("  ({}) pair (none)",i+1);
    } else if (rq->fix) {
      out += fmt::format("  ({}) fix {}",i+1,((Fix *) rq->requestor)->style);
    } else if (rq->compute) {
      out += fmt::format("  ({}) compute {}",i+1,((Compute *) rq->requestor)->style);
    } else if (rq->command) {
      out += fmt::format("  ({}) command {}",i+1,rq->command_style);
    } else if (rq->neigh) {
      out += fmt::format("  ({}) neighbor class addition",i+1);
    }

    if (rq->occasional) out += ", occasional";
    else out += ", perpetual";

    // order these to get single output of most relevant

    if (rq->copy) {
      if (rq->trim)
        out += fmt::format(", trim from ({})",rq->copylist+1);
      else
        out += fmt::format(", copy from ({})",rq->copylist+1);
    } else if (rq->halffull)
      if (rq->trim)
        out += fmt::format(", half/full trim from ({})",rq->halffulllist+1);
      else
        out += fmt::format(", half/full from ({})",rq->halffulllist+1);
    else if (rq->skip)
      out += fmt::format(", skip from ({})",rq->skiplist+1);
    out += "\n";

    // list of neigh list attributes

    out += "      attributes: ";
    if (rq->half) out += "half";
    else if (rq->full) out += "full";

    if (rq->newton == 0) {
      if (force->newton_pair) out += ", newton on";
      else out += ", newton off";
    } else if (rq->newton == 1) out += ", newton on";
    else if (rq->newton == 2) out += ", newton off";

    if (rq->ghost) out += ", ghost";
    if (rq->size) out += ", size";
    if (rq->history) out += ", history";
    if (rq->granonesided) out += ", onesided";
    if (rq->respamiddle) out += ", respa outer/middle/inner";
    else if (rq->respainner) out += ", respa outer/inner";
    if (rq->bond) out += ", bond";
    if (rq->omp) out += ", omp";
    if (rq->intel) out += ", intel";
    if (rq->kokkos_device) out += ", kokkos_device";
    if (rq->kokkos_host) out += ", kokkos_host";
    if (rq->ssa) out += ", ssa";
    if (rq->cut) out += fmt::format(", cut {}",rq->cutoff);
    if (rq->off2on) out += ", off2on";
    out += "\n";

    out += "      ";
    if (lists[i]->pair_method == 0) out += "pair build: none\n";
    else out += fmt::format("pair build: {}\n",pairnames[lists[i]->pair_method-1]);

    out += "      ";
    if (lists[i]->stencil_method == 0) out += "stencil: none\n";
    else out += fmt::format("stencil: {}\n",stencilnames[lists[i]->stencil_method-1]);

    out += "      ";
    if (lists[i]->bin_method == 0) out += "bin: none\n";
    else out += fmt::format("bin: {}\n",binnames[lists[i]->bin_method-1]);
  }
  utils::logmesg(lmp,out);
}

/* ----------------------------------------------------------------------
   make copy of current requests and Neighbor params
   used to compare to when next run occurs
------------------------------------------------------------------------- */

void Neighbor::requests_new2old()
{
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);

  old_nrequest = nrequest;
  old_requests = (NeighRequest **)
    memory->smalloc(old_nrequest*sizeof(NeighRequest *),"neighbor:old_requests");

  for (int i = 0; i < old_nrequest; i++)
    old_requests[i] = new NeighRequest(requests[i]);

  old_style = style;
  old_triclinic = triclinic;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
}

/* ----------------------------------------------------------------------
   find and return request made by classptr
   if not found or classptr = nullptr, return nullptr
   id is optional and defaults to 0, which is the request id value unless set explicitly
------------------------------------------------------------------------- */

NeighRequest *Neighbor::find_request(void *classptr, const int id) const
{
  if (classptr == nullptr) return nullptr;

  for (int i = 0; i < nrequest; i++)
    if ((requests[i]->requestor == classptr) && (requests[i]->id == id)) return requests[i];

  return nullptr;
}

/* ----------------------------------------------------------------------
   return vector with neighbor list requests from pair styles
------------------------------------------------------------------------- */

const std::vector<NeighRequest *> Neighbor::get_pair_requests() const
{
  std::vector<NeighRequest *> matches;
  for (int i=0; i < nrequest; ++i)
    if (requests[i]->pair) matches.push_back(requests[i]);
  return matches;
}

/* ----------------------------------------------------------------------
   find and return list requested by classptr
   if not found or classptr = nullptr, return nullptr
   id is optional and defaults to 0, which is the request id value unless set explicitly
------------------------------------------------------------------------- */

NeighList *Neighbor::find_list(void *classptr, const int id) const
{
  if (classptr == nullptr) return nullptr;

  for (int i = 0; i < nlist; i++)
    if ((lists[i]->requestor == classptr) && (lists[i]->id == id)) return lists[i];

  return nullptr;
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

  if (style == Neighbor::NSQ) return 0;
  if (rq->skip || rq->copy || rq->halffull) return 0;

  // use request settings to match exactly one NBin class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < nbclass; i++) {
    mask = binmasks[i];

    // require match of these request flags and mask bits
    // (!A != !B) is effectively a logical xor

    if (!rq->intel != !(mask & NB_INTEL)) continue;
    if (!rq->ssa != !(mask & NB_SSA)) continue;
    if (!rq->kokkos_device != !(mask & NB_KOKKOS_DEVICE)) continue;
    if (!rq->kokkos_host != !(mask & NB_KOKKOS_HOST)) continue;

    // multi neighbor style require multi bin style
    if (style == Neighbor::MULTI) {
      if (!(mask & NB_MULTI)) continue;
    } else {
      if (!(mask & NB_STANDARD)) continue;
    }

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

  if (style == Neighbor::NSQ) return 0;
  if (rq->skip || rq->copy || rq->halffull) return 0;

  // convert newton request to newtflag = on or off

  int newtflag = 1;
  if (rq->newton == 0 && newton_pair) newtflag = 1;
  else if (rq->newton == 0 && !newton_pair) newtflag = 0;
  else if (rq->newton == 1) newtflag = 1;
  else if (rq->newton == 2) newtflag = 0;

  // request a full stencil if building full neighbor list or newton is off
  int fullflag = 0;
  if (rq->full) fullflag = 1;
  if (!newtflag) fullflag = 1;

  //printf("STENCIL RQ FLAGS: hff %d %d n %d g %d s %d newtflag %d fullflag %d\n",
  //       rq->half,rq->full,rq->newton,rq->ghost,rq->ssa,
  //       newtflag, fullflag);

  // use request and system settings to match exactly one NStencil class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < nsclass; i++) {
    mask = stencilmasks[i];

    //printf("III %d: half %d full %d ghost %d ssa %d\n",
    //       i,mask & NS_HALF,mask & NS_FULL,mask & NS_GHOST,mask & NS_SSA);

    // exactly one of half or full is set and must match

    if (fullflag) {
      if (!(mask & NS_FULL)) continue;
    } else {
      if (!(mask & NS_HALF)) continue;
    }

    // require match of these request flags and mask bits
    // (!A != !B) is effectively a logical xor

    if (!rq->ghost != !(mask & NS_GHOST)) continue;
    if (!rq->ssa != !(mask & NS_SSA)) continue;

    // neighbor style is one of BIN, MULTI_OLD, or MULTI and must match

    if (style == Neighbor::BIN) {
      if (!(mask & NS_BIN)) continue;
    } else if (style == Neighbor::MULTI_OLD) {
      if (!(mask & NS_MULTI_OLD)) continue;
    } else if (style == Neighbor::MULTI) {
      if (!(mask & NS_MULTI)) continue;
    }

    // dimension is 2 or 3 and must match

    if (dimension == 2) {
      if (!(mask & NS_2D)) continue;
    } else if (dimension == 3) {
      if (!(mask & NS_3D)) continue;
    }

    // domain triclinic flag is on or off and must match

    if (triclinic) {
      if (!(mask & NS_TRI)) continue;
    } else if (!triclinic) {
      if (!(mask & NS_ORTHO)) continue;
    }

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
  // error check for includegroup with ghost neighbor request

  if (includegroup && rq->ghost)
    error->all(FLERR,"Neighbor include group not allowed with ghost neighbors");

  // convert newton request to newtflag = on or off

  bool newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = true;
  else if (rq->newton == 0 && !newton_pair) newtflag = false;
  else if (rq->newton == 1) newtflag = true;
  else if (rq->newton == 2) newtflag = false;
  else error->all(FLERR,"Illegal 'newton' flag in neighbor list request");

  int molecular = atom->molecular;

  //printf("PAIR RQ FLAGS: hf %d %d n %d g %d sz %d gos %d r %d b %d o %d i %d "
  //       "kk %d %d ss %d dn %d sk %d cp %d hf %d oo %d\n",
  //        rq->half,rq->full,rq->newton,rq->ghost,rq->size,
  //        rq->granonesided,rq->respaouter,rq->bond,rq->omp,rq->intel,
  //       rq->kokkos_host,rq->kokkos_device,rq->ssa,rq->dnum,
  //      rq->skip,rq->copy,rq->halffull,rq->off2on);

  // use request and system settings to match exactly one NPair class mask
  // checks are bitwise using NeighConst bit masks

  int mask;

  for (int i = 0; i < npclass; i++) {
    mask = pairmasks[i];

    //printf("  PAIR NAMES i %d %d name %s mask %d\n",i,nrequest,
    //       pairnames[i],pairmasks[i]);

    // if copy request, no further checks needed, just return or continue
    // trim and Kokkos device/host flags must also match in order to copy
    // intel and omp flags must match to trim

    if (rq->copy) {
      if (!(mask & NP_COPY)) continue;
      if (rq->trim) {
        if (!rq->trim != !(mask & NP_TRIM)) continue;
        if (!rq->omp != !(mask & NP_OMP)) continue;
        if (!rq->intel != !(mask & NP_INTEL)) continue;
      }
      if (rq->kokkos_device || rq->kokkos_host) {
        if (!rq->kokkos_device != !(mask & NP_KOKKOS_DEVICE)) continue;
        if (!rq->kokkos_host != !(mask & NP_KOKKOS_HOST)) continue;
      }
      if (!requests[rq->copylist]->kokkos_device != !(mask & NP_KOKKOS_DEVICE)) continue;
      if (!requests[rq->copylist]->kokkos_host != !(mask & NP_KOKKOS_HOST)) continue;
      return i+1;
    }

    // exactly one of half or full is set and must match

    if (rq->half) {
      if (!(mask & NP_HALF)) continue;
    } else if (rq->full) {
      if (!(mask & NP_FULL)) continue;
    }

    // newtflag is on or off and must match

    if (newtflag) {
      if (!(mask & NP_NEWTON)) continue;
    } else if (!newtflag) {
      if (!(mask & NP_NEWTOFF)) continue;
    }

    // if molecular on, do not match ATOMONLY (b/c a MOLONLY Npair exists)
    // if molecular off, do not match MOLONLY (b/c an ATOMONLY Npair exists)

    if (molecular != Atom::ATOMIC) {
      if (mask & NP_ATOMONLY) continue;
    } else if (molecular == Atom::ATOMIC) {
      if (mask & NP_MOLONLY) continue;
    }

    // require match of these request flags and mask bits
    // (!A != !B) is effectively a logical xor

    if (!rq->ghost != !(mask & NP_GHOST)) continue;
    if (!rq->size != !(mask & NP_SIZE)) continue;
    if (!rq->respaouter != !(mask & NP_RESPA)) continue;
    if (!rq->granonesided != !(mask & NP_ONESIDE)) continue;
    if (!rq->respaouter != !(mask & NP_RESPA)) continue;
    if (!rq->bond != !(mask & NP_BOND)) continue;
    if (!rq->omp != !(mask & NP_OMP)) continue;
    if (!rq->intel != !(mask & NP_INTEL)) continue;
    if (!rq->kokkos_device != !(mask & NP_KOKKOS_DEVICE)) continue;
    if (!rq->kokkos_host != !(mask & NP_KOKKOS_HOST)) continue;
    if (!rq->ssa != !(mask & NP_SSA)) continue;

    if (!rq->skip != !(mask & NP_SKIP)) continue;

    if (!rq->trim != !(mask & NP_TRIM)) continue;

    if (!rq->halffull != !(mask & NP_HALF_FULL)) continue;
    if (!rq->off2on != !(mask & NP_OFF2ON)) continue;

    // neighbor style is one of NSQ, BIN, MULTI_OLD, or MULTI and must match

    if (style == Neighbor::NSQ) {
      if (!(mask & NP_NSQ)) continue;
    } else if (style == Neighbor::BIN) {
      if (!(mask & NP_BIN)) continue;
    } else if (style == Neighbor::MULTI_OLD) {
      if (!(mask & NP_MULTI_OLD)) continue;
    } else if (style == Neighbor::MULTI) {
      if (!(mask & NP_MULTI)) continue;
    }

    // domain triclinic flag is on or off and must match

    if (triclinic) {
      if (!(mask & NP_TRI)) continue;
    } else if (!triclinic) {
      if (!(mask & NP_ORTHO)) continue;
    }

    return i+1;
  }

  // error return if matched none

  return -1;
}

/* ----------------------------------------------------------------------
   called internally to request a pairwise neighbor list
------------------------------------------------------------------------- */

int Neighbor::request(void *requestor, int instance)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)
      memory->srealloc(requests,maxrequest*sizeof(NeighRequest *), "neighbor:requests");
  }

  requests[nrequest] = new NeighRequest(lmp, requestor, instance);
  nrequest++;
  return nrequest-1;
}

/* ----------------------------------------------------------------------
  add_request() is called by other classes to request a pairwise neighbor list
------------------------------------------------------------------------- */

NeighRequest *Neighbor::add_request(Pair *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->apply_flags(flags);
  // apply intel flag. omp flag is set globally via set_omp_neighbor()
  if (requestor->suffix_flag & Suffix::INTEL) {
    req->intel = 1;
    req->omp = 0;
  }
  return req;
}

NeighRequest *Neighbor::add_request(Fix *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->pair = 0;
  req->fix = 1;
  req->apply_flags(flags);
  return req;
}

NeighRequest *Neighbor::add_request(Compute *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->pair = 0;
  req->compute = 1;
  req->apply_flags(flags);
  return req;
}

NeighRequest *Neighbor::add_request(Command *requestor, const char *style, int flags)
{
  int irequest = request(requestor, 0);
  auto req = requests[irequest];
  req->pair = 0;
  req->command = 1;
  req->occasional = 1;
  req->command_style = style;
  req->apply_flags(flags);
  return req;
}

// set neighbor list request OpenMP flag

void Neighbor::set_omp_neighbor(int flag)
{
  // flag *all* neighbor list requests as OPENMP threaded,
  // but skip lists already flagged as INTEL threaded
  for (int i = 0; i < nrequest; ++i)
    if (!requests[i]->intel) requests[i]->omp = flag;
}

/* report if there is a neighbor list with the intel flag set */

bool Neighbor::has_intel_request() const
{
  return (((nrequest > 0) && (requests[0]->intel > 0))
          ||  ((old_nrequest > 0) && (old_requests[0]->intel > 0)));
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
     if reduced skin distance is negative, set to zero
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
      if (delta < 0.0) delta = 0.0;
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
      if (delta < 0.0) delta = 0.0;
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
  // rebuild collection array from scratch
  if (style == Neighbor::MULTI) build_collection(0);

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
  // not just NBin associated with perpetual lists, also occasional lists
  // b/c cannot wait to bin occasional lists in build_one() call
  // if bin then, atoms may have moved outside of proc domain & bin extent,
  //   leading to errors or even a crash

  if (style != Neighbor::NSQ) {
    if (last_setup_bins < 0) setup_bins();
    for (i = 0; i < nbin; i++) {
      neigh_bin[i]->bin_atoms_setup(nall);
      neigh_bin[i]->bin_atoms();
    }
  }

  // build pairwise lists for all perpetual NPair/NeighList
  // grow() with nlocal/nall args so that only realloc if have to

  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->copy || lists[m]->trim || lists[m]->kk2cpu)
      lists[m]->grow(nlocal,nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }

  // build topology lists for bonds/angles/etc

  if ((atom->molecular != Atom::ATOMIC) && topoflag) build_topology();
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

  if (mylist == nullptr)
    error->all(FLERR,"Trying to build an occasional neighbor list before initialization complete");

  // build_one() should never be invoked on a perpetual list

  if (!mylist->occasional) error->all(FLERR,"Neighbor::build_one() invoked on perpetual list");

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
  // or this is halffull and parent is occasional list,
  // or this is skip list and parent is occasional list,
  // ensure parent is current

  if (mylist->listcopy && mylist->listcopy->occasional)
    build_one(mylist->listcopy,preflag);
  if (mylist->listfull && mylist->listfull->occasional)
    build_one(mylist->listfull,preflag);
  if (mylist->listskip && mylist->listskip->occasional)
    build_one(mylist->listskip,preflag);

  // create stencil if hasn't been created since last setup_bins() call

  NStencil *ns = np->ns;
  if (ns && ns->last_stencil < last_setup_bins) {
    ns->create_setup();
    ns->create();
  }

  // build the list

  if (!mylist->copy || mylist->trim || mylist->kk2cpu)
    mylist->grow(atom->nlocal,atom->nlocal+atom->nghost);
  np->build_setup();
  np->build(mylist);
}

/* ----------------------------------------------------------------------
   set neighbor style and skin distance
------------------------------------------------------------------------- */

void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command: expected 2 arguments but found {}", narg);

  skin = utils::numeric(FLERR,arg[0],false,lmp);
  if (skin < 0.0) error->all(FLERR, "Invalid neighbor argument: {}", arg[0]);

  if (strcmp(arg[1],"nsq") == 0) style = Neighbor::NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = Neighbor::BIN;
  else if (strcmp(arg[1],"multi") == 0) {
    style = Neighbor::MULTI;
    ncollections = atom->ntypes;
  } else if (strcmp(arg[1],"multi/old") == 0) style = Neighbor::MULTI_OLD;
  else error->all(FLERR,"Unknown neighbor {} argument: {}", arg[0], arg[1]);

  if (style == Neighbor::MULTI_OLD && lmp->citeme) lmp->citeme->add(cite_neigh_multi_old);
  if (style == Neighbor::MULTI && lmp->citeme) lmp->citeme->add(cite_neigh_multi);
}

/* ----------------------------------------------------------------------
   reset timestamps in all NeignBin, NStencil, NPair classes
   so that neighbor lists will rebuild properly with timestep change
   ditto for lastcall and last_setup_bins
------------------------------------------------------------------------- */

void Neighbor::reset_timestep(bigint /*ntimestep*/)
{
  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->last_bin = -1;
  for (int i = 0; i < nstencil; i++)
    neigh_stencil[i]->last_stencil = -1;
  for (int i = 0; i < nlist; i++) {
    if (!neigh_pair[i]) continue;
    neigh_pair[i]->last_build = -1;
  }

  lastcall = -1;
  last_setup_bins = -1;
}

/* ----------------------------------------------------------------------
   modify parameters of the pair-wise neighbor build
------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify every", error);
      every = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (every <= 0) error->all(FLERR, "Invalid neigh_modify every argument: {}", every);
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify delay", error);
      delay = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (delay < 0) error->all(FLERR, "Invalid neigh_modify delay argument: {}", delay);
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify check", error);
      dist_check = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"once") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify once", error);
      build_once = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify page", error);
      old_pgsize = pgsize;
      pgsize = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify one", error);
      old_oneatom = oneatom;
      oneatom = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify binsize", error);
      binsize_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cluster") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify cluster", error);
      cluster_check = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"include") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify include", error);
      includegroup = group->find(arg[iarg+1]);
      if (includegroup < 0)
        error->all(FLERR, "Invalid include keyword: group {} not found", arg[iarg+1]);
      if (atom->firstgroupname == nullptr)
          error->all(FLERR, "Invalid include keyword: atom_modify first command must be used");
      if (strcmp(arg[iarg+1],atom->firstgroupname) != 0)
        error->all(FLERR, "Neigh_modify include group != atom_modify first group: {}", atom->firstgroupname);
      iarg += 2;
    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude", error);

      if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude type", error);
        if (nex_type == maxex_type) {
          maxex_type += EXDELTA;
          memory->grow(ex1_type,maxex_type,"neigh:ex1_type");
          memory->grow(ex2_type,maxex_type,"neigh:ex2_type");
        }
        ex1_type[nex_type] = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        ex2_type[nex_type] = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
        nex_type++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"group") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude group", error);
        if (nex_group == maxex_group) {
          maxex_group += EXDELTA;
          memory->grow(ex1_group,maxex_group,"neigh:ex1_group");
          memory->grow(ex2_group,maxex_group,"neigh:ex2_group");
        }
        ex1_group[nex_group] = group->find(arg[iarg+2]);
        ex2_group[nex_group] = group->find(arg[iarg+3]);
        if (ex1_group[nex_group] == -1)
          error->all(FLERR, "Invalid exclude group keyword: group {} not found", arg[iarg+2]);
        if (ex2_group[nex_group] == -1)
            error->all(FLERR, "Invalid exclude group keyword: group {} not found", arg[iarg+3]);
        nex_group++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"molecule/inter") == 0 ||
                 strcmp(arg[iarg+1],"molecule/intra") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude molecule", error);
        if (atom->molecule_flag == 0)
          error->all(FLERR,"Neigh_modify exclude molecule "
                     "requires atom attribute molecule");
        if (nex_mol == maxex_mol) {
          maxex_mol += EXDELTA;
          memory->grow(ex_mol_group,maxex_mol,"neigh:ex_mol_group");
          if (lmp->kokkos)
            grow_ex_mol_intra_kokkos();
          else
            memory->grow(ex_mol_intra,maxex_mol,"neigh:ex_mol_intra");
        }
        ex_mol_group[nex_mol] = group->find(arg[iarg+2]);
        if (ex_mol_group[nex_mol] == -1)
          error->all(FLERR, "Invalid exclude keyword:group {} not found", arg[iarg+2]);
        if (strcmp(arg[iarg+1],"molecule/intra") == 0)
          ex_mol_intra[nex_mol] = 1;
        else
          ex_mol_intra[nex_mol] = 0;
        nex_mol++;
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"none") == 0) {
        nex_type = nex_group = nex_mol = 0;
        iarg += 2;
      } else error->all(FLERR,"Unknown neigh_modify exclude keyword: {}", arg[iarg+1]);
    } else if (strcmp(arg[iarg],"collection/interval") == 0) {
      if (style != Neighbor::MULTI)
        error->all(FLERR,"Cannot use collection/interval command without multi setting");

      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "neigh_modify collection/interval", error);
      ncollections = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ncollections < 1)
        error->all(FLERR, "Invalid collection/interval keyword: illegal number of custom collections: {}", ncollections);
      if (iarg+2+ncollections > narg)
        error->all(FLERR, "Invalid collection/interval keyword: expected {} separate lists of types", ncollections);

      int i;

      // Invalidate old user cutoffs

      comm->ncollections_cutoff = 0;
      interval_collection_flag = 1;
      custom_collection_flag = 1;
      memory->grow(collection2cut,ncollections,"neigh:collection2cut");

      // Set upper cutoff for each collection

      double cut_interval;
      for (i = 0; i < ncollections; i++){
        cut_interval = utils::numeric(FLERR,arg[iarg+2+i],false,lmp);
        collection2cut[i] = cut_interval;

        if (i != 0)
          if (collection2cut[i-1] >= collection2cut[i])
            error->all(FLERR,"Nonsequential interval cutoffs in collection/interval setting");
      }

      iarg += 2 + ncollections;
    } else if (strcmp(arg[iarg],"collection/type") == 0) {
      if (style != Neighbor::MULTI)
        error->all(FLERR,"Cannot use collection/type command without multi setting");

      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "neigh_modify collection/type", error);
      ncollections = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ncollections < 1)
        error->all(FLERR, "Invalid collection/type keyword: illegal number of custom collections: {}", ncollections);
      if (iarg+2+ncollections > narg)
        error->all(FLERR, "Invalid collection/type keyword: expected {} separate lists of types", ncollections);

      int ntypes = atom->ntypes;
      int nlo, nhi, i, k;

      // Invalidate old user cutoffs

      comm->ncollections_cutoff = 0;
      interval_collection_flag = 0;
      custom_collection_flag = 1;
      if (!type2collection)
        memory->create(type2collection,ntypes+1,"neigh:type2collection");

      // Erase previous mapping

      for (i = 1; i <= ntypes; i++)
        type2collection[i] = -1;

      // For each custom range, define mapping for types in interval

      for (i = 0; i < ncollections; i++){
        std::vector<std::string> words = Tokenizer(arg[iarg+2+i], ",").as_vector();
        for (const auto &word : words) {
          utils::bounds(FLERR,word,1,ntypes,nlo,nhi,error);
          for (k = nlo; k <= nhi; k++) {
            if (type2collection[k] != -1)
              error->all(FLERR,"Type specified more than once in collection/type commnd");
            type2collection[k] = i;
          }
        }
      }

      // Check for undefined atom type

      for (i = 1; i <= ntypes; i++){
        if (type2collection[i] == -1) {
          error->all(FLERR,"Type missing in collection/type commnd");
        }
      }

      iarg += 2 + ncollections;
    } else error->all(FLERR,"Unknown neigh_modify keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   convenience function to allow modifying parameters from a single string
------------------------------------------------------------------------- */

void Neighbor::modify_params(const std::string &modcmd)
{
  auto args = utils::split_words(modcmd);
  auto newarg = new char*[args.size()];
  int i=0;
  for (const auto &arg : args) {
    newarg[i++] = (char *)arg.c_str();
  }
  modify_params(args.size(),newarg);
  delete[] newarg;
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
   check if any of the old requested neighbor lists are full
------------------------------------------------------------------------- */

int Neighbor::any_full()
{
  int any_full = 0;
  for (int i = 0; i < old_nrequest; i++) {
    if (old_requests[i]->full) any_full = 1;
  }
  return any_full;
}

/* ----------------------------------------------------------------------
   populate collection array for multi starting at the index istart
------------------------------------------------------------------------- */

void Neighbor::build_collection(int istart)
{
  if (style != Neighbor::MULTI)
    error->all(FLERR, "Cannot define atom collections without neighbor style multi");

  int nmax = atom->nlocal+atom->nghost;
  if (nmax > nmax_collection) {
    nmax_collection = nmax+DELTA_PERATOM;
    memory->grow(collection, nmax_collection, "neigh:collection");
  }

  if (finite_cut_flag) {
    double cut;
    int icollection;
    for (int i = istart; i < nmax; i++){
      cut = force->pair->atom2cut(i);
      collection[i] = -1;

      for (icollection = 0; icollection < ncollections; icollection++){
        if (collection2cut[icollection] >= cut) {
          collection[i] = icollection;
          break;
        }
      }

      if (collection[i] == -1)
        error->one(FLERR, "Atom cutoff exceeds interval cutoffs for multi");
    }
  } else {
    int *type = atom->type;
    for (int i = istart; i < nmax; i++){
      collection[i] = type2collection[type[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   for neighbor list statistics in Finish class
------------------------------------------------------------------------- */

bigint Neighbor::get_nneigh_full()
{
  // find a non-skip neighbor list containing full pairwise interactions
  // count neighbors in that list for stats purposes
  // allow it to be Kokkos neigh list as well

  int m;
  for (m = 0; m < old_nrequest; m++)
    if (old_requests[m]->full && !old_requests[m]->skip) break;

  bigint nneighfull = -1;
  if (m < old_nrequest) {
    nneighfull = 0;
    if (!lists[m]->kokkos && lists[m]->numneigh) {
      int inum = neighbor->lists[m]->inum;
      int *ilist = neighbor->lists[m]->ilist;
      int *numneigh = neighbor->lists[m]->numneigh;
      for (int i = 0; i < inum; i++)
        nneighfull += numneigh[ilist[i]];
    } else if (lmp->kokkos) nneighfull = lmp->kokkos->neigh_count(m);
  }
  return nneighfull;
}

bigint Neighbor::get_nneigh_half()
{
  // find a non-skip neighbor list containing half pairwise interactions
  // count neighbors in that list for stats purposes
  // allow it to be Kokkos neigh list as well

  int m;
  for (m = 0; m < old_nrequest; m++)
    if (old_requests[m]->half && !old_requests[m]->skip && lists[m] && lists[m]->numneigh) break;

  bigint nneighhalf = -1;
  if (m < old_nrequest) {
    nneighhalf = 0;
    if (!lists[m]->kokkos) {
      int inum = neighbor->lists[m]->inum;
      int *ilist = neighbor->lists[m]->ilist;
      int *numneigh = neighbor->lists[m]->numneigh;
      for (int i = 0; i < inum; i++)
        nneighhalf += numneigh[ilist[i]];
    } else if (lmp->kokkos) nneighhalf = lmp->kokkos->neigh_count(m);
  }
  return nneighhalf;
}

/* ----------------------------------------------------------------------
   add pair of atoms to bondlist array
   will only persist until the next neighbor build
------------------------------------------------------------------------- */

void Neighbor::add_temporary_bond(int i1, int i2, int btype)
{
  neigh_bond->add_temporary_bond(i1, i2, btype);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double Neighbor::memory_usage()
{
  double bytes = 0;
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
