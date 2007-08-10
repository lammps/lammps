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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "limits.h"
#include "neighbor.h"
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

#define PGDELTA 1
#define LB_FACTOR 1.5
#define SMALL 1.0e-6
#define EXDELTA 1
#define BIG 1.0e20
#define CUT2BIN_RATIO 100

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{NSQ,BIN,MULTI};     // also in neigh_stencil.cpp

/* ---------------------------------------------------------------------- */

Neighbor::Neighbor(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  style = BIN;
  every = 1;
  delay = 10;
  dist_check = 1;
  pgsize = 10000;
  oneatom = 2000;

  maxlocal = 0;
    
  cutneighsq = NULL;
  cuttype = NULL;
  cuttypesq = NULL;
  fixchecklist = NULL;

  // last neighbor info

  maxhold = 0;
  xhold = NULL;

  // pair exclusion list info

  nex_type = maxex_type = 0;
  ex1_type = ex2_type = NULL;
  ex_type = NULL;

  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = NULL;

  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = NULL;

  // bin info

  maxhead = 0;
  binhead = NULL;
  maxbin = 0;
  bins = NULL;

  nstencil = maxstencil = 0;
  stencil = NULL;
  nstencil_full = maxstencil_full = 0;
  stencil_full = NULL;

  maxstencil_multi = 0;
  nstencil_multi = NULL;
  stencil_multi = NULL;
  distsq_multi = NULL;

  maxstencil_full_multi = 0;
  nstencil_full_multi = NULL;
  stencil_full_multi = NULL;
  distsq_full_multi = NULL;

  // half neighbor list info

  half = half_command = 0;
  numneigh = NULL;
  firstneigh = NULL;
  maxpage = 0;
  pages = NULL;

  // full neighbor list info

  full = 0;
  numneigh_full = NULL;
  firstneigh_full = NULL;
  maxpage_full = 0;
  pages_full = NULL;

  // shear history neighbor list info

  fix_history = NULL;
  firsttouch = NULL;
  firstshear = NULL;
  maxpage_history = 0;
  pages_touch = NULL;
  pages_shear = NULL;

  // multiple respa neighbor list info

  respa = 0;
  numneigh_inner = NULL;
  firstneigh_inner = NULL;
  numneigh_middle = NULL;
  firstneigh_middle = NULL;
  maxpage_inner = 0;
  maxpage_middle = 0;
  pages_inner = NULL;
  pages_middle = NULL;

  // bond list info

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
  memory->destroy_2d_double_array(cutneighsq);
  delete [] cuttype;
  delete [] cuttypesq;

  delete [] fixchecklist;
  memory->destroy_2d_double_array(xhold);

  memory->sfree(ex1_type);
  memory->sfree(ex2_type);
  memory->destroy_2d_int_array(ex_type);

  memory->sfree(ex1_group);
  memory->sfree(ex2_group);
  delete [] ex1_bit;
  delete [] ex2_bit;

  memory->sfree(ex_mol_group);
  delete [] ex_mol_bit;

  memory->sfree(binhead);
  memory->sfree(bins);

  memory->sfree(stencil);
  memory->sfree(stencil_full);

  if (nstencil_multi) {
    for (int i = 1; i <= atom->ntypes; i++) {
      memory->sfree(stencil_multi[i]);
      memory->sfree(distsq_multi[i]);
    }
    delete [] nstencil_multi;
    delete [] stencil_multi;
    delete [] distsq_multi;
  }
  if (nstencil_full_multi) {
    for (int i = 1; i <= atom->ntypes; i++) {
      memory->sfree(stencil_full_multi[i]);
      memory->sfree(distsq_full_multi[i]);
    }
    delete [] nstencil_full_multi;
    delete [] stencil_full_multi;
    delete [] distsq_full_multi;
  }

  memory->destroy_2d_int_array(bondlist);
  memory->destroy_2d_int_array(anglelist);
  memory->destroy_2d_int_array(dihedrallist);
  memory->destroy_2d_int_array(improperlist);

  memory->sfree(numneigh);
  memory->sfree(firstneigh);
  clear_pages();

  memory->sfree(numneigh_full);
  memory->sfree(firstneigh_full);
  clear_pages_full();

  memory->sfree(firsttouch);
  memory->sfree(firstshear);
  clear_pages_history();

  memory->sfree(numneigh_inner);
  memory->sfree(firstneigh_inner);
  memory->sfree(numneigh_middle);
  memory->sfree(firstneigh_middle);
  clear_pages_inner();
  clear_pages_middle();
}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
  int i,j,m,n;

  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;

  // error check

  if (delay > 0 && (delay % every) != 0)
    error->all("Neighbor delay must be 0 or multiple of every setting");

  // ------------------------------------------------------------------
  // settings

  // set neighbor cutoffs (force cutoff + skin)
  // trigger determines when atoms migrate and neighbor lists are rebuilt
  // cutneigh and cutneighsq determine what pairs go into neighbor list
  //   set to 0 if cutforce = 0
  // cutneighmin/max used for neighbor bin sizes
  // cutghost determines comm distance = max of cutneigh & skin
  //   may need ghosts for bonds even if all cutneigh = 0 (pair = NULL)

  triggersq = 0.25*skin*skin;

  n = atom->ntypes;
  if (cutneighsq == NULL) {
    cutneighsq = memory->create_2d_double_array(n+1,n+1,"neigh:cutneighsq");
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
    }
  }
  cutghost = MAX(cutneighmax,skin);
  cutneighmaxsq = cutneighmax * cutneighmax;

  // check other classes that can induce reneighboring in decide()

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

  // ------------------------------------------------------------------
  // memory management

  // free xhold and bins if not needed for this run

  if (dist_check == 0) {
    memory->destroy_2d_double_array(xhold);
    maxhold = 0;
    xhold = NULL;
  }

  if (style == NSQ) {
    memory->sfree(bins);
    memory->sfree(binhead);
    memory->sfree(stencil);
    memory->sfree(stencil_full);
    maxbin = maxhead = maxstencil = maxstencil_full = 0;
    binhead = NULL;
    bins = NULL;
  }

  // 1st time allocation of xhold and bins

  if (dist_check) {
    if (maxhold == 0) {
      maxhold = atom->nmax;
      xhold = memory->create_2d_double_array(maxhold,3,"neigh:xhold");
    }
  }

  if (style != NSQ) {
    if (maxbin == 0) {
      maxbin = atom->nmax;
      bins = (int *) memory->smalloc(maxbin*sizeof(int),"bins");
    }
  }
    
  // exclusion lists for type, group, molecule settings from neigh_modify

  n = atom->ntypes;

  if (nex_type == 0 && nex_group == 0 && nex_mol == 0) exclude = 0;
  else exclude = 1;

  if (nex_type) {
    memory->destroy_2d_int_array(ex_type);
    ex_type = (int **) memory->create_2d_int_array(n+1,n+1,"neigh:ex_type");

    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
	ex_type[i][j] = 0;

    for (i = 0; i < nex_type; i++) {
      if (ex1_type[i] <= 0 || ex1_type[i] > n || 
	  ex2_type[i] <= 0 || ex2_type[i] > n)
	error->all("Invalid atom type in neighbor exclusion list");
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
  // neighbor list flags and memory allocation/deallocation

  // determine whether to build half and full lists
  // query pair,fix,compute for their requirements

  half_once = full_once = 0;
  half_every = full_every = 0;
  if (force->pair) half_every = force->pair->neigh_half_every;
  if (force->pair) full_every = force->pair->neigh_full_every;

  for (i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->neigh_half_every) half_every = 1;
    if (modify->fix[i]->neigh_full_every) full_every = 1;
    if (modify->fix[i]->neigh_half_once) half_once = 1;
    if (modify->fix[i]->neigh_full_once) full_once = 1;
  }

  for (i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->neigh_half_once) half_once = 1;
    if (modify->compute[i]->neigh_full_once) full_once = 1;
  }

  half = full = 0;
  if (half_every || half_once || half_command) half = 1;
  if (full_every || full_once) full = 1;

  // determine whether to build granular history lists
  // fix_history = granular shear history fix
  
  fix_history = NULL;
  if (force->pair_match("gran/history") || force->pair_match("gran/hertzian"))
    for (i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"SHEAR_HISTORY") == 0) 
	fix_history = (FixShearHistory *) modify->fix[i];

  // determine whether to build extra rRESPA lists
  // respa = 1,2 if rRESPA requires inner,middle neighbor lists
  // set neighbor cutoffs for multiple lists

  respa = 0;
  if (update->whichflag == 0 && strcmp(update->integrate_style,"respa") == 0) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  if (respa && half_every == 0)
    error->all("Cannot use rRESPA with full neighbor lists");

  if (respa) {
    double *cut_respa = ((Respa *) update->integrate)->cutoff;
    cut_inner_sq = (cut_respa[1] + skin) * (cut_respa[1] + skin);
    cut_middle_sq = (cut_respa[3] + skin) * (cut_respa[3] + skin);
    cut_middle_inside_sq = (cut_respa[0] - skin) * (cut_respa[0] - skin);
  }

  // zero atom-length numneigh and firstneigh arrays
  // will be (re)allocated on first build()

  maxlocal = 0;

  memory->sfree(numneigh);
  memory->sfree(firstneigh);
  memory->sfree(numneigh_full);
  memory->sfree(firstneigh_full);
  memory->sfree(firsttouch);
  memory->sfree(firstshear);
  memory->sfree(numneigh_inner);
  memory->sfree(firstneigh_inner);
  memory->sfree(numneigh_middle);
  memory->sfree(firstneigh_middle);

  numneigh = numneigh_full = numneigh_inner = numneigh_middle = NULL;
  firstneigh = firstneigh_full = NULL;
  firsttouch = NULL;
  firstshear = NULL;
  firstneigh_inner = firstneigh_middle = NULL;

  // clear old neighbor lists if no longer needed (whether exist or not)

  if (half == 0) clear_pages();
  if (full == 0) clear_pages_full();
  if (fix_history == NULL) clear_pages_history();
  if (respa == 0) clear_pages_inner();
  if (respa == 0 || respa == 1) clear_pages_middle();

  // setup new neighbor lists
  // only if they don't exist so memory will persist from run to run

  if (half && pages == NULL) add_pages(0);
  if (full && pages_full == NULL) add_pages_full(0);
  if (fix_history && pages_touch == NULL) add_pages_history(0);
  if (respa >= 1 && pages_inner == NULL) add_pages_inner(0);
  if (respa == 2 && pages_middle == NULL) add_pages_middle(0);

  // set ptrs to half/full/multi/triclinic build & stencil functions

  if (half) {
    if (fix_history) {
      if (style == NSQ) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::granular_nsq_no_newton;
	else half_build = &Neighbor::granular_nsq_newton;
      } else if (style == BIN) {
	if (force->newton_pair == 0) {
	  half_build = &Neighbor::granular_bin_no_newton;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_no_newton;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_no_newton;
	} else if (triclinic) {
	  half_build = &Neighbor::granular_bin_newton_tri;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_newton_tri;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_newton_tri;
	} else {
	  half_build = &Neighbor::granular_bin_newton;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_newton;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_newton;
	}
      } else error->all("Neighbor multi not allowed with granular");

    } else if (respa) {
      if (style == NSQ) {
	if (force->newton_pair == 0) 
	  half_build = &Neighbor::respa_nsq_no_newton;
	else half_build = &Neighbor::respa_nsq_newton;
      } else if (style == BIN) {
	if (force->newton_pair == 0) {
	  half_build = &Neighbor::respa_bin_no_newton;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_no_newton;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_no_newton;
	} else if (triclinic) {
	  half_build = &Neighbor::respa_bin_newton_tri;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_newton_tri;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_newton_tri;
	} else {
	  half_build = &Neighbor::respa_bin_newton;
	  if (dimension == 3)
	    half_stencil = &Neighbor::stencil_half_3d_newton;
	  else
	    half_stencil = &Neighbor::stencil_half_2d_newton;
	}
      } else error->all("Neighbor multi not allowed with rRESPA");

    } else {
      if (style == NSQ) {
	if (force->newton_pair == 0) {
	  if (full_every) half_build = &Neighbor::half_full_no_newton;
	  else half_build = &Neighbor::half_nsq_no_newton;
	} else {
	  if (full_every) half_build = &Neighbor::half_full_newton;
	  else half_build = &Neighbor::half_nsq_newton;
	}
      } else if (style == BIN) {
	if (force->newton_pair == 0) {
	  if (full_every) {
	    half_build = &Neighbor::half_full_no_newton;
	    half_stencil = &Neighbor::stencil_none;
	  } else {
	    half_build = &Neighbor::half_bin_no_newton;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_no_newton;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_no_newton;
	  }
	} else {
	  if (full_every) {
	    half_build = &Neighbor::half_full_newton;
	    half_stencil = &Neighbor::stencil_none;
	  } else if (triclinic) {
	    half_build = &Neighbor::half_bin_newton_tri;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_newton_tri;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_newton_tri;
	  } else {
	    half_build = &Neighbor::half_bin_newton;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_newton;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_newton;
	  }
	}
      } else if (style == MULTI) {
	if (force->newton_pair == 0) {
	  if (full_every) {
	    half_build = &Neighbor::half_full_no_newton;
	    half_stencil = &Neighbor::stencil_none;
	  } else {
	    half_build = &Neighbor::half_bin_no_newton_multi;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_no_newton_multi;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_no_newton_multi;
	  }
	} else {
	  if (full_every) {
	    half_build = &Neighbor::half_full_newton;
	    half_stencil = &Neighbor::stencil_none;
	  } else if (triclinic) {
	    half_build = &Neighbor::half_bin_newton_multi_tri;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_newton_multi_tri;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_newton_multi_tri;
	  } else {
	    half_build = &Neighbor::half_bin_newton_multi;
	    if (dimension == 3)
	      half_stencil = &Neighbor::stencil_half_3d_newton_multi;
	    else
	      half_stencil = &Neighbor::stencil_half_2d_newton_multi;
	  }
	}
      }
    }

  } else half_build = NULL;

  if (full) {
    if (style == NSQ) full_build = &Neighbor::full_nsq;
    else if (style == BIN) {
      full_build = &Neighbor::full_bin;
      if (dimension == 3)
	full_stencil = &Neighbor::stencil_full_3d;
      else
	full_stencil = &Neighbor::stencil_full_2d;
    } else {
      full_build = &Neighbor::full_bin_multi;
      if (dimension == 3)
	full_stencil = &Neighbor::stencil_full_3d_multi;
      else
	full_stencil = &Neighbor::stencil_full_2d_multi;
    }
  } else full_build = NULL;

  // ------------------------------------------------------------------
  // bond neighbor lists

  // 1st time allocation of bond lists

  if (atom->molecular && atom->nbonds && maxbond == 0) {
    if (nprocs == 1) maxbond = atom->nbonds;
    else maxbond = static_cast<int> (LB_FACTOR * atom->nbonds / nprocs);
    bondlist = memory->create_2d_int_array(maxbond,3,"neigh:bondlist");
  }

  if (atom->molecular && atom->nangles && maxangle == 0) {
    if (nprocs == 1) maxangle = atom->nangles;
    else maxangle = static_cast<int> (LB_FACTOR * atom->nangles / nprocs);
    anglelist =  memory->create_2d_int_array(maxangle,4,"neigh:anglelist");
  }

  if (atom->molecular && atom->ndihedrals && maxdihedral == 0) {
    if (nprocs == 1) maxdihedral = atom->ndihedrals;
    else maxdihedral = static_cast<int> 
	   (LB_FACTOR * atom->ndihedrals / nprocs);
    dihedrallist = 
      memory->create_2d_int_array(maxdihedral,5,"neigh:dihedrallist");
  }

  if (atom->molecular && atom->nimpropers && maximproper == 0) {
    if (nprocs == 1) maximproper = atom->nimpropers;
    else maximproper = static_cast<int>
	   (LB_FACTOR * atom->nimpropers / nprocs);
    improperlist = 
      memory->create_2d_int_array(maximproper,5,"neigh:improperlist");
  }

  // set flags that determine which bond neighboring routines to use
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

  // set ptrs to intra-molecular build functions

  if (bond_off) bond_build = &Neighbor::bond_partial;
  else bond_build = &Neighbor::bond_all;

  if (angle_off) angle_build = &Neighbor::angle_partial;
  else angle_build = &Neighbor::angle_all;

  if (dihedral_off) dihedral_build = &Neighbor::dihedral_partial;
  else dihedral_build = &Neighbor::dihedral_all;

  if (improper_off) improper_build = &Neighbor::improper_partial;
  else improper_build = &Neighbor::improper_all;

  // set intra-molecular neighbor list counts to 0
  // in case all are turned off but potential is still defined

  nbondlist = nanglelist = ndihedrallist = nimproperlist = 0;
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
    if (dist_check == 0) return 1;
    else return check_distance();
  } else return 0;
}

/* ---------------------------------------------------------------------- */

int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;

  int nlocal = atom->nlocal;
  double **x = atom->x;
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - xhold[i][0];
    dely = x[i][1] - xhold[i][1];
    delz = x[i][2] - xhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > triggersq) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}

/* ----------------------------------------------------------------------
   build all needed neighbor lists every few timesteps
   half, full, bond lists are created as needed
------------------------------------------------------------------------- */

void Neighbor::build()
{
  ago = 0;
  ncalls++;

  // store current nlocal used on this build (used by fix shear/history)

  nlocal_neighbor = atom->nlocal;

  // store current atom positions if needed

  if (dist_check) {
    double **x = atom->x;
    int nlocal = atom->nlocal;
    if (nlocal > maxhold) {
      maxhold = atom->nmax;
      memory->destroy_2d_double_array(xhold);
      xhold = memory->create_2d_double_array(maxhold,3,"neigh:xhold");
    }
    for (int i = 0; i < nlocal; i++) {
      xhold[i][0] = x[i][0];
      xhold[i][1] = x[i][1];
      xhold[i][2] = x[i][2];
    }
  }

  // extend atom arrays if necessary
  // check half/full instead of half_every/full_every so memory will be
  //   allocated correctly whenever build_half() and build_full() are called

  if (atom->nlocal > maxlocal) {
    maxlocal = atom->nmax;

    if (half) {
      memory->sfree(numneigh);
      memory->sfree(firstneigh);
      numneigh = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh");
      firstneigh = (int **)
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh");
    }

    if (full) {
      memory->sfree(numneigh_full);
      memory->sfree(firstneigh_full);
      numneigh_full = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_full");
      firstneigh_full = (int **)
      memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_full");
    }

    if (fix_history) {
      memory->sfree(firsttouch);
      memory->sfree(firstshear);
      firsttouch = (int **) 
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firsttouch");
      firstshear = (double **)
	memory->smalloc(maxlocal*sizeof(double *),"neigh:firstshear");
    }

    if (respa) {
      memory->sfree(numneigh_inner);
      memory->sfree(firstneigh_inner);
      numneigh_inner = (int *)
	memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_inner");
      firstneigh_inner = (int **)
	memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_inner");
      if (respa == 2) {
	memory->sfree(numneigh_middle);
	memory->sfree(firstneigh_middle);
	numneigh_middle = (int *)
	  memory->smalloc(maxlocal*sizeof(int),"neigh:numneigh_middle");
	firstneigh_middle = (int **)
	  memory->smalloc(maxlocal*sizeof(int *),"neigh:firstneigh_middle");
      }
    }
  }

  // extend bin list if necessary

  if (style != NSQ && atom->nmax > maxbin) {
    maxbin = atom->nmax;
    memory->sfree(bins);
    bins = (int *) memory->smalloc(maxbin*sizeof(int),"bins");
  }

  // list construction for pairs and bonds
  // full comes first in case half is built from full

  if (full_every) (this->*full_build)();
  if (half_every) (this->*half_build)();

  if (atom->molecular) {
    if (atom->nbonds) (this->*bond_build)();
    if (atom->nangles) (this->*angle_build)();
    if (atom->ndihedrals) (this->*dihedral_build)();
    if (atom->nimpropers) (this->*improper_build)();
  }
}

/* ----------------------------------------------------------------------
   one-time call to build a half neighbor list made by other classes
------------------------------------------------------------------------- */

void Neighbor::build_half()
{
  (this->*half_build)();
}

/* ----------------------------------------------------------------------
   one-time call to build a full neighbor list made by other classes
------------------------------------------------------------------------- */

void Neighbor::build_full()
{
  (this->*full_build)();
}

/* ----------------------------------------------------------------------
   setup neighbor binning parameters
   bin numbering in each dimension is global:
     0 = 0.0 to binsize, 1 = binsize to 2*binsize, etc
     nbin-1,nbin,etc = bbox-binsize to binsize, bbox to bbox+binsize, etc
     -1,-2,etc = -binsize to 0.0, -2*size to -size, etc
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
   stencil() = bin offsets in 1d sense for stencil of surrounding bins
   stencil_full() = bin offsets in 1d sense for stencil for full neighbor list
------------------------------------------------------------------------- */

void Neighbor::setup_bins()
{
  // bbox lo/hi = bounding box of entire domain
  // bbox = size of bbox of entire domain
  // bsubbox lo/hi = bounding box of my subdomain extended by ghost atoms
  // for triclinic:
  //   bbox bounds all 8 corners of tilted box
  //   subdomain is in lamda coords
  //   include dimension-dependent extension via comm->cutghost
  //   domain->bbox() converts lamda extent to box coords and computes bbox

  double bbox[3],bsubboxlo[3],bsubboxhi[3];

  if (triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
    bsubboxlo[0] = domain->sublo[0] - cutghost;
    bsubboxlo[1] = domain->sublo[1] - cutghost;
    bsubboxlo[2] = domain->sublo[2] - cutghost;
    bsubboxhi[0] = domain->subhi[0] + cutghost;
    bsubboxhi[1] = domain->subhi[1] + cutghost;
    bsubboxhi[2] = domain->subhi[2] + cutghost;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
    double lo[3],hi[3];
    lo[0] = domain->sublo_lamda[0] - comm->cutghost[0];
    lo[1] = domain->sublo_lamda[1] - comm->cutghost[1];
    lo[2] = domain->sublo_lamda[2] - comm->cutghost[2];
    hi[0] = domain->subhi_lamda[0] + comm->cutghost[0];
    hi[1] = domain->subhi_lamda[1] + comm->cutghost[1];
    hi[2] = domain->subhi_lamda[2] + comm->cutghost[2];
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
  if (style == BIN) binsize_optimal = 0.5*cutneighmax;
  else binsize_optimal = 0.5*cutneighmin;
  if (binsize_optimal == 0.0) binsize_optimal = bbox[0];
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain

  if (bbox[0]*binsizeinv > INT_MAX || bbox[1]*binsizeinv > INT_MAX ||
      bbox[2]*binsizeinv > INT_MAX)
    error->all("Domain too large for neighbor bins");

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
    error->all("Cannot use neighbor bins - box size << cutoff");
  
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

  mbins = mbinx*mbiny*mbinz;
  if (mbins > maxhead) {
    maxhead = mbins;
    memory->sfree(binhead);
    binhead = (int *) memory->smalloc(maxhead*sizeof(int),"neigh:binhead");
  }

  // create stencil of bins to search over in neighbor list construction
  // sx,sy,sz = max range of stencil extent
  // stencil is empty if cutneighmax = 0.0

  int sx = static_cast<int> (cutneighmax*bininvx);
  if (sx*binsizex < cutneighmax) sx++;
  int sy = static_cast<int> (cutneighmax*bininvy);
  if (sy*binsizey < cutneighmax) sy++;
  int sz = static_cast<int> (cutneighmax*bininvz);
  if (sz*binsizez < cutneighmax) sz++;
  if (dimension == 2) sz = 0;

  // allocate stencil memory and create stencil(s)
  // check half/full instead of half_every/full_every so stencils will be
  //   allocated correctly whenever build_half() and build_full() are called

  stencil_allocate(sx,sy,sz);
  if (half) (this->*half_stencil)(sx,sy,sz);
  if (full) (this->*full_stencil)(sx,sy,sz);
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
  if (narg != 2) error->all("Illegal neighbor command");

  skin = atof(arg[0]);
  if (skin < 0.0) error->all("Illegal neighbor command");

  if (strcmp(arg[1],"nsq") == 0) style = NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = BIN;
  else if (strcmp(arg[1],"multi") == 0) style = MULTI;
  else error->all("Illegal neighbor command");
}

/* ----------------------------------------------------------------------
   modify parameters of the pair-wise neighbor build
------------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      every = atoi(arg[iarg+1]);
      if (every <= 0) error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      delay = atoi(arg[iarg+1]);
      if (delay < 0) error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) dist_check = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) dist_check = 0;
      else error->all("Illegal neigh_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      pgsize = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");
      oneatom = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) error->all("Illegal neigh_modify command");

      if (strcmp(arg[iarg+1],"type") == 0) {
	if (iarg+4 > narg) error->all("Illegal neigh_modify command");
	if (nex_type == maxex_type) {
	  maxex_type += EXDELTA;
	  ex1_type = (int *) memory->srealloc(ex1_type,maxex_type*sizeof(int),
					      "neigh:ex1_type");
	  ex2_type = (int *) memory->srealloc(ex2_type,maxex_type*sizeof(int),
					      "neigh:ex2_type");
	}
	ex1_type[nex_type] = atoi(arg[iarg+2]);
	ex2_type[nex_type] = atoi(arg[iarg+3]);
	nex_type++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"group") == 0) {
	if (iarg+4 > narg) error->all("Illegal neigh_modify command");
	if (nex_group == maxex_group) {
	  maxex_group += EXDELTA;
	  ex1_group = 
	    (int *) memory->srealloc(ex1_group,maxex_group*sizeof(int),
				     "neigh:ex1_group");
	  ex2_group = 
	    (int *) memory->srealloc(ex2_group,maxex_group*sizeof(int),
				     "neigh:ex2_group");
	}
	ex1_group[nex_group] = group->find(arg[iarg+2]);
	ex2_group[nex_group] = group->find(arg[iarg+3]);
	if (ex1_group[nex_group] == -1 || ex2_group[nex_group] == -1)
	  error->all("Invalid group ID in neigh_modify command");
	nex_group++;
	iarg += 4;

      } else if (strcmp(arg[iarg+1],"molecule") == 0) {
	if (iarg+3 > narg) error->all("Illegal neigh_modify command");
	if (atom->molecular == 0) {
	  char *str =
	    "Must use molecular atom style with neigh_modify exclude molecule";
	  error->all(str);
	}
	if (nex_mol == maxex_mol) {
	  maxex_mol += EXDELTA;
	  ex_mol_group = 
	    (int *) memory->srealloc(ex_mol_group,maxex_mol*sizeof(int),
				     "neigh:ex_mol_group");
	}
	ex_mol_group[nex_mol] = group->find(arg[iarg+2]);
	if (ex_mol_group[nex_mol] == -1)
	  error->all("Invalid group ID in neigh_modify command");
	nex_mol++;
	iarg += 3;

      } else if (strcmp(arg[iarg+1],"none") == 0) {
	nex_type = nex_group = nex_mol = 0;
	iarg += 2;
      } else error->all("Illegal neigh_modify command");

    } else error->all("Illegal neigh_modify command");
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

int Neighbor::memory_usage()
{
  int bytes = 0;

  bytes += maxhold*3 * sizeof(double);

  if (style == NSQ) {
    bytes += maxbin * sizeof(int);
    bytes += maxhead * sizeof(int);
    bytes += maxstencil * sizeof(int);
    bytes += maxstencil_full * sizeof(int);
  }

  if (half) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage*pgsize * sizeof(int);
  }

  if (full) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage_full*pgsize * sizeof(int);
  }

  if (fix_history) {
    bytes += maxlocal * sizeof(int *);
    bytes += maxlocal * sizeof(double *);
    bytes += maxpage_history*pgsize * sizeof(int);
    bytes += maxpage_history*pgsize*3 * sizeof(double);
  }

  if (respa) {
    bytes += maxlocal * sizeof(int);
    bytes += maxlocal * sizeof(int *);
    bytes += maxpage_inner*pgsize * sizeof(int);
    if (respa == 2) {
      bytes += maxlocal * sizeof(int);
      bytes += maxlocal * sizeof(int *);
      bytes += maxpage_middle*pgsize * sizeof(int);
    }
  }

  bytes += maxbond*3 * sizeof(int);
  bytes += maxangle*4 * sizeof(int);
  bytes += maxdihedral*5 * sizeof(int);
  bytes += maximproper*5 * sizeof(int);

  return bytes;
}

/* ----------------------------------------------------------------------
   add pages to half/full/granular/rRESPA neighbor lists, starting at npage
------------------------------------------------------------------------- */

void Neighbor::add_pages(int npage)
{
  maxpage += PGDELTA;
  pages = (int **) 
    memory->srealloc(pages,maxpage*sizeof(int *),"neigh:pages");
  for (int i = npage; i < maxpage; i++)
    pages[i] = (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages[i]");
}

void Neighbor::add_pages_full(int npage)
{
  maxpage_full += PGDELTA;
  pages_full = (int **) 
    memory->srealloc(pages_full,maxpage_full*sizeof(int *),"neigh:pages_full");
  for (int i = npage; i < maxpage_full; i++)
    pages_full[i] =
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_full[i]");
}

void Neighbor::add_pages_history(int npage)
{
  maxpage_history += PGDELTA;
  pages_touch = (int **)
    memory->srealloc(pages_touch,maxpage_history*sizeof(int *),
		     "neigh:pages_touch");
  pages_shear = (double **)
    memory->srealloc(pages_shear,maxpage_history*sizeof(double *),
		     "neigh:pages_shear");
  for (int i = npage; i < maxpage_history; i++) {
    pages_touch[i] = (int *)
      memory->smalloc(pgsize*sizeof(int),"neigh:pages_touch[i]");
    pages_shear[i] = (double *)
      memory->smalloc(3*pgsize*sizeof(double),"neigh:pages_shear[i]");
  }
}

void Neighbor::add_pages_inner(int npage_inner)
{
  maxpage_inner += PGDELTA;
  pages_inner = (int **) 
    memory->srealloc(pages_inner,maxpage_inner*sizeof(int *),
		     "neigh:pages_inner");
  for (int i = npage_inner; i < maxpage_inner; i++)
    pages_inner[i] = 
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_inner[i]");
}

void Neighbor::add_pages_middle(int npage_middle)
{
  maxpage_middle += PGDELTA;
  pages_middle = (int **) 
    memory->srealloc(pages_middle,maxpage_middle*sizeof(int *),
		     "neigh:pages_middle");
  for (int i = npage_middle; i < maxpage_middle; i++)
    pages_middle[i] = 
      (int *) memory->smalloc(pgsize*sizeof(int),"neigh:pages_middle[i]");
}

/* ----------------------------------------------------------------------
   clear half/full/granular/rRESPA neighbor lists
------------------------------------------------------------------------- */

void Neighbor::clear_pages()
{
  for (int i = 0; i < maxpage; i++) memory->sfree(pages[i]);
  memory->sfree(pages);
  pages = NULL;
  maxpage = 0;
}

void Neighbor::clear_pages_full()
{
  for (int i = 0; i < maxpage_full; i++) memory->sfree(pages_full[i]);
  memory->sfree(pages_full);
  pages_full = NULL;
  maxpage_full = 0;
}

void Neighbor::clear_pages_history()
{
  for (int i = 0; i < maxpage_history; i++) memory->sfree(pages_touch[i]);
  for (int i = 0; i < maxpage_history; i++) memory->sfree(pages_shear[i]);
  memory->sfree(pages_touch);
  memory->sfree(pages_shear);
  pages_touch = NULL;
  pages_shear = NULL;
  maxpage_history = 0;
}

void Neighbor::clear_pages_inner()
{
  for (int i = 0; i < maxpage_inner; i++) memory->sfree(pages_inner[i]);
  memory->sfree(pages_inner);
  pages_inner = NULL;
  maxpage_inner = 0;
}

void Neighbor::clear_pages_middle()
{
  for (int i = 0; i < maxpage_middle; i++) memory->sfree(pages_middle[i]);
  memory->sfree(pages_middle);
  pages_middle = NULL;
  maxpage_middle = 0;
}

/* ----------------------------------------------------------------------
   determine if atom j is in special list of atom i
   if it is not, return 0
   if it is and special flag is 0 (both coeffs are 0.0), return -1
   if it is and special flag is 1 (both coeffs are 1.0), return 0
   if it is and special flag is 2 (otherwise), return 1,2,3
     for which neighbor it is (and which coeff it maps to)
------------------------------------------------------------------------- */

int Neighbor::find_special(int i, int j)
{
  int *list = atom->special[i];
  int n1 = atom->nspecial[i][0];
  int n2 = atom->nspecial[i][1];
  int n3 = atom->nspecial[i][2];
  int tag = atom->tag[j];

  for (int i = 0; i < n3; i++) {
    if (list[i] == tag) {
      if (i < n1) {
	if (special_flag[1] == 0) return -1;
	else if (special_flag[1] == 1) return 0;
	else return 1;
      } else if (i < n2) {
	if (special_flag[2] == 0) return -1;
	else if (special_flag[2] == 1) return 0;
	else return 2;
      } else {
	if (special_flag[3] == 0) return -1;
	else if (special_flag[3] == 1) return 0;
	else return 3;
      }
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void Neighbor::bin_atoms()
{
  int i,ibin,nlocal,nall;
  double **x;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  x = atom->x;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin ghost atoms 1st, so will be at end of linked list
  // then bin owned atoms

  for (i = nlocal; i < nall; i++) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  for (i = 0; i < nlocal; i++) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  /*
  for (i = nall-1; i >= 0; i--) {
    ibin = coord2bin(x[i]);
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
  }
  */
}

/* ----------------------------------------------------------------------
   convert atom coords into local bin #
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are put in correct bins
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int Neighbor::coord2bin(double *x)
{
  int ix,iy,iz;

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx - mbinxlo;
  else if (x[0] >= bboxlo[0])
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - mbinxlo;
  else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - mbinxlo - 1;
  
  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny - mbinylo;
  else if (x[1] >= bboxlo[1])
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - mbinylo;
  else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - mbinylo - 1;
  
  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz - mbinzlo;
  else if (x[2] >= bboxlo[2])
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - mbinzlo;
  else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - mbinzlo - 1;

  return (iz*mbiny*mbinx + iy*mbinx + ix + 1);
}

/* ----------------------------------------------------------------------
   test if atom pair i,j is excluded from neighbor list
   due to type, group, molecule settings from neigh_modify command
   return 1 if should be excluded, 0 if included
------------------------------------------------------------------------- */

int Neighbor::exclusion(int i, int j, int *type, int *mask, int *molecule)
{
  int m;

  if (nex_type && ex_type[type[i]][type[j]]) return 1;

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
