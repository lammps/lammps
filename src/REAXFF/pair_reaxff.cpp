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
   Contributing author:
   Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

   Per-atom energy/virial added by Ray Shan (Sandia)
   Fix reaxff/bonds and fix reaxff/species for pair_style reaxff added by
        Ray Shan (Sandia)
   Hybrid and hybrid/overlay compatibility added by Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "pair_reaxff.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "fix_reaxff.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "fix_acks2_reaxff.h"

#include <cmath>
#include <cstring>

#include "reaxff_api.h"

using namespace LAMMPS_NS;
using namespace ReaxFF;

static const char cite_pair_reax_c[] =
  "pair reaxff command: doi:10.1016/j.parco.2011.08.005\n\n"
  "@Article{Aktulga12,\n"
  " author = {H. M. Aktulga and J. C. Fogarty and S. A. Pandit and A. Y. Grama},\n"
  " title = {Parallel Reactive Molecular Dynamics: {N}umerical Methods and Algorithmic Techniques},\n"
  " journal = {Parallel Computing},\n"
  " year =    2012,\n"
  " volume =  38,\n"
  " number =  {4--5},\n"
  " pages =   {245--259}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairReaxFF::PairReaxFF(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_reax_c);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  ghostneigh = 1;

  fix_id = utils::strdup("REAXFF_" + std::to_string(instance_me));

  api = new API;

  api->system = new reax_system;
  memset(api->system,0,sizeof(reax_system));
  api->control = new control_params;
  memset(api->control,0,sizeof(control_params));
  api->data = new simulation_data;
  memset(api->data,0,sizeof(simulation_data));
  api->workspace = new storage;
  memset(api->workspace,0,sizeof(storage));
  memory->create(api->lists, LIST_N,"reaxff:lists");
  memset(api->lists,0,LIST_N * sizeof(reax_list));

  api->control->me = api->system->my_rank = comm->me;

  api->system->num_nbrs = 0;
  api->system->n = 0;                // my atoms
  api->system->N = 0;                // mine + ghosts
  api->system->local_cap = 0;
  api->system->total_cap = 0;
  api->system->my_atoms = nullptr;
  api->system->pair_ptr = this;
  api->system->mem_ptr = memory;
  api->system->error_ptr = error;
  api->control->error_ptr = error;
  api->control->lmp_ptr = lmp;

  api->system->omp_active = 0;

  fix_reaxff = nullptr;
  tmpid = nullptr;
  tmpbo = nullptr;

  nextra = 14;
  pvector = new double[nextra];

  setup_flag = 0;
  fixspecies_flag = 0;
  nmax = 0;
  list_blocking_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairReaxFF::~PairReaxFF()
{
  if (copymode) return;

  if (fix_reaxff) modify->delete_fix(fix_id);
  delete[] fix_id;

  if (setup_flag) {

    // deallocate reax data-structures

    if (api->control->tabulate) Deallocate_Lookup_Tables(api->system);
    if (api->control->hbond_cut > 0) Delete_List(api->lists+HBONDS);

    Delete_List(api->lists+BONDS);
    Delete_List(api->lists+THREE_BODIES);
    Delete_List(api->lists+FAR_NBRS);

    DeAllocate_Workspace(api->workspace);
    DeAllocate_System(api->system);
  }

  delete api->system;
  delete api->control;
  delete api->data;
  delete api->workspace;
  memory->destroy(api->lists);
  delete api;

  // deallocate interface storage
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);

    delete[] chi;
    delete[] eta;
    delete[] gamma;
    delete[] bcut_acks2;
  }

  memory->destroy(tmpid);
  memory->destroy(tmpbo);

  delete[] pvector;
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");
  map = new int[n+1];

  chi = new double[n+1];
  eta = new double[n+1];
  gamma = new double[n+1];
  bcut_acks2 = new double[n+1];
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");

  if (comm->me == 0) {
    // read name of control file or use default controls

    if (strcmp(arg[0],"NULL") == 0) {
      api->control->tabulate = 0;

      api->control->bond_cut = 5.;
      api->control->hbond_cut = 7.50;
      api->control->thb_cut = 0.001;
      api->control->thb_cutsq = 0.00001;
      api->control->bg_cut = 0.3;

      api->control->nthreads = 1;

    } else Read_Control_File(arg[0], api->control);
  }
  MPI_Bcast(api->control,sizeof(control_params),MPI_CHAR,0,world);

  // must reset these to local values after broadcast
  api->control->me = comm->me;
  api->control->error_ptr = error;
  api->control->lmp_ptr = lmp;

  // default values

  qeqflag = 1;
  api->control->lgflag = 0;
  api->control->enobondsflag = 1;
  api->system->mincap = REAX_MIN_CAP;
  api->system->minhbonds = REAX_MIN_HBONDS;
  api->system->safezone = REAX_SAFE_ZONE;
  api->system->saferzone = REAX_SAFER_ZONE;

  // process optional keywords

  int iarg = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"checkqeq") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      qeqflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"enobonds") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->control->enobondsflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"lgvdw") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->control->lgflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"safezone") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->system->safezone = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (api->system->safezone < 0.0)
        error->all(FLERR,"Illegal pair_style reaxff safezone command");
      api->system->saferzone = api->system->safezone*1.2 + 0.2;
      iarg += 2;
    } else if (strcmp(arg[iarg],"mincap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->system->mincap = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (api->system->mincap < 0)
        error->all(FLERR,"Illegal pair_style reaxff mincap command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"minhbonds") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->system->minhbonds = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (api->system->minhbonds < 0)
        error->all(FLERR,"Illegal pair_style reaxff minhbonds command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"list/blocking") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      list_blocking_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tabulate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reaxff command");
      api->control->tabulate = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (api->control->tabulate < 0)
        error->all(FLERR,"Illegal pair_style reaxff tabulate command");
      iarg += 2;
    } else error->all(FLERR,"Illegal pair_style reaxff command");
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::coeff(int nargs, char **args)
{
  if (!allocated) allocate();

  if (nargs != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read ffield file

  Read_Force_Field(args[2], &(api->system->reax_param), api->control, world);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"

  int itmp = 0;
  int nreax_types = api->system->reax_param.num_atom_types;
  for (int i = 3; i < nargs; i++) {
    if (strcmp(args[i],"NULL") == 0) {
      map[i-2] = -1;
      itmp ++;
      continue;
    }
  }

  int n = atom->ntypes;

  // pair_coeff element map
  for (int i = 3; i < nargs; i++)
    for (int j = 0; j < nreax_types; j++)
      if (utils::lowercase(args[i]) == utils::lowercase(api->system->reax_param.sbp[j].name)) {
        map[i-2] = j;
        itmp ++;
      }

  // error check
  if (itmp != n)
    error->all(FLERR,"Non-existent ReaxFF type");

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

}

/* ---------------------------------------------------------------------- */

void PairReaxFF::init_style()
{
  if (!atom->q_flag) error->all(FLERR,"Pair style reaxff requires atom attribute q");

  auto acks2_fixes = modify->get_fix_by_style("^acks2/reax");
  int have_qeq = modify->get_fix_by_style("^qeq/reax").size()
    + modify->get_fix_by_style("^qeq/shielded").size() + acks2_fixes.size();

  if (qeqflag && (have_qeq != 1))
    error->all(FLERR,"Pair style reaxff requires use of exactly one of the "
               "fix qeq/reaxff or fix qeq/shielded or fix acks2/reaxff commands");

  api->system->acks2_flag = acks2_fixes.size();
  if (api->system->acks2_flag)
    api->workspace->s = (dynamic_cast<FixACKS2ReaxFF *>(acks2_fixes.front()))->get_s();

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts
  api->system->wsize = comm->nprocs;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style reaxff requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style reaxff requires newton pair on");

  // need a half neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  neighbor->add_request(this, NeighConst::REQ_GHOST | NeighConst::REQ_NEWTON_OFF);

  cutmax = MAX3(api->control->nonb_cut, api->control->hbond_cut, api->control->bond_cut);
  if ((cutmax < 2.0*api->control->bond_cut) && (comm->me == 0))
    error->warning(FLERR,"Total cutoff < 2*bond cutoff. May need to use an "
                   "increased neighbor list skin.");

  if (fix_reaxff == nullptr)
    fix_reaxff = dynamic_cast<FixReaxFF *>(modify->add_fix(fmt::format("{} all REAXFF",fix_id)));
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::setup()
{
  int oldN;
  int mincap = api->system->mincap;
  double safezone = api->system->safezone;

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts
  oldN = api->system->N;

  if (setup_flag == 0) {

    setup_flag = 1;

    int *num_bonds = fix_reaxff->num_bonds;
    int *num_hbonds = fix_reaxff->num_hbonds;

    // determine the local and total capacity

    api->system->local_cap = MAX((int)(api->system->n * safezone), mincap);
    api->system->total_cap = MAX((int)(api->system->N * safezone), mincap);

    // initialize my data structures

    PreAllocate_Space(api->system, api->workspace);
    write_reax_atoms();

    api->system->wsize = comm->nprocs;

    int num_nbrs = estimate_reax_lists();
    if (num_nbrs < 0)
      error->all(FLERR,"Too many neighbors for pair style reaxff");

    Make_List(api->system->total_cap,num_nbrs,TYP_FAR_NEIGHBOR,api->lists+FAR_NBRS);
    (api->lists+FAR_NBRS)->error_ptr=error;

    write_reax_lists();

    Initialize(api->system,api->control,api->data,api->workspace,&api->lists,world);
    for (int k = 0; k < api->system->N; ++k) {
      num_bonds[k] = api->system->my_atoms[k].num_bonds;
      num_hbonds[k] = api->system->my_atoms[k].num_hbonds;
    }

  } else {

    // fill in reax datastructures

    write_reax_atoms();

    // reset the bond list info for new atoms

    for (int k = oldN; k < api->system->N; ++k)
      Set_End_Index(k, Start_Index(k, api->lists+BONDS), api->lists+BONDS);

    // check if I need to shrink/extend my data-structs

    ReAllocate(api->system, api->control, api->data, api->workspace, &api->lists);
  }
}

/* ---------------------------------------------------------------------- */

double PairReaxFF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cutghost[i][j] = cutghost[j][i] = cutmax;
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::compute(int eflag, int vflag)
{
  // communicate num_bonds once every reneighboring
  // 2 num arrays stored by fix, grab ptr to them

  if (neighbor->ago == 0) comm->forward_comm(fix_reaxff);
  int *num_bonds = fix_reaxff->num_bonds;
  int *num_hbonds = fix_reaxff->num_hbonds;

  ev_init(eflag,vflag);

  api->system->n = atom->nlocal; // my atoms
  api->system->N = atom->nlocal + atom->nghost; // mine + ghosts

  if (api->system->acks2_flag) {
    auto ifix = modify->get_fix_by_style("^acks2/reax").front();
    api->workspace->s = (dynamic_cast<FixACKS2ReaxFF*>(ifix))->get_s();
  }

  // setup data structures

  setup();

  Reset(api->system, api->control, api->data, api->workspace, &api->lists);
  api->workspace->realloc.num_far = write_reax_lists();

  // forces

  Compute_Forces(api->system,api->control,api->data,api->workspace,&api->lists);
  read_reax_forces(vflag);

  for (int k = 0; k < api->system->N; ++k) {
    num_bonds[k] = api->system->my_atoms[k].num_bonds;
    num_hbonds[k] = api->system->my_atoms[k].num_hbonds;
  }

  // energies and pressure

  if (eflag_global) {

    // Store the different parts of the energy
    // in a list for output by compute pair command

    pvector[0] = api->data->my_en.e_bond;
    pvector[1] = api->data->my_en.e_ov + api->data->my_en.e_un;
    pvector[2] = api->data->my_en.e_lp;
    pvector[3] = 0.0;
    pvector[4] = api->data->my_en.e_ang;
    pvector[5] = api->data->my_en.e_pen;
    pvector[6] = api->data->my_en.e_coa;
    pvector[7] = api->data->my_en.e_hb;
    pvector[8] = api->data->my_en.e_tor;
    pvector[9] = api->data->my_en.e_con;
    pvector[10] = api->data->my_en.e_vdW;
    pvector[11] = api->data->my_en.e_ele;
    pvector[12] = 0.0;
    pvector[13] = api->data->my_en.e_pol;
  }

  if (vflag_fdotr) virial_fdotr_compute();

// Set internal timestep counter to that of LAMMPS

  api->data->step = update->ntimestep;

  // populate tmpid and tmpbo arrays for fix reaxff/species
  int i, j;

  if (fixspecies_flag) {
    if (api->system->N > nmax) {
      memory->destroy(tmpid);
      memory->destroy(tmpbo);
      nmax = api->system->N;
      memory->create(tmpid,nmax,MAXSPECBOND,"pair:tmpid");
      memory->create(tmpbo,nmax,MAXSPECBOND,"pair:tmpbo");
    }

    for (i = 0; i < api->system->N; i ++)
      for (j = 0; j < MAXSPECBOND; j ++) {
        tmpbo[i][j] = 0.0;
        tmpid[i][j] = 0;
      }
    FindBond();
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::write_reax_atoms()
{
  int *num_bonds = fix_reaxff->num_bonds;
  int *num_hbonds = fix_reaxff->num_hbonds;

  if (api->system->N > api->system->total_cap)
    error->all(FLERR,"Too many ghost atoms");

  for (int i = 0; i < api->system->N; ++i) {
    api->system->my_atoms[i].orig_id = atom->tag[i];
    api->system->my_atoms[i].type = map[atom->type[i]];
    api->system->my_atoms[i].x[0] = atom->x[i][0];
    api->system->my_atoms[i].x[1] = atom->x[i][1];
    api->system->my_atoms[i].x[2] = atom->x[i][2];
    api->system->my_atoms[i].q = atom->q[i];
    api->system->my_atoms[i].num_bonds = num_bonds[i];
    api->system->my_atoms[i].num_hbonds = num_hbonds[i];
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::get_distance(rvec xj, rvec xi, double *d_sqr, rvec *dvec)
{
  (*dvec)[0] = xj[0] - xi[0];
  (*dvec)[1] = xj[1] - xi[1];
  (*dvec)[2] = xj[2] - xi[2];
  *d_sqr = SQR((*dvec)[0]) + SQR((*dvec)[1]) + SQR((*dvec)[2]);
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::set_far_nbr(far_neighbor_data *fdest,
                            int j, double d, rvec dvec)
{
  fdest->nbr = j;
  fdest->d = d;
  rvec_Copy(fdest->dvec, dvec);
  ivec_MakeZero(fdest->rel_box);
}

/* ---------------------------------------------------------------------- */

int PairReaxFF::estimate_reax_lists()
{
  int itr_i, itr_j, i, j;
  int num_nbrs;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double d_sqr;
  rvec dvec;
  double **x;

  int mincap = api->system->mincap;
  double safezone = api->system->safezone;

  x = atom->x;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  num_nbrs = 0;

  int numall = list->inum + list->gnum;

  for (itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    jlist = firstneigh[i];

    for (itr_j = 0; itr_j < numneigh[i]; ++itr_j) {
      j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance(x[j], x[i], &d_sqr, &dvec);

      if (d_sqr <= SQR(api->control->nonb_cut))
        ++num_nbrs;
    }
  }

  return static_cast<int> (MAX(num_nbrs*safezone, mincap*REAX_MIN_NBRS));
}

/* ---------------------------------------------------------------------- */

int PairReaxFF::write_reax_lists()
{
  int itr_i, itr_j, i, j;
  int num_nbrs;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double d_sqr, cutoff_sqr;
  rvec dvec;
  double *dist, **x;
  reax_list *far_nbrs;
  far_neighbor_data *far_list;

  x = atom->x;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  far_nbrs = api->lists + FAR_NBRS;
  far_list = far_nbrs->select.far_nbr_list;

  num_nbrs = 0;
  int inum = list->inum;
  dist = (double*) calloc(api->system->N, sizeof(double));

  int numall = list->inum + list->gnum;

  for (itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    jlist = firstneigh[i];
    Set_Start_Index(i, num_nbrs, far_nbrs);

    if (itr_i < inum)
      cutoff_sqr = SQR(api->control->nonb_cut);
    else
      cutoff_sqr = SQR(api->control->bond_cut);

    for (itr_j = 0; itr_j < numneigh[i]; ++itr_j) {
      j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance(x[j], x[i], &d_sqr, &dvec);

      if (d_sqr <= (cutoff_sqr)) {
        dist[j] = sqrt(d_sqr);
        set_far_nbr(&far_list[num_nbrs], j, dist[j], dvec);
        ++num_nbrs;
      }
    }
    Set_End_Index(i, num_nbrs, far_nbrs);
  }

  free(dist);

  return num_nbrs;
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::read_reax_forces(int /*vflag*/)
{
  for (int i = 0; i < api->system->N; ++i) {
    api->system->my_atoms[i].f[0] = api->workspace->f[i][0];
    api->system->my_atoms[i].f[1] = api->workspace->f[i][1];
    api->system->my_atoms[i].f[2] = api->workspace->f[i][2];

    atom->f[i][0] += -api->workspace->f[i][0];
    atom->f[i][1] += -api->workspace->f[i][1];
    atom->f[i][2] += -api->workspace->f[i][2];
  }

}

/* ---------------------------------------------------------------------- */

void *PairReaxFF::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"chi") == 0 && chi) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) chi[i] = api->system->reax_param.sbp[map[i]].chi;
      else chi[i] = 0.0;
    return (void *) chi;
  }
  if (strcmp(str,"eta") == 0 && eta) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) eta[i] = api->system->reax_param.sbp[map[i]].eta;
      else eta[i] = 0.0;
    return (void *) eta;
  }
  if (strcmp(str,"gamma") == 0 && gamma) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) gamma[i] = api->system->reax_param.sbp[map[i]].gamma;
      else gamma[i] = 0.0;
    return (void *) gamma;
   }
   if (strcmp(str,"bcut_acks2") == 0 && bcut_acks2) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) bcut_acks2[i] = api->system->reax_param.sbp[map[i]].bcut_acks2;
      else bcut_acks2[i] = 0.0;
    return (void *) bcut_acks2;
  }
  if (strcmp(str,"bond_softness") == 0) {
      double* bond_softness = &api->system->reax_param.gp.l[34];
    return (void *) bond_softness;
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

double PairReaxFF::memory_usage()
{
  double bytes = 0.0;

  // From pair_reax_c
  bytes += (double)1.0 * api->system->N * sizeof(int);
  bytes += (double)1.0 * api->system->N * sizeof(double);

  // From reaxff_allocate: BO
  bytes += (double)1.0 * api->system->total_cap * sizeof(reax_atom);
  bytes += (double)19.0 * api->system->total_cap * sizeof(double);
  bytes += (double)3.0 * api->system->total_cap * sizeof(int);

  // From reaxff_lists
  bytes += (double)2.0 * api->lists->n * sizeof(int);
  bytes += (double)api->lists->num_intrs * sizeof(three_body_interaction_data);
  bytes += (double)api->lists->num_intrs * sizeof(bond_data);
  bytes += (double)api->lists->num_intrs * sizeof(far_neighbor_data);
  bytes += (double)api->lists->num_intrs * sizeof(hbond_data);

  if (fixspecies_flag)
    bytes += (double)2 * nmax * MAXSPECBOND * sizeof(double);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void PairReaxFF::FindBond()
{
  int i, j, pj, nj;
  double bo_tmp, bo_cut;

  bond_data *bo_ij;
  bo_cut = 0.10;

  for (i = 0; i < api->system->n; i++) {
    nj = 0;
    for (pj = Start_Index(i, api->lists); pj < End_Index(i, api->lists); ++pj) {
      bo_ij = &(api->lists->select.bond_list[pj]);
      j = bo_ij->nbr;
      if (j < i) continue;

      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp >= bo_cut) {
        tmpid[i][nj] = j;
        tmpbo[i][nj] = bo_tmp;
        nj ++;
        if (nj > MAXSPECBOND) error->all(FLERR,"Increase MAXSPECBOND in reaxff_defs.h");
      }
    }
  }
}
