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

/* ----------------------------------------------------------------------
   Contributing authors: Koenraad Janssens and David Olmsted (SNL)
   Modification for bcc provided by: Tegar Wicaksono (UBC)
   For a tutorial, please see "Order parameters of crystals in LAMMPS"
            (https://doi.org/10.6084/m9.figshare.1488628.v1
------------------------------------------------------------------------- */

#include "fix_orient_bcc.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define BIG 1000000000

static const char cite_fix_orient_bcc[] =
  "fix orient/bcc command:\n\n"
  "@Article{Wicaksono16,\n"
  " author = {A. T. Wicaksono, C. W. Sinclair, M. Militzer},\n"
  " title = {An atomistic study of the correlation between the migration of planar and curved grain boundaries},\n"
  " journal = {Computational Materials Science},\n"
  " year =    2016,\n"
  " volume =  117,\n"
  " pages =   {397--405}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixOrientBCC::FixOrientBCC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xifilename(nullptr), chifilename(nullptr), order(nullptr), nbr(nullptr), sort(nullptr), list(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_orient_bcc);

  MPI_Comm_rank(world,&me);

  if (narg != 11) error->all(FLERR,"Illegal fix orient/bcc command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  peratom_flag = 1;
  size_peratom_cols = 2;
  peratom_freq = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  nstats = utils::inumeric(FLERR,arg[3],false,lmp);
  direction_of_motion = utils::inumeric(FLERR,arg[4],false,lmp);
  a = utils::numeric(FLERR,arg[5],false,lmp);
  Vxi = utils::numeric(FLERR,arg[6],false,lmp);
  uxif_low = utils::numeric(FLERR,arg[7],false,lmp);
  uxif_high = utils::numeric(FLERR,arg[8],false,lmp);

  if (direction_of_motion == 0) {
    chifilename = utils::strdup(arg[9]);
    xifilename = utils::strdup(arg[10]);
  } else if (direction_of_motion == 1) {
    xifilename = utils::strdup(arg[9]);
    chifilename = utils::strdup(arg[10]);
  } else error->all(FLERR,"Illegal fix orient/bcc command");

  // initializations

  half_bcc_nn = 4;
  use_xismooth = false;
  double xicutoff = 1.57;
  xicutoffsq = xicutoff * xicutoff;
  cutsq = 0.75 * a*a*xicutoffsq;
  nmax = 0;

  // read xi and chi reference orientations from files

  if (me == 0) {
    char line[IMGMAX];
    char *result;
    int count;

    FILE *inpfile = fopen(xifilename,"r");
    if (inpfile == nullptr) error->one(FLERR,"Fix orient/bcc file open failed");
    for (int i = 0; i < 4; i++) {
      result = fgets(line,IMGMAX,inpfile);
      if (!result) error->one(FLERR,"Fix orient/bcc file read failed");
      count = sscanf(line,"%lg %lg %lg",&Rxi[i][0],&Rxi[i][1],&Rxi[i][2]);
      if (count != 3) error->one(FLERR,"Fix orient/bcc file read failed");
    }
    fclose(inpfile);

    inpfile = fopen(chifilename,"r");
    if (inpfile == nullptr) error->one(FLERR,"Fix orient/bcc file open failed");
    for (int i = 0; i < 4; i++) {
      result = fgets(line,IMGMAX,inpfile);
      if (!result) error->one(FLERR,"Fix orient/bcc file read failed");
      count = sscanf(line,"%lg %lg %lg",&Rchi[i][0],&Rchi[i][1],&Rchi[i][2]);
      if (count != 3) error->one(FLERR,"Fix orient/bcc file read failed");
    }
    fclose(inpfile);
  }

  MPI_Bcast(&Rxi[0][0],18,MPI_DOUBLE,0,world);
  MPI_Bcast(&Rchi[0][0],18,MPI_DOUBLE,0,world);

  // make copy of the reference vectors

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      half_xi_chi_vec[0][i][j] = Rxi[i][j];
      half_xi_chi_vec[1][i][j] = Rchi[i][j];
    }

  // compute xiid,xi0,xi1 for all 8 neighbors
  // xi is the favored crystal
  // want order parameter when actual is Rchi

  double xi_sq,dxi[3],rchi[3];

  xiid = 0.0;
  for (int i = 0; i < 4; i++) {
    rchi[0] = Rchi[i][0];
    rchi[1] = Rchi[i][1];
    rchi[2] = Rchi[i][2];
    find_best_ref(rchi,0,xi_sq,dxi);
    xiid += sqrt(xi_sq);
    for (int j = 0; j < 3; j++) rchi[j] = -rchi[j];
    find_best_ref(rchi,0,xi_sq,dxi);
    xiid += sqrt(xi_sq);
  }

  xiid /= 8.0;
  xi0 = uxif_low * xiid;
  xi1 = uxif_high * xiid;

  // set comm size needed by this Fix
  // NOTE: doesn't seem that use_xismooth is ever true

  if (use_xismooth) comm_forward = 62;
  else comm_forward = 50;

  added_energy = 0.0;

  nmax = atom->nmax;
  nbr = (Nbr *) memory->smalloc(nmax*sizeof(Nbr),"orient/bcc:nbr");
  memory->create(order,nmax,2,"orient/bcc:order");
  array_atom = order;

  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) order[i][0] = order[i][1] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixOrientBCC::~FixOrientBCC()
{
  delete [] xifilename;
  delete [] chifilename;
  memory->sfree(nbr);
  memory->destroy(order);
}

/* ---------------------------------------------------------------------- */

int FixOrientBCC::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::init()
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  // need a full neighbor list
  // perpetual list, built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::post_force(int /*vflag*/)
{
  int i,j,k,ii,jj,inum,jnum,m,n,nn,nsort;
  tagint id_self;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double edelta,omega;
  double dx,dy,dz,rsq,xismooth,xi_sq,duxi,duxi_other;
  double dxi[3];
  double *dxiptr;
  bool found_myself;

  // set local ptrs

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // insure nbr and order data structures are adequate size

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(nbr);
    memory->destroy(order);
    nbr = (Nbr *) memory->smalloc(nmax*sizeof(Nbr),"orient/bcc:nbr");
    memory->create(order,nmax,2,"orient/bcc:order");
    array_atom = order;
  }

  // loop over owned atoms and build Nbr data structure of neighbors
  // use full neighbor list

  added_energy = 0.0;
  int count = 0;
  int mincount = BIG;
  int maxcount = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (jnum < mincount) mincount = jnum;
    if (jnum > maxcount) {
      if (maxcount) delete [] sort;
      sort = new Sort[jnum];
      maxcount = jnum;
    }

    // loop over all neighbors of atom i
    // for those within cutsq, build sort data structure
    // store local id, rsq, delta vector, xismooth (if included)

    nsort = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      count++;

      dx = x[i][0] - x[j][0];
      dy = x[i][1] - x[j][1];
      dz = x[i][2] - x[j][2];
      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < cutsq) {
        sort[nsort].id = j;
        sort[nsort].rsq = rsq;
        sort[nsort].delta[0] = dx;
        sort[nsort].delta[1] = dy;
        sort[nsort].delta[2] = dz;
        if (use_xismooth) {
          xismooth = (xicutoffsq - 4.0*rsq/(3.0*a*a)) / (xicutoffsq - 1.0);
          sort[nsort].xismooth = 1.0 - fabs(1.0-xismooth);
        }
        nsort++;
      }
    }

    // sort neighbors by rsq distance
    // no need to sort if nsort <= 8

    if (nsort > 8) qsort(sort,nsort,sizeof(Sort),compare);

    // copy up to 8 nearest neighbors into nbr data structure
    // operate on delta vector via find_best_ref() to compute dxi

    n = MIN(8,nsort);
    nbr[i].n = n;
    if (n == 0) continue;

    double xi_total = 0.0;
    for (j = 0; j < n; j++) {
      find_best_ref(sort[j].delta,0,xi_sq,dxi);
      xi_total += sqrt(xi_sq);
      nbr[i].id[j] = sort[j].id;
      nbr[i].dxi[j][0] = dxi[0]/n;
      nbr[i].dxi[j][1] = dxi[1]/n;
      nbr[i].dxi[j][2] = dxi[2]/n;
      if (use_xismooth) nbr[i].xismooth[j] = sort[j].xismooth;
    }
    xi_total /= n;
    order[i][0] = xi_total;

    // compute potential derivative to xi

    if (xi_total < xi0) {
      nbr[i].duxi = 0.0;
      edelta = 0.0;
      order[i][1] = 0.0;
    } else if (xi_total > xi1) {
      nbr[i].duxi = 0.0;
      edelta = Vxi;
      order[i][1] = 1.0;
    } else {
      omega = MY_PI2*(xi_total-xi0) / (xi1-xi0);
      nbr[i].duxi = MY_PI*Vxi*sin(2.0*omega) / (2.0*(xi1-xi0));
      edelta = Vxi*(1 - cos(2.0*omega)) / 2.0;
      order[i][1] = omega / MY_PI2;
    }
    added_energy += edelta;
  }

  if (maxcount) delete [] sort;

  // communicate to acquire nbr data for ghost atoms

  comm->forward_comm_fix(this);

  // compute grain boundary force on each owned atom
  // skip atoms not in group

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    n = nbr[i].n;
    duxi = nbr[i].duxi;

    for (j = 0; j < n; j++) {
      dxiptr = &nbr[i].dxi[j][0];
      if (use_xismooth) {
        xismooth = nbr[i].xismooth[j];
        f[i][0] += duxi * dxiptr[0] * xismooth;
        f[i][1] += duxi * dxiptr[1] * xismooth;
        f[i][2] += duxi * dxiptr[2] * xismooth;
      } else {
        f[i][0] += duxi * dxiptr[0];
        f[i][1] += duxi * dxiptr[1];
        f[i][2] += duxi * dxiptr[2];
      }

      // m = local index of neighbor
      // id_self = ID for atom I in atom M's neighbor list
      // if M is local atom, id_self will be local ID of atom I
      // if M is ghost atom, id_self will be global ID of atom I

      m = nbr[i].id[j];
      if (m < nlocal) id_self = i;
      else id_self = tag[i];
      found_myself = false;
      nn = nbr[m].n;

      for (k = 0; k < nn; k++) {
        if (id_self == nbr[m].id[k]) {
          if (found_myself) error->one(FLERR,"Fix orient/bcc found self twice");
          found_myself = true;
          duxi_other = nbr[m].duxi;
          dxiptr = &nbr[m].dxi[k][0];
          if (use_xismooth) {
            xismooth = nbr[m].xismooth[k];
            f[i][0] -= duxi_other * dxiptr[0] * xismooth;
            f[i][1] -= duxi_other * dxiptr[1] * xismooth;
            f[i][2] -= duxi_other * dxiptr[2] * xismooth;
          } else {
            f[i][0] -= duxi_other * dxiptr[0];
            f[i][1] -= duxi_other * dxiptr[1];
            f[i][2] -= duxi_other * dxiptr[2];
          }
        }
      }
    }
  }

  // print statistics every nstats timesteps

  if (nstats && update->ntimestep % nstats == 0) {
    int total;
    MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,world);
    double ave = static_cast<double>(total)/atom->natoms;

    int min,max;
    MPI_Allreduce(&mincount,&min,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&maxcount,&max,1,MPI_INT,MPI_MAX,world);

    if (me == 0) {
      std::string mesg = fmt::format("orient step {}: {} atoms have {} "
                                     "neighbors\n", update->ntimestep,
                                     atom->natoms,total);
      mesg += fmt::format("  neighs: min = {}, max ={}, ave = {}\n",
                          min,max,ave);
      utils::logmesg(lmp,mesg);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixOrientBCC::compute_scalar()
{
  double added_energy_total;
  MPI_Allreduce(&added_energy,&added_energy_total,1,MPI_DOUBLE,MPI_SUM,world);
  return added_energy_total;
}

/* ---------------------------------------------------------------------- */

int FixOrientBCC::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,num;
  tagint id;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int m = 0;

  for (i = 0; i < n; i++) {
    k = list[i];
    num = nbr[k].n;
    buf[m++] = num;
    buf[m++] = nbr[k].duxi;

    for (j = 0; j < num; j++) {
      if (use_xismooth) buf[m++] = nbr[k].xismooth[j];
      buf[m++] = nbr[k].dxi[j][0];
      buf[m++] = nbr[k].dxi[j][1];
      buf[m++] = nbr[k].dxi[j][2];

      // id stored in buf needs to be global ID
      // if k is a local atom, it stores local IDs, so convert to global
      // if k is a ghost atom (already comm'd), its IDs are already global

      id = nbr[k].id[j];
      if (k < nlocal) id = tag[id];
      buf[m++] = id;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,num;
  int last = first + n;
  int m = 0;

  for (i = first; i < last; i++) {
    nbr[i].n = num = static_cast<int> (buf[m++]);
    nbr[i].duxi = buf[m++];

    for (j = 0; j < num; j++) {
      if (use_xismooth) nbr[i].xismooth[j] = buf[m++];
      nbr[i].dxi[j][0] = buf[m++];
      nbr[i].dxi[j][1] = buf[m++];
      nbr[i].dxi[j][2] = buf[m++];
      nbr[i].id[j] = static_cast<tagint> (buf[m++]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixOrientBCC::find_best_ref(double *displs, int which_crystal,
                                 double &xi_sq, double *dxi)
{
  int i;
  double dot,tmp;

  double  best_dot  = -1.0;         // best is biggest (smallest angle)
  int     best_i    = -1;
  int     best_sign = 0;

  for (i = 0; i < half_bcc_nn; i++) {
    dot = displs[0] * half_xi_chi_vec[which_crystal][i][0] +
      displs[1] * half_xi_chi_vec[which_crystal][i][1] +
      displs[2] * half_xi_chi_vec[which_crystal][i][2];
    if (fabs(dot) > best_dot) {
      best_dot = fabs(dot);
      best_i = i;
      if (dot < 0.0) best_sign = -1;
      else best_sign = 1;
    }
  }

  xi_sq = 0.0;
  for (i = 0; i < 3; i++) {
    tmp = displs[i] - best_sign * half_xi_chi_vec[which_crystal][best_i][i];
    xi_sq += tmp*tmp;
  }

  if (xi_sq > 0.0) {
    double xi = sqrt(xi_sq);
    for (i = 0; i < 3; i++)
      dxi[i] = (best_sign * half_xi_chi_vec[which_crystal][best_i][i] -
                displs[i]) / xi;
  } else dxi[0] = dxi[1] = dxi[2] = 0.0;
}

/* ----------------------------------------------------------------------
   compare two neighbors I and J in sort data structure
   called via qsort in post_force() method
   is a static method so can't access sort data structure directly
   return -1 if I < J, 0 if I = J, 1 if I > J
   do comparison based on rsq distance
------------------------------------------------------------------------- */

int FixOrientBCC::compare(const void *pi, const void *pj)
{
  FixOrientBCC::Sort *ineigh = (FixOrientBCC::Sort *) pi;
  FixOrientBCC::Sort *jneigh = (FixOrientBCC::Sort *) pj;

  if (ineigh->rsq < jneigh->rsq) return -1;
  else if (ineigh->rsq > jneigh->rsq) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixOrientBCC::memory_usage()
{
  double bytes = (double)nmax * sizeof(Nbr);
  bytes += (double)2*nmax * sizeof(double);
  return bytes;
}
