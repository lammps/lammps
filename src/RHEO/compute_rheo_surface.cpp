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
   Contributing authors:
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "compute_rheo_surface.h"

#include "fix_rheo.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_solids.h"
#include "atom.h"
#include "memory.h"
#include "atom.h"
#include "comm.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "error.h"
#include "force.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::ComputeRHEOSurface(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix rheo/surface command");

  cut = utils::numeric(FLERR,arg[3],false,lmp);
  divR_limit = utils::numeric(FLERR,arg[4],false,lmp);
  coord_limit = utils::inumeric(FLERR,arg[5],false,lmp);

  divr_flag = 1;
  if (narg == 7) {
    divr_flag = 0;
  }

  int dim = domain->dimension;

  peratom_flag = 1;
  size_peratom_cols = dim;
  peratom_freq = 1;


  comm_forward = 2;
  comm_reverse = dim*dim + 1;

  cutsq = cut*cut;

  B = nullptr;
  gradC = nullptr;
  n_surface = nullptr;

  int nall = atom->nlocal + atom->nghost;
  nmax = nall;
  memory->create(B,nmax,dim*dim,"fix/rheo/surface:B");
  memory->create(gradC,nmax,dim*dim,"fix/rheo/surface:gradC");
  memory->create(n_surface,nmax,dim,"fix/rheo/surface:B");
  array_atom = n_surface;

  compute_kernel = nullptr;
  compute_solids = NULL;
  fix_rheo = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOSurface::~ComputeRHEOSurface()
{
  if (modify->nfix) modify->delete_fix("PROPERTY_ATOM_RHEO_SURFACE");

  memory->destroy(B);
  memory->destroy(gradC);
  memory->destroy(n_surface);
}

void ComputeRHEOSurface::post_constructor()
{
  //Store persistent per atom quantities
  char **fixarg = new char*[5];
  fixarg[0] = (char *) "PROPERTY_ATOM_RHEO_SURFACE";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "property/atom";
  fixarg[3] = (char *) "d_divr";
  fixarg[4] = (char *) "d_rsurf";
  modify->add_fix(5,fixarg,1);

  int temp_flag;
  index_divr = atom->find_custom("divr", temp_flag);
  if ((index_divr < 0) || (temp_flag != 1))
      error->all(FLERR, "Pair rheo/surface can't find fix property/atom divr");

  index_rsurf = atom->find_custom("rsurf", temp_flag);
  if ((index_rsurf < 0) || (temp_flag != 1))
      error->all(FLERR, "Pair rheo/surface can't find fix property/atom rsurf");

  delete [] fixarg;
  divr = atom->dvector[index_divr];
}
/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init()
{
  // need an occasional full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;

  int flag;
  int ifix = modify->find_fix_by_style("rheo");
  if (ifix == -1) error->all(FLERR, "Need to define fix rheo to use fix rheo/surface");
  fix_rheo = ((FixRHEO *) modify->fix[ifix]);
  compute_kernel = fix_rheo->compute_kernel;
  compute_solids = fix_rheo->compute_solids;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::setup_pre_force(int /*vflag*/)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::pre_force(int /*vflag*/)
{
  int i, j, ii, jj, jnum, a, b, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, r, wp, Voli, Volj, rhoi, rhoj;
  int *jlist;
  double *dWij, *dWji;
  double dx[3];

  divr = atom->dvector[index_divr];

  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double **x = atom->x;
  int *surface = atom->surface;
  int *phase = atom->phase;
  double *rsurf = atom->dvector[index_rsurf];
  int newton = force->newton;
  int dim = domain->dimension;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rho = atom->rho;
  double *temp = atom->temp;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (nmax <= nall) {
    nmax = nall;
    memory->destroy(B);
    memory->destroy(gradC);
    memory->destroy(n_surface);

    memory->create(B,nmax,dim*dim,"fix/rheo/surface:B");
    memory->create(gradC,nmax,dim*dim,"fix/rheo/surface:gradC");
    memory->create(n_surface,nmax,dim,"fix/rheo/surface:n_surface");
    array_atom = n_surface;
  }

  for (i = 0; i < nall; i++) {
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        B[i][a*dim + b] = 0.0;
        gradC[i][a*dim + b] = 0.0;
      }
      n_surface[i][a] = 0.0;
    }
    divr[i] = 0.0;
    surface[i] = 0;
  }

  // loop over neighbors to calculate the average orientation of neighbors
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    itype = type[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      dx[0] = delx;
      dx[1] = dely;
      dx[2] = delz;
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq) {
        jtype = type[j];

        rhoi = rho[i];
        rhoj = rho[j];

        // Add corrections for walls
        if (phase[i] <= FixRHEO::FLUID_MAX && phase[j] > FixRHEO::FLUID_MAX) {
          rhoj = compute_solids->correct_rho(j,i);
        } else if (phase[i] > FixRHEO::FLUID_MAX && phase[j] <= FixRHEO::FLUID_MAX) {
          rhoi = compute_solids->correct_rho(i,j);
        } else if (phase[i] > FixRHEO::FLUID_MAX && phase[j] > FixRHEO::FLUID_MAX) {
          rhoi = 1.0;
          rhoj = 1.0;
        }

        Voli = mass[itype]/rhoi;
        Volj = mass[jtype]/rhoj;

        //compute kernel gradient
        wp = compute_kernel->calc_dw_quintic(i, j, delx, dely, delz, sqrt(rsq),compute_kernel->dWij,compute_kernel->dWji);
        //wp = compute_kernel->calc_dw(i, j, delx, dely, delz, sqrt(rsq));//,compute_kernel->dWij,compute_kernel->dWji);

        dWij = compute_kernel->dWij;
        dWji = compute_kernel->dWji;

        for (a=0; a<dim; a++){
          divr[i] -= dWij[a]*dx[a]*Volj; // dx = xi-xj = xji = -xij
          gradC[i][a] +=  dWij[a]*Volj;
        }

        if (j < nlocal || newton) {
          for (a=0; a<dim; a++){
            divr[j] += dWji[a]*dx[a]*Voli;
            gradC[j][a] +=  dWji[a]*Voli;
          }
        }
      }
    }
  }

  comm_stage = 0;
  comm_reverse = dim*dim + 1;   // gradC and divr
  comm_forward = 1;             // divr
  if (newton) comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);

  int *coordination = compute_kernel->coordination;
  // Find the free-surface
  //0-bulk 1-surf vicinity  2-surface 3-splash
  if (divr_flag) {
    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        surface[i] = 0;
        rsurf[i] = cut; //Maximum range that can be seen
        if (divr[i] < divR_limit) {
          surface[i] = 2;
          rsurf[i] = 0.0;
          if (coordination[i] < coord_limit) surface[i] = 3;
        }
      }
    }
  } else {
    for (i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        surface[i] = 0;
        rsurf[i] = cut; //Maximum range that can be seen
        if (coordination[i] < divR_limit) {
          surface[i] = 2;
          rsurf[i] = 0.0;
          if (coordination[i] < coord_limit) surface[i] = 3;
        }
      }
    }
  }

  //comm_stage = 1;
  //comm_forward = 1;
  //comm->forward_comm_fix(this);  // communicate free surface particles

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq) {
        r = sqrt(rsq);
        if (surface[i] == 0 && surface[j] == 2) surface[i] = 1;
        if (surface[j] == 0 && surface[i] == 2) surface[j] = 1;
        if (surface[j] == 2) rsurf[i] = MIN(rsurf[i], r);
        if (surface[i] == 2) rsurf[j] = MIN(rsurf[j], r);
      }
    }
  }

  comm_stage = 1;
  comm_reverse = 2;
  comm_forward = 2;
  if (newton) comm->reverse_comm_fix(this);
  comm->forward_comm_fix(this);

  //Now loop again and for each surface particle (2)
  // find its neighbors that are bulk (0) and convert to surface vicinity (1)
  // if the surface particle has no (0) or (1) neighbors then it is a spash (3)

  //for (ii = 0; ii < inum; ii++) {  // is this the right i and j loop for this?
  //  i = ilist[ii];
  //
  //  if (surface[i]!=2) continue; //Only consider surface particles
  //
  //  bool nobulkneigh = true; // whether we have no bulk neighbors
  //  xtmp = x[i][0];
  //  ytmp = x[i][1];
  //  ztmp = x[i][2];
  //  jlist = firstneigh[i];
  //  jnum = numneigh[i];
  //
  //  for (jj = 0; jj < jnum; jj++) {
  //    j = jlist[jj];
  //    j &= NEIGHMASK;
  //
  //    //other surface or splash neighbors do not need labeling
  //    if (surface[j]>=2){
  //      continue;
  //    }
  //
  //    //check distance criterion rij < h = cutsq/9 for quintic kernel
  //    delx = xtmp - x[j][0];
  //    dely = ytmp - x[j][1];
  //    delz = ztmp - x[j][2];
  //    dx[0] = 3.0*delx;   // multiplied by three here to make criterion r<h instead of r<3*h
  //    dx[1] = 3.0*dely;
  //    dx[2] = 3.0*delz;
  //    rsq = delx * delx + dely * dely + delz * delz;
  //    if (rsq < cutsq) {
  //      //We have identified 1 bulk fluid neighbor
  //      nobulkneigh = false;
  //      //that bulk fluid neighbor is in the vicinity of hte surface
  //      surface[j] = 1;
  //    }
  //  }
  //  if (nobulkneigh){
  //    surface[i] = 3;
  //  }
  //}
//
//  //Reverse comm surface?
//
//  // loop over neighbors to calculate the average orientation
//  // skip for bulk or splash
//  for (ii = 0; ii < inum; ii++) {
//    i = ilist[ii];
//    if ((surface[i]==0)||(surface[i]==3)){
//      continue;
//    }
//
//    itype = type[i];
//    rhoi = rho[i];
//    Voli = mass[itype]/rhoi;
//
//    xtmp = x[i][0];
//    ytmp = x[i][1];
//    ztmp = x[i][2];
//
//    jlist = firstneigh[i];
//    jnum = numneigh[i];
//    for (jj = 0; jj < jnum; jj++) {
//      j = jlist[jj];
//      j &= NEIGHMASK;
//
//      delx = xtmp - x[j][0];
//      dely = ytmp - x[j][1];
//      delz = ztmp - x[j][2];
//      dx[0] = delx;   // multiplied by three here to make criterion r<h instead of r<3*h
//      dx[1] = dely;
//      dx[2] = delz;
//      rsq = delx * delx + dely * dely + delz * delz;
//      if (rsq < cutsq) {
//
//        jtype = type[j];
//        rhoj = rho[j];
//        Volj = mass[jtype]/rhoj;
//
//        for (a=0; a<dim; a++){
//          for (b=0; b<dim; b++){
//            B[i][a*dim+b] -= dx[a]*dWij[b]*Volj;
//          }
//        }
//
//        if (j < nlocal || newton) {
//          for (a=0; a<dim; a++){
//            for (b=0; b<dim; b++){
//              B[j][a*dim+b] += dx[a]*dWji[b]*Voli;
//            }
//          }
//        }
//      }
//    }
//  }
//  //reverse comm to populate B[j] if Newton is on
//  comm_stage = 2;
//  comm_reverse = dim*dim;        // B
//  if (newton) comm->reverse_comm_fix(this);
//
//
//  // Now need to invert each B
//  int status, s;
//  //LU requires a permuation matrix
//  gsl_permutation * p = gsl_permutation_alloc(dim);
//  for (ii = 0; ii < inum; ii++) {
//    i = ilist[ii];
//    if ((surface[i]==0)||(surface[i]==3)){
//      continue;
//    }
//
//    //Use gsl to get Binv
//    //B is not symmteric so we will use a LU decomp
//    gsl_matrix_view gB = gsl_matrix_view_array(B[i],dim,dim);
//    status = 0;
//    status = gsl_linalg_LU_decomp(&gB.matrix,p,&s);  //B[i] is now the LU decomp
//    // check if decomposition failure
//    if (status) {
//      fprintf(stderr, "failed, gsl_errno=%d.n", status);
//      continue;
//    } else {
//      gsl_linalg_LU_invx(&gB.matrix,p); //B[i] is now inv(B[i])
//    }
//  }
//  gsl_permutation_free(p);
  double maggC = 0.0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      maggC=0;
      for (a=0;a<dim;a++){
        maggC += gradC[i][a]*gradC[i][a];
      }
      maggC = sqrt(maggC) + 1e-10;
      for (a=0;a<dim;a++){
        n_surface[i][a] = -gradC[i][a]/maggC;
      }//dr can then be calculated by fix vshift
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::pack_reverse_comm(int n, int first, double *buf)
{
  int i,a,b,k,m,last;
  int dim = domain->dimension;
  int *surface = atom->surface;
  double *rsurf = atom->dvector[index_rsurf];

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      buf[m++] = divr[i];
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          buf[m++] = gradC[i][a*dim + b];
    } else if (comm_stage == 1) {
      buf[m++] = (double) surface[i];
      buf[m++] = rsurf[i];
    } else if (comm_stage == 2) {
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          buf[m++] = B[i][a*dim + b];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,a,b,k,j,m;
  int dim = domain->dimension;
  int *surface = atom->surface;
  double *rsurf = atom->dvector[index_rsurf];

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      divr[j] += buf[m++];
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          gradC[j][a*dim + b] += buf[m++];
    } else if (comm_stage == 1) {
      int temp = (int) buf[m++];
      surface[j] = MAX(surface[j], temp);
      double temp2 = buf[m++];
      rsurf[j] = MIN(rsurf[j], temp2);
    } else if (comm_stage == 2) {
      for (a = 0; a < dim; a ++ )
        for (b = 0; b < dim; b ++)
          B[j][a*dim + b] += buf[m++];
    }
  }
}


/* ---------------------------------------------------------------------- */

int ComputeRHEOSurface::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,a,b,k,m;
  int *surface = atom->surface;
  double *rsurf = atom->dvector[index_rsurf];
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      buf[m++] = divr[j];
    } else if (comm_stage == 1) {
      buf[m++] = (double) surface[j];
      buf[m++] = rsurf[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOSurface::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, a, b, m, last;
  int *surface = atom->surface;
  double *rsurf = atom->dvector[index_rsurf];

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      divr[i] = buf[m++];
    } else if (comm_stage == 1) {
      surface[i] = (int) buf[m++];
      rsurf[i] = buf[m++];
    }
  }
}
