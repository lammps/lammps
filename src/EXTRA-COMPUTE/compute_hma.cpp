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
This compute implements harmonically-mapped averaging for crystalline solids.
The current implementation handles atomic crystals.

Computing the heat capacity relies on being able to compute the second
derivative of the energy with respect to atom positions.  This capability is
provided by the single2 method in Pair, but is currently only implemented for
the shifted-force LJ potential (lj/smooth/linear).  Pair::single2 takes a single
pair and (like Pair::single) returns the energy and sets the force as an out
parameter, but also sets the elements of 6-element double array out parameter,
which are the unique components of the atomic Hessian tensor for the pair.  A
helper method exists (Pair::pairTensor), which will compute the tensor from
linear derivatives and the vector between the positions.  HMA Heat capacity can
be computed for other models by implementing single2 in those Pair classes.

More information about HMA is available in these publications:

A. J. Schultz, D. A. Kofke, “Comprehensive high-precision high-accuracy
equation of state and coexistence properties for classical Lennard-Jones
crystals and low-temperature fluid phases”, J. Chem. Phys. 149, 204508 (2018)
https://doi.org/10.1063/1.5053714

S. G. Moustafa, A. J. Schultz, D. A. Kofke, “Harmonically Assisted Methods for
Computing the Free Energy of Classical Crystals by Molecular Simulation: A
Comparative Study”, J. Chem. Theory Comput. 13, 825-834 (2017)
https://doi.org/10.1021/acs.jctc.6b01082

S. G. Moustafa, A. J. Schultz, D. A. Kofke, “Very fast averaging of thermal
properties of crystals by molecular simulation”, Phys. Rev. E 92, 043303 (2015)
https://doi.org/10.1103/PhysRevE.92.043303
------------------------------------------------------------------------- */

#include "compute_hma.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_store.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHMA::ComputeHMA(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), id_temp(nullptr), deltaR(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute hma command");
  if (igroup) error->all(FLERR,"Compute hma must use group all");
  if (strcmp(arg[3],"NULL") == 0)
    error->all(FLERR,"fix ID specifying the set temperature of canonical simulation is required");
  else id_temp = utils::strdup(arg[3]);

  create_attribute = 1;
  extscalar = 1;
  timeflag = 1;

  // (from compute displace/atom) create a new fix STORE style
  // our new fix's id (id_fix)= compute-ID + COMPUTE_STORE
  // our new fix's group = same as compute group

  id_fix = utils::strdup(std::string(id)+"_COMPUTE_STORE");
  fix = (FixStore *)modify->add_fix(fmt::format("{} {} STORE peratom 1 3",
                                                id_fix, group->names[igroup]));

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;
    double **x = atom->x;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      domain->unmap(x[i],image[i],xoriginal[i]);
  }

  vector_flag = 1;
  extvector = -1;
  comm_forward = 0; // 3 if 2nd derivative needed

  computeU = computeP = computeCv = -1;
  returnAnharmonic = 0;
  size_vector = 0;
  extlist = new int[3];
  for (int iarg=4; iarg<narg; iarg++) {
    if (!strcmp(arg[iarg], "u")) {
      if (computeU>-1) continue;
      computeU = size_vector;
      extlist[size_vector] = 1;
      size_vector++;
    } else if (!strcmp(arg[iarg], "p")) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute hma command");
      if (computeP>-1) continue;
      computeP = size_vector;
      deltaPcap = utils::numeric(FLERR, arg[iarg+1],false,lmp);
      extlist[size_vector] = 0;
      size_vector++;
      iarg++;
    } else if (!strcmp(arg[iarg], "cv")) {
      if (computeCv>-1) continue;
      computeCv = size_vector;
      comm_forward = 3;
      extlist[size_vector] = 1;
      size_vector++;
    } else if (!strcmp(arg[iarg], "anharmonic")) {
      // the first time we're called, we'll grab lattice pressure and energy
      returnAnharmonic = -1;
    } else {
      error->all(FLERR,"Illegal compute hma command");
    }
  }

  if (size_vector == 0) error->all(FLERR,"Illegal compute hma command");
  vector = new double[size_vector];

  if (computeU>-1 || computeCv>-1) peflag = 1;
  if (computeP>-1) pressflag = 1;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeHMA::~ComputeHMA()
{
  // check nfix in case all fixes have already been deleted
  if (modify->nfix) modify->delete_fix(id_fix);

  delete[] id_fix;
  delete[] id_temp;
  delete[] extlist;
  delete[] vector;

  memory->destroy(deltaR);
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::init() {
  if (computeCv>-1) {
    if (force->pair == nullptr)
      error->all(FLERR,"No pair style is defined for compute hma cv");
    if (force->pair->single_enable == 0)
      error->all(FLERR,"Pair style does not support compute hma cv");
  }

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

void ComputeHMA::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

void ComputeHMA::setup()
{
  int dummy=0;
  int ifix = modify->find_fix(id_temp);
  if (ifix < 0) error->all(FLERR,"Could not find compute hma temperature ID");
  double * temperat = (double *) modify->fix[ifix]->extract("t_target",dummy);
  if (temperat==nullptr) error->all(FLERR,"Could not find compute hma temperature ID");
  finaltemp = * temperat;

  // set fix which stores original atom coords

  int ifix2 = modify->find_fix(id_fix);
  if (ifix2 < 0) error->all(FLERR,"Could not find hma store fix ID");
  fix = (FixStore *) modify->fix[ifix2];
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::compute_vector()
{
  invoked_vector = update->ntimestep;

  // grow deltaR array if necessary
  if (comm_forward>0 && atom->nmax > nmax) {
    memory->destroy(deltaR);
    nmax = atom->nmax;
    memory->create(deltaR,nmax,3,"hma:deltaR");
  }

  double **xoriginal = fix->astore;
  double fdr = 0.0;
  double **x = atom->x;
  double **f = atom->f;

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double u = 0.0;
  if (computeU>-1 || computeCv>-1) {
    if (force->pair) u += force->pair->eng_vdwl + force->pair->eng_coul;
    if (force->bond) u += force->bond->energy;
    if (force->angle) u += force->angle->energy;
    if (force->dihedral) u += force->dihedral->energy;
    if (force->improper) u += force->improper->energy;
  }

  int dimension = domain->dimension;
  double p = 0, vol = 0;
  if (computeP>-1) {
    p = virial_compute(dimension);
    vol = xprd * yprd;
    if (dimension == 3) vol *= zprd;
    p *= force->nktv2p / (dimension*vol);
    if (returnAnharmonic == -1) {
      pLat = p;
    }
  }

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++) {
      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      double dx = x[i][0] + xbox*xprd - xoriginal[i][0];
      double dy = x[i][1] + ybox*yprd - xoriginal[i][1];
      double dz = x[i][2] + zbox*zprd - xoriginal[i][2];
      if (comm_forward>0) {
        deltaR[i][0] = dx;
        deltaR[i][1] = dy;
        deltaR[i][2] = dz;
      }
      fdr += dx*f[i][0] + dy*f[i][1] + dz*f[i][2];
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      int xbox = (image[i] & IMGMASK) - IMGMAX;
      int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image[i] >> IMG2BITS) - IMGMAX;
      double dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
      double dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
      double dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
      if (comm_forward>0) {
        deltaR[i][0] = dx;
        deltaR[i][1] = dy;
        deltaR[i][2] = dz;
      }
      fdr += ((dx*f[i][0])+(dy*f[i][1])+(dz*f[i][2]));
    }
  }

  double phiSum = 0.0;
  if (computeCv>-1) {
    comm->forward_comm_compute(this);
    double** cutsq = force->pair->cutsq;
    if (force->pair) {
      double **x = atom->x;
      int *type = atom->type;
      int nlocal = atom->nlocal;
      double *special_lj = force->special_lj;
      double *special_coul = force->special_coul;
      int newton_pair = force->newton_pair;

      if (update->firststep == update->ntimestep) neighbor->build_one(list,1);
      else neighbor->build_one(list);
      int inum = list->inum;
      int *ilist = list->ilist;
      int *numneigh = list->numneigh;
      int **firstneigh = list->firstneigh;

      for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii];
        double fac = (newton_pair || i < nlocal) ? 1.0 : 0.5;
        double* ix = x[i];
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        double *idr = deltaR[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          if (!newton_pair && j>=nlocal) fac -= 0.5;
          double factor_lj = special_lj[sbmask(j)];
          double factor_coul = special_coul[sbmask(j)];
          j &= NEIGHMASK;
          double* jx = x[j];
          double delr[3];
          delr[0] = ix[0] - jx[0];
          delr[1] = ix[1] - jx[1];
          delr[2] = ix[2] - jx[2];
          double rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];
          int jtype = type[j];
          if (rsq < cutsq[itype][jtype]) {
            double* jdr = deltaR[j];
            double fforce, d2u[6];
            force->pair->single_hessian(i, j, itype, jtype, rsq, delr, factor_coul, factor_lj, fforce, d2u);
            int m = 0;
            for (int k=0; k<3; k++) {
              double a = fac;
              for (int l=k; l<3; l++) {
                phiSum += a*(idr[k]*jdr[l]+jdr[k]*idr[l])*d2u[m];
                phiSum -= a*(idr[k]*idr[l]*d2u[m] + jdr[k]*jdr[l]*d2u[m]);
                m++;
                if (k==l) a *= 2;
              }
            }
          }
        }
      }
    }
  }

  // compute and sum up properties across processors

  double fdrTotal;
  MPI_Allreduce(&fdr,&fdrTotal,1,MPI_DOUBLE,MPI_SUM,world);
  double uTotal;
  if (computeU>-1 || computeCv>-1) {
    MPI_Allreduce(&u,&uTotal,1,MPI_DOUBLE,MPI_SUM,world);
    if (returnAnharmonic == -1) {
      uLat = uTotal;
    }
    if (computeU>-1) {
      if (returnAnharmonic) {
        vector[computeU] = uTotal - uLat + 0.5*fdrTotal;
      } else {
        vector[computeU] = uTotal + 0.5*fdrTotal + 0.5*dimension*(atom->natoms - 1)*force->boltz*finaltemp;
      }
    }
  }

  if (computeP>-1) {
    double fv = ((deltaPcap)-(force->boltz*finaltemp*force->nktv2p*atom->natoms/vol))/(force->boltz*finaltemp*dimension*(atom->natoms - 1));
    if (returnAnharmonic) {
      vector[computeP] = p - pLat + (fv*fdrTotal);
    } else {
      vector[computeP] = p + (fv*fdrTotal) + deltaPcap;
    }
  }

  if (computeCv>-1) {
    if (computeU==-1) MPI_Allreduce(&u,&uTotal,1,MPI_DOUBLE,MPI_SUM,world);
    double buTot;
    if (returnAnharmonic) {
      buTot = (uTotal - uLat + 0.5*fdrTotal)/finaltemp;
    } else {
      buTot = (uTotal + 0.5*fdrTotal)/finaltemp + 0.5*dimension*(atom->natoms - 1)*force->boltz;
    }
    double one = -0.25*(fdr + phiSum)/finaltemp;
    double Cv;
    MPI_Allreduce(&one,&Cv,1,MPI_DOUBLE,MPI_SUM,world);
    vector[computeCv] = Cv + buTot*buTot;
    if (!returnAnharmonic) {
      vector[computeCv] += 0.5*dimension*(atom->natoms-1);
    }
  }
  if (returnAnharmonic == -1) {
    // we have our lattice properties
    returnAnharmonic = 1;
  }
}

double ComputeHMA::virial_compute(int n)
{
  double v = 0;

  // sum contributions to virial from forces and fixes

  if (force->pair) v += sumVirial(n, force->pair->virial);
  if (force->bond) v += sumVirial(n, force->bond->virial);
  if (force->angle) v += sumVirial(n, force->angle->virial);
  if (force->dihedral) v += sumVirial(n, force->dihedral->virial);
  if (force->improper) v += sumVirial(n, force->improper->virial);
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->thermo_virial) v += sumVirial(n, modify->fix[i]->virial);

  // sum virial across procs

  double virial;
  MPI_Allreduce(&v,&virial,1,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (force->kspace)
    for (int i = 0; i < n; i++) virial += force->kspace->virial[i];
  return virial;
}

/* ---------------------------------------------------------------------- */

int ComputeHMA::pack_forward_comm(int n, int *list, double *buf,
                                  int /* pbc_flag */, int * /* pbc */)
{
  int m = 0;
  for (int ii = 0; ii < n; ii++) {
    int i = list[ii];
    buf[m++] = deltaR[i][0];
    buf[m++] = deltaR[i][1];
    buf[m++] = deltaR[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeHMA::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    deltaR[i][0] = buf[m++];
    deltaR[i][1] = buf[m++];
    deltaR[i][2] = buf[m++];
  }
}


/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeHMA::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

double ComputeHMA::memory_usage()
{
  double bytes = (double)nmax * 3 * sizeof(double);
  return bytes;
}
