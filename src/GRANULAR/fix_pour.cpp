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
#include "fix_pour.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "molecule.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "region.h"
#include "region_block.h"
#include "region_cylinder.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001

enum{ATOM,MOLECULE};
enum{ONE,RANGE,POLY};

/* ---------------------------------------------------------------------- */

FixPour::FixPour(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix pour command");

  time_depend = 1;

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR,"Fix pour requires atom attributes radius, rmass");

  // required args

  ninsert = force->inumeric(FLERR,arg[3]);
  ntype = force->inumeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]);

  if (seed <= 0) error->all(FLERR,"Illegal fix pour command");

  // read options from end of input line

  options(narg-6,&arg[6]);

  // error check on type

  if (mode == ATOM && (ntype <= 0 || ntype > atom->ntypes))
    error->all(FLERR,"Invalid atom type in fix pour command");

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix pour");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix pour region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix pour region cannot be dynamic");

  if (strcmp(domain->regions[iregion]->style,"block") == 0) {
    region_style = 1;
    xlo = ((RegBlock *) domain->regions[iregion])->xlo;
    xhi = ((RegBlock *) domain->regions[iregion])->xhi;
    ylo = ((RegBlock *) domain->regions[iregion])->ylo;
    yhi = ((RegBlock *) domain->regions[iregion])->yhi;
    zlo = ((RegBlock *) domain->regions[iregion])->zlo;
    zhi = ((RegBlock *) domain->regions[iregion])->zhi;
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
        ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  } else if (strcmp(domain->regions[iregion]->style,"cylinder") == 0) {
    region_style = 2;
    char axis = ((RegCylinder *) domain->regions[iregion])->axis;
    xc = ((RegCylinder *) domain->regions[iregion])->c1;
    yc = ((RegCylinder *) domain->regions[iregion])->c2;
    rc = ((RegCylinder *) domain->regions[iregion])->radius;
    zlo = ((RegCylinder *) domain->regions[iregion])->lo;
    zhi = ((RegCylinder *) domain->regions[iregion])->hi;
    if (axis != 'z')
      error->all(FLERR,"Must use a z-axis cylinder with fix pour");
    if (xc-rc < domain->boxlo[0] || xc+rc > domain->boxhi[0] ||
        yc-rc < domain->boxlo[1] || yc+rc > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  } else error->all(FLERR,"Must use a block or cylinder region with fix pour");

  if (region_style == 2 && domain->dimension == 2)
    error->all(FLERR,
               "Must use a block region with fix pour for 2d simulations");

  // error check and further setup for mode = MOLECULE

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix_pour unless atoms have IDs");

  if (mode == MOLECULE) {
    if (onemol->xflag == 0)
      error->all(FLERR,"Fix pour molecule must have coordinates");
    if (onemol->typeflag == 0)
      error->all(FLERR,"Fix pour molecule must have atom types");
    if (ntype+onemol->maxtype <= 0 || ntype+onemol->maxtype > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix pour mol command");

    // fix pour uses geoemetric center of molecule for insertion

    onemol->compute_center();
  }

  if (rigidflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix pour rigid and not molecule");
  if (shakeflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix pour shake and not molecule");
  if (rigidflag && shakeflag)
    error->all(FLERR,"Cannot use fix pour rigid and shake");

  // setup of coords and imageflags array

  if (mode == ATOM) natom = 1;
  else natom = onemol->natoms;
  memory->create(coords,natom,4,"pour:coords");
  memory->create(imageflags,natom,"pour:imageflags");

  // find current max atom and molecule IDs if necessary

  if (idnext) find_maxid();

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // allgather arrays

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // grav = gravity in distance/time^2 units
  // assume grav = -magnitude at this point, enforce in init()

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
    if (strcmp(modify->fix[ifix]->style,"gravity/omp") == 0) break;
  }
  if (ifix == modify->nfix)
    error->all(FLERR,"No fix gravity defined for fix pour");
  grav = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;

  // nfreq = timesteps between insertions
  // should be time for a particle to fall from top of insertion region
  //   to bottom, taking into account that the region may be moving
  // set these 2 eqs equal to each other, solve for smallest positive t
  //   x = zhi + vz*t + 1/2 grav t^2
  //   x = zlo + rate*t
  //   gives t = [-(vz-rate) - sqrt((vz-rate)^2 - 2*grav*(zhi-zlo))] / grav
  //   where zhi-zlo > 0, grav < 0, and vz & rate can be either > 0 or < 0

  double v_relative,delta;
  if (domain->dimension == 3) {
    v_relative = vz - rate;
    delta = zhi - zlo;
  } else {
    v_relative = vy - rate;
    delta = yhi - ylo;
  }
  double t =
    (-v_relative - sqrt(v_relative*v_relative - 2.0*grav*delta)) / grav;
  nfreq = static_cast<int> (t/update->dt + 0.5);

  // 1st insertion on next timestep

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;

  // nper = # to insert each time
  // depends on specified volume fraction
  // volume = volume of insertion region
  // volume_one = volume of inserted particle (with max possible radius)
  // in 3d, insure dy >= 1, for quasi-2d simulations

  double volume,volume_one;
  if (domain->dimension == 3) {
    if (region_style == 1) {
      double dy = yhi - ylo;
      if (dy < 1.0) dy = 1.0;
      volume = (xhi-xlo) * dy * (zhi-zlo);
    } else volume = MY_PI*rc*rc * (zhi-zlo);
    if (mode == MOLECULE) {
      double molradius = onemol->molradius;
      volume_one = 4.0/3.0 * MY_PI * molradius*molradius*molradius;
    } else if (dstyle == ONE || dstyle == RANGE) {
      volume_one = 4.0/3.0 * MY_PI * radius_max*radius_max*radius_max;
    } else if (dstyle == POLY) {
      volume_one = 0.0;
      for (int i = 0; i < npoly; i++)
        volume_one += (4.0/3.0 * MY_PI * 
          radius_poly[i]*radius_poly[i]*radius_poly[i]) * frac_poly[i];
    }
  } else {
    volume = (xhi-xlo) * (yhi-ylo);
    if (mode == MOLECULE) {
      double molradius = onemol->molradius;
      volume_one = MY_PI * molradius*molradius;
    } else if (dstyle == ONE || dstyle == RANGE) {
      volume_one = MY_PI * radius_max*radius_max;
    } else if (dstyle == POLY) {
      volume_one = 0.0;
      for (int i = 0; i < npoly; i++)
        volume_one += (MY_PI * radius_poly[i]*radius_poly[i]) * frac_poly[i];
    }
  }

  nper = static_cast<int> (volfrac*volume/volume_one);
  int nfinal = update->ntimestep + 1 + (ninsert-1)/nper * nfreq;

  // print stats

  if (me == 0) {
    if (screen)
      fprintf(screen,
              "Particle insertion: %d every %d steps, %d by step %d\n",
              nper,nfreq,ninsert,nfinal);
    if (logfile)
      fprintf(logfile,
              "Particle insertion: %d every %d steps, %d by step %d\n",
              nper,nfreq,ninsert,nfinal);
  }
}

/* ---------------------------------------------------------------------- */

FixPour::~FixPour()
{
  delete random;
  delete [] idrigid;
  delete [] idshake;
  delete [] radius_poly;
  delete [] frac_poly;
  memory->destroy(coords);
  memory->destroy(imageflags);
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

int FixPour::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPour::init()
{
  if (domain->triclinic) 
    error->all(FLERR,"Cannot use fix pour with triclinic box");

  // insure gravity fix exists
  // for 3d must point in -z, for 2d must point in -y
  // else insertion cannot work

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
    if (strcmp(modify->fix[ifix]->style,"gravity/omp") == 0) break;
  }
  if (ifix == modify->nfix)
    error->all(FLERR,"No fix gravity defined for fix pour");

  double xgrav = ((FixGravity *) modify->fix[ifix])->xgrav;
  double ygrav = ((FixGravity *) modify->fix[ifix])->ygrav;
  double zgrav = ((FixGravity *) modify->fix[ifix])->zgrav;

  if (domain->dimension == 3) {
    if (fabs(xgrav) > EPSILON || fabs(ygrav) > EPSILON ||
        fabs(zgrav+1.0) > EPSILON)
      error->all(FLERR,"Gravity must point in -z to use with fix pour in 3d");
  } else {
    if (fabs(xgrav) > EPSILON || fabs(ygrav+1.0) > EPSILON ||
        fabs(zgrav) > EPSILON)
      error->all(FLERR,"Gravity must point in -y to use with fix pour in 2d");
  }

  double gnew = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;
  if (gnew != grav)
    error->all(FLERR,"Gravity changed since fix pour was created");

  // if rigidflag defined, check for rigid/small fix
  // its molecule template must be same as this one

  fixrigid = NULL;
  if (rigidflag) {
    int ifix = modify->find_fix(idrigid);
    if (ifix < 0) error->all(FLERR,"Fix pour rigid fix does not exist");
    fixrigid = modify->fix[ifix];
    int tmp;
    if (onemol != (Molecule *) fixrigid->extract("onemol",tmp))
      error->all(FLERR,
                 "Fix pour and fix rigid/small not using same molecule ID");
  }

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one

  fixshake = NULL;
  if (shakeflag) {
    int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix pour shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    if (onemol != (Molecule *) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix pour and fix shake not using same molecule ID");
  }
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixPour::pre_exchange()
{
  int i,j,m,flag,nlocalprev;
  double r[3],rotmat[3][3],quat[4],vnew[3];
  double *newcoord;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // find current max atom and molecule IDs if necessary

  if (!idnext) find_maxid();

  // nnew = # of particles (atoms or molecules) to insert this timestep

  int nnew = nper;
  if (ninserted + nnew > ninsert) nnew = ninsert - ninserted;

  // lo/hi current = z (or y) bounds of insertion region this timestep

  int dimension = domain->dimension;
  if (dimension == 3) {
    lo_current = zlo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = zhi + (update->ntimestep - nfirst) * update->dt * rate;
  } else {
    lo_current = ylo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = yhi + (update->ntimestep - nfirst) * update->dt * rate;
  }

  // ncount = # of my atoms that overlap the insertion region
  // nprevious = total of ncount across all procs

  int ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) ncount++;

  int nprevious;
  MPI_Allreduce(&ncount,&nprevious,1,MPI_INT,MPI_SUM,world);

  // xmine is for my atoms
  // xnear is for atoms from all procs + atoms to be inserted

  double **xmine,**xnear;
  memory->create(xmine,ncount,4,"fix_pour:xmine");
  memory->create(xnear,nprevious+nnew*natom,4,"fix_pour:xnear");
  int nnear = nprevious;

  // setup for allgatherv

  int n = 4*ncount;
  MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // load up xmine array

  double **x = atom->x;
  double *radius = atom->radius;

  ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) {
      xmine[ncount][0] = x[i][0];
      xmine[ncount][1] = x[i][1];
      xmine[ncount][2] = x[i][2];
      xmine[ncount][3] = radius[i];
      ncount++;
    }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = NULL;
  if (ncount) ptr = xmine[0];
  MPI_Allgatherv(ptr,4*ncount,MPI_DOUBLE,
                 xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  // insert new particles into xnear list, one by one
  // check against all nearby atoms and previously inserted ones
  // if there is an overlap then try again at same z (3d) or y (2d) coord
  // else insert by adding to xnear list
  // max = maximum # of insertion attempts for all particles
  // h = height, biased to give uniform distribution in time of insertion
  // for MOLECULE mode:
  //   coords = coords of all atoms in particle
  //   perform random rotation around center pt
  //   apply PBC so final coords are inside box
  //   store image flag modified due to PBC

  int success;
  double radtmp,delx,dely,delz,rsq,radsum,rn,h;
  double coord[3],xcm[3];

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  AtomVec *avec = atom->avec;
  double denstmp,vxtmp,vytmp,vztmp;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  int attempt = 0;
  int maxiter = nnew * maxattempt;
  int ntotal = nprevious + nnew*natom;

  while (nnear < ntotal) {
    rn = random->uniform();
    h = hi_current - rn*rn * (hi_current-lo_current);
    if (mode == ATOM) radtmp = radius_sample();

    success = 0;
    while (attempt < maxiter) {
      attempt++;
      xyz_random(h,coord);

      if (mode == ATOM) {
        coords[0][0] = coord[0];
        coords[0][1] = coord[1];
        coords[0][2] = coord[2];
        coords[0][3] = radtmp;
        imageflags[0] = ((imageint) IMGMAX << IMG2BITS) |
            ((imageint) IMGMAX << IMGBITS) | IMGMAX;
      } else {
        if (dimension == 3) {
          r[0] = random->uniform() - 0.5;
          r[1] = random->uniform() - 0.5;
          r[2] = random->uniform() - 0.5;
        } else {
          r[0] = r[1] = 0.0;
          r[2] = 1.0;
        }
        double theta = random->uniform() * MY_2PI;
        MathExtra::norm3(r);
        MathExtra::axisangle_to_quat(r,theta,quat);
        MathExtra::quat_to_mat(quat,rotmat);
        for (i = 0; i < natom; i++) {
          MathExtra::matvec(rotmat,onemol->dx[i],coords[i]);
          coords[i][0] += coord[0];
          coords[i][1] += coord[1];
          coords[i][2] += coord[2];

          // coords[3] = particle radius
          // default to 0.5, if radii not defined in Molecule
          //   same as atom->avec->create_atom(), invoked below

          if (onemol->radiusflag) coords[i][3] = onemol->radius[i];
          else coords[i][3] = 0.5;

          imageflags[i] = ((imageint) IMGMAX << IMG2BITS) |
            ((imageint) IMGMAX << IMGBITS) | IMGMAX;
          domain->remap(coords[i],imageflags[i]);
        }
      }

      // if any pair of atoms overlap, try again
      // use minimum_image() to account for PBC

      for (m = 0; m < natom; m++) {
        for (i = 0; i < nnear; i++) {
          delx = coords[m][0] - xnear[i][0];
          dely = coords[m][1] - xnear[i][1];
          delz = coords[m][2] - xnear[i][2];
	  domain->minimum_image(delx,dely,delz);
          rsq = delx*delx + dely*dely + delz*delz;
          radsum = coords[m][3] + xnear[i][3];
          if (rsq <= radsum*radsum) break;
        }
        if (i < nnear) break;
      }
      if (m == natom) {
        success = 1;
        break;
      }
    }

    if (!success) break;

    // proceed with insertion

    nlocalprev = atom->nlocal;

    // add all atoms in particle to xnear

    for (m = 0; m < natom; m++) {
      xnear[nnear][0] = coords[m][0];
      xnear[nnear][1] = coords[m][1];
      xnear[nnear][2] = coords[m][2];
      xnear[nnear][3] = coords[m][3];
      nnear++;
    }

    // choose random velocity for new particle
    // used for every atom in molecule
    // z velocity set to what velocity would be if particle
    //   had fallen from top of insertion region
    //   this gives continuous stream of atoms
    //   solution for v from these 2 eqs, after eliminate t:
    //     v = vz + grav*t
    //     coord[2] = hi_current + vz*t + 1/2 grav t^2

    if (dimension == 3) {
      vnew[0] = vxlo + random->uniform() * (vxhi-vxlo);
      vnew[1] = vylo + random->uniform() * (vyhi-vylo);
      vnew[2] = -sqrt(vz*vz + 2.0*grav*(coord[2]-hi_current));
    } else {
      vnew[0] = vxlo + random->uniform() * (vxhi-vxlo);
      vnew[1] = -sqrt(vy*vy + 2.0*grav*(coord[1]-hi_current));
      vnew[2] = 0.0;
    }

    // check if new atoms are in my sub-box or above it if I am highest proc
    // if so, add atom to my list via create_atom()
    // initialize additional info about the atoms
    // set group mask to "all" plus fix group

    for (m = 0; m < natom; m++) {
      if (mode == ATOM)
        denstmp = density_lo + random->uniform() * (density_hi-density_lo);
      newcoord = coords[m];

      flag = 0;
      if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
          newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
      else if (dimension == 3 && newcoord[2] >= domain->boxhi[2] &&
               comm->myloc[2] == comm->procgrid[2]-1 &&
               newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
               newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
      else if (dimension == 2 && newcoord[1] >= domain->boxhi[1] &&
               comm->myloc[1] == comm->procgrid[1]-1 &&
               newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;

      if (flag) {
        if (mode == ATOM) atom->avec->create_atom(ntype,coords[m]);
        else atom->avec->create_atom(ntype+onemol->type[m],coords[m]);
        int n = atom->nlocal - 1;
        atom->tag[n] = maxtag_all + m+1;
        if (mode == MOLECULE && atom->molecule_flag) 
          atom->molecule[n] = maxmol_all+1;
        atom->mask[n] = 1 | groupbit;
        atom->image[n] = imageflags[m];
        atom->v[n][0] = vnew[0];
        atom->v[n][1] = vnew[1];
        atom->v[n][2] = vnew[2];
        if (mode == ATOM) {
          radtmp = newcoord[3];
          atom->radius[n] = radtmp;
          atom->rmass[n] = 4.0*MY_PI/3.0 * radtmp*radtmp*radtmp * denstmp;
        } else atom->add_molecule_atom(onemol,m,n,maxtag_all);
        for (j = 0; j < nfix; j++)
          if (fix[j]->create_attribute) fix[j]->set_arrays(n);
      }
    }

    // FixRigidSmall::set_molecule stores rigid body attributes
    //   coord is new position of geometric center of mol, not COM
    // FixShake::set_molecule stores shake info for molecule

    if (rigidflag)
      fixrigid->set_molecule(nlocalprev,maxtag_all,coord,vnew,quat);
    else if (shakeflag)
      fixshake->set_molecule(nlocalprev,maxtag_all,coord,vnew,quat);

    maxtag_all += natom;
    if (mode == MOLECULE && atom->molecule_flag) maxmol_all++;
  }

  // warn if not successful with all insertions b/c too many attempts

  int ninserted_atoms = nnear - nprevious;
  int ninserted_mols = ninserted_atoms / natom;
  ninserted += ninserted_mols;
  if (ninserted_mols < nnew && me == 0)
    error->warning(FLERR,"Less insertions than requested",0);

  // reset global natoms,nbonds,etc
  // increment maxtag_all and maxmol_all if necessary
  // if global map exists, reset it now instead of waiting for comm
  // since adding atoms messes up ghosts

  if (ninserted_atoms) {
    atom->natoms += ninserted_atoms;
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
      error->all(FLERR,"Too many total atoms");
    if (mode == MOLECULE) {
      atom->nbonds += onemol->nbonds * ninserted_mols;
      atom->nangles += onemol->nangles * ninserted_mols;
      atom->ndihedrals += onemol->ndihedrals * ninserted_mols;
      atom->nimpropers += onemol->nimpropers * ninserted_mols;
    }
    if (maxtag_all >= MAXTAGINT)
      error->all(FLERR,"New atom IDs exceed maximum allowed ID");
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  // free local memory

  memory->destroy(xmine);
  memory->destroy(xnear);

  // next timestep to insert

  if (ninserted < ninsert) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   maxtag_all = current max atom ID for all atoms
   maxmol_all = current max molecule ID for all atoms
------------------------------------------------------------------------- */

void FixPour::find_maxid()
{
  tagint *tag = atom->tag;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  if (mode == MOLECULE && molecule) {
    max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,molecule[i]);
    MPI_Allreduce(&max,&maxmol_all,1,MPI_INT,MPI_MAX,world);
  }
}

/* ----------------------------------------------------------------------
   check if particle i could overlap with a particle inserted into region
   return 1 if yes, 0 if no
   for ATOM mode, use delta with maximum size for inserted atoms
   for MOLECULE mode, use delta with radius of inserted molecules
   account for PBC in overlap decision via outside() and minimum_image()
------------------------------------------------------------------------- */

int FixPour::overlap(int i)
{
  double delta;
  if (mode == ATOM) delta = atom->radius[i] + radius_max;
  else delta = atom->radius[i] + onemol->molradius;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;
  int *periodicity = domain->periodicity;

  double *x = atom->x[i];

  if (domain->dimension == 3) {
    if (region_style == 1) {
      if (outside(0,x[0],xlo-delta,xhi+delta)) return 0;
      if (outside(1,x[1],ylo-delta,yhi+delta)) return 0;
      if (outside(2,x[2],lo_current-delta,hi_current+delta)) return 0;
    } else {
      double delx = x[0] - xc;
      double dely = x[1] - yc;
      double delz = 0.0;
      domain->minimum_image(delx,dely,delz);
      double rsq = delx*delx + dely*dely;
      double r = rc + delta;
      if (rsq > r*r) return 0;
      if (outside(2,x[2],lo_current-delta,hi_current+delta)) return 0;
    }
  } else {
    if (outside(0,x[0],xlo-delta,xhi+delta)) return 0;
    if (outside(1,x[1],lo_current-delta,hi_current+delta)) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   check if value is inside/outside lo/hi bounds in dimension
   account for PBC if needed
   return 1 if value is outside, 0 if inside
------------------------------------------------------------------------- */

int FixPour::outside(int dim, double value, double lo, double hi)
{
  double boxlo = domain->boxlo[dim];
  double boxhi = domain->boxhi[dim];

  if (domain->periodicity[dim]) {
    if (lo < boxlo && hi > boxhi) {
      return 0;
    } else if (lo < boxlo) {
      if (value > hi && value < lo + domain->prd[dim]) return 1;
    } else if (hi > boxhi) {
      if (value > hi - domain->prd[dim] && value < lo) return 1;
    } else {
      if (value < lo || value > hi) return 1;
    }
  } 

  if (value < lo || value > hi) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixPour::xyz_random(double h, double *coord)
{
  if (domain->dimension == 3) {
    if (region_style == 1) {
      coord[0] = xlo + random->uniform() * (xhi-xlo);
      coord[1] = ylo + random->uniform() * (yhi-ylo);
      coord[2] = h;
    } else {
      double r1,r2;
      while (1) {
        r1 = random->uniform() - 0.5;
        r2 = random->uniform() - 0.5;
        if (r1*r1 + r2*r2 < 0.25) break;
      }
      coord[0] = xc + 2.0*r1*rc;
      coord[1] = yc + 2.0*r2*rc;
      coord[2] = h;
    }
  } else {
    coord[0] = xlo + random->uniform() * (xhi-xlo);
    coord[1] = h;
    coord[2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

double FixPour::radius_sample()
{
  if (dstyle == ONE) return radius_one;
  if (dstyle == RANGE) return radius_lo + 
                         random->uniform()*(radius_hi-radius_lo);

  double value = random->uniform();

  int i = 0;
  double sum = 0.0;
  while (sum < value) {
    sum += frac_poly[i];
    i++;
  }
  return radius_poly[i-1];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixPour::options(int narg, char **arg)
{
  // defaults

  iregion = -1;
  mode = ATOM;
  rigidflag = 0;
  idrigid = NULL;
  shakeflag = 0;
  idshake = NULL;
  idnext = 0;
  dstyle = ONE;
  radius_max = radius_one = 0.5;
  radius_poly = frac_poly = NULL;
  density_lo = density_hi = 1.0;
  volfrac = 0.25;
  maxattempt = 50;
  rate = 0.0;
  vxlo = vxhi = vylo = vyhi = vy = vz = 0.0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all(FLERR,"Fix pour region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      int imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule ID for fix pour does not exist");
      mode = MOLECULE;
      onemol = atom->molecules[imol];
      iarg += 2;
    } else if (strcmp(arg[iarg],"rigid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idrigid;
      idrigid = new char[n];
      strcpy(idrigid,arg[iarg+1]);
      rigidflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"shake") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idshake;
      idshake = new char[n];
      strcpy(idshake,arg[iarg+1]);
      shakeflag = 1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      if (strcmp(arg[iarg+1],"max") == 0) idnext = 0;
      else if (strcmp(arg[iarg+1],"next") == 0) idnext = 1;
      else error->all(FLERR,"Illegal fix pour command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diam") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      if (strcmp(arg[iarg+1],"one") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
        dstyle = ONE;
        radius_one = 0.5 * force->numeric(FLERR,arg[iarg+2]);
        radius_max = radius_one;
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"range") == 0) {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix pour command");
        dstyle = RANGE;
        radius_lo = 0.5 * force->numeric(FLERR,arg[iarg+2]);
        radius_hi = 0.5 * force->numeric(FLERR,arg[iarg+3]);
        radius_max = radius_hi;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"poly") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
        dstyle = POLY;
        npoly = force->inumeric(FLERR,arg[iarg+2]);
        if (npoly <= 0) error->all(FLERR,"Illegal fix pour command");
        if (iarg+3 + 2*npoly > narg) 
          error->all(FLERR,"Illegal fix pour command");
        radius_poly = new double[npoly];
        frac_poly = new double[npoly];
        iarg += 3;
        radius_max = 0.0;
        for (int i = 0; i < npoly; i++) {
          radius_poly[i] = 0.5 * force->numeric(FLERR,arg[iarg++]);
          frac_poly[i] = force->numeric(FLERR,arg[iarg++]);
          if (radius_poly[i] <= 0.0 || frac_poly[i] < 0.0)
            error->all(FLERR,"Illegal fix pour command");
          radius_max = MAX(radius_max,radius_poly[i]);
        }
        double sum = 0.0;
        for (int i = 0; i < npoly; i++) sum += frac_poly[i];
        if (sum != 1.0)
          error->all(FLERR,"Fix pour polydisperse fractions do not sum to 1.0");
      } else error->all(FLERR,"Illegal fix pour command");

    } else if (strcmp(arg[iarg],"dens") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
      density_lo = force->numeric(FLERR,arg[iarg+1]);
      density_hi = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pour command");
      volfrac = force->numeric(FLERR,arg[iarg+1]);
      maxattempt = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pour command");
      rate = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (domain->dimension == 3) {
        if (iarg+6 > narg) error->all(FLERR,"Illegal fix pour command");
        vxlo = force->numeric(FLERR,arg[iarg+1]);
        vxhi = force->numeric(FLERR,arg[iarg+2]);
        vylo = force->numeric(FLERR,arg[iarg+3]);
        vyhi = force->numeric(FLERR,arg[iarg+4]);
        vz = force->numeric(FLERR,arg[iarg+5]);
        iarg += 6;
      } else {
        if (iarg+4 > narg) error->all(FLERR,"Illegal fix pour command");
        vxlo = force->numeric(FLERR,arg[iarg+1]);
        vxhi = force->numeric(FLERR,arg[iarg+2]);
        vy = force->numeric(FLERR,arg[iarg+3]);
        vz = 0.0;
        iarg += 4;
      }
    } else error->all(FLERR,"Illegal fix pour command");
  }
}

/* ---------------------------------------------------------------------- */

void FixPour::reset_dt()
{
  error->all(FLERR,"Cannot change timestep with fix pour");
}

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixPour::extract(const char *str, int &itype)
{
  if (strcmp(str,"radius") == 0) {
    if (mode == ATOM) {
      if (itype == ntype) oneradius = radius_max;
      else oneradius = 0.0;
    } else {
      double *radius = onemol->radius;
      int *type = onemol->type;
      int natoms = onemol->natoms;

      // check radii of matching types in Molecule
      // default to 0.5, if radii not defined in Molecule
      //   same as atom->avec->create_atom(), invoked in pre_exchange()

      oneradius = 0.0;
      for (int i = 0; i < natoms; i++)
        if (type[i] == itype-ntype) {
          if (radius) oneradius = MAX(oneradius,radius[i]);
          else oneradius = 0.5;
        }
    }
    itype = 0;
    return &oneradius;
  }
  return NULL;
}
