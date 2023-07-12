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
   Contributing author: Charlie Sievers (UC Davis), charliesievers at cox.net
------------------------------------------------------------------------- */

#include "dynamical_matrix.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <algorithm>

using namespace LAMMPS_NS;
enum{REGULAR,ESKM};

/* ---------------------------------------------------------------------- */

DynamicalMatrix::DynamicalMatrix(LAMMPS *lmp) : Command(lmp), fp(nullptr)
{
  external_force_clear = 0;
}

/* ---------------------------------------------------------------------- */

DynamicalMatrix::~DynamicalMatrix()
{
  if (fp && comm->me == 0) {
    if (compressed) platform::pclose(fp);
    else fclose(fp);
    memory->destroy(groupmap);
    fp = nullptr;
  }
}

/* ----------------------------------------------------------------------
   setup without output or one-time post-init setup
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void DynamicalMatrix::setup()
{
  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  domain->image_check();
  domain->box_too_small_check();
  neighbor->build(1);

  // compute all forces
  eflag=0;
  vflag=0;
  if (force->kspace) {
    force->kspace->setup();
  }
  update_force();


  update->setupflag = 0;

  //if all then skip communication groupmap population
  if (gcount == atom->natoms)
    for (bigint i=0; i<atom->natoms; i++)
      groupmap[i] = i;
  else
    create_groupmap();
}

/* ---------------------------------------------------------------------- */

void DynamicalMatrix::command(int narg, char **arg)
{
  MPI_Comm_rank(world,&me);

  if (domain->box_exist == 0)
    error->all(FLERR,"Dynamical_matrix command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal dynamical_matrix command");

  lmp->init();

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;

  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find dynamical matrix group ID");
  groupbit = group->bitmask[igroup];
  gcount = group->count(igroup);
  dynlen = (gcount)*3;
  memory->create(groupmap,atom->natoms,"total_group_map:totalgm");
  update->setupflag = 1;

  int style = -1;
  if (strcmp(arg[1],"regular") == 0) style = REGULAR;
  else if (strcmp(arg[1],"eskm") == 0) style = ESKM;
  else error->all(FLERR,"Illegal Dynamical_matrix command");
  del = utils::numeric(FLERR, arg[2],false,lmp);

  // set option defaults

  binaryflag = 0;
  scaleflag = 0;
  compressed = 0;
  file_flag = 0;
  file_opened = 0;
  folded = 0;
  conversion = 1;

  // read options from end of input line
  if (style == REGULAR) options(narg-3,&arg[3]);
  else if (style == ESKM) options(narg-3,&arg[3]);
  else if (me == 0 && screen) fprintf(screen,"Illegal Dynamical Matrix command\n");

  if (!folded) dynlenb = dynlen;
  else dynlenb = (atom->natoms)*3;

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Dynamical_matrix command requires an atom map");

  // move atoms by 3-vector or specified variable(s)

  if (style == REGULAR) {
    setup();
    timer->init();
    timer->barrier_start();
    calculateMatrix();
    timer->barrier_stop();
  }

  if (style == ESKM) {
    setup();
    convert_units(update->unit_style);
    conversion = conv_energy/conv_distance/conv_mass;
    timer->init();
    timer->barrier_start();
    calculateMatrix();
    timer->barrier_stop();
  }

  Finish finish(lmp);
  finish.end(1);
}

/* ----------------------------------------------------------------------
   parse optional parameters
------------------------------------------------------------------------- */

void DynamicalMatrix::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal dynamical_matrix command");
  int iarg = 0;
  const char* filename = "dynmat.dyn";

  while (iarg < narg) {
    if (strcmp(arg[iarg],"binary") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal dynamical_matrix command");
      if (strcmp(arg[iarg+1],"gzip") == 0) {
          compressed = 1;
      } else {
          binaryflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dynamical_matrix command");
      filename = arg[iarg + 1];
      file_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"fold") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal dynamical_matrix command");
      if (strcmp(arg[iarg+1],"yes") == 0) {
        folded = 1;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
        folded = 0;
      } else error->all(FLERR,"Illegal input for dynamical_matrix fold option");
      iarg += 2;
    } else error->all(FLERR,"Illegal dynamical_matrix command");
  }
  if (file_flag == 1) {
    openfile(filename);
  }
}

/* ----------------------------------------------------------------------
   generic opening of a file
   ASCII or binary or compressed
   some derived classes override this function
------------------------------------------------------------------------- */

void DynamicalMatrix::openfile(const char *filename)
{
  // if file already opened, return
  if (file_opened) return;
  fp = nullptr;

  if (me == 0) {
    if (compressed) {
      fp = platform::compressed_write(std::string(filename)+".gz");
      if (!fp) error->one(FLERR,"Cannot open compressed file");
    } else if (binaryflag) {
      fp = fopen(filename,"wb");
    } else {
      fp = fopen(filename,"w");
    }
    if (!fp) error->one(FLERR,"Cannot open dynmat file: {}", utils::getsyserror());
  }

  file_opened = 1;
}

/* ----------------------------------------------------------------------
   create dynamical matrix
------------------------------------------------------------------------- */

void DynamicalMatrix::calculateMatrix()
{
  int local_idx; // local index
  int local_jdx; // second local index
  int nlocal = atom->nlocal;
  bigint natoms = atom->natoms;
  int *type = atom->type;
  bigint *gm = groupmap;
  double imass; // dynamical matrix element
  double *m = atom->mass;
  double **f = atom->f;

  auto dynmat = new double*[3];
  for (int i=0; i<3; i++)
    dynmat[i] = new double[dynlenb];

  auto fdynmat = new double*[3];
  for (int i=0; i<3; i++)
    fdynmat[i] = new double[dynlenb];

  //initialize dynmat to all zeros
  dynmat_clear(dynmat);

  if (me == 0 && screen) {
    fputs("Calculating Dynamical Matrix ...\n", screen);
    fmt::print(screen,"  Total # of atoms = {}\n"
                      "  Atoms in group = {}\n"
                      "  Total dynamical matrix elements = {}\n",
               natoms, gcount, dynlen*dynlen);
  }

  // emit dynlen rows of dimalpha*dynlen*dimbeta elements

  update->nsteps = 0;
  int prog = 0;
  for (bigint i=1; i<=natoms; i++) {
    local_idx = atom->map(i);
    if (gm[i-1] < 0)
      continue;
    for (int alpha=0; alpha<3; alpha++) {
      displace_atom(local_idx, alpha, 1);
      update_force();
      for (bigint j=1; j<=natoms; j++) {
        local_jdx = atom->map(j);
        if (local_idx >= 0 && local_jdx >= 0 && local_jdx < nlocal
          && (gm[j-1] >= 0 || folded)){
          if (folded) {
            for (int beta=0; beta<3; beta++){
              dynmat[alpha][(j-1)*3+beta] -= f[local_jdx][beta];
            }
          } else {
            for (int beta=0; beta<3; beta++){
              dynmat[alpha][gm[j-1]*3+beta] -= f[local_jdx][beta];
            }
          }
        }
      }
      displace_atom(local_idx,alpha,-2);
      update_force();
      for (bigint j=1; j<=natoms; j++) {
        local_jdx = atom->map(j);
        if (local_idx >= 0 && local_jdx >= 0 && local_jdx < nlocal
          && (gm[j-1] >= 0 || folded)){
          if (atom->rmass_flag == 1)
            imass = sqrt(m[local_idx] * m[local_jdx]);
          else
            imass = sqrt(m[type[local_idx]] * m[type[local_jdx]]);
          if (folded){
            for (int beta=0; beta<3; beta++){
              dynmat[alpha][(j-1)*3+beta] -= -f[local_jdx][beta];
              dynmat[alpha][(j-1)*3+beta] /= (2 * del * imass);
              dynmat[alpha][(j-1)*3+beta] *= conversion;
            }
          } else {
            for (int beta=0; beta<3; beta++){
              dynmat[alpha][gm[j-1]*3+beta] -= -f[local_jdx][beta];
              dynmat[alpha][gm[j-1]*3+beta] /= (2 * del * imass);
              dynmat[alpha][gm[j-1]*3+beta] *= conversion;
            }
          }
        }
      }
      displace_atom(local_idx,alpha,1);
    }
    for (int k=0; k<3; k++)
      MPI_Reduce(dynmat[k],fdynmat[k],dynlenb,MPI_DOUBLE,MPI_SUM,0,world);
    if (me == 0)
      writeMatrix(fdynmat);
    dynmat_clear(dynmat);
    if (me == 0 && screen) {
      int p = 10 * gm[i-1] / gcount;
      if (p > prog) {
        prog = p;
        fprintf(screen," %d%%",p*10);
        fflush(screen);
      }
    }
  }
  if (me == 0 && screen) fprintf(screen,"\n");

  for (int i=0; i < 3; i++)
    delete [] dynmat[i];
  delete [] dynmat;

  for (int i=0; i < 3; i++)
    delete [] fdynmat[i];
  delete [] fdynmat;

  if (screen && me == 0) fprintf(screen,"Finished Calculating Dynamical Matrix\n");
}

/* ----------------------------------------------------------------------
   write dynamical matrix
------------------------------------------------------------------------- */

void DynamicalMatrix::writeMatrix(double **dynmat)
{
  if (me != 0 || !fp)
    return;

  clearerr(fp);
  if (binaryflag) {
    for (int i=0; i<3; i++)
      fwrite(dynmat[i], sizeof(double), dynlenb, fp);
    if (ferror(fp))
      error->one(FLERR, "Error writing to binary file");
  } else {
    for (int i = 0; i < 3; i++) {
      for (bigint j = 0; j < dynlenb; j++) {
        if ((j+1)%3==0) fprintf(fp, "%4.8f\n", dynmat[i][j]);
        else fprintf(fp, "%4.8f ",dynmat[i][j]);
      }
    }
    if (ferror(fp))
      error->one(FLERR,"Error writing to file");
  }
}

/* ----------------------------------------------------------------------
  Displace atoms
   ---------------------------------------------------------------------- */

void DynamicalMatrix::displace_atom(int local_idx, int direction, int magnitude)
{
  if (local_idx < 0) return;

  double **x = atom->x;
  int *sametag = atom->sametag;
  int j = local_idx;

  x[local_idx][direction] += del*magnitude;

  while (sametag[j] >= 0) {
    j = sametag[j];
    x[j][direction] += del*magnitude;
  }
}


/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

void DynamicalMatrix::update_force()
{
  neighbor->ago = 0;
  if (modify->get_fix_by_id("package_intel")) neighbor->decide();
  force_clear();
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force_any;

  if (n_pre_force) {
    modify->pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
    timer->stamp(Timer::BOND);
  }
  if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }
  if (n_pre_reverse) {
    modify->pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }
  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }
  // force modifications

  if (n_post_force) {
    modify->post_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  ++ update->nsteps;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void DynamicalMatrix::force_clear()
{
  if (external_force_clear) return;

  // clear global force array
  // if either newton flag is set, also include ghosts

  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&atom->f[0][0],0,3*nbytes);
  }
}

/* ----------------------------------------------------------------------
   clear dynmat needed
------------------------------------------------------------------------- */

void DynamicalMatrix::dynmat_clear(double **dynmat)
{
  size_t nbytes = sizeof(double) * dynlenb;

  if (nbytes) {
    for (int i=0; i<3; i++)
      memset(&dynmat[i][0],0,nbytes);
  }
}

/* ---------------------------------------------------------------------- */

void DynamicalMatrix::convert_units(const char *style)
{
  // physical constants from:
  // https://physics.nist.gov/cuu/Constants/Table/allascii.txt
  // using thermochemical calorie = 4.184 J

  if (strcmp(style,"lj") == 0) {
    error->all(FLERR,"Conversion Not Set");
    //conversion = 1; // lj -> 10 J/mol

  } else if (strcmp(style,"real") == 0) {
    conv_energy = 418.4; // kcal/mol -> 10 J/mol
    conv_mass = 1; // g/mol -> g/mol
    conv_distance = 1; // angstrom -> angstrom

  } else if (strcmp(style,"metal") == 0) {
    conv_energy = 9648.5; // eV -> 10 J/mol
    conv_mass = 1; // g/mol -> g/mol
    conv_distance = 1; // angstrom -> angstrom

  } else if (strcmp(style,"si") == 0) {
    if (me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
    conv_energy = 6.022E22; // J -> 10 J/mol
    conv_mass = 6.022E26; // kg -> g/mol
    conv_distance = 1E-10; // meter -> angstrom

  } else if (strcmp(style,"cgs") == 0) {
    if (me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
    conv_energy = 6.022E12; // Erg -> 10 J/mol
    conv_mass = 6.022E23; // g -> g/mol
    conv_distance = 1E-7; // centimeter -> angstrom

  } else if (strcmp(style,"electron") == 0) {
    conv_energy = 262550; // Hartree -> 10 J/mol
    conv_mass = 1; // amu -> g/mol
    conv_distance = 0.529177249; // bohr -> angstrom

  } else if (strcmp(style,"micro") == 0) {
    if (me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
    conv_energy = 6.022E10; // picogram-micrometer^2/microsecond^2 -> 10 J/mol
    conv_mass = 6.022E11; // pg -> g/mol
    conv_distance = 1E-4; // micrometer -> angstrom

  } else if (strcmp(style,"nano") == 0) {
    if (me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
    conv_energy = 6.022E4; // attogram-nanometer^2/nanosecond^2 -> 10 J/mol
    conv_mass = 6.022E5; // ag -> g/mol
    conv_distance = 0.1; // angstrom -> angstrom

  } else error->all(FLERR,"Units Type Conversion Not Found");

}

/* ---------------------------------------------------------------------- */

void DynamicalMatrix::create_groupmap()
{
  //Create a group map which maps atom order onto group
  // groupmap[global atom index-1] = output column/row

  int local_idx; // local index
  int gid = 0; //group index
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  bigint natoms = atom->natoms;
  int *recv = new int[comm->nprocs];
  int *displs = new int[comm->nprocs];
  auto temp_groupmap = new bigint[natoms];

  //find number of local atoms in the group (final_gid)
  for (bigint i=1; i<=natoms; i++) {
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal) && mask[local_idx] & groupbit)
      gid += 1; // gid at the end of loop is final_Gid
  }
  //create an array of length final_gid
  auto sub_groupmap = new bigint[gid];

  gid = 0;
  //create a map between global atom id and group atom id for each proc
  for (bigint i=1; i<=natoms; i++) {
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal) && mask[local_idx] & groupbit) {
      sub_groupmap[gid] = i;
      gid += 1;
    }
  }

  //populate arrays for Allgatherv
  for (int i=0; i < comm->nprocs; i++) {
    recv[i] = 0;
  }
  recv[me] = gid;
  MPI_Allreduce(recv,displs,comm->nprocs,MPI_INT,MPI_SUM,world);
  for (int i=0; i<comm->nprocs; i++) {
    recv[i]=displs[i];
    if (i>0) displs[i] = displs[i-1]+recv[i-1];
    else displs[i] = 0;
  }

  //combine subgroup maps into total temporary groupmap
  MPI_Allgatherv(sub_groupmap,gid,MPI_LMP_BIGINT,temp_groupmap,recv,displs,MPI_LMP_BIGINT,world);
  std::sort(temp_groupmap,temp_groupmap+gcount);

  //populate member groupmap based on temp groupmap
  bigint j = 0;
  for (bigint i=1; i<=natoms; i++) {
    // flag groupmap contents that are in temp_groupmap
    if (j < gcount && i == temp_groupmap[j])
      groupmap[i-1] = j++;
    else
      groupmap[i-1] = -1;
  }

  //free that memory!
  delete[] recv;
  delete[] displs;
  delete[] sub_groupmap;
  delete[] temp_groupmap;
}
