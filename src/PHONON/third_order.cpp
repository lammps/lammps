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
   Contributing author: Charlie Sievers (UC Davis), charliesievers at cox.net
------------------------------------------------------------------------- */

#include "third_order.h"

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
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

#include <cstring>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathSpecial;
enum{REGULAR,ESKM};

/* ---------------------------------------------------------------------- */

ThirdOrder::ThirdOrder(LAMMPS *lmp) : Command(lmp), fp(nullptr)
{
  external_force_clear = 1;
}

/* ---------------------------------------------------------------------- */

ThirdOrder::~ThirdOrder()
{
  if (fp && me == 0) {
    if (compressed) platform::pclose(fp);
    else fclose(fp);
  }
  fp = nullptr;
  memory->destroy(groupmap);
}

/* ----------------------------------------------------------------------
   setup without output or one-time post-init setup
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void ThirdOrder::setup()
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

  // build neighbor list this command needs based on earlier request

  neighbor->build_one(list);

  // compute all forces
  external_force_clear = 0;
  eflag=0;
  vflag=0;
  if (force->kspace) {
    force->kspace->setup();
  }
  update_force();

  modify->setup(vflag);
  update->setupflag = 0;

  if (gcount == atom->natoms)
    for (bigint i=0; i<atom->natoms; i++)
      groupmap[i] = i;
  else
    create_groupmap();
}

/* ---------------------------------------------------------------------- */

void ThirdOrder::command(int narg, char **arg)
{
  MPI_Comm_rank(world,&me);

  if (domain->box_exist == 0)
    error->all(FLERR,"third_order command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal third_order command");

  // request a full neighbor list for use by this command

  neighbor->add_request(this, "third_order", NeighConst::REQ_FULL);

  lmp->init();
  list = neighbor->find_list(this);

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;

  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find third_order group ID");
  groupbit = group->bitmask[igroup];
  gcount = group->count(igroup);
  dynlen = (gcount)*3;
  memory->create(groupmap,atom->natoms,"total_group_map:totalgm");
  update->setupflag = 1;

  int style = -1;
  if (strcmp(arg[1],"regular") == 0) style = REGULAR;
  else if (strcmp(arg[1],"eskm") == 0) style = ESKM;
  else error->all(FLERR,"Illegal Dynamical Matrix command");

  // set option defaults

  binaryflag = 0;
  scaleflag = 0;
  compressed = 0;
  file_flag = 0;
  file_opened = 0;
  conversion = 1;
  folded = 0;

  // set Neigborlist attributes to NULL
  ijnum = nullptr;
  neighbortags = nullptr;

  // read options from end of input line
  if (style == REGULAR) options(narg-3,&arg[3]);
  else if (style == ESKM) options(narg-3,&arg[3]);
  else error->all(FLERR,"Illegal Third Order command");
  del = utils::numeric(FLERR, arg[2],false,lmp);

  if (!folded) dynlenb = dynlen;
  else dynlenb = (atom->natoms)*3;

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Third Order command requires an atom map, see atom_modify");

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
    conversion = conv_energy/conv_distance/conv_distance;
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

void ThirdOrder::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal Third Order command");
  int iarg = 0;
  const char *filename = "Third Order.dat";

  while (iarg < narg) {
    if (strcmp(arg[iarg],"binary") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal Third Order command");
      if (strcmp(arg[iarg+1],"gzip") == 0) {
        compressed = 1;
      } else {
        binaryflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal third_order command");
      filename = arg[iarg + 1];
      file_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"fold") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal Third Order command");
      if (strcmp(arg[iarg+1],"yes") == 0) {
        folded = 1;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
        folded = 0;
      } else error->all(FLERR,"Illegal input for Third Order fold option");
      iarg += 2;
    } else error->all(FLERR,"Illegal Third Order command");
  }
  if (file_flag == 1 && me == 0) {
    openfile(filename);
  }
}

/* ----------------------------------------------------------------------
   generic opening of a file
   ASCII or binary or compressed
   some derived classes override this function
------------------------------------------------------------------------- */

void ThirdOrder::openfile(const char* filename)
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
    if (!fp) error->one(FLERR,"Cannot open third_order file: {}", utils::getsyserror());
  }


  file_opened = 1;
}

/* ----------------------------------------------------------------------
   create third order tensor
------------------------------------------------------------------------- */

void ThirdOrder::calculateMatrix()
{
  int local_idx; // local index
  int local_jdx; // second local index
  int local_kdx; // third local index
  int nlocal = atom->nlocal;
  bigint natoms = atom->natoms;
  bigint *gm = groupmap;
  double **f = atom->f;
  int inum;
  bigint j;
  bigint *firstneigh;

  auto dynmat = new double[dynlenb];
  auto fdynmat = new double[dynlenb];
  memset(&dynmat[0],0,dynlenb*sizeof(double));
  memset(&fdynmat[0],0,dynlenb*sizeof(double));

  getNeighbortags();

  if (comm->me == 0 && screen) {
    fputs("Calculating Third Order ...\n", screen);
    fmt::print(screen,"  Total # of atoms = {}\n"
                      "  Atoms in group = {}\n"
                      "  Total third order elements = {}\n",
                      natoms, gcount, dynlen*dynlen*dynlen);
  }

  update->nsteps = 0;
  int prog = 0;
  for (bigint i=1; i<=natoms; i++){
    if (gm[i-1] < 0)
      continue;
    inum = ijnum[i-1];
    firstneigh = neighbortags[i-1];
    local_idx = atom->map(i);
    for (int alpha=0; alpha<3; alpha++){
      for (int jj=0; jj<inum; jj++){
        j = firstneigh[jj];
        if (gm[j] < 0 && !folded)
          continue;
        local_jdx = atom->map(j+1);
        for (int beta=0; beta<3; beta++){
          displace_atom(local_idx, alpha, 1);
          displace_atom(local_jdx, beta, 1);
          update_force();
          for (bigint k=1; k<=natoms; k++) {
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++) {
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && ((gm[j] >= 0 && gm[k-1] >= 0) || folded)
                  && local_kdx < nlocal) {
                if (folded) {
                  dynmat[(k-1)*3+gamma] += f[local_kdx][gamma];
                } else {
                  dynmat[gm[k-1]*3+gamma] += f[local_kdx][gamma];
                }
              }
            }
          }
          displace_atom(local_jdx, beta, -2);
          update_force();
          for (bigint k=1; k<=natoms; k++) {
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++) {
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && ((gm[j] >= 0 && gm[k-1] >= 0) || folded)
                  && local_kdx < nlocal) {
                if (folded) {
                  dynmat[(k-1)*3+gamma] -= f[local_kdx][gamma];
                } else {
                  dynmat[gm[k-1]*3+gamma] -= f[local_kdx][gamma];
                }
              }
            }
          }
          displace_atom(local_jdx, beta, 1);
          displace_atom(local_idx,alpha,-2);
          displace_atom(local_jdx, beta, 1);
          update_force();
          for (bigint k=1; k<=natoms; k++) {
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++) {
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && ((gm[j] >= 0 && gm[k-1] >= 0) || folded)
                  && local_kdx < nlocal) {
                if (folded) {
                  dynmat[(k-1)*3+gamma] -= f[local_kdx][gamma];
                } else {
                  dynmat[gm[k-1]*3+gamma] -= f[local_kdx][gamma];
                }
              }
            }
          }
          displace_atom(local_jdx, beta, -2);
          update_force();
          for (bigint k=1; k<=natoms; k++) {
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++) {
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && ((gm[j] >= 0 && gm[k-1] >= 0) || folded)
                  && local_kdx < nlocal) {
                if (folded) {
                  dynmat[(k-1)*3+gamma] += f[local_kdx][gamma];
                  dynmat[(k-1)*3+gamma] /= (4 * del * del);
                } else {
                  dynmat[gm[k-1]*3+gamma] += f[local_kdx][gamma];
                  dynmat[gm[k-1]*3+gamma] /= (4 * del * del);
                }
              }
            }
          }
          displace_atom(local_jdx, beta, 1);
          displace_atom(local_idx, alpha, 1);
          MPI_Reduce(dynmat,fdynmat,dynlenb,MPI_DOUBLE,MPI_SUM,0,world);
          if (me == 0){
            if (folded) {
              writeMatrix(fdynmat, gm[i-1], alpha, j, beta);
            } else {
              writeMatrix(fdynmat, gm[i-1], alpha, gm[j], beta);
            }
          }
          memset(&dynmat[0],0,dynlenb*sizeof(double));
        }
      }
    }
    if (me == 0 && screen) {
      int p = 10 * gm[i-1] / gcount;
      if (p > prog) {
        prog = p;
        fprintf(screen," %d%%",p*10);
        fflush(screen);
      }
    }
  }

  delete [] dynmat;
  delete [] fdynmat;

  if (screen && me ==0)
    fprintf(screen,"Finished Calculating Third Order Tensor\n");
}

/* ----------------------------------------------------------------------
   write third order tensor
------------------------------------------------------------------------- */

void ThirdOrder::writeMatrix(double *dynmat, bigint i, int a, bigint j, int b)
{
  if (me != 0)
    return;

  double norm;
  if (!binaryflag && fp) {
    clearerr(fp);
    if (folded){
      for (int k = 0; k < atom->natoms; k++){
        norm = square(dynmat[k*3])+square(dynmat[k*3+1])+square(dynmat[k*3+2]);
        if (norm > 1.0e-16)
          fmt::print(fp, "{} {} {} {} {} {:17.8f} {:17.8f} {:17.8f}\n",
                     i+1, a+1, j+1, b+1, k+1, dynmat[k*3] * conversion,
                     dynmat[k*3+1] * conversion, dynmat[k*3+2] * conversion);
      }
    } else {
      for (int k = 0; k < gcount; k++){
        norm = square(dynmat[k*3])+square(dynmat[k*3+1])+square(dynmat[k*3+2]);
        if (norm > 1.0e-16)
          fmt::print(fp, "{} {} {} {} {} {:17.8f} {:17.8f} {:17.8f}\n",
                     i+1, a+1, j+1, b+1, groupmap[k]+1, dynmat[k*3] * conversion,
                     dynmat[k*3+1] * conversion, dynmat[k*3+2] * conversion);
      }
    }
  } else if (binaryflag && fp) {
    clearerr(fp);
    fwrite(&dynmat[0], sizeof(double), dynlen, fp);
  }
  if (ferror(fp)) error->one(FLERR,"Error writing to file");

}

/* ----------------------------------------------------------------------
  Displace atoms
   ---------------------------------------------------------------------- */

void ThirdOrder::displace_atom(int local_idx, int direction, int magnitude)
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

void ThirdOrder::update_force()
{
  neighbor->ago = 0;
  if (modify->get_fix_by_id("package_intel")) neighbor->decide();
  force_clear();
  int n_post_force = modify->n_post_force_any;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;

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

void ThirdOrder::force_clear()
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

/* ---------------------------------------------------------------------- */

void ThirdOrder::convert_units(const char *style)
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

void ThirdOrder::create_groupmap()
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
    if ((local_idx >= 0) && (local_idx < nlocal)
        && (mask[local_idx] & groupbit)) {
      sub_groupmap[gid] = i;
      gid += 1;
    }
  }

  //populate arrays for Allgatherv
  for (int i=0; i<comm->nprocs; i++) {
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

/* ---------------------------------------------------------------------- */

void ThirdOrder::getNeighbortags() {
  // Create an extended neighbor list which is indexed by atom tag and yields atom tags
  // groupmap[global atom index-1] = global atom indices (-1) of extended neighbors

  bigint natoms = atom->natoms;
  int *ilist,*jlist,*numneigh,**firstneigh;
  bigint *Jlist,*klist;
  int ii,jj,kk,inum,jnum,knum,sum;
  int *temptags = (int*) malloc(natoms*sizeof(int));
  int *ijnumproc = (int*) malloc(natoms*sizeof(int));
  memory->create(ijnum, natoms, "thirdorder:ijnum");
  bigint **firsttags;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  memset(&ijnum[0],0,natoms*sizeof(int));
  for (ii = 0; ii < inum; ii++) {
    sum = 0;
    memset(&temptags[0],0,natoms*sizeof(int));
    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    temptags[atom->tag[ilist[ii] & NEIGHMASK]-1] = 1;
    for (jj = 0; jj < jnum; jj++) {
      temptags[atom->tag[jlist[jj] & NEIGHMASK]-1] = 1;
    }
    for (bigint i=0; i<=natoms-1; i++) {
      sum += temptags[i];
    }
    ijnum[atom->tag[ilist[ii] & NEIGHMASK]-1] = sum;
  }
  MPI_Allreduce(ijnum,ijnumproc,natoms,MPI_INT,MPI_SUM,world);
  memset(&ijnum[0],0,natoms*sizeof(int));
  sum = 0;
  for (bigint i=0; i<=natoms-1; i++) {
    sum += ijnumproc[i];
  }

  bigint nbytes = ((bigint) sizeof(bigint)) * sum;
  auto data = (bigint *) memory->smalloc(nbytes, "thirdorder:firsttags");
  auto datarecv = (bigint *) memory->smalloc(nbytes, "thirdorder:neighbortags");
  nbytes = ((bigint) sizeof(bigint *)) * natoms;
  firsttags = (bigint **) memory->smalloc(nbytes, "thirdorder:firsttags");
  neighbortags = (bigint **) memory->smalloc(nbytes, "thirdorder:neighbortags");
  memset(&data[0],0,sum*sizeof(bigint));
  memset(&datarecv[0],0,sum*sizeof(bigint));

  bigint n = 0;
  for (bigint i = 0; i < natoms; i++) {
    firsttags[i] = &data[n];
    neighbortags[i] = &datarecv[n];
    n += ijnumproc[i];
  }

  for (ii = 0; ii < inum; ii++) {
    int m = 0;
    memset(&temptags[0],0,natoms*sizeof(int));
    jnum = numneigh[ii];
    jlist = firstneigh[ii];
    temptags[atom->tag[ilist[ii] & NEIGHMASK]-1] = 1;
    for (jj = 0; jj < jnum; jj++) {
      temptags[atom->tag[jlist[jj] & NEIGHMASK]-1] = 1;
    }
    for (int j=0; j < natoms; j++) {
      if (temptags[j] == 1) {
        neighbortags[atom->tag[ilist[ii] & NEIGHMASK]-1][m] = j;
        m += 1;
      }
    }
  }
  MPI_Allreduce(datarecv,data,sum,MPI_LMP_BIGINT,MPI_SUM,world);

  for (bigint i = 0; i < natoms; i++) {
    ijnum[i] = 0;
    sum = 0;
    memset(&temptags[0],0,natoms*sizeof(int));
    jnum = ijnumproc[i];
    Jlist = firsttags[i];
    temptags[i] = 1;
    for (jj = 0; jj < jnum; jj++) {
      temptags[Jlist[jj]] = 1;
      klist = firsttags[Jlist[jj]];
      knum = ijnumproc[Jlist[jj]];
      for (kk = 0; kk < knum; kk++) {
        temptags[klist[kk]] = 1;
      }
    }
    for (bigint j=0; j<natoms; j++)
      sum += temptags[j];

    ijnum[i] = sum;
  }

  sum = 0;
  for (bigint i=0; i<=natoms-1; i++) {
    sum += ijnum[i];
  }

  free (neighbortags);
  nbytes = ((bigint) sizeof(bigint)) * sum;
  datarecv = (bigint *) memory->smalloc(nbytes, "thirdorder:firsttags");
  nbytes = ((bigint) sizeof(bigint *)) * natoms;
  neighbortags = (bigint **) memory->smalloc(nbytes, "thirdorder:neighbortags");
  memset(&datarecv[0],0,sum*sizeof(bigint));

  n = 0;
  for (bigint i = 0; i < natoms; i++) {
    neighbortags[i] = &datarecv[n];
    n += ijnum[i];
  }

  for (bigint i = 0; i < natoms; i++) {
    int m = 0;
    memset(&temptags[0],0,natoms*sizeof(int));
    jnum = ijnumproc[i];
    Jlist = firsttags[i];
    temptags[i] = 1;
    for (int j = 0; j < jnum; j++) {
      temptags[Jlist[j]] = 1;
      klist = firsttags[Jlist[j]];
      knum = ijnumproc[Jlist[j]];
      for (kk = 0; kk < knum; kk++) {
        temptags[klist[kk]] = 1;
      }
    }
    for (bigint j=0; j < natoms; j++) {
      if (temptags[j] == 1) {
        neighbortags[i][m] = j;
        m += 1;
      }
    }
  }

  free (firsttags);
  free (ijnumproc);
  free (temptags);
}
