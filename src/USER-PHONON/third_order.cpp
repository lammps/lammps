//
// Created by charlie sievers on 7/5/18.
//

#include "third_order.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "force.h"
#include "memory.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "finish.h"
#include "math_special.h"
#include <algorithm>
#include <complex>

using namespace LAMMPS_NS;
using namespace MathSpecial;
enum{REGULAR,BALLISTICO};

/* ---------------------------------------------------------------------- */

ThirdOrder::ThirdOrder(LAMMPS *lmp) : Pointers(lmp), fp(NULL)
{
  external_force_clear = 1;
}

/* ---------------------------------------------------------------------- */

ThirdOrder::~ThirdOrder()
{
  if (fp && me == 0) fclose(fp);
  fp = NULL;
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

  // compute all forces
  external_force_clear = 0;
  eflag=0;
  vflag=0;
  update_force();

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
  else if (strcmp(arg[1],"eskm") == 0) style = BALLISTICO;
  else error->all(FLERR,"Illegal Dynamical Matrix command");

  // set option defaults

  binaryflag = 0;
  scaleflag = 0;
  compressed = 0;
  file_flag = 0;
  file_opened = 0;
  conversion = 1;

  // read options from end of input line
  if (style == REGULAR) options(narg-3,&arg[3]);  //COME BACK
  else if (style == BALLISTICO) options(narg-3,&arg[3]); //COME BACK
  else if (comm->me == 0 && screen) fprintf(screen,"Illegal Dynamical Matrix command\n");
  del = force->numeric(FLERR, arg[2]);

  if (atom->map_style == 0)
    error->all(FLERR,"third_order command requires an atom map, see atom_modify");

  // move atoms by 3-vector or specified variable(s)

  if (style == REGULAR) {
    setup();
    timer->init();
    timer->barrier_start();
    calculateMatrix();
    timer->barrier_stop();
  }

  if (style == BALLISTICO) {
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
  if (narg < 0) error->all(FLERR,"Illegal third_order command");
  int iarg = 0;
  const char *filename = "third_order.dat";
  std::stringstream fss;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal third_order command");
      fss << arg[iarg + 1];
      filename = fss.str().c_str();
      file_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"binary") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal third_order command");
      if (strcmp(arg[iarg+1],"gzip") == 0) {
        compressed = 1;
      } else if (strcmp(arg[iarg+1],"yes") == 0) {
        binaryflag = 1;
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal third_order command");
  }
  if (file_flag == 1 and me == 0) {
    openfile(filename);
  }
}

/* ----------------------------------------------------------------------
   generic opening of a file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void ThirdOrder::openfile(const char* filename)
{
  // if file already opened, return
  if (file_opened) return;

  if (compressed) {
#ifdef LAMMPS_GZIP
    char gzip[128];
    sprintf(gzip,"gzip -6 > %s",filename);
#ifdef _WIN32
    fp = _popen(gzip,"wb");
#else
    fp = popen(gzip,"w");
#endif
#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  } else if (binaryflag) {
    fp = fopen(filename,"wb");
  } else {
    fp = fopen(filename,"w");
  }

  if (fp == NULL) error->one(FLERR,"Cannot open dump file");

  file_opened = 1;
}

/* ----------------------------------------------------------------------
   create dynamical matrix
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

  double *dynmat = new double[3*dynlen];
  double *fdynmat = new double[3*dynlen];
  memset(&dynmat[0],0,dynlen*sizeof(double));
  memset(&fdynmat[0],0,dynlen*sizeof(double));

  if (comm->me == 0 && screen) {
    fprintf(screen,"Calculating Third Order ...\n");
    fprintf(screen,"  Total # of atoms = " BIGINT_FORMAT "\n", natoms);
    fprintf(screen,"  Atoms in group = " BIGINT_FORMAT "\n", gcount);
    fprintf(screen,"  Total third order elements = "
            BIGINT_FORMAT "\n", (dynlen*dynlen*dynlen) );
  }

  update->nsteps = 0;
  int prog = 0;
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    for (int alpha=0; alpha<3; alpha++){
      for (bigint j=1; j<=natoms; j++){
        local_jdx = atom->map(j);
        for (int beta=0; beta<3; beta++){
          displace_atom(local_idx, alpha, 1);
          displace_atom(local_jdx, beta, 1);
          update_force();
          for (bigint k=1; k<=natoms; k++){
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++){
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && gm[i-1] >= 0 && gm[j-1] >= 0 && gm[k-1] >= 0
                  && local_kdx < nlocal) {
                dynmat[gm[k-1]*3+gamma] += f[local_kdx][gamma];
              }
            }
          }
          displace_atom(local_jdx, beta, -2);
          update_force();
          for (bigint k=1; k<=natoms; k++){
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++){
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && gm[i-1] >= 0 && gm[j-1] >= 0 && gm[k-1] >= 0
                  && local_kdx < nlocal) {
                dynmat[gm[k-1]*3+gamma] -= f[local_kdx][gamma];
              }
            }
          }
          displace_atom(local_jdx, beta, 1);
          displace_atom(local_idx,alpha,-2);
          displace_atom(local_jdx, beta, 1);
          update_force();
          for (bigint k=1; k<=natoms; k++){
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++){
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && gm[i-1] >= 0 && gm[j-1] >= 0 && gm[k-1] >= 0
                  && local_kdx < nlocal) {
                dynmat[gm[k-1]*3+gamma] -= f[local_kdx][gamma];
              }
            }
          }
          displace_atom(local_jdx, beta, -2);
          update_force();
          for (bigint k=1; k<=natoms; k++){
            local_kdx = atom->map(k);
            for (int gamma=0; gamma<3; gamma++){
              if (local_idx >= 0 && local_jdx >= 0 && local_kdx >= 0
                  && gm[i-1] >= 0 && gm[j-1] >= 0 && gm[k-1] >= 0
                  && local_kdx < nlocal) {
                dynmat[gm[k-1]*3+gamma] += f[local_kdx][gamma];
                dynmat[gm[k-1]*3+gamma] /= (4 * del * del);
              }
            }
          }
          displace_atom(local_jdx, beta, 1);
          displace_atom(local_idx, alpha, 1);
          MPI_Reduce(dynmat,fdynmat,3*dynlen,MPI_DOUBLE,MPI_SUM,0,world);
          if (me == 0){
            writeMatrix(fdynmat, gm[i-1], alpha, gm[j-1], beta);
          }
          memset(&dynmat[0],0,dynlen*sizeof(double));
        }
      }
    }
    if (comm->me == 0 && screen) {
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

  if (screen && me ==0 )
    fprintf(screen,"Finished Calculating Third Order Tensor\n");
}

/* ----------------------------------------------------------------------
   write dynamical matrix
------------------------------------------------------------------------- */

void ThirdOrder::writeMatrix(double *dynmat, bigint i, int a, bigint j, int b)
{
  if (me != 0)
    return;

  double norm;
  if (!binaryflag && fp) {
    clearerr(fp);
    for (int k = 0; k < gcount; k++){
      norm = square(dynmat[k*3])+
        square(dynmat[k*3+1])+
        square(dynmat[k*3+2]);
      if (norm > 1.0e-16)
        fprintf(fp,
                BIGINT_FORMAT " %d " BIGINT_FORMAT " %d " BIGINT_FORMAT
                " %7.8f %7.8f %7.8f\n",
                i+1, a + 1, j+1, b + 1, groupmap[k]+1,
                dynmat[k*3] * conversion,
                dynmat[k*3+1] * conversion,
                dynmat[k*3+2] * conversion);
    }
  } else if (binaryflag && fp){
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

  while (sametag[j] >= 0){
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
  force_clear();

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (atom->molecular) {
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
  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
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
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
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
    if (comm->me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
    conv_energy = 6.022E22; // J -> 10 J/mol
    conv_mass = 6.022E26; // kg -> g/mol
    conv_distance = 1E-10; // meter -> angstrom

  } else if (strcmp(style,"cgs") == 0) {
    if (comm->me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
    conv_energy = 6.022E12; // Erg -> 10 J/mol
    conv_mass = 6.022E23; // g -> g/mol
    conv_distance = 1E-7; // centimeter -> angstrom

  } else if (strcmp(style,"electron") == 0) {
    conv_energy = 262550; // Hartree -> 10 J/mol
    conv_mass = 1; // amu -> g/mol
    conv_distance = 0.529177249; // bohr -> angstrom

  } else if (strcmp(style,"micro") == 0) {
    if (comm->me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
    conv_energy = 6.022E10; // picogram-micrometer^2/microsecond^2 -> 10 J/mol
    conv_mass = 6.022E11; // pg -> g/mol
    conv_distance = 1E-4; // micrometer -> angstrom

  } else if (strcmp(style,"nano") == 0) {
    if (comm->me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
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
  bigint *temp_groupmap = new bigint[natoms];

  //find number of local atoms in the group (final_gid)
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal) && mask[local_idx] & groupbit)
      gid += 1; // gid at the end of loop is final_Gid
  }
  //create an array of length final_gid
  bigint *sub_groupmap = new bigint[gid];

  gid = 0;
  //create a map between global atom id and group atom id for each proc
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal)
        && (mask[local_idx] & groupbit)){
      sub_groupmap[gid] = i;
      gid += 1;
    }
  }

  //populate arrays for Allgatherv
  for (int i=0; i<comm->nprocs; i++){
    recv[i] = 0;
  }
  recv[comm->me] = gid;
  MPI_Allreduce(recv,displs,comm->nprocs,MPI_INT,MPI_SUM,world);
  for (int i=0; i<comm->nprocs; i++){
    recv[i]=displs[i];
    if (i>0) displs[i] = displs[i-1]+recv[i-1];
    else displs[i] = 0;
  }

  //combine subgroup maps into total temporary groupmap
  MPI_Allgatherv(sub_groupmap,gid,MPI_LMP_BIGINT,
                 temp_groupmap,recv,displs,MPI_LMP_BIGINT,world);
  std::sort(temp_groupmap,temp_groupmap+gcount);

  //populate member groupmap based on temp groupmap
  bigint j = 0;
  for (bigint i=1; i<=natoms; i++){
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
