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

#include "fix_rigid.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "random_mars.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "rigid_const.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

/* ---------------------------------------------------------------------- */

FixRigid::FixRigid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), step_respa(NULL),
  inpfile(NULL), nrigid(NULL), mol2body(NULL), body2mol(NULL),
  body(NULL), displace(NULL), masstotal(NULL), xcm(NULL),
  vcm(NULL), fcm(NULL), inertia(NULL), ex_space(NULL),
  ey_space(NULL), ez_space(NULL), angmom(NULL), omega(NULL),
  torque(NULL), quat(NULL), imagebody(NULL), fflag(NULL),
  tflag(NULL), langextra(NULL), sum(NULL), all(NULL),
  remapflag(NULL), xcmimage(NULL), eflags(NULL), orient(NULL),
  dorient(NULL), id_dilate(NULL), id_gravity(NULL), random(NULL),
  avec_ellipsoid(NULL), avec_line(NULL), avec_tri(NULL)
{
  int i,ibody;

  scalar_flag = 1;
  extscalar = 0;
  time_integrate = 1;
  rigid_flag = 1;
  virial_flag = 1;
  thermo_virial = 1;
  create_attribute = 1;
  dof_flag = 1;
  enforce2d_flag = 1;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  extended = orientflag = dorientflag = 0;
  body = NULL;
  xcmimage = NULL;
  displace = NULL;
  eflags = NULL;
  orient = NULL;
  dorient = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // parse args for rigid body specification
  // set nbody and body[i] for each atom

  if (narg < 4) error->all(FLERR,"Illegal fix rigid command");
  int iarg;

  mol2body = NULL;
  body2mol = NULL;

  // single rigid body
  // nbody = 1
  // all atoms in fix group are part of body

  if (strcmp(arg[3],"single") == 0) {
    rstyle = SINGLE;
    iarg = 4;
    nbody = 1;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = 0;
    }

  // each molecule in fix group is a rigid body
  // maxmol = largest molecule ID
  // ncount = # of atoms in each molecule (have to sum across procs)
  // nbody = # of non-zero ncount values
  // use nall as incremented ptr to set body[] values for each atom

  } else if (strcmp(arg[3],"molecule") == 0 || strcmp(arg[3],"custom") == 0) {
    rstyle = MOLECULE;
    tagint *molecule;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int custom_flag = strcmp(arg[3],"custom") == 0;
    if (custom_flag) {
      if (narg < 5) error->all(FLERR,"Illegal fix rigid command");

      // determine whether atom-style variable or atom property is used

      if (strstr(arg[4],"i_") == arg[4]) {
        int is_double=0;
        int custom_index = atom->find_custom(arg[4]+2,is_double);
        if (custom_index == -1)
          error->all(FLERR,"Fix rigid custom requires "
                     "previously defined property/atom");
        else if (is_double)
          error->all(FLERR,"Fix rigid custom requires "
                     "integer-valued property/atom");
        int minval = INT_MAX;
        int *value = atom->ivector[custom_index];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) minval = MIN(minval,value[i]);
        int vmin = minval;
        MPI_Allreduce(&vmin,&minval,1,MPI_INT,MPI_MIN,world);
        molecule = new tagint[nlocal];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit)
            molecule[i] = (tagint)(value[i] - minval + 1);
          else
            molecule[i] = 0;

      } else if (strstr(arg[4],"v_") == arg[4]) {
        int ivariable = input->variable->find(arg[4]+2);
        if (ivariable < 0)
          error->all(FLERR,"Variable name for fix rigid custom does not exist");
        if (input->variable->atomstyle(ivariable) == 0)
          error->all(FLERR,"Fix rigid custom variable is no atom-style variable");
        double *value = new double[nlocal];
        input->variable->compute_atom(ivariable,0,value,1,0);
        int minval = INT_MAX;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) minval = MIN(minval,(int)value[i]);
        int vmin = minval;
        MPI_Allreduce(&vmin,&minval,1,MPI_INT,MPI_MIN,world);
        molecule = new tagint[nlocal];
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit)
            molecule[i] = (tagint)((tagint)value[i] - minval + 1);
        delete[] value;

      } else error->all(FLERR,"Unsupported fix rigid custom property");

    } else {
      if (atom->molecule_flag == 0)
        error->all(FLERR,"Fix rigid molecule requires atom attribute molecule");
      molecule = atom->molecule;
    }

    iarg = 4 + custom_flag;

    tagint maxmol_tag = -1;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) maxmol_tag = MAX(maxmol_tag,molecule[i]);

    tagint itmp;
    MPI_Allreduce(&maxmol_tag,&itmp,1,MPI_LMP_TAGINT,MPI_MAX,world);
    if (itmp+1 > MAXSMALLINT)
      error->all(FLERR,"Too many molecules for fix rigid");
    maxmol = (int) itmp;

    int *ncount;
    memory->create(ncount,maxmol+1,"rigid:ncount");
    for (i = 0; i <= maxmol; i++) ncount[i] = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) ncount[molecule[i]]++;

    memory->create(mol2body,maxmol+1,"rigid:mol2body");
    MPI_Allreduce(ncount,mol2body,maxmol+1,MPI_INT,MPI_SUM,world);

    nbody = 0;
    for (i = 0; i <= maxmol; i++)
      if (mol2body[i]) mol2body[i] = nbody++;
      else mol2body[i] = -1;

    memory->create(body2mol,nbody,"rigid:body2mol");

    nbody = 0;
    for (i = 0; i <= maxmol; i++)
      if (mol2body[i] >= 0) body2mol[nbody++] = i;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = mol2body[molecule[i]];
    }

    memory->destroy(ncount);
    if (custom_flag) delete [] molecule;

  // each listed group is a rigid body
  // check if all listed groups exist
  // an atom must belong to fix group and listed group to be in rigid body
  // error if atom belongs to more than 1 rigid body

  } else if (strcmp(arg[3],"group") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal fix rigid command");
    rstyle = GROUP;
    nbody = force->inumeric(FLERR,arg[4]);
    if (nbody <= 0) error->all(FLERR,"Illegal fix rigid command");
    if (narg < 5+nbody) error->all(FLERR,"Illegal fix rigid command");
    iarg = 5+nbody;

    int *igroups = new int[nbody];
    for (ibody = 0; ibody < nbody; ibody++) {
      igroups[ibody] = group->find(arg[5+ibody]);
      if (igroups[ibody] == -1)
        error->all(FLERR,"Could not find fix rigid group ID");
    }

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit)
        for (ibody = 0; ibody < nbody; ibody++)
          if (mask[i] & group->bitmask[igroups[ibody]]) {
            if (body[i] >= 0) flag = 1;
            body[i] = ibody;
          }
    }

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall)
      error->all(FLERR,"One or more atoms belong to multiple rigid bodies");

    delete [] igroups;

  } else error->all(FLERR,"Illegal fix rigid command");

  // error check on nbody

  if (nbody == 0) error->all(FLERR,"No rigid bodies defined");

  // create all nbody-length arrays

  memory->create(nrigid,nbody,"rigid:nrigid");
  memory->create(masstotal,nbody,"rigid:masstotal");
  memory->create(xcm,nbody,3,"rigid:xcm");
  memory->create(vcm,nbody,3,"rigid:vcm");
  memory->create(fcm,nbody,3,"rigid:fcm");
  memory->create(inertia,nbody,3,"rigid:inertia");
  memory->create(ex_space,nbody,3,"rigid:ex_space");
  memory->create(ey_space,nbody,3,"rigid:ey_space");
  memory->create(ez_space,nbody,3,"rigid:ez_space");
  memory->create(angmom,nbody,3,"rigid:angmom");
  memory->create(omega,nbody,3,"rigid:omega");
  memory->create(torque,nbody,3,"rigid:torque");
  memory->create(quat,nbody,4,"rigid:quat");
  memory->create(imagebody,nbody,"rigid:imagebody");
  memory->create(fflag,nbody,3,"rigid:fflag");
  memory->create(tflag,nbody,3,"rigid:tflag");
  memory->create(langextra,nbody,6,"rigid:langextra");

  memory->create(sum,nbody,6,"rigid:sum");
  memory->create(all,nbody,6,"rigid:all");
  memory->create(remapflag,nbody,4,"rigid:remapflag");

  // initialize force/torque flags to default = 1.0
  // for 2d: fz, tx, ty = 0.0

  array_flag = 1;
  size_array_rows = nbody;
  size_array_cols = 15;
  global_freq = 1;
  extarray = 0;

  for (i = 0; i < nbody; i++) {
    fflag[i][0] = fflag[i][1] = fflag[i][2] = 1.0;
    tflag[i][0] = tflag[i][1] = tflag[i][2] = 1.0;
    if (domain->dimension == 2) fflag[i][2] = tflag[i][0] = tflag[i][1] = 0.0;
  }

  // number of linear rigid bodies is counted later

  nlinear = 0;

  // parse optional args

  int seed;
  langflag = 0;
  reinitflag = 1;

  tstat_flag = 0;
  pstat_flag = 0;
  allremap = 1;
  t_chain = 10;
  t_iter = 1;
  t_order = 3;
  p_chain = 10;

  inpfile = NULL;
  id_gravity = NULL;
  id_dilate = NULL;

  pcouple = NONE;
  pstyle = ANISO;
  dimension = domain->dimension;

  for (int i = 0; i < 3; i++) {
    p_start[i] = p_stop[i] = p_period[i] = 0.0;
    p_flag[i] = 0;
  }

  while (iarg < narg) {
    if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix rigid command");

      int mlo,mhi;
      force->bounds(FLERR,arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(arg[iarg+3],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");

      if (domain->dimension == 2 && zflag == 1.0)
        error->all(FLERR,"Fix rigid z force cannot be on for 2d simulation");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
        fflag[m-1][0] = xflag;
        fflag[m-1][1] = yflag;
        fflag[m-1][2] = zflag;
        count++;
      }
      if (count == 0) error->all(FLERR,"Illegal fix rigid command");

      iarg += 5;

    } else if (strcmp(arg[iarg],"torque") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix rigid command");

      int mlo,mhi;
      force->bounds(FLERR,arg[iarg+1],nbody,mlo,mhi);

      double xflag,yflag,zflag;
      if (strcmp(arg[iarg+2],"off") == 0) xflag = 0.0;
      else if (strcmp(arg[iarg+2],"on") == 0) xflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(arg[iarg+3],"off") == 0) yflag = 0.0;
      else if (strcmp(arg[iarg+3],"on") == 0) yflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(arg[iarg+4],"off") == 0) zflag = 0.0;
      else if (strcmp(arg[iarg+4],"on") == 0) zflag = 1.0;
      else error->all(FLERR,"Illegal fix rigid command");

      if (domain->dimension == 2 && (xflag == 1.0 || yflag == 1.0))
        error->all(FLERR,"Fix rigid xy torque cannot be on for 2d simulation");

      int count = 0;
      for (int m = mlo; m <= mhi; m++) {
        tflag[m-1][0] = xflag;
        tflag[m-1][1] = yflag;
        tflag[m-1][2] = zflag;
        count++;
      }
      if (count == 0) error->all(FLERR,"Illegal fix rigid command");

      iarg += 5;

    } else if (strcmp(arg[iarg],"langevin") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid") != 0 && strcmp(style,"rigid/nve") != 0 &&
          strcmp(style,"rigid/omp") != 0 && strcmp(style,"rigid/nve/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      langflag = 1;
      t_start = force->numeric(FLERR,arg[iarg+1]);
      t_stop = force->numeric(FLERR,arg[iarg+2]);
      t_period = force->numeric(FLERR,arg[iarg+3]);
      seed = force->inumeric(FLERR,arg[iarg+4]);
      if (t_period <= 0.0)
        error->all(FLERR,"Fix rigid langevin period must be > 0.0");
      if (seed <= 0) error->all(FLERR,"Illegal fix rigid command");
      iarg += 5;

    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/nvt") != 0 && strcmp(style,"rigid/npt") != 0 &&
          strcmp(style,"rigid/nvt/omp") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      tstat_flag = 1;
      t_start = force->numeric(FLERR,arg[iarg+1]);
      t_stop = force->numeric(FLERR,arg[iarg+2]);
      t_period = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"iso") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      pcouple = XYZ;
      p_start[0] = p_start[1] = p_start[2] = force->numeric(FLERR,arg[iarg+1]);
      p_stop[0] = p_stop[1] = p_stop[2] = force->numeric(FLERR,arg[iarg+2]);
      p_period[0] = p_period[1] = p_period[2] =
        force->numeric(FLERR,arg[iarg+3]);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
              p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg],"aniso") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      p_start[0] = p_start[1] = p_start[2] = force->numeric(FLERR,arg[iarg+1]);
      p_stop[0] = p_stop[1] = p_stop[2] = force->numeric(FLERR,arg[iarg+2]);
      p_period[0] = p_period[1] = p_period[2] =
        force->numeric(FLERR,arg[iarg+3]);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
              p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      p_start[0] = force->numeric(FLERR,arg[iarg+1]);
      p_stop[0] = force->numeric(FLERR,arg[iarg+2]);
      p_period[0] = force->numeric(FLERR,arg[iarg+3]);
      p_flag[0] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      p_start[1] = force->numeric(FLERR,arg[iarg+1]);
      p_stop[1] = force->numeric(FLERR,arg[iarg+2]);
      p_period[1] = force->numeric(FLERR,arg[iarg+3]);
      p_flag[1] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      p_start[2] = force->numeric(FLERR,arg[iarg+1]);
      p_stop[2] = force->numeric(FLERR,arg[iarg+2]);
      p_period[2] = force->numeric(FLERR,arg[iarg+3]);
      p_flag[2] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(arg[iarg+1],"xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1],"xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1],"yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1],"xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1],"none") == 0) pcouple = NONE;
      else error->all(FLERR,"Illegal fix rigid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix rigid npt/nph command");
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else {
        allremap = 0;
        delete [] id_dilate;
        int n = strlen(arg[iarg+1]) + 1;
        id_dilate = new char[n];
        strcpy(id_dilate,arg[iarg+1]);
        int idilate = group->find(id_dilate);
        if (idilate == -1)
          error->all(FLERR,
                     "Fix rigid npt/nph dilate group ID does not exist");
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"tparam") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/nvt") != 0 && strcmp(style,"rigid/npt") != 0 &&
          strcmp(style,"rigid/nvt/omp") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      t_chain = force->inumeric(FLERR,arg[iarg+1]);
      t_iter = force->inumeric(FLERR,arg[iarg+2]);
      t_order = force->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"pchain") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp(style,"rigid/npt") != 0 && strcmp(style,"rigid/nph") != 0 &&
          strcmp(style,"rigid/npt/omp") != 0 &&
          strcmp(style,"rigid/nph/omp") != 0)
        error->all(FLERR,"Illegal fix rigid command");
      p_chain = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"infile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid command");
      delete [] inpfile;
      int n = strlen(arg[iarg+1]) + 1;
      inpfile = new char[n];
      strcpy(inpfile,arg[iarg+1]);
      restart_file = 1;
      reinitflag = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"reinit") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid command");
      if (strcmp("yes",arg[iarg+1]) == 0) reinitflag = 1;
      else if  (strcmp("no",arg[iarg+1]) == 0) reinitflag = 0;
      else error->all(FLERR,"Illegal fix rigid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"gravity") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix rigid command");
      delete [] id_gravity;
      int n = strlen(arg[iarg+1]) + 1;
      id_gravity = new char[n];
      strcpy(id_gravity,arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal fix rigid command");
  }

  // set pstat_flag

  pstat_flag = 0;
  for (int i = 0; i < 3; i++)
    if (p_flag[i]) pstat_flag = 1;

  if (pcouple == XYZ || (dimension == 2 && pcouple == XY)) pstyle = ISO;
  else pstyle = ANISO;

  // initialize Marsaglia RNG with processor-unique seed

  if (langflag) random = new RanMars(lmp,seed + me);
  else random = NULL;

  // initialize vector output quantities in case accessed before run

  for (i = 0; i < nbody; i++) {
    xcm[i][0] = xcm[i][1] = xcm[i][2] = 0.0;
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    fcm[i][0] = fcm[i][1] = fcm[i][2] = 0.0;
    torque[i][0] = torque[i][1] = torque[i][2] = 0.0;
  }

  // nrigid[n] = # of atoms in Nth rigid body
  // error if one or zero atoms

  int *ncount = new int[nbody];
  for (ibody = 0; ibody < nbody; ibody++) ncount[ibody] = 0;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) ncount[body[i]]++;

  MPI_Allreduce(ncount,nrigid,nbody,MPI_INT,MPI_SUM,world);
  delete [] ncount;

  for (ibody = 0; ibody < nbody; ibody++)
    if (nrigid[ibody] <= 1) error->all(FLERR,"One or zero atoms in rigid body");

  // wait to setup bodies until first init() using current atom properties

  setupflag = 0;

  // compute per body forces and torques at final_integrate() by default

  earlyflag = 0;

  // print statistics

  int nsum = 0;
  for (ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];

  if (me == 0) {
    if (screen) fprintf(screen,"%d rigid bodies with %d atoms\n",nbody,nsum);
    if (logfile) fprintf(logfile,"%d rigid bodies with %d atoms\n",nbody,nsum);
  }
}

/* ---------------------------------------------------------------------- */

FixRigid::~FixRigid()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete random;
  delete [] inpfile;
  delete [] id_dilate;
  delete [] id_gravity;

  memory->destroy(mol2body);
  memory->destroy(body2mol);

  // delete locally stored per-atom arrays

  memory->destroy(body);
  memory->destroy(xcmimage);
  memory->destroy(displace);
  memory->destroy(eflags);
  memory->destroy(orient);
  memory->destroy(dorient);

  // delete nbody-length arrays

  memory->destroy(nrigid);
  memory->destroy(masstotal);
  memory->destroy(xcm);
  memory->destroy(vcm);
  memory->destroy(fcm);
  memory->destroy(inertia);
  memory->destroy(ex_space);
  memory->destroy(ey_space);
  memory->destroy(ez_space);
  memory->destroy(angmom);
  memory->destroy(omega);
  memory->destroy(torque);
  memory->destroy(quat);
  memory->destroy(imagebody);
  memory->destroy(fflag);
  memory->destroy(tflag);
  memory->destroy(langextra);

  memory->destroy(sum);
  memory->destroy(all);
  memory->destroy(remapflag);
}

/* ---------------------------------------------------------------------- */

int FixRigid::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  if (langflag) mask |= POST_FORCE;
  mask |= PRE_NEIGHBOR;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigid::init()
{
  int i,ibody;

  triclinic = domain->triclinic;

  // atom style pointers to particles that store extra info

  avec_ellipsoid = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  avec_line = (AtomVecLine *) atom->style_match("line");
  avec_tri = (AtomVecTri *) atom->style_match("tri");

  // warn if more than one rigid fix
  // if earlyflag, warn if any post-force fixes come after a rigid fix

  int count = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) count++;
  if (count > 1 && me == 0) error->warning(FLERR,"More than one fix rigid");

  if (earlyflag) {
    int rflag = 0;
    for (i = 0; i < modify->nfix; i++) {
      if (modify->fix[i]->rigid_flag) rflag = 1;
      if (rflag && (modify->fmask[i] & POST_FORCE) &&
          !modify->fix[i]->rigid_flag) {
        char str[128];
        snprintf(str,128,"Fix %s alters forces after fix rigid",
                 modify->fix[i]->id);
        error->warning(FLERR,str);
      }
    }
  }

  // warn if body properties are read from inpfile
  //   and the gravity keyword is not set and a gravity fix exists
  // this could mean body particles are overlapped
  //   and gravity is not applied correctly

  if (inpfile && !id_gravity) {
    for (i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style,"gravity") == 0) {
        if (comm->me == 0)
          error->warning(FLERR,"Gravity may not be correctly applied "
                         "to rigid bodies if they consist of "
                         "overlapped particles");
        break;
      }
    }
  }

  // error if npt,nph fix comes before rigid fix

  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) break;
    if (strcmp(modify->fix[i]->style,"nph") == 0) break;
  }

  if (i < modify->nfix) {
    for (int j = i; j < modify->nfix; j++)
      if (strcmp(modify->fix[j]->style,"rigid") == 0)
        error->all(FLERR,"Rigid fix must come before NPT/NPH fix");
  }

  // add gravity forces based on gravity vector from fix

  if (id_gravity) {
    int ifix = modify->find_fix(id_gravity);
    if (ifix < 0) error->all(FLERR,"Fix rigid cannot find fix gravity ID");
    if (strcmp(modify->fix[ifix]->style,"gravity") != 0)
      error->all(FLERR,"Fix rigid gravity fix is invalid");
    int tmp;
    gvec = (double *) modify->fix[ifix]->extract("gvec",tmp);
  }

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;

  // setup rigid bodies, using current atom info. if reinitflag is not set,
  // do the initialization only once, b/c properties may not be re-computable
  // especially if overlapping particles.
  //   do not do dynamic init if read body properties from inpfile.
  // this is b/c the inpfile defines the static and dynamic properties and may
  // not be computable if contain overlapping particles.
  //   setup_bodies_static() reads inpfile itself

  if (reinitflag || !setupflag) {
    setup_bodies_static();
    if (!inpfile) setup_bodies_dynamic();
    setupflag = 1;
  }

  // temperature scale factor

  double ndof = 0.0;
  for (ibody = 0; ibody < nbody; ibody++) {
    ndof += fflag[ibody][0] + fflag[ibody][1] + fflag[ibody][2];
    ndof += tflag[ibody][0] + tflag[ibody][1] + tflag[ibody][2];
  }
  ndof -= nlinear;
  if (ndof > 0.0) tfactor = force->mvv2e / (ndof * force->boltz);
  else tfactor = 0.0;
}

/* ----------------------------------------------------------------------
   invoke pre_neighbor() to insure body xcmimage flags are reset
     needed if Verlet::setup::pbc() has remapped/migrated atoms for 2nd run
------------------------------------------------------------------------- */

void FixRigid::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   compute initial fcm and torque on bodies, also initial virial
   reset all particle velocities to be consistent with vcm and omega
------------------------------------------------------------------------- */

void FixRigid::setup(int vflag)
{
  int i,n,ibody;

  // fcm = force on center-of-mass of each rigid body

  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
  }

  // torque = torque on each rigid body

  double **x = atom->x;

  double dx,dy,dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    domain->unmap(x[i],xcmimage[i],unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][0] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][1] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][2] += dx * f[i][1] - dy * f[i][0];
  }

  // extended particles add their torque to torque of body

  if (extended) {
    double **torque_one = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];
      if (eflags[i] & TORQUE) {
        sum[ibody][0] += torque_one[i][0];
        sum[ibody][1] += torque_one[i][1];
        sum[ibody][2] += torque_one[i][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    torque[ibody][0] = all[ibody][0];
    torque[ibody][1] = all[ibody][1];
    torque[ibody][2] = all[ibody][2];
  }

  // zero langextra in case Langevin thermostat not used
  // no point to calling post_force() here since langextra
  // is only added to fcm/torque in final_integrate()

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) langextra[ibody][i] = 0.0;

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set velocities from angmom & omega

  for (ibody = 0; ibody < nbody; ibody++)
    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);

  set_v();

  // guesstimate virial as 2x the set_v contribution

  if (vflag_global)
    for (n = 0; n < 6; n++) virial[n] *= 2.0;
  if (vflag_atom) {
    for (i = 0; i < nlocal; i++)
      for (n = 0; n < 6; n++)
        vatom[i][n] *= 2.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate(int vflag)
{
  double dtfm;

  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // update xcm by full step

    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
    MathExtra::richardson(quat[ibody],angmom[ibody],omega[ibody],
                          inertia[ibody],dtq);
    MathExtra::q_to_exyz(quat[ibody],
                         ex_space[ibody],ey_space[ibody],ez_space[ibody]);
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  set_xv();
}

/* ----------------------------------------------------------------------
   apply Langevin thermostat to all 6 DOF of rigid bodies
   computed by proc 0, broadcast to other procs
   unlike fix langevin, this stores extra force in extra arrays,
     which are added in when final_integrate() calculates a new fcm/torque
------------------------------------------------------------------------- */

void FixRigid::apply_langevin_thermostat()
{
  if (me == 0) {
    double gamma1,gamma2;

    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    t_target = t_start + delta * (t_stop-t_start);
    double tsqrt = sqrt(t_target);

    double boltz = force->boltz;
    double dt = update->dt;
    double mvv2e = force->mvv2e;
    double ftm2v = force->ftm2v;

    for (int i = 0; i < nbody; i++) {
      gamma1 = -masstotal[i] / t_period / ftm2v;
      gamma2 = sqrt(masstotal[i]) * tsqrt *
        sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      langextra[i][0] = gamma1*vcm[i][0] + gamma2*(random->uniform()-0.5);
      langextra[i][1] = gamma1*vcm[i][1] + gamma2*(random->uniform()-0.5);
      langextra[i][2] = gamma1*vcm[i][2] + gamma2*(random->uniform()-0.5);

      gamma1 = -1.0 / t_period / ftm2v;
      gamma2 = tsqrt * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      langextra[i][3] = inertia[i][0]*gamma1*omega[i][0] +
        sqrt(inertia[i][0])*gamma2*(random->uniform()-0.5);
      langextra[i][4] = inertia[i][1]*gamma1*omega[i][1] +
        sqrt(inertia[i][1])*gamma2*(random->uniform()-0.5);
      langextra[i][5] = inertia[i][2]*gamma1*omega[i][2] +
        sqrt(inertia[i][2])*gamma2*(random->uniform()-0.5);
    }
  }

  MPI_Bcast(&langextra[0][0],6*nbody,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   called from FixEnforce2d post_force() for 2d problems
   zero all body values that should be zero for 2d model
------------------------------------------------------------------------- */

void FixRigid::enforce2d()
{
  for (int ibody = 0; ibody < nbody; ibody++) {
    xcm[ibody][2] = 0.0;
    vcm[ibody][2] = 0.0;
    fcm[ibody][2] = 0.0;
    torque[ibody][0] = 0.0;
    torque[ibody][1] = 0.0;
    angmom[ibody][0] = 0.0;
    angmom[ibody][1] = 0.0;
    omega[ibody][0] = 0.0;
    omega[ibody][1] = 0.0;
    if (langflag && langextra) {
      langextra[ibody][2] = 0.0;
      langextra[ibody][3] = 0.0;
      langextra[ibody][4] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::compute_forces_and_torques()
{
  int i,ibody;

  // sum over atoms to get force and torque on rigid body

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];

    domain->unmap(x[i],xcmimage[i],unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][3] += dy*f[i][2] - dz*f[i][1];
    sum[ibody][4] += dz*f[i][0] - dx*f[i][2];
    sum[ibody][5] += dx*f[i][1] - dy*f[i][0];
  }

  // extended particles add their torque to torque of body

  if (extended) {
    double **torque_one = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & TORQUE) {
        sum[ibody][3] += torque_one[i][0];
        sum[ibody][4] += torque_one[i][1];
        sum[ibody][5] += torque_one[i][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  // include Langevin thermostat forces

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0] + langextra[ibody][0];
    fcm[ibody][1] = all[ibody][1] + langextra[ibody][1];
    fcm[ibody][2] = all[ibody][2] + langextra[ibody][2];
    torque[ibody][0] = all[ibody][3] + langextra[ibody][3];
    torque[ibody][1] = all[ibody][4] + langextra[ibody][4];
    torque[ibody][2] = all[ibody][5] + langextra[ibody][5];
  }

  // add gravity force to COM of each body

  if (id_gravity) {
    for (ibody = 0; ibody < nbody; ibody++) {
      fcm[ibody][0] += gvec[0]*masstotal[ibody];
      fcm[ibody][1] += gvec[1]*masstotal[ibody];
      fcm[ibody][2] += gvec[2]*masstotal[ibody];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::post_force(int /*vflag*/)
{
  if (langflag) apply_langevin_thermostat();
  if (earlyflag) compute_forces_and_torques();
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate()
{
  int ibody;
  double dtfm;

  if (!earlyflag) compute_forces_and_torques();

  // update vcm and angmom
  // fflag,tflag = 0 for some dimensions in 2d

  for (ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}

/* ---------------------------------------------------------------------- */

void FixRigid::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtq = 0.5 * step_respa[ilevel];

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixRigid::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ----------------------------------------------------------------------
   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to first-time definition of rigid body in setup_bodies_static()
     or due to box flip
   also adjust imagebody = rigid body image flags, due to xcm remap
   also reset body xcmimage flags of all atoms in bodies
   xcmimage flags are relative to xcm so that body can be unwrapped
   if don't do this, would need xcm to move with true image flags
     then a body could end up very far away from box
     set_xv() will then compute huge displacements every step to
       reset coords of all body atoms to be back inside the box,
       ditto for triclinic box flip, which causes numeric problems
------------------------------------------------------------------------- */

void FixRigid::pre_neighbor()
{
  for (int ibody = 0; ibody < nbody; ibody++)
    domain->remap(xcm[ibody],imagebody[ibody]);
  image_shift();
}

/* ----------------------------------------------------------------------
   reset body xcmimage flags of atoms in bodies
   xcmimage flags are relative to xcm so that body can be unwrapped
   xcmimage = true image flag - imagebody flag
------------------------------------------------------------------------- */

void FixRigid::image_shift()
{
  int ibody;
  imageint tdim,bdim,xdim[3];

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    tdim = image[i] & IMGMASK;
    bdim = imagebody[ibody] & IMGMASK;
    xdim[0] = IMGMAX + tdim - bdim;
    tdim = (image[i] >> IMGBITS) & IMGMASK;
    bdim = (imagebody[ibody] >> IMGBITS) & IMGMASK;
    xdim[1] = IMGMAX + tdim - bdim;
    tdim = image[i] >> IMG2BITS;
    bdim = imagebody[ibody] >> IMG2BITS;
    xdim[2] = IMGMAX + tdim - bdim;

    xcmimage[i] = (xdim[2] << IMG2BITS) | (xdim[1] << IMGBITS) | xdim[0];
  }
}

/* ----------------------------------------------------------------------
   count # of DOF removed by rigid bodies for atoms in igroup
   return total count of DOF
------------------------------------------------------------------------- */

int FixRigid::dof(int tgroup)
{
  // cannot count DOF correctly unless setup_bodies_static() has been called

  if (!setupflag) {
    if (comm->me == 0)
      error->warning(FLERR,"Cannot count rigid body degrees-of-freedom "
                     "before bodies are initialized");
    return 0;
  }

  int tgroupbit = group->bitmask[tgroup];

  // nall = # of point particles in each rigid body
  // mall = # of finite-size particles in each rigid body
  // particles must also be in temperature group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *ncount = new int[nbody];
  int *mcount = new int[nbody];
  for (int ibody = 0; ibody < nbody; ibody++)
    ncount[ibody] = mcount[ibody] = 0;

  for (int i = 0; i < nlocal; i++)
    if (body[i] >= 0 && mask[i] & tgroupbit) {
      // do not count point particles or point dipoles as extended particles
      // a spheroid dipole will be counted as extended
      if (extended && (eflags[i] & ~(POINT | DIPOLE))) mcount[body[i]]++;
      else ncount[body[i]]++;
    }

  int *nall = new int[nbody];
  int *mall = new int[nbody];
  MPI_Allreduce(ncount,nall,nbody,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(mcount,mall,nbody,MPI_INT,MPI_SUM,world);

  // warn if nall+mall != nrigid for any body included in temperature group

  int flag = 0;
  for (int ibody = 0; ibody < nbody; ibody++) {
    if (nall[ibody]+mall[ibody] > 0 &&
        nall[ibody]+mall[ibody] != nrigid[ibody]) flag = 1;
  }
  if (flag && me == 0)
    error->warning(FLERR,"Computing temperature of portions of rigid bodies");

  // remove appropriate DOFs for each rigid body wholly in temperature group
  // N = # of point particles in body
  // M = # of finite-size particles in body
  // 3d body has 3N + 6M dof to start with
  // 2d body has 2N + 3M dof to start with
  // 3d point-particle body with all non-zero I should have 6 dof, remove 3N-6
  // 3d point-particle body (linear) with a 0 I should have 5 dof, remove 3N-5
  // 2d point-particle body should have 3 dof, remove 2N-3
  // 3d body with any finite-size M should have 6 dof, remove (3N+6M) - 6
  // 2d body with any finite-size M should have 3 dof, remove (2N+3M) - 3

  int n = 0;
  nlinear = 0;
  if (domain->dimension == 3) {
    for (int ibody = 0; ibody < nbody; ibody++)
      if (nall[ibody]+mall[ibody] == nrigid[ibody]) {
        n += 3*nall[ibody] + 6*mall[ibody] - 6;
        if (inertia[ibody][0] == 0.0 || inertia[ibody][1] == 0.0 ||
            inertia[ibody][2] == 0.0) {
          n++;
          nlinear++;
        }
      }
  } else if (domain->dimension == 2) {
    for (int ibody = 0; ibody < nbody; ibody++)
      if (nall[ibody]+mall[ibody] == nrigid[ibody])
        n += 2*nall[ibody] + 3*mall[ibody] - 3;
  }

  delete [] ncount;
  delete [] mcount;
  delete [] nall;
  delete [] mall;

  return n;
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

void FixRigid::deform(int flag)
{
  if (flag == 0)
    for (int ibody = 0; ibody < nbody; ibody++)
      domain->x2lamda(xcm[ibody],xcm[ibody]);
  else
    for (int ibody = 0; ibody < nbody; ibody++)
      domain->lamda2x(xcm[ibody],xcm[ibody]);
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_xv()
{
  int ibody;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double ione[3],exone[3],eyone[3],ezone[3],vr[6],p[3][3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set x and v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    // save old positions and velocities for virial

    if (evflag) {
      if (triclinic == 0) {
        x0 = x[i][0] + xbox*xprd;
        x1 = x[i][1] + ybox*yprd;
        x2 = x[i][2] + zbox*zprd;
      } else {
        x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
        x1 = x[i][1] + ybox*yprd + zbox*yz;
        x2 = x[i][2] + zbox*zprd;
      }
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],x[i]);

    v[i][0] = omega[ibody][1]*x[i][2] - omega[ibody][2]*x[i][1] +
      vcm[ibody][0];
    v[i][1] = omega[ibody][2]*x[i][0] - omega[ibody][0]*x[i][2] +
      vcm[ibody][1];
    v[i][2] = omega[ibody][0]*x[i][1] - omega[ibody][1]*x[i][0] +
      vcm[ibody][2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    if (triclinic == 0) {
      x[i][0] += xcm[ibody][0] - xbox*xprd;
      x[i][1] += xcm[ibody][1] - ybox*yprd;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    } else {
      x[i][0] += xcm[ibody][0] - xbox*xprd - ybox*xy - zbox*xz;
      x[i][1] += xcm[ibody][1] - ybox*yprd - zbox*yz;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    }

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set orientation, omega, angmom of each extended particle

  if (extended) {
    double theta_body,theta;
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double **mu = atom->mu;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & SPHERE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::quatquat(quat[ibody],orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,ione,
                                   angmom_one[i]);
      } else if (eflags[i] & LINE) {
        if (quat[ibody][3] >= 0.0) theta_body = 2.0*acos(quat[ibody][0]);
        else theta_body = -2.0*acos(quat[ibody][0]);
        theta = orient[i][0] + theta_body;
        while (theta <= -MY_PI) theta += MY_2PI;
        while (theta > MY_PI) theta -= MY_2PI;
        lbonus[line[i]].theta = theta;
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::quatquat(quat[ibody],orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,
                                   inertiaatom,angmom_one[i]);
      }
      if (eflags[i] & DIPOLE) {
        MathExtra::quat_to_mat(quat[ibody],p);
        MathExtra::matvec(p,dorient[i],mu[i]);
        MathExtra::snormalize3(mu[i][3],mu[i],mu[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigid::set_v()
{
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double ione[3],exone[3],eyone[3],ezone[3],delta[3],vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    const int ibody = body[i];

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],delta);

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[ibody][1]*delta[2] - omega[ibody][2]*delta[1] +
      vcm[ibody][0];
    v[i][1] = omega[ibody][2]*delta[0] - omega[ibody][0]*delta[2] +
      vcm[ibody][1];
    v[i][2] = omega[ibody][0]*delta[1] - omega[ibody][1]*delta[0] +
      vcm[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
      ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

      if (triclinic == 0) {
        x0 = x[i][0] + xbox*xprd;
        x1 = x[i][1] + ybox*yprd;
        x2 = x[i][2] + zbox*zprd;
      } else {
        x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
        x1 = x[i][1] + ybox*yprd + zbox*yz;
        x2 = x[i][2] + zbox*zprd;
      }

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set omega, angmom of each extended particle

  if (extended) {
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      const int ibody = body[i];

      if (eflags[i] & SPHERE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,ione,
                                   angmom_one[i]);
      } else if (eflags[i] & LINE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,
                                   inertiaatom,angmom_one[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   one-time initialization of static rigid body attributes
   sets extended flags, masstotal, center-of-mass
   sets Cartesian and diagonalized inertia tensor
   sets body image flags
   may read some properties from inpfile
------------------------------------------------------------------------- */

void FixRigid::setup_bodies_static()
{
  int i,ibody;

  // extended = 1 if any particle in a rigid body is finite size
  //              or has a dipole moment

  extended = orientflag = dorientflag = 0;

  AtomVecEllipsoid::Bonus *ebonus;
  if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
  AtomVecLine::Bonus *lbonus;
  if (avec_line) lbonus = avec_line->bonus;
  AtomVecTri::Bonus *tbonus;
  if (avec_tri) tbonus = avec_tri->bonus;
  double **mu = atom->mu;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (atom->radius_flag || atom->ellipsoid_flag || atom->line_flag ||
      atom->tri_flag || atom->mu_flag) {
    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      if (radius && radius[i] > 0.0) flag = 1;
      if (ellipsoid && ellipsoid[i] >= 0) flag = 1;
      if (line && line[i] >= 0) flag = 1;
      if (tri && tri[i] >= 0) flag = 1;
      if (mu && mu[i][3] > 0.0) flag = 1;
    }

    MPI_Allreduce(&flag,&extended,1,MPI_INT,MPI_MAX,world);
  }

  // grow extended arrays and set extended flags for each particle
  // orientflag = 4 if any particle stores ellipsoid or tri orientation
  // orientflag = 1 if any particle stores line orientation
  // dorientflag = 1 if any particle stores dipole orientation

  if (extended) {
    if (atom->ellipsoid_flag) orientflag = 4;
    if (atom->line_flag) orientflag = 1;
    if (atom->tri_flag) orientflag = 4;
    if (atom->mu_flag) dorientflag = 1;
    grow_arrays(atom->nmax);

    for (i = 0; i < nlocal; i++) {
      eflags[i] = 0;
      if (body[i] < 0) continue;

      // set to POINT or SPHERE or ELLIPSOID or LINE

      if (radius && radius[i] > 0.0) {
        eflags[i] |= SPHERE;
        eflags[i] |= OMEGA;
        eflags[i] |= TORQUE;
      } else if (ellipsoid && ellipsoid[i] >= 0) {
        eflags[i] |= ELLIPSOID;
        eflags[i] |= ANGMOM;
        eflags[i] |= TORQUE;
      } else if (line && line[i] >= 0) {
        eflags[i] |= LINE;
        eflags[i] |= OMEGA;
        eflags[i] |= TORQUE;
      } else if (tri && tri[i] >= 0) {
        eflags[i] |= TRIANGLE;
        eflags[i] |= ANGMOM;
        eflags[i] |= TORQUE;
      } else eflags[i] |= POINT;

      // set DIPOLE if atom->mu and mu[3] > 0.0

      if (atom->mu_flag && mu[i][3] > 0.0)
        eflags[i] |= DIPOLE;
    }
  }

  // set body xcmimage flags = true image flags

  imageint *image = atom->image;
  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) xcmimage[i] = image[i];
    else xcmimage[i] = 0;

  // compute masstotal & center-of-mass of each rigid body
  // error if image flag is not 0 in a non-periodic dim

  double **x = atom->x;

  int *periodicity = domain->periodicity;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox,ybox,zbox;
  double massone,xunwrap,yunwrap,zunwrap;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    if ((xbox && !periodicity[0]) || (ybox && !periodicity[1]) ||
        (zbox && !periodicity[2]))
      error->one(FLERR,"Fix rigid atom has non-zero image flag "
                 "in a non-periodic dimension");

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    sum[ibody][0] += xunwrap * massone;
    sum[ibody][1] += yunwrap * massone;
    sum[ibody][2] += zunwrap * massone;
    sum[ibody][3] += massone;
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  for (ibody = 0; ibody < nbody; ibody++) {
    masstotal[ibody] = all[ibody][3];
    xcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    xcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    xcm[ibody][2] = all[ibody][2]/masstotal[ibody];
  }

  // set vcm, angmom = 0.0 in case inpfile is used
  // and doesn't overwrite all body's values
  // since setup_bodies_dynamic() will not be called

  for (ibody = 0; ibody < nbody; ibody++) {
    vcm[ibody][0] = vcm[ibody][1] = vcm[ibody][2] = 0.0;
    angmom[ibody][0] = angmom[ibody][1] = angmom[ibody][2] = 0.0;
  }

  // set rigid body image flags to default values

  for (ibody = 0; ibody < nbody; ibody++)
    imagebody[ibody] = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // overwrite masstotal, center-of-mass, image flags with file values
  // inbody[i] = 0/1 if Ith rigid body is initialized by file

  int *inbody;
  if (inpfile) {
    memory->create(inbody,nbody,"rigid:inbody");
    for (ibody = 0; ibody < nbody; ibody++) inbody[ibody] = 0;
    readfile(0,masstotal,xcm,vcm,angmom,imagebody,inbody);
  }

  // remap the xcm of each body back into simulation box
  //   and reset body and atom xcmimage flags via pre_neighbor()

  pre_neighbor();

  // compute 6 moments of inertia of each body in Cartesian reference frame
  // dx,dy,dz = coords relative to center-of-mass
  // symmetric 3x3 inertia tensor stored in Voigt notation as 6-vector

  double dx,dy,dz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += massone * (dy*dy + dz*dz);
    sum[ibody][1] += massone * (dx*dx + dz*dz);
    sum[ibody][2] += massone * (dx*dx + dy*dy);
    sum[ibody][3] -= massone * dy*dz;
    sum[ibody][4] -= massone * dx*dz;
    sum[ibody][5] -= massone * dx*dy;
  }

  // extended particles may contribute extra terms to moments of inertia

  if (extended) {
    double ivec[6];
    double *shape,*quatatom,*inertiaatom;
    double length,theta;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      if (eflags[i] & SPHERE) {
        sum[ibody][0] += SINERTIA*massone * radius[i]*radius[i];
        sum[ibody][1] += SINERTIA*massone * radius[i]*radius[i];
        sum[ibody][2] += SINERTIA*massone * radius[i]*radius[i];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::inertia_ellipsoid(shape,quatatom,massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      } else if (eflags[i] & LINE) {
        length = lbonus[line[i]].length;
        theta = lbonus[line[i]].theta;
        MathExtra::inertia_line(length,theta,massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::inertia_triangle(inertiaatom,quatatom,massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  // overwrite Cartesian inertia tensor with file values

  if (inpfile) readfile(1,NULL,all,NULL,NULL,NULL,inbody);

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy_space = 3 evectors = principal axes of rigid body

  int ierror;
  double cross[3];
  double tensor[3][3],evectors[3][3];

  for (ibody = 0; ibody < nbody; ibody++) {
    tensor[0][0] = all[ibody][0];
    tensor[1][1] = all[ibody][1];
    tensor[2][2] = all[ibody][2];
    tensor[1][2] = tensor[2][1] = all[ibody][3];
    tensor[0][2] = tensor[2][0] = all[ibody][4];
    tensor[0][1] = tensor[1][0] = all[ibody][5];

    ierror = MathExtra::jacobi(tensor,inertia[ibody],evectors);
    if (ierror) error->all(FLERR,
                           "Insufficient Jacobi rotations for rigid body");

    ex_space[ibody][0] = evectors[0][0];
    ex_space[ibody][1] = evectors[1][0];
    ex_space[ibody][2] = evectors[2][0];
    ey_space[ibody][0] = evectors[0][1];
    ey_space[ibody][1] = evectors[1][1];
    ey_space[ibody][2] = evectors[2][1];
    ez_space[ibody][0] = evectors[0][2];
    ez_space[ibody][1] = evectors[1][2];
    ez_space[ibody][2] = evectors[2][2];

    // if any principal moment < scaled EPSILON, set to 0.0

    double max;
    max = MAX(inertia[ibody][0],inertia[ibody][1]);
    max = MAX(max,inertia[ibody][2]);

    if (inertia[ibody][0] < EPSILON*max) inertia[ibody][0] = 0.0;
    if (inertia[ibody][1] < EPSILON*max) inertia[ibody][1] = 0.0;
    if (inertia[ibody][2] < EPSILON*max) inertia[ibody][2] = 0.0;

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed

    MathExtra::cross3(ex_space[ibody],ey_space[ibody],cross);
    if (MathExtra::dot3(cross,ez_space[ibody]) < 0.0)
      MathExtra::negate3(ez_space[ibody]);

    // create initial quaternion

    MathExtra::exyz_to_q(ex_space[ibody],ey_space[ibody],ez_space[ibody],
                         quat[ibody]);
  }

  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  double qc[4],delta[3];
  double *quatatom;
  double theta_body;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) {
      displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
      continue;
    }

    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    delta[0] = xunwrap - xcm[ibody][0];
    delta[1] = yunwrap - xcm[ibody][1];
    delta[2] = zunwrap - xcm[ibody][2];
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],
                                ez_space[ibody],delta,displace[i]);

    if (extended) {
      if (eflags[i] & ELLIPSOID) {
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::qconjugate(quat[ibody],qc);
        MathExtra::quatquat(qc,quatatom,orient[i]);
        MathExtra::qnormalize(orient[i]);
      } else if (eflags[i] & LINE) {
        if (quat[ibody][3] >= 0.0) theta_body = 2.0*acos(quat[ibody][0]);
        else theta_body = -2.0*acos(quat[ibody][0]);
        orient[i][0] = lbonus[line[i]].theta - theta_body;
        while (orient[i][0] <= -MY_PI) orient[i][0] += MY_2PI;
        while (orient[i][0] > MY_PI) orient[i][0] -= MY_2PI;
        if (orientflag == 4) orient[i][1] = orient[i][2] = orient[i][3] = 0.0;
      } else if (eflags[i] & TRIANGLE) {
        quatatom = tbonus[tri[i]].quat;
        MathExtra::qconjugate(quat[ibody],qc);
        MathExtra::quatquat(qc,quatatom,orient[i]);
        MathExtra::qnormalize(orient[i]);
      } else if (orientflag == 4) {
        orient[i][0] = orient[i][1] = orient[i][2] = orient[i][3] = 0.0;
      } else if (orientflag == 1)
        orient[i][0] = 0.0;

      if (eflags[i] & DIPOLE) {
        MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],
                                    ez_space[ibody],mu[i],dorient[i]);
        MathExtra::snormalize3(mu[i][3],dorient[i],dorient[i]);
      } else if (dorientflag)
        dorient[i][0] = dorient[i][1] = dorient[i][2] = 0.0;
    }
  }

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles may contribute extra terms to moments of inertia

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += massone *
      (displace[i][1]*displace[i][1] + displace[i][2]*displace[i][2]);
    sum[ibody][1] += massone *
      (displace[i][0]*displace[i][0] + displace[i][2]*displace[i][2]);
    sum[ibody][2] += massone *
      (displace[i][0]*displace[i][0] + displace[i][1]*displace[i][1]);
    sum[ibody][3] -= massone * displace[i][1]*displace[i][2];
    sum[ibody][4] -= massone * displace[i][0]*displace[i][2];
    sum[ibody][5] -= massone * displace[i][0]*displace[i][1];
  }

  if (extended) {
    double ivec[6];
    double *shape,*inertiaatom;
    double length;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      if (eflags[i] & SPHERE) {
        sum[ibody][0] += SINERTIA*massone * radius[i]*radius[i];
        sum[ibody][1] += SINERTIA*massone * radius[i]*radius[i];
        sum[ibody][2] += SINERTIA*massone * radius[i]*radius[i];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        MathExtra::inertia_ellipsoid(shape,orient[i],massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      } else if (eflags[i] & LINE) {
        length = lbonus[line[i]].length;
        MathExtra::inertia_line(length,orient[i][0],massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        MathExtra::inertia_triangle(inertiaatom,orient[i],massone,ivec);
        sum[ibody][0] += ivec[0];
        sum[ibody][1] += ivec[1];
        sum[ibody][2] += ivec[2];
        sum[ibody][3] += ivec[3];
        sum[ibody][4] += ivec[4];
        sum[ibody][5] += ivec[5];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  // error check that re-computed moments of inertia match diagonalized ones
  // do not do test for bodies with params read from inpfile

  double norm;
  for (ibody = 0; ibody < nbody; ibody++) {
    if (inpfile && inbody[ibody]) continue;
    if (inertia[ibody][0] == 0.0) {
      if (fabs(all[ibody][0]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][0]-inertia[ibody][0])/inertia[ibody][0]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    if (inertia[ibody][1] == 0.0) {
      if (fabs(all[ibody][1]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][1]-inertia[ibody][1])/inertia[ibody][1]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    if (inertia[ibody][2] == 0.0) {
      if (fabs(all[ibody][2]) > TOLERANCE)
        error->all(FLERR,"Fix rigid: Bad principal moments");
    } else {
      if (fabs((all[ibody][2]-inertia[ibody][2])/inertia[ibody][2]) >
          TOLERANCE) error->all(FLERR,"Fix rigid: Bad principal moments");
    }
    norm = (inertia[ibody][0] + inertia[ibody][1] + inertia[ibody][2]) / 3.0;
    if (fabs(all[ibody][3]/norm) > TOLERANCE ||
        fabs(all[ibody][4]/norm) > TOLERANCE ||
        fabs(all[ibody][5]/norm) > TOLERANCE)
      error->all(FLERR,"Fix rigid: Bad principal moments");
  }

  if (inpfile) memory->destroy(inbody);
}

/* ----------------------------------------------------------------------
   one-time initialization of dynamic rigid body attributes
   set vcm and angmom, computed explicitly from constituent particles
   not done if body properites read from file, e.g. for overlapping particles
------------------------------------------------------------------------- */

void FixRigid::setup_bodies_dynamic()
{
  int i,ibody;
  double massone,radone;

  // vcm = velocity of center-of-mass of each rigid body
  // angmom = angular momentum of each rigid body

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double dx,dy,dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    if (rmass) massone = rmass[i];
    else massone = mass[type[i]];

    sum[ibody][0] += v[i][0] * massone;
    sum[ibody][1] += v[i][1] * massone;
    sum[ibody][2] += v[i][2] * massone;

    domain->unmap(x[i],xcmimage[i],unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][3] += dy * massone*v[i][2] - dz * massone*v[i][1];
    sum[ibody][4] += dz * massone*v[i][0] - dx * massone*v[i][2];
    sum[ibody][5] += dx * massone*v[i][1] - dy * massone*v[i][0];
  }

  // extended particles add their rotation to angmom of body

  if (extended) {
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double *radius = atom->radius;
    int *line = atom->line;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & OMEGA) {
        if (eflags[i] & SPHERE) {
          radone = radius[i];
          sum[ibody][3] += SINERTIA*rmass[i] * radone*radone * omega_one[i][0];
          sum[ibody][4] += SINERTIA*rmass[i] * radone*radone * omega_one[i][1];
          sum[ibody][5] += SINERTIA*rmass[i] * radone*radone * omega_one[i][2];
        } else if (eflags[i] & LINE) {
          radone = lbonus[line[i]].length;
          sum[ibody][5] += LINERTIA*rmass[i] * radone*radone * omega_one[i][2];
        }
      }
      if (eflags[i] & ANGMOM) {
        sum[ibody][3] += angmom_one[i][0];
        sum[ibody][4] += angmom_one[i][1];
        sum[ibody][5] += angmom_one[i][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  // normalize velocity of COM

  for (ibody = 0; ibody < nbody; ibody++) {
    vcm[ibody][0] = all[ibody][0]/masstotal[ibody];
    vcm[ibody][1] = all[ibody][1]/masstotal[ibody];
    vcm[ibody][2] = all[ibody][2]/masstotal[ibody];
    angmom[ibody][0] = all[ibody][3];
    angmom[ibody][1] = all[ibody][4];
    angmom[ibody][2] = all[ibody][5];
  }
}

/* ----------------------------------------------------------------------
   read per rigid body info from user-provided file
   which = 0 to read everthing except 6 moments of inertia
   which = 1 to read 6 moments of inertia
   flag inbody = 0 for bodies whose info is read from file
   nlines = # of lines of rigid body info
   one line = rigid-ID mass xcm ycm zcm ixx iyy izz ixy ixz iyz
              vxcm vycm vzcm lx ly lz ix iy iz
------------------------------------------------------------------------- */

void FixRigid::readfile(int which, double *vec,
                        double **array1, double **array2, double **array3,
                        imageint *ivec, int *inbody)
{
  int j,nchunk,id,eofflag,xbox,ybox,zbox;
  int nlines;
  FILE *fp;
  char *eof,*start,*next,*buf;
  char line[MAXLINE];
  char *r_token;

  if (me == 0) {
    fp = fopen(inpfile,"r");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open fix rigid inpfile %s",inpfile);
      error->one(FLERR,str);
    }

    while (1) {
      eof = fgets(line,MAXLINE,fp);
      if (eof == NULL) error->one(FLERR,"Unexpected end of fix rigid file");
      start = &line[strspn(line," \t\n\v\f\r")];
      if (*start != '\0' && *start != '#') break;
    }

    sscanf(line,"%d",&nlines);
  }

  MPI_Bcast(&nlines,1,MPI_INT,0,world);
  if (nlines == 0) error->all(FLERR,"Fix rigid file has no lines");

  char *buffer = new char[CHUNK*MAXLINE];
  char **values = new char*[ATTRIBUTE_PERBODY];

  int nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines-nread,CHUNK);
    eofflag = comm->read_lines_from_file(fp,nchunk,MAXLINE,buffer);
    if (eofflag) error->all(FLERR,"Unexpected end of fix rigid file");

    buf = buffer;
    next = strchr(buf,'\n');
    *next = '\0';
    int nwords = atom->count_words(buf);
    *next = '\n';

    if (nwords != ATTRIBUTE_PERBODY)
      error->all(FLERR,"Incorrect rigid body format in fix rigid file");

    // loop over lines of rigid body attributes
    // tokenize the line into values
    // id = rigid body ID
    // use ID as-is for SINGLE, as mol-ID for MOLECULE, as-is for GROUP
    // for which = 0, store all but inertia in vecs and arrays
    // for which = 1, store inertia tensor array, invert 3,4,5 values to Voigt

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf,'\n');
      r_token = buf;

      values[0] = utils::strtok_r(r_token," \t\n\r\f",&r_token);
      for (j = 1; j < nwords; j++)
        values[j] = utils::strtok_r(NULL," \t\n\r\f",&r_token);

      id = atoi(values[0]);
      if (rstyle == MOLECULE) {
        if (id <= 0 || id > maxmol)
          error->all(FLERR,"Invalid rigid body ID in fix rigid file");
        id = mol2body[id];
      } else id--;

      if (id < 0 || id >= nbody)
        error->all(FLERR,"Invalid rigid body ID in fix rigid file");
      inbody[id] = 1;

      if (which == 0) {
        vec[id] = atof(values[1]);
        array1[id][0] = atof(values[2]);
        array1[id][1] = atof(values[3]);
        array1[id][2] = atof(values[4]);
        array2[id][0] = atof(values[11]);
        array2[id][1] = atof(values[12]);
        array2[id][2] = atof(values[13]);
        array3[id][0] = atof(values[14]);
        array3[id][1] = atof(values[15]);
        array3[id][2] = atof(values[16]);
        xbox = atoi(values[17]);
        ybox = atoi(values[18]);
        zbox = atoi(values[19]);
        ivec[id] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
          (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
          (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
      } else {
        array1[id][0] = atof(values[5]);
        array1[id][1] = atof(values[6]);
        array1[id][2] = atof(values[7]);
        array1[id][3] = atof(values[10]);
        array1[id][4] = atof(values[9]);
        array1[id][5] = atof(values[8]);
      }

      buf = next + 1;
    }

    nread += nchunk;
  }

  if (me == 0) fclose(fp);

  delete [] buffer;
  delete [] values;
}

/* ----------------------------------------------------------------------
   write out restart info for mass, COM, inertia tensor, image flags to file
   identical format to inpfile option, so info can be read in when restarting
   only proc 0 writes list of global bodies to file
------------------------------------------------------------------------- */

void FixRigid::write_restart_file(char *file)
{
  if (me) return;

  char outfile[128];
  snprintf(outfile,128,"%s.rigid",file);
  FILE *fp = fopen(outfile,"w");
  if (fp == NULL) {
    char str[192];
    snprintf(str,192,"Cannot open fix rigid restart file %s",outfile);
    error->one(FLERR,str);
  }

  fprintf(fp,"# fix rigid mass, COM, inertia tensor info for "
          "%d bodies on timestep " BIGINT_FORMAT "\n\n",
          nbody,update->ntimestep);
  fprintf(fp,"%d\n",nbody);

  // compute I tensor against xyz axes from diagonalized I and current quat
  // Ispace = P Idiag P_transpose
  // P is stored column-wise in exyz_space

  int xbox,ybox,zbox;
  double p[3][3],pdiag[3][3],ispace[3][3];

  int id;
  for (int i = 0; i < nbody; i++) {
    if (rstyle == SINGLE || rstyle == GROUP) id = i;
    else id = body2mol[i];

    MathExtra::col2mat(ex_space[i],ey_space[i],ez_space[i],p);
    MathExtra::times3_diag(p,inertia[i],pdiag);
    MathExtra::times3_transpose(pdiag,p,ispace);

    xbox = (imagebody[i] & IMGMASK) - IMGMAX;
    ybox = (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (imagebody[i] >> IMG2BITS) - IMGMAX;

    fprintf(fp,"%d %-1.16e %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e "
            "%d %d %d\n",
            id,masstotal[i],xcm[i][0],xcm[i][1],xcm[i][2],
            ispace[0][0],ispace[1][1],ispace[2][2],
            ispace[0][1],ispace[0][2],ispace[1][2],
            vcm[i][0],vcm[i][1],vcm[i][2],
            angmom[i][0],angmom[i][1],angmom[i][2],
            xbox,ybox,zbox);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixRigid::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(imageint);
  bytes += nmax*3 * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);    // vatom
  if (extended) {
    bytes += nmax * sizeof(int);
    if (orientflag) bytes = nmax*orientflag * sizeof(double);
    if (dorientflag) bytes = nmax*3 * sizeof(double);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixRigid::grow_arrays(int nmax)
{
  memory->grow(body,nmax,"rigid:body");
  memory->grow(xcmimage,nmax,"rigid:xcmimage");
  memory->grow(displace,nmax,3,"rigid:displace");
  if (extended) {
    memory->grow(eflags,nmax,"rigid:eflags");
    if (orientflag) memory->grow(orient,nmax,orientflag,"rigid:orient");
    if (dorientflag) memory->grow(dorient,nmax,3,"rigid:dorient");
  }

  // check for regrow of vatom
  // must be done whether per-atom virial is accumulated on this step or not
  //   b/c this is only time grow_array() may be called
  // need to regrow b/c vatom is calculated before and after atom migration

  if (nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->grow(vatom,maxvatom,6,"fix:vatom");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixRigid::copy_arrays(int i, int j, int /*delflag*/)
{
  body[j] = body[i];
  xcmimage[j] = xcmimage[i];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];
  if (extended) {
    eflags[j] = eflags[i];
    for (int k = 0; k < orientflag; k++)
      orient[j][k] = orient[i][k];
    if (dorientflag) {
      dorient[j][0] = dorient[i][0];
      dorient[j][1] = dorient[i][1];
      dorient[j][2] = dorient[i][2];
    }
  }

  // must also copy vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++)
      vatom[j][k] = vatom[i][k];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixRigid::set_arrays(int i)
{
  body[i] = -1;
  xcmimage[i] = 0;
  displace[i][0] = 0.0;
  displace[i][1] = 0.0;
  displace[i][2] = 0.0;

  // must also zero vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++)
      vatom[i][k] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixRigid::pack_exchange(int i, double *buf)
{
  buf[0] = ubuf(body[i]).d;
  buf[1] = ubuf(xcmimage[i]).d;
  buf[2] = displace[i][0];
  buf[3] = displace[i][1];
  buf[4] = displace[i][2];
  if (!extended) return 5;

  int m = 5;
  buf[m++] = eflags[i];
  for (int j = 0; j < orientflag; j++)
    buf[m++] = orient[i][j];
  if (dorientflag) {
    buf[m++] = dorient[i][0];
    buf[m++] = dorient[i][1];
    buf[m++] = dorient[i][2];
  }

  // must also pack vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++)
      buf[m++] = vatom[i][k];

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixRigid::unpack_exchange(int nlocal, double *buf)
{
  body[nlocal] = (int) ubuf(buf[0]).i;
  xcmimage[nlocal] = (imageint) ubuf(buf[1]).i;
  displace[nlocal][0] = buf[2];
  displace[nlocal][1] = buf[3];
  displace[nlocal][2] = buf[4];
  if (!extended) return 5;

  int m = 5;
  eflags[nlocal] = static_cast<int> (buf[m++]);
  for (int j = 0; j < orientflag; j++)
    orient[nlocal][j] = buf[m++];
  if (dorientflag) {
    dorient[nlocal][0] = buf[m++];
    dorient[nlocal][1] = buf[m++];
    dorient[nlocal][2] = buf[m++];
  }

  // must also unpack vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++)
      vatom[nlocal][k] = buf[m++];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixRigid::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
}

/* ----------------------------------------------------------------------
   zero linear momentum of each rigid body
   set Vcm to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixRigid::zero_momentum()
{
  for (int ibody = 0; ibody < nbody; ibody++)
    vcm[ibody][0] = vcm[ibody][1] = vcm[ibody][2] = 0.0;

  evflag = 0;
  set_v();
}

/* ----------------------------------------------------------------------
   zero angular momentum of each rigid body
   set angmom/omega to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixRigid::zero_rotation()
{
  for (int ibody = 0; ibody < nbody; ibody++) {
    angmom[ibody][0] = angmom[ibody][1] = angmom[ibody][2] = 0.0;
    omega[ibody][0] = omega[ibody][1] = omega[ibody][2] = 0.0;
  }

  evflag = 0;
  set_v();
}

/* ---------------------------------------------------------------------- */

int FixRigid::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"bodyforces") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (strcmp(arg[1],"early") == 0) earlyflag = 1;
    else if (strcmp(arg[1],"late") == 0) earlyflag = 0;
    else error->all(FLERR,"Illegal fix_modify command");

    // reset fix mask
    // must do here and not in init,
    // since modify.cpp::init() uses fix masks before calling fix::init()

    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->id,id) == 0) {
        if (earlyflag) modify->fmask[i] |= POST_FORCE;
        else if (!langflag) modify->fmask[i] &= ~POST_FORCE;
        break;
      }
    return 2;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
   non-active DOF are removed by fflag/tflag and in tfactor
------------------------------------------------------------------------- */

double FixRigid::compute_scalar()
{
  double wbody[3],rot[3][3];

  double t = 0.0;
  for (int i = 0; i < nbody; i++) {
    t += masstotal[i] * (fflag[i][0]*vcm[i][0]*vcm[i][0] +
                         fflag[i][1]*vcm[i][1]*vcm[i][1] +
                         fflag[i][2]*vcm[i][2]*vcm[i][2]);

    // wbody = angular velocity in body frame

    MathExtra::quat_to_mat(quat[i],rot);
    MathExtra::transpose_matvec(rot,angmom[i],wbody);
    if (inertia[i][0] == 0.0) wbody[0] = 0.0;
    else wbody[0] /= inertia[i][0];
    if (inertia[i][1] == 0.0) wbody[1] = 0.0;
    else wbody[1] /= inertia[i][1];
    if (inertia[i][2] == 0.0) wbody[2] = 0.0;
    else wbody[2] /= inertia[i][2];

    t += tflag[i][0]*inertia[i][0]*wbody[0]*wbody[0] +
      tflag[i][1]*inertia[i][1]*wbody[1]*wbody[1] +
      tflag[i][2]*inertia[i][2]*wbody[2]*wbody[2];
  }

  t *= tfactor;
  return t;
}

/* ---------------------------------------------------------------------- */

void *FixRigid::extract(const char *str, int &dim)
{
  if (strcmp(str,"body") == 0) {
    dim = 1;
    return body;
  }
  if (strcmp(str,"masstotal") == 0) {
    dim = 1;
    return masstotal;
  }
  if (strcmp(str,"t_target") == 0) {
    dim = 0;
    return &t_target;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
------------------------------------------------------------------------- */

double FixRigid::extract_ke()
{
  double ke = 0.0;
  for (int i = 0; i < nbody; i++)
    ke += masstotal[i] *
      (vcm[i][0]*vcm[i][0] + vcm[i][1]*vcm[i][1] + vcm[i][2]*vcm[i][2]);

  return 0.5*ke;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

double FixRigid::extract_erotational()
{
  double wbody[3],rot[3][3];

  double erotate = 0.0;
  for (int i = 0; i < nbody; i++) {

    // wbody = angular velocity in body frame

    MathExtra::quat_to_mat(quat[i],rot);
    MathExtra::transpose_matvec(rot,angmom[i],wbody);
    if (inertia[i][0] == 0.0) wbody[0] = 0.0;
    else wbody[0] /= inertia[i][0];
    if (inertia[i][1] == 0.0) wbody[1] = 0.0;
    else wbody[1] /= inertia[i][1];
    if (inertia[i][2] == 0.0) wbody[2] = 0.0;
    else wbody[2] /= inertia[i][2];

    erotate += inertia[i][0]*wbody[0]*wbody[0] +
      inertia[i][1]*wbody[1]*wbody[1] + inertia[i][2]*wbody[2]*wbody[2];
  }

  return 0.5*erotate;
}

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   15 values per body
   xcm = 0,1,2; vcm = 3,4,5; fcm = 6,7,8; torque = 9,10,11; image = 12,13,14
------------------------------------------------------------------------- */

double FixRigid::compute_array(int i, int j)
{
  if (j < 3) return xcm[i][j];
  if (j < 6) return vcm[i][j-3];
  if (j < 9) return fcm[i][j-6];
  if (j < 12) return torque[i][j-9];
  if (j == 12) return (imagebody[i] & IMGMASK) - IMGMAX;
  if (j == 13) return (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
  return (imagebody[i] >> IMG2BITS) - IMGMAX;
}
