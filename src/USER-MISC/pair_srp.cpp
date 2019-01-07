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
   Contributing authors: Timothy Sirk (ARL), Pieter in't Veld (BASF)

This pair style srp command calculates a segmental repulsive force
between bonds. This is useful for preventing the crossing of bonds if
soft non-bonded potentials are used, such as DPD polymer chains.

See the doc page for pair_style srp command for usage instructions.

There is an example script for this package in examples/USER/srp.

Please contact Timothy Sirk for questions (tim.sirk@us.army.mil).
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "pair_srp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_srp.h"
#include "thermo.h"
#include "output.h"
#include "citeme.h"

using namespace LAMMPS_NS;

#define SMALL 1.0e-10
#define BIG 1e10
#define ONETWOBIT 0x40000000

static const char cite_srp[] =
  "@Article{Sirk2012\n"
  " author = {T. Sirk and Y. Sliozberg and J. Brennan and M. Lisal and J. Andzelm},\n"
  " title = {An enhanced entangled polymer model for dissipative particle dynamics},\n"
  " journal = {J.~Chem.~Phys.},\n"
  " year =    2012,\n"
  " volume =  136,\n"
  " pages =   {134903}\n"
  "}\n\n";

static int srp_instance = 0;

/* ----------------------------------------------------------------------
 set size of pair comms in constructor
 ---------------------------------------------------------------------- */

PairSRP::PairSRP(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;

  if (lmp->citeme) lmp->citeme->add(cite_srp);

  nextra = 1;
  segment = NULL;

  // generate unique fix-id for this pair style instance
  fix_id = strdup("XX_FIX_SRP");
  fix_id[0] = '0' + srp_instance / 10;
  fix_id[1] = '0' + srp_instance % 10;
  ++srp_instance;

  // create fix SRP instance here, as it has to
  // be executed before all other fixes
  char **fixarg = new char*[3];
  fixarg[0] = fix_id;
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "SRP";
  modify->add_fix(3,fixarg);
  f_srp = (FixSRP *) modify->fix[modify->nfix-1];
  delete [] fixarg;
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSRP::allocate()
{
    allocated = 1;
    // particles of bptype inserted by fix srp
    // bptype is the highest numbered atom type
    int n = bptype;
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(cut, n + 1, n + 1, "pair:cut");
    memory->create(a0, n + 1, n + 1, "pair:a0");

    // setflag for atom types
    memory->create(setflag,n+1,n+1,"pair:setflag");
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            setflag[i][j] = 0;

    maxcount = 0;
}

/* ----------------------------------------------------------------------
 free
 ------------------------------------------------------------------------- */

PairSRP::~PairSRP()
{
    if (allocated)
    {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(cut);
        memory->destroy(a0);
        memory->destroy(segment);
    }

  // check nfix in case all fixes have already been deleted
  if (modify->nfix) modify->delete_fix(fix_id);
  free(fix_id);
}

/* ----------------------------------------------------------------------
 compute bond-bond repulsions
 ------------------------------------------------------------------------- */

void PairSRP::compute(int eflag, int vflag)

{
    // setup energy and virial
    if (eflag || vflag)
        ev_setup(eflag, vflag);
    else
        evflag = vflag_fdotr = 0;

    double **x = atom->x;
    double **f = atom->f;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    int i0, i1, j0, j1;
    int i,j,ii,jj,inum,jnum;
    double dijsq, dij;

    int *ilist,*jlist,*numneigh,**firstneigh;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    double dx,dy,dz,ti,tj;
    double wd, lever0, lever1, evdwl, fpair;
    double fxlever0, fylever0, fzlever0, fxlever1, fylever1, fzlever1;
    double fx, fy, fz;
    evdwl = 0.0;

    // mapping global to local for atoms inside bond particles
    // exclude 1-2 neighs if requested
    if (neighbor->ago == 0){
      remapBonds(nall);
      if(exclude) onetwoexclude(ilist, inum, jlist, numneigh, firstneigh);
    }

  // this pair style only used with hybrid
  // due to exclusions
  // each atom i is type bptype
  // each neigh j is type bptype

  // using midpoint distance option
  if(midpoint){

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];
      jnum = numneigh[i];
      // two atoms inside bond particle
      i0 = segment[i][0];
      j0 = segment[i][1];

      for (jj = 0; jj < jnum; jj++) {

        jlist = firstneigh[i];
        j = jlist[jj];

        // enforce 1-2 exclusions
        if ((sbmask(j) & exclude))
          continue;

        j &= NEIGHMASK;
        //retrieve atoms from bond particle
        i1 = segment[j][0];
        j1 = segment[j][1];

        // midpt dist bond 0 and 1
        dx = 0.5*(x[i0][0] - x[i1][0] + x[j0][0] - x[j1][0]);
        dy = 0.5*(x[i0][1] - x[i1][1] + x[j0][1] - x[j1][1]);
        dz = 0.5*(x[i0][2] - x[i1][2] + x[j0][2] - x[j1][2]);
        dijsq = dx*dx + dy*dy + dz*dz;

        if (dijsq < cutsq[bptype][bptype]){
        dij = sqrt(dijsq);

        if (dij < SMALL)
          continue;     // dij can be 0.0 with soft potentials

        wd = 1.0 - dij / cut[bptype][bptype];
        fpair = 0.5 * a0[bptype][bptype] * wd / dij; // 0.5 factor for lever rule

        // force for bond 0, beads 0,1
        //force between bonds
        fx = fpair * dx;
        fy = fpair * dy;
        fz = fpair * dz;

        f[i0][0] += fx; //keep force sign for bond 0
        f[i0][1] += fy;
        f[i0][2] += fz;

        f[j0][0] += fx;
        f[j0][1] += fy;
        f[j0][2] += fz;

        f[i1][0] -= fx; //flip force sign for bond 1
        f[i1][1] -= fy;
        f[i1][2] -= fz;

        f[j1][0] -= fx;
        f[j1][1] -= fy;
        f[j1][2] -= fz;

        // ************************************************* //

        if (eflag){
          evdwl = 0.5 * a0[bptype][bptype] * cut[bptype][bptype] * wd * wd;
        }

        if (evflag){
          ev_tally(i0,i1,nlocal,1,0.5*evdwl,0.0,fpair,dx,dy,dz);
          ev_tally(j0,j1,nlocal,1,0.5*evdwl,0.0,fpair,dx,dy,dz);
        }

        if (vflag_fdotr) virial_fdotr_compute();

        }
      }
   }
 } else {
  // using min distance option

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];
      jnum = numneigh[i];
      i0 = segment[i][0];
      j0 = segment[i][1];

      for (jj = 0; jj < jnum; jj++) {

        jlist = firstneigh[i];
        j = jlist[jj];

        // enforce 1-2 exclusions
        if ((sbmask(j) & exclude))
          continue;

        j &= NEIGHMASK;

        i1 = segment[j][0];
        j1 = segment[j][1];

        getMinDist(x, dx, dy, dz, ti, tj, i0, j0, i1, j1);
        dijsq = dx*dx + dy*dy + dz*dz;

        if (dijsq < cutsq[bptype][bptype]){

        dij = sqrt(dijsq);

        if (dij < SMALL)
      continue;     // dij can be 0.0 with soft potentials

        wd = 1.0 - dij / cut[bptype][bptype];
        fpair = a0[bptype][bptype] * wd / dij;

        // force for bond 0, beads 0,1
        lever0 = 0.5 + ti; // assign force according to lever rule
        lever1 = 0.5 + tj; // assign force according to lever rule
        //force between bonds
        fx = fpair * dx;
        fy = fpair * dy;
        fz = fpair * dz;

        //decompose onto atoms
        fxlever0 = fx * lever0;
        fylever0 = fy * lever0;
        fzlever0 = fz * lever0;
        fxlever1 = fx * lever1;
        fylever1 = fy * lever1;
        fzlever1 = fz * lever1;

        f[i0][0] += fxlever0; //keep force sign for bond 0
        f[i0][1] += fylever0;
        f[i0][2] += fzlever0;

        f[j0][0] += (fx - fxlever0);
        f[j0][1] += (fy - fylever0);
        f[j0][2] += (fz - fzlever0);

        f[i1][0] -= fxlever1; //flip force sign for bond 1
        f[i1][1] -= fylever1;
        f[i1][2] -= fzlever1;

        f[j1][0] -= (fx - fxlever1);
        f[j1][1] -= (fy - fylever1);
        f[j1][2] -= (fz - fzlever1);

        // ************************************************* //

        if (eflag){
          evdwl = 0.5 * a0[bptype][bptype] * cut[bptype][bptype] * wd * wd;
        }

        if (evflag){
          ev_tally(i0,i1,nlocal,1,0.5*evdwl,0.0,0.5*fpair,dx,dy,dz);
          ev_tally(j0,j1,nlocal,1,0.5*evdwl,0.0,0.5*fpair,dx,dy,dz);
        }

       if (vflag_fdotr) virial_fdotr_compute();

      }
    }
  }
 }
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSRP::settings(int narg, char **arg)
{
  if (narg < 3 || narg > 7)
    error->all(FLERR,"Illegal pair_style command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair_style srp requires atom IDs");

  cut_global = force->numeric(FLERR,arg[0]);
  // wildcard
  if (strcmp(arg[1],"*") == 0) {
    btype = 0;
  } else {
    btype = force->inumeric(FLERR,arg[1]);
    if ((btype > atom->nbondtypes) || (btype <= 0))
      error->all(FLERR,"Illegal pair_style command");
  }

  // settings
  midpoint = 0;
  min = 0;

  if (strcmp(arg[2],"min") == 0) min = 1;
  else if (strcmp(arg[2],"mid") == 0) midpoint = 1;
  else
    error->all(FLERR,"Illegal pair_style command");

  int iarg = 3;
  // default exclude 1-2
  // scaling for 1-2, etc not supported
  exclude = 1;

  // use last atom type by default for bond particles
  bptype = atom->ntypes;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair srp command");
      if (strcmp(arg[iarg+1],"yes") == 0)
        exclude = 1;
      if (strcmp(arg[iarg+1],"no") == 0){
        if (min) error->all(FLERR,"Illegal exclude option in pair srp command");
        exclude = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"bptype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair srp command");
      bptype = force->inumeric(FLERR,arg[iarg+1]);
      if ((bptype < 1) || (bptype > atom->ntypes))
        error->all(FLERR,"Illegal bond particle type for srp");
      iarg += 2;
    } else error->all(FLERR,"Illegal pair srp command");
  }

  // reset cutoffs if explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= bptype; i++)
      for (j = i; j <= bptype; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
 set coeffs
 ------------------------------------------------------------------------- */

void PairSRP::coeff(int narg, char **arg)
{
    if (narg < 3 || narg > 4)
        error->all(FLERR,"PairSRP: Incorrect args for pair coeff");
    if (!allocated) allocate();

    // set ij bond-bond cutoffs
    int ilo, ihi, jlo, jhi;
    force->bounds(FLERR,arg[0], bptype, ilo, ihi);
    force->bounds(FLERR,arg[1], bptype, jlo, jhi);

    double a0_one = force->numeric(FLERR,arg[2]);
    double cut_one = cut_global;
    if (narg == 4)  cut_one = force->numeric(FLERR,arg[3]);

    int count = 0;
    for (int i = ilo; i <= ihi; i++)
    {
        for (int j = MAX(jlo,i); j <= jhi; j++)
        {
            a0[i][j] = a0_one;
            cut[i][j] = cut_one;
            cutsq[i][j] = cut_one * cut_one;
            setflag[i][j] = 1;
            count++;
        }
    }

    if (count == 0) error->warning(FLERR,"PairSRP: No pair coefficients were set");
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSRP::init_style()
{
  if (!force->newton_pair)
    error->all(FLERR,"PairSRP: Pair srp requires newton pair on");

  // verify that fix SRP is still defined and has not been changed.
  int ifix = modify->find_fix(fix_id);
  if (f_srp != (FixSRP *)modify->fix[ifix])
    error->all(FLERR,"Fix SRP has been changed unexpectedly");

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Using type %d for bond particles\n",bptype);
    if (logfile) fprintf(logfile,"Using type %d for bond particles\n",bptype);
  }

  // set bond and bond particle types in fix srp
  // bonds of this type will be represented by bond particles
  // if bond type is 0, then all bonds have bond particles
  // btype = bond type
  char c0[20];
  char* arg0[2];
  sprintf(c0, "%d", btype);
  arg0[0] = (char *) "btype";
  arg0[1] = c0;
  f_srp->modify_params(2, arg0);

  // bptype = bond particle type
  sprintf(c0, "%d", bptype);
  arg0[0] = (char *) "bptype";
  arg0[1] = c0;
  f_srp->modify_params(2, arg0);

  // bond particles do not contribute to energy or virial
  // bond particles do not belong to group all
  // but thermo normalization is by nall
  // therefore should turn off normalization
  int me;
  MPI_Comm_rank(world,&me);
  char *arg1[2];
  arg1[0] = (char *) "norm";
  arg1[1] = (char *) "no";
  output->thermo->modify_params(2, arg1);
  if (me == 0)
    error->message(FLERR,"Thermo normalization turned off by pair srp");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSRP::init_one(int i, int j)
{

 if (setflag[i][j] == 0) error->all(FLERR,"PairSRP: All pair coeffs are not set");

  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
 find min distance for bonds i0/j0 and i1/j1
 ------------------------------------------------------------------------- */
inline void PairSRP::getMinDist(double** &x, double &dx, double &dy, double &dz, double &ti, double &tj, int &i0, int &j0, int &i1, int &j1)
{
    // move these outside the loop
    double diffx0, diffy0, diffz0, diffx1, diffy1, diffz1, dPx, dPy, dPz, RiRi, RiRj, RjRj;
    double denom, termx0, termy0, termz0, num0, termx1, termy1, termz1, num1;

    // compute midpt dist from 1st atom, 1st bond
    diffx0 = x[j0][0] - x[i0][0]; // x,y,z from bond 0
    diffy0 = x[j0][1] - x[i0][1];
    diffz0 = x[j0][2] - x[i0][2];

    // compute midpt dist from 1st atom, 2nd bond
    diffx1 = x[j1][0] - x[i1][0];
    diffy1 = x[j1][1] - x[i1][1];
    diffz1 = x[j1][2] - x[i1][2];

    // midpoint distance
    dPx = 0.5*(diffx0-diffx1) + x[i0][0]-x[i1][0];
    dPy = 0.5*(diffy0-diffy1) + x[i0][1]-x[i1][1];
    dPz = 0.5*(diffz0-diffz1) + x[i0][2]-x[i1][2];

    // Ri^2 Rj^2
    RiRi = diffx0*diffx0 + diffy0*diffy0 + diffz0*diffz0;
    RiRj = diffx0*diffx1 + diffy0*diffy1 + diffz0*diffz1;
    RjRj = diffx1*diffx1 + diffy1*diffy1 + diffz1*diffz1;
    denom = RiRj*RiRj - RiRi*RjRj;

    // handle case of parallel lines
    // reduce to midpt distance
    if (fabs(denom) < SMALL){
        if(denom < 0) denom = -BIG;
        else denom = BIG;
    }

    // calc ti
    termx0 = RiRj*diffx1 - RjRj*diffx0;
    termy0 = RiRj*diffy1 - RjRj*diffy0;
    termz0 = RiRj*diffz1 - RjRj*diffz0;
    num0 = dPx*termx0 + dPy*termy0 + dPz*termz0;
    ti = num0 / denom;
    if (ti > 0.5) ti = 0.5;
    if (ti < -0.5) ti = -0.5;

    // calc tj
    termx1 = RiRj*diffx0 - RiRi*diffx1;
    termy1 = RiRj*diffy0 - RiRi*diffy1;
    termz1 = RiRj*diffz0 - RiRi*diffz1;
    num1 = dPx*termx1 + dPy*termy1 + dPz*termz1;
    tj = -num1/ denom;
    if (tj > 0.5)  tj = 0.5;
    if (tj < -0.5) tj = -0.5;

    // min dist
    dx = dPx - ti*diffx0 + tj*diffx1;
    dy = dPy - ti*diffy0 + tj*diffy1;
    dz = dPz - ti*diffz0 + tj*diffz1;
}

/* --------------------------------------------------------
map global id of atoms in stored by each bond particle
 ------------------------------------------------------- */
inline void PairSRP::remapBonds(int &nall)
{
  if(nall > maxcount){
    memory->grow(segment, nall, 2, "pair:segment");
    maxcount = nall;
  }

  // loop over all bond particles
  // each bond paricle holds two bond atoms
  // map global ids of bond atoms to local ids
  // might not be able to map both bond atoms of j, if j is outside neighcut
  // these are not on neighlist, so are not used
  int tmp;
  srp = f_srp->array_atom;

    for (int i = 0; i < nall; i++) {
      if(atom->type[i] == bptype){
        // tmp is local id
        // tmp == -1 is ok
        tmp = atom->map((int)srp[i][0]);
        segment[i][0] = domain->closest_image(i,tmp);
        // repeat with other id
        tmp = atom->map((int)srp[i][1]);
        segment[i][1] = domain->closest_image(i,tmp);
      }
    }
}

/* --------------------------------------------------------
add exclusions for 1-2 neighs, if requested
more complex exclusions or scaling probably not needed
 ------------------------------------------------------- */
inline void PairSRP::onetwoexclude(int* &ilist, int &inum, int* &jlist, int* &numneigh, int** &firstneigh)
{
    int i0, i1, j0, j1;
    int i,j,ii,jj,jnum;

    // encode neighs with exclusions
    // only need 1-2 info for normal uses of srp
    // add 1-3, etc later if ever needed

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];
      jnum = numneigh[i];
      // two atoms inside bond particle
      i0 = segment[i][0];
      j0 = segment[i][1];

      for (jj = 0; jj < jnum; jj++) {

        jlist = firstneigh[i];
        j = jlist[jj];
        j &= NEIGHMASK;
        //two atoms inside bond particle
        i1 = segment[j][0];
        j1 = segment[j][1];

        // check for a 1-2 neigh
        if(i0 == i1 || i0 == j1 || i1 == j0 || j0 == j1){
          j |= ONETWOBIT;
          jlist[jj] = j;
        }
      }
    }
}

/* ----------------------------------------------------------------------
proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSRP::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,a0[i][i]);
}

/* ----------------------------------------------------------------------
proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSRP::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g\n",i,j,a0[i][j],cut[i][j]);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSRP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSRP::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          printf(" i %d j %d \n",i,j);
          fread(&a0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSRP::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&bptype,sizeof(int),1,fp);
  fwrite(&btype,sizeof(int),1,fp);
  fwrite(&min,sizeof(int),1,fp);
  fwrite(&midpoint,sizeof(int),1,fp);
  fwrite(&exclude,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSRP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&bptype,sizeof(int),1,fp);
    fread(&btype,sizeof(int),1,fp);
    fread(&min,sizeof(int),1,fp);
    fread(&midpoint,sizeof(int),1,fp);
    fread(&exclude,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
}
