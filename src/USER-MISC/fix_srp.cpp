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
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_srp.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "atom_vec.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSRP::FixSRP(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // settings
  nevery=1;
  peratom_freq = 1;
  time_integrate = 0;
  create_attribute = 0;
  comm_border = 2;

  // restart settings
  restart_global = 1;
  restart_peratom = 1;
  restart_pbc = 1;

  // per-atom array width 2
  peratom_flag = 1;
  size_peratom_cols = 2;                                                  

  // initial allocation of atom-based array                      
  // register with Atom class
  array = NULL;
  grow_arrays(atom->nmax);                                                

  // extends pack_exchange() 
  atom->add_callback(0);
  atom->add_callback(1); // restart
  atom->add_callback(2); 

  // zero 
  for (int i = 0; i < atom->nmax; i++) 
    for (int m = 0; m < 3; m++) 
      array[i][m] = 0.0; 

  // assume setup of fix is needed to insert particles
  // yes if reading datafile
  // maybe if reading restart
  // a datafile cannot contain BPs
  // a restart file written during the run has BPs as per atom data
  // a restart file written after the run does not have BPs
  setup = 1; 
}

/* ---------------------------------------------------------------------- */

FixSRP::~FixSRP()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1); // for restart
  atom->delete_callback(id,2);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixSRP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_EXCHANGE;
  mask |= POST_RUN;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSRP::init()
{
  if (force->pair_match("hybrid",1) == NULL)
    error->all(FLERR,"Cannot use pair srp without pair_style hybrid");

  // reserve the highest numbered atom type for bond particles (BPs)
  // set bptype if not read from restart 
  if(!restart_reset) bptype = atom->ntypes;

  // check for prescense of extra atom type for bond particle
  // if reading a rst file written in the middle of run
  // then bond particles are already present
  // otherwise need to insert particles in setup
   for(int i=0; i < atom->nlocal; i++)
      if(atom->type[i] == bptype){
        // error if missing extra atom type in either
        // data file or rst file without this fix
        if(!restart_reset)
          error->all(FLERR,"Fix SRP requires an extra atom type");
        else
          setup = 0; 
      }
 
  // setup neigh exclusions for diff atom types
  // BPs do not interact with other types
  // type bptype only interacts with itself
  char* arg1[4];
  arg1[0] = "exclude"; 
  arg1[1] = "type"; 
  char c0[20];
  char c1[20];

  for(int z = 1; z < bptype; z++)
  {
   sprintf(c0, "%d", z);
   arg1[2] = c0; 
   sprintf(c1, "%d", bptype);
   arg1[3] = c1; 
   neighbor->modify_params(4, arg1);
  }
}

/* ----------------------------------------------------------------------
   insert bond particles
------------------------------------------------------------------------- */

void FixSRP::setup_pre_force(int zz)
{
  if(!setup) return;

  double rsqold = 0.0;
  double delx, dely, delz, rmax, rsq, rsqmax;
  double **x = atom->x;
  bigint nall = atom->nlocal + atom->nghost;
  double xold[nall][3];

  // make a copy of coordinates and tags
  // create_atom overwrites ghost atoms
  for(int i = 0; i < nall; i++){
    xold[i][0] = x[i][0];
    xold[i][1] = x[i][1];
    xold[i][2] = x[i][2];
  }

  int *tag = atom->tag;
  int tagold[nall];

  for(int i = 0; i < nall; i++){
    tagold[i]=tag[i];
  }

  int nlocal = atom->nlocal;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nadd = 0;
  int i,j;

  for (int n = 0; n < nbondlist; n++) {

   // consider only the user defined bond type
   if(bondlist[n][2] != btype) continue;

   i = bondlist[n][0];
   j = bondlist[n][1];

   // position of bond i
   xone[0] = (xold[i][0] + xold[j][0])*0.5;
   xone[1] = (xold[i][1] + xold[j][1])*0.5;
   xone[2] = (xold[i][2] + xold[j][2])*0.5;

   // record longest bond for ghostcut
    delx = xold[j][0] - xold[i][0];
    dely = xold[j][1] - xold[i][1];
    delz = xold[j][2] - xold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if(rsq > rsqold) rsqold = rsq; 

   // make one particle for each bond
   // i is local
   // if newton bond, always make particle
   // if j is local, always make particle
   // if j is ghost, decide from tag  

   if( force->newton_bond || j < nlocal || tagold[i] > tagold[j] ){
        atom->natoms += 1;
        atom->avec->create_atom(bptype,xone);
        // pack tag i/j into buffer for comm
        array[atom->nlocal-1][0] =  (double)tagold[i];
        array[atom->nlocal-1][1] =  (double)tagold[j];
        nadd++;
    }
  } 

  // new total # of atoms
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // assign tags for new atoms, update map
  if (atom->tag_enable) atom->tag_extend();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  char str[128];
  int Nadd = 0;
  MPI_Allreduce(&nadd,&Nadd,1,MPI_INT,MPI_SUM,world);
  if(comm->me == 0){
    sprintf(str, "Inserted %d bond particles.", Nadd);
    error->message(FLERR,str);
  }

  // check ghost comm distances
  // warn and change if shorter from estimate
  // ghost atoms must be present for bonds on edge of neighbor cutoff
  // extend cutghost slightly more than half of the longest bond
  MPI_Allreduce(&rsqold,&rsqmax,1,MPI_DOUBLE,MPI_MAX,world);
  rmax = sqrt(rsqmax); 
  double cutneighmax_srp = neighbor->cutneighmax + 0.51*rmax;

  // find smallest cutghost
  double cutghostmin = comm->cutghost[0];
  if (cutghostmin > comm->cutghost[1])
    cutghostmin = comm->cutghost[1];
  if (cutghostmin > comm->cutghost[2])
    cutghostmin = comm->cutghost[2];

  // reset cutghost if needed
  if (cutneighmax_srp > cutghostmin){
    if(comm->me == 0){
      sprintf(str, "Extending ghost comm cutoff. New %f, old %f.", cutneighmax_srp, cutghostmin);
      error->message(FLERR,str);
    }
    // cutghost updated by comm->setup 
    comm->cutghostuser = cutneighmax_srp;
  }

  // put new particles in the box before exchange
  // move owned to new procs
  // get ghosts
  // build neigh lists again
  domain->pbc();
  comm->setup();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  neighbor->build();
  neighbor->ncalls = 0;

  // new atom counts
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  // zero all forces 
  for(int i = 0; i < nall; i++)
    for(int n = 0; n < 3; n++)
      atom->f[i][n] = 0.0;

  // do not include bond particles in thermo output
  // remove them from all groups
  for(int i=0; i< nlocal; i++)
    if(atom->type[i] == bptype)
      atom->mask[i] = 0;
}

/* ----------------------------------------------------------------------
   set position of bond particles
------------------------------------------------------------------------- */

void FixSRP::pre_exchange()
{
  // update ghosts 
  comm->forward_comm();

  // reassign bond particle coordinates to midpoint of bonds
  // only needed before neigh rebuild
  double **x=atom->x;
  int i,j;
  int nlocal = atom->nlocal;

  for(int ii = 0; ii < nlocal; ii++){
    if(atom->type[ii] != bptype) continue;

    i = atom->map((int)array[ii][0]);
    if(i < 0) error->all(FLERR,"Fix SRP failed to map atom");
    i = domain->closest_image(ii,i);

    j = atom->map((int)array[ii][1]);
    if(j < 0) error->all(FLERR,"Fix SRP failed to map atom");
    j = domain->closest_image(ii,j);

    // position of bond particle ii 
    atom->x[ii][0] = (x[i][0] + x[j][0])*0.5;
    atom->x[ii][1] = (x[i][1] + x[j][1])*0.5;
    atom->x[ii][2] = (x[i][2] + x[j][2])*0.5;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSRP::memory_usage()
{
  double bytes = atom->nmax*2 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSRP::grow_arrays(int nmax)
{
  memory->grow(array,nmax,2,"fix_srp:array");
  array_atom = array;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
   called when move to new proc
------------------------------------------------------------------------- */

void FixSRP::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < 2; m++)
    array[j][m] = array[i][m];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values
   called when atom is created
------------------------------------------------------------------------- */

void FixSRP::set_arrays(int i)
{
  array[i][0] = -1;
  array[i][1] = -1; 
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixSRP::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < 2; m++) buf[m] = array[i][m];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixSRP::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < 2; m++) array[nlocal][m] = buf[m];
  return 2;
}
/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSRP::pack_border(int n, int *list, double *buf)
{
  // pack buf for border com
  int i,j;
  int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = array[j][0];
      buf[m++] = array[j][1];
    }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSRP::unpack_border(int n, int first, double *buf)
{
  // unpack buf into array
  int i,last;
  int m = 0;
  last = first + n;

      for (i = first; i < last; i++){
        array[i][0] = static_cast<int> (buf[m++]);
        array[i][1] = static_cast<int> (buf[m++]);
      }
  return m;
}

/* ----------------------------------------------------------------------
   remove particles after run
------------------------------------------------------------------------- */

void FixSRP::post_run()
{
  // all bond particles are removed after each run
  // useful for write_data and write_restart commands 
  // since it occurs between runs
 
  bigint natoms_previous = atom->natoms;
  int nlocal = atom->nlocal;
  int* dlist;
  memory->create(dlist,nlocal,"fix_srp:dlist");

  for (int i = 0; i < nlocal; i++){
    if(atom->type[i] == bptype)
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  // delete local atoms flagged in dlist
  // reset nlocal

  AtomVec *avec = atom->avec;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  // reset atom->natoms
  // reset atom->map if it exists
  // set nghost to 0 so old ghosts won't be mapped

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // print before and after atom count

  bigint ndelete = natoms_previous - atom->natoms;

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Deleted " BIGINT_FORMAT
                        " atoms, new total = " BIGINT_FORMAT "\n",
                        ndelete,atom->natoms);
    if (logfile) fprintf(logfile,"Deleted " BIGINT_FORMAT
                         " atoms, new total = " BIGINT_FORMAT "\n",
                         ndelete,atom->natoms);
  }

  // verlet calls box_too_small_check() in post_run
  // this check maps all bond partners
  // therefore need ghosts
  domain->pbc();
  comm->setup();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
} 

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixSRP::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 3;
  buf[m++] = array[i][0]; 
  buf[m++] = array[i][1];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixSRP::unpack_restart(int nlocal, int nth) 
{
  double **extra = atom->extra;

// skip to Nth set of extra values
  int m = 0; 
  for (int i = 0; i < nth; i++){
    m += static_cast<int> (extra[nlocal][m]);
  }

  m++; 
  array[nlocal][0] = extra[nlocal][m++];
  array[nlocal][1] = extra[nlocal][m++];

}
/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixSRP::maxsize_restart()
{
  return 3;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixSRP::size_restart(int nlocal)
{
  return 3;
}

/* ----------------------------------------------------------------------
   pack global state of Fix 
------------------------------------------------------------------------- */

void FixSRP::write_restart(FILE *fp)
{
  int n = 0;
  double list[3];
  list[n++] = comm->cutghostuser;
  list[n++] = btype;
  list[n++] = bptype;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixSRP::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  comm->cutghostuser = static_cast<double> (list[n++]);
  btype = static_cast<int> (list[n++]);
  bptype = static_cast<int> (list[n++]);
}

/* ----------------------------------------------------------------------
   interface with other classes
   pair srp sets the bond type for this fix
------------------------------------------------------------------------- */

int FixSRP::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"btype") == 0) {
    btype = atoi(arg[1]);
    return 2;
  }
  return 0;
}

