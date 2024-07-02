/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_surface_local.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"
#include "molecule.h"

#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;

// NOTE: epsilon should be defined as fraction of max surf size

#define EPSILON 0.001
#define NBIN 100
#define BIG 1.0e20
#define MAXLINE 256
#define MAXTRIPOINT 24

#define DELTA_BONUS 4            // NOTE: make it larger when done testing

enum{DATAFILE,MOLFILE};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

// allocate space for static class variable

FixSurfaceLocal *FixSurfaceLocal::fptr;

/* ---------------------------------------------------------------------- */

FixSurfaceLocal::FixSurfaceLocal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix surface/local command");

  create_attribute = 1;

  dimension = domain->dimension;

  index = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(2);

  nlocal_connect = nghost_connect = nmax_connect = 0;
  connect2d = NULL;
  connect3d = NULL;
  
  // NOTE: what is optimal way to initialize MyPoolChunk?

  tcp = NULL;
  if (dimension == 3) tcp = new MyPoolChunk<tagint>(1,MAXTRIPOINT,6);

  if (sizeof(double) == sizeof(tagint)) tagintdoubleratio = 1;
  else if (sizeof(double) == 2*sizeof(tagint)) tagintdoubleratio = 2;
  else error->all(FLERR,"Internal error in fix surface/local");

  // process args

  idmol = NULL;

  if (strcmp(arg[3],"NULL") == 0) mode = DATAFILE;
  else {
    mode = MOLFILE;
    int n = strlen(arg[3]) + 1;
    idmol = new char[n];
    strcpy(idmol,arg[3]);
  }
}

/* ---------------------------------------------------------------------- */

FixSurfaceLocal::~FixSurfaceLocal()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete [] idmol;
  memory->destroy(index);

  if (dimension == 2) memory->sfree(connect2d);
  else {
    int nall = nlocal_connect + nghost_connect;
    for (int i = 0; i < nall; i++) {
      tcp->put(connect3d[i].indexc1);
      tcp->put(connect3d[i].indexc2);
      tcp->put(connect3d[i].indexc3);
    }
    memory->sfree(connect3d);
  }

  delete tcp;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceLocal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;     // NOTE: just for debugging
  return mask;
}

/* ----------------------------------------------------------------------
   process and setup lines/tris
   NOTE: doc why done here
------------------------------------------------------------------------- */

void FixSurfaceLocal::post_constructor()
{
  // for DATAFILE mode:
  // lines or tri are already defined from reading data file(s)
  // already distributed to procs
  // task is to infer line/tri connectivity in distributed manner

  if (mode == DATAFILE) {
    if (comm->me == 0 && screen)
      fprintf(screen,"Connecting line/tri particles ...\n");

    if (dimension == 2) connectivity2d_local();
    else connectivity3d_local();
  }

  // for SURFFILE mode:
  // each proc builds global data struct from one or more molecules
  // all procs infer global connectivity from global data struct
  // assign lines/surfs to each proc, based on center pt coords

  if (mode == MOLFILE) {
    if (comm->me == 0 && screen)
      fprintf(screen,"Converting molecule file to line/tri particles ...\n");

    extract_from_molecules(idmol);

    if (dimension == 2) {
      connectivity2d_global();
      assign2d();
    } else {
      connectivity3d_global();
      assign3d();
    }

    memory->sfree(points);
    memory->sfree(lines);
    memory->sfree(tris);
    if (dimension == 3) memory->destroy(clist);
  }

  // set max size of connection info
  // NOTE: set 3d based on # max of corner connections (not 12), is 10 correct?

  if (dimension == 2) comm_border = 4;
  else comm_border = 10 + 3*12;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_arrays(int nmax)
{
  memory->grow(index,nmax,"surface/local:index");
}

/* ----------------------------------------------------------------------
   DEBUG: see pre_neighbor
------------------------------------------------------------------------- */

void FixSurfaceLocal::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   DEBUG: check that index vector is identical to line/tri
   check that index vector is consistent with connect->ilocal
------------------------------------------------------------------------- */

void FixSurfaceLocal::pre_neighbor()
{
  int n = atom->nlocal + atom->nghost;

  int *surf = atom->line;
  if (dimension == 3) surf = atom->tri;

  int count1 = 0;
  int count2 = 0;

  // DEBUG
  //if (comm->me == 0) printf("RENEIGH %ld\n",update->ntimestep);

  for (int i = 0; i < n; i++) {
    if (surf[i] < 0) continue;

    if (index[i] != surf[i]) {
      count1++;
      continue;
    }
    if (dimension == 2 && connect2d[index[i]].ilocal != i) count2++;
    if (dimension == 3 && connect3d[index[i]].ilocal != i) count2++;
  }

  int all1,all2;
  MPI_Allreduce(&count1,&all1,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&count2,&all2,1,MPI_INT,MPI_SUM,world);

  if (all1 || all2 && comm->me == 0) {
    char str[128];
    sprintf(str,"FSL index vector mis-match: %d %d\n",all1,all2);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   grow connect array by DELTA
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_connect() 
{
  nmax_connect = nmax_connect/DELTA_BONUS * DELTA_BONUS;
  nmax_connect += DELTA_BONUS;
  if (nmax_connect < 0) error->one(FLERR,"Per-processor system is too big");

  if (dimension == 2) 
    connect2d = (Connect2d *) 
      memory->srealloc(connect2d,nmax_connect*sizeof(Connect2d),
                       "surface/local:connect2d");
  else
    connect3d = (Connect3d *) 
      memory->srealloc(connect3d,nmax_connect*sizeof(Connect3d),
                       "surface/local:connect3d");
}

/* ----------------------------------------------------------------------
   copy values within local connect array
------------------------------------------------------------------------- */

void FixSurfaceLocal::copy_arrays(int i, int j, int delflag)
{
  // if deleting atom J via delflag and J has connect data, then delete it
  
  if (dimension == 2) {
    if (delflag && index[j] >= 0) {
      int k = index[j];
      copy_connect(nlocal_connect-1,k);
      nlocal_connect--;
    }
  } else {
    if (delflag && index[j] >= 0) {
      int k = index[j];
      tcp->put(connect3d[k].indexc1);
      tcp->put(connect3d[k].indexc2);
      tcp->put(connect3d[k].indexc3);
      copy_connect(nlocal_connect-1,k);
      nlocal_connect--;
    }
  }

  // if atom I has connect data, reset I's connect.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's connection data 
  //   is already deleted

  if (index[i] >= 0 && i != j) {
    if (dimension == 2) connect2d[index[i]].ilocal = j;
    else connect3d[index[i]].ilocal = j;
  }
  index[j] = index[i];
}

/* ----------------------------------------------------------------------
   copy bonus data from I to J, effectively deleting the J entry
   also reset index that points to I to now point to J
------------------------------------------------------------------------- */

void FixSurfaceLocal::copy_connect(int i, int j)
{
  if (dimension == 2) {
    index[connect2d[i].ilocal] = j;
    memcpy(&connect2d[j],&connect2d[i],sizeof(Connect2d));
  } else {
    index[connect3d[i].ilocal] = j;
    memcpy(&connect3d[j],&connect3d[i],sizeof(Connect3d));
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixSurfaceLocal::set_arrays(int i)
{
  index[i] = -1;
}

/* ----------------------------------------------------------------------
   clear ghost info in connect array
   called before ghosts are recommunicated in comm and irregular
------------------------------------------------------------------------- */

void FixSurfaceLocal::clear_bonus()
{
  if (dimension == 3) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = nlocal_connect; i < nall; i++) {
      tcp->put(connect3d[i].indexc1);
      tcp->put(connect3d[i].indexc2);
      tcp->put(connect3d[i].indexc3);
    }
  }
  nghost_connect = 0;
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSurfaceLocal::pack_border(int n, int *list, double *buf)
{
  int i,j,k,m,ic,nc;

  m = 0;

  if (dimension == 2) {

    for (i = 0; i < n; i++) {
      j = list[i];
      if (index[j] < 0) buf[m++] = ubuf(0).d;
      else {
        ic = index[j];
        buf[m++] = ubuf(1).d;
        buf[m++] = ubuf(connect2d[ic].neigh_p1).d;
        buf[m++] = ubuf(connect2d[ic].neigh_p2).d;
        buf[m++] = ubuf(connect2d[ic].flags).d;
      }
    }

  } else {

    for (i = 0; i < n; i++) {
      j = list[i];
      if (index[j] < 0) buf[m++] = ubuf(0).d;
      else {
        ic = index[j];
        buf[m++] = ubuf(1).d;
        buf[m++] = ubuf(connect3d[ic].neigh_e1).d;
        buf[m++] = ubuf(connect3d[ic].neigh_e2).d;
        buf[m++] = ubuf(connect3d[ic].neigh_e3).d;
        buf[m++] = ubuf(connect3d[ic].nc1).d;
        buf[m++] = ubuf(connect3d[ic].nc2).d;
        buf[m++] = ubuf(connect3d[ic].nc3).d;
        buf[m++] = ubuf(connect3d[ic].flags).d;
        
        // NOTE: which pack method to use
        
        nc = connect3d[ic].nc1;
        for (k = 0; k < nc; k++)
          buf[m++] = ubuf(connect3d[ic].neigh_c1[k]).d;
        nc = connect3d[ic].nc2;
        for (k = 0; k < nc; k++)
          buf[m++] = ubuf(connect3d[ic].neigh_c2[k]).d;
        nc = connect3d[ic].nc3;
        for (k = 0; k < nc; k++)
          buf[m++] = ubuf(connect3d[ic].neigh_c3[k]).d;
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSurfaceLocal::unpack_border(int n, int first, double *buf)
{
  int i,j,k,m,last,nc;

  m = 0;
  last = first + n;

  if (dimension == 2) {

    for (i = first; i < last; i++) {
      index[i] = (int) ubuf(buf[m++]).i;
      if (index[i] == 0) index[i] = -1;
      else {
        j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();
        connect2d[j].neigh_p1 = (tagint) ubuf(buf[m++]).i;
        connect2d[j].neigh_p2 = (tagint) ubuf(buf[m++]).i;
        connect2d[j].flags = (int) ubuf(buf[m++]).i;
        
        connect2d[j].ilocal = i;
        index[i] = j;
        nghost_connect++;
      }
    }

  } else {

    for (i = first; i < last; i++) {
      index[i] = (int) ubuf(buf[m++]).i;
      if (index[i] == 0) index[i] = -1;
      else {
        j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();
        
        connect3d[j].neigh_e1 = (tagint) ubuf(buf[m++]).i;
        connect3d[j].neigh_e2 = (tagint) ubuf(buf[m++]).i;
        connect3d[j].neigh_e3 = (tagint) ubuf(buf[m++]).i;
        connect3d[j].nc1 = (int) ubuf(buf[m++]).i;
        connect3d[j].nc2 = (int) ubuf(buf[m++]).i;
        connect3d[j].nc3 = (int) ubuf(buf[m++]).i;
        connect3d[j].flags = (int) ubuf(buf[m++]).i;
      
        connect3d[j].neigh_c1 = tcp->get(connect3d[j].nc1,connect3d[j].indexc1);
        connect3d[j].neigh_c2 = tcp->get(connect3d[j].nc2,connect3d[j].indexc2);
        connect3d[j].neigh_c3 = tcp->get(connect3d[j].nc3,connect3d[j].indexc3);
        
        // NOTE: which unpack method to use
      
        nc = connect3d[j].nc1;
        for (k = 0; k < nc; k++)
          connect3d[j].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
        nc = connect3d[j].nc2;
        for (k = 0; k < nc; k++)
          connect3d[j].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
        nc = connect3d[j].nc3;
        for (k = 0; k < nc; k++)
          connect3d[j].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;
        
        connect3d[j].ilocal = i;
        index[i] = j;
        nghost_connect++;
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local connect array for exchange with another proc
------------------------------------------------------------------------- */

int FixSurfaceLocal::pack_exchange(int i, double *buf)
{
  int j,k,n;
  int m = 0;

  if (dimension == 2) {
    if (index[i] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = index[i];

      buf[m++] = ubuf(connect2d[j].neigh_p1).d;
      buf[m++] = ubuf(connect2d[j].neigh_p2).d;
      buf[m++] = ubuf(connect2d[j].flags).d;
    }

  } else {
    if (index[i] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = index[i];

      buf[m++] = ubuf(connect3d[j].neigh_e1).d;
      buf[m++] = ubuf(connect3d[j].neigh_e2).d;
      buf[m++] = ubuf(connect3d[j].neigh_e3).d;
      buf[m++] = ubuf(connect3d[j].nc1).d;
      buf[m++] = ubuf(connect3d[j].nc2).d;
      buf[m++] = ubuf(connect3d[j].nc3).d;
      buf[m++] = ubuf(connect3d[j].flags).d;
      
      // NOTE: which method to use

      n = connect3d[j].nc1;
      for (k = 0; k < n; k++)
        buf[m++] = ubuf(connect3d[j].neigh_c1[k]).d;
      n = connect3d[j].nc2;
      for (k = 0; k < n; k++)
        buf[m++] = ubuf(connect3d[j].neigh_c2[k]).d;
      n = connect3d[j].nc3;
      for (k = 0; k < n; k++)
        buf[m++] = ubuf(connect3d[j].neigh_c3[k]).d;
      
      /*
        memcpy(&buf[m],connect3d[i].neigh_c1,connect3d[i].nc1*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[i].nc1;
        else m += (connect3d[i].nc1+1)/2;
        memcpy(&buf[m],connect3d[i].neigh_c2,connect3d[i].nc2*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[i].nc1;
        else m += (connect3d[i].nc2+1)/2;
        memcpy(&buf[m],connect3d[i].neigh_c3,connect3d[i].nc3*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[i].nc1;
        else m += (connect3d[i].nc3+1)/2;
      */
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local connect array from exchange with another proc
------------------------------------------------------------------------- */

int FixSurfaceLocal::unpack_exchange(int nlocal, double *buf)
{
  int j,k,n;
  int m = 0;

  if (dimension == 2) {
    index[nlocal] = (int) ubuf(buf[m++]).i;
    if (index[nlocal] == 0) index[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      connect2d[nlocal_connect].neigh_p1 = (tagint) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].neigh_p2 = (tagint) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].flags = (int) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].ilocal = nlocal;

      index[nlocal] = nlocal_connect++;
    }

  } else {
    index[nlocal] = (int) ubuf(buf[m++]).i;
    if (index[nlocal] == 0) index[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      connect3d[nlocal_connect].neigh_e1 = (tagint) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].neigh_e2 = (tagint) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].neigh_e3 = (tagint) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].nc1 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].nc2 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].nc3 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].flags = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].ilocal = nlocal;

      connect3d[nlocal_connect].neigh_c1 = 
        tcp->get(connect3d[nlocal_connect].nc1,
                 connect3d[nlocal_connect].indexc1);
      connect3d[nlocal_connect].neigh_c2 = 
        tcp->get(connect3d[nlocal_connect].nc2,
                 connect3d[nlocal_connect].indexc2);
      connect3d[nlocal_connect].neigh_c3 = 
        tcp->get(connect3d[nlocal_connect].nc3,
                 connect3d[nlocal_connect].indexc3);

      // NOTE: which method to use

      n = connect3d[nlocal_connect].nc1;
      for (k = 0; k < n; k++)
        connect3d[nlocal_connect].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
      n = connect3d[nlocal_connect].nc2;
      for (k = 0; k < n; k++)
        connect3d[nlocal_connect].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
      n = connect3d[nlocal_connect].nc3;
      for (k = 0; k < n; k++)
        connect3d[nlocal_connect].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;

      /*
        memcpy(connect3d[nlocal].neigh_c1,&buf[m],
        connect3d[nlocal].nc1*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[nlocal].nc1;
        else m += (connect3d[nlocal].nc1+1)/2;
        memcpy(connect3d[nlocal].neigh_c2,&buf[m],
        connect3d[nlocal].nc2*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[nlocal].nc2;
        else m += (connect3d[nlocal].nc2+1)/2;
        memcpy(connect3d[nlocal].neigh_c3,&buf[m],
        connect3d[nlocal].nc3*sizeof(tagint));
        if (tagintdoubleratio == 1) m += connect3d[nlocal].nc3;
        else m += (connect3d[nlocal].nc3+1)/2;
      */

      index[nlocal] = nlocal_connect++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSurfaceLocal::memory_usage()
{
  double bytes = 0.0;
  if (dimension == 2) bytes = nmax_connect*sizeof(Connect2d);
  else bytes = nmax_connect*sizeof(Connect3d);
  bytes += atom->nmax * sizeof(int);               // index array
  return bytes;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for distributed 2d connectivity build
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for owned lines
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity2d_local()
{
  int i,j,k,m,n,ibin;

  // error check

  avec_line = (AtomVecLine *) atom->style_match("line");
  if (!avec_line) 
    error->all(FLERR,"Fix surface/local requires atom style line");

  // calculate current endpts of owned lines

  int nlocal = atom->nlocal;
  memory->create(endpts,nlocal,4,"surface/local:endpts");
  calculate_endpts(nlocal);

  // nline = # of line segments I own
  // maxline = max of nline across all procs
  // ebuf = array of per-line values to circulate to all procs

  int *line = atom->line;

  int nline = 0;
  for (i = 0; i < nlocal; i++)
    if (line[i] >= 0) nline++;

  int maxline;
  MPI_Allreduce(&nline,&maxline,1,MPI_INT,MPI_MAX,world);

  // initialize ebuf for lines I own
  // 5 values = id, 2 end pt coords (x,y)

  double **ebuf;
  memory->create(ebuf,maxline,5,"surface/local:ebuf");

  AtomVecLine::Bonus *bonus = avec_line->bonus;
  double **x = atom->x;
  tagint *tag = atom->tag;

  nline = 0;
  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    ebuf[nline][0] = ubuf(tag[i]).d;
    memcpy(&ebuf[nline][1],endpts[i],4*sizeof(double));
    nline++;
  }

  // epssq = square of EPSILON fraction of min line length across all procs

  double minlen = BIG;
  for (i = 0; i < nlocal; i++)
    if (line[i] >= 0) minlen = MIN(minlen,bonus[line[i]].length);

  double eps;
  MPI_Allreduce(&minlen,&eps,1,MPI_DOUBLE,MPI_MIN,world);
  eps *= EPSILON;
  epssq = eps*eps;

  // compute bbox for my line endpts, add 2*EPS to all 4 bounds
  // set bbox size to 0.0 +/- EPS if no lines on this proc

  binlo[0] = binlo[1] = BIG;
  binhi[0] = binhi[1] = -BIG;

  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    binlo[0] = MIN(binlo[0],endpts[i][0]);
    binhi[0] = MAX(binhi[0],endpts[i][0]);
    binlo[1] = MIN(binlo[1],endpts[i][1]);
    binhi[1] = MAX(binhi[1],endpts[i][1]);
    binlo[0] = MIN(binlo[0],endpts[i][2]);
    binhi[0] = MAX(binhi[0],endpts[i][2]);
    binlo[1] = MIN(binlo[1],endpts[i][3]);
    binhi[1] = MAX(binhi[1],endpts[i][3]);
  }

  if (binlo[0] > binhi[0]) {
    binlo[0] = binhi[0] = 0.0;
    binlo[1] = binhi[1] = 0.0;
  }

  binlo[0] -= 2.0*eps;
  binlo[1] -= 2.0*eps;
  binhi[0] += 2.0*eps;
  binhi[1] += 2.0*eps;

  // create upto NBIN x NBIN bins that tile bbox
  // insure bin size is not <= 2*EPS so that pt +/- EPS cannot overlap > 4 bins

  nbinx = static_cast<int> ((binhi[0]-binlo[0])/(4.0*eps));
  nbinx = MIN(nbinx,NBIN);
  nbinx = MAX(nbinx,1);
  nbiny = static_cast<int> ((binhi[1]-binlo[1])/(4.0*eps));
  nbiny = MIN(nbiny,NBIN);
  nbiny = MAX(nbiny,1);

  invbinx = nbinx / (binhi[0] - binlo[0]);
  invbiny = nbiny / (binhi[1] - binlo[1]);

  nbins = nbinx*nbiny;
  memory->create(bincount,nbins,"surface/local:bincount");
  memory->create(binfirst,nbins,"surface/local:binfirst");

  // count # of end pts in each bin, including overlaps by eps

  int indices[4];
  memset(bincount,0,nbins*sizeof(int));

  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    n = overlap2bin2d(&endpts[i][0],eps,indices);
    for (m = 0; m < n; m++) bincount[indices[m]]++;
    n = overlap2bin2d(&endpts[i][2],eps,indices);
    for (m = 0; m < n; m++) bincount[indices[m]]++;
  }

  // setup binfirst = index to first point in bin
  // allocate pts = list of pts in all bins

  binfirst[0] = 0;
  for (m = 1; m < nbins; m++)
    binfirst[m] = binfirst[m-1] + bincount[m-1];
  int ntotal = binfirst[nbins-1] + bincount[nbins-1];

  pts = (OnePt *) memory->smalloc(ntotal*sizeof(OnePt),"surface/local:bins");

  // add each of my line endpts to bins, including overlaps by eps

  memset(bincount,0,nbins*sizeof(int));

  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    for (j = 0; j < 2; j++) {
      if (j == 0) n = overlap2bin2d(&endpts[i][0],eps,indices);
      else n = overlap2bin2d(&endpts[i][2],eps,indices);
      for (m = 0; m < n; m++) {
        ibin = indices[m];
        k = binfirst[ibin] + bincount[ibin];
        pts[k].iatom = i;
        pts[k].iconnect = line[i];
        pts[k].ptwhich = j+1;
        bincount[ibin]++;
      }
    }
  }

  // allocate & initialize my line info as if no connections
  // comm->ring() operation will fill it in with actual connections

  nlocal_connect = nmax_connect = nline;
  nghost_connect = 0;
  grow_connect();

  for (i = 0; i < nlocal; i++) {
    index[i] = line[i];
    if (line[i] < 0) continue;
    j = line[i];
    connect2d[j].neigh_p1 = connect2d[j].neigh_p2 = 0;
    connect2d[j].flags = 0;
    connect2d[j].ilocal = i;
  }

  // circulate my line info to all other procs, including self

  fptr = this;
  nmatch = errormatch = 0;
  
  /*
  if (ebuf)
    comm->ring(nline,5*sizeof(double),(void *) &ebuf[0][0],0,linematch,NULL);
  else 
    comm->ring(nline,5*sizeof(double),NULL,0,linematch,NULL);
  */
  
  // check for errors = matches with more than 2 points
  // print stats on # of matches

  int all = 0;
  MPI_Allreduce(&errormatch,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Fix surface/local found %g matching end pts "
            "with more than 2 lines",0.5*all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nmatch,&all,1,MPI_INT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"  matched %g line end point pairs\n",0.5*all);
    if (logfile) 
      fprintf(logfile,"  matched %g line end point pairs\n",0.5*all);
  }

  // clean-up

  memory->destroy(ebuf);
  memory->destroy(bincount);
  memory->destroy(binfirst);
  memory->destroy(pts);
  memory->destroy(endpts);
}

/* ----------------------------------------------------------------------
   callback from comm->ring() for 2d
------------------------------------------------------------------------- */

void FixSurfaceLocal::linematch(int n, char *cbuf)
{
  int i,j,m,ibin,npt,jfirst,jlast,iatom,iconnect,ptwhich,optwhich;
  tagint tagother;
  double ptx,pty,dx,dy,rsq;
  double *pt,*opt;

  double *buf = (double *) cbuf;

  // access class data for neighbors and bins

  Connect2d *connect2d = fptr->connect2d;
  double **endpts = fptr->endpts;

  OnePt *pts = fptr->pts;
  int *bincount = fptr->bincount;
  int *binfirst = fptr->binfirst;

  double epssq = fptr->epssq;
  int *nmatch = &fptr->nmatch;
  int *errormatch = &fptr->errormatch;

  tagint *tag = fptr->atom->tag;

  // loop over each entry in buf
  // buf = end pts for each line owned by another proc (including self)
  // 0 = ID of line, 1/2 = 1st end pt, 3/4 = 2nd end pt

  m = 0;
  for (i = 0; i < n; i++) {
    tagother = (tagint) ubuf(buf[m]).i;

    // loop over 2 pts in buf line, owned by another proc (or me)
    // search for matching owned line pt in ibin
    // same owned pt may appear many times in bins as part of multiple lines

    for (optwhich = 1; optwhich <= 2; optwhich++) {
      if (optwhich == 1) opt = &buf[m+1];
      else opt = &buf[m+3];

      ibin = fptr->pt2bin2d(opt);

      if (ibin >= 0) {
        jfirst = binfirst[ibin];
        jlast = jfirst + bincount[ibin];

        // loop over all line end pts in bins[jfirst:jlast], owned by me
        // skip bin pt if from same line ID as buf pt is from

        for (j = jfirst; j < jlast; j++) {
          iatom = pts[j].iatom;
          iconnect = pts[j].iconnect;
          ptwhich = pts[j].ptwhich;

          if (tagother == tag[iatom]) continue;

          if (ptwhich == 1) pt = &endpts[iconnect][0];
          else pt = &endpts[iconnect][2];

          dx = pt[0] - opt[0];
          dy = pt[1] - opt[1];
          rsq = dx*dx + dy*dy;

          // two points are close enough to match
          // set neigh1 or neigh2 of my line = tagother, based on ptwhich
          // error if this point is already common to 2 lines
          // NOTE: this is where to set coupling flags

          if (rsq < epssq) {
            (*nmatch)++;

            if (ptwhich == 1) {
              if (connect2d[iconnect].neigh_p1 > 0) (*errormatch)++;
              else connect2d[iconnect].neigh_p1 = tagother;
            } else if (ptwhich == 2) {
              if (connect2d[iconnect].neigh_p2 > 0) (*errormatch)++;
              else connect2d[iconnect].neigh_p2 = tagother;
            }
          }
        }
      }
    }

    m += 5;
  }
}

/* ----------------------------------------------------------------------
   compute current end points of N lines
   nothing computed for particles that are not lines
------------------------------------------------------------------------- */

void FixSurfaceLocal::calculate_endpts(int n)
{
  double length,theta,dx,dy;
  double *endpt;

  AtomVecLine::Bonus *bonus = avec_line->bonus;
  double **x = atom->x;
  int *line = atom->line;

  for (int i = 0; i < n; i++) {
    if (line[i] < 0) continue;
    endpt = endpts[i];
    length = bonus[line[i]].length;
    theta = bonus[line[i]].theta;
    dx = 0.5*length*cos(theta);
    dy = 0.5*length*sin(theta);
    endpt[0] = x[i][0] - dx;
    endpt[1] = x[i][1] - dy;
    endpt[2] = x[i][0] + dx;
    endpt[3] = x[i][1] + dy;
  }
}

/* ----------------------------------------------------------------------
   map point pt to a bin
   return ibin from 0 to Nbins-1
   return -1 if outside bin bounds
------------------------------------------------------------------------- */

int FixSurfaceLocal::pt2bin2d(double *pt)
{
  if (pt[0] < binlo[0] || pt[0] >= binhi[0] || 
      pt[1] < binlo[1] || pt[1] >= binhi[1])
    return -1;

  int ix = static_cast<int> ((pt[0]-binlo[0]) * invbinx);
  ix = MIN(ix,nbinx-1);
  int iy = static_cast<int> ((pt[1]-binlo[1]) * invbiny);
  iy = MIN(iy,nbiny-1);

  int ibin = iy*nbinx + ix;
  return ibin;
}

/* ----------------------------------------------------------------------
   map point +/ EPS in both dims to all bins it overlaps with
   return N = # of overlapped bins
   return indices = list of bin IDs, each from 0 to Nbins-1
------------------------------------------------------------------------- */

int FixSurfaceLocal::overlap2bin2d(double *pt, double eps, int *indices)
{
  int ilo = static_cast<int> ((pt[0]-eps-binlo[0]) * invbinx);
  int ihi = static_cast<int> ((pt[0]+eps-binlo[0]) * invbinx);
  int jlo = static_cast<int> ((pt[1]-eps-binlo[1]) * invbiny);
  int jhi = static_cast<int> ((pt[1]+eps-binlo[1]) * invbiny);

  ilo = MAX(ilo,0);
  ihi = MIN(ihi,nbinx-1);
  jlo = MAX(jlo,0);
  jhi = MIN(jhi,nbiny-1);

  if (ilo == ihi && jlo == jhi) {
    indices[0] = jlo*nbinx + ilo;
    return 1;
  }

  int i,j;
  int n = 0;
  for (j = jlo; j <= jhi; j++) 
    for (i = ilo; i <= ihi; i++)
      indices[n++] = j*nbinx + i;
  return n;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for distributed 3d connectivity build
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all owned tris
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity3d_local()
{
  int i,j,k,m,n,ibin;

  // error check

  avec_tri = (AtomVecTri *) atom->style_match("tri");
  if (!avec_tri) 
    error->all(FLERR,"Fix surface/local requires atom style tri");

  // calculate current corner pts of owned triangles

  int nlocal = atom->nlocal;
  memory->create(corners,nlocal,9,"surface/local:corners");
  calculate_corners(nlocal);

  // ntri = # of triangles I own
  // maxtri = max of ntri across all procs
  // tbuf = array of per-tri values to circulate to all procs

  int *tri = atom->tri;

  int ntri = 0;
  for (i = 0; i < nlocal; i++)
    if (tri[i] >= 0) ntri++;

  int maxtri;
  MPI_Allreduce(&ntri,&maxtri,1,MPI_INT,MPI_MAX,world);

  // initialize tbuf for tris I own
  // 10 values = id, 3 corner pt coords (x,y,z)

  double **tbuf;
  memory->create(tbuf,maxtri,10,"surface/local:tbuf");

  AtomVecTri::Bonus *bonus = avec_tri->bonus;
  double **x = atom->x;
  tagint *tag = atom->tag;

  ntri = 0;
  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    tbuf[ntri][0] = ubuf(tag[i]).d;
    memcpy(&tbuf[ntri][1],corners[i],9*sizeof(double));
    ntri++;
  }

  // epssq = square of EPSILON fraction of min tri diameter across all procs

  double *radius = atom->radius;

  double minlen = BIG;
  for (i = 0; i < nlocal; i++)
    if (tri[i] >= 0) minlen = MIN(minlen,radius[i]);
  minlen *= 2.0;

  double eps;
  MPI_Allreduce(&minlen,&eps,1,MPI_DOUBLE,MPI_MIN,world);
  eps *= EPSILON;
  epssq = eps*eps;

  // compute bbox for my tri corner pts, add 2*EPS to all 6 bounds
  // set bbox size to 0.0 +/- EPS if no tris on this proc

  binlo[0] = binlo[1] = binlo[2] = BIG;
  binhi[0] = binhi[1] = binhi[2] = -BIG;

  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    m = 0;
    for (j = 0; j < 3; j++) {
      binlo[0] = MIN(binlo[0],corners[i][m]);
      binhi[0] = MAX(binhi[0],corners[i][m]);
      binlo[1] = MIN(binlo[1],corners[i][m+1]);
      binhi[1] = MAX(binhi[1],corners[i][m+1]);
      binlo[2] = MIN(binlo[2],corners[i][m+2]);
      binhi[2] = MAX(binhi[2],corners[i][m+2]);
      m += 3;
    }
  }

  if (binlo[0] > binhi[0]) {
    binlo[0] = binhi[0] = 0.0;
    binlo[1] = binhi[1] = 0.0;
    binlo[2] = binhi[2] = 0.0;
  }

  binlo[0] -= 2.0*eps;
  binlo[1] -= 2.0*eps;
  binlo[2] -= 2.0*eps;
  binhi[0] += 2.0*eps;
  binhi[1] += 2.0*eps;
  binhi[2] += 2.0*eps;

  // create upto NBIN x NBIN x NBIN bins that tile bbox
  // insure bin size is not <= 2*EPS so that pt +/- EPS cannot overlap > 8 bins

  nbinx = static_cast<int> ((binhi[0]-binlo[0])/(4.0*eps));
  nbinx = MIN(nbinx,NBIN);
  nbinx = MAX(nbinx,1);
  nbiny = static_cast<int> ((binhi[1]-binlo[1])/(4.0*eps));
  nbiny = MIN(nbiny,NBIN);
  nbiny = MAX(nbiny,1);
  nbinz = static_cast<int> ((binhi[2]-binlo[2])/(4.0*eps));
  nbinz = MIN(nbinz,NBIN);
  nbinz = MAX(nbinz,1);

  invbinx = nbinx / (binhi[0] - binlo[0]);
  invbiny = nbiny / (binhi[1] - binlo[1]);
  invbinz = nbinz / (binhi[2] - binlo[2]);

  nbins = nbinx*nbiny*nbinz;
  memory->create(bincount,nbins,"surface/local:bincount");
  memory->create(binfirst,nbins,"surface/local:binfirst");

  // count # of corner pts in each bin, including overlaps by eps

  int indices[8];
  memset(bincount,0,nbins*sizeof(int));

  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    n = overlap2bin3d(&corners[i][0],eps,indices);
    for (m = 0; m < n; m++) bincount[indices[m]]++;
    n = overlap2bin3d(&corners[i][3],eps,indices);
    for (m = 0; m < n; m++) bincount[indices[m]]++;
    n = overlap2bin3d(&corners[i][6],eps,indices);
    for (m = 0; m < n; m++) bincount[indices[m]]++;
  }

  // setup binfirst = index to first point in bin
  // allocate pts = list of pts in all bins

  binfirst[0] = 0;
  for (m = 1; m < nbins; m++)
    binfirst[m] = binfirst[m-1] + bincount[m-1];
  int ntotal = binfirst[nbins-1] + bincount[nbins-1];

  pts = (OnePt *) memory->smalloc(ntotal*sizeof(OnePt),"surface/local:bins");

  // add each of my tri corners to bins, including overlaps by eps

  memset(bincount,0,nbins*sizeof(int));

  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    for (j = 0; j < 3; j++) {
      if (j == 0) n = overlap2bin3d(&corners[i][0],eps,indices);
      else if (j == 1) n = overlap2bin3d(&corners[i][3],eps,indices);
      else n = overlap2bin3d(&corners[i][6],eps,indices);
      for (m = 0; m < n; m++) {
        ibin = indices[m];
        k = binfirst[ibin] + bincount[ibin];
        pts[k].iatom = i;
        pts[k].iconnect = tri[i];
        pts[k].ptwhich = j+1;
        bincount[ibin]++;
      }
    }
  }

  // allocate & initialize my tri info as if no connections
  // circulate my tri info to all other procs, including self
  // do this twice:
  //   first time to accumulate corner counts
  //   allocate corner connection lists
  //   second time to fill lists (and redo everything else)
  // initialize my tri info as if no connections both times
  // comm->ring() operation will fill it in with actual connections

  nlocal_connect = nmax_connect = ntri;
  nghost_connect = 0;
  grow_connect();

  for (i = 0; i < nlocal; i++) {
    index[i] = tri[i];
    if (tri[i] < 0) continue;
    j = tri[i];
    connect3d[j].neigh_e1 = connect3d[j].neigh_e2 = connect3d[j].neigh_e3 = 0;
    connect3d[j].nc1 = connect3d[j].nc2 = connect3d[j].nc3 = 0;
    connect3d[j].flags = 0;
  }

  fptr = this;
  nmatch1 = nmatch2 = errormatch1 = errormatch2 = 0;
  vecflag = 0;
  
  /*
  if (tbuf) 
    comm->ring(ntri,10*sizeof(double),(void *) &tbuf[0][0],0,trimatch,NULL);
  else 
    comm->ring(ntri,10*sizeof(double),NULL,0,trimatch,NULL);
  */
  
  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    j = tri[i];
    connect3d[j].indexc1 = connect3d[j].indexc2 = connect3d[j].indexc3 = -1;
    if (connect3d[j].nc1) 
      connect3d[j].neigh_c1 = tcp->get(connect3d[j].nc1,connect3d[j].indexc1);
    else connect3d[j].neigh_c1 = NULL;
    if (connect3d[j].nc2) 
      connect3d[j].neigh_c2 = tcp->get(connect3d[j].nc2,connect3d[j].indexc2);
    else connect3d[j].neigh_c2 = NULL;
    if (connect3d[j].nc3) 
      connect3d[j].neigh_c3 = tcp->get(connect3d[j].nc3,connect3d[j].indexc3);
    else connect3d[j].neigh_c3 = NULL;
    connect3d[j].neigh_e1 = connect3d[j].neigh_e2 = connect3d[j].neigh_e3 = 0;
    connect3d[j].nc1 = connect3d[j].nc2 = connect3d[j].nc3 = 0;
    connect3d[j].flags = 0;
    connect3d[j].ilocal = i;
  }

  fptr = this;
  nmatch1 = nmatch2 = errormatch1 = errormatch2 = 0;
  vecflag = 1;
  
  /*
  if (tbuf) 
    comm->ring(ntri,10*sizeof(double),(void *) &tbuf[0][0],0,trimatch,NULL);
  else 
    comm->ring(ntri,10*sizeof(double),NULL,0,trimatch,NULL);
  */
  
  // check for errors = matches with misaligned or more than 2 edges
  // print stats on # of matches

  int all = 0;
  MPI_Allreduce(&errormatch1,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Fix surface/local found %g matching edges in same direction",
            0.5*all);
    error->all(FLERR,str);
  }

  all = 0;
  MPI_Allreduce(&errormatch2,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Fix surface/local found %g matching edges "
            "with more than 2 tris",0.5*all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nmatch1,&all,1,MPI_INT,MPI_SUM,world);
  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,
              "  matched %g tri edge pairs in fix surface/local\n",0.5*all);
    if (logfile) 
      fprintf(logfile,
              "  matched %g tri edge pairs in fix surface/local\n",0.5*all);
  }

  MPI_Allreduce(&nmatch2,&all,1,MPI_INT,MPI_SUM,world);
  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,
              "  matched %g tri corners in fix surface/local\n",0.5*all);
    if (logfile) 
      fprintf(logfile,
              "  matched %g tri corners in fix surface/local\n",0.5*all);
  }

  // clean-up

  memory->destroy(tbuf);
  memory->destroy(bincount);
  memory->destroy(binfirst);
  memory->destroy(pts);
  memory->destroy(corners);
}

/* ----------------------------------------------------------------------
   callback from comm->ring() for 3d
------------------------------------------------------------------------- */

void FixSurfaceLocal::trimatch(int n, char *cbuf)
{
  int i,j,m,ibin,npt,ipoint,jfirst,jlast,iatom,iconnect;
  int ptwhich,optwhich,match;
  tagint tagother;
  double ptx,pty,dx,dy,dz,rsq,rsq11,rsq12,rsq21,rsq22;
  double *pt,*ptnext1,*ptnext2,*opt,*optnext1,*optnext2;
  double *cornersother;

  double *buf = (double *) cbuf;

  // access class data for neighbors and bins

  Connect3d *connect3d = fptr->connect3d;
  double **corners = fptr->corners;

  OnePt *pts = fptr->pts;
  int *bincount = fptr->bincount;
  int *binfirst = fptr->binfirst;

  double epssq = fptr->epssq;
  int *nmatch1 = &fptr->nmatch1;
  int *nmatch2 = &fptr->nmatch2;
  int *errormatch1 = &fptr->errormatch1;
  int *errormatch2 = &fptr->errormatch2;

  tagint *tag = fptr->atom->tag;
  int vecflag = fptr->vecflag;

  // loop over each entry in buf
  // buf = corner pts for each tri owned by another proc (including self)
  // 0 = ID of tri, 1-9 = 3 corner pts

  m = 0;
  for (i = 0; i < n; i++) {
    tagother = (tagint) ubuf(buf[m]).i;
    cornersother = &buf[m+1];

    // loop over 3 pts in buf tri, owned by another proc (or me)
    // search for matching owned tri pt in ibin
    // same owned pt may appear many times in bins as part of multiple tris
    // opnext 1/2 are next points in order around tri perimeter

    for (optwhich = 1; optwhich <= 3; optwhich++) {
      if (optwhich == 1) {
        opt = &cornersother[0];
        optnext1 = &cornersother[3];
        optnext2 = &cornersother[6];
      } else if (optwhich == 2) {
        opt = &cornersother[3];
        optnext1 = &cornersother[6];
        optnext2 = &cornersother[0];
      } else {
        opt = &cornersother[6];
        optnext1 = &cornersother[0];
        optnext2 = &cornersother[3];
      }

      ibin = fptr->pt2bin3d(opt);

      if (ibin >= 0) {
        jfirst = binfirst[ibin];
        jlast = jfirst + bincount[ibin];

        // loop over all tri corner pts in bins[jfirst:jlast], owned by me
        // skip bin pt if from same tri ID as buf pt is from
        
        for (j = jfirst; j < jlast; j++) {
          iatom = pts[j].iatom;
          iconnect = pts[j].iconnect;
          ptwhich = pts[j].ptwhich;

          if (tagother == tag[iatom]) continue;

          if (ptwhich == 1) {
            pt = &corners[iconnect][0];
            ptnext1 = &corners[iconnect][3];
            ptnext2 = &corners[iconnect][6];
          } else if (ptwhich == 2) {
            pt = &corners[iconnect][3];
            ptnext1 = &corners[iconnect][6];
            ptnext2 = &corners[iconnect][0];
          } else {
            pt = &corners[iconnect][6];
            ptnext1 = &corners[iconnect][0];
            ptnext2 = &corners[iconnect][3];
          }

          dx = pt[0] - opt[0];
          dy = pt[1] - opt[1];
          dz = pt[2] - opt[2];
          rsq = dx*dx + dy*dy + dz*dz;

          // two points are close enough to match
          // use ptnext 1/2 and optnext 1/2 to determine if edge vs corner match
          // NOTE: this is where to set coupling flags

          if (rsq < epssq) {
            dx = ptnext1[0] - optnext1[0];
            dy = ptnext1[1] - optnext1[1];
            dz = ptnext1[2] - optnext1[2];
            rsq11 = dx*dx + dy*dy + dz*dz;
            dx = ptnext1[0] - optnext2[0];
            dy = ptnext1[1] - optnext2[1];
            dz = ptnext1[2] - optnext2[2];
            rsq12 = dx*dx + dy*dy + dz*dz;
            dx = ptnext2[0] - optnext1[0];
            dy = ptnext2[1] - optnext1[1];
            dz = ptnext2[2] - optnext1[2];
            rsq21 = dx*dx + dy*dy + dz*dz;
            dx = ptnext2[0] - optnext2[0];
            dy = ptnext2[1] - optnext2[1];
            dz = ptnext2[2] - optnext2[2];
            rsq22 = dx*dx + dy*dy + dz*dz;

            // edge match in same direction, should never happen

            if (rsq11 < epssq || rsq22 < epssq) {
              (*errormatch2)++;

            // edge match at first endpt, error if already marked
            // also add tagother to two corner pt lists

            } else if (rsq12 < epssq) {
              (*nmatch1)++;
              if (ptwhich == 1) {
                if (connect3d[iconnect].neigh_e1 > 0) (*errormatch1)++;
                else {
                  connect3d[iconnect].neigh_e1 = tagother;
                  if (vecflag) {
                    connect3d[iconnect].neigh_c1[connect3d[iconnect].nc1] = 
                      tagother;
                    connect3d[iconnect].neigh_c2[connect3d[iconnect].nc2] = 
                      tagother;
                  }
                  connect3d[iconnect].nc1++;
                  connect3d[iconnect].nc2++;
                }
              } else if (ptwhich == 2) {
                if (connect3d[iconnect].neigh_e2 > 0) (*errormatch1)++;
                else {
                  connect3d[iconnect].neigh_e2 = tagother;
                  if (vecflag) {
                    connect3d[iconnect].neigh_c2[connect3d[iconnect].nc2] = 
                      tagother;
                    connect3d[iconnect].neigh_c3[connect3d[iconnect].nc3] = 
                      tagother;
                  }
                  connect3d[iconnect].nc2++;
                  connect3d[iconnect].nc3++;
                }
              } else if (ptwhich == 3) {
                if (connect3d[iconnect].neigh_e3 > 0) (*errormatch1)++;
                else {
                  connect3d[iconnect].neigh_e3 = tagother;
                  if (vecflag) {
                    connect3d[iconnect].neigh_c3[connect3d[iconnect].nc3] = 
                      tagother;
                    connect3d[iconnect].neigh_c1[connect3d[iconnect].nc1] = 
                      tagother;
                  }
                  connect3d[iconnect].nc3++;
                  connect3d[iconnect].nc1++;
                }
              }

            // edge match at second endpt, other endpt will mark it

            } else if (rsq21 < epssq) {
              continue;

            // match at corner pt

            } else {
              (*nmatch2)++;
              if (ptwhich == 1) {
                if (vecflag) {
                  connect3d[iconnect].neigh_c1[connect3d[iconnect].nc1] = 
                    tagother;
                }
                connect3d[iconnect].nc1++;
              } else if (ptwhich == 2) {
                if (vecflag) {
                  connect3d[iconnect].neigh_c2[connect3d[iconnect].nc2] = 
                    tagother;
                }
                connect3d[iconnect].nc2++;
              } else if (ptwhich == 3) {
                if (vecflag) {
                  connect3d[iconnect].neigh_c3[connect3d[iconnect].nc3] = 
                    tagother;
                }
                connect3d[iconnect].nc3++;
              }
            }
          }
        }
      }
    }

    m += 10;
  }
}

/* ----------------------------------------------------------------------
   compute current corner points of N triangles
   nothing computed for particles that are not tris
------------------------------------------------------------------------- */

void FixSurfaceLocal::calculate_corners(int n)
{
  int ibonus;
  double p[3][3];
  double *corner;

  AtomVecTri::Bonus *bonus = avec_tri->bonus;
  double **x = atom->x;
  int *tri = atom->tri;

  for (int i = 0; i < n; i++) {
    if (tri[i] < 0) continue;
    ibonus = tri[i];
    corner = corners[i];
    MathExtra::quat_to_mat(bonus[ibonus].quat,p);
    MathExtra::matvec(p,bonus[ibonus].c1,&corner[0]);
    MathExtra::add3(x[i],&corner[0],&corner[0]);
    MathExtra::matvec(p,bonus[ibonus].c2,&corner[3]);
    MathExtra::add3(x[i],&corner[3],&corner[3]);
    MathExtra::matvec(p,bonus[ibonus].c3,&corner[6]);
    MathExtra::add3(x[i],&corner[6],&corner[6]);
  }
}

/* ----------------------------------------------------------------------
   map point pt to a bin
   return ibin from 0 to Nbins-1
   return -1 if outside bin bounds
------------------------------------------------------------------------- */

int FixSurfaceLocal::pt2bin3d(double *pt)
{
  if (pt[0] < binlo[0] || pt[0] >= binhi[0] || 
      pt[1] < binlo[1] || pt[1] >= binhi[1] ||
      pt[2] < binlo[2] || pt[2] >= binhi[2])
    return -1;

  int ix = static_cast<int> ((pt[0]-binlo[0]) * invbinx);
  ix = MIN(ix,nbinx-1);
  int iy = static_cast<int> ((pt[1]-binlo[1]) * invbiny);
  iy = MIN(iy,nbiny-1);
  int iz = static_cast<int> ((pt[2]-binlo[2]) * invbinz);
  iz = MIN(iz,nbinz-1);

  int ibin = iz*nbiny*nbinx + iy*nbinx + ix;
  return ibin;
}

/* ----------------------------------------------------------------------
   map point +/ EPS in both dims to all bins it overlaps with
   return N = # of overlapped bins
   return indices = list of bin IDs, each from 0 to Nbins-1
------------------------------------------------------------------------- */

int FixSurfaceLocal::overlap2bin3d(double *pt, double eps, int *indices)
{
  int ilo = static_cast<int> ((pt[0]-eps-binlo[0]) * invbinx);
  int ihi = static_cast<int> ((pt[0]+eps-binlo[0]) * invbinx);
  int jlo = static_cast<int> ((pt[1]-eps-binlo[1]) * invbiny);
  int jhi = static_cast<int> ((pt[1]+eps-binlo[1]) * invbiny);
  int klo = static_cast<int> ((pt[2]-eps-binlo[2]) * invbinz);
  int khi = static_cast<int> ((pt[2]+eps-binlo[2]) * invbinz);
  
  ilo = MAX(ilo,0);
  ihi = MIN(ihi,nbinx-1);
  jlo = MAX(jlo,0);
  jhi = MIN(jhi,nbiny-1);
  klo = MAX(klo,0);
  khi = MIN(khi,nbinz-1);

  if (ilo == ihi && jlo == jhi && klo == khi) {
    indices[0] = klo*nbiny*nbinx + jlo*nbinx + ilo;
    return 1;
  }

  int i,j,k;
  int n = 0;
  for (k = klo; k <= khi; k++) 
    for (j = jlo; j <= jhi; j++) 
      for (i = ilo; i <= ihi; i++)
	indices[n++] = k*nbiny*nbinx + j*nbinx + i;
  return n;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for global surf connectivity build, assignment to procs
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   extract points,line,surfs from molecule ID for one or more mol files
   concatenate into single list of points,lines,tris
------------------------------------------------------------------------- */

void FixSurfaceLocal::extract_from_molecules(char *molID)
{
  // check that no line/tri particles already exist
  // no connectivity would be produced for them

  int *surf = atom->line;
  if (dimension == 3) surf = atom->tri;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (surf[i] >= 0) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall) 
    error->all(FLERR,"Fix surface local molecule file "
               "when surfaces already exist");

  // populate global point/line/tri data structs

  points = NULL;
  lines = NULL;
  tris = NULL;
  npoints = nlines = ntris = 0;

  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for "
               "fix surface/local does not exist");
  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;
  for (int m = 0; m < nmol; m++) {
    if (onemols[m]->pointflag == 0)
      error->all(FLERR,"Fix surface/local molecule must have points");
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface/local molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface/local molecule must have triangles");

    // NOTE: anything else about molfile surf to check?
    //       e.g. are types within bounds? 
    //       or did Molecule check at read?

    // NOTE: for nmol = 1, could just set points,lines,tris to pt
    //       to chunk of 2d data in molfile arrays

    int np = onemols[m]->npoints;
    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    npoints += np;
    nlines += nl;
    ntris += nt;
    points = (Point *) memory->srealloc(points,npoints*sizeof(Point),
                                        "surface/local:points");
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface/local:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface/local:tris");

    double **pts = onemols[m]->points;
    int j = npoints - np;
    for (int i = 0; i < np; i++) {
      points[j].x[0] = pts[i][0];
      points[j].x[1] = pts[i][1];
      points[j].x[2] = pts[i][2];
      j++;
    }

    // need to offset line/tri index lists by previous npoints & subtract one

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      int **epts = onemols[m]->lines;
      int j = nlines - nl;
      for (int i = 0; i < nl; i++) {
        lines[j].mol = molline[i];
        lines[j].type = typeline[i];
        lines[j].p1 = epts[i][0] + npoints-np - 1;
        lines[j].p2 = epts[i][1] + npoints-np - 1;
        j++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      int **cpts = onemols[m]->tris;
      int j = ntris - nt;
      for (int i = 0; i < nt; i++) {
        tris[j].mol = moltri[i];
        tris[j].type = typetri[i];
        tris[j].p1 = cpts[i][0] + npoints-np - 1;
        tris[j].p2 = cpts[i][1] + npoints-np - 1;
        tris[j].p3 = cpts[i][2] + npoints-np - 1;
        j++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for all lines
   stored in connect2dall
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity2d_global()
{
  int j,p1,p2;

  // allocate and initilize global connectivity list
  // ilocal will be initialized when assign to procs

  connect2dall = (Connect2d *) memory->smalloc(nlines*sizeof(Connect2d),
                                               "surface/local:connect2dall");

  for (int i = 0; i < nlines; i++) {
    connect2dall[i].neigh_p1 = connect2dall[i].neigh_p2 = 0;
    connect2dall[i].flags = connect2dall[i].ilocal = 0;
  }

  // create ptflag = pair of flags for each point
  //   1st index:
  //     0 if pt not in line
  //     1 if 1st point in one line, 2 if 2nd point in one line
  //    -1 if already in 2 lines
  //   2nd index: which line it is in (for one line)
  // as create it, check for errors:
  //   no point can be part of more than 2 lines
  //   point that is part of 2 lines must be 1st in one, 2nd in other
  // as create it, set connections in connect2dall
  //   set neigh_p1,neigh_p2 = 1 to Nlines for now
  //   will offset by existing particle IDs when assign to procs

  int **ptflag;
  memory->create(ptflag,npoints,2,"surface/local:ptflag");
  for (int i = 0; i < npoints; i++) ptflag[i][0] = 0;

  for (int i = 0; i < nlines; i++) {
    p1 = lines[i].p1;
    p2 = lines[i].p2;

    if (ptflag[p1][0] < 0) 
      error->all(FLERR,"Fix surface/local point part of more than 2 lines");
    else if (ptflag[p1][0] == 1)
      error->all(FLERR,"Fix surface/local pair of lines are misoriented");
    else if (ptflag[p1][0] == 0) {
      ptflag[p1][0] = 1;
      ptflag[p1][1] = i;
    } else if (ptflag[p1][0] == 2) {
      ptflag[p1][0] = -1;
      j = ptflag[p1][1];
      connect2dall[i].neigh_p1 = j+1;
      connect2dall[j].neigh_p2 = i+1;
    }

    if (ptflag[p2][0] < 0) 
      error->all(FLERR,"Fix surface/local point part of more than 2 lines");
    else if (ptflag[p2][0] == 2)
      error->all(FLERR,"Fix surface/local pair of lines are misoriented");
    else if (ptflag[p2][0] == 0) {
      ptflag[p2][0] = 2;
      ptflag[p2][1] = i;
    } else if (ptflag[p2][0] == 1) {
      ptflag[p2][0] = -1;
      j = ptflag[p2][1];
      connect2dall[i].neigh_p2 = j+1;
      connect2dall[j].neigh_p1 = i+1;
    }
  }

  memory->destroy(ptflag);

  // NOTE: here is where to set connect2d flags
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all triangles
   stored in connect3dall
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity3d_global()
{
  int itri,jtri,iedge,jedge;
  bigint p1,p2,p3;

  // allocate and initilize global connectivity list
  // ilocal and indexc123 will be initialized when assign to procs

  connect3dall = (Connect3d *) memory->smalloc(ntris*sizeof(Connect3d),
                                               "surface/local:connect3dall");

  for (int i = 0; i < ntris; i++) {
    connect3dall[i].neigh_e1 = connect3dall[i].neigh_e2 = 
      connect3dall[i].neigh_e3 = 0;
    connect3dall[i].flags = 0;
  }

  // hash = STL map of ordered edges
  // key = (p1,p2) via bit-shifting by 32-bits into bigint
  // value = (itri,iedge) also bit-shifted, itri = 0 to Ntri-1, iedge = 1,2,3

  if (sizeof(bigint) != 2*sizeof(int))
    error->all(FLERR,"Fix surface/local triangle connections cannot be formed");

  bigint key,value;
  std::map<bigint,bigint> hash;
  std::map<bigint,bigint>::iterator it;

  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    key = (p1 << 32) | p2;
    value = (((bigint) i) << 32) | 1;
    if (hash.find(key) != hash.end())
      error->all(FLERR,"Fix surface/local edge is part of more than 2 tris "
                 "or its triangle is misoriented");
    else hash[key] = value;

    key = (p2 << 32) | p3;
    value = (((bigint) i) << 32) | 2;
    if (hash.find(key) != hash.end())
      error->all(FLERR,"Fix surface/local edge is part of more than 2 tris "
                 "or its triangle is misoriented");
    else hash[key] = value;

    key = (p3 << 32) | p1;
    value = (((bigint) i) << 32) | 3;
    if (hash.find(key) != hash.end())
      error->all(FLERR,"Fix surface/local edge is part of more than 2 tris "
                 "or its triangle is misoriented");
    else hash[key] = value;
  }

  // set edge connections via hash
  // test itri < jtri to avoid resetting (with identical values)

  for (it = hash.begin(); it != hash.end(); ++it) {
    p1 = it->first >> 32;
    p2 = it->first & MAXSMALLINT;
    key = (p2 << 32) | p1;
    if (hash.find(key) != hash.end()) {
      itri = it->second >> 32;
      iedge = it->second & MAXSMALLINT;
      value = hash[key];
      jtri = value >> 32;
      jedge = value & MAXSMALLINT;
      if (itri < jtri) {
        if (iedge == 1) connect3dall[itri].neigh_e1 = jtri+1;
        else if (iedge == 2) connect3dall[itri].neigh_e2 = jtri+1;
        else if (iedge == 3) connect3dall[itri].neigh_e3 = jtri+1;
        if (jedge == 1) connect3dall[jtri].neigh_e1 = itri+1;
        else if (jedge == 2) connect3dall[jtri].neigh_e2 = itri+1;
        else if (jedge == 3) connect3dall[jtri].neigh_e3 = itri+1;
      }
    }
  }

  // setup corner point connectivity lists
  // count # of tris containing each point
  // create ragged 2d array to contain all tri indices, then fill it
  // set neigh_c123 vector ptrs in connect3dall to rows of ragged array
  //   nc123 counts are set correctly (exclude self tri)
  //   neigh_c123 vectors are one too long (include self, unless count = 0)

  int *counts;
  memory->create(counts,npoints,"surface/local:count");
  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tris[i].p1]++;
    counts[tris[i].p2]++;
    counts[tris[i].p3]++;
  }

  memory->create_ragged(clist,npoints,counts,"surface/local:clist");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    clist[tris[i].p1][counts[tris[i].p1]++] = i+1;
    clist[tris[i].p2][counts[tris[i].p2]++] = i+1;
    clist[tris[i].p3][counts[tris[i].p3]++] = i+1;
  }

  for (int i = 0; i < ntris; i++)
    connect3dall[i].nc1 = connect3dall[i].nc2 = connect3dall[i].nc3 = 0;

  for (int i = 0; i < ntris; i++) {
    connect3dall[i].nc1 = counts[tris[i].p1] - 1;
    if (connect3dall[i].nc1 == 0) connect3dall[i].neigh_c1 = NULL;
    else connect3dall[i].neigh_c1 = clist[tris[i].p1];
    connect3dall[i].nc2 = counts[tris[i].p2] - 1;
    if (connect3dall[i].nc2 == 0) connect3dall[i].neigh_c2 = NULL;
    else connect3dall[i].neigh_c2 = clist[tris[i].p2];
    connect3dall[i].nc3 = counts[tris[i].p3] - 1;
    if (connect3dall[i].nc3 == 0) connect3dall[i].neigh_c3 = NULL;
    else connect3dall[i].neigh_c3 = clist[tris[i].p3];
  }

  memory->destroy(counts);

  // NOTE: here is where to set connect3d flags
}

/* ----------------------------------------------------------------------
   assign lines and their connectivity from global data structs to each proc
   logic for proc assign follows Atom::data_atoms(), as if read from data file
------------------------------------------------------------------------- */

void FixSurfaceLocal::assign2d()
{
  int i;

  AtomVec *avec = atom->avec;
  avec_line = (AtomVecLine *) atom->style_match("line");
  if (!avec_line) 
    error->all(FLERR,"Fix surface/local requires atom style line");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }

  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // idmaxall = largest existing atom ID

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (i = 0; i < nlocal; i++) idmax = MAX(tag[i],idmax);
  tagint idmaxall;
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // loop over all global lines
  // compute their midpoints
  // keep lines belonging to this proc

  imageint imagedata;
  double xmid[3],lamda[3];
  double *coord,*x1,*x2;

  char pts[4][32];
  char *values[4];
  for (i = 0; i < 4; i++) values[i] = &pts[i][0];
  std::vector<std::string> svalues(4);
  
  for (i = 0; i < nlines; i++) {

    // midpt of line

    x1 = points[lines[i].p1].x;
    x2 = points[lines[i].p2].x;
    xmid[0] = 0.5*(x1[0]+x2[0]);
    xmid[1] = 0.5*(x1[1]+x2[1]);
    xmid[2] = 0.5*(x1[2]+x2[2]);

    imagedata = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    // NOTE: this will wrap line center back into periodic box

    domain->remap(xmid,imagedata);
    if (triclinic) {
      domain->x2lamda(xmid,lamda);
      coord = lamda;
    } else coord = xmid;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {

      // create a default particle of correct style

      avec->create_atom(lines[i].type,xmid);
      
      // change it to be a line

      int n = atom->nlocal - 1;
      atom->tag[n] = idmaxall + i + 1;
      atom->image[n] = imagedata;
      atom->molecule[n] = lines[i].mol;
      atom->line[n] = 0;
      atom->rmass[n] = 1.0;              // set line density = 1.0
      atom->radius[n] = 0.0;

      // set the bonus data for a line = end points

      sprintf(values[0],"%-1.16e",x1[0]);
      sprintf(values[1],"%-1.16e",x1[1]);
      sprintf(values[2],"%-1.16e",x2[0]);
      sprintf(values[3],"%-1.16e",x2[1]);

      svalues[0] = (const char *) values[0];
      svalues[1] = (const char *) values[1];
      svalues[2] = (const char *) values[2];
      svalues[3] = (const char *) values[3];
      
      avec_line->data_atom_bonus(n,svalues);

      // copy global connectivity enty for line I to local list
      // reset neigh p1/p2 to new particle IDs beyond idmaxall

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect2d[nlocal_connect],&connect2dall[i],sizeof(Connect2d));
      connect2d[nlocal_connect].neigh_p1 += idmaxall;
      connect2d[nlocal_connect].neigh_p2 += idmaxall;
      connect2d[nlocal_connect].neigh_p2 += idmaxall;
      connect2d[nlocal_connect].ilocal = n;
      index[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect2dall);

  // create global mapping that includes added lines

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}

/* ----------------------------------------------------------------------
   assign tris and their connectivity from global data structs to each proc
------------------------------------------------------------------------- */

void FixSurfaceLocal::assign3d()
{
  int i,j,m;

  AtomVec *avec = atom->avec;
  avec_tri = (AtomVecTri *) atom->style_match("tri");
  if (!avec_tri) 
    error->all(FLERR,"Fix surface/local requires atom style tri");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
    }

  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
  }

  // idmaxall = largest existing atom ID

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (i = 0; i < nlocal; i++) idmax = MAX(tag[i],idmax);
  tagint idmaxall;
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // loop over all global lines
  // compute their midpoints
  // keep lines belonging to this proc

  int nc;
  imageint imagedata;
  tagint idself;
  double xmid[3],lamda[3];
  double *coord,*x1,*x2,*x3;
  tagint *self,*neigh;

  char pts[9][32];
  char *values[9];
  for (i = 0; i < 9; i++) values[i] = &pts[i][0];
  std::vector<std::string> svalues(9);

  for (i = 0; i < ntris; i++) {

    // center pt of triangle

    x1 = points[tris[i].p1].x;
    x2 = points[tris[i].p2].x;
    x3 = points[tris[i].p3].x;
    xmid[0] = (x1[0]+x2[0]+x3[0]) / 3.0;
    xmid[1] = (x1[1]+x2[1]+x3[1]) / 3.0;
    xmid[2] = (x1[2]+x2[2]+x3[2]) / 3.0;

    imagedata = ((imageint) IMGMAX << IMG2BITS) |
      ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    // NOTE: this will wrap triangle center back into periodic box

    domain->remap(xmid,imagedata);
    if (triclinic) {
      domain->x2lamda(xmid,lamda);
      coord = lamda;
    } else coord = xmid;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {

      // create a default particle of correct style

      avec->create_atom(tris[i].type,xmid);
      
      // change it to be a triangle

      int n = atom->nlocal - 1;
      atom->tag[n] = idmaxall + i + 1;
      atom->image[n] = imagedata;
      atom->molecule[n] = tris[i].mol;
      atom->tri[n] = 0;
      atom->rmass[n] = 1.0;              // set tri density = 1.0
      atom->radius[n] = 0.0;

      // set the bonus data for a tri = corner pts

      sprintf(values[0],"%-1.16e",x1[0]);
      sprintf(values[1],"%-1.16e",x1[1]);
      sprintf(values[2],"%-1.16e",x1[2]);
      sprintf(values[3],"%-1.16e",x2[0]);
      sprintf(values[4],"%-1.16e",x2[1]);
      sprintf(values[5],"%-1.16e",x2[2]);
      sprintf(values[6],"%-1.16e",x3[0]);
      sprintf(values[7],"%-1.16e",x3[1]);
      sprintf(values[8],"%-1.16e",x3[2]);

      svalues[0] = (const char *) values[0];
      svalues[1] = (const char *) values[1];
      svalues[2] = (const char *) values[2];
      svalues[3] = (const char *) values[3];
      svalues[4] = (const char *) values[4];
      svalues[5] = (const char *) values[5];
      svalues[6] = (const char *) values[6];
      svalues[7] = (const char *) values[7];
      svalues[8] = (const char *) values[8];

      avec_tri->data_atom_bonus(n,svalues);

      // copy global connectivity enty for tri I to local list
      // reset neigh e1/e2/e3 to new particle IDs beyond idmaxall
      // replace neigh c1/c2/c3 ptrs into clist with tcp pool allocs
      // nc1,nc2,nc3 are correct reduced lengths, but clist ptrs have
      //   one extra self value that must be removed

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect3d[nlocal_connect],&connect3dall[i],sizeof(Connect3d));
      connect3d[nlocal_connect].neigh_e1 += idmaxall;
      connect3d[nlocal_connect].neigh_e2 += idmaxall;
      connect3d[nlocal_connect].neigh_e3 += idmaxall;

      idself = atom->tag[n];

      nc = connect3d[nlocal_connect].nc1;
      self = connect3d[nlocal_connect].neigh_c1;
      for (j = 0; j < nc; j++) self[j] += idmaxall;
      if (nc) {
        neigh = connect3d[nlocal_connect].neigh_c1 = 
          tcp->get(nc,connect3d[nlocal_connect].indexc1);
        m = 0;
        for (j = 0; j <= nc; j++)
          if (self[j] != idself) neigh[m++] = self[j];
      } else connect3d[nlocal_connect].neigh_c1 = NULL;

      nc = connect3d[nlocal_connect].nc2;
      self = connect3d[nlocal_connect].neigh_c2;
      for (j = 0; j < nc; j++) self[j] += idmaxall;
      if (nc) {
        neigh = connect3d[nlocal_connect].neigh_c2 = 
          tcp->get(nc,connect3d[nlocal_connect].indexc2);
        m = 0;
        for (j = 0; j <= nc; j++)
          if (self[j] != idself) neigh[m++] = self[j];
      } else connect3d[nlocal_connect].neigh_c2 = NULL;

      nc = connect3d[nlocal_connect].nc3;
      self = connect3d[nlocal_connect].neigh_c3;
      for (j = 0; j < nc; j++) self[j] += idmaxall;
      if (nc) {
        neigh = connect3d[nlocal_connect].neigh_c3 = 
          tcp->get(nc,connect3d[nlocal_connect].indexc3);
        m = 0;
        for (j = 0; j <= nc; j++)
          if (self[j] != idself) neigh[m++] = self[j];
      } else connect3d[nlocal_connect].neigh_c3 = NULL;

      connect3d[nlocal_connect].ilocal = n;
      index[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect3dall);
  memory->destroy(clist);

  // create global mapping that includes added tris

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}
