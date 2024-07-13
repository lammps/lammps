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
#include "stl_reader.h"

#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;

// NOTE: epsilon should be defined as fraction of max surf size

#define EPSILON 0.001
#define NBIN 100
#define BIG 1.0e20
#define MAXLINE 256
#define MAXTRIPOINT 24

#define DELTA 128
#define DELTA_CONNECT 4            // NOTE: make it larger when done testing

enum{DATAFILE,MOLTEMPLATE,STLFILE};
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

  cindex = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(2);

  nlocal_connect = nghost_connect = nmax_connect = 0;
  connect2d = nullptr;
  connect3d = nullptr;

  // NOTE: what is optimal way to initialize MyPoolChunk?

  tcp = new MyPoolChunk<tagint>(1,MAXTRIPOINT,6);

  // process args

  sourceID = nullptr;
  
  if (strcmp(arg[3],"NULL") == 0) mode = DATAFILE;
  else {
    int n = strlen(arg[3]) + 1;
    sourceID = new char[n];
    strcpy(sourceID,arg[3]);
    int imol = atom->find_molecule(sourceID);
    if (imol >= 0) mode = MOLTEMPLATE;
    else mode = STLFILE;
  }
}

/* ---------------------------------------------------------------------- */

FixSurfaceLocal::~FixSurfaceLocal()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete [] sourceID;
  memory->destroy(cindex);

  // return all connection vector memory to TCP
  
  if (dimension == 2) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = 0; i < nall; i++) {
      if (connect2d[i].neigh_p1) tcp->put(connect2d[i].indexp1);
      if (connect2d[i].neigh_p2) tcp->put(connect2d[i].indexp2);
    }
    memory->sfree(connect2d);
  } else {
    int nall = nlocal_connect + nghost_connect;
    for (int i = 0; i < nall; i++) {
      if (connect3d[i].neigh_e1) tcp->put(connect3d[i].indexe1);
      if (connect3d[i].neigh_e2) tcp->put(connect3d[i].indexe2);
      if (connect3d[i].neigh_e3) tcp->put(connect3d[i].indexe3);
      if (connect3d[i].neigh_c1) tcp->put(connect3d[i].indexc1);
      if (connect3d[i].neigh_c2) tcp->put(connect3d[i].indexc2);
      if (connect3d[i].neigh_c3) tcp->put(connect3d[i].indexc3);
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
   define distributed lines/tris and their connectivity in one of 3 ways
   (1) mode = DATAFILE, lines/tris were already read from data file
   (2) mode = MOLTEMPLATE, lines/tris are in molecule template ID
   {3} mode = STLFILE, tris are in an STL file
   must be done in post_constructor() b/c MOLTEMPLATE and STLFILE
     add owned lines/tris to AtomVec class
     its grow() method makes a callback to grow_arrays() in this fix
     callback can't be invoked unless fix is fully instantiated
------------------------------------------------------------------------- */

void FixSurfaceLocal::post_constructor()
{
  // for DATAFILE:
  // lines or tri are already distributed across procs
  // task is to infer line/tri connectivity in distributed manner

  if (mode == DATAFILE) {
    if (comm->me == 0 && screen)
      fprintf(screen,"Connecting line/tri particles ...\n");

    if (dimension == 2) connectivity2d_local();
    else connectivity3d_local();
  }

  // for MOLTEMPLATE or STLFILE
  // error check that no line/tri particles already exist
  //   b/c no connectivity would be produced for them
  // read in lines/tris from appropriate sourrce
  // each proc builds global data structs of points/lines/tris
  // each proc infers global connectivity from global data structs
  // distribute lines/surfs across procs, based on center pt coords
  // delete global data structs
  
  if (mode == MOLTEMPLATE || mode == STLFILE) {
    if (check_exist())
      error->all(FLERR,"Fix surface/local with non-NULL source when lines/tris already exist");
    
    if (comm->me == 0 && screen) {
      if (mode == MOLTEMPLATE)
        fprintf(screen,"Converting molecule file to line/tri particles ...\n");
      if (mode == STLFILE)
        fprintf(screen,"Reading STL file for triangle particles ...\n");
    }

    points = nullptr;
    lines = nullptr;
    tris = nullptr;

    plist = nullptr;
    elist = nullptr;
    clist = nullptr;
    
    if (mode == MOLTEMPLATE) extract_from_molecules(sourceID);
    if (mode == STLFILE) extract_from_stlfile(sourceID);

    if (dimension == 2) {
      connectivity2d_global();
      assign2d();
    } else {
      connectivity3d_global();
      assign3d();
    }

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

    memory->sfree(points);
    memory->sfree(lines);
    memory->sfree(tris);
    
    memory->destroy(plist);
    memory->destroy(elist);
    memory->destroy(clist);
  }
  
  // set max size of connection info
  // NOTE: set 2d/3d based on # max of endpt/corner connections (not 12) ??
  
  if (dimension == 2) comm_border = 4 + 2*12;
  else comm_border = 8 + 6*12;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_arrays(int nmax)
{
  memory->grow(cindex,nmax,"surface/local:cindex");
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

  for (int i = 0; i < n; i++) {
    if (surf[i] < 0) continue;

    if (cindex[i] != surf[i]) {
      count1++;
      continue;
    }
    if (dimension == 2 && connect2d[cindex[i]].ilocal != i) count2++;
    if (dimension == 3 && connect3d[cindex[i]].ilocal != i) count2++;
  }

  int all1,all2;
  MPI_Allreduce(&count1,&all1,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&count2,&all2,1,MPI_INT,MPI_SUM,world);

  if (all1 || all2 && comm->me == 0) {
    char str[128];
    sprintf(str,"FSL cindex vector mis-match: %d %d\n",all1,all2);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   grow connect array by DELTA_CONNECT
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_connect()
{
  nmax_connect = nmax_connect/DELTA_CONNECT * DELTA_CONNECT;
  bigint newmax = (bigint) nmax_connect + DELTA_CONNECT;
  if (newmax > MAXSMALLINT) error->one(FLERR,"Too many fix surface/local connections");
  nmax_connect = newmax;

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
   copy atom I to atom J with optional delflag
   I or J or both can be a line/triangle with connection info
------------------------------------------------------------------------- */

void FixSurfaceLocal::copy_arrays(int i, int j, int delflag)
{
  // if overwriting atom J via delflag and J is a line/tri with connection data:
  //   return J's vector memory to TCP before memcpy()
  //   shrink list of owned connections by copying last element to Kth location

  if (dimension == 2) {
    if (delflag && cindex[j] >= 0) {
      int k = cindex[j];
      if (connect2d[k].neigh_p1) tcp->put(connect2d[k].indexp1);
      if (connect2d[k].neigh_p2) tcp->put(connect2d[k].indexp2);
      cindex[connect2d[nlocal_connect-1].ilocal] = k;
      memcpy(&connect2d[k],&connect2d[nlocal_connect-1],sizeof(Connect2d));
      nlocal_connect--;
    }
  } else {
    if (delflag && cindex[j] >= 0) {
      int k = cindex[j];
      if (connect3d[k].neigh_e1) tcp->put(connect3d[k].indexe1);
      if (connect3d[k].neigh_e2) tcp->put(connect3d[k].indexe2);
      if (connect3d[k].neigh_e3) tcp->put(connect3d[k].indexe3);
      if (connect3d[k].neigh_c1) tcp->put(connect3d[k].indexc1);
      if (connect3d[k].neigh_c2) tcp->put(connect3d[k].indexc2);
      if (connect3d[k].neigh_c3) tcp->put(connect3d[k].indexc3);
      cindex[connect3d[nlocal_connect-1].ilocal] = k;
      memcpy(&connect3d[k],&connect3d[nlocal_connect-1],sizeof(Connect3d));
      nlocal_connect--;
    }
  }

  // if atom I has connection data, reset I's connect.ilocal to loc J
  // do NOT do this if self-copy (I=J) since I's connection data is already deleted
  
  if (cindex[i] >= 0 && i != j) {
    if (dimension == 2) connect2d[cindex[i]].ilocal = j;
    else connect3d[cindex[i]].ilocal = j;
  }
  cindex[j] = cindex[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixSurfaceLocal::set_arrays(int i)
{
  cindex[i] = -1;
}

/* ----------------------------------------------------------------------
   clear ghost info in connect array
   called before ghosts are recommunicated in comm and irregular
   return all vector memory to TCP
------------------------------------------------------------------------- */

void FixSurfaceLocal::clear_bonus()
{
  if (dimension == 2) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = nlocal_connect; i < nall; i++) {
      if (connect2d[i].neigh_p1) tcp->put(connect2d[i].indexp1);
      if (connect2d[i].neigh_p2) tcp->put(connect2d[i].indexp2);
    }
  } else if (dimension == 3) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = nlocal_connect; i < nall; i++) {
      if (connect3d[i].neigh_e1) tcp->put(connect3d[i].indexe1);
      if (connect3d[i].neigh_e2) tcp->put(connect3d[i].indexe2);
      if (connect3d[i].neigh_e3) tcp->put(connect3d[i].indexe3);
      if (connect3d[i].neigh_c1) tcp->put(connect3d[i].indexc1);
      if (connect3d[i].neigh_c2) tcp->put(connect3d[i].indexc2);
      if (connect3d[i].neigh_c3) tcp->put(connect3d[i].indexc3);
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
    int np1,np2;

    for (i = 0; i < n; i++) {
      j = list[i];
      if (cindex[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        ic = cindex[j];

        np1 = connect2d[ic].np1;
        np2 = connect2d[ic].np2;
        buf[m++] = ubuf(np1).d;
        buf[m++] = ubuf(np2).d;
        if (np1 > 1)
          for (k = 0; k < np1; k++)
            buf[m++] = ubuf(connect2d[ic].neigh_p1[k]).d;
        if (np2 > 1)
          for (k = 0; k < np2; k++)
            buf[m++] = ubuf(connect2d[ic].neigh_p2[k]).d;
        
        buf[m++] = ubuf(connect2d[ic].flags).d;
      }
    }

  } else {
    int ne1,ne2,ne3,nc1,nc2,nc3;

    for (i = 0; i < n; i++) {
      j = list[i];
      if (cindex[j] < 0) buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        ic = cindex[j]; 

        ne1 = connect3d[ic].ne1;
        ne2 = connect3d[ic].ne2;
        ne3 = connect3d[ic].ne3;
        buf[m++] = ubuf(ne1).d;
        buf[m++] = ubuf(ne2).d;
        buf[m++] = ubuf(ne3).d;
        if (ne1 > 1)
          for (k = 0; k < ne1; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_e1[k]).d;
        if (ne2 > 1)
          for (k = 0; k < ne2; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_e2[k]).d;
        if (ne3 > 1)
          for (k = 0; k < ne3; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_e3[k]).d;

        nc1 = connect3d[ic].nc1;
        nc2 = connect3d[ic].nc2;
        nc3 = connect3d[ic].nc3;
        buf[m++] = ubuf(nc1).d;
        buf[m++] = ubuf(nc2).d;
        buf[m++] = ubuf(nc3).d;
        if (nc1 > 1)
          for (k = 0; k < nc1; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_c1[k]).d;
        if (nc2 > 1)
          for (k = 0; k < nc2; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_c2[k]).d;
        if (nc3 > 1)
          for (k = 0; k < nc3; k++)
            buf[m++] = ubuf(connect3d[ic].neigh_c3[k]).d;

        buf[m++] = ubuf(connect3d[ic].flags).d;
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
  int i,j,k,m,last,flag;

  m = 0;
  last = first + n;

  if (dimension == 2) {
    int np1,np2;
    
    for (i = first; i < last; i++) {
      flag = (int) ubuf(buf[m++]).i;
      if (flag == 0) cindex[i] = -1;
      else {
        j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();

        np1 = (int) ubuf(buf[m++]).i;
        np2 = (int) ubuf(buf[m++]).i;
        connect2d[j].np1 = np1;
        connect2d[j].np2 = np2;
        if (np1 > 1) {
          connect2d[j].neigh_p1 = tcp->get(np1,connect2d[j].indexp1);
          for (k = 0; k < np1; k++)
            connect2d[j].neigh_p1[k] = (tagint) ubuf(buf[m++]).i;
        } else connect2d[j].neigh_p1 = nullptr;
        if (np2 > 1) {
          connect2d[j].neigh_p2 = tcp->get(np2,connect2d[j].indexp2);
          for (k = 0; k < np2; k++)
            connect2d[j].neigh_p2[k] = (tagint) ubuf(buf[m++]).i;
        } else connect2d[j].neigh_p2 = nullptr;

        connect2d[j].flags = (int) ubuf(buf[m++]).i;

        connect2d[j].ilocal = i;
        cindex[i] = j;
        nghost_connect++;
      }
    }

  } else {
    int ne1,ne2,ne3,nc1,nc2,nc3;

    for (i = first; i < last; i++) {
      flag = (int) ubuf(buf[m++]).i;
      if (flag == 0) cindex[i] = -1;
      else {
        j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();

        ne1 = (int) ubuf(buf[m++]).i;
        ne2 = (int) ubuf(buf[m++]).i;
        ne3 = (int) ubuf(buf[m++]).i;
        connect3d[j].ne1 = ne1;
        connect3d[j].ne2 = ne2;
        connect3d[j].ne3 = ne3;
        if (ne1 > 1) {
          connect3d[j].neigh_e1 = tcp->get(ne1,connect3d[j].indexe1);
          for (k = 0; k < ne1; k++)
            connect3d[j].neigh_e1[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_e1 = nullptr;
        if (ne2 > 1) {
          connect3d[j].neigh_e2 = tcp->get(ne2,connect3d[j].indexe2);
          for (k = 0; k < ne2; k++)
            connect3d[j].neigh_e2[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_e2 = nullptr;
        if (ne3 > 1) {
          connect3d[j].neigh_e3 = tcp->get(ne3,connect3d[j].indexe3);
          for (k = 0; k < ne3; k++)
            connect3d[j].neigh_e3[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_e3 = nullptr;

        nc1 = (int) ubuf(buf[m++]).i;
        nc2 = (int) ubuf(buf[m++]).i;
        nc3 = (int) ubuf(buf[m++]).i;
        connect3d[j].nc1 = nc1;
        connect3d[j].nc2 = nc2;
        connect3d[j].nc3 = nc3;
        if (nc1 > 1) {
          connect3d[j].neigh_c1 = tcp->get(nc1,connect3d[j].indexc1);
          for (k = 0; k < nc1; k++)
            connect3d[j].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_c1 = nullptr;
        if (nc2 > 1) {
          connect3d[j].neigh_c2 = tcp->get(nc2,connect3d[j].indexc2);
          for (k = 0; k < nc2; k++)
            connect3d[j].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_c2 = nullptr;
        if (nc3 > 1) {
          connect3d[j].neigh_c3 = tcp->get(nc3,connect3d[j].indexc3);
          for (k = 0; k < nc3; k++)
            connect3d[j].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;
        } else connect3d[j].neigh_c3 = nullptr;
        
        connect3d[j].flags = (int) ubuf(buf[m++]).i;

        connect3d[j].ilocal = i;
        cindex[i] = j;
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
    int np1,np2;

    if (cindex[i] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = cindex[i];

      np1 = connect2d[j].np1;
      np2 = connect2d[j].np2;
      buf[m++] = ubuf(np1).d;
      buf[m++] = ubuf(np2).d;
      if (np1 > 1)
        for (k = 0; k < np1; k++)
          buf[m++] = ubuf(connect2d[j].neigh_p1[k]).d;
      if (np2 > 1)
        for (k = 0; k < np2; k++)
          buf[m++] = ubuf(connect2d[j].neigh_p2[k]).d;

      buf[m++] = ubuf(connect2d[j].flags).d;
    }

  } else {
    int ne1,ne2,ne3,nc1,nc2,nc3;

    if (cindex[i] < 0) buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = cindex[i];

      ne1 = connect3d[j].ne1;
      ne2 = connect3d[j].ne2;
      ne3 = connect3d[j].ne3;
      buf[m++] = ubuf(ne1).d;
      buf[m++] = ubuf(ne2).d;
      buf[m++] = ubuf(ne3).d;
      if (ne1 > 1)
        for (k = 0; k < ne1; k++)
          buf[m++] = ubuf(connect3d[j].neigh_e1[k]).d;
      if (ne2 > 1)
        for (k = 0; k < ne2; k++)
          buf[m++] = ubuf(connect3d[j].neigh_e2[k]).d;
      if (ne3 > 1)
        for (k = 0; k < ne3; k++)
          buf[m++] = ubuf(connect3d[j].neigh_e3[k]).d;

      nc1 = connect3d[j].nc1;
      nc2 = connect3d[j].nc2;
      nc3 = connect3d[j].nc3;
      buf[m++] = ubuf(nc1).d;
      buf[m++] = ubuf(nc2).d;
      buf[m++] = ubuf(nc3).d;
      if (nc1 > 1)
        for (k = 0; k < nc1; k++)
          buf[m++] = ubuf(connect3d[j].neigh_c1[k]).d;
      if (nc2 > 1)
        for (k = 0; k < nc2; k++)
          buf[m++] = ubuf(connect3d[j].neigh_c2[k]).d;
      if (nc3 > 1)
        for (k = 0; k < nc3; k++)
          buf[m++] = ubuf(connect3d[j].neigh_c3[k]).d;

      buf[m++] = ubuf(connect3d[j].flags).d;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local connect array from exchange with another proc
------------------------------------------------------------------------- */

int FixSurfaceLocal::unpack_exchange(int nlocal, double *buf)
{
  int j,k,n,flag;
  int m = 0;

  if (dimension == 2) {
    int np1,np2;

    flag = (int) ubuf(buf[m++]).i;
    if (flag == 0) cindex[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      np1 = (int) ubuf(buf[m++]).i;
      np2 = (int) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].np1 = np1;
      connect2d[nlocal_connect].np2 = np2;
      if (np1 > 1) {
        connect2d[nlocal_connect].neigh_p1 = tcp->get(np1,connect2d[nlocal_connect].indexp1);
        for (k = 0; k < np1; k++)
          connect2d[nlocal_connect].neigh_p1[k] = (tagint) ubuf(buf[m++]).i;
      } else connect2d[nlocal_connect].neigh_p1 = nullptr;
      if (np2 > 1) {
        connect2d[nlocal_connect].neigh_p2 = tcp->get(np2,connect2d[nlocal_connect].indexp2);
        for (k = 0; k < np2; k++)
          connect2d[nlocal_connect].neigh_p2[k] = (tagint) ubuf(buf[m++]).i;
      } else connect2d[nlocal_connect].neigh_p2 = nullptr;

      connect2d[nlocal_connect].flags = (int) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].ilocal = nlocal;
      cindex[nlocal] = nlocal_connect++;
    }

  } else {
    int ne1,ne2,ne3,nc1,nc2,nc3;

    flag = (int) ubuf(buf[m++]).i;
    if (flag == 0) cindex[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      ne1 = (int) ubuf(buf[m++]).i;
      ne2 = (int) ubuf(buf[m++]).i;
      ne3 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].ne1 = ne1;
      connect3d[nlocal_connect].ne2 = ne2;
      connect3d[nlocal_connect].ne3 = ne3;
      if (ne1 > 1) {
        connect3d[nlocal_connect].neigh_e1 = tcp->get(ne1,connect3d[nlocal_connect].indexe1);
        for (k = 0; k < ne1; k++)
          connect3d[nlocal_connect].neigh_e1[k] = (tagint) ubuf(buf[m++]).i;
      } else connect3d[nlocal_connect].neigh_e1 = nullptr;
      if (ne2 > 1) {
        connect3d[nlocal_connect].neigh_e2 = tcp->get(ne2,connect3d[nlocal_connect].indexe2);
        for (k = 0; k < ne2; k++)
          connect3d[nlocal_connect].neigh_e2[k] = (tagint) ubuf(buf[m++]).i; 
      } else connect3d[nlocal_connect].neigh_e2 = nullptr;
      if (ne3 > 1) {
        connect3d[nlocal_connect].neigh_e3 = tcp->get(ne3,connect3d[nlocal_connect].indexe3);
        for (k = 0; k < ne3; k++)
          connect3d[nlocal_connect].neigh_e3[k] = (tagint) ubuf(buf[m++]).i;
      } else connect3d[nlocal_connect].neigh_e3 = nullptr;

      nc1 = (int) ubuf(buf[m++]).i;
      nc2 = (int) ubuf(buf[m++]).i;
      nc3 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].nc1 = nc1;
      connect3d[nlocal_connect].nc2 = nc2;
      connect3d[nlocal_connect].nc3 = nc3;
      if (nc1 > 1) {
        connect3d[nlocal_connect].neigh_c1 = tcp->get(nc1,connect3d[nlocal_connect].indexc1);
        for (k = 0; k < nc1; k++)
          connect3d[nlocal_connect].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
      } else connect3d[nlocal_connect].neigh_c1 = nullptr;
      if (nc2 > 1) {
        connect3d[nlocal_connect].neigh_c2 = tcp->get(nc2,connect3d[nlocal_connect].indexc2);
        for (k = 0; k < nc2; k++)
          connect3d[nlocal_connect].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
      } else connect3d[nlocal_connect].neigh_c2 = nullptr;
      if (nc3 > 1) {
        connect3d[nlocal_connect].neigh_c3 = tcp->get(nc3,connect3d[nlocal_connect].indexc3);
        for (k = 0; k < nc3; k++)
          connect3d[nlocal_connect].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;
      } else connect3d[nlocal_connect].neigh_c3 = nullptr;
        
      connect3d[nlocal_connect].flags = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].ilocal = nlocal;
      cindex[nlocal] = nlocal_connect++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   // NOTE: need to also tally pages in tcp
------------------------------------------------------------------------- */

double FixSurfaceLocal::memory_usage()
{
  double bytes = 0.0;
  if (dimension == 2) bytes = nmax_connect*sizeof(Connect2d);
  else bytes = nmax_connect*sizeof(Connect3d);
  bytes += atom->nmax * sizeof(int);               // cindex vector
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

  // allocate connection info for my owned lines
  // initialize cindex for all my particles, including lines
  // comm->ring() operation will populate connection info
  
  nlocal_connect = nmax_connect = nline;
  nghost_connect = 0;
  grow_connect();

  for (i = 0; i < nlocal; i++) {
    cindex[i] = line[i];
    if (line[i] < 0) continue;
    j = line[i];
    connect2d[j].ilocal = i;
  }

  // NOTE: have to initialize rest of values in Connect2d for each line I own
  
  // circulate my line info to all other procs, including self

  fptr = this;
  nmatch = errormatch = 0;

  /*
  if (ebuf)
    comm->ring(nline,5*sizeof(double),(void *) &ebuf[0][0],0,linematch,nullptr);
  else
    comm->ring(nline,5*sizeof(double),nullptr,0,linematch,nullptr);
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

          if (rsq < epssq) {
            (*nmatch)++;

            /*
            if (ptwhich == 1) {
              if (connect2d[iconnect].neigh_p1 > 0) (*errormatch)++;
              else connect2d[iconnect].neigh_p1 = tagother;
            } else if (ptwhich == 2) {
              if (connect2d[iconnect].neigh_p2 > 0) (*errormatch)++;
              else connect2d[iconnect].neigh_p2 = tagother;
            }
            */
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

  // allocate connection info for my owned triangles
  // initialize cindex for all my particles, including triangles
  // comm->ring() operations will populate connection info
  // do this twice:
  //   first time to accumulate corner counts
  //   allocate corner connection lists
  //   second time to fill lists (and redo everything else)

  nlocal_connect = nmax_connect = ntri;
  nghost_connect = 0;
  grow_connect();

  for (i = 0; i < nlocal; i++) {
    cindex[i] = tri[i];
    if (tri[i] < 0) continue;
    j = tri[i];
    connect3d[j].flags = 0;
  }

  // NOTE: have to initialize rest of values in Connect3d for each tri I own

  fptr = this;
  nmatch1 = nmatch2 = errormatch1 = errormatch2 = 0;
  vecflag = 0;

  /*
  if (tbuf)
    comm->ring(ntri,10*sizeof(double),(void *) &tbuf[0][0],0,trimatch,nullptr);
  else
    comm->ring(ntri,10*sizeof(double),nullptr,0,trimatch,nullptr);
  */

  for (i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    j = tri[i];
    connect3d[j].indexc1 = connect3d[j].indexc2 = connect3d[j].indexc3 = -1;
    if (connect3d[j].nc1)
      connect3d[j].neigh_c1 = tcp->get(connect3d[j].nc1,connect3d[j].indexc1);
    else connect3d[j].neigh_c1 = nullptr;
    if (connect3d[j].nc2)
      connect3d[j].neigh_c2 = tcp->get(connect3d[j].nc2,connect3d[j].indexc2);
    else connect3d[j].neigh_c2 = nullptr;
    if (connect3d[j].nc3)
      connect3d[j].neigh_c3 = tcp->get(connect3d[j].nc3,connect3d[j].indexc3);
    else connect3d[j].neigh_c3 = nullptr;
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
    comm->ring(ntri,10*sizeof(double),(void *) &tbuf[0][0],0,trimatch,nullptr);
  else
    comm->ring(ntri,10*sizeof(double),nullptr,0,trimatch,nullptr);
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
              /*
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
              */
              
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
   check if distributed line/tri particles already exist
------------------------------------------------------------------------- */

int FixSurfaceLocal::check_exist()
{
  int *surf = atom->line;
  if (dimension == 3) surf = atom->tri;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (surf[i] >= 0) flag = 1;
  
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);

  if (flagall) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   extract lines or surfs from molecule template ID for one or more molecules
   concatenate into single list of lines and tris
   create list of unique points using hash
------------------------------------------------------------------------- */

void FixSurfaceLocal::extract_from_molecules(char *molID)
{
  // populate global point/line/tri data structs

  npoints = nlines = ntris = 0;
  int maxpoints = 0;

  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for fix surface/local does not exist");

  // loop over one or more molecules in molID

  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;

  for (int m = 0; m < nmol; m++) {
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface/local molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface/local molecule must have triangles");

    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    nlines += nl;
    ntris += nt;
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface/global:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface/global:tris");

    // create a map
    // key = xyz coords of a point
    // value = index into unique points vector

    std::map<std::tuple<double,double,double>,int> hash;

    // offset line/tri index lists by previous npoints
    // pi,p2,p3 are C-style indices into points vector

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      double **epts = onemols[m]->lines;
      int iline = nlines - nl;

      for (int i = 0; i < nl; i++) {
        lines[iline].mol = molline[i];
        lines[iline].type = typeline[i];

        auto key = std::make_tuple(epts[i][0],epts[i][1],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/local:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][0];
          points[npoints].x[1] = epts[i][1];
          points[npoints].x[2] = 0.0;
          lines[iline].p1 = npoints;
          npoints++;
        } else lines[iline].p1 = hash[key];

        key = std::make_tuple(epts[i][2],epts[i][3],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/local:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][2];
          points[npoints].x[1] = epts[i][3];
          points[npoints].x[2] = 0.0;
          lines[iline].p2 = npoints;
          npoints++;
        } else lines[iline].p2 = hash[key];

        iline++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      double **cpts = onemols[m]->tris;
      int itri = ntris - nt;

      for (int i = 0; i < nt; i++) {
        tris[itri].mol = moltri[i];
        tris[itri].type = typetri[i];

        auto key = std::make_tuple(cpts[i][0],cpts[i][1],cpts[i][2]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/local:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][0];
          points[npoints].x[1] = cpts[i][1];
          points[npoints].x[2] = cpts[i][2];
          tris[itri].p1 = npoints;
          npoints++;
        } else tris[itri].p1 = hash[key];

        key = std::make_tuple(cpts[i][3],cpts[i][4],cpts[i][5]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/local:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][3];
          points[npoints].x[1] = cpts[i][4];
          points[npoints].x[2] = cpts[i][5];
          tris[itri].p2 = npoints;
          npoints++;
        } else tris[itri].p2 = hash[key];

        key = std::make_tuple(cpts[i][6],cpts[i][7],cpts[i][8]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/local:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][6];
          points[npoints].x[1] = cpts[i][7];
          points[npoints].x[2] = cpts[i][8];
          tris[itri].p3 = npoints;
          npoints++;
        } else tris[itri].p3 = hash[key];

        itri++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   extract triangles from an STL file, can be text or binary
   create list of unique points using hash
------------------------------------------------------------------------- */

void FixSurfaceLocal::extract_from_stlfile(char *filename)
{
  if (dimension == 2)
    error->all(FLERR,"Fix surface/global cannot use an STL file for 2d simulations");

  // read tris from STL file
  // stltris = tri coords internal to STL reader
  
  STLReader *stl = new STLReader(lmp);
  double **stltris;
  ntris = stl->read_file(filename,stltris);

  // create points and tris data structs

  npoints = 0;
  int maxpoints = 0;

  tris = (Tri *) memory->smalloc(ntris*sizeof(Tri),"surface/global:tris");
  
  // create a map
  // key = xyz coords of a point
  // value = index into unique points vector
  
  std::map<std::tuple<double,double,double>,int> hash;

  // loop over STL tris
  // populate points and tris data structs
  // set molecule and type of tri = 1
  
  for (int itri = 0; itri < ntris; itri++) {
    tris[itri].mol = 1;
    tris[itri].type = 1;

    auto key = std::make_tuple(stltris[itri][0],stltris[itri][1],stltris[itri][2]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][0];
      points[npoints].x[1] = stltris[itri][1];
      points[npoints].x[2] = stltris[itri][2];
      tris[itri].p1 = npoints;
      npoints++;
    } else tris[itri].p1 = hash[key];

    key = std::make_tuple(stltris[itri][3],stltris[itri][4],stltris[itri][5]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][3];
      points[npoints].x[1] = stltris[itri][4];
      points[npoints].x[2] = stltris[itri][5];
      tris[itri].p2 = npoints;
      npoints++;
    } else tris[itri].p2 = hash[key];
    
    key = std::make_tuple(stltris[itri][6],stltris[itri][7],stltris[itri][8]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][6];
      points[npoints].x[1] = stltris[itri][7];
      points[npoints].x[2] = stltris[itri][8];
      tris[itri].p3 = npoints;
      npoints++;
    } else tris[itri].p3 = hash[key];
  }
  
  // delete STL reader
  
  delete stl;
}

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for all lines
   stored in connect2dall
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity2d_global()
{
  // allocate and initilize global connectivity list
  // ilocal will be initialized when assign to procs

  connect2dall = (Connect2d *) memory->smalloc(nlines*sizeof(Connect2d),
                                               "surface/local:connect2dall");

  // setup end point connectivity lists
  // count # of lines containing each end point
  // create ragged 2d array to contain all point indices, then fill it
  // set neigh_p12 vector ptrs in connect2d to rows of ragged array
  // np12 counts and vectors include self line
  
  int *counts;
  memory->create(counts,npoints,"surface/local:count");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    counts[lines[i].p1]++;
    counts[lines[i].p2]++;
  }

  memory->create_ragged(plist,npoints,counts,"surface/local:plist");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    plist[lines[i].p1][counts[lines[i].p1]++] = i;
    plist[lines[i].p2][counts[lines[i].p2]++] = i;
  }

  // note that in connect2dall, neigh_p12 stores global line indices (0 to Nlines-1)
  // in final connect2d, neigh_p12 will store line IDs

  for (int i = 0; i < nlines; i++) {
    connect2dall[i].np1 = counts[lines[i].p1];
    if (connect2dall[i].np1 == 1) connect2dall[i].neigh_p1 = nullptr;
    else connect2dall[i].neigh_p1 = plist[lines[i].p1];

    connect2dall[i].np2 = counts[lines[i].p2];
    if (connect2dall[i].np2 == 1) connect2dall[i].neigh_p2 = nullptr;
    else connect2dall[i].neigh_p2 = plist[lines[i].p2];
  }

  memory->destroy(counts);
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all triangles
   stored in connect3dall
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity3d_global()
{
  int p1,p2,p3;

  // allocate and initilize global connectivity list
  // ilocal will be initialized when assign to procs

  connect3dall = (Connect3d *) memory->smalloc(ntris*sizeof(Connect3d),
                                               "surface/local:connect3d");
  int **tri2edge;
  memory->create(tri2edge,ntris,3,"surfface/global::tri2edge");

  // create a map
  // key = <p1,p2> indices of 2 points, in either order
  // value = index into count of unique edges

  std::map<std::tuple<int,int>,int> hash;
  int nedges = 0;

  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    auto key1 = std::make_tuple(p1,p2);
    auto key2 = std::make_tuple(p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][0] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][0] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][0] = hash[key2];

    key1 = std::make_tuple(p2,p3);
    key2 = std::make_tuple(p3,p2);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][1] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][1] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][1] = hash[key2];

    key1 = std::make_tuple(p3,p1);
    key2 = std::make_tuple(p1,p3);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][2] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][2] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][2] = hash[key2];
  }

  // setup tri edge connectivity lists
  // count # of tris containing each edge
  // create ragged 2d array to contain all edge indices, then fill it
  // set neigh_e123 vector ptrs in connect3d to rows of ragged array
  // ne123 counts and vectors include self tri

  int *counts;
  memory->create(counts,nedges,"surface/local:count");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tri2edge[i][0]]++;
    counts[tri2edge[i][1]]++;
    counts[tri2edge[i][2]]++;
  }

  memory->create_ragged(elist,nedges,counts,"surface/local:elist");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    elist[tri2edge[i][0]][counts[tri2edge[i][0]]++] = i;
    elist[tri2edge[i][1]][counts[tri2edge[i][1]]++] = i;
    elist[tri2edge[i][2]][counts[tri2edge[i][2]]++] = i;
  }

  // note that in connect3dall, neigh_e123 stores global tri indices (0 to Ntri-1)
  // in final connect3d, neigh_e123 will store tri IDs

  for (int i = 0; i < ntris; i++) {
    connect3dall[i].ne1 = counts[tri2edge[i][0]];
    if (connect3dall[i].ne1 == 1) connect3dall[i].neigh_e1 = nullptr;
    else connect3dall[i].neigh_e1 = elist[tri2edge[i][0]];

    connect3dall[i].ne2 = counts[tri2edge[i][1]];
    if (connect3dall[i].ne2 == 1) connect3dall[i].neigh_e2 = nullptr;
    else connect3dall[i].neigh_e2 = elist[tri2edge[i][1]];

    connect3dall[i].ne3 = counts[tri2edge[i][2]];
    if (connect3dall[i].ne3 == 1) connect3dall[i].neigh_e3 = nullptr;
    else connect3dall[i].neigh_e3 = elist[tri2edge[i][2]];
  }

  memory->destroy(counts);
  memory->destroy(tri2edge);

  // setup corner point connectivity lists
  // count # of tris containing each point
  // create ragged 2d array to contain all tri indices, then fill it
  // set neigh_c123 vector ptrs in connect3d to rows of ragged array
  // nc123 counts and vectors include self tri

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
    clist[tris[i].p1][counts[tris[i].p1]++] = i;
    clist[tris[i].p2][counts[tris[i].p2]++] = i;
    clist[tris[i].p3][counts[tris[i].p3]++] = i;
  }

  // note that in connect3dall, neigh_c123 stores global tri indices (0 to Ntri-1)
  // in final connect3d, neigh_c123 will store tri IDs

  for (int i = 0; i < ntris; i++) {
    connect3dall[i].nc1 = counts[tris[i].p1];
    if (connect3dall[i].nc1 == 1) connect3dall[i].neigh_c1 = nullptr;
    else connect3dall[i].neigh_c1 = clist[tris[i].p1];

    connect3dall[i].nc2 = counts[tris[i].p2];
    if (connect3dall[i].nc2 == 1) connect3dall[i].neigh_c2 = nullptr;
    else connect3dall[i].neigh_c2 = clist[tris[i].p2];

    connect3dall[i].nc3 = counts[tris[i].p3];
    if (connect3dall[i].nc3 == 1) connect3dall[i].neigh_c3 = nullptr;
    else connect3dall[i].neigh_c3 = clist[tris[i].p3];
  }

  memory->destroy(counts);
}

/* ----------------------------------------------------------------------
   assign lines and their connectivity from global data structs to each proc
   logic for proc assign follows Atom::data_atoms(), as if read from data file
------------------------------------------------------------------------- */

void FixSurfaceLocal::assign2d()
{
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
  // set cindex = -1 for all current particles
  
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (int i = 0; i < nlocal; i++) {
    idmax = MAX(tag[i],idmax);
    cindex[i] = -1;
  }
  tagint idmaxall;
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // loop over all global lines
  // compute their midpoints
  // keep lines belonging to this proc

  int j,n,num;
  imageint imagedata;
  double xmid[3],lamda[3];
  double *coord,*x1,*x2;
  tagint *global,*local;

  char pts[4][32];
  char *values[4];
  for (int i = 0; i < 4; i++) values[i] = &pts[i][0];
  std::vector<std::string> svalues(5);

  for (int i = 0; i < nlines; i++) {

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

      // svalues[0] is not used here
      // used when caller of data_atom_bonus() is via read_data command

      svalues[1] = (const char *) values[0];
      svalues[2] = (const char *) values[1];
      svalues[3] = (const char *) values[2];
      svalues[4] = (const char *) values[3];

      avec_line->data_atom_bonus(n,svalues);

      // copy connectivity for global line I to local list
      // replace neigh p1/p2 ptrs into plist with tcp pool allocs
      // reset endpt line IDs to new IDs beyond idmaxall

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect2d[nlocal_connect],&connect2dall[i],sizeof(Connect2d));

      num = connect2d[nlocal_connect].np1;
      global = connect2d[nlocal_connect].neigh_p1;
      if (num > 1) {
        local = connect2d[nlocal_connect].neigh_p1 =
          tcp->get(num,connect2d[nlocal_connect].indexp1);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect2d[nlocal_connect].neigh_p1 = nullptr;

      num = connect2d[nlocal_connect].np2;
      global = connect2d[nlocal_connect].neigh_p2;
      if (num > 1) {
        local = connect2d[nlocal_connect].neigh_p2 =
          tcp->get(num,connect2d[nlocal_connect].indexp2);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect2d[nlocal_connect].neigh_p2 = nullptr;

      connect2d[nlocal_connect].ilocal = n;
      cindex[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect2dall);

  // recreate atom map that includes added lines

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
  // set cindex = -1 for all current particles

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (int i = 0; i < nlocal; i++) {
    idmax = MAX(tag[i],idmax);
    cindex[i] = -1;
  }
  tagint idmaxall;
  MPI_Allreduce(&idmax,&idmaxall,1,MPI_LMP_TAGINT,MPI_MAX,world);

  // loop over all global triangles
  // compute their center points
  // keep lines belonging to this proc

  int j,n,num;
  imageint imagedata;
  double xmid[3],lamda[3];
  double *coord,*x1,*x2,*x3;
  tagint *global,*local;

  char pts[9][32];
  char *values[9];
  for (int i = 0; i < 9; i++) values[i] = &pts[i][0];
  std::vector<std::string> svalues(10);

  for (int i = 0; i < ntris; i++) {

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

      // svalues[0] is not used here
      // used when caller of data_atom_bonus() is via read_data command

      svalues[1] = (const char *) values[0];
      svalues[2] = (const char *) values[1];
      svalues[3] = (const char *) values[2];
      svalues[4] = (const char *) values[3];
      svalues[5] = (const char *) values[4];
      svalues[6] = (const char *) values[5];
      svalues[7] = (const char *) values[6];
      svalues[8] = (const char *) values[7];
      svalues[9] = (const char *) values[8];

      avec_tri->data_atom_bonus(n,svalues);

      // copy connectivity for globa triangle I to local list
      // replace neigh e1/e2/e3 ptrs into elist with tcp pool allocs
      // replace neigh c1/c2/c3 ptrs into clist with tcp pool allocs
      // reset edge and corner tri IDs to new IDs beyond idmaxall

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect3d[nlocal_connect],&connect3dall[i],sizeof(Connect3d));
      
      num = connect3d[nlocal_connect].ne1;
      global = connect3d[nlocal_connect].neigh_e1;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_e1 =
          tcp->get(num,connect3d[nlocal_connect].indexe1);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_e1 = nullptr;

      num = connect3d[nlocal_connect].ne2;
      global = connect3d[nlocal_connect].neigh_e2;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_e2 =
          tcp->get(num,connect3d[nlocal_connect].indexe2);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_e2 = nullptr;

      num = connect3d[nlocal_connect].ne3;
      global = connect3d[nlocal_connect].neigh_e3;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_e3 =
          tcp->get(num,connect3d[nlocal_connect].indexe3);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_e3 = nullptr;

      num = connect3d[nlocal_connect].nc1;
      global = connect3d[nlocal_connect].neigh_c1;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_c1 =
          tcp->get(num,connect3d[nlocal_connect].indexc1);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_c1 = nullptr;

      num = connect3d[nlocal_connect].nc2;
      global = connect3d[nlocal_connect].neigh_c2;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_c2 =
          tcp->get(num,connect3d[nlocal_connect].indexc2);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_c2 = nullptr;

      num = connect3d[nlocal_connect].nc3;
      global = connect3d[nlocal_connect].neigh_c3;
      if (num > 1) {
        local = connect3d[nlocal_connect].neigh_c3 =
          tcp->get(num,connect3d[nlocal_connect].indexc3);
        for (j = 0; j < num; j++)
          local[j] = global[j] + idmaxall + 1;
      } else connect3d[nlocal_connect].neigh_c3 = nullptr;

      connect3d[nlocal_connect].ilocal = n;
      cindex[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect3dall);

  // recreate atom map that includes added tris

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}
