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
   Contributing author: Daniel Schwen
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_voronoi_atom.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "variable.h"
#include "input.h"

#include <vector>

using namespace LAMMPS_NS;
using namespace voro;

/* ---------------------------------------------------------------------- */

ComputeVoronoi::ComputeVoronoi(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  int sgroup;
  
  size_peratom_cols = 2;
  peratom_flag = 1;

  surface = VOROSURF_NONE;
  maxedge = 0;
  fthresh = ethresh = 0.0;
  radstr = NULL;
  onlyGroup = false;

  int iarg = 3;
  while ( iarg<narg ) {
    if (strcmp(arg[iarg],"only_group") == 0) {
      onlyGroup = true;
      iarg++;
    } 
    else if (strcmp(arg[iarg],"radius") == 0) {
      if (iarg + 2 > narg || strstr(arg[iarg+1],"v_") != arg[iarg+1] ) 
        error->all(FLERR,"Illegal compute voronoi/atom command");
      int n = strlen(&arg[iarg+1][2]) + 1;
      radstr = new char[n];
      strcpy(radstr,&arg[iarg+1][2]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg],"surface") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR,"Illegal compute voronoi/atom command");
      // group all is a special case where we just skip group testing
      if(strcmp(arg[iarg+1], "all") == 0) {
        surface = VOROSURF_ALL;
      } else {
        sgroup = group->find(arg[iarg+1]);
        if (sgroup == -1) 
          error->all(FLERR,"Could not find compute/voronoi surface group ID");
        sgroupbit = group->bitmask[sgroup]; 
        surface = VOROSURF_GROUP;
      }
      size_peratom_cols = 3;
      iarg += 2;
    } else if (strcmp(arg[iarg], "edge_histo") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR,"Illegal compute voronoi/atom command");
      maxedge = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "face_threshold") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR,"Illegal compute voronoi/atom command");
      fthresh = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "edge_threshold") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR,"Illegal compute voronoi/atom command");
      ethresh = atof(arg[iarg+1]);
      iarg += 2;
    } 
    else 
      error->all(FLERR,"Illegal compute voronoi/atom command");
  }

  nmax = rmax = 0;
  edge = rfield = sendvector = NULL;
  voro = NULL;

  if ( maxedge > 0 ) {
    vector_flag = 1;
    size_vector = maxedge+1;
    memory->create(edge,maxedge+1,"voronoi/atom:edge");
    memory->create(sendvector,maxedge+1,"voronoi/atom:sendvector");
    vector = edge;
  }
}

/* ---------------------------------------------------------------------- */

ComputeVoronoi::~ComputeVoronoi()
{
  memory->destroy(edge);
  memory->destroy(rfield);
  memory->destroy(sendvector);
  memory->destroy(voro);
  delete[] radstr;
}

/* ---------------------------------------------------------------------- */

void ComputeVoronoi::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"voronoi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute voronoi/atom command");
}

/* ----------------------------------------------------------------------
   gather compute vector data from other nodes
------------------------------------------------------------------------- */

void ComputeVoronoi::compute_peratom()
{
  int i, j;
  const double e = 0.01;

  invoked_peratom = update->ntimestep;

  // grow per atom array if necessary

  int nlocal = atom->nlocal;
  if (nlocal > nmax) {
    memory->destroy(voro);
    nmax = atom->nmax;
    memory->create(voro,nmax,size_peratom_cols,"voronoi/atom:voro");
    array_atom = voro;
  }

  // in the onlyGroup mode we are not setting values for all atoms later in the voro loop
  // initialize everything to zero here
  if (onlyGroup) {
    if (surface == VOROSURF_NONE) 
      for (i = 0; i < nlocal; i++) voro[i][0] = voro[i][1] = 0.0;
    else
      for (i = 0; i < nlocal; i++) voro[i][0] = voro[i][1] = voro[i][2] = 0.0;
  }

  double *sublo = domain->sublo, *sublo_lamda = domain->sublo_lamda, *boxlo = domain->boxlo;
  double *subhi = domain->subhi, *subhi_lamda = domain->subhi_lamda, *boxhi = domain->boxhi;
  double *cut = comm->cutghost;
  double sublo_bound[3], subhi_bound[3], cut_bound[3];
  double **x = atom->x;
   
  // setup bounds for voro++ domain for orthogonal and triclinic simulation boxes
  if( domain->triclinic ) {
    // triclinic box: embed parallelepiped into orthogonal voro++ domain
    double mx, my, sxy,sxz,syz;
    mx = (boxhi[0]-boxlo[0])/(subhi[0]-sublo[0]);
    my = (boxhi[1]-boxlo[1])/(subhi[1]-sublo[1]);
    sxy = domain->xy/mx;
    sxz = domain->xz/mx;
    syz = domain->yz/my;

    // cutghost is in lamda coordinates for triclinic boxes, use subxx_lamda
    double *h = domain->h, cuttri[3];
    sublo_bound[0] = h[0]*sublo_lamda[0] + h[5]*sublo_lamda[1] + h[4]*sublo_lamda[2] + boxlo[0];
    sublo_bound[1] = h[1]*sublo_lamda[1] + h[3]*sublo_lamda[2] + boxlo[1];
    sublo_bound[2] = h[2]*sublo_lamda[2] + boxlo[2];
    subhi_bound[0] = h[0]*subhi_lamda[0] + h[5]*subhi_lamda[1] + h[4]*subhi_lamda[2] + boxlo[0];
    subhi_bound[1] = h[1]*subhi_lamda[1] + h[3]*subhi_lamda[2] + boxlo[1];
    subhi_bound[2] = h[2]*subhi_lamda[2] + boxlo[2];
    cut_bound[0] = h[0]*cut[0] + h[5]*cut[1] + h[4]*cut[2];
    cut_bound[1] = h[1]*cut[1] + h[3]*cut[2];
    cut_bound[2] = h[2]*cut[2];

  } else {
    // orthogonal box
    for( i=0; i<3; ++i ) {
      sublo_bound[i] = sublo[i];
      subhi_bound[i] = subhi[i];
      cut_bound[i] = cut[i];
    }
  }

  // n = # of voro++ spatial hash cells (with approximately cubic cells)
  int nall = nlocal + atom->nghost;
  double n[3], V;
  for( i=0; i<3; ++i ) n[i] = subhi_bound[i] - sublo_bound[i];
  V = n[0]*n[1]*n[2];
  for( i=0; i<3; ++i ) {
    n[i] = round( n[i]*pow( double(nall)/(V*8.0), 0.333333 ) );
    n[i] = n[i]==0 ? 1 : n[i];
  } 

  // clear edge statistics
  for (i = 0; i < maxedge; ++i) edge[i]=0;

  // initialize voro++ container
  // preallocates 8 atoms per cell
  // voro++ allocates more memory if needed
  voronoicell_neighbor c;
  int *mask = atom->mask;
  if(radstr) {
    // check and fetch atom style variable data
    int radvar = input->variable->find(radstr);
    if (radvar < 0)
      error->all(FLERR,"Variable name for voronoi radius does not exist");
    if (!input->variable->atomstyle(radvar))
      error->all(FLERR,"Variable for voronoi radius is not atom style");
    // prepare destination buffer for variable evaluation
    if (nlocal > rmax) {
      memory->destroy(rfield);
      rmax = atom->nmax;
      memory->create(rfield,rmax,"voronoi/atom:rfield");
    }
    // compute atom style radius variable
    input->variable->compute_atom(radvar,0,rfield,1,0);

    // communicate values to ghost atoms of neighboring nodes
    comm->forward_comm_compute(this);

    // polydisperse voro++ container
    container_poly con(sublo_bound[0]-cut_bound[0]-e,subhi_bound[0]+cut_bound[0]+e,
                       sublo_bound[1]-cut_bound[1]-e,subhi_bound[1]+cut_bound[1]+e,
                       sublo_bound[2]-cut_bound[2]-e,subhi_bound[2]+cut_bound[2]+e,
                       int(n[0]),int(n[1]),int(n[2]),false,false,false,8); 

    // pass coordinates for local and ghost atoms to voro++
    for (i = 0; i < nall; i++)
      if( !onlyGroup || (mask[i] & groupbit) )
        con.put(i,x[i][0],x[i][1],x[i][2],rfield[i]);

    // invoke voro++ and fetch results for owned atoms in group
    c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c,cl)) {
      i = cl.pid();
      processCell(c,i);
    } while (cl.inc());
  } else {
    // monodisperse voro++ container
    container con(sublo_bound[0]-cut_bound[0]-e,subhi_bound[0]+cut_bound[0]+e,
                  sublo_bound[1]-cut_bound[1]-e,subhi_bound[1]+cut_bound[1]+e,
                  sublo_bound[2]-cut_bound[2]-e,subhi_bound[2]+cut_bound[2]+e,
                  int(n[0]),int(n[1]),int(n[2]),false,false,false,8); 

    // pass coordinates for local and ghost atoms to voro++
    for (i = 0; i < nall; i++)
      if( !onlyGroup || (mask[i] & groupbit) )
        con.put(i,x[i][0],x[i][1],x[i][2]);

    // invoke voro++ and fetch results for owned atoms in group
    c_loop_all cl(con);
    if (cl.start()) do if (con.compute_cell(c,cl)) {
      i = cl.pid();
      processCell(c,i);
    } while (cl.inc());
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
void ComputeVoronoi::processCell(voronoicell_neighbor &c, int i) 
{
  int j,k, *mask = atom->mask;
  std::vector<int> neigh, norder, vlist;
  std::vector<double> narea, vcell;
  bool have_narea = false;

  // zero out surface area if surface computation was requested
  if (surface != VOROSURF_NONE && !onlyGroup) voro[i][2] = 0.0;

  if (i < atom->nlocal && (mask[i] & groupbit)) {
    // cell volume
    voro[i][0] = c.volume();

    // number of cell faces
    c.neighbors(neigh);
    if (fthresh > 0) {
      // count only faces above area threshold
      c.face_areas(narea);
      have_narea = true;
      voro[i][1] = 0.0;
      for (j=0; j<narea.size(); ++j)  
        if (narea[j] > fthresh) voro[i][1] += 1.0;
    } else {
      // unthresholded face count
      voro[i][1] = neigh.size();
    }

    // cell surface area
    if (surface == VOROSURF_ALL) {
      voro[i][2] = c.surface_area();
    } else if (surface == VOROSURF_GROUP) {
      if (!have_narea) c.face_areas(narea);
      voro[i][2] = 0.0;
      // loop over all faces (neighbors) and check if they are in the surface group
      for (j=0; j<voro[i][1]; ++j)  
        if (mask[neigh[j]] & sgroupbit) voro[i][2] += narea[j];
    }

    // histogram of number of face edges
    if (maxedge>0) {
      if (ethresh > 0) {
        // count only edges above length threshold 
        c.vertices(vcell);
        c.face_vertices(vlist); // for each face: vertex count followed list of vertex indices (n_1,v1_1,v2_1,v3_1,..,vn_1,n_2,v2_1,...)
        double dx, dy, dz, r2, t2 = ethresh*ethresh;
        for( j=0; j<vlist.size(); j+=vlist[j]+1 ) { 
          int a, b, nedge = 0;
          // vlist[j] contains number of vertex indices for the current face
          for( k=0; k<vlist[j]; ++k ) { 
            a = vlist[j+1+k];              // first vertex in edge 
            b = vlist[j+1+(k+1)%vlist[j]]; // second vertex in edge (possible wrap around to first vertex in list)
            dx = vcell[a*3]   - vcell[b*3];
            dy = vcell[a*3+1] - vcell[b*3+1];
            dz = vcell[a*3+2] - vcell[b*3+2];
            r2 = dx*dx+dy*dy+dz*dz;
            if (r2 > t2) nedge++;
          }
          // counted edges above threshold, now put into the correct bin
          if (nedge>0) {
            if (nedge<=maxedge) 
              edge[nedge-1]++;
            else
              edge[maxedge]++;
          }
        }
      } else {
        // unthresholded edge counts
        c.face_orders(norder);
        for (j=0; j<voro[i][1]; ++j)
          if (norder[j]>0) {
            if (norder[j]<=maxedge) 
              edge[norder[j]-1]++;
            else
              edge[maxedge]++;
          }
      }
    }
  } else if (i < atom->nlocal) voro[i][0] = voro[i][1] = 0.0;
}

double ComputeVoronoi::memory_usage()
{
  double bytes = size_peratom_cols * nmax * sizeof(double);
  return bytes;
}

void ComputeVoronoi::compute_vector()
{
  invoked_vector = update->ntimestep;
  if( invoked_peratom < invoked_vector ) compute_peratom();

  for( int i=0; i<size_vector; ++i ) sendvector[i] = edge[i];
  MPI_Allreduce(sendvector,edge,size_vector,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

int ComputeVoronoi::pack_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int i,m=0;
  for (i = 0; i < n; i++) buf[m++] = rfield[list[i]];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeVoronoi::unpack_comm(int n, int first, double *buf)
{
  int i,last,m=0;
  last = first + n;
  for (i = first; i < last; i++) rfield[i] = buf[m++];
}

