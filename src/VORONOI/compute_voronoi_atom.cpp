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
#include "force.h"

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
  occupation = false;

  con_mono = NULL;
  con_poly = NULL;
  tags = occvec = sendocc = lroot = lnext = NULL;

  int iarg = 3;
  while ( iarg<narg ) {
    if (strcmp(arg[iarg], "occupation") == 0) {
      occupation = true;
      iarg++;
    } 
    else if (strcmp(arg[iarg], "only_group") == 0) {
      onlyGroup = true;
      iarg++;
    } 
    else if (strcmp(arg[iarg], "radius") == 0) {
      if (iarg + 2 > narg || strstr(arg[iarg+1],"v_") != arg[iarg+1] ) error->all(FLERR,"Missing atom style variable for radical voronoi tesselation radius.");
      int n = strlen(&arg[iarg+1][2]) + 1;
      radstr = new char[n];
      strcpy(radstr,&arg[iarg+1][2]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg], "surface") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Missing group name after keyword 'surface'.");
      // group all is a special case where we just skip group testing
      if(strcmp(arg[iarg+1], "all") == 0) {
        surface = VOROSURF_ALL;
      } else {
        sgroup = group->find(arg[iarg+1]);
        if (sgroup == -1) error->all(FLERR,"Could not find compute/voronoi surface group ID");
        sgroupbit = group->bitmask[sgroup]; 
        surface = VOROSURF_GROUP;
      }
      size_peratom_cols = 3;
      iarg += 2;
    } else if (strcmp(arg[iarg], "edge_histo") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Missing maximum edge count after keyword 'edge_histo'.");
      maxedge = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "face_threshold") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Missing minimum face area after keyword 'face_threshold'.");
      fthresh = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "edge_threshold") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Missing minimum edge length after keyword 'edge_threshold'.");
      ethresh = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } 
    else 
      error->all(FLERR,"Illegal compute voronoi/atom command");
  }

  if (occupation && ( surface!=VOROSURF_NONE || maxedge>0 ) ) 
    error->all(FLERR,"Illegal compute voronoi/atom command (occupation and (surface or edges))");

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

  // voro++ container classes
  delete con_mono;
  delete con_poly;

  // occupation analysis stuff
  memory->destroy(lroot);
  memory->destroy(lnext);
  memory->destroy(occvec);
#ifdef NOTINPLACE
  memory->destroy(sendocc);
#endif
  memory->destroy(tags);
}

/* ---------------------------------------------------------------------- */

void ComputeVoronoi::init()
{
}

/* ----------------------------------------------------------------------
   gather compute vector data from other nodes
------------------------------------------------------------------------- */

void ComputeVoronoi::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow per atom array if necessary
  int nlocal = atom->nlocal;
  if (nlocal > nmax) {
    memory->destroy(voro);
    nmax = atom->nmax;
    memory->create(voro,nmax,size_peratom_cols,"voronoi/atom:voro");
    array_atom = voro;
  }

  // decide between occupation or per-frame tesselation modes
  if (occupation) {
    // build cells only once
    int i, j, 
        nall = nlocal + atom->nghost;
    if (con_mono==NULL && con_poly==NULL) {
      // generate the voronoi cell network for the initial structure
      buildCells();

      // save tags of atoms (i.e. of each voronoi cell)
      memory->create(tags,nall,"voronoi/atom:tags");
      for (i=0; i<nall; i++) tags[i] = atom->tag[i];

      // linked list structure for cell occupation count on the atoms
      oldnall= nall;
      memory->create(lroot,nall,"voronoi/atom:lroot"); // point to first atom index in cell (or -1 for empty cell)
      lnext = NULL; 
      lmax = 0;
  
      // build the occupation buffer
      oldnatoms = atom->natoms;
      memory->create(occvec,oldnatoms,"voronoi/atom:occvec");
#ifdef NOTINPLACE
      memory->create(sendocc,oldnatoms,"voronoi/atom:sendocc");
#endif
    }

    // get the occupation of each original voronoi cell
    checkOccupation();
  } else {
    // build cells for each output
    buildCells();
    loopCells();
  }
}

void ComputeVoronoi::buildCells()
{
  int i, j;
  const double e = 0.01;
  int nlocal = atom->nlocal;

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
  int *mask = atom->mask;
  if (radstr) {
    // check and fetch atom style variable data
    int radvar = input->variable->find(radstr);
    if (radvar < 0)
            error->all(FLERR,"Variable name for voronoi radius set does not exist");
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
    delete con_poly;
    con_poly = new container_poly(sublo_bound[0]-cut_bound[0]-e,subhi_bound[0]+cut_bound[0]+e,
                       sublo_bound[1]-cut_bound[1]-e,subhi_bound[1]+cut_bound[1]+e,
                       sublo_bound[2]-cut_bound[2]-e,subhi_bound[2]+cut_bound[2]+e,
                       int(n[0]),int(n[1]),int(n[2]),false,false,false,8); 

    // pass coordinates for local and ghost atoms to voro++
    for (i = 0; i < nall; i++)
      if( !onlyGroup || (mask[i] & groupbit) )
        con_poly->put(i,x[i][0],x[i][1],x[i][2],rfield[i]);
  } else {
    // monodisperse voro++ container
    delete con_mono;
    con_mono = new container(sublo_bound[0]-cut_bound[0]-e,subhi_bound[0]+cut_bound[0]+e,
                  sublo_bound[1]-cut_bound[1]-e,subhi_bound[1]+cut_bound[1]+e,
                  sublo_bound[2]-cut_bound[2]-e,subhi_bound[2]+cut_bound[2]+e,
                  int(n[0]),int(n[1]),int(n[2]),false,false,false,8); 

    // pass coordinates for local and ghost atoms to voro++
    for (i = 0; i < nall; i++)
      if( !onlyGroup || (mask[i] & groupbit) )
        con_mono->put(i,x[i][0],x[i][1],x[i][2]);
  }
}

void ComputeVoronoi::checkOccupation()
{
  // clear occupation vector
  memset(occvec, 0, oldnatoms*sizeof(*occvec));

  int i, j, k,
      nlocal = atom->nlocal,
      nall = atom->nghost + nlocal;
  double rx, ry, rz,
         **x = atom->x;

  // prepare destination buffer for variable evaluation
  if (nall > lmax) {
    memory->destroy(lnext);
    lmax = atom->nmax;
    memory->create(lnext,lmax,"voronoi/atom:lnext");
  }

  // clear lroot
  for (i=0; i<oldnall; ++i) lroot[i] = -1;

  // clear lnext
  for (i=0; i<nall; ++i) lnext[i] = -1;

  // loop over all local atoms and find out in which of the local first frame voronoi cells the are in
  // (need to loop over ghosts, too, to get correct occupation numbers for the second column)
  for (i=0; i<nall; ++i) {
    // again: find_voronoi_cell() should be in the common base class. Why it is not, I don't know. Ask the voro++ author.
    if ((  radstr && con_poly->find_voronoi_cell(x[i][0], x[i][1], x[i][2], rx, ry, rz, k)) ||
        ( !radstr && con_mono->find_voronoi_cell(x[i][0], x[i][1], x[i][2], rx, ry, rz, k) )) {
      // increase occupation count of this particular cell
      // only for local atoms, as we do an MPI reduce sum later
      if (i<nlocal) occvec[tags[k]-1]++;

      // add this atom to the linked list of cell j
      if (lroot[k]<0) 
        lroot[k]=i;
      else {
        j = lroot[k];
        while (lnext[j]>=0) j=lnext[j];
        lnext[j] = i; 
      }
    }
  }

  // MPI sum occupation
#ifdef NOTINPLACE
  memcpy(sendocc, occvec, oldnatoms*sizeof(*occvec));
  MPI_Allreduce(sendocc, occvec, oldnatoms, MPI_INT, MPI_SUM, world);
#else
  MPI_Allreduce(MPI_IN_PLACE, occvec, oldnatoms, MPI_INT, MPI_SUM, world);
#endif

  // determine the total number of atoms in this atom's currently occupied cell
  int c;
  for (i=0; i<oldnall; i++) { // loop over lroot (old voronoi cells)
    // count
    c = 0;
    j = lroot[i];
    while (j>=0) {
      c++;
      j = lnext[j];
    }
    // set
    j = lroot[i];
    while (j>=0) {
      voro[j][1] = c;
      j = lnext[j];
    }
  }

  // cherry pick currently owned atoms 
  for (i=0; i<nlocal; i++) {
    // set the new atom count in the atom's first frame voronoi cell
    voro[i][0] = occvec[atom->tag[i]-1];
  }
}

void ComputeVoronoi::loopCells()
{
  // invoke voro++ and fetch results for owned atoms in group
  voronoicell_neighbor c;
  int i;
  if (radstr) {
    c_loop_all cl(*con_poly);
    if (cl.start()) do if (con_poly->compute_cell(c,cl)) {
      i = cl.pid();
      processCell(c,i);
    } while (cl.inc());
  } else {
    c_loop_all cl(*con_mono);
    if (cl.start()) do if (con_mono->compute_cell(c,cl)) {
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

