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
   Contributing author (NUMA option) : Mike Brown (ORNL)
------------------------------------------------------------------------- */

#include "procmap.h"
#include "universe.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

#include <map>
#include <string>

using namespace LAMMPS_NS;

enum{MULTIPLE};                   // same as in Comm

/* ---------------------------------------------------------------------- */

ProcMap::ProcMap(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   create a one-level 3d grid of procs via procs2box()
------------------------------------------------------------------------- */

int ProcMap::onelevel_grid(int nprocs, int *user_procgrid, int *procgrid,
			   int otherflag, int other_style_caller,
			   int *other_procgrid_caller)
{
  other_style = other_style_caller;
  other_procgrid[0] = other_procgrid_caller[0];
  other_procgrid[1] = other_procgrid_caller[1];
  other_procgrid[2] = other_procgrid_caller[2];

  int flag = procs2box(nprocs,user_procgrid,procgrid,1,1,1,otherflag);
  return flag;
}

/* ----------------------------------------------------------------------
   create a two-level 3d grid of procs and cores via procs2box()
------------------------------------------------------------------------- */

int ProcMap::twolevel_grid(int nprocs, int *user_procgrid, int *procgrid,
			   int ncores, int *user_coregrid, int *coregrid,
			   int otherflag, int other_style_caller,
			   int *other_procgrid_caller)
{
  error->all(FLERR,
	     "The twolevel option is not yet supported, but will be soon");
  return 1;
}

/* ----------------------------------------------------------------------
   create a 3d grid of procs that does a 2-level hierarchy within a node
   auto-detects NUMA sockets within a multi-core node
   return 1 if successful, 0 if not
------------------------------------------------------------------------- */

int ProcMap::numa_grid(int nprocs, int *user_procgrid, int *procgrid,
		       int *numagrid)
{
  // hardwire this for now

  int numa_nodes = 1;

  // get names of all nodes

  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char node_names[MPI_MAX_PROCESSOR_NAME*nprocs];
  MPI_Get_processor_name(node_name,&name_length);
  MPI_Allgather(&node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,&node_names,
                MPI_MAX_PROCESSOR_NAME,MPI_CHAR,world);
  std::string node_string = std::string(node_name);
  
  // get number of procs per node

  std::map<std::string,int> name_map;
  std::map<std::string,int>::iterator np;
  for (int i = 0; i < nprocs; i++) {
    std::string i_string = std::string(&node_names[i*MPI_MAX_PROCESSOR_NAME]);
    np = name_map.find(i_string);
    if (np == name_map.end()) name_map[i_string] = 1;
    else np->second++;
  }
  procs_per_node = name_map.begin()->second;
  procs_per_numa = procs_per_node / numa_nodes;
  
  // error return if any of these conditions met
  
  if (procs_per_numa < 4 ||            // less than 4 procs per numa node
      procs_per_node % numa_nodes ||   // no-op since numa_nodes = 1 for now
      nprocs % procs_per_numa ||       // total procs not a multiple of node
      nprocs == procs_per_numa ||      // only 1 node used
      user_procgrid[0] > 1 ||          // user specified grid > 1 in any dim
      user_procgrid[1] > 1 ||
      user_procgrid[2] > 1)
    return 0;
  
  // user settings for the factorization per numa node
  // currently not user settable

  int user_numagrid[3];
  user_numagrid[0] = user_numagrid[1] = user_numagrid[2] = 0;
  
  // if user specifies 1 for a proc grid dimension,
  // also use 1 for the numa grid dimension

  if (user_procgrid[0] == 1) user_numagrid[0] = 1;
  if (user_procgrid[1] == 1) user_numagrid[1] = 1;
  if (user_procgrid[2] == 1) user_numagrid[2] = 1;
  
  // initial factorization within NUMA node

  procs2box(procs_per_numa,user_numagrid,numagrid,1,1,1,0);
  if (numagrid[0]*numagrid[1]*numagrid[2] != procs_per_numa)
    error->all(FLERR,"Bad grid of processors");
  
  // factorization for the grid of NUMA nodes

  int node_count = nprocs / procs_per_numa;
  procs2box(node_count,user_procgrid,nodegrid,
	    numagrid[0],numagrid[1],numagrid[2],0);
  if (procgrid[0]*procgrid[1]*procgrid[2] != node_count)
    error->all(FLERR,"Bad grid of processors");
  
  // repeat NUMA node factorization using subdomain sizes
  // refines the factorization if the user specified the node layout

  procs2box(procs_per_numa,user_numagrid,numagrid,
	    procgrid[0],procgrid[1],procgrid[2],0);

  // assign a unique id to each node

  node_id = 0;
  int node_num = 0;
  for (np = name_map.begin(); np != name_map.end(); ++np) {
    if (np->first == node_string) node_id = node_num;
    node_num++;
  }

  // return the proc-level factorization

  procgrid[0] = nodegrid[0] * numagrid[0];
  procgrid[1] = nodegrid[1] * numagrid[1];
  procgrid[2] = nodegrid[2] * numagrid[2];

  return 1;
}

/* ----------------------------------------------------------------------
   create a 1-level 3d grid of procs via procs2box()
------------------------------------------------------------------------- */

void ProcMap::custom_grid(char *cfile, int nprocs,
			  int *user_procgrid, int *procgrid)
{
  error->all(FLERR,
	     "The grid custom option is not yet supported, but will be soon");
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d box so as to minimize surface area
   area = surface area of each of 3 faces of simulation box divided by sx,sy,sz
   for triclinic, area = cross product of 2 edge vectors stored in h matrix
   valid assignment will be factorization of nprocs = Px by Py by Pz
   user_factors = if non-zero, factors are specified by user
   sx,sy,sz = scale box xyz dimension vy dividing by sx,sy,sz
   other = 1 to enforce compatability with other partition's layout
   return factors = # of procs assigned to each dimension
   return 1 if factor successfully, 0 if not
------------------------------------------------------------------------- */

int ProcMap::procs2box(int nprocs, int *user_factors, int *factors, 
		       const int sx, const int sy, const int sz, int other)
{
  factors[0] = user_factors[0];
  factors[1] = user_factors[1];
  factors[2] = user_factors[2];

  // all 3 proc counts are specified

  if (factors[0] && factors[1] && factors[2]) return 1;

  // 2 out of 3 proc counts are specified

  if (factors[0] > 0 && factors[1] > 0) {
    factors[2] = nprocs/(factors[0]*factors[1]);
    return 1;
  } else if (factors[0] > 0 && factors[2] > 0) {
    factors[1] = nprocs/(factors[0]*factors[2]);
    return 1;
  } else if (factors[1] > 0 && factors[2] > 0) {
    factors[0] = nprocs/(factors[1]*factors[2]);
    return 1;
  } 

  // determine cross-sectional areas for orthogonal and triclinic boxes
  // area[0] = xy, area[1] = xz, area[2] = yz

  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx * sy);
    area[1] = domain->xprd * domain->zprd / (sx * sz);
    area[2] = domain->yprd * domain->zprd / (sy * sz);
  } else {
    double *h = domain->h;
    double a[3],b[3],c[3];
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[5]; b[1] = h[1]; b[2] = 0.0;
    MathExtra::cross3(a,b,c);
    area[0] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx * sy);
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[1] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx * sz);
    a[0] = h[5]; a[1] = h[1]; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[2] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sy * sz);
  }

  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of nprocs
  // only consider valid cases that match procgrid settings
  // surf = surface area of a proc sub-domain
  // only consider cases that match user_factors & other_procgrid settings
  // success = 1 if valid factoriztion is found
  // may not be if other constraint is enforced

  int ipx,ipy,ipz,valid;
  double surf;
  
  int success = 0;
  ipx = 1;
  while (ipx <= nprocs) {
    valid = 1;
    if (user_factors[0] && ipx != user_factors[0]) valid = 0;
    if (other) {
      if (other_style == MULTIPLE && other_procgrid[0] % ipx) valid = 0;
    }
    if (nprocs % ipx) valid = 0;

    if (!valid) {
      ipx++;
      continue;
    }

    ipy = 1;
    while (ipy <= nprocs/ipx) {
      valid = 1;
      if (user_factors[1] && ipy != user_factors[1]) valid = 0;
      if (other) {
	if (other_style == MULTIPLE && other_procgrid[1] % ipy) valid = 0;
      }
      if ((nprocs/ipx) % ipy) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      ipz = nprocs/ipx/ipy;
      valid = 1;
      if (user_factors[2] && ipz != user_factors[2]) valid = 0;
      if (other) {
	if (other_style == MULTIPLE && other_procgrid[2] % ipz) valid = 0;
      }
      if (domain->dimension == 2 && ipz != 1) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
      if (surf < bestsurf) {
	success = 1;
	bestsurf = surf;
	factors[0] = ipx;
	factors[1] = ipy;
	factors[2] = ipz;
      }
      ipy++;
    }

    ipx++;
  }

  return success;
}

/* ----------------------------------------------------------------------
   map processors to 3d grid via MPI_Cart routines
   MPI may do layout in machine-optimized fashion
------------------------------------------------------------------------- */

void ProcMap::cart_map(int reorder, int *procgrid,
		       int *myloc, int procneigh[3][2], int ***grid2proc)
{
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;
      
  MPI_Cart_create(world,3,procgrid,periods,reorder,&cartesian);
  MPI_Cart_get(cartesian,3,procgrid,periods,myloc);
  MPI_Cart_shift(cartesian,0,1,&procneigh[0][0],&procneigh[0][1]);
  MPI_Cart_shift(cartesian,1,1,&procneigh[1][0],&procneigh[1][1]);
  MPI_Cart_shift(cartesian,2,1,&procneigh[2][0],&procneigh[2][1]);

  int coords[3];
  int i,j,k;
  for (i = 0; i < procgrid[0]; i++)
    for (j = 0; j < procgrid[1]; j++)
      for (k = 0; k < procgrid[2]; k++) {
	coords[0] = i; coords[1] = j; coords[2] = k;
	MPI_Cart_rank(cartesian,coords,&grid2proc[i][j][k]);
      }
  
  MPI_Comm_free(&cartesian);
}

/* ----------------------------------------------------------------------
   map processors to 3d grid via MPI_Cart routines
   respect sub-grid of cores within each node
   MPI may do layout in machine-optimized fashion
------------------------------------------------------------------------- */

void ProcMap::cart_map(int reorder, int *procgrid, int *coregrid,
		       int *myloc, int procneigh[3][2], int ***grid2proc)
{
}

/* ----------------------------------------------------------------------
   map processors to 3d grid in XYZ order
------------------------------------------------------------------------- */

void ProcMap::xyz_map(char *xyz, int *procgrid,
		      int *myloc, int procneigh[3][2], int ***grid2proc)
{
  int me;
  MPI_Comm_rank(world,&me);

  int i,j,k;
  for (i = 0; i < procgrid[0]; i++)
    for (j = 0; j < procgrid[1]; j++)
      for (k = 0; k < procgrid[2]; k++) {
	grid2proc[i][j][k] = k*procgrid[1]*procgrid[0] + j*procgrid[0] + i;
	if (xyz[0] == 'x' && xyz[1] == 'y' && xyz[2] == 'z')
	  grid2proc[i][j][k] = k*procgrid[1]*procgrid[0] + j*procgrid[0] + i;
	else if (xyz[0] == 'x' && xyz[1] == 'z' && xyz[2] == 'y')
	  grid2proc[i][j][k] = j*procgrid[2]*procgrid[0] + k*procgrid[0] + i;
	else if (xyz[0] == 'y' && xyz[1] == 'x' && xyz[2] == 'z')
	  grid2proc[i][j][k] = k*procgrid[0]*procgrid[1] + i*procgrid[1] + j;
	else if (xyz[0] == 'y' && xyz[1] == 'z' && xyz[2] == 'x')
	  grid2proc[i][j][k] = i*procgrid[2]*procgrid[1] + k*procgrid[1] + j;
	else if (xyz[0] == 'z' && xyz[1] == 'x' && xyz[2] == 'y')
	  grid2proc[i][j][k] = j*procgrid[0]*procgrid[2] + i*procgrid[2] + k;
	else if (xyz[0] == 'z' && xyz[1] == 'y' && xyz[2] == 'x')
	  grid2proc[i][j][k] = i*procgrid[1]*procgrid[2] + j*procgrid[2] + k;

	if (grid2proc[i][j][k] == me) {
	  myloc[0] = i; myloc[1] = j, myloc[2] = k;
	}
      }

  int minus,plus;
  grid_shift(myloc[0],procgrid[0],minus,plus);
  procneigh[0][0] = grid2proc[minus][myloc[1]][myloc[2]];
  procneigh[0][1] = grid2proc[plus][myloc[1]][myloc[2]];

  grid_shift(myloc[1],procgrid[1],minus,plus);
  procneigh[1][0] = grid2proc[myloc[0]][minus][myloc[2]];
  procneigh[1][1] = grid2proc[myloc[0]][plus][myloc[2]];

  grid_shift(myloc[2],procgrid[2],minus,plus);
  procneigh[2][0] = grid2proc[myloc[0]][myloc[1]][minus];
  procneigh[2][1] = grid2proc[myloc[0]][myloc[1]][plus];
}

/* ----------------------------------------------------------------------
   map processors to 3d grid in XYZ order
   respect sub-grid of cores within each node
------------------------------------------------------------------------- */

void ProcMap::xyz_map(char *xyz, int *procgrid, int *coregrid,
		      int *myloc, int procneigh[3][2], int ***grid2proc)
{
}

/* ----------------------------------------------------------------------
   map processors to 3d grid in 2-level NUMA ordering
------------------------------------------------------------------------- */

void ProcMap::numa_map(int *numagrid,
		       int *myloc, int procneigh[3][2], int ***grid2proc)
{
  // setup a per node communicator and find rank within

  MPI_Comm node_comm;
  MPI_Comm_split(world,node_id,0,&node_comm);  
  int node_rank;
  MPI_Comm_rank(node_comm,&node_rank);
  
  // setup a per numa communicator and find rank within

  MPI_Comm numa_comm;
  int local_numa = node_rank / procs_per_numa;
  MPI_Comm_split(node_comm,local_numa,0,&numa_comm);     
  int numa_rank;
  MPI_Comm_rank(numa_comm,&numa_rank);
  
  // setup a communicator with the rank 0 procs from each numa node

  MPI_Comm numa_leaders;
  MPI_Comm_split(world,numa_rank,0,&numa_leaders);
  
  // use the MPI Cartesian routines to map the nodes to the grid
  // could implement xyz mapflag as in non-NUMA case?

  int reorder = 0;
  int periods[3];
  periods[0] = periods[1] = periods[2] = 1;
  MPI_Comm cartesian;
  if (numa_rank == 0) {
    MPI_Cart_create(numa_leaders,3,nodegrid,periods,reorder,&cartesian);
    MPI_Cart_get(cartesian,3,nodegrid,periods,myloc);
  }
  
  // broadcast numa node location in grid to other procs in numa node

  MPI_Bcast(myloc,3,MPI_INT,0,numa_comm);
  
  // compute my location within the node grid

  int z_offset = numa_rank / (numagrid[0] * numagrid[1]);
  int y_offset = (numa_rank % (numagrid[0] * numagrid[1]))/numagrid[0];
  int x_offset = numa_rank % numagrid[0];
  myloc[0] = myloc[0] * numagrid[0] + x_offset;
  myloc[1] = myloc[1] * numagrid[1] + y_offset;
  myloc[2] = myloc[2] * numagrid[2] + z_offset;
  
  // allgather of locations to fill grid2proc

  int nprocs;
  MPI_Comm_size(world,&nprocs);

  int **gridi;
  memory->create(gridi,nprocs,3,"comm:gridi");
  MPI_Allgather(&myloc,3,MPI_INT,gridi[0],3,MPI_INT,world);
  for (int i = 0; i < nprocs; i++)
    grid2proc[gridi[i][0]][gridi[i][1]][gridi[i][2]] = i;
  memory->destroy(gridi);
  
  // proc IDs of neighbors

  int minus,plus;
  grid_shift(myloc[0],nodegrid[0]*numagrid[0],minus,plus);
  procneigh[0][0] = grid2proc[minus][myloc[1]][myloc[2]];
  procneigh[0][1] = grid2proc[plus][myloc[1]][myloc[2]];

  grid_shift(myloc[1],nodegrid[1]*numagrid[1],minus,plus);
  procneigh[1][0] = grid2proc[myloc[0]][minus][myloc[2]];
  procneigh[1][1] = grid2proc[myloc[0]][plus][myloc[2]];

  grid_shift(myloc[2],nodegrid[2]*numagrid[2],minus,plus);
  procneigh[2][0] = grid2proc[myloc[0]][myloc[1]][minus];
  procneigh[2][1] = grid2proc[myloc[0]][myloc[1]][plus];

  // clean-up

  if (numa_rank == 0) MPI_Comm_free(&cartesian);
  MPI_Comm_free(&numa_leaders);
  MPI_Comm_free(&numa_comm);
  MPI_Comm_free(&node_comm);
}

/* ----------------------------------------------------------------------
   map processors to 3d grid in custom ordering
------------------------------------------------------------------------- */

void ProcMap::custom_map(int *myloc, int procneigh[3][2], int ***grid2proc)
{
}

/* ----------------------------------------------------------------------
   minus,plus = indices of neighboring processors in a dimension
------------------------------------------------------------------------- */

void ProcMap::grid_shift(int myloc, int nprocs, int &minus, int &plus) 
{
  minus = myloc - 1;
  if (minus < 0) minus = nprocs - 1;
  plus = myloc + 1;
  if (plus == nprocs) plus = 0;
}

/* ----------------------------------------------------------------------
   output mapping of processors to 3d grid to file
------------------------------------------------------------------------- */

void ProcMap::output(char *file, int *procgrid, int ***grid2proc)
{
  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  FILE *fp;
  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) error->one(FLERR,"Cannot open processors custom file");
    fprintf(fp,"LAMMPS mapping of processors to 3d grid\n");
    fprintf(fp,"partition = %d\n",universe->iworld+1);
    fprintf(fp,"Px Py Pz = %d %d %d\n",procgrid[0],procgrid[1],procgrid[2]);
    fprintf(fp,"world-ID universe-ID original-ID: I J K: name\n\n");
  }

  // find me in the grid

  int ime,jme,kme;
  for (int i = 0; i < procgrid[0]; i++)
    for (int j = 0; j < procgrid[1]; j++)
      for (int k = 0; k < procgrid[2]; k++)
	if (grid2proc[i][j][k] == me) {
	  ime = i; jme = j; kme = k;
	}

  // polled comm of grid mapping info from each proc to proc 0

  int tmp;
  int vec[6];
  char procname[MPI_MAX_PROCESSOR_NAME+1];
  MPI_Status status;

  vec[0] = me;
  vec[1] = universe->me;
  MPI_Comm_rank(universe->uorig,&vec[2]);
  vec[3] = ime + 1;
  vec[4] = jme + 1;
  vec[5] = kme + 1;

  int len;
  MPI_Get_processor_name(procname,&len);
  procname[len] = '\0';

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Recv(vec,6,MPI_INT,iproc,0,world,&status);
	MPI_Recv(procname,MPI_MAX_PROCESSOR_NAME+1,MPI_CHAR,
		 iproc,0,world,&status);
      }

      fprintf(fp,"%d %d %d: %d %d %d: %s\n",
	      vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],procname);
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Send(vec,6,MPI_INT,0,0,world);
    MPI_Send(procname,strlen(procname)+1,MPI_CHAR,0,0,world);
  }

  // close output file

  if (me == 0) fclose(fp);
}
