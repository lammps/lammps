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

#define MAXLINE 128

enum{MULTIPLE};                   // same as in Comm

/* ---------------------------------------------------------------------- */

ProcMap::ProcMap(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   create a one-level 3d grid of procs
------------------------------------------------------------------------- */

void ProcMap::onelevel_grid(int nprocs, int *user_procgrid, int *procgrid,
			    int otherflag, int other_style,
			    int *other_procgrid, int *other_coregrid)
{
  int **factors;

  // factors = list of all possible 3 factors of processor count

  int npossible = factor(nprocs,NULL);
  memory->create(factors,npossible,3,"procmap:factors");
  npossible = factor(nprocs,factors);

  // constrain by 2d, user request, other partition

  if (domain->dimension == 2) npossible = cull_2d(npossible,factors,3);
  npossible = cull_user(npossible,factors,3,user_procgrid);
  if (otherflag) npossible = cull_other(npossible,factors,3,
					other_style,other_procgrid,
					other_coregrid);

  // user/other constraints make failure possible

  if (npossible == 0)
    error->all(FLERR,"Could not create 3d grid of processors");

  // select best set of 3 factors based on surface area of proc sub-domains

  best_factors(npossible,factors,procgrid,1,1,1);

  // clean-up

  memory->destroy(factors);
}

/* ----------------------------------------------------------------------
   create a two-level 3d grid of procs
------------------------------------------------------------------------- */

void ProcMap::twolevel_grid(int nprocs, int *user_procgrid, int *procgrid,
			    int ncores, int *user_coregrid, int *coregrid,
			    int otherflag, int other_style, 
			    int *other_procgrid, int *other_coregrid)
{
  int **nfactors,**cfactors,**factors;

  if (nprocs % ncores) 
    error->all(FLERR,"Processors twogrid requires proc count "
	       "be a multiple of core count");

  // nfactors = list of all possible 3 factors of node count
  // constrain by 2d

  int nnpossible = factor(nprocs/ncores,NULL);
  memory->create(nfactors,nnpossible,3,"procmap:nfactors");
  nnpossible = factor(nprocs/ncores,nfactors);

  if (domain->dimension == 2) nnpossible = cull_2d(nnpossible,nfactors,3);

  // cfactors = list of all possible 3 factors of core count
  // constrain by 2d

  int ncpossible = factor(ncores,NULL);
  memory->create(cfactors,ncpossible,3,"procmap:cfactors");
  ncpossible = factor(ncores,cfactors);

  if (domain->dimension == 2) ncpossible = cull_2d(ncpossible,cfactors,3);
  ncpossible = cull_user(ncpossible,cfactors,3,user_coregrid);

  // factors = all combinations of nfactors and cfactors
  // factors stores additional index pointing to corresponding cfactors
  // constrain by user request, other partition

  int npossible = nnpossible * ncpossible;
  memory->create(factors,npossible,4,"procmap:factors");
  npossible = combine_factors(nnpossible,nfactors,ncpossible,cfactors,factors);

  npossible = cull_user(npossible,factors,4,user_procgrid);
  if (otherflag) npossible = cull_other(npossible,factors,4,
					other_style,other_procgrid,
					other_coregrid);

  // user/other constraints make failure possible

  if (npossible == 0)
    error->all(FLERR,"Could not create twolevel 3d grid of processors");

  // select best set of 3 factors based on surface area of proc sub-domains
  // index points to corresponding core factorization

  int index = best_factors(npossible,factors,procgrid,1,1,1);

  coregrid[0] = cfactors[factors[index][3]][0];
  coregrid[1] = cfactors[factors[index][3]][1];
  coregrid[2] = cfactors[factors[index][3]][2];

  // clean-up

  memory->destroy(nfactors);
  memory->destroy(cfactors);
  memory->destroy(factors);
}

/* ----------------------------------------------------------------------
   create a 3d grid of procs that does a 2-level hierarchy within a node
   auto-detects NUMA sockets within a multi-core node
------------------------------------------------------------------------- */

void ProcMap::numa_grid(int nprocs, int *user_procgrid, int *procgrid,
			int *numagrid)
{
  // hardwire this for now

  int numa_nodes = 1;

  // get names of all nodes

  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(node_name,&name_length);
  node_name[name_length] = '\0';
  char *node_names = new char[MPI_MAX_PROCESSOR_NAME*nprocs];
  MPI_Allgather(node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,node_names,
                MPI_MAX_PROCESSOR_NAME,MPI_CHAR,world);
  std::string node_string = std::string(node_name);
  
  // get number of procs per node
  // NOTE: could do this without STL map

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
  
  delete [] node_names;

  // error if any of these conditions met
  
  if (nprocs % procs_per_numa ||       // total procs not a multiple of node
      user_procgrid[0] > 1 ||          // user specified grid > 1 in any dim
      user_procgrid[1] > 1 ||
      user_procgrid[2] > 1)
    error->all(FLERR,"Could not create numa grid of processors");
  
  // user settings for the factorization per numa node
  // currently not user settable
  // if user specifies 1 for a proc grid dimension,
  // also use 1 for the numa grid dimension

  int user_numagrid[3];
  user_numagrid[0] = user_numagrid[1] = user_numagrid[2] = 0;

  if (user_procgrid[0] == 1) user_numagrid[0] = 1;
  if (user_procgrid[1] == 1) user_numagrid[1] = 1;
  if (user_procgrid[2] == 1) user_numagrid[2] = 1;
  
  // initial factorization within NUMA node

  int **numafactors;
  int numapossible = factor(procs_per_numa,NULL);
  memory->create(numafactors,numapossible,3,"procmap:numafactors");
  numapossible = factor(procs_per_numa,numafactors);

  if (domain->dimension == 2) 
    numapossible = cull_2d(numapossible,numafactors,3);
  numapossible = cull_user(numapossible,numafactors,3,user_numagrid);

  if (numapossible == 0)
    error->all(FLERR,"Could not create numa grid of processors");

  best_factors(numapossible,numafactors,numagrid,1,1,1);

  // user_nodegrid = implied user contraints on nodes

  int user_nodegrid[3];
  user_nodegrid[0] = user_procgrid[0] / numagrid[0];
  user_nodegrid[1] = user_procgrid[1] / numagrid[1];
  user_nodegrid[2] = user_procgrid[2] / numagrid[2];

  // factorization for the grid of NUMA nodes

  int node_count = nprocs / procs_per_numa;

  int **nodefactors;
  int nodepossible = factor(node_count,NULL);
  memory->create(nodefactors,nodepossible,3,"procmap:nodefactors");
  nodepossible = factor(node_count,nodefactors);

  if (domain->dimension == 2) 
    nodepossible = cull_2d(nodepossible,nodefactors,3);
  nodepossible = cull_user(nodepossible,nodefactors,3,user_nodegrid);

  if (nodepossible == 0)
    error->all(FLERR,"Could not create numa grid of processors");

  best_factors(nodepossible,nodefactors,nodegrid,
	       numagrid[0],numagrid[1],numagrid[2]);
  
  // repeat NUMA node factorization using subdomain sizes
  // refines the factorization if the user specified the node layout
  // NOTE: this will not re-enforce user-procgrid constraint will it?

  best_factors(numapossible,numafactors,numagrid,
	       nodegrid[0],nodegrid[1],nodegrid[2]);

  memory->destroy(numafactors);
  memory->destroy(nodefactors);

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
}

/* ----------------------------------------------------------------------
   define a 3d grid from a custom file
------------------------------------------------------------------------- */

void ProcMap::custom_grid(char *cfile, int nprocs,
			  int *user_procgrid, int *procgrid)
{
  FILE *fp;
  char line[MAXLINE];

  int me;
  MPI_Comm_rank(world,&me);

  if (me == 0) {
    FILE *fp = fopen(cfile,"r");
    if (fp == NULL) error->one(FLERR,"Cannot open custom file");

    // skip header = blank and comment lines
    
    char *ptr;
    if (!fgets(line,MAXLINE,fp))
      error->one(FLERR,"Unexpected end of custom file");
    while (1) {
      if (ptr = strchr(line,'#')) *ptr = '\0';
      if (strspn(line," \t\n\r") != strlen(line)) break;
      if (!fgets(line,MAXLINE,fp))
	error->one(FLERR,"Unexpected end of custom file");
    }
  }

  int n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d %d %d",&procgrid[0],&procgrid[1],&procgrid[2]);

  int flag = 0;
  if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs) flag = 1;
  if (user_procgrid[0] && procgrid[0] != user_procgrid[0]) flag = 1;
  if (user_procgrid[1] && procgrid[1] != user_procgrid[1]) flag = 1;
  if (user_procgrid[2] && procgrid[2] != user_procgrid[2]) flag = 1;
  if (flag) error->all(FLERR,"Processors custom grid file is inconsistent");

  // cmap = map of procs to grid
  // store for use in custom_map()

  memory->create(cmap,nprocs,4,"procmap:cmap");
  for (int i = 0; i < nprocs; i++) cmap[i][0] = -1;

  if (me == 0) {
    for (int i = 0; i < nprocs; i++) {
      if (!fgets(line,MAXLINE,fp))
	error->one(FLERR,"Unexpected end of custom file");
      sscanf(line,"%d %d %d %d",
	     &cmap[i][0],&cmap[i][1],&cmap[i][2],&cmap[i][3]);
    }
    fclose(fp);
  }

  MPI_Bcast(&cmap[0][0],nprocs*4,MPI_INT,0,world);

  // error check on cmap values

  flag = 0;
  for (int i = 0; i < nprocs; i++) {
    if (cmap[i][0] == -1) flag = 1;
    else {
      if (cmap[i][1] <= 0 || cmap[i][1] > procgrid[0]) flag = 1;
      if (cmap[i][2] <= 0 || cmap[i][2] > procgrid[1]) flag = 1;
      if (cmap[i][3] <= 0 || cmap[i][3] > procgrid[2]) flag = 1;
    }
  }
  if (flag) error->all(FLERR,"Processors custom grid file is inconsistent");
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

void ProcMap::cart_map(int reorder, int *procgrid, int ncores, int *coregrid,
		       int *myloc, int procneigh[3][2], int ***grid2proc)
{
  // setup NUMA params that numa_grid() sets up

  int me;
  MPI_Comm_rank(world,&me);

  procs_per_node = ncores;
  procs_per_numa = ncores;
  node_id = me/ncores;
  nodegrid[0] = procgrid[0] / coregrid[0];
  nodegrid[1] = procgrid[1] / coregrid[1];
  nodegrid[2] = procgrid[2] / coregrid[2];

  // now can use numa_map() to perform mapping

  numa_map(reorder,coregrid,myloc,procneigh,grid2proc);
}

/* ----------------------------------------------------------------------
   map processors to 3d grid in XYZ order
   called by onelevel
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

  // proc IDs of neighbors

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
   called by twolevel
------------------------------------------------------------------------- */

void ProcMap::xyz_map(char *xyz, int *procgrid, int ncores, int *coregrid,
		      int *myloc, int procneigh[3][2], int ***grid2proc)
{
  int me;
  MPI_Comm_rank(world,&me);

  nodegrid[0] = procgrid[0] / coregrid[0];
  nodegrid[1] = procgrid[1] / coregrid[1];
  nodegrid[2] = procgrid[2] / coregrid[2];

  int i,j,k,inode,jnode,knode,icore,jcore,kcore;
  for (i = 0; i < procgrid[0]; i++)
    for (j = 0; j < procgrid[1]; j++)
      for (k = 0; k < procgrid[2]; k++) {
	inode = i/coregrid[0];
	jnode = j/coregrid[1];
	knode = k/coregrid[2];
	icore = i % coregrid[0];
	jcore = j % coregrid[1];
	kcore = k % coregrid[2];

	if (xyz[0] == 'x' && xyz[1] == 'y' && xyz[2] == 'z') {
	  grid2proc[i][j][k] = ncores * 
	    (knode*nodegrid[1]*nodegrid[0] + jnode*nodegrid[0] + inode) +
	    (kcore*coregrid[1]*coregrid[0] + jcore*coregrid[0] + icore);
	} else if (xyz[0] == 'x' && xyz[1] == 'z' && xyz[2] == 'y')
	  grid2proc[i][j][k] = ncores * 
	    (jnode*nodegrid[2]*nodegrid[0] + knode*nodegrid[0] + inode) + 
	    (jcore*coregrid[2]*coregrid[0] + kcore*coregrid[0] + icore);
	else if (xyz[0] == 'y' && xyz[1] == 'x' && xyz[2] == 'z')
	  grid2proc[i][j][k] = ncores *
	    (knode*nodegrid[0]*nodegrid[1] + inode*nodegrid[1] + jnode) +
	    (kcore*coregrid[0]*coregrid[1] + icore*coregrid[1] + jcore);
	else if (xyz[0] == 'y' && xyz[1] == 'z' && xyz[2] == 'x')
	  grid2proc[i][j][k] = ncores *
	    (inode*nodegrid[2]*nodegrid[1] + knode*nodegrid[1] + jnode) +
	    (icore*coregrid[2]*coregrid[1] + kcore*coregrid[1] + jcore);
	else if (xyz[0] == 'z' && xyz[1] == 'x' && xyz[2] == 'y')
	  grid2proc[i][j][k] = ncores *
	    (jnode*nodegrid[0]*nodegrid[2] + inode*nodegrid[2] + knode) +
	    (jcore*coregrid[0]*coregrid[2] + icore*coregrid[2] + kcore);
	else if (xyz[0] == 'z' && xyz[1] == 'y' && xyz[2] == 'x')
	  grid2proc[i][j][k] = ncores *
	    (inode*nodegrid[1]*nodegrid[2] + jnode*nodegrid[2] + knode) +
	    (icore*coregrid[1]*coregrid[2] + jcore*coregrid[2] + kcore);

	if (grid2proc[i][j][k] == me) {
	  myloc[0] = i; myloc[1] = j, myloc[2] = k;
	}
      }

  // proc IDs of neighbors

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
   map processors to 3d grid in 2-level NUMA ordering
------------------------------------------------------------------------- */

void ProcMap::numa_map(int reorder, int *numagrid,
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
  
  // allgather of myloc into gridi to fill grid2proc

  int nprocs;
  MPI_Comm_size(world,&nprocs);

  int **gridi;
  memory->create(gridi,nprocs,3,"comm:gridi");
  MPI_Allgather(myloc,3,MPI_INT,gridi[0],3,MPI_INT,world);
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

void ProcMap::custom_map(int *procgrid,
			 int *myloc, int procneigh[3][2], int ***grid2proc)
{
  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  for (int i = 0; i < nprocs; i++) {
    grid2proc[cmap[i][1]-1][cmap[i][2]-1][cmap[i][3]-1] = cmap[i][0];
    if (cmap[i][0] == me) {
      myloc[0] = cmap[i][1] - 1;
      myloc[1] = cmap[i][2] - 1;
      myloc[2] = cmap[i][3] - 1;
    }
  }

  // proc IDs of neighbors

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

  memory->destroy(cmap);
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
    if (fp == NULL) error->one(FLERR,"Cannot open processors output file");
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

/* ----------------------------------------------------------------------
   generate all possible 3-integer factorizations of N
   store them in factors if non-NULL
   return # of factorizations
------------------------------------------------------------------------- */

int ProcMap::factor(int n, int **factors)
{
  int i,j,nyz;

  int m = 0;
  for (i = 1; i <= n; i++) {
    if (n % i) continue;
    nyz = n/i;
    for (j = 1; j <= nyz; j++) {
      if (nyz % j) continue;
      if (factors) {
	factors[m][0] = i;
	factors[m][1] = j;
	factors[m][2] = nyz/j;
      }
      m++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   create N1*N2 new factors (procs) from factors1 (nodes) and factors2 (cores)
   store index of corresponding core factors in factors[][3]
------------------------------------------------------------------------- */

int ProcMap::combine_factors(int n1, int **factors1, int n2, int **factors2,
			     int **factors)
{
  int m = 0;
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      factors[m][0] = factors1[i][0]*factors2[j][0];
      factors[m][1] = factors1[i][1]*factors2[j][1];
      factors[m][2] = factors1[i][2]*factors2[j][2];
      factors[m][3] = j;
      m++;
    }
  return n1*n2;
}

/* ----------------------------------------------------------------------
   remove any factors where Pz != 1 for 2d
------------------------------------------------------------------------- */

int ProcMap::cull_2d(int n, int **factors, int m)
{
  int i = 0;
  while (i < n) {
    if (factors[i][2] != 1) {
      for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
      n--;
    } else i++;
  }
  return n;
}

/* ----------------------------------------------------------------------
   remove any factors that do not match non-zero user_factors Px,Py,Pz
------------------------------------------------------------------------- */

int ProcMap::cull_user(int n, int **factors, int m, int *user_factors)
{
  int i = 0;
  while (i < n) {
    int flag = 0;
    if (user_factors[0] && factors[i][0] != user_factors[0]) flag = 1;
    if (user_factors[1] && factors[i][1] != user_factors[1]) flag = 1;
    if (user_factors[2] && factors[i][2] != user_factors[2]) flag = 1;
    if (flag) {
      for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
      n--;
    } else i++;
  }
  return n;
}

/* ----------------------------------------------------------------------
   remove any factors that do not match settings from other partition
   MULTIPLE = other Nx,Ny,Nz must be multiple of my Px,Py,Pz
              where Nx,Ny,Nz = node grid = procgrid/coregrid
------------------------------------------------------------------------- */

int ProcMap::cull_other(int n, int **factors, int m, 
			int other_style, int *other_procgrid,
			int *other_coregrid)
{
  int i = 0;
  while (i < n) {
    if (other_style == MULTIPLE) {
      int flag = 0;
      if ((other_procgrid[0]/other_coregrid[0]) % factors[i][0]) flag = 1;
      if ((other_procgrid[1]/other_coregrid[1]) % factors[i][1]) flag = 1;
      if ((other_procgrid[2]/other_coregrid[2]) % factors[i][2]) flag = 1;
      if (flag) {
	for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
	n--;
      } else i++;
    }
  }
  return n;
}

/* ----------------------------------------------------------------------
   choose best factors from list of Npossible factors
   best = minimal surface area of sub-domain
   return best = 3 factors
   return index of best factors in factors
------------------------------------------------------------------------- */

int ProcMap::best_factors(int npossible, int **factors, int *best,
			  const int sx, const int sy, const int sz)
{
  // determine cross-sectional areas for orthogonal and triclinic boxes
  // for triclinic, area = cross product of 2 edge vectors stored in h matrix
  // area[3] = surface area 3 box faces divided by sx,sy,sz
  // area[0] = xy, area[1] = xz, area[2] = yz
  
  double area[3];
  if (domain->triclinic == 0) {
    area[0] = domain->xprd * domain->yprd / (sx*sy);
    area[1] = domain->xprd * domain->zprd / (sx*sz);
    area[2] = domain->yprd * domain->zprd / (sy*sz);
  } else {
    double *h = domain->h;
    double a[3],b[3],c[3];
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[5]; b[1] = h[1]; b[2] = 0.0;
    MathExtra::cross3(a,b,c);
    area[0] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx*sy);
    a[0] = h[0]; a[1] = 0.0; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[1] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sx*sz);
    a[0] = h[5]; a[1] = h[1]; a[2] = 0.0;
    b[0] = h[4]; b[1] = h[3]; b[2] = h[2];
    MathExtra::cross3(a,b,c);
    area[2] = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]) / (sy*sz);
  }

  int index;
  double surf;
  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  for (int m = 0; m < npossible; m++) {
    surf = area[0]/factors[m][0]/factors[m][1] + 
      area[1]/factors[m][0]/factors[m][2] + 
      area[2]/factors[m][1]/factors[m][2];
    if (surf < bestsurf) {
      bestsurf = surf;
      best[0] = factors[m][0];
      best[1] = factors[m][1];
      best[2] = factors[m][2];
      index = m;
    }
  }

  return index;
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

