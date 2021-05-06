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
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)
   Trinayan Baruah, Northeastern University(baruah.t@northeastern.edu)
   Nicholas Curtis, AMD(nicholas.curtis@amd.com)
   David Kaeli,     Northeastern University(kaeli@ece.neu.edu)
   Per-atom energy/virial added by Ray Shan (Sandia)
   Hybrid and hybrid/overlay compatibility added by Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "pair_reaxc_hip.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "fix.h"
#include "fix_reaxc_hip.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#include "reaxc_defs_hip.h"
#include "reaxc_types_hip.h"
#include "reaxc_allocate_hip.h"
#include "reaxc_control_hip.h"
#include "reaxc_ffield_hip.h"
#include "reaxc_init_md_hip.h"
#include "reaxc_io_tools_hip.h"
#include "reaxc_list_hip.h"
#include "reaxc_reset_tools_hip.h"
#include "reaxc_forces_hip.h"
#include "reaxc_vector_hip.h"

extern "C" void Cuda_Init_Block_Sizes( reax_system *system, control_params *control );
extern "C" void Setup_Cuda_Environment( int, int, int );
extern "C" void Cuda_Initialize( reax_system*, control_params*, simulation_data*,
		storage*,reax_list**, reax_list*,output_controls*, mpi_datatypes* );


extern "C" void Cuda_Adjust_End_Index_Before_ReAllocation(int oldN, int systemN, reax_list **gpu_lists);


extern "C" void Cuda_ReAllocate( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace, reax_list **lists);


extern "C" void Cuda_Reset( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace, reax_list **lists);


extern "C" void Cuda_Write_Reax_Lists(reax_system *system, reax_list**, reax_list*);


extern "C" int Cuda_Compute_Forces( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace, reax_list **lists,
		output_controls *out_control, mpi_datatypes *mpi_data );

extern "C" void  Cuda_Allocate_Atoms(reax_system *system);

extern "C" void Cuda_Update_Atoms_On_Device(reax_system *system);

extern "C" void Cuda_Make_List( int n, int num_intrs, int type, reax_list *l);

extern "C" void Output_Sync_Forces(storage *workspace, int total_cap);

extern "C" void  Output_Sync_Atoms( reax_system *sys );

extern "C" void Output_Sync_Simulation_Data( simulation_data *host, simulation_data *dev );



using namespace LAMMPS_NS;

static const char cite_pair_reax_c[] =
		"pair reax/c command:\n\n"
		"@Article{Aktulga12,\n"
		" author = {H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama},\n"
		" title = {Parallel reactive molecular dynamics: Numerical methods and algorithmic techniques},\n"
		" journal = {Parallel Computing},\n"
		" year =    2012,\n"
		" volume =  38,\n"
		" pages =   {245--259}\n"
		"}\n\n";

/* ---------------------------------------------------------------------- */

PairReaxCGPU::PairReaxCGPU(LAMMPS *lmp) : Pair(lmp)
{
	if (lmp->citeme) lmp->citeme->add(cite_pair_reax_c);

	single_enable = 0;
	restartinfo = 0;
	one_coeff = 1;
	manybody_flag = 1;
	ghostneigh = 1;

	fix_id = new char[24];
	snprintf(fix_id,24,"REAXC_%d",instance_me);

	system = (reax_system *)
    																																										memory->smalloc(sizeof(reax_system),"reax:system");
	memset(system,0,sizeof(reax_system));
	control = (control_params *)
    																																										memory->smalloc(sizeof(control_params),"reax:control");
	memset(control,0,sizeof(control_params));
	data = (simulation_data *)
    																																										memory->smalloc(sizeof(simulation_data),"reax:data");
	workspace = (storage *)
    																																										memory->smalloc(sizeof(storage),"reax:storage");

	workspace->d_workspace = (storage *)memory->smalloc(sizeof(storage),"reax:gpu_storage");


	grid * const g = &system->my_grid;

	gpu_lists = (reax_list **)memory->smalloc(sizeof(reax_list*) * LIST_N ,"reax:gpu_lists");
	for ( int i = 0; i < LIST_N; ++i )
	{
		gpu_lists[i] = (reax_list *)memory->smalloc( sizeof(reax_list),
				"Setup::pmd_handle->lists[i]" );
		gpu_lists[i]->allocated = FALSE;
	}

	cpu_lists = (reax_list *)memory->smalloc(LIST_N * sizeof(reax_list),"reax:lists");
	memset(cpu_lists,0,LIST_N * sizeof(reax_list));



	out_control = (output_controls *)memory->smalloc(sizeof(output_controls),"reax:out_control");
	memset(out_control,0,sizeof(output_controls));
	mpi_data = (mpi_datatypes *)memory->smalloc(sizeof(mpi_datatypes),"reax:mpi");
	control->me = system->my_rank = comm->me;


	system->my_coords[0] = 0;
	system->my_coords[1] = 0;
	system->my_coords[2] = 0;
	system->num_nbrs = 0;
	system->n = 0; // my atoms
	system->N = 0; // mine + ghosts
	system->bigN = 0;  // all atoms in the system
	system->local_cap = 0;
	system->total_cap = 0;
	system->gcell_cap = 0;
	system->bndry_cuts.ghost_nonb = 0;
	system->bndry_cuts.ghost_hbond = 0;
	system->bndry_cuts.ghost_bond = 0;
	system->bndry_cuts.ghost_cutoff = 0;
	system->my_atoms = NULL;
	system->pair_ptr = this;
	system->error_ptr = error;
	control->error_ptr = error;

	system->omp_active = 0;

	fix_reax = NULL;
	tmpid = NULL;
	tmpbo = NULL;

	nextra = 14;
	pvector = new double[nextra];

	setup_flag = 0;
	fixspecies_flag = 0;

	nmax = 0;

}

/* ---------------------------------------------------------------------- */

PairReaxCGPU::~PairReaxCGPU()
{

	if (copymode) return;

	if (fix_reax) modify->delete_fix(fix_id);
	delete[] fix_id;

	if (setup_flag) {
		Close_Output_Files( system, control, out_control, mpi_data );

		// deallocate reax data-structures

		if (control->tabulate )
		{
		  printf("Tabulate option is not supported\n");
		}

		if (control->hbond_cut > 0 )  Delete_List( cpu_lists+HBONDS );
		Delete_List( cpu_lists+BONDS );
		Delete_List( cpu_lists+THREE_BODIES );
		Delete_List( cpu_lists+FAR_NBRS );

		DeAllocate_Workspace( control, workspace );
		DeAllocate_System( system );
	}

	memory->destroy( system );
	memory->destroy( control );
	memory->destroy( data );
	memory->destroy( workspace );
	memory->destroy( cpu_lists );
	memory->destroy( out_control );
	memory->destroy( mpi_data );

	// deallocate interface storage
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(cutghost);
		delete [] map;

		delete [] chi;
		delete [] eta;
		delete [] gamma;
	}

	memory->destroy(tmpid);
	memory->destroy(tmpbo);

	delete [] pvector;

}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::allocate( )
{
	allocated = 1;
	int n = atom->ntypes;

	memory->create(setflag,n+1,n+1,"pair:setflag");
	memory->create(cutsq,n+1,n+1,"pair:cutsq");
	memory->create(cutghost,n+1,n+1,"pair:cutghost");
	map = new int[n+1];

	chi = new double[n+1];
	eta = new double[n+1];
	gamma = new double[n+1];
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::settings(int narg, char **arg)
{
	if (narg < 1) error->all(FLERR,"Illegal pair_style command");

	// read name of control file or use default controls

	if (strcmp(arg[0],"NULL") == 0) {
		strcpy( control->sim_name, "simulate" );
		control->ensemble = 0;
		out_control->energy_update_freq = 0;
		control->tabulate = 0;

		control->reneighbor = 1;
		control->vlist_cut = control->nonb_cut;
		control->bond_cut = 5.;
		control->hbond_cut = 7.50;
		control->thb_cut = 0.001;
		control->thb_cutsq = 0.00001;
		control->bg_cut = 0.3;

		// Initialize for when omp style included
		control->nthreads = 1;

		out_control->write_steps = 0;
		out_control->traj_method = 0;
		strcpy( out_control->traj_title, "default_title" );
		out_control->atom_info = 0;
		out_control->bond_info = 0;
		out_control->angle_info = 0;
	} else Read_Control_File(arg[0], control, out_control);

	// default values

	qeqflag = 1;
	control->lgflag = 0;
	control->enobondsflag = 1;
	system->mincap = MIN_CAP;
	system->safezone = SAFE_ZONE;
	system->saferzone = SAFER_ZONE;

	// process optional keywords

	int iarg = 1;

	while (iarg < narg) {
		if (strcmp(arg[iarg],"checkqeq") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
			if (strcmp(arg[iarg+1],"yes") == 0) qeqflag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) qeqflag = 0;
			else error->all(FLERR,"Illegal pair_style reax/c command");
			iarg += 2;
		} else if (strcmp(arg[iarg],"enobonds") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
			if (strcmp(arg[iarg+1],"yes") == 0) control->enobondsflag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) control->enobondsflag = 0;
			else error->all(FLERR,"Illegal pair_style reax/c command");
			iarg += 2;
		} else if (strcmp(arg[iarg],"lgvdw") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
			if (strcmp(arg[iarg+1],"yes") == 0) control->lgflag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) control->lgflag = 0;
			else error->all(FLERR,"Illegal pair_style reax/c command");
			iarg += 2;
		} else if (strcmp(arg[iarg],"safezone") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
			system->safezone = utils::numeric(FLERR,arg[iarg+1],false,lmp);
			if (system->safezone < 0.0)
				error->all(FLERR,"Illegal pair_style reax/c safezone command");
			system->saferzone = system->safezone*1.2 + 0.2;
			iarg += 2;
		} else if (strcmp(arg[iarg],"mincap") == 0) {
			if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
			system->mincap = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
			if (system->mincap < 0)
				error->all(FLERR,"Illegal pair_style reax/c mincap command");
			iarg += 2;
		} else error->all(FLERR,"Illegal pair_style reax/c command");
	}

	// LAMMPS is responsible for generating nbrs

	control->reneighbor = 1;

	Setup_Cuda_Environment(system->my_rank,
			control->nprocs, control->gpus_per_node );


}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::coeff( int nargs, char **args )
{
	if (!allocated) allocate();

	if (nargs != 3 + atom->ntypes)
		error->all(FLERR,"Incorrect args for pair coefficients");

	// insure I,J args are * *

	if (strcmp(args[0],"*") != 0 || strcmp(args[1],"*") != 0)
		error->all(FLERR,"Incorrect args for pair coefficients");

	// read ffield file

	char *file = args[2];
	FILE *fp;
	fp = utils::open_potential(file,lmp,nullptr);
	if (fp != NULL)
		Read_Force_Field(fp, &(system->reax_param), control);
	else {
		char str[128];
		snprintf(str,128,"Cannot open ReaxFF potential file %s",file);
		error->all(FLERR,str);
	}

	// read args that map atom types to elements in potential file
	// map[i] = which element the Ith atom type is, -1 if NULL

	int itmp = 0;
	int nreax_types = system->reax_param.num_atom_types;
	for (int i = 3; i < nargs; i++) {
		if (strcmp(args[i],"NULL") == 0) {
			map[i-2] = -1;
			itmp ++;
			continue;
		}
	}

	int n = atom->ntypes;

	// pair_coeff element map
	for (int i = 3; i < nargs; i++)
		for (int j = 0; j < nreax_types; j++)
			if (strcasecmp(args[i],system->reax_param.sbp[j].name) == 0) {
				map[i-2] = j;
				itmp ++;
			}

	// error check
	if (itmp != n)
		error->all(FLERR,"Non-existent ReaxFF type");

	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;

	// set setflag i,j for type pairs where both are mapped to elements

	int count = 0;
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			if (map[i] >= 0 && map[j] >= 0) {
				setflag[i][j] = 1;
				count++;
			}

	if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");


}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::init_style( )
{
	if (!atom->q_flag)
		error->all(FLERR,"Pair style reax/c requires atom attribute q");

	// firstwarn = 1;

	bool have_qeq = ((modify->find_fix_by_style("^qeq/reax") != -1)
			|| (modify->find_fix_by_style("^qeq/shielded") != -1));
	if (!have_qeq && qeqflag == 1)
		error->all(FLERR,"Pair reax/c requires use of fix qeq/reax or qeq/shielded");

	system->n = atom->nlocal; // my atoms
	system->N = atom->nlocal + atom->nghost; // mine + ghosts
	system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system
	system->wsize = comm->nprocs;

	system->big_box.V = 0;
	system->big_box.box_norms[0] = 0;
	system->big_box.box_norms[1] = 0;
	system->big_box.box_norms[2] = 0;

	if (atom->tag_enable == 0)
		error->all(FLERR,"Pair style reax/c requires atom IDs");
	if (force->newton_pair == 0)
		error->all(FLERR,"Pair style reax/c requires newton pair on");
	if ((atom->map_tag_max > 99999999) && (comm->me == 0))
		error->warning(FLERR,"Some Atom-IDs are too large. Pair style reax/c "
				"native output files may get misformatted or corrupted");

	// because system->bigN is an int, we cannot have more atoms than MAXSMALLINT

	if (atom->natoms > MAXSMALLINT)
		error->all(FLERR,"Too many atoms for pair style reax/c");

	// need a half neighbor list w/ Newton off and ghost neighbors
	// built whenever re-neighboring occurs

	int irequest = neighbor->request(this,instance_me);
	//neighbor->requests[irequest]->newton = 2;
	neighbor->requests[irequest]->ghost = 1;
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;

	cutmax = MAX3(control->nonb_cut, control->hbond_cut, control->bond_cut);
	if ((cutmax < 2.0*control->bond_cut) && (comm->me == 0))
		error->warning(FLERR,"Total cutoff < 2*bond cutoff. May need to use an "
				"increased neighbor list skin.");

	/*  for( int i = 0; i < LIST_N; ++i )
    if (lists[i].allocated != 1)
      lists[i].allocated = 0;
	 */
	if (fix_reax == NULL) {
		char **fixarg = new char*[3];
		fixarg[0] = (char *) fix_id;
		fixarg[1] = (char *) "all";
		fixarg[2] = (char *) "REAXC";
		modify->add_fix(3,fixarg);
		delete [] fixarg;
		fix_reax = (FixReaxC *) modify->fix[modify->nfix-1];
	}
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::setup( )
{
	int oldN;
	int mincap = system->mincap;
	double safezone = system->safezone;

	system->n = atom->nlocal; // my atoms
	system->N = atom->nlocal + atom->nghost; // mine + ghosts
	oldN = system->N;
	system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system
	Cuda_Init_Block_Sizes(system, control);

	if (setup_flag == 0) {

		setup_flag = 1;

		int *num_bonds = fix_reax->num_bonds;
		int *num_hbonds = fix_reax->num_hbonds;

		control->vlist_cut = neighbor->cutneighmax;

		// determine the local and total capacity

		system->local_cap = MAX( (int)(system->n * safezone), mincap );
		system->total_cap = MAX( (int)(system->N * safezone), mincap );

		// initialize my data structures

		PreAllocate_Space( system, control, workspace);
		Cuda_Allocate_Atoms(system);
		update_and_copy_reax_atoms_to_device();

		int num_nbrs = estimate_reax_lists(); //TB:: Should this be moved to GPU?
		system->total_far_nbrs = num_nbrs;

		//printf("Num nbrs %d \n", num_nbrs);

		if(!Make_List(system->total_cap, num_nbrs, TYP_FAR_NEIGHBOR,
				(cpu_lists+FAR_NBRS)))
			error->one(FLERR,"Pair reax/c problem in far neighbor list");
		(cpu_lists+FAR_NBRS)->error_ptr=error;
		Cuda_Make_List(system->total_cap, system->total_far_nbrs,
				TYP_FAR_NEIGHBOR, gpu_lists[FAR_NBRS]);
		update_and_write_reax_lists_to_device();

		//exit(0);

		Initialize( system, control, data, workspace, &cpu_lists, out_control,
				mpi_data, world );

		Cuda_Initialize(system, control, data, workspace, gpu_lists,cpu_lists, out_control,
				mpi_data);


		for( int k = 0; k < system->N; ++k )
		{
			num_bonds[k] = system->my_atoms[k].num_bonds;
			num_hbonds[k] = system->my_atoms[k].num_hbonds;
		}
	}
	else
	{
		//printf("Realloc setup \n");
		// fill in reax datastructures
		update_and_copy_reax_atoms_to_device();

		// reset the bond list info for new atoms
		//printf("Initial setup done. far numbers gpu %d \n", gpu_lists[FAR_NBRS]->num_intrs);


		Cuda_Adjust_End_Index_Before_ReAllocation(oldN, system->N, gpu_lists);

		//printf("Initial setup done. far numbers gpu %d \n", gpu_lists[FAR_NBRS]->num_intrs);

		// check if I need to shrink/extend my data-structs

		Cuda_ReAllocate(system, control,
				data, workspace, gpu_lists);
	}

	bigint local_ngroup = list->inum;
	MPI_Allreduce( &local_ngroup, &ngroup, 1, MPI_LMP_BIGINT, MPI_SUM, world );
}

/* ---------------------------------------------------------------------- */

double PairReaxCGPU::init_one(int i, int j)
{
	if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

	cutghost[i][j] = cutghost[j][i] = cutmax;
	return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::compute(int eflag, int vflag)
{

	double evdwl,ecoul;
	double t_start, t_end;

	// communicate num_bonds once every reneighboring
	// 2 num arrays stored by fix, grab ptr to them

	if (neighbor->ago == 0) comm->forward_comm_fix(fix_reax);
	int *num_bonds = fix_reax->num_bonds;
	int *num_hbonds = fix_reax->num_hbonds;

	evdwl = ecoul = 0.0;
	ev_init(eflag,vflag);

	if (vflag_global) control->virial = 1;
	else control->virial = 0;

	system->n = atom->nlocal; // my atoms
	system->N = atom->nlocal + atom->nghost; // mine + ghosts
	system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system

	system->big_box.V = 0;
	system->big_box.box_norms[0] = 0;
	system->big_box.box_norms[1] = 0;
	system->big_box.box_norms[2] = 0;
	if (comm->me == 0 ) t_start = MPI_Wtime();

	// setup data structures

	setup();



	Cuda_Reset(system, control, data, workspace, gpu_lists);
	//printf("fixxxing \n");
	workspace->realloc.num_far = update_and_write_reax_lists_to_device();

	// timing for filling in the reax lists
	if (comm->me == 0) {
		t_end = MPI_Wtime();
		data->timing.nbrs = t_end - t_start;
	}


	// forces
	Cuda_Compute_Forces(system, control, data, workspace, gpu_lists, out_control, mpi_data);


	read_reax_forces_from_device(vflag);
	Output_Sync_Atoms(system);


	for(int k = 0; k < system->N; ++k)
	{
		num_bonds[k] = system->my_atoms[k].num_bonds;
		num_hbonds[k] = system->my_atoms[k].num_hbonds;
	}


	Output_Sync_Simulation_Data( data, (simulation_data *)data->d_simulation_data);



	// energies and pressure

	if (eflag_global) {

		evdwl += data->my_en.e_bond;
		evdwl += data->my_en.e_ov;
		evdwl += data->my_en.e_un;
		evdwl += data->my_en.e_lp;
		evdwl += data->my_en.e_ang;
		evdwl += data->my_en.e_pen;
		evdwl += data->my_en.e_coa;
		evdwl += data->my_en.e_hb;
		evdwl += data->my_en.e_tor;
		evdwl += data->my_en.e_con;
		evdwl += data->my_en.e_vdW;

		ecoul += data->my_en.e_ele;
		ecoul += data->my_en.e_pol;

		// eng_vdwl += evdwl;
		// eng_coul += ecoul;

		// Store the different parts of the energy
		// in a list for output by compute pair command

		pvector[0] = data->my_en.e_bond;
		pvector[1] = data->my_en.e_ov + data->my_en.e_un;
		pvector[2] = data->my_en.e_lp;
		pvector[3] = 0.0;
		pvector[4] = data->my_en.e_ang;
		pvector[5] = data->my_en.e_pen;
		pvector[6] = data->my_en.e_coa;
		pvector[7] = data->my_en.e_hb;
		pvector[8] = data->my_en.e_tor;
		pvector[9] = data->my_en.e_con;
		pvector[10] = data->my_en.e_vdW;
		pvector[11] = data->my_en.e_ele;
		pvector[12] = 0.0;
		pvector[13] = data->my_en.e_pol;
	}

	if (vflag_fdotr)
	{
		virial_fdotr_compute();
	}



	// Set internal timestep counter to that of LAMMPS

	data->step = update->ntimestep;



	Output_Results( system, control, data, gpu_lists, out_control, mpi_data );

	// populate tmpid and tmpbo arrays for fix reax/c/species
	int i, j;

	if(fixspecies_flag) {
		printf("fix species not implemented for GPU  %d\n",fixspecies_flag);

	}

	//printf("Finished loop \n");
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::update_and_copy_reax_atoms_to_device()
{
	int *num_bonds = fix_reax->num_bonds;
	int *num_hbonds = fix_reax->num_hbonds;

	if (system->N > system->total_cap)
		error->all(FLERR,"Too many ghost atoms");

	for( int i = 0; i < system->N; ++i ){
		system->my_atoms[i].orig_id = atom->tag[i];
		system->my_atoms[i].type = map[atom->type[i]];
		system->my_atoms[i].x[0] = atom->x[i][0];
		system->my_atoms[i].x[1] = atom->x[i][1];
		system->my_atoms[i].x[2] = atom->x[i][2];
		system->my_atoms[i].q = atom->q[i];
		system->my_atoms[i].num_bonds = num_bonds[i];
		system->my_atoms[i].num_hbonds = num_hbonds[i];
	}

	Cuda_Update_Atoms_On_Device(system);
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::get_distance( rvec xj, rvec xi, double *d_sqr, rvec *dvec )
{
	(*dvec)[0] = xj[0] - xi[0];
	(*dvec)[1] = xj[1] - xi[1];
	(*dvec)[2] = xj[2] - xi[2];
	*d_sqr = SQR((*dvec)[0]) + SQR((*dvec)[1]) + SQR((*dvec)[2]);
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::set_far_nbr( far_neighbor_data *fdest,
		int j, double d, rvec dvec )
{

	fdest->nbr = j;
	fdest->d = d;
	rvec_Copy( fdest->dvec, dvec );
	ivec_MakeZero( fdest->rel_box );
}

/* ---------------------------------------------------------------------- */

int PairReaxCGPU::estimate_reax_lists()
{
	int itr_i, itr_j, i, j;
	int num_nbrs, num_marked;
	int *ilist, *jlist, *numneigh, **firstneigh, *marked;
	double d_sqr;
	rvec dvec;
	double **x;

	int mincap = system->mincap;
	double safezone = system->safezone;

	x = atom->x;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	num_nbrs = 0;
	num_marked = 0;
	marked = (int*) calloc( system->N, sizeof(int) );

	int numall = list->inum + list->gnum;

	for( itr_i = 0; itr_i < numall; ++itr_i ){
		i = ilist[itr_i];
		marked[i] = 1;
		++num_marked;
		jlist = firstneigh[i];

		for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ){
			j = jlist[itr_j];
			if( i < j)
			{
				j &= NEIGHMASK;
				get_distance( x[j], x[i], &d_sqr, &dvec );

				if (d_sqr <= SQR(control->nonb_cut))
					++num_nbrs;
			}
		}

		for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ){
			j = jlist[itr_j];
			if( i > j)
			{
				j &= NEIGHMASK;
				get_distance( x[i], x[j], &d_sqr, &dvec );

				if (d_sqr <= SQR(control->nonb_cut))
				{
					++num_nbrs;
				}
			}
		}

	}


	free(marked);


	return static_cast<int> (MAX( num_nbrs*safezone, mincap*MIN_NBRS ));
}

/* ---------------------------------------------------------------------- */

int PairReaxCGPU::update_and_write_reax_lists_to_device()
{


	int itr_i, itr_j, i, j;
	int num_nbrs;
	int *ilist, *jlist, *numneigh, **firstneigh;
	double d_sqr, cutoff_sqr;
	rvec dvec;
	double *dist, **x;
	reax_list *far_nbrs;
	far_neighbor_data *far_list;

	x = atom->x;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	far_nbrs = (cpu_lists +FAR_NBRS);
	far_list = far_nbrs->select.far_nbr_list;



	num_nbrs = 0;
	int inum = list->inum;
	dist = (double*) calloc( system->N, sizeof(double) );



	int numall = list->inum + list->gnum;

	for( itr_i = 0; itr_i < numall; ++itr_i ){
		i = ilist[itr_i];
		jlist = firstneigh[i];

		Set_Start_Index( i, num_nbrs, far_nbrs );

		if (i < inum)
			cutoff_sqr = control->nonb_cut*control->nonb_cut;
		else
			cutoff_sqr = control->bond_cut*control->bond_cut;

		for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ) {
			j = jlist[itr_j];
			if ( i <  j)
			{
				j &= NEIGHMASK;

				get_distance( x[j], x[i], &d_sqr, &dvec );

				if (d_sqr <= (cutoff_sqr)) {
					dist[j] = sqrt( d_sqr );
					set_far_nbr( &far_list[num_nbrs], j, dist[j], dvec );
					++num_nbrs;
				}
			}
		}

		for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ) {
			j = jlist[itr_j];
			if ( i >  j)
			{
				j &= NEIGHMASK;

				get_distance( x[i], x[j], &d_sqr, &dvec );

				if (d_sqr <= (cutoff_sqr)) {
					dist[j] = sqrt( d_sqr );
					set_far_nbr( &far_list[num_nbrs], j, dist[j], dvec );
					++num_nbrs;
				}
			}
		}





		Set_End_Index( i, num_nbrs, far_nbrs );
	}

	free( dist );

	Cuda_Write_Reax_Lists(system,  gpu_lists, cpu_lists);


	return num_nbrs;
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::read_reax_forces_from_device(int /*vflag*/)
{

	Output_Sync_Forces(workspace,system->total_cap);


	int world_rank;
	MPI_Comm_rank(world, &world_rank);


	for( int i = 0; i < system->N; ++i ) {
		system->my_atoms[i].f[0] = workspace->f[i][0];
		system->my_atoms[i].f[1] = workspace->f[i][1];
		system->my_atoms[i].f[2] = workspace->f[i][2];



		//if(i < 20)
			//printf("%d,%f,%f,%f\n",system->my_atoms[i].orig_id, system->my_atoms[i].f[0],system->my_atoms[i].f[1],system->my_atoms[i].f[2]);


		atom->f[i][0] += -workspace->f[i][0];
		atom->f[i][1] += -workspace->f[i][1];
		atom->f[i][2] += -workspace->f[i][2];
	}

	//exit(0);

	//printf("Computation done\n");


}

/* ---------------------------------------------------------------------- */

void *PairReaxCGPU::extract(const char *str, int &dim)
{
	dim = 1;
	if (strcmp(str,"chi") == 0 && chi) {
		for (int i = 1; i <= atom->ntypes; i++)
			if (map[i] >= 0) chi[i] = system->reax_param.sbp[map[i]].chi;
			else chi[i] = 0.0;
		return (void *) chi;
	}
	if (strcmp(str,"eta") == 0 && eta) {
		for (int i = 1; i <= atom->ntypes; i++)
			if (map[i] >= 0) eta[i] = system->reax_param.sbp[map[i]].eta;
			else eta[i] = 0.0;
		return (void *) eta;
	}
	if (strcmp(str,"gamma") == 0 && gamma) {
		for (int i = 1; i <= atom->ntypes; i++)
			if (map[i] >= 0) gamma[i] = system->reax_param.sbp[map[i]].gamma;
			else gamma[i] = 0.0;
		return (void *) gamma;
	}
	return NULL;
}

/* ---------------------------------------------------------------------- */

double PairReaxCGPU::memory_usage()
{

	double bytes = 0.0;

	// From pair_reax_c
	bytes += 1.0 * system->N * sizeof(int);
	bytes += 1.0 * system->N * sizeof(double);

	// From reaxc_allocate: BO
	bytes += 1.0 * system->total_cap * sizeof(reax_atom);
	bytes += 19.0 * system->total_cap * sizeof(double);
	bytes += 3.0 * system->total_cap * sizeof(int);

	// From reaxc_lists
	/*bytes += 2.0 * lists->n * sizeof(int);
  bytes += lists->num_intrs * sizeof(three_body_interaction_data);
  bytes += lists->num_intrs * sizeof(bond_data);
  bytes += lists->num_intrs * sizeof(dbond_data);
  bytes += lists->num_intrs * sizeof(dDelta_data);
  bytes += lists->num_intrs * sizeof(far_neighbor_data);
  bytes += lists->num_intrs * sizeof(hbond_data);

  if(fixspecies_flag)
    bytes += 2 * nmax * MAXSPECBOND * sizeof(double);*/

	return bytes;
}

/* ---------------------------------------------------------------------- */

void PairReaxCGPU::FindBond()
{
}
