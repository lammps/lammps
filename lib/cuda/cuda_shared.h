/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef _CUDA_SHARED_H_
#define _CUDA_SHARED_H_
#include "cuda_precision.h"

#define CUDA_MAX_DEBUG_SIZE 1000 //size of debugdata array (allows for so many doubles or twice as many int)

struct dev_array {
  void* dev_data;			// pointer to memory address on cuda device
  unsigned dim[3];		// array dimensions
};

struct cuda_shared_atom {	// relevent data from atom class
  dev_array dx; 			// cumulated distance for binning settings
  dev_array x;			// position
  dev_array v;			// velocity
  dev_array f;			// force
  dev_array tag;
  dev_array type; 		// global ID number, there are ghosttype = ntypes  (ntypescuda=ntypes+1)
  dev_array mask;
  dev_array image;
  dev_array q;			// charges
  dev_array mass;			// per-type masses
  dev_array rmass;		// per-atom masses
  dev_array radius;		// per-atom radius
  dev_array density;
  dev_array omega;
  dev_array torque;
  dev_array molecule;

  dev_array special;
  int maxspecial;
  dev_array nspecial;
  int* special_flag;
  int molecular;

  dev_array eatom;		// per-atom energy
  dev_array vatom;		// per-atom virial
  int need_eatom;
  int need_vatom;

  dev_array x_type;		// position + type in X_CFLOAT4 struct
  dev_array v_radius;		// velociyt + radius in V_CFLOAT4 struct currently only used for granular atom_style
  dev_array omega_rmass;		// velociyt + radius in V_CFLOAT4 struct currently only used for granular atom_style

  double* mass_host;		// remember per-type host pointer to masses
  //int natoms;				// total # of atoms in system, could be 0
  int nghost;				// and ghost atoms on this proc
  int nlocal;				// # of owned
  int nall;			    // total # of atoms in this proc
  int nmax;				// max # of owned+ghost in arrays on this proc
  int ntypes;
  int q_flag;				// do we have charges?
  int rmass_flag;			// do we have per-atom masses?
  int firstgroup;
  int nfirst;

  int update_nlocal;
  int update_nmax;
  int update_neigh;

  dev_array xhold;	    // position at last neighboring
  X_CFLOAT triggerneighsq;		// maximum square movement before reneighboring
  int reneigh_flag;		// is reneighboring necessary
  int maxhold;			// size of xhold
  int dist_check; 		//perform distance check for reneighboring
  dev_array binned_id;    //id of each binned atom (not tag!!)
  dev_array binned_idnew; //new id of each binned atom for sorting basically setting atom[binned_id[k]] at atom[binned_newid[k]]
  float bin_extraspace;
  int bin_dim[3];
  int bin_nmax;
  dev_array map_array;
};

struct cuda_shared_pair {	// relevent data from pair class
  char cudable_force;		// check for (cudable_force!=0)
  X_CFLOAT cut_global;
  X_CFLOAT cut_inner_global;
  X_CFLOAT cut_coul_global;
  double** cut;			// type-type cutoff
  double** cutsq;			// type-type cutoff
  double** cut_inner;			// type-type cutoff for coul
  double** cut_coul;			// type-type cutoff for coul
  double** coeff1;		// tpye-type pair parameters
  double** coeff2;
  double** coeff3;
  double** coeff4;
  double** coeff5;
  double** coeff6;
  double** coeff7;
  double** coeff8;
  double** coeff9;
  double** coeff10;
  double** offset;
  double* special_lj;
  double* special_coul;
  dev_array virial; // ENERGY_CFLOAT
  dev_array eng_vdwl; // ENERGY_CFLOAT
  dev_array eng_coul; // ENERGY_CFLOAT
  X_CFLOAT cut_coulsq_global;
  F_CFLOAT g_ewald, kappa;
  int freeze_group_bit;

  dev_array coeff1_gm;
  dev_array coeff2_gm;
  dev_array coeff3_gm;
  dev_array coeff4_gm;
  dev_array coeff5_gm;
  dev_array coeff6_gm;
  dev_array coeff7_gm;
  dev_array coeff8_gm;
  dev_array coeff9_gm;
  dev_array coeff10_gm;

  int lastgridsize;
  int n_energy_virial;
  int collect_forces_later;
  int use_block_per_atom;
  int override_block_per_atom;
  bool neighall;

};

struct cuda_shared_domain {	// relevent data from domain class
  X_CFLOAT sublo[3];			// orthogonal box -> sub-box bounds on this proc
  X_CFLOAT subhi[3];
  X_CFLOAT boxlo[3];
  X_CFLOAT boxhi[3];
  X_CFLOAT prd[3];
  int periodicity[3];		// xyz periodicity as array

  int triclinic;
  X_CFLOAT xy;
  X_CFLOAT xz;
  X_CFLOAT yz;
  X_CFLOAT boxlo_lamda[3];
  X_CFLOAT boxhi_lamda[3];
  X_CFLOAT prd_lamda[3];
  X_CFLOAT h[6];
  X_CFLOAT h_inv[6];
  V_CFLOAT h_rate[6];
  int update;
};

struct cuda_shared_pppm {
  char cudable_force;
#ifdef FFT_CUFFT
  FFT_CFLOAT* work1;
  FFT_CFLOAT* work2;
  FFT_CFLOAT* work3;
  PPPM_CFLOAT* greensfn;
  PPPM_CFLOAT* fkx;
  PPPM_CFLOAT* fky;
  PPPM_CFLOAT* fkz;
  PPPM_CFLOAT* vg;
#endif
  int* part2grid;
  PPPM_CFLOAT* density_brick;
  int* density_brick_int;
  PPPM_CFLOAT density_intScale;
  PPPM_CFLOAT* vdx_brick;
  PPPM_CFLOAT* vdy_brick;
  PPPM_CFLOAT* vdz_brick;
  PPPM_CFLOAT* density_fft;
  ENERGY_CFLOAT* energy;
  ENERGY_CFLOAT* virial;
  int nxlo_in;
  int nxhi_in;
  int nxlo_out;
  int nxhi_out;
  int nylo_in;
  int nyhi_in;
  int nylo_out;
  int nyhi_out;
  int nzlo_in;
  int nzhi_in;
  int nzlo_out;
  int nzhi_out;
  int nx_pppm;
  int ny_pppm;
  int nz_pppm;
  PPPM_CFLOAT qqrd2e;
  int order;
  // float3 sublo;
  PPPM_CFLOAT* rho_coeff;
  int nmax;
  int nlocal;
  PPPM_CFLOAT* debugdata;
  PPPM_CFLOAT delxinv;
  PPPM_CFLOAT delyinv;
  PPPM_CFLOAT delzinv;
  int nlower;
  int nupper;
  PPPM_CFLOAT shiftone;
  PPPM_CFLOAT3* fH;
};

struct cuda_shared_comm {
  int maxswap;
  int maxlistlength;
  dev_array pbc;
  dev_array slablo;
  dev_array slabhi;
  dev_array multilo;
  dev_array multihi;
  dev_array sendlist;
  int grow_flag;
  int comm_phase;

  int nsend;
  int* nsend_swap;
  int* send_size;
  int* recv_size;
  double** buf_send;
  void** buf_send_dev;
  double** buf_recv;
  void** buf_recv_dev;
  void* buffer;
  int buffer_size;
  double overlap_split_ratio;
};

struct cuda_shared_neighlist { // member of CudaNeighList, has no instance in cuda_shared_data
  int maxlocal;
  int inum;                // # of I atoms neighbors are stored for local indices of I atoms
  int inum_border2;
  dev_array inum_border;         // # of atoms which interact with border atoms
  dev_array ilist;
  dev_array ilist_border;
  dev_array numneigh;
  dev_array numneigh_inner;
  dev_array numneigh_border;
  dev_array firstneigh;
  dev_array neighbors;
  dev_array neighbors_border;
  dev_array neighbors_inner;
  int maxpage;
  dev_array page_pointers;
  dev_array* pages;
  int maxneighbors;
  int neigh_lists_per_page;
  double** cutneighsq;
  CUDA_CFLOAT* cu_cutneighsq;
  int* binned_id;
  int* bin_dim;
  int bin_nmax;
  float bin_extraspace;
  double maxcut;
  dev_array ex_type;
  int nex_type;
  dev_array ex1_bit;
  dev_array ex2_bit;
  int nex_group;
  dev_array ex_mol_bit;
  int nex_mol;

};

struct cuda_compile_settings {	// this is used to compare compile settings (i.e. precision) of the cu files, and the cpp files
  int prec_glob;
  int prec_x;
  int prec_v;
  int prec_f;
  int prec_pppm;
  int prec_fft;
  int cufft;
  int arch;
};

struct cuda_timings_struct {
  //Debug:
  double test1;
  double test2;
  //transfers
  double transfer_upload_tmp_constr;
  double transfer_download_tmp_deconstr;

  //communication
  double comm_forward_total;
  double comm_forward_mpi_upper;
  double comm_forward_mpi_lower;
  double comm_forward_kernel_pack;
  double comm_forward_kernel_unpack;
  double comm_forward_kernel_self;
  double comm_forward_upload;
  double comm_forward_download;

  double comm_exchange_total;
  double comm_exchange_mpi;
  double comm_exchange_kernel_pack;
  double comm_exchange_kernel_unpack;
  double comm_exchange_kernel_fill;
  double comm_exchange_cpu_pack;
  double comm_exchange_upload;
  double comm_exchange_download;

  double comm_border_total;
  double comm_border_mpi;
  double comm_border_kernel_pack;
  double comm_border_kernel_unpack;
  double comm_border_kernel_self;
  double comm_border_kernel_buildlist;
  double comm_border_upload;
  double comm_border_download;

  //pair forces
  double pair_xtype_conversion;
  double pair_kernel;
  double pair_virial;
  double pair_force_collection;

  //neighbor
  double neigh_bin;
  double neigh_build;
  double neigh_special;

  //PPPM
  double pppm_particle_map;
  double pppm_make_rho;
  double pppm_brick2fft;
  double pppm_poisson;
  double pppm_fillbrick;
  double pppm_fieldforce;
  double pppm_compute;

};

struct cuda_shared_data {	// holds space for all relevent data from the different classes
  void* buffer; //holds temporary GPU data [data used in subroutines, which has not to be consistend outside of that routine]
  int buffersize; //maxsize of buffer
  int buffer_new; //should be 1 if the pointer to buffer has changed
  void* flag;
  void* debugdata;  //array for easily collecting debugdata from device class cuda contains the corresponding cu_debugdata and host array
  cuda_shared_atom atom;
  cuda_shared_pair pair;
  cuda_shared_domain domain;
  cuda_shared_pppm pppm;
  cuda_shared_comm comm;
  cuda_compile_settings compile_settings;
  cuda_timings_struct cuda_timings;
  int exchange_dim;
  int me; //mpi rank
  unsigned int datamask;
  int overlap_comm;
};


#endif // #ifndef _CUDA_SHARED_H_
