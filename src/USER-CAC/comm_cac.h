/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMM_CLASS

CommStyle(cac,CommCAC)

#else

#ifndef LMP_COMM_CAC_H
#define LMP_COMM_CAC_H

#include "comm_tiled.h"

namespace LAMMPS_NS {

class CommCAC : public CommTiled {
 public:
  
  CommCAC(class LAMMPS *);
  CommCAC(class LAMMPS *, class Comm *);
  virtual ~CommCAC();

  virtual void post_constructor();
  virtual void init();
  virtual void setup();                        // setup comm pattern
  virtual void forward_comm(int dummy = 0);    // forward comm of atom coords
  
  virtual void exchange();                     // move atoms to new procs
  virtual void borders();                      // setup list of atoms to comm

  virtual void coord2proc_setup();
  virtual int coord2proc(double *, int &, int &, int &);

  bigint memory_usage();

 private:
  int nswap;                    // # of swaps to perform = 2*dim
  int nswap_border;                    // # of swaps to perform = 2*dim

  // forward/reverse comm info, proc lists include self

  int *nsendproc,*nrecvproc;    // # of procs to send/recv to/from per swap
  int *sendother,*recvother;    // 1 if send/recv to/from other proc per swap
  int *sendself;                // 1 if send to self per swap
  int *nprocmax;                // current max # of send procs per swap
  int *nrecv_procmax;                // current max # of send procs per swap
  int **sendproc,**recvproc;    // procs to send/recv to/from per swap
  int **sendbox_flag;           //decide whether this send should add data to sendlist
  int **repeatsend_flag;           //decide whether this send should add data to sendlist
  int **sendsize,**recvsize; // size of buffer sent by each overlap proc per swap
  int **sendoffset,**recvoffset; // offset of buffer for each overlap proc per swap
  int **overlap_sendsize,**overlap_recvsize; // size of buffer sent by each overlap proc per swap for overlapping elements
  int **overlap_sendoffset,**overlap_recvoffset; // offset of buffer for each overlap proc per swap for overlapping elements
  int **sendnum,**recvnum,**overlap_sendnum,**overlap_recvnum;      // # of atoms to send/recv per swap/proc
  int **size_forward_recv;      // # of values to recv in each forward swap/proc
  int **firstrecv, **overlap_firstrecv;              // where to put 1st recv atom per swap/proc
  int **size_reverse_send;      // # of values to send in each reverse swap/proc
  int **size_reverse_recv;      // # of values to recv in each reverse swap/proc
  int **forward_recv_offset;  // forward comm offsets in buf_recv per swap/proc
  int **reverse_recv_offset;  // reverse comm offsets in buf_recv per swap/proc
  double **cutghostmulti;           // cutghost on a per-type basis
  double **cutghostCAC;           // cutghost on a per-element scale basis for CAC package
  int ***sendlist, ***overlap_sendlist, ***tag_sendlist, ***repeat_list;  // list of atoms to send per swap/proc
  int **maxsendlist, **overlap_maxsendlist;            // max size of send list per swap/proc
  double ***aug_oboxes;         // other task boxes augmented to include all their local element bounding boxes
  int **pbc_flag;               // general flag for sending atoms thru PBC
  int ***pbc;                   // dimension flags for PBC adjustments
  double **eboxes;              // element bounding boxes
  double **foreign_eboxes;      //eboxes of other tasks communicated due to element overlap in my box
  int *ebox_ref;                 //local element index for this ebox ranging from 1:nlocal+nghost
  int neboxes;                  //number of element bounding boxes in me
  int local_neboxes;            //number of element bounding boxes in me that don't bound ghosts
  int maxebox;                  //maximum size of ebox array
  int maxforeign_ebox;
  int nforeign_eboxes;          //number of ghost eboxes in me
  int *foreign_eprocs;          //set of ghost eboxes in me declared by local element index 1:nghost+number of eboxes sent in overlap-comm step
  int **foreign_image;          //set of ghost eboxes in me declared by local element index 1:nghost+number of eboxes sent in overlap-comm step
  int *foreign_swaps;          //set of ghost eboxes in me declared by local element index 1:nghost+number of eboxes sent in overlap-comm step
  int nebox_considered;          //number of local and ghost eboxes that have already performed their role in previous swaps
  int ebox_limit;               //total number of eboxes nlocal+nghost+nforeign
  int ebox_limit_recheck;       //total number of eboxes nlocal+nforeign
  int elimit;                   //element limit 
  int recheck_flag;             //determines if ghost vs. nforeign box overlaps have been rechecked and communicated for corner cases
  int ***sent_flag;             //determines if the current element was already sent for the given swap and overlap proc in the first border comm round
  int *recv_flag;              //determines in which swap the current element was received
  int **maxsent;               
  int maxall;
  int *overlap_repeat;           //stores flags for each proc to determine if overlap array has repeats in O(P)
  int *work1,*work2;                // work vectors
  int foreign_swap;             //stores swap index in which a foreign ebox was sent
  double element_overlap_range[6]; //upper bound on range than an element can overlap into another task's subbox
  double aug_box[6];             //subbox of me expanded by element overlap of local elements
  int pbc_overlap;               //flag that guides the behavior of box_drop
  int current_pbc[3];           //grid offset for this proc need search (corner, face, edge)
  double **lo2_set, **hi2_set;             //stores overlap boxes to compute send/recv need for pbc images
  double **current_pbc_set;             //stores overlap boxes to compute send/recv need for pbc images
  int max_image;                //max image count a single overlap check has incurred
  int nstencil;                 //local variables to bin and stencil info
  int *stencil;
  int *bin_ncontent;
  int **bin_content;
  int *nbin_element_overlap;  //array storing the number of bins this element overlaps
  int **bin_element_overlap;  //set of bins this element overlaps
 
  double ***sendbox;            // bounding box of atoms to send per swap/proc
  double ***overlap_sendbox;            // bounding box of atoms to send per swap/proc
  double ****sendbox_multi;     // bounding box of atoms to send per swap/proc for multi comm
  double max_search_range;       //largest bin scale range
  // exchange comm info, proc lists do not include self

  int *nexchproc;               // # of procs to send/recv to/from in each dim
  int *nexchprocmax;            // current max # of exch procs for each dim
  int **exchproc;               // procs to exchange with per dim
  int **exchnum;                // # of values received per dim/proc

 
  double *buf_send;             // send buffer for all comm
  double *buf_recv;             // recv buffer for all comm
  int maxsend,maxrecv;          // current size of send/recv buffer
  
  int bufextra;                 // extra space beyond maxsend in send buffer
  int smaxone,rmaxone;          // max size in atoms of single borders send/recv
  int smaxall,rmaxall;          // max size in atoms of any borders send/recv
                                //   for comm to all procs in one swap

  int maxreqstat;               // max size of Request and Status vectors
  MPI_Request *requests;

  struct RCBinfo {
    double mysplit[3][2];      // fractional RCB bounding box for one proc
    double cutfrac;            // fractional position of cut this proc owns
    int dim;                   // dimension = 0/1/2 of cut
  };

  RCBinfo *rcbinfo;            // list of RCB info for all procs

  int noverlap;                // # of overlapping procs
  int maxoverlap;              // current max length of overlap
  int maxoverlap_box;              // current max length of overlap
  int *overlap;                // list of overlapping procs
  int **overlap_pbc;           // list of image offsets in each dim for each overlap proc
  double **proc2box;           // list of proc boxes setup by overlap calculation for brick    
  int overlap_counter;          // global counter used in box other to retain the typedef for box other

  double *prd;                 // local ptrs to Domain attributes
  double *boxlo,*boxhi;
  double *sublo,*subhi;
  int dimension;

  // NOTE: init_buffers is called from a constructor and must not be made virtual
  void init_buffers();

  // box drop and other functions

  typedef void (CommCAC::*BoxDropPtr)(int, double *, double *, int &);
  BoxDropPtr box_drop;
  typedef void (CommCAC::*BoxDropFullPtr)(int, double *, double *, int &);
  BoxDropFullPtr box_drop_full;
  void box_drop_brick(int, double *, double *, int &);
  void box_drop_brick_full(int, double *, double *, int &);
  void box_drop_tiled(int, double *, double *, int &);
  void box_drop_tiled_recurse(double *, double *, int, int, int &);
  void box_drop_tiled_full(int, double *, double *, int &);
  void box_drop_tiled_recurse_full(double *, double *, int, int, int &);

  typedef void (CommCAC::*BoxOtherPtr)(int, int, int, double *, double *);
  BoxOtherPtr box_other;
  typedef void (CommCAC::*BoxOtherFullPtr)(int, int, int, double *, double *);
  BoxOtherFullPtr box_other_full;
  void box_other_brick(int, int, int, double *, double *);
  void box_other_brick_full(int, int, int, double *, double *);
  void box_other_tiled(int, int, int, double *, double *);

  typedef int (CommCAC::*BoxTouchPtr)(int, int, int);
  BoxTouchPtr box_touch;
  int box_touch_brick(int, int, int);
  int box_touch_tiled(int, int, int);

  typedef int (CommCAC::*PointDropPtr)(int, double *);
  PointDropPtr point_drop;
  int point_drop_brick(int, double *);
  int point_drop_tiled(int, double *);
  int point_drop_tiled_recurse(double *, int, int);
  int closer_subbox_edge(int, double *);
 
  void overlap_element_comm(int);       // communicate element overlaps
  void get_aug_oboxes(int);       // communicate element overlaps
  void compute_eboxes(int);                //compute set of local and ghost eboxes as needed
  int  sendbox_include(int,int,int);            // decide if elements or atoms should be in the ghost sendbox
  int  pack_eboxes(int, int *, double *, int, int *,int); //pack element bounding boxes to send
  void unpack_eboxes(int, int, double *); //unpack element bounding boxes on recv
  void grow_send(int, int);            // reallocate send buffer
  void grow_eboxes(int);               // storage for element bounding boxes
  void grow_recv(int);                 // free/allocate recv buffer
  void grow_list(int, int, int);       // reallocate sendlist for one swap/proc
  void grow_sent_list(int, int, int);       // reallocate sent_flag for one swap/proc
  void overlap_grow_list(int, int, int);       // reallocate element overlap sendlist for one swap/proc
  void allocate_swap(int);             // allocate swap arrays
  void grow_swap_send(int, int, int);  // grow swap arrays for send and recv
  void grow_swap_recv(int, int, int);
  void deallocate_swap(int);           // deallocate swap arrays

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot yet use comm_style cac with triclinic box

Self-explanatory.

E: Cannot yet use comm_style cac with multi-mode comm

Self-explanatory

E: Cannot use comm_style cac with non CAC atom style

Self-explanatory.

E: Can only use comm style CAC with brick and rcb decompositions

Self-explanatory

E: Cannot use the CAC comm style without a CAC pair style

Self-explanatory

E: excessive/negative bin index

internal error check that should not occur unless a bug is present. Contact the Author.

E: Communication cutoff for comm_style cac cannot exceed periodic box length

Self-explanatory.

E: Reverse comm fix variable not yet supported by CommTiled

UNDOCUMENTED

E: Comm cac mis-match in box drop brick

Internal error check in comm_style cac which should not occur.
Contact the Author.

E: Comm cac invalid index in box drop brick

Internal error check in comm_style cac which should not occur.
Contact the Author.

U: KOKKOS package does not yet support comm_style cac

Self-explanatory.

*/
