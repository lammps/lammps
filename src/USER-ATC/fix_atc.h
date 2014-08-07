#ifdef FIX_CLASS

FixStyle(atc,FixATC)

#else

#ifndef FIX_ATC_H
#define FIX_ATC_H

#include "fix.h"
#include "pointers.h" // access to lammps pointers

namespace ATC {
  class ATC_Method;
}
namespace LAMMPS_NS {
  class  NeighList;

  /**
   *  @class FixATC 
   *  @brief Class for an atom-to-continuum (ATC) Lammps fix. 
   */

  class FixATC : public Fix {
  public: 
    /** constructor & destructor */
    FixATC(class LAMMPS *, int, char **);
    ~FixATC();

    /** initialization functions */
    void init();
    void init_list(int id, NeighList *ptr) ; 
    void setup(int vflag);
    void min_setup(int vflag);

    /** setmask: tell LAMMPS which fix methods to call */
    int setmask();

    /** initial_integrate */
    void initial_integrate(int vflag);

    /** after first integrate phase */
    void post_integrate();

    /** final_integrate */
    void final_integrate();

    /** end of step for run or minimize */
    void end_of_step();

    /** pre_exchange is used to modify fix-specific data
       and is called before domain->pbc() and comm->exchange().  */
    void setup_pre_exchange();
    void pre_exchange();
    void min_setup_pre_exchange();
    void min_pre_exchange();

    double memory_usage(); 
    void grow_arrays(int);
    void copy_arrays(int, int, int);

    /** pack_exchange called from atom_vec->pack_exchange()
       and packs fix-specific data for a given real (local) 
       atom being sent to another processor.  */
    int pack_exchange(int, double *);

    /** unpack_exchange called from atom_vec->unpack_exchange()
       and unpacks fix-specific data for a given real (local) 
       atom received from another processor. */
    int unpack_exchange(int, double *);

    /** pack_forward_comm called from comm->forward_comm_fix and
       packs fix-specific data for a given ghost atom
       from exchange with another proc */
    int pack_forward_comm(int , int *, double *, int, int *);  
 
    /** unpack_comm called from comm->forward_comm_fix and
       unpacks fix-specific data for a given ghost atom
       from exchange with another proc */
    void unpack_forward_comm(int, int, double *);

    /** pre_neighbor is used to modify fix-specific data
       and is called before neighbor list is built in 
       neighbor->build().  */
    void pre_neighbor();
    void setup_pre_neighbor();
    void min_setup_pre_neighbor();

    /** pre/post_force is used to modify fix-specific data
        and is before/after the various force computations. */
    void pre_force(int vflag);
    void post_force(int vflag);

    /** post_run is called after a run completes */
    void post_run();

    /** min_pre_force is called before forces are calculated in minimize */
    void min_pre_force(int vflag);

    /** min_post_force is called after forces are calculated in minimize */
    void min_post_force(int vflag);

    /** modify atc parameters (parser) */
    int modify_param(int narg, char** arg);

    /** calls ATC_Method to handle restarting/checkpointing */
    /** these four methods are for writing per-atom quantities */
    int pack_restart(int, double *);
    void unpack_restart(int, int);
    int size_restart(int);
    int maxsize_restart();
    /** these two methods are for writing all other quantities */
    void write_restart(FILE *);
    void restart(char *);

    /** accessor function for ATC_Method class pointer */
    const ATC::ATC_Method* atc() { return atc_; }

  protected:
    LAMMPS * lammps_;

    /** functions for "thermo" output */
    virtual double compute_scalar() ;
    virtual double compute_vector(int n) ;
    virtual double compute_array(int irow, int icol) ;
    double dtv,dtf;
    ATC::ATC_Method *atc_;
  };
}

#endif
#endif
