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

#ifdef FIX_CLASS

FixStyle(atc,FixATC)

#else

#ifndef LMP_FIX_ATC_H
#define LMP_FIX_ATC_H

#include "fix.h"

#include "pointers.h" // access to lammps pointers

// NOTE what is the correct way?
//class ATC::ATC_Transfer;
#include "ATC_Transfer.h"
#include "LammpsInterface.h"

namespace LAMMPS_NS {
  // fwd decl
  class  NeighList;

  class FixATC : public Fix {
  public: 
    /** constructor & destructor */
    FixATC(class LAMMPS *, int, char **);
    ~FixATC();

    /** calls ATC_Transfer */ 
    void init();
    void init_list(int id, NeighList *ptr) {
      ATC::LammpsInterface::instance()->init_list(id,ptr);
    }
    void setup(int vflag);

    /** setmask: tell LAMMPS which fix methods to call */
    int setmask();

    /** initial_integrate: calls ATC_Transfer */
    void initial_integrate(int vflag);

    /** final_integrate: calls ATC_Transfer */
    void final_integrate();

    /** calls ATC_Transfer */
    void pre_exchange();
    void min_pre_exchange();
    double memory_usage();
    void grow_arrays(int);
    void copy_arrays(int, int);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int pack_comm(int , int *, double *, int, int *);  
    void unpack_comm(int, int, double *);

    /** modify atc parameters (parser) */
    int modify_param(int narg, char** arg);

    /** calls ATC_Transfer to handle restarting/checkpointing */
    /** these four methods are for writing per-atom quantities */
    int pack_restart(int, double *);
    void unpack_restart(int, int);
    int size_restart(int);
    int maxsize_restart();
    /** these two methods are for writing all other quantities */
    void write_restart(FILE *);
    void restart(char *);

  protected:
    /** functions for "thermo" output */
    virtual double compute_scalar() {return atcTransfer_->compute_scalar();}
    virtual double compute_vector(int n) {return atcTransfer_->compute_vector(n);}
    double dtv,dtf;
    ATC::ATC_Transfer *atcTransfer_;
  };

}

#endif
#endif
