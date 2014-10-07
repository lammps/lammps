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

#ifndef CUDA_H
#define CUDA_H

#include "pointers.h"
#include "cuda_shared.h"
#include "cuda_data.h"
#include "cuda_precision.h"
#include <map>

#ifdef _DEBUG
#define MYDBG(a) a
#else
#define MYDBG(a)
#endif

namespace LAMMPS_NS
{
class Cuda : protected Pointers
{
  public:
    Cuda(class LAMMPS*);
    ~Cuda();
    //static void setDevice(class LAMMPS*);
    void allocate();

    void accelerator(int, char**);
    void activate();

    void setSharedDataZero();
    void setSystemParams();

    void setDomainParams();

    void checkResize();
    void evsetup_eatom_vatom(int eflag_atom, int vflag_atom);
    void uploadAll();
    void downloadAll();
    void upload(int datamask);
    void download(int datamask);
    void downloadX();

    class CudaNeighList* registerNeighborList(class NeighList* neigh_list);
    void uploadAllNeighborLists();
    void downloadAllNeighborLists();
    void set_neighinit(int dist_check, double triggerneighsq) {
      shared_data.atom.dist_check = dist_check;
      shared_data.atom.triggerneighsq = triggerneighsq;
    }
    bool decide_by_integrator() {
      return neighbor_decide_by_integrator  && cu_xhold && finished_setup;
    }
    void update_xhold(int &maxhold, double* xhold);

    void setTimingsZero();
    void print_timings();

    void cu_x_download() {
      cu_x->download();
    }
    bool device_set;
    bool dotiming;
    bool dotestatom;
    int testatom;

    double uploadtime, downloadtime;
    bool finished_setup, begin_setup;
    bool oncpu;
    bool finished_run;

    int self_comm;

    int cuda_exists;

    double extent[6];
    int* debugdata;
    // data shared between host code and device code
    // (number of atoms, device pointers for up- & download)
    cuda_shared_data shared_data;

    cCudaData<double  , F_CFLOAT , x >* cu_q;
    cCudaData<double  , F_CFLOAT , yx>* cu_f;
    cCudaData<double  , V_CFLOAT , x >* cu_mass;
    cCudaData<double  , V_CFLOAT , x >* cu_rmass;
    cCudaData<double  , V_CFLOAT , yx>* cu_v;
    cCudaData<double  , X_CFLOAT , yx>* cu_x;
    cCudaData<double  , X_CFLOAT , yx>* cu_xhold;
    cCudaData<int     , int     , x >* cu_mask;
    cCudaData<int     , int     , x >* cu_tag;
    cCudaData<int     , int     , x >* cu_type;
    cCudaData<int     , int     , x >* cu_image;
    cCudaData<double  , ENERGY_CFLOAT, x >* cu_eatom;
    cCudaData<double  , ENERGY_CFLOAT, yx>* cu_vatom;
    cCudaData<double  , ENERGY_CFLOAT, x >* cu_virial;
    cCudaData<double  , ENERGY_CFLOAT, x >* cu_eng_vdwl;
    cCudaData<double  , ENERGY_CFLOAT, x >* cu_eng_coul;
    cCudaData<double  , double  , x >* cu_extent;
    int* binned_id;
    cCudaData<int           , int            , xx >* cu_binned_id;
    int* binned_idnew;
    cCudaData<int           , int            , xx >* cu_binned_idnew;
    cCudaData<int           , int            , x >* cu_debugdata;
    cCudaData<double  , X_CFLOAT , x>* cu_radius;
    cCudaData<double  , F_CFLOAT , x>* cu_density;
    cCudaData<double  , V_CFLOAT , yx>* cu_omega;
    cCudaData<double  , F_CFLOAT , yx>* cu_torque;
    cCudaData<int           , int            , yx >* cu_special;
    cCudaData<int           , int            , yx >* cu_nspecial;
    cCudaData<int     , int     , x >* cu_molecule;


    cCudaData<X_CFLOAT  , X_CFLOAT , x>* cu_x_type;
    X_CFLOAT* x_type;

    cCudaData<V_CFLOAT  , V_CFLOAT , x>* cu_v_radius;
    V_CFLOAT* v_radius;

    cCudaData<V_CFLOAT  , V_CFLOAT , x>* cu_omega_rmass;
    V_CFLOAT* omega_rmass;

    cCudaData<int     , int     , x >* cu_map_array;
    int neighbor_decide_by_integrator;

    bool pinned;

    void* copy_buffer;
    int copy_buffersize;

  private:
    int pppn;                  // number of GPUs/node
    int *devicelist;           // IDs of GPUs

    std::map<class NeighList*, class CudaNeighList*> neigh_lists;
};
}

#endif // CUDA_H
