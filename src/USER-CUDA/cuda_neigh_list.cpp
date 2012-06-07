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

#include "cuda_neigh_list.h"
#include "neigh_list.h"
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include "cuda.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

CudaNeighList::CudaNeighList(LAMMPS *lmp, class NeighList* neigh_list) : Pointers(lmp)
{
        cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        MYDBG(printf("# CUDA: CudaNeighList::cudaNeighList() ... start\n");)
        this->neigh_list = neigh_list;
        neigh_list->cuda_list=this;
        sneighlist.maxlocal = neigh_list->get_maxlocal();
        sneighlist.maxneighbors = 32;
        sneighlist.maxcut = 0.0;
        sneighlist.cutneighsq = NULL;
        cu_neighbors = NULL;
        cu_neighbors_border = NULL;
        cu_neighbors_inner = NULL;
        cu_numneigh_border = NULL;
        cu_numneigh_inner = NULL;
        cu_numneigh = NULL;
        cu_ilist = NULL;
        cu_ilist_border = NULL;
        cu_inum_border = NULL;
        inum_border = 0;
        neighbors = NULL;
        neighbors_inner = NULL;
        neighbors_border = NULL;
        numneigh_border = NULL;
        numneigh_inner = NULL;
        ilist_border = NULL;

        build_cuda = false;
        sneighlist.binned_id=NULL;
        sneighlist.bin_dim=new int[3];
        sneighlist.bin_dim[0]=0;
        sneighlist.bin_dim[1]=0;
        sneighlist.bin_dim[2]=0;

        cu_ex_type = NULL;
        cu_ex1_bit = NULL;
        cu_ex2_bit = NULL;
        cu_ex_mol_bit = NULL;
        sneighlist.nex_type=0;
        sneighlist.nex_group=0;
        sneighlist.nex_mol=0;

        sneighlist.bin_nmax=0;
        sneighlist.bin_extraspace=0.05;
        MYDBG(printf("# CUDA: CudaNeighList::cudaNeighList() ... end\n");)

}

CudaNeighList::~CudaNeighList()
{
        dev_free();
}

void CudaNeighList::dev_alloc()
{
        MYDBG( printf("# CUDA: CudaNeighList::dev_alloc() ... start\n"); )
        cu_ilist         = new cCudaData<int , int , x> (neigh_list->ilist   , & sneighlist.ilist     , sneighlist.maxlocal );
        cu_numneigh      = new cCudaData<int , int , x> (neigh_list->numneigh, & sneighlist.numneigh  , sneighlist.maxlocal );
        neighbors = new int[atom->nmax*sneighlist.maxneighbors];
        cu_neighbors= new cCudaData<int, int, x> (neighbors                                          , & sneighlist.neighbors, atom->nmax*sneighlist.maxneighbors );

        if(cuda->shared_data.overlap_comm)
        {
        ilist_border  = new int[sneighlist.maxlocal];
        numneigh_border        = new int[sneighlist.maxlocal];
        numneigh_inner        = new int[sneighlist.maxlocal];
        cu_inum_border  = new cCudaData<int , int , x> (&inum_border                 , & sneighlist.inum_border      , 1 );
        cu_ilist_border  = new cCudaData<int , int , x> (ilist_border                 , & sneighlist.ilist_border     , sneighlist.maxlocal );
        cu_numneigh_border        = new cCudaData<int , int , x> (numneigh_border  , & sneighlist.numneigh_border  , sneighlist.maxlocal );
        cu_numneigh_inner         = new cCudaData<int , int , x> (numneigh_inner   , & sneighlist.numneigh_inner   , sneighlist.maxlocal );
        neighbors_border = new int[sneighlist.maxlocal*sneighlist.maxneighbors];
        cu_neighbors_border= new cCudaData<int, int, x> (neighbors_border         , & sneighlist.neighbors_border, sneighlist.maxlocal*sneighlist.maxneighbors );
        neighbors_inner = new int[sneighlist.maxlocal*sneighlist.maxneighbors];
        cu_neighbors_inner = new cCudaData<int, int, x> (neighbors_inner         , & sneighlist.neighbors_inner , sneighlist.maxlocal*sneighlist.maxneighbors );
        }
        cuda->shared_data.atom.update_neigh=2;
        MYDBG( printf("# CUDA: CudaNeighList::dev_alloc() ... end\n"); )
}

void CudaNeighList::dev_free()
{
        MYDBG( printf("# CUDA: CudaNeighList::dev_free() ... start\n"); )
        delete cu_numneigh;
        delete cu_ilist;
        delete [] neighbors;
        delete cu_neighbors;

        if(cuda->shared_data.overlap_comm)
        {
        delete [] ilist_border;
        delete [] numneigh_border;
        delete [] numneigh_inner;
        delete [] neighbors_border;
        delete [] neighbors_inner;
        delete cu_inum_border;
        delete cu_neighbors_border;
        delete cu_neighbors_inner;
        delete cu_numneigh_border;
        delete cu_numneigh_inner;
        delete cu_ilist_border;
        }
        MYDBG( printf("# CUDA: CudaNeighList::dev_free() ... end\n"); )
}

void CudaNeighList::grow_device()
{
        MYDBG(printf("# CUDA: CudaNeighList::grow_device() ... start\n");)
        // if host has allocated more memory for atom arrays than device has, then allocate more memory on device
        int new_maxlocal = neigh_list->get_maxlocal();
        if(sneighlist.maxlocal < new_maxlocal)
        {
                sneighlist.maxlocal = new_maxlocal;
                dev_free();
                dev_alloc();
        }

        if(!cu_ilist || !cu_numneigh) dev_alloc();

        // check, if hosts data has been allocated somewhere else
        if(cu_ilist   ->get_host_data() != neigh_list->ilist)    cu_ilist   ->set_host_data(neigh_list->ilist);
        if(cu_numneigh->get_host_data() != neigh_list->numneigh) cu_numneigh->set_host_data(neigh_list->numneigh);

        MYDBG(printf("# CUDA: CudaNeighList::grow_device() ... end\n");)
}


void CudaNeighList::nl_upload(bool will_be_changed)
{
        //return;
        MYDBG(printf("# CUDA: CudaNeighList::nl_upload() ... start\n");)
        if(cu_ilist)
        cu_ilist->upload();
        if(cu_numneigh)
        cu_numneigh->upload();
        MYDBG(printf("# CUDA: CudaNeighList::nl_upload() ... end\n");)
}

void CudaNeighList::nl_download(bool will_be_changed)
{
        MYDBG(printf("# CUDA: CudaNeighList::nl_download() ... start\n");)
        if(cu_ilist)
        cu_ilist->download();
        if(cu_numneigh)
        cu_numneigh->download();
        MYDBG(printf("# CUDA: CudaNeighList::nl_download() ... end\n");)
}
