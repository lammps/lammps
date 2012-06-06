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

#include "neighbor_cuda.h"
#include "neigh_list.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "cuda_neigh_list.h"
#include "cuda.h"
#include "neighbor_cu.h"
#include <cmath>
using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   N^2 search for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */
void NeighborCuda::full_bin_cuda(NeighList *list)
{
  MYDBG(printf(" # CUDA::NeighFullBinCuda ... start\n");)
  if(includegroup) error->warning(FLERR,"Warning using inlcudegroup neighborbuild. This is not yet supported by CUDA neighborbuild styles.\n");
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if(nlocal==0) return;
  CudaNeighList* clist=list->cuda_list;
  cuda_shared_neighlist* slist=&clist->sneighlist;

  if(not clist) cuda->registerNeighborList(list);

  clist->build_cuda=true;

  if(slist->bin_extraspace<0.09)
  {
    for(int i=1;i<=atom->ntypes;i++)
    for(int j=1;j<=atom->ntypes;j++)
    {
            if(slist->maxcut<cutneighsq[i][j]) slist->maxcut=cutneighsq[i][j];
    }
    slist->maxcut=sqrt(slist->maxcut);
  }
  int bin_dim_tmp[3];
  int bin_nmax_tmp;
//printf("Hallo\n");
  timespec starttime,endtime;
  do
  {
    do
    {
      bin_dim_tmp[0]=static_cast <int> ((domain->subhi[0]-domain->sublo[0])/slist->maxcut);
      bin_dim_tmp[1]=static_cast <int> ((domain->subhi[1]-domain->sublo[1])/slist->maxcut);
      bin_dim_tmp[2]=static_cast <int> ((domain->subhi[2]-domain->sublo[2])/slist->maxcut);
      if(bin_dim_tmp[0]==0) bin_dim_tmp[0]+=1;
      if(bin_dim_tmp[1]==0) bin_dim_tmp[1]+=1;
      if(bin_dim_tmp[2]==0) bin_dim_tmp[2]+=1;
      bin_nmax_tmp=static_cast <int> ((1.0+slist->bin_extraspace)*nlocal/(bin_dim_tmp[0]*bin_dim_tmp[1]*bin_dim_tmp[2]));
      bin_dim_tmp[0]+=4;
      bin_dim_tmp[1]+=4;
      bin_dim_tmp[2]+=4;
           if(bin_nmax_tmp<32) slist->maxcut*=1.2;
          // printf("slist->maxcut: %lf\n", slist->maxcut);
    } while(bin_nmax_tmp<32);
    if((slist->bin_dim[0]!=bin_dim_tmp[0])||(slist->bin_dim[1]!=bin_dim_tmp[1])||(slist->bin_dim[2]!=bin_dim_tmp[2])||(slist->bin_nmax!=bin_nmax_tmp))
    {
            if(slist->binned_id!=NULL)
            CudaWrapper_FreeCudaData(slist->binned_id,slist->bin_dim[0]*slist->bin_dim[1]*slist->bin_dim[2]*slist->bin_nmax*sizeof(int));
            slist->bin_dim[0] = bin_dim_tmp[0];
            slist->bin_dim[1] = bin_dim_tmp[1];
            slist->bin_dim[2] = bin_dim_tmp[2];
            slist->bin_nmax = bin_nmax_tmp;
            slist->binned_id=(int*) CudaWrapper_AllocCudaData(slist->bin_dim[0]*slist->bin_dim[1]*slist->bin_dim[2]*slist->bin_nmax*sizeof(int));
           //printf("slist->bin: %i %i %i %i \n", bin_dim_tmp[0],bin_dim_tmp[1],bin_dim_tmp[2],bin_nmax_tmp);
    }
    //if(list->cuda_list->sneighlist.bin_nmax>512) error->all(FLERR,"To many atoms per bin. Likely cause is very long pair cutoff. This needs major rewrite of code and is not yet scheduled to be done.\n");
  }while(Cuda_BinAtoms(&cuda->shared_data, &list->cuda_list->sneighlist));

 // cuda->cu_debugdata->memset_device(0);
  int maxneighbors=slist->maxneighbors;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;

  if((nex_type!=slist->nex_type)||
  (nex_group!=slist->nex_group)||
  (nex_mol!=slist->nex_mol))
  {
          slist->nex_type=nex_type;
          slist->nex_group=nex_group;
          slist->nex_mol=nex_mol;
          //printf("%i %i %i\n",nex_type,nex_group,nex_mol);
          if(nex_type)
          {
          delete clist->cu_ex_type;
          clist->cu_ex_type=new cCudaData<int , int , x> (&ex_type[0][0]   , & slist->ex_type     , (atom->ntypes+1)*(atom->ntypes+1) );
          clist->cu_ex_type->upload();
          }
         //printf("AA %i %i %i\n",nex_type,nex_group,nex_mol);
          if(nex_group)
          {
           delete clist->cu_ex1_bit;
          clist->cu_ex1_bit=new cCudaData<int , int , x> (ex1_bit   , & slist->ex1_bit     , nex_group );
          clist->cu_ex1_bit->upload();
          //printf("A %i %i %i\n",nex_type,nex_group,nex_mol);
          delete clist->cu_ex2_bit;
          clist->cu_ex2_bit=new cCudaData<int , int , x> (ex2_bit   , & slist->ex2_bit     , nex_group );
          clist->cu_ex2_bit->upload();
          }
          //printf("B %i %i %i\n",nex_type,nex_group,nex_mol);
          if(nex_mol)
          {
          delete clist->cu_ex_mol_bit;
          clist->cu_ex_mol_bit=new cCudaData<int , int , x> (ex_mol_bit   , & slist->ex_mol_bit     , nex_mol );
          clist->cu_ex_mol_bit->upload();
          }
          //printf("C %i %i %i\n",nex_type,nex_group,nex_mol);
  }
  int overflow = 0;
  int inum = 0;
  int npnt = 0;
  do
  {
          npnt=0;
          inum=0;
    overflow=0;
    clist->grow_device();
    slist->cutneighsq=cutneighsq;
    slist->maxneighbors=maxneighbors;
    slist->inum = list->inum = nlocal;
    //list->cuda_list->grow_device();
    if(cuda->shared_data.overlap_comm)
    {
          list->cuda_list->inum_border=0;
          list->cuda_list->cu_inum_border->upload();
    }

    cuda->shared_data.atom.nall=nall;
    //Cuda_NeighborReBuildFirstneigh(&cuda->shared_data, &list->cuda_list->sneighlist);
    overflow= Cuda_NeighborBuildFullBin(&cuda->shared_data, &list->cuda_list->sneighlist);

        /*cuda->cu_debugdata->download();
        printf("Debugdata: %i ",cuda->debugdata[0]);
        for(int i=0;i<cuda->debugdata[0];i+=3) printf("// %i %i %i",cuda->debugdata[i+1],cuda->debugdata[i+2],cuda->debugdata[i+3]);
        printf("\n");*/
        //printf("maxneighborsA: %i %i %i %i\n",maxneighbors,pgsize,oneatom,atom->nmax);

    if(overflow<0)
    {
            maxneighbors+=32;
            if(-overflow>maxneighbors) maxneighbors=((-overflow+37)/32)*32;
            delete list->cuda_list->cu_neighbors;
            delete [] list->cuda_list->neighbors;
            list->cuda_list->neighbors= new int[slist->maxlocal*maxneighbors];
            list->cuda_list->sneighlist.maxneighbors=maxneighbors;
        //printf("maxneighborsA1: %i %i %i %i %i\n",maxneighbors,pgsize,oneatom,atom->nmax,slist->maxlocal);
            list->cuda_list->cu_neighbors= new cCudaData<int, int, x> (list->cuda_list->neighbors                          , & list->cuda_list->sneighlist.neighbors, slist->maxlocal*maxneighbors );
        //printf("maxneighborsA2: %i %i %i %i\n",maxneighbors,pgsize,oneatom,atom->nmax);

            if(cuda->shared_data.overlap_comm)
            {
              list->cuda_list->sneighlist.maxneighbors=maxneighbors;
              list->cuda_list->dev_free();
              list->cuda_list->dev_alloc();
            }
        //printf("maxneighborsA3: %i %i %i %i\n",maxneighbors,pgsize,oneatom,atom->nmax);
    }
        //printf("maxneighborsB: %i %i %i %i\n",maxneighbors,pgsize,oneatom,atom->nmax);
    if(cuda->shared_data.overlap_comm)
    {
                  list->cuda_list->cu_inum_border->download();
                  list->cuda_list->sneighlist.inum_border2=list->cuda_list->inum_border;
    }
  }
  while(overflow<0);

  //cuda->cu_debugdata->download();
 // printf("Differences in: %i\n",cuda->debugdata[0]);
 // for(int i=0;i<20;i++) printf("%i %i %i %i// ",cuda->debugdata[4*i+1],cuda->debugdata[4*i+2],cuda->debugdata[4*i+3],cuda->debugdata[4*i+4]);
//  printf("\n");
/*for(int i=0;i<10;i++)
{
        printf("%i %i // ",i,numneigh[i]);
        for(int j=0;j<numneigh[i];j++)
         printf("%i ",list->cuda_list->neighbors[i+j*nlocal]);
        printf("\n");
}*/
/*  int count=0;
  if(cuda->shared_data.overlap_comm)
  {
  list->cuda_list->cu_inum_border->download();
  list->cuda_list->cu_ilist_border->download();
  list->cuda_list->cu_numneigh_border->download();
  list->cuda_list->cu_numneigh_inner->download();
  list->cuda_list->cu_neighbors->download();
  list->cuda_list->cu_neighbors_inner->download();
  list->cuda_list->cu_neighbors_border->download();

  //list->cuda_list->cu_firstneigh->download();
 // list->cuda_list->nl_download();
  list->cuda_list->cu_numneigh->download();
  int diff=0;
  //for(int i=0;i<nlocal;i++)*/
 /* int i=123;
  {
          int k=-1;
          //printf("inum_border: %i\n",list->cuda_list->inum_border);
          //for(int j=0;j<list->numneigh[i];j++) printf("%i ",list->firstneigh[i][j]);printf("\n");
          for(int j=0;j<list->cuda_list->inum_border;j++)
          if(list->cuda_list->ilist_border[j]==i) k=j;
          int d=numneigh[i]-list->cuda_list->numneigh_inner[i];
          if(k>-1) d-=list->cuda_list->numneigh_border[k];
          if(d!=0) {printf("Error at %i %i %i %i %i\n",i,k,d,numneigh[i],list->cuda_list->numneigh_inner[i]); diff++;}
          if(k>-1 && count<10)
          {
                  printf("Numneighs: %i %i %i  Border_i: %i %i\n",numneigh[i],list->cuda_list->numneigh_inner[i],list->cuda_list->numneigh_border[k],k,(int)list->cuda_list->cu_ilist_border->dev_data());
        cuda->shared_data.me=k;
        for(int j=0;j<numneigh[i];j++)
         printf("%i ",list->cuda_list->neighbors[i+j*nlocal]);
           printf("\n");
        for(int j=0;j<list->cuda_list->numneigh_inner[i];j++)
         printf("%i ",list->cuda_list->neighbors_inner[i+j*nlocal]);
         printf(" // ");
        for(int j=0;j<list->cuda_list->numneigh_border[k];j++)
         printf("%i ",list->cuda_list->neighbors_border[k+j*nlocal]);
           printf("\n");
           count++;
          }
  }
  printf("%i\n",diff);
  }*/
  list->cuda_list->cu_numneigh->download();
  list->cuda_list->cu_ilist->download();
  cuda->shared_data.atom.update_neigh=2;
        //printf("Done\n");

  MYDBG(printf(" # CUDA::NeighFullBinCuda ... end\n");)

}


void NeighborCuda::full_nsq_cuda(NeighList *list)
{
        printf("Full_Nsq cuda neighbor list build is not implemented anymore.\n");
return;
/*
  MYDBG(printf(" # CUDA::NeighFullNSQCuda ... start\n");)
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if(cuda->cu_xhold) cuda->cu_xhold->upload();


  if(not list->cuda_list) cuda->registerNeighborList(list);
  list->cuda_list->build_cuda=true;
  int maxneighbors=list->cuda_list->sneighlist.maxneighbors;
  int neigh_lists_per_page=pgsize/maxneighbors;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;

  int overflow = 0;
  int inum = 0;
  int npage = 0;
  int npnt = 0;
  do
  {
          npage=0;
          npnt=0;
          inum=0;
    overflow=0;
          neigh_lists_per_page=pgsize/maxneighbors;
    npage=(2*nlocal*maxneighbors-1)/pgsize;
    while(npage>list->maxpage) list->add_pages();
    pages = list->pages;
    npage=0;
          list->cuda_list->sneighlist.neigh_lists_per_page=pgsize/maxneighbors;
    list->cuda_list->grow_device();
    list->cuda_list->sneighlist.cutneighsq=cutneighsq;
    list->cuda_list->sneighlist.maxneighbors=maxneighbors;
    list->cuda_list->sneighlist.inum = list->inum = nlocal;

    cuda->shared_data.atom.nall=nall;
    Cuda_NeighborReBuildFirstneigh(&cuda->shared_data, &list->cuda_list->sneighlist);
    overflow= not Cuda_NeighborBuildFullNsq(&cuda->shared_data, &list->cuda_list->sneighlist);



     if(overflow) maxneighbors+=32;
  }
  while(overflow);
   if(not cudable) list->cuda_list->nl_download();
  MYDBG(printf(" # CUDA::NeighFullNSQCuda ... end\n");)
  */
}
