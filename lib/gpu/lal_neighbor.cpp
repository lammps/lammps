/***************************************************************************
                                 neighbor.cpp
                             -------------------
                            W. Michael Brown (ORNL)
                              Peng Wang (Nvidia)

  Class for handling neighbor lists

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov, penwang@nvidia.com
 ***************************************************************************/

#include "lal_precision.h"
#include "lal_neighbor.h"
#include "lal_device.h"
#include <cmath>
using namespace LAMMPS_AL;

int Neighbor::bytes_per_atom(const int max_nbors) const {
  if (_gpu_nbor==1)
    return (max_nbors+2)*sizeof(int);
  else if (_gpu_nbor==2)
    return (max_nbors+3)*sizeof(int);
  else if (_use_packing)
    return ((max_nbors+2)*2)*sizeof(int);
  else
    return (max_nbors+3)*sizeof(int);
}

bool Neighbor::init(NeighborShared *shared, const int inum,
                       const int host_inum, const int max_nbors,
                       const int maxspecial, UCL_Device &devi,
                       const int gpu_nbor, const int gpu_host,
                       const bool pre_cut, const int block_cell_2d,
                       const int block_cell_id, const int block_nbor_build,
                       const int threads_per_atom, const int warp_size,
                       const bool time_device,
                       const std::string compile_flags) {
  clear();

  _threads_per_atom=threads_per_atom;
  _block_cell_2d=block_cell_2d;
  _block_cell_id=block_cell_id;
  _max_block_nbor_build=block_nbor_build;
  _block_nbor_build=block_nbor_build;
  _warp_size=warp_size;
  _shared=shared;
  dev=&devi;
  _gpu_nbor=gpu_nbor;
  _time_device=time_device;
  if (gpu_host==0)
    _gpu_host=false;
  else if (gpu_host==1)
    _gpu_host=true;
  else
    // Not yet implemented
    assert(0==1);

  if (pre_cut || gpu_nbor==0)
    _alloc_packed=true;
  else
    _alloc_packed=false;

  if (pre_cut)
    _packed_permissions=UCL_READ_WRITE;
  else
    _packed_permissions=UCL_READ_ONLY;

  bool success=true;

  // Initialize timers for the selected GPU
  _nbor_time_avail=false;
  time_nbor.init(*dev);
  time_kernel.init(*dev);
  time_hybrid1.init(*dev);
  time_hybrid2.init(*dev);
  time_transpose.init(*dev);
  time_nbor.zero();
  time_kernel.zero();
  time_hybrid1.zero();
  time_hybrid2.zero();
  time_transpose.zero();

  _max_atoms=static_cast<int>(static_cast<double>(inum)*1.10);
  if (_max_atoms==0)
    _max_atoms=1000;

  _max_host=static_cast<int>(static_cast<double>(host_inum)*1.10);
  _max_nbors=(max_nbors/threads_per_atom+1)*threads_per_atom;

  _maxspecial=maxspecial;
  if (gpu_nbor==0)
    _maxspecial=0;

  if (gpu_nbor==0)
    success=success && (host_packed.alloc(2*IJ_SIZE,*dev,
                                          UCL_WRITE_ONLY)==UCL_SUCCESS);
  alloc(success);
  if (!success)
    return false;

  if (_use_packing==false)
    _shared->compile_kernels(devi,gpu_nbor,compile_flags);

  return success;
}

void Neighbor::alloc(bool &success) {
  dev_nbor.clear();
  host_acc.clear();
  int nt=_max_atoms+_max_host;
  if (_use_packing==false || _gpu_nbor>0)
    success=success &&
            (dev_nbor.alloc((_max_nbors+2)*_max_atoms,*dev)==UCL_SUCCESS);
  else
    success=success && (dev_nbor.alloc(3*_max_atoms,*dev,
                                       UCL_READ_ONLY)==UCL_SUCCESS);
  success=success && (host_acc.alloc(nt*2,*dev,
                                     UCL_READ_WRITE)==UCL_SUCCESS);

  _c_bytes=dev_nbor.row_bytes();
  if (_alloc_packed) {
    dev_packed.clear();
    success=success && (dev_packed.alloc((_max_nbors+2)*_max_atoms,*dev,
                                         _packed_permissions)==UCL_SUCCESS);
    dev_ilist.clear();
    success=success && (dev_ilist.alloc(_max_atoms,*dev,
                                      UCL_READ_WRITE)==UCL_SUCCESS);
    _c_bytes+=dev_packed.row_bytes()+dev_ilist.row_bytes();
  }
  if (_max_host>0) {
    nbor_host.clear();
    dev_numj_host.clear();
    host_ilist.clear();
    host_jlist.clear();

    success=(nbor_host.alloc(_max_nbors*_max_host,*dev,UCL_READ_WRITE,
                             UCL_READ_WRITE)==UCL_SUCCESS) && success;
    success=success && (dev_numj_host.alloc(_max_host,*dev,
                                            UCL_READ_WRITE)==UCL_SUCCESS);
    success=success && (host_ilist.alloc(nt,*dev,UCL_NOT_PINNED)==UCL_SUCCESS);
    if (!success)
      return;
    for (int i=0; i<nt; i++)
      host_ilist[i]=i;
    success=success && (host_jlist.alloc(_max_host,*dev,
                                         UCL_NOT_PINNED)==UCL_SUCCESS);
    if (!success)
      return;
    int *ptr=nbor_host.host.begin();
    for (int i=0; i<_max_host; i++) {
      host_jlist[i]=ptr;
      ptr+=_max_nbors;
    }
    _c_bytes+=nbor_host.device.row_bytes()+dev_numj_host.row_bytes();
  } else {
    // Some OpenCL implementations return errors for NULL pointers as args
    nbor_host.device.view(dev_nbor);
    dev_numj_host.view(dev_nbor);
  }
  if (_maxspecial>0) {
    dev_nspecial.clear();
    dev_special.clear();
    dev_special_t.clear();
    int at=_max_atoms+_max_host;
    success=success && (dev_nspecial.alloc(3*at,*dev,
                                           UCL_READ_ONLY)==UCL_SUCCESS);
    success=success && (dev_special.alloc(_maxspecial*at,*dev,
                                          UCL_READ_WRITE)==UCL_SUCCESS);
    success=success && (dev_special_t.alloc(_maxspecial*at,*dev,
                                            UCL_READ_ONLY)==UCL_SUCCESS);
    _gpu_bytes+=dev_nspecial.row_bytes()+dev_special.row_bytes()+
                dev_special_t.row_bytes();
  }

  _allocated=true;
}

void Neighbor::clear() {
  _gpu_bytes=0.0;
  _cell_bytes=0.0;
  _c_bytes=0.0;
  _bin_time=0.0;
  if (_ncells>0) {
    _ncells=0;
    cell_counts.clear();
    if (_gpu_nbor==2)
      delete [] cell_iter;
  }
  if (_allocated) {
    _allocated=false;
    _nbor_time_avail=false;

    host_packed.clear();
    host_acc.clear();
    dev_ilist.clear();
    dev_nbor.clear();
    nbor_host.clear();
    dev_packed.clear();
    dev_numj_host.clear();
    host_ilist.clear();
    host_jlist.clear();
    dev_nspecial.clear();
    dev_special.clear();
    dev_special_t.clear();

    time_kernel.clear();
    time_nbor.clear();
    time_hybrid1.clear();
    time_hybrid2.clear();
    time_transpose.clear();
  }
}

double Neighbor::host_memory_usage() const {
  if (_gpu_nbor>0) {
    if (_gpu_host)
      return nbor_host.device.row_bytes()*nbor_host.rows()+
             host_ilist.row_bytes()+host_jlist.row_bytes();
    else
      return 0;
  } else
    return host_packed.row_bytes()*host_packed.rows()+host_acc.row_bytes()+
           sizeof(Neighbor);
}

void Neighbor::get_host(const int inum, int *ilist, int *numj,
                        int **firstneigh, const int block_size) {
  _nbor_time_avail=true;
  time_nbor.start();

  UCL_H_Vec<int> ilist_view;
  ilist_view.view(ilist,inum,*dev);
  ucl_copy(dev_nbor,ilist_view,false);

  UCL_D_Vec<int> nbor_offset;
  UCL_H_Vec<int> host_offset;

  int copy_count=0;
  int ij_count=0;
  int acc_count=0;
  int dev_count=0;
  int *h_ptr=host_packed.begin();
  _nbor_pitch=inum;

  for (int ii=0; ii<inum; ii++) {
    int i=ilist[ii];
    int nj=numj[i];
    host_acc[ii]=nj;
    host_acc[ii+inum]=acc_count;

    acc_count+=nj;

    int *jlist=firstneigh[i];
    for (int jj=0; jj<nj; jj++) {
      *h_ptr=jlist[jj];
      h_ptr++;
      ij_count++;

      if (ij_count==IJ_SIZE) {
        dev_nbor.sync();
        host_offset.view_offset(IJ_SIZE*(copy_count%2),host_packed,IJ_SIZE);
        nbor_offset.view_offset(dev_count,dev_packed,IJ_SIZE);
        ucl_copy(nbor_offset,host_offset,true);
        copy_count++;
        ij_count=0;
        dev_count+=IJ_SIZE;
        h_ptr=host_packed.begin()+(IJ_SIZE*(copy_count%2));
      }
    }
  }
  if (ij_count!=0) {
    dev_nbor.sync();
    host_offset.view_offset(IJ_SIZE*(copy_count%2),host_packed,ij_count);
    nbor_offset.view_offset(dev_count,dev_packed,ij_count);
    ucl_copy(nbor_offset,host_offset,true);
  }
  UCL_D_Vec<int> acc_view;
  acc_view.view_offset(inum,dev_nbor,inum*2);
  ucl_copy(acc_view,host_acc,inum*2,true);

  UCL_H_Vec<int> host_view;
  host_view.alloc(_max_atoms,*dev,UCL_READ_WRITE);
  for (int ii=0; ii<inum; ii++) {
    int i=ilist[ii];
    host_view[i] = ii;
  }
  ucl_copy(dev_ilist,host_view,true);

  time_nbor.stop();

  if (_use_packing==false) {
    time_kernel.start();
    int GX=static_cast<int>(ceil(static_cast<double>(inum)*_threads_per_atom/
                                 block_size));
    _shared->k_nbor.set_size(GX,block_size);
    _shared->k_nbor.run(&dev_nbor, &dev_packed, &inum, &_threads_per_atom);
    time_kernel.stop();
  }
}

// This is the same as get host, but the requirement that ilist[i]=i and
// inum=nlocal is forced to be true to allow direct indexing of neighbors of
// neighbors
void Neighbor::get_host3(const int inum, const int nlist, int *ilist, int *numj,
                         int **firstneigh, const int block_size) {
  _nbor_time_avail=true;
  time_nbor.start();

  UCL_H_Vec<int> ilist_view;
  ilist_view.view(ilist,inum,*dev);
  ucl_copy(dev_nbor,ilist_view,false);

  UCL_D_Vec<int> nbor_offset;
  UCL_H_Vec<int> host_offset;

  int copy_count=0;
  int ij_count=0;
  int acc_count=0;
  int dev_count=0;
  int *h_ptr=host_packed.begin();
  _nbor_pitch=inum;

  if (nlist!=inum)
    host_acc.zero(inum);

  for (int ii=0; ii<nlist; ii++) {
    int i=ilist[ii];
    int nj=numj[i];
    host_acc[i]=nj;
    host_acc[i+inum]=acc_count;
    acc_count+=nj;
  }

  for (int i=0; i<inum; i++) {
    int nj=host_acc[i];
    int *jlist=firstneigh[i];
    for (int jj=0; jj<nj; jj++) {
      *h_ptr=jlist[jj];
      h_ptr++;
      ij_count++;

      if (ij_count==IJ_SIZE) {
        dev_nbor.sync();
        host_offset.view_offset(IJ_SIZE*(copy_count%2),host_packed,IJ_SIZE);
        nbor_offset.view_offset(dev_count,dev_packed,IJ_SIZE);
        ucl_copy(nbor_offset,host_offset,true);
        copy_count++;
        ij_count=0;
        dev_count+=IJ_SIZE;
        h_ptr=host_packed.begin()+(IJ_SIZE*(copy_count%2));
      }
    }
  }
  if (ij_count!=0) {
    dev_nbor.sync();
    host_offset.view_offset(IJ_SIZE*(copy_count%2),host_packed,ij_count);
    nbor_offset.view_offset(dev_count,dev_packed,ij_count);
    ucl_copy(nbor_offset,host_offset,true);
  }
  UCL_D_Vec<int> acc_view;
  acc_view.view_offset(inum,dev_nbor,inum*2);
  ucl_copy(acc_view,host_acc,inum*2,true);
  time_nbor.stop();

  if (_use_packing==false) {
    time_kernel.start();
    int GX=static_cast<int>(ceil(static_cast<double>(inum)*_threads_per_atom/
                                 block_size));
    _shared->k_nbor.set_size(GX,block_size);
    _shared->k_nbor.run(&dev_nbor, &dev_packed, &inum, &_threads_per_atom);
    time_kernel.stop();
  }
}

template <class numtyp, class acctyp>
void Neighbor::resize_max_neighbors(const int maxn, bool &success) {
  if (maxn>_max_nbors) {
    int mn=static_cast<int>(static_cast<double>(maxn)*1.10);
    mn=(mn/_threads_per_atom+1)*_threads_per_atom;
    success=success && (dev_nbor.resize((mn+1)*_max_atoms)==UCL_SUCCESS);
    _gpu_bytes=dev_nbor.row_bytes();
    if (_max_host>0) {
      success=success && (nbor_host.resize(mn*_max_host)==UCL_SUCCESS);
      int *ptr=nbor_host.host.begin();
      for (int i=0; i<_max_host; i++) {
        host_jlist[i]=ptr;
        ptr+=mn;
      }
      _gpu_bytes+=nbor_host.row_bytes();
    } else {
      nbor_host.device.view(dev_nbor);
      dev_numj_host.view(dev_nbor);
    }
    if (_alloc_packed) {
      success=success && (dev_packed.resize((mn+2)*_max_atoms)==UCL_SUCCESS);
      _gpu_bytes+=dev_packed.row_bytes();
    }
    _max_nbors=mn;
  }
}

template <class numtyp, class acctyp>
void Neighbor::build_nbor_list(double **x, const int inum, const int host_inum,
                               const int nall, Atom<numtyp,acctyp> &atom,
                               double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, bool &success,
                               int &mn) {
  _nbor_time_avail=true;
  const int nt=inum+host_inum;

  // Calculate number of cells and allocate storage for binning as necessary
  int ncellx, ncelly, ncellz, ncell_3d;
  int ghost_cells=2*_cells_in_cutoff;
  ncellx = static_cast<int>(ceil((subhi[0]-sublo[0])/_cell_size))+ghost_cells;
  ncelly = static_cast<int>(ceil((subhi[1]-sublo[1])/_cell_size))+ghost_cells;
  ncellz = static_cast<int>(ceil((subhi[2]-sublo[2])/_cell_size))+ghost_cells;
  ncell_3d = ncellx * ncelly * ncellz;
  if (ncell_3d+1>_ncells) {
    cell_counts.clear();

    if (_gpu_nbor==2) {
      if (_ncells>0)
        delete [] cell_iter;
      cell_iter = new int[ncell_3d+1];
      cell_counts.alloc(ncell_3d+1,dev_nbor,UCL_READ_WRITE,UCL_READ_ONLY);
    } else {
      cell_counts.device.clear();
      cell_counts.device.alloc(ncell_3d+1,dev_nbor);
    }

    _ncells=ncell_3d+1;
    _cell_bytes=cell_counts.device.row_bytes();
  }

  const numtyp cutoff_cast=static_cast<numtyp>(_cutoff);

  if (_maxspecial>0) {
    time_nbor.start();
    UCL_H_Vec<int> view_nspecial;
    UCL_H_Vec<tagint> view_special, view_tag;
    view_nspecial.view(nspecial[0],nt*3,*dev);
    view_special.view(special[0],nt*_maxspecial,*dev);
    view_tag.view(tag,nall,*dev);
    ucl_copy(dev_nspecial,view_nspecial,nt*3,false);
    ucl_copy(dev_special_t,view_special,nt*_maxspecial,false);
    ucl_copy(atom.dev_tag,view_tag,nall,false);
    time_nbor.stop();
    if (_time_device)
      time_nbor.add_to_total();
    time_transpose.start();
    const int b2x=_block_cell_2d;
    const int b2y=_block_cell_2d;
    const int g2x=static_cast<int>(ceil(static_cast<double>(_maxspecial)/b2x));
    const int g2y=static_cast<int>(ceil(static_cast<double>(nt)/b2y));
    _shared->k_transpose.set_size(g2x,g2y,b2x,b2y);
    _shared->k_transpose.run(&dev_special,&dev_special_t,&_maxspecial,&nt);
    time_transpose.stop();
  }

  // If binning on CPU, do this now
  if (_gpu_nbor==2) {
    double stime = MPI_Wtime();
    int *cell_id=atom.host_cell_id.begin();
    int *particle_id=atom.host_particle_id.begin();

    // Build cell list on CPU
    cell_counts.host.zero();
    double i_cell_size=1.0/_cell_size;

    int offset_hi=_cells_in_cutoff+1;
    for (int i=0; i<nt; i++) {
      double px, py, pz;
      px=x[i][0]-sublo[0];
      py=x[i][1]-sublo[1];
      pz=x[i][2]-sublo[2];

      int ix = static_cast<int>(px*i_cell_size+1);
      ix = std::max(ix,_cells_in_cutoff);
      ix = std::min(ix,ncellx-offset_hi);
      int iy = static_cast<int>(py*i_cell_size+1);
      iy = std::max(iy,_cells_in_cutoff);
      iy = std::min(iy,ncelly-offset_hi);
      int iz = static_cast<int>(pz*i_cell_size+1);
      iz = std::max(iz,_cells_in_cutoff);
      iz = std::min(iz,ncellz-offset_hi);

      int id = ix+iy*ncellx+iz*ncellx*ncelly;
      cell_id[i] = id;
      cell_counts[id+1]++;
    }

    for (int i=nt; i<nall; i++) {
      double px, py, pz;
      px=x[i][0]-sublo[0];
      py=x[i][1]-sublo[1];
      pz=x[i][2]-sublo[2];

      int ix = static_cast<int>(px*i_cell_size+1);
      ix = std::max(ix,0);
      ix = std::min(ix,ncellx-1);
      int iy = static_cast<int>(py*i_cell_size+1);
      iy = std::max(iy,0);
      iy = std::min(iy,ncelly-1);
      int iz = static_cast<int>(pz*i_cell_size+1);
      iz = std::max(iz,0);
      iz = std::min(iz,ncellz-1);

      int id = ix+iy*ncellx+iz*ncellx*ncelly;
      cell_id[i] = id;
      cell_counts[id+1]++;
    }

    mn=0;
    for (int i=0; i<_ncells; i++)
      mn=std::max(mn,cell_counts[i]);
    mn*=8;
    set_nbor_block_size(mn/2);

    resize_max_neighbors<numtyp,acctyp>(mn,success);
    if (!success)
      return;
    _total_atoms=nt;

    cell_iter[0]=0;
    for (int i=1; i<_ncells; i++) {
      cell_counts[i]+=cell_counts[i-1];
      cell_iter[i]=cell_counts[i];
    }
    time_hybrid1.start();
    cell_counts.update_device(true);
    time_hybrid1.stop();
    for (int i=0; i<nall; i++) {
      int celli=cell_id[i];
      int ploc=cell_iter[celli];
      cell_iter[celli]++;
      particle_id[ploc]=i;
    }
    time_hybrid2.start();
    ucl_copy(atom.dev_particle_id,atom.host_particle_id,true);
    time_hybrid2.stop();
    _bin_time+=MPI_Wtime()-stime;
  }

  time_kernel.start();

  _nbor_pitch=inum;
  _shared->neigh_tex.bind_float(atom.x,4);

  // If binning on GPU, do this now
  if (_gpu_nbor==1) {
    const numtyp i_cell_size=static_cast<numtyp>(1.0/_cell_size);
    const int neigh_block=_block_cell_id;
    const int GX=(int)ceil((float)nall/neigh_block);
    const numtyp sublo0=static_cast<numtyp>(sublo[0]);
    const numtyp sublo1=static_cast<numtyp>(sublo[1]);
    const numtyp sublo2=static_cast<numtyp>(sublo[2]);
    _shared->k_cell_id.set_size(GX,neigh_block);
    _shared->k_cell_id.run(&atom.x, &atom.dev_cell_id,
                           &atom.dev_particle_id, &sublo0, &sublo1,
                           &sublo2, &i_cell_size, &ncellx, &ncelly, &ncellz,
                           &nt, &nall, &_cells_in_cutoff);

    atom.sort_neighbor(nall);

    /* calculate cell count */
    _shared->k_cell_counts.set_size(GX,neigh_block);
    _shared->k_cell_counts.run(&atom.dev_cell_id, &cell_counts, &nall,
                               &ncell_3d);
  }

  /* build the neighbor list */
  const int cell_block=_block_nbor_build;
  _shared->k_build_nbor.set_size(ncellx-ghost_cells,(ncelly-ghost_cells)*
                                 (ncellz-ghost_cells),cell_block,1);
  _shared->k_build_nbor.run(&atom.x, &atom.dev_particle_id,
                            &cell_counts, &dev_nbor, &nbor_host,
                            &dev_numj_host, &_max_nbors, &cutoff_cast, &ncellx,
                            &ncelly, &ncellz, &inum, &nt, &nall,
                            &_threads_per_atom, &_cells_in_cutoff);

  /* Get the maximum number of nbors and realloc if necessary */
  UCL_D_Vec<int> numj;
  numj.view_offset(inum,dev_nbor,inum);
  ucl_copy(host_acc,numj,inum,true);
  if (nt>inum) {
    UCL_H_Vec<int> host_offset;
    host_offset.view_offset(inum,host_acc,nt-inum);
    ucl_copy(host_offset,dev_numj_host,nt-inum,true);
  }

  if (_gpu_nbor!=2) {
    host_acc.sync();
    mn=host_acc[0];
    for (int i=1; i<nt; i++)
      mn=std::max(mn,host_acc[i]);
    set_nbor_block_size(mn);

    if (mn>_max_nbors) {
      resize_max_neighbors<numtyp,acctyp>(mn,success);
      if (!success)
        return;
      time_kernel.stop();
      if (_time_device)
        time_kernel.add_to_total();
      build_nbor_list(x, inum, host_inum, nall, atom, sublo, subhi, tag,
                      nspecial, special, success, mn);
      return;
    }
  }

  if (_maxspecial>0) {
    const int GX2=static_cast<int>(ceil(static_cast<double>
                                          (nt*_threads_per_atom)/cell_block));
    _shared->k_special.set_size(GX2,cell_block);
    _shared->k_special.run(&dev_nbor, &nbor_host, &dev_numj_host,
                           &atom.dev_tag, &dev_nspecial, &dev_special,
                           &inum, &nt, &_max_nbors, &_threads_per_atom);
  }
  time_kernel.stop();

  time_nbor.start();
  if (inum<nt) {
    nbor_host.update_host(true);
    nbor_host.sync();
  }
  time_nbor.stop();
}

template void Neighbor::build_nbor_list<PRECISION,ACC_PRECISION>
     (double **x, const int inum, const int host_inum, const int nall,
      Atom<PRECISION,ACC_PRECISION> &atom, double *sublo, double *subhi,
      tagint *, int **, tagint **, bool &success, int &mn);

