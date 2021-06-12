// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MEMORY_KOKKOS_H
#define LMP_MEMORY_KOKKOS_H

#include "memory.h"             // IWYU pragma: export
#include "kokkos_type.h"

namespace LAMMPS_NS {

class MemoryKokkos : public Memory {
 public:
  MemoryKokkos(class LAMMPS *lmp) : Memory(lmp) {}

/* ----------------------------------------------------------------------
   Kokkos versions of create/grow/destroy multi-dimensional arrays
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type *&array,
                   int n1, const char *name)
{
  data = TYPE(name,n1);
  array = data.h_view.data();
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     typename TYPE::value_type *&array, int n1,
                     const char *name)
{
  data = TYPE(std::string(name),n1);
  h_data = Kokkos::create_mirror_view(data);
  array = h_data.data();
  return data;
}


template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
  h_data = Kokkos::create_mirror_view(data);
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 1d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type *&array,
                 int n1, const char *name)
{
  if (array == nullptr) return create_kokkos(data,array,n1,name);

  data.resize(n1);
  array = data.h_view.data();
  return data;
}

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type* &array)
{
  if (array == nullptr) return;
  data = TYPE();
  array = nullptr;
}

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE destroy_kokkos(TYPE &data)
{
  /*if (data.data()!=nullptr)
    free(data.data());*/
  data = TYPE();
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4,n5);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 , int n6 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*n6*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name) ,n1,n2,n3,n4,n5,n6);
  return data;
}



template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, int n1, int n2,
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  h_data = Kokkos::create_mirror_view(data);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array,
                   int n1, int n2, const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    if (n2==0)
      array[i] = nullptr;
    else
      array[i] = &data.h_view(i,0);
    n += n2;
  }
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     typename TYPE::value_type **&array, int n1, int n2,
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  h_data = Kokkos::create_mirror_view(data);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    if (n2==0)
      array[i] = nullptr;
    else
      array[i] = &h_data(i,0);
    n += n2;
  }
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array,
                 int n1, int n2, const char *name)
{
  if (array == nullptr) return create_kokkos(data,array,n1,n2,name);
  data.resize(n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type**) srealloc(array,nbytes,name);

  for (int i = 0; i < n1; i++)
    if (n2==0)
      array[i] = nullptr;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array,
                   int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  for (int i = 0; i < n1; i++)
    if (data.h_view.extent(1)==0)
      array[i] = nullptr;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array,
                 int n1, const char *name)
{
  if (array == nullptr) return create_kokkos(data,array,n1,name);

  data.resize(n1);

  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) srealloc(array,nbytes,name);

  for (int i = 0; i < n1; i++)
    if (data.h_view.extent(1)==0)
      array[i] = nullptr;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type** &array)
{
  if (array == nullptr) return;
  data = TYPE();
  sfree(array);
  array = nullptr;
}

};

}

#endif

