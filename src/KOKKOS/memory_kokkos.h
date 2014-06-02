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
  array = data.h_view.ptr_on_device();
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, 
                     typename TYPE::value_type *&array, int n1, 
                     const char *name)
{
  data = TYPE(std::string(name),n1);
#ifndef KOKKOS_USE_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  array = h_data.ptr_on_device();
  return data;
}


template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
#ifndef KOKKOS_USE_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
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
  if (array == NULL) return create_kokkos(data,array,n1,name);
  
  data.resize(n1);
  array = data.h_view.ptr_on_device();
  return data;
}

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type* &array)
{
  if (array == NULL) return;
  data = TYPE();
  array = NULL;
}

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, const char *name)
{
  data = TYPE(name,n1,n2);
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, int n1, int n2, 
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
#ifndef KOKKOS_USE_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
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
#ifndef KOKKOS_USE_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  bigint n = 0;
  for (int i = 0; i < n1; i++) {
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
  if (array == NULL) return create_kokkos(data,array,n1,n2,name);
  data.resize(n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type**) srealloc(array,nbytes,name);
  
  for (int i = 0; i < n1; i++)
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
    array[i] = &data.h_view(i,0);
  
  return data;
}

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                 int n1, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,name);
  
  data.resize(n1);
  
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  for (int i = 0; i < n1; i++)
    array[i] = &data.h_view(i,0);
  
  return data;
}

/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type** &array)
{
  if (array == NULL) return;
  data = TYPE();
  sfree(array);
  array = NULL;
}
