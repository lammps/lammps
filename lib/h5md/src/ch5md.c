/* This file is part of ch5md
 * Copyright (C) 2013-2015 Pierre de Buyl
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD license. See the LICENSE file for details.
 */

#include "hdf5.h"
#include "ch5md.h"
#include <string.h>
#include <stdlib.h>

#define MIN_CHUNK 10
#define MAX_CHUNK 256
#define MAX_RANK 5

h5md_file h5md_create_file (const char *filename, const char *author, const char *author_email, const char *creator, const char *creator_version)
{
  h5md_file file;
  hid_t g, g1;
  hid_t a, s, t;
  hsize_t dims[1];

  file.version[0] = 1;
  file.version[1] = 0;

  file.id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  g = H5Gcreate(file.id, "h5md", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  dims[0] = 2;
  s = H5Screate_simple(1, dims, NULL);
  a = H5Acreate(g, "version", H5T_NATIVE_INT, s, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(a, H5T_NATIVE_INT, file.version);
  H5Aclose(a);
  H5Sclose(s);

  s = H5Screate(H5S_SCALAR);

  g1 = H5Gcreate(g, "author", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  t = H5Tcopy(H5T_C_S1);
  H5Tset_size(t, strlen(author));
  a = H5Acreate(g1, "name", t, s, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(a, t, author);
  H5Aclose(a);
  H5Tclose(t);
  if (NULL!=author_email) {
    t = H5Tcopy(H5T_C_S1);
    H5Tset_size(t, strlen(author_email));
    a = H5Acreate(g1, "author_email", t, s, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, t, author_email);
    H5Aclose(a);
    H5Tclose(t);
  }
  H5Gclose(g1);

  g1 = H5Gcreate(g, "creator", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  t = H5Tcopy(H5T_C_S1);
  H5Tset_size(t, strlen(creator));
  a = H5Acreate(g1, "name", t, s, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(a, t, creator);
  H5Aclose(a);
  H5Tclose(t);
  t = H5Tcopy(H5T_C_S1);
  H5Tset_size(t, strlen(creator_version));
  a = H5Acreate(g1, "version", t, s, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(a, t, creator_version);
  H5Aclose(a);
  H5Tclose(t);
  H5Gclose(g1);

  file.particles = H5Gcreate(file.id, "particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  file.observables = H5Gcreate(file.id, "observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  file.parameters = H5Gcreate(file.id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  return file;
}

int h5md_close_file(h5md_file file) {
  H5Gclose(file.particles);
  H5Gclose(file.observables);
  H5Gclose(file.parameters);
  H5Fclose(file.id);

  return 0;
}

hid_t h5md_open_file (const char *filename)
{
  hid_t file;
  hid_t g;
  hid_t a, s;
  hsize_t dims[1];
  int version[2];
  int version_ok;

  file = H5Fopen(filename,  H5F_ACC_RDWR, H5P_DEFAULT);

  g = H5Gopen(file, "h5md", H5P_DEFAULT);

  version_ok = false;
  dims[0] = 2;
  a = H5Aopen(g, "version", H5P_DEFAULT);
  s = H5Aget_space(a);
  if (!(H5Sis_simple(s)>0)) {
    printf("H5MD version is not a simple dataspace");
    H5Sclose(s);
    H5Aclose(a);
    H5Gclose(g);
    H5Fclose(file);
  } else {
    if (H5Sget_simple_extent_ndims(s)==1) {
      H5Sget_simple_extent_dims(s, dims, NULL);
      if (dims[0]==2) {
	H5Aread(a, H5T_NATIVE_INT, version);
	if ( (version[0]==1) && (version[1]==0) ) {version_ok = true;}
      }
    }
  }

  H5Aclose(a);
  H5Sclose(s);
  H5Gclose(g);

  if (!version_ok) H5Fclose(file);

  return file;

}

h5md_particles_group h5md_create_particles_group(h5md_file file, const char *name)
{
  h5md_particles_group group;

  memset(&group,0,sizeof(group));
  group.group = H5Gcreate(file.particles, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  return group;
}

hid_t h5md_open_particles_group(hid_t particles, const char *name)
{
  hid_t group;

  group = H5Gopen(particles, name, H5P_DEFAULT);

  return group;
}

h5md_element h5md_create_time_data(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, h5md_element *link)
{

  h5md_element td;

  hid_t spc, plist;
  hsize_t dims[MAX_RANK], max_dims[MAX_RANK], chunks[MAX_RANK];

  int i;

  memset(&td,0,sizeof(td));
  dims[0] = 0 ;
  max_dims[0] = H5S_UNLIMITED ;
  for (i=0; i<rank; i++) {
    dims[i+1] = int_dims[i];
    max_dims[i+1] = int_dims[i];
  }
  chunks[0] = 1 ;
  if (MAX_CHUNK<int_dims[0]/4) {
    chunks[1] = MAX_CHUNK;
  } else {
    chunks[1] = int_dims[0];
  }
  for (i=1; i<rank; i++) {
    chunks[i+1]=int_dims[i];
  }

  td.group = H5Gcreate(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (NULL==link) {
    spc = H5Screate_simple( 1 , dims, max_dims);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 1, chunks);
    td.time = H5Dcreate(td.group, "time", H5T_NATIVE_DOUBLE, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
    td.step = H5Dcreate(td.group, "step", H5T_NATIVE_INT, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(spc);
    td.link=NULL;
  } else {
    td.link = link;
    H5Lcreate_hard(link->group, "step", td.group, "step", H5P_DEFAULT, H5P_DEFAULT);
    H5Lcreate_hard(link->group, "time", td.group, "time", H5P_DEFAULT, H5P_DEFAULT);
  }

  spc = H5Screate_simple( rank+1 , dims, max_dims) ;
  plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist, rank+1, chunks);
  td.value = H5Dcreate(td.group, "value", datatype, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Pclose(plist);
  H5Sclose(spc);

  td.datatype = datatype;
  td.is_time = true;

  return td;

}

int h5md_close_element(h5md_element e)
{
  if (!e.is_time) return 0;

  if (e.link==NULL) {
    H5Dclose(e.step);
    H5Dclose(e.time);
  }
  H5Dclose(e.value);
  H5Gclose(e.group);

  return 0;

}

h5md_element h5md_create_fixed_data_simple(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, void *data)
{

  h5md_element fd;

  hid_t spc;
  hsize_t dims[H5S_MAX_RANK];
  int i;

  memset(&fd,0,sizeof(fd));
  for (i=0; i<rank; i++) {
    dims[i] = int_dims[i];
  }

  spc = H5Screate_simple(rank, dims, NULL);
  fd.value = H5Dcreate(loc, name, datatype, spc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(spc);
  H5Dwrite(fd.value, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(fd.value);

  fd.is_time = false;

  return fd;

}

h5md_element h5md_create_fixed_data_scalar(hid_t loc, const char *name, hid_t datatype, void *data)
{

  h5md_element fd;
  hid_t spc;

  memset(&fd,0,sizeof(fd));
  spc = H5Screate(H5S_SCALAR);
  fd.value = H5Dcreate(loc, name, datatype, spc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(spc);
  H5Dwrite(fd.value, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(fd.value);

  fd.is_time = false;

  return fd;

}

int h5md_extend_by_one(hid_t dset, hsize_t *dims) {

  hid_t file_space;
  int rank;
  hsize_t maxdims[H5S_MAX_RANK];

  // Get dataset information
  file_space = H5Dget_space(dset);
  rank = H5Sget_simple_extent_ndims(file_space);
  if (rank > H5S_MAX_RANK) {
    return CH5MD_RANK_ERROR;
  }
  H5Sget_simple_extent_dims(file_space, dims, maxdims);
  H5Sclose(file_space);

  // Extend dimensions by one
  dims[0] = dims[0]+1;
  H5Dset_extent(dset, dims);

  return 0;

}

int h5md_append(h5md_element e, void *data, int step, double time) {

  hid_t mem_space, file_space;
  int i, rank;
  hsize_t dims[H5S_MAX_RANK];
  hsize_t start[H5S_MAX_RANK], count[H5S_MAX_RANK];

  // If not a time-dependent H5MD element, do nothing
  if (!e.is_time) return 0;

  if (NULL==e.link) {
    h5md_extend_by_one(e.step, dims);

    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;

    // Select the hyperslab
    file_space = H5Dget_space(e.step);
    rank = H5Sget_simple_extent_ndims(file_space);
    mem_space = H5Screate_simple(rank-1, dims+1, NULL);
    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;
    for (i=1 ; i<rank ; i++) {
      start[i] = 0;
      count[i] = dims[i];
    }
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(e.step, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, (void *)&step);
    H5Sclose(file_space);
    H5Sclose(mem_space);

    h5md_extend_by_one(e.time, dims);

    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;

    // Select the hyperslab
    file_space = H5Dget_space(e.time);
    rank = H5Sget_simple_extent_ndims(file_space);
    mem_space = H5Screate_simple(rank-1, dims+1, NULL);
    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;
    for (i=1 ; i<rank ; i++) {
      start[i] = 0;
      count[i] = dims[i];
    }
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(e.time, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, (void *)&time);
    H5Sclose(file_space);
    H5Sclose(mem_space);
  }

  h5md_extend_by_one(e.value, dims);

  // Select the hyperslab
  file_space = H5Dget_space(e.value);
  rank = H5Sget_simple_extent_ndims(file_space);
  mem_space = H5Screate_simple(rank-1, dims+1, NULL);

  // Define hyperslab selection
  start[0] = dims[0]-1;
  count[0] = 1;
  for (i=1 ; i<rank ; i++) {
    start[i] = 0;
    count[i] = dims[i];
  }
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
  H5Dwrite(e.value, e.datatype, mem_space, file_space, H5P_DEFAULT, data);
  H5Sclose(file_space);
  H5Sclose(mem_space);

  return 0;
}

int h5md_create_box(h5md_particles_group *group, int dim, char *boundary[], bool is_time, double value[], h5md_element *link)
{
  hid_t spc, att, t;
  hsize_t dims[1];
  int int_dims[1];
  int i;
  size_t boundary_length, tmp;
  //char *tmp_boundary;

  // Create box group
  group->box = H5Gcreate(group->group, "box", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create dimension attribute
  spc = H5Screate(H5S_SCALAR);
  att = H5Acreate(group->box, "dimension", H5T_NATIVE_INT, spc, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att, H5T_NATIVE_INT, &dim);
  H5Aclose(att);
  H5Sclose(spc);

  // Compute the size of the string type for boundary
  dims[0] = dim;
  boundary_length=0;
  for (i=0; i<dim; i++) {
    tmp = strlen(boundary[i])+1;
    if (tmp>boundary_length) {
      boundary_length=tmp;
    }
  }
  char *tmp_boundary = malloc(dim*sizeof(char)*boundary_length);
  for (i=0; i<dim; i++) {
    strcpy((tmp_boundary+i*boundary_length), boundary[i]);
  }
  // Create boundary attribute
  t = H5Tcopy(H5T_C_S1);
  H5Tset_size(t, boundary_length);
  spc = H5Screate_simple(1, dims, NULL);
  att = H5Acreate(group->box, "boundary", t, spc, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att, t, tmp_boundary);
  free(tmp_boundary);
  H5Aclose(att);
  H5Sclose(spc);
  H5Tclose(t);

  // Create edges
  // Check if the box is time-dependent or not
  int_dims[0]=dim;
  if (is_time) {
    group->box_edges = h5md_create_time_data(group->box, "edges", 1, int_dims, H5T_NATIVE_DOUBLE, link);
  } else {
    if (NULL!=value) {
      group->box_edges = h5md_create_fixed_data_simple(group->box, "edges", 1, int_dims, H5T_NATIVE_DOUBLE, value);
    }
  }

  H5Gclose(group->box);

  return 0;
}

int h5md_write_string_attribute(hid_t loc, const char *obj_name,
    const char *att_name, const char *value)
{
  hid_t obj;
  hid_t s, t, a;

  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);

  t = H5Tcopy(H5T_C_S1);
  H5Tset_size(t, strlen(value));
  s = H5Screate(H5S_SCALAR);
  a = H5Acreate(obj, att_name, t, s, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(a, t, value);
  H5Aclose(a);
  H5Sclose(s);
  H5Tclose(t);

  H5Oclose(obj);

  return 0;
}
