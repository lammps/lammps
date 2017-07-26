/* This file is part of ch5md
 * Copyright (C) 2013-2015 Pierre de Buyl
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD license. See the LICENSE file for details.
 */

#ifndef CH5MD_H
#define CH5MD_H

#include "hdf5.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CH5MD_RANK_ERROR -10

typedef struct h5md_element_struct {
  hid_t group;
  hid_t step;
  hid_t time;
  hid_t value;
  hid_t datatype;
  int is_time;
  int current_step;
  struct h5md_element_struct *link;
  struct h5md_particles_group_struct *particles_group;
} h5md_element;

typedef struct h5md_particles_group_struct {
  hid_t group;
  h5md_element position;
  hid_t box;
  h5md_element box_edges;
  h5md_element image;
  h5md_element velocity;
  h5md_element force;
  h5md_element mass;
  h5md_element species;
  h5md_element id;
  h5md_element charge;
  int local_size_max;
} h5md_particles_group;

typedef struct {
  hid_t id;
  int version[2];
  hid_t particles;
  hid_t observables;
  hid_t parameters;
} h5md_file;

h5md_file h5md_create_file (const char *filename, const char *author, const char *author_email, const char *creator, const char *creator_version);
int h5md_close_file(h5md_file file);
hid_t h5md_open_file (const char *filename);
h5md_particles_group h5md_create_particles_group(h5md_file file, const char *name);
hid_t h5md_open_particles_group(hid_t particles, const char *name);
h5md_element h5md_create_time_data(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, h5md_element *link);
int h5md_close_element(h5md_element e);
h5md_element h5md_create_fixed_data_simple(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, void *data);
h5md_element h5md_create_fixed_data_scalar(hid_t loc, const char *name, hid_t datatype, void *data);
int h5md_append(h5md_element e, void *data, int step, double time);
int h5md_create_box(h5md_particles_group *group, int dim, char *boundary[], bool is_time, double value[], h5md_element *link);
int h5md_write_string_attribute(hid_t loc, const char *obj_name,
    const char *att_name, const char *value);

#ifdef __cplusplus
}
#endif

#endif
