/* Copyright (C) 2013-2014 Pierre de Buyl
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
  a.  Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the <organization> nor the
     names of its contributors may be used to endorse or promote products
     derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CH5MD_H
#define CH5MD_H

#ifdef __cplusplus
extern "C" {
#endif
#include "hdf5.h"
#include <stdbool.h>

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
h5md_element h5md_create_time_data(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, h5md_element *link);
h5md_element h5md_create_fixed_data_simple(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, void *data);
h5md_element h5md_create_fixed_data_scalar(hid_t loc, const char *name, hid_t datatype, void *data);
int h5md_append(h5md_element e, void *data, int step, double time);
int h5md_create_box(h5md_particles_group *group, int dim, char *boundary[], bool is_time, double value[]);
#ifdef __cplusplus
}
#endif
#endif
