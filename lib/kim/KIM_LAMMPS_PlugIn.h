// Copyright 2018 Regents of the University of Minnesota
// All rights reserved.
//
// Contributing authors: Ryan S. Elliott

#define KIM_STATUS_OK 1
#define KIM_STATUS_NEIGH_ITER_PAST_END 2
#define KIM_STATUS_NEIGH_ITER_INIT_OK 3
#define KIM_STATUS_NEIGH_INVALID_REQUEST -11
#define KIM_STATUS_NEIGH_INVALID_MODE -6
#define KIM_COMPUTE_FALSE 0
#define KIM_COMPUTE_TRUE 1
#define KIM_STATUS_FAIL 0

#include "dlfcn.h"
#include "error.h"

namespace LAMMPS_NS {

class KIM_LAMMPS_PlugIn {
 public:
  ~KIM_LAMMPS_PlugIn() {};
  KIM_LAMMPS_PlugIn() : library(NULL) {};

  bool setup_kim_api_library(Error * const error);
  bool is_strictly_between_1_5_and_2_0;

  // KIM symbols
  void* library;
  int (*report_error)(int line, const char *, const char *, int);
  void (*setm_compute_by_index)(void *, int *, ...);
  int (*model_compute)(void *);
  int (*model_init)(void *);
  int (*model_reinit)(void *);
  void * (*get_sim_buffer)(void *, int *);
  void * (*get_data_by_index)(void *, int, int *);
  int (*model_destroy)(void *);
  void (*free)(void *, int *);
  int (*string_init)(void *, const char *, const char *);
  int (*is_half_neighbors)(void *, int *);
  int (*get_NBC_method)(void *, const char **);
  int (*get_index)(void *, const char *, int *);
  int (*getm_index)(void *, int *, int, ...);
  int (*get_species_code)(void *, const char *, int *);
  void (*setm_data_by_index)(void *, int *, int, ...);
  int (*set_method_by_index)(void *, int, intptr_t, void (*)());
  void (*set_sim_buffer)(void *, void *, int *);
  int (*set_data_by_index)(void *, int, intptr_t, void *);
  void (*set_compute_by_index)(void *, int, int, int*);
  int (*get_num_params)(void *, int *, int *);
  int (*get_free_parameter)(void *, const int, const char **);
  void * (*get_data)(void *, const char *, int *);
  intptr_t (*get_rank)(void *, const char *, int *);
  intptr_t (*get_shape)(void *, const char *, int *, int *);
  int (*model_info)(void *, const char *);
};


bool KIM_LAMMPS_PlugIn::setup_kim_api_library(Error * error)
{
  if (library == NULL)
  {
#ifdef KIM_INSTALL_DIR
    library = dlopen(KIM_INSTALL_DIR "/lib/libkim-api-v1.so", RTLD_NOW);
#endif

    if (library == NULL)
      library = dlopen("kim-api-v1.so", RTLD_NOW);

    if (library == NULL)
    {
      error->all(FLERR,"KIM API library cannot be found");
      return false;
    }
    else
    {
      // check for version and set is_strictly_between_1_5_and_2_0
      void * ver = dlsym(library, "KIM_API_get_version");
      if (ver == NULL)
        is_strictly_between_1_5_and_2_0 = false;
      else
        is_strictly_between_1_5_and_2_0 = true;

      std::string error_prefix("KIM API library error: ");
      // get symbols
      if (NULL == (report_error =
                   (int (*)(int, const char *, const char *, int))
                   dlsym(library, "KIM_API_report_error")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (setm_compute_by_index =
                   (void (*)(void *, int *, ...))
                   dlsym(library, "KIM_API_setm_compute_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (model_compute = (int (*)(void *))
                   dlsym(library, "KIM_API_model_compute")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (model_init = (int (*)(void *))
                   dlsym(library, "KIM_API_model_init")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (model_reinit = (int (*)(void *))
                   dlsym(library, "KIM_API_model_reinit")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_sim_buffer = (void * (*)(void *, int *))
                   dlsym(library, "KIM_API_get_sim_buffer")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_data_by_index = (void * (*)(void *, int, int *))
                   dlsym(library, "KIM_API_get_data_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (model_destroy = (int (*)(void *))
                   dlsym(library, "KIM_API_model_destroy")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (free = (void (*)(void *, int *))
                   dlsym(library, "KIM_API_free")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (string_init =
                   (int (*)(void *, const char *, const char *))
                   dlsym(library, "KIM_API_string_init")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (is_half_neighbors = (int (*)(void *, int *))
                   dlsym(library, "KIM_API_is_half_neighbors")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_NBC_method = (int (*)(void *, const char **))
                   dlsym(library, "KIM_API_get_NBC_method")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_index = (int (*)(void *, const char *, int *))
                   dlsym(library, "KIM_API_get_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (getm_index = (int (*)(void *, int *, int, ...))
                   dlsym(library, "KIM_API_getm_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_species_code =
                   (int (*)(void *, const char *, int *))
                   dlsym(library, "KIM_API_get_species_code")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (setm_data_by_index =
                   (void (*)(void *, int *, int, ...))
                   dlsym(library, "KIM_API_setm_data_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (set_method_by_index =
                   (int (*)(void *, int, intptr_t, void (*)()))
                   dlsym(library, "KIM_API_set_method_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (set_sim_buffer = (void (*)(void *, void *, int *))
                   dlsym(library, "KIM_API_set_sim_buffer")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (set_data_by_index =
                   (int (*)(void *, int, intptr_t, void *))
                   dlsym(library, "KIM_API_set_data_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (set_compute_by_index =
                   (void (*)(void *, int, int, int*))
                   dlsym(library, "KIM_API_set_compute_by_index")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_num_params = (int (*)(void *, int *, int *))
                   dlsym(library, "KIM_API_get_num_params")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_free_parameter =
                   (int (*)(void *, const int, const char **))
                   dlsym(library, "KIM_API_get_free_parameter")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_data = (void * (*)(void *, const char *, int *))
                   dlsym(library, "KIM_API_get_data")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_rank =
                   (intptr_t (*)(void *, const char *, int *))
                   dlsym(library, "KIM_API_get_rank")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (get_shape =
                   (intptr_t (*)(void *, const char *, int *, int *))
                   dlsym(library, "KIM_API_get_shape")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

      if (NULL == (model_info = (int (*)(void *, const char *))
                   dlsym(library, "KIM_API_model_info")))
        error->all(FLERR,(error_prefix + std::string(dlerror())).c_str());

    return true;
    }
  }
  else
  {
    return true;
  }
}

}  // namespace LAMMPS_NS
