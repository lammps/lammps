%{
/* extern void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr); */
extern void *lammps_open_no_mpi(int argc, char **argv, void **ptr);
extern void *lammps_open_fortran(int argc, char **argv, int f_comm);
extern void  lammps_close(void *handle);
extern void  lammps_mpi_init();
extern void  lammps_mpi_finalize();
extern void  lammps_free(void *ptr);
extern void  lammps_file(void *handle, const char *file);
extern char *lammps_command(void *handle, const char *cmd);
extern void  lammps_commands_list(void *handle, int ncmd, const char **cmds);
extern void  lammps_commands_string(void *handle, const char *str);
extern int    lammps_version(void *handle);
extern void   lammps_get_os_info(char *buffer, int buf_size);
extern void   lammps_memory_usage(void *handle, double *meminfo);
extern int    lammps_get_mpi_comm(void *handle);
extern double lammps_get_natoms(void *handle);
extern double lammps_get_thermo(void *handle, const char *keyword);
extern void   lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                          double *xy, double *yz, double *xz,
                          int *pflags, int *boxflag);
extern void   lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                        double xy, double yz, double xz);
extern int    lammps_extract_setting(void *handle, const char *keyword);
extern void  *lammps_extract_global(void *handle, const char *name);
extern void  *lammps_extract_atom(void *handle, const char *name);
extern int lammps_extract_global_datatype(void *handle, const char *name);
extern int lammps_extract_atom_datatype(void *handle, const char *name);
/* extern int    lammps_create_atoms(void *handle, int n, int *id, int *type,
 extern int    lammps_create_atoms(void *handle, int n, int64_t *id, int *type, */
extern void *lammps_extract_compute(void *handle, char *id, int, int);
extern void *lammps_extract_fix(void *handle, char *, int, int, int, int);
extern void *lammps_extract_variable(void *handle, char *, char *);
extern int   lammps_set_variable(void *, char *, char *);
extern void lammps_gather(void *, char *, int, int, void *);
extern void lammps_gather_concat(void *, char *, int, int, void *);
extern void lammps_gather_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_scatter(void *, char *, int, int, void *);
extern void lammps_scatter_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_gather_atoms(void *, char *, int, int, void *);
extern void lammps_gather_atoms_concat(void *, char *, int, int, void *);
extern void lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_scatter_atoms(void *, char *, int, int, void *);
extern void lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);
extern int lammps_config_has_mpi_support();
extern int lammps_config_has_package(const char *);
extern int lammps_config_package_count();
extern int lammps_config_package_name(int, char *, int);
extern int lammps_config_has_gzip_support();
extern int lammps_config_has_png_support();
extern int lammps_config_has_jpeg_support();
extern int lammps_config_has_ffmpeg_support();
extern int lammps_config_has_exceptions();
extern int lammps_has_style(void *, const char *, const char *);
extern int lammps_style_count(void *, const char *);
extern int lammps_style_name(void *, const char *, int, char *, int);
extern int lammps_has_id(void *, const char *, const char *);
extern int lammps_id_count(void *, const char *);
extern int lammps_id_name(void *, const char *, int, char *, int);
extern int lammps_find_pair_neighlist(void*, char *, int, int, int);
extern int lammps_find_fix_neighlist(void*, char *, int);
extern int lammps_find_compute_neighlist(void*, char *, int);
extern int lammps_neighlist_num_elements(void*, int);
extern void lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );
/*
extern int lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int image, int *flags);
extern int64_t lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int64_t image, int *flags);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_fix_external_set_energy_global(void *, char *, double);
extern void lammps_fix_external_set_virial_global(void *, char *, double *);
*/
extern int lammps_is_running(void *handle);
extern void lammps_force_timeout(void *handle);
extern int lammps_has_error(void *handle);
extern int lammps_get_last_error_message(void *handle, char *buffer, int buf_size);
%}

/* extern void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr); */
extern void *lammps_open_no_mpi(int argc, char **argv, void **ptr);
extern void *lammps_open_fortran(int argc, char **argv, int f_comm);
extern void  lammps_close(void *handle);
extern void  lammps_mpi_init();
extern void  lammps_mpi_finalize();
extern void  lammps_free(void *ptr);
extern void  lammps_file(void *handle, const char *file);
extern char *lammps_command(void *handle, const char *cmd);
extern void  lammps_commands_list(void *handle, int ncmd, const char **cmds);
extern void  lammps_commands_string(void *handle, const char *str);
extern int    lammps_version(void *handle);
extern void   lammps_get_os_info(char *buffer, int buf_size);
extern void   lammps_memory_usage(void *handle, double *meminfo);
extern int    lammps_get_mpi_comm(void *handle);
extern double lammps_get_natoms(void *handle);
extern double lammps_get_thermo(void *handle, const char *keyword);
extern void   lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                          double *xy, double *yz, double *xz,
                          int *pflags, int *boxflag);
extern void   lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                        double xy, double yz, double xz);
extern int    lammps_extract_setting(void *handle, const char *keyword);
extern void  *lammps_extract_global(void *handle, const char *name);
extern void  *lammps_extract_atom(void *handle, const char *name);
extern int lammps_extract_global_datatype(void *handle, const char *name);
extern int lammps_extract_atom_datatype(void *handle, const char *name);
/*
extern int    lammps_create_atoms(void *handle, int n, int *id, int *type,
extern int    lammps_create_atoms(void *handle, int n, int64_t *id, int *type, */
extern void *lammps_extract_compute(void *handle, char *id, int, int);
extern void *lammps_extract_fix(void *handle, char *, int, int, int, int);
extern void *lammps_extract_variable(void *handle, char *, char *);
extern int   lammps_set_variable(void *, char *, char *);
extern void lammps_gather(void *, char *, int, int, void *);
extern void lammps_gather_concat(void *, char *, int, int, void *);
extern void lammps_gather_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_scatter(void *, char *, int, int, void *);
extern void lammps_scatter_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_gather_atoms(void *, char *, int, int, void *);
extern void lammps_gather_atoms_concat(void *, char *, int, int, void *);
extern void lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void lammps_scatter_atoms(void *, char *, int, int, void *);
extern void lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);
extern int lammps_config_has_mpi_support();
extern int lammps_config_has_package(const char *);
extern int lammps_config_package_count();
extern int lammps_config_package_name(int, char *, int);
extern int lammps_config_has_gzip_support();
extern int lammps_config_has_png_support();
extern int lammps_config_has_jpeg_support();
extern int lammps_config_has_ffmpeg_support();
extern int lammps_config_has_exceptions();
extern int lammps_has_style(void *, const char *, const char *);
extern int lammps_style_count(void *, const char *);
extern int lammps_style_name(void *, const char *, int, char *, int);
extern int lammps_has_id(void *, const char *, const char *);
extern int lammps_id_count(void *, const char *);
extern int lammps_id_name(void *, const char *, int, char *, int);
extern int lammps_find_pair_neighlist(void*, char *, int, int, int);
extern int lammps_find_fix_neighlist(void*, char *, int);
extern int lammps_find_compute_neighlist(void*, char *, int);
extern int lammps_neighlist_num_elements(void*, int);
extern void lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );
/*
extern int lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int image, int *flags);
extern int64_t lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int64_t image, int *flags);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
extern void lammps_fix_external_set_energy_global(void *, char *, double);
extern void lammps_fix_external_set_virial_global(void *, char *, double *);
*/
extern int lammps_is_running(void *handle);
extern void lammps_force_timeout(void *handle);
extern int lammps_has_error(void *handle);
extern int lammps_get_last_error_message(void *handle, char *buffer, int buf_size);

