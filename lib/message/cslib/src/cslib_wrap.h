/* ----------------------------------------------------------------------
   CSlib - Client/server library for code coupling
   http://cslib.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level CSlib directory.
------------------------------------------------------------------------- */

/* C style library interface to CSlib class
   ifdefs allow this file to be included in a C program
*/

#ifdef __cplusplus
extern "C" {
#endif

void cslib_open(int, const char *, const void *, const void *, void **);
void cslib_open_fortran(int, const char *, const char *, const void *, void **);
void cslib_open_fortran_mpi_one(int, const char *, const void *, 
                                const void *, void **);
void cslib_close(void *);

void cslib_send(void *, int, int);

void cslib_pack_int(void *, int, int);
void cslib_pack_int64(void *, int, int64_t);
void cslib_pack_float(void *, int, float);
void cslib_pack_double(void *, int, double);
void cslib_pack_string(void *, int, char *);
void cslib_pack(void *, int, int, int, void *);
void cslib_pack_parallel(void *, int, int, int, int *, int, void *);

int cslib_recv(void *, int *, int **, int **, int **);

int cslib_unpack_int(void *, int);
int64_t cslib_unpack_int64(void *, int);
float cslib_unpack_float(void *, int);
double cslib_unpack_double(void *, int);
char *cslib_unpack_string(void *, int);
void *cslib_unpack(void *, int);
void cslib_unpack_data(void *, int, void *);
void cslib_unpack_parallel(void *, int, int, int *, int, void *);

int cslib_extract(void *, int);
  
#ifdef __cplusplus
}
#endif
