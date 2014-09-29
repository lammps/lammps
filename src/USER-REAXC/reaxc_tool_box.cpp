/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reax_c.h"
#include "reaxc_tool_box.h"

void Transform( rvec x1, simulation_box *box, char flag, rvec x2 )
{
  int i, j;
  real tmp;

  if (flag > 0) {
    for (i=0; i < 3; i++) {
      tmp = 0.0;
      for (j=0; j < 3; j++)
        tmp += box->trans[i][j]*x1[j];
      x2[i] = tmp;
    }
  }
  else {
    for (i=0; i < 3; i++) {
      tmp = 0.0;
      for (j=0; j < 3; j++)
        tmp += box->trans_inv[i][j]*x1[j];
      x2[i] = tmp;
    }
  }
}


void Transform_to_UnitBox( rvec x1, simulation_box *box, char flag, rvec x2 )
{
  Transform( x1, box, flag, x2 );

  x2[0] /= box->box_norms[0];
  x2[1] /= box->box_norms[1];
  x2[2] /= box->box_norms[2];
}

void Fit_to_Periodic_Box( simulation_box *box, rvec *p )
{
  int i;

  for( i = 0; i < 3; ++i ) {
    if( (*p)[i] < box->min[i] ) {
      /* handle lower coords */
      while( (*p)[i] < box->min[i] )
        (*p)[i] += box->box_norms[i];
    }
    else if( (*p)[i] >= box->max[i] ) {
      /* handle higher coords */
      while( (*p)[i] >= box->max[i] )
        (*p)[i] -= box->box_norms[i];
    }
  }
}

struct timeval tim;
real t_end;

real Get_Time( )
{
  gettimeofday(&tim, NULL );
  return( tim.tv_sec + (tim.tv_usec / 1000000.0) );
}


real Get_Timing_Info( real t_start )
{
  gettimeofday(&tim, NULL );
  t_end = tim.tv_sec + (tim.tv_usec / 1000000.0);
  return (t_end - t_start);
}


void Update_Timing_Info( real *t_start, real *timing )
{
  gettimeofday(&tim, NULL );
  t_end = tim.tv_sec + (tim.tv_usec / 1000000.0);
  *timing += (t_end - *t_start);
  *t_start = t_end;
}

int Get_Atom_Type( reax_interaction *reax_param, char *s, MPI_Comm comm )
{
  int i;

  for( i = 0; i < reax_param->num_atom_types; ++i )
    if( !strcmp( reax_param->sbp[i].name, s ) )
      return i;

  fprintf( stderr, "Unknown atom type %s. Terminating...\n", s );
  MPI_Abort( comm, UNKNOWN_ATOM_TYPE );

  return -1;
}



char *Get_Element( reax_system *system, int i )
{
  return &( system->reax_param.sbp[system->my_atoms[i].type].name[0] );
}



char *Get_Atom_Name( reax_system *system, int i )
{
  return &(system->my_atoms[i].name[0]);
}



int Allocate_Tokenizer_Space( char **line, char **backup, char ***tokens )
{
  int i;

  if( (*line = (char*) malloc( sizeof(char) * MAX_LINE )) == NULL )
    return FAILURE;

  if( (*backup = (char*) malloc( sizeof(char) * MAX_LINE )) == NULL )
    return FAILURE;

  if( (*tokens = (char**) malloc( sizeof(char*) * MAX_TOKENS )) == NULL )
    return FAILURE;

  for( i = 0; i < MAX_TOKENS; i++ )
    if( ((*tokens)[i] = (char*) malloc(sizeof(char) * MAX_TOKEN_LEN)) == NULL )
      return FAILURE;

  return SUCCESS;
}



int Tokenize( char* s, char*** tok )
{
  char test[MAX_LINE];
  const char *sep = (const char *)"\t \n\r\f!=";
  char *word;
  int count=0;

  strncpy( test, s, MAX_LINE );

  for( word = strtok(test, sep); word; word = strtok(NULL, sep) ) {
    strncpy( (*tok)[count], word, MAX_LINE );
    count++;
  }

  return count;
}

/* safe malloc */
void *smalloc( long n, const char *name, MPI_Comm comm )
{
  void *ptr;

  if( n <= 0 ) {
    fprintf( stderr, "WARNING: trying to allocate %ld bytes for array %s. ",
             n, name );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  ptr = malloc( n );
  if( ptr == NULL ) {
    fprintf( stderr, "ERROR: failed to allocate %ld bytes for array %s",
             n, name );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return ptr;
}


/* safe calloc */
void *scalloc( int n, int size, const char *name, MPI_Comm comm )
{
  void *ptr;

  if( n <= 0 ) {
    fprintf( stderr, "WARNING: trying to allocate %d elements for array %s. ",
             n, name );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  if( size <= 0 ) {
    fprintf( stderr, "WARNING: elements size for array %s is %d. ",
             name, size );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  ptr = calloc( n, size );
  if( ptr == NULL ) {
    fprintf( stderr, "ERROR: failed to allocate %d bytes for array %s",
             n*size, name );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return ptr;
}


/* safe free */
void sfree( void *ptr, const char *name )
{
  if( ptr == NULL ) {
    fprintf( stderr, "WARNING: trying to free the already NULL pointer %s!\n",
             name );
    return;
  }

  free( ptr );
  ptr = NULL;
}

/* ----------------------------------------------------------------------
   strip off leading part of path, return just the filename
------------------------------------------------------------------------- */

static const char *potname(const char *path)
{
  const char *pot;

  if (path == NULL) return NULL;

#if defined(_WIN32)
  // skip over the disk drive part of windows pathnames
  if (isalpha(path[0]) && path[1] == ':') 
    path += 2;
#endif

  for (pot = path; *path != '\0'; ++path) {
#if defined(_WIN32)
    if ((*path == '\\') || (*path == '/')) pot = path + 1;
#else
    if (*path == '/') pot = path + 1;
#endif
  }

  return pot;
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name; failing that,
   search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */

FILE *lmp_open_potential(const char *name)
{
  FILE *fp;

  if (name == NULL) return NULL;

  // attempt to open file directly
  // if successful, return ptr

  fp = fopen(name,"r");
  if (fp) return fp;

  // try the environment variable directory

  const char *path = getenv("LAMMPS_POTENTIALS");
  if (path == NULL) return NULL;

  const char *pot = potname(name);
  if (pot == NULL) return NULL;

  size_t len1 = strlen(path);
  size_t len2 = strlen(pot);
  char *newpath = new char[len1+len2+2];

  strcpy(newpath,path);
#if defined(_WIN32)
  newpath[len1] = '\\';
  newpath[len1+1] = 0;
#else
  newpath[len1] = '/';
  newpath[len1+1] = 0;
#endif
  strcat(newpath,pot);

  fp = fopen(newpath,"r");
  delete[] newpath;
  return fp;
}


