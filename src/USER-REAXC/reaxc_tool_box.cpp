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

#include "pair_reaxc.h"
#include "reaxc_tool_box.h"

struct timeval tim;
double t_end;

double Get_Time( )
{
  gettimeofday(&tim, NULL );
  return( tim.tv_sec + (tim.tv_usec / 1000000.0) );
}

int Tokenize( char* s, char*** tok )
{
  char test[MAX_LINE];
  const char *sep = (const char *)"\t \n\r\f!=";
  char *word;
  int count=0;

  strncpy( test, s, MAX_LINE-1);

  for( word = strtok(test, sep); word; word = strtok(NULL, sep) ) {
    strncpy( (*tok)[count], word, MAX_LINE );
    count++;
  }

  return count;
}

/* safe malloc */
void *smalloc( rc_bigint n, const char *name, MPI_Comm comm )
{
  void *ptr;

  if (n <= 0) {
    fprintf( stderr, "WARNING: trying to allocate %ld bytes for array %s. ",
             n, name );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  ptr = malloc( n );
  if (ptr == NULL) {
    fprintf( stderr, "ERROR: failed to allocate %ld bytes for array %s",
             n, name );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return ptr;
}


/* safe calloc */
void *scalloc( rc_bigint n, rc_bigint size, const char *name, MPI_Comm comm )
{
  void *ptr;

  if (n <= 0) {
    fprintf( stderr, "WARNING: trying to allocate %ld elements for array %s. ",
             n, name );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  if (size <= 0) {
    fprintf( stderr, "WARNING: elements size for array %s is %ld. ",
             name, size );
    fprintf( stderr, "returning NULL.\n" );
    return NULL;
  }

  ptr = calloc( n, size );
  if (ptr == NULL) {
    fprintf( stderr, "ERROR: failed to allocate %ld bytes for array %s",
             n*size, name );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return ptr;
}


/* safe free */
void sfree( void *ptr, const char *name )
{
  if (ptr == NULL) {
    fprintf( stderr, "WARNING: trying to free the already NULL pointer %s!\n",
             name );
    return;
  }

  free( ptr );
  ptr = NULL;
}

