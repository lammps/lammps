/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

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

#include "reaxc_types.h"

#include "tool_box.h"



/* Safe wrapper around libc malloc
 *
 * n: num. of bytes to allocated
 * name: message with details about pointer, used for warnings/errors
 *
 * returns: ptr to allocated memory
 * */
void * smalloc( size_t n, const char *name )
{
    void *ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array %s.\n",
                n, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] requesting %zu bytes for %s\n", n, name );
    fflush( stderr );
#endif

    ptr = malloc( n );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array %s.\n",
                n, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] granted memory at address: %p\n", ptr );
    fflush( stderr );
#endif

    return ptr;
}


/* Safe wrapper around libc realloc
 *
 * n: num. of bytes to reallocated
 * name: message with details about pointer, used for warnings/errors
 *
 * returns: ptr to reallocated memory
 * */
void * srealloc( void *ptr, size_t n, const char *name )
{
    void *new_ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to reallocate %zu bytes for array %s.\n",
                n, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    if ( ptr == NULL )
    {
        fprintf( stderr, "[INFO] trying to allocate %zu NEW bytes for array %s.\n",
                n, name );
    }

    fprintf( stderr, "[INFO] requesting %zu bytes for %s\n", n, name );
    fflush( stderr );
#endif

    new_ptr = realloc( ptr, n );

    /* technically, ptr may still be allocated and valid,
     * but we needed more memory, so abort */
    if ( new_ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to reallocate %zu bytes for array %s.\n",
                n, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] granted memory at address: %p\n", new_ptr );
    fflush( stderr );
#endif

    return new_ptr;
}


/* Safe wrapper around libc calloc
 *
 * n: num. of elements to allocated (each of size bytes)
 * size: num. of bytes per element
 * name: message with details about pointer, used for warnings/errors
 *
 * returns: ptr to allocated memory, all bits initialized to zeros
 * */
void * scalloc( size_t n, size_t size, const char *name )
{
    void *ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array %s.\n",
                n * size, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] requesting %zu bytes for %s\n", n * size, name );
    fflush( stderr );
#endif

    ptr = calloc( n, size );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array %s.\n",
                n * size, name );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] granted memory at address: %p\n", ptr );
    fflush( stderr );
#endif

    return ptr;
}


/* Safe wrapper around libc free
 *
 * ptr: pointer to dynamically allocated memory which will be deallocated
 * name: message with details about pointer, used for warnings/errors
 * */
void sfree( void *ptr, const char *name )
{
    if ( ptr == NULL )
    {
        fprintf( stderr, "[WARNING] trying to free the already NULL pointer %s!\n",
                name );
        return;
    }

//    printf("Freeing pointer %s \n", name);

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] trying to free pointer (%s), address: %p\n",
            name, (void *) ptr );
    fflush( stderr );
#endif

    free( ptr );
}


/* Safe wrapper around libc fopen
 *
 * fname: name of file to be opened
 * mode: mode in which to open file
 * msg: message to be printed in case of error
 * */
FILE * sfopen( const char * fname, const char * mode, const char * msg )
{
    FILE * ptr;

    if ( fname == NULL )
    {
        fprintf( stderr, "[ERROR] trying to open file: NULL file name (%s). Terminating...\n",
                msg );
        exit( INVALID_INPUT );
    }
    if ( mode == NULL )
    {
        fprintf( stderr, "[ERROR] trying to open file: NULL mode (%s). Terminating...\n",
                msg );
        exit( INVALID_INPUT );
    }

    ptr = fopen( fname, mode );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to open file %s with mode %s (%s)\n",
              fname, mode, msg );
        exit( INVALID_INPUT );
    }

    return ptr;
}


/* Safe wrapper around libc fclose
 *
 * fname: name of file to be opened
 * mode: mode in which to open file
 * msg: message to be printed in case of error
 * */
void sfclose( FILE * fp, const char * msg )
{
    int ret;

    if ( fp == NULL )
    {
        fprintf( stderr, "[WARNING] trying to close NULL file pointer (%s). Returning...\n", msg );
        return;
    }

    ret = fclose( fp );

    if ( ret != 0 )
    {
        fprintf( stderr, "[ERROR] error detected when closing file (%s). Terminating...\n", msg );
        exit( INVALID_INPUT );
    }
}
