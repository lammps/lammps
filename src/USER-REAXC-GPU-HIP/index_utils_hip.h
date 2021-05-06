#ifndef __INDEX_UTILS_H_
#define __INDEX_UTILS_H_

#include "reaxc_types_hip.h"


/* Indexing routine for grid cells */
static inline CUDA_HOST_DEVICE int index_grid_3d( int i, int j, int k,
        const grid * const g )
{
    return (i * g->ncells[1] * g->ncells[2]) + (j * g->ncells[2]) + k;
}

/* Indexing routine for grid cells, identical in functionality to above function */
static inline CUDA_HOST_DEVICE int index_grid_3d_v( ivec x,
        const grid * const g )
{
    return (x[0] * g->ncells[1] * g->ncells[2]) + (x[1] * g->ncells[2]) + x[2];
}

/* Indexing routine for neighbors of binned atoms within grid cells */
static inline CUDA_HOST_DEVICE int index_grid_nbrs( int i, int j, int k, int l,
        const grid * const g )
{
    return (i * g->ncells[1] * g->ncells[2] * g->max_nbrs) +
        (j * g->ncells[2] * g->max_nbrs) + (k * g->max_nbrs) + l;
}

/* Indexing routine for binned atoms within grid cells */
static inline CUDA_HOST_DEVICE int index_grid_atoms( int i, int j, int k, int l,
        const grid * const g )
{
    return (i * g->ncells[1] * g->ncells[2] * g->max_atoms) +
        (j * g->ncells[2] * g->max_atoms) + (k * g->max_atoms) + l;
}

/* Indexing routine for workspace system structures */
static inline CUDA_HOST_DEVICE int index_wkspace_sys( int i, int j, int N )
{
    return (i * N) + j;
}

/* Indexing routine for workspace res structures */
static inline CUDA_HOST_DEVICE int index_wkspace_res( int i, int j )
{
    return (i * (RESTART_REAX + 1)) + j;
}

/* Indexing routine for two body parameters */
static inline CUDA_HOST_DEVICE int index_tbp( int i, int j, int num_atom_types )
{
    return (i * num_atom_types) + j;
}

/* Indexing routine for three body parameters */
static inline CUDA_HOST_DEVICE int index_thbp( int i, int j, int k, int num_atom_types )
{
    return (i * num_atom_types * num_atom_types) + (j * num_atom_types) + k;
}

/* Indexing routine for hydrogen bonding parameters */
static inline CUDA_HOST_DEVICE int index_hbp( int i, int j, int k, int num_atom_types )
{
    return (i * num_atom_types * num_atom_types) + (j * num_atom_types) + k;
}

/* Indexing routine for four body parameters */
static inline CUDA_HOST_DEVICE int index_fbp( int i, int j, int k, int l, int num_atom_types )
{
    return (i * num_atom_types * num_atom_types * num_atom_types) +
        (j * num_atom_types * num_atom_types) + (k * num_atom_types) + l;
}

/* Indexing routine for LR table (force tabulation) */
static inline CUDA_HOST_DEVICE int index_lr( int i, int j, int num_atom_types )
{
	//printf("index %d \n", (i * num_atom_types) + j);
    return (i * num_atom_types) + j;
}

#endif
