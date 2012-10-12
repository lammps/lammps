/*
 * pair_manybody_const.h
 *
 *  Created on: Oct 11, 2011
 *      Author: chmu-tph
 */

#define MANYBODY_NPAIR 3

__device__ __constant__ int elem2param[(MANYBODY_NPAIR + 1) * (MANYBODY_NPAIR + 1) * (MANYBODY_NPAIR + 1)];
__device__ __constant__ int nelements;
__device__ __constant__ int map[MANYBODY_NPAIR + 2];
__device__ __constant__ int* _glob_numneigh_red;  //number of neighbors within force cutoff (as opposed to neighbor cutoff)
__device__ __constant__ int* _glob_neighbors_red; //indices of neighbors within force cutoff
__device__ __constant__ int* _glob_neightype_red; //type of neighbors within force cutoff

