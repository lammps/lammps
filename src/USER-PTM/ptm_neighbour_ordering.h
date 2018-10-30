#ifndef PTM_NEIGHBOUR_ORDERING_H
#define PTM_NEIGHBOUR_ORDERING_H

#include <inttypes.h>

namespace ptm {

int calculate_neighbour_ordering(void* voronoi_handle, int num_points, const double (*_points)[3], int8_t* ordering);

int calculate_diamond_neighbour_ordering(	int num_points, double (*unpermuted_points)[3], int32_t* unpermuted_numbers,
						int8_t* ordering, double (*points)[3], int32_t* numbers);

void* voronoi_initialize_local();
void voronoi_uninitialize_local(void* ptr);

}

#endif

