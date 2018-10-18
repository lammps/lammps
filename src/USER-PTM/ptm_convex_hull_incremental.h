#ifndef PTM_CONVEX_HULL_INCREMENTAL_H
#define PTM_CONVEX_HULL_INCREMENTAL_H


#include <stdint.h>
#include <stdbool.h>
#include "ptm_constants.h"

namespace ptm {

typedef struct
{
	int8_t facets[PTM_MAX_FACETS][3];
	double plane_normal[PTM_MAX_FACETS][3];
	bool processed[PTM_MAX_POINTS];
	int initial_vertices[4];
	double barycentre[3];
	int num_facets;
	int num_prev;
	bool ok;

} convexhull_t;

void add_facet(const double (*points)[3], int a, int b, int c, int8_t* facet, double* plane_normal, double* barycentre);
int get_convex_hull(int num_points, const double (*points)[3], convexhull_t* ch, int8_t simplex[][3]);

}

#endif

