#ifndef PTM_GRAPH_TOOLS_H
#define PTM_GRAPH_TOOLS_H

#include <stdint.h>
#include "ptm_constants.h"

namespace ptm {

bool build_facet_map(int num_facets, int8_t facets[][3], int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS]);
int graph_degree(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree);

}

#endif

