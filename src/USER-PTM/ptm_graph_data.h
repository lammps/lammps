#ifndef PTM_GRAPH_DATA_H
#define PTM_GRAPH_DATA_H

#include <stdint.h>
#include "ptm_constants.h"

namespace ptm {

typedef struct
{
	int id;
	uint64_t hash;
	int automorphism_index;
	int num_automorphisms;
	int8_t canonical_labelling[PTM_MAX_POINTS];
	int8_t facets[PTM_MAX_FACETS][3];
} graph_t;

#define NUM_SC_GRAPHS 1
#define NUM_ICO_GRAPHS 1
#define NUM_FCC_GRAPHS 8
#define NUM_HCP_GRAPHS 16
#define NUM_BCC_GRAPHS 218
#define NUM_DCUB_GRAPHS 12
#define NUM_DHEX_GRAPHS 24

extern int8_t automorphisms[][PTM_MAX_POINTS];

extern graph_t graphs_sc[NUM_SC_GRAPHS];
extern graph_t graphs_fcc[NUM_FCC_GRAPHS];
extern graph_t graphs_hcp[NUM_HCP_GRAPHS];
extern graph_t graphs_ico[NUM_ICO_GRAPHS];
extern graph_t graphs_bcc[NUM_BCC_GRAPHS];
extern graph_t graphs_dcub[NUM_DCUB_GRAPHS];
extern graph_t graphs_dhex[NUM_DHEX_GRAPHS];

}

#endif

