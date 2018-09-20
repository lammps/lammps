#ifndef PTM_INITIALIZE_DATA_H
#define PTM_INITIALIZE_DATA_H


#include "ptm_graph_data.h"
#include "ptm_graph_tools.h"
#include "ptm_deformation_gradient.h"
#include "ptm_fundamental_mappings.h"
#include "ptm_neighbour_ordering.h"
#include "ptm_canonical_coloured.h"
#include "ptm_convex_hull_incremental.h"


namespace ptm {

typedef struct
{
	int type;
	int num_nbrs;
	int num_facets;
	int max_degree;
	int num_graphs;
	int num_mappings;
	graph_t* graphs;
	const double (*points)[3];
	const double (*penrose)[3];
	const int8_t (*mapping)[PTM_MAX_POINTS];
} refdata_t;


//refdata_t structure_sc =  { .type = PTM_MATCH_SC,  .num_nbrs =  6, .num_facets =  8, .max_degree = 4, .num_graphs = NUM_SC_GRAPHS,  .graphs = graphs_sc,  .points = ptm_template_sc,  .penrose = penrose_sc , .mapping = mapping_sc };
const refdata_t structure_sc =   { PTM_MATCH_SC,    6,  8, 4, NUM_SC_GRAPHS,   NUM_CUBIC_MAPPINGS, graphs_sc,   ptm_template_sc,   penrose_sc,   mapping_sc   };
const refdata_t structure_fcc =  { PTM_MATCH_FCC,  12, 20, 6, NUM_FCC_GRAPHS,  NUM_CUBIC_MAPPINGS, graphs_fcc,  ptm_template_fcc,  penrose_fcc,  mapping_fcc  };
const refdata_t structure_hcp =  { PTM_MATCH_HCP,  12, 20, 6, NUM_HCP_GRAPHS,  NUM_HEX_MAPPINGS,   graphs_hcp,  ptm_template_hcp,  penrose_hcp,  mapping_hcp  };
const refdata_t structure_ico =  { PTM_MATCH_ICO,  12, 20, 6, NUM_ICO_GRAPHS,  NUM_ICO_MAPPINGS,   graphs_ico,  ptm_template_ico,  penrose_ico,  mapping_ico  };
const refdata_t structure_bcc =  { PTM_MATCH_BCC,  14, 24, 8, NUM_BCC_GRAPHS,  NUM_CUBIC_MAPPINGS, graphs_bcc,  ptm_template_bcc,  penrose_bcc,  mapping_bcc  };
const refdata_t structure_dcub = { PTM_MATCH_DCUB, 16, 28, 8, NUM_DCUB_GRAPHS, NUM_DCUB_MAPPINGS,  graphs_dcub, ptm_template_dcub, penrose_dcub, mapping_dcub };
const refdata_t structure_dhex = { PTM_MATCH_DHEX, 16, 28, 8, NUM_DHEX_GRAPHS, NUM_DHEX_MAPPINGS,  graphs_dhex, ptm_template_dhex, penrose_dhex, mapping_dhex };

}

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ptm_local_handle* ptm_local_handle_t;
ptm_local_handle_t ptm_initialize_local();
void ptm_uninitialize_local(ptm_local_handle_t ptr);

int ptm_initialize_global();

//------------------------------------
//    global initialization switch
//------------------------------------
extern bool ptm_initialized;


#ifdef __cplusplus
}
#endif


#endif

