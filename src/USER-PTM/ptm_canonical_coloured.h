#ifndef PTM_CANONICAL_COLOURED_H
#define PTM_CANONICAL_COLOURED_H

#include <stdint.h>

namespace ptm {

int canonical_form_coloured(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* colours, int8_t* canonical_labelling, int8_t* best_code, uint64_t* p_hash);
}

#endif

