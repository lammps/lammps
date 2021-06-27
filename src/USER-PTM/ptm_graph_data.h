/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef PTM_GRAPH_DATA_H
#define PTM_GRAPH_DATA_H

#include "ptm_constants.h"
#include <cstdint>

namespace ptm {

typedef struct {
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

}    // namespace ptm

#endif
