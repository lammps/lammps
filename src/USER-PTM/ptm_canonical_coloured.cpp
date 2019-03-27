/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstring>
#include <climits>
#include <algorithm>
#include "ptm_graph_tools.h"
#include "ptm_constants.h"

namespace ptm {

static bool weinberg_coloured(int num_nodes, int num_edges, int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS], int8_t* colours, int8_t* best_code, int8_t* canonical_labelling, int a, int b)
{
        bool m[PTM_MAX_NBRS][PTM_MAX_NBRS];
        memset(m, 0, sizeof(bool) * PTM_MAX_NBRS * PTM_MAX_NBRS);

        int8_t index[PTM_MAX_NBRS];
        memset(index, -1, sizeof(int8_t) * PTM_MAX_NBRS);


        int n = 0;
        index[a] = colours[a] * num_nodes + n++;
        if (index[a] > best_code[0])
                return false;

        bool winning = false;
        if (index[a] < best_code[0])
        {
                best_code[0] = index[a];
                winning = true;
        }

        int c = -1;
        for (int it=1;it<2*num_edges;it++)
        {
                bool newvertex = index[b] == -1;

                if (newvertex)
                        index[b] = colours[b] * num_nodes + n++;

                if (!winning && index[b] > best_code[it])
                        return false;

                if (winning || index[b] < best_code[it])
                {
                        winning = true;
                        best_code[it] = index[b];
                }

                if (newvertex)
                {
                        //When a new vertex is reached, take the right-most edge
                        //relative to the edge on which the vertex is reached.

                        c = common[a][b];
                }
                else if (m[b][a] == false)
                {
                        //When an old vertex is reached on a new path, go back
                        //in the opposite direction.

                        c = a;
                }
                else
                {
                        //When an old vertex is reached on an old path, leave the
                        //vertex on the right-most edge that has not previously
                        //been traversed in that direction.

                        c = common[a][b];
                        while (m[b][c] == true)
                                c = common[c][b];
                }

                m[a][b] = true;
                a = b;
                b = c;
        }

        if (winning)
        {
                memcpy(canonical_labelling, index, sizeof(int8_t) * num_nodes);
                return true;
        }

        return false;
}

int canonical_form_coloured(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree, int8_t* colours, int8_t* canonical_labelling, int8_t* best_code, uint64_t* p_hash)
{
        int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS] = {{0}};
        int num_edges = 3 * num_facets / 2;
        if (!build_facet_map(num_facets, facets, common))
                return -1;

        memset(best_code, SCHAR_MAX, sizeof(int8_t) * 2 * PTM_MAX_EDGES);

        bool equal = true;
        for (int i = 1;i<num_nodes;i++)
                if (degree[i] != degree[0] || colours[i] != colours[0])
                        equal = false;

        if (equal)
        {
                weinberg_coloured(num_nodes, num_edges, common, colours, best_code, canonical_labelling, facets[0][0], facets[0][1]);
        }
        else
        {
                uint32_t best_degree = 0;
                for (int i = 0;i<num_facets;i++)
                {
                        int a = facets[i][0];
                        int b = facets[i][1];
                        int c = facets[i][2];

                        //int da = colours[a] * num_nodes + degree[a];
                        //int db = colours[b] * num_nodes + degree[b];
                        //int dc = colours[c] * num_nodes + degree[c];

                        int da = degree[a];
                        int db = degree[b];
                        int dc = degree[c];

                        best_degree = std::max(best_degree, ((uint32_t)da << 16) | ((uint32_t)db << 8) | ((uint32_t)dc << 0));
                        best_degree = std::max(best_degree, ((uint32_t)da << 0) | ((uint32_t)db << 16) | ((uint32_t)dc << 8));
                        best_degree = std::max(best_degree, ((uint32_t)da << 8) | ((uint32_t)db << 0) | ((uint32_t)dc << 16));
                }

                for (int i = 0;i<num_facets;i++)
                {
                        int a = facets[i][0];
                        int b = facets[i][1];
                        int c = facets[i][2];

                        //int da = colours[a] * num_nodes + degree[a];
                        //int db = colours[b] * num_nodes + degree[b];
                        //int dc = colours[c] * num_nodes + degree[c];

                        int da = degree[a];
                        int db = degree[b];
                        int dc = degree[c];

                        if (best_degree == (((uint32_t)da << 16) | ((uint32_t)db << 8) | ((uint32_t)dc << 0)))
                                weinberg_coloured(num_nodes, num_edges, common, colours, best_code, canonical_labelling, a, b);

                        if (best_degree == (((uint32_t)da << 0) | ((uint32_t)db << 16) | ((uint32_t)dc << 8)))
                                weinberg_coloured(num_nodes, num_edges, common, colours, best_code, canonical_labelling, b, c);

                        if (best_degree == (((uint32_t)da << 8) | ((uint32_t)db << 0) | ((uint32_t)dc << 16)))
                                weinberg_coloured(num_nodes, num_edges, common, colours, best_code, canonical_labelling, c, a);
                }
        }

        for (int i = num_nodes-1;i>=0;i--)
                canonical_labelling[i+1] = (canonical_labelling[i] % num_nodes) + 1;
        canonical_labelling[0] = 0;

        uint64_t hash = 0;
        for (int i = 0;i<2 * num_edges;i++)
        {
                uint64_t e = best_code[i];
                e += i % 8;
                e &= 0xF;
                e <<= (4 * i) % 64;
                hash ^= e;
        }

        *p_hash = hash;
        return PTM_NO_ERROR;
}

}

