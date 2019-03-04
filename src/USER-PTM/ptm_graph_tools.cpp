/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstring>
#include <algorithm>
#include "ptm_graph_tools.h"
#include "ptm_constants.h"


namespace ptm {

bool build_facet_map(int num_facets, int8_t facets[][3], int8_t common[PTM_MAX_NBRS][PTM_MAX_NBRS])
{
        memset(common, -1, sizeof(int8_t) * PTM_MAX_NBRS * PTM_MAX_NBRS);

        for (int i = 0;i<num_facets;i++)
        {
                int a = facets[i][0];
                int b = facets[i][1];
                int c = facets[i][2];

                //assert(common[a][b] == -1);
                //assert(common[b][c] == -1);
                //assert(common[c][a] == -1);
                if (common[a][b] != -1 || common[b][c] != -1 || common[c][a] != -1)
                        return false;

                common[a][b] = c;
                common[b][c] = a;
                common[c][a] = b;
        }

        return true;
}

int graph_degree(int num_facets, int8_t facets[][3], int num_nodes, int8_t* degree)
{
        memset(degree, 0, sizeof(int8_t) * num_nodes);

        for (int i = 0;i<num_facets;i++)
        {
                int a = facets[i][0];
                int b = facets[i][1];
                int c = facets[i][2];

                degree[a]++;
                degree[b]++;
                degree[c]++;
        }

        int8_t max_degree = 0;
        for (int i = 0;i<num_nodes;i++)
                max_degree = std::max(max_degree, degree[i]);

        return max_degree;
}

}

