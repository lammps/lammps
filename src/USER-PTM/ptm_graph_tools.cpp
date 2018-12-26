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

