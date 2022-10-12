// clang-format off
/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "ptm_structure_matcher.h"
#include <cstring>
#include <cmath>
#include <utility>
#include "ptm_convex_hull_incremental.h"
#include "ptm_canonical_coloured.h"
#include "ptm_graph_data.h"
#include "ptm_graph_tools.h"
#include "ptm_normalize_vertices.h"
#include "ptm_polar.h"
#include "ptm_constants.h"


namespace ptm {

static double calc_rmsd(int num_points, const double (*ideal_points)[3], double (*normalized)[3], int8_t* mapping,
                        double G1, double G2, double E0, double* q, double* p_scale)
{
        double A0[9];
        InnerProduct(A0, num_points, ideal_points, normalized, mapping);

        double nrmsdsq, rot[9];
        FastCalcRMSDAndRotation(A0, E0, &nrmsdsq, q, rot);

        double k0 = 0;
        for (int i=0;i<num_points;i++)
        {
                for (int j=0;j<3;j++)
                {
                        double v = 0.0;
                        for (int k=0;k<3;k++)
                                v += rot[j*3+k] * ideal_points[i][k];

                        k0 += v * normalized[mapping[i]][j];
                }
        }

        double scale = k0 / G2;
        *p_scale = scale;
        return sqrt(fabs(G1 - scale*k0) / num_points);
}

static void check_graphs(        const refdata_t* s,
                                uint64_t hash,
                                int8_t* canonical_labelling,
                                double (*normalized)[3],
                                result_t* res)
{
        int num_points = s->num_nbrs + 1;
        const double (*ideal_points)[3] = s->points;
        int8_t inverse_labelling[PTM_MAX_POINTS];
        int8_t mapping[PTM_MAX_POINTS];

        for (int i=0; i<num_points; i++)
                inverse_labelling[ canonical_labelling[i] ] = i;

        double G1 = 0, G2 = 0;
        for (int i=0;i<num_points;i++)
        {
                double x1 = ideal_points[i][0];
                double y1 = ideal_points[i][1];
                double z1 = ideal_points[i][2];

                double x2 = normalized[i][0];
                double y2 = normalized[i][1];
                double z2 = normalized[i][2];

                G1 += x1 * x1 + y1 * y1 + z1 * z1;
                G2 += x2 * x2 + y2 * y2 + z2 * z2;
        }
        double E0 = (G1 + G2) / 2;

        for (int i = 0;i<s->num_graphs;i++)
        {
                if (hash != s->graphs[i].hash)
                        continue;

                graph_t* gref = &s->graphs[i];
                for (int j = 0;j<gref->num_automorphisms;j++)
                {
                        for (int k=0;k<num_points;k++)
                                mapping[automorphisms[gref->automorphism_index + j][k]] = inverse_labelling[ gref->canonical_labelling[k] ];

                        double q[4], scale = 0;
                        double rmsd = calc_rmsd(num_points, ideal_points, normalized, mapping, G1, G2, E0, q, &scale);
                        if (rmsd < res->rmsd)
                        {
                                res->rmsd = rmsd;
                                res->scale = scale;
                                res->ref_struct = s;
                                memcpy(res->q, q, 4 * sizeof(double));
                                memcpy(res->mapping, mapping, sizeof(int8_t) * num_points);
                        }
                }
        }
}

int match_general(const refdata_t* s, double (*ch_points)[3], double (*points)[3], convexhull_t* ch, result_t* res)
{
        int8_t degree[PTM_MAX_NBRS];
        int8_t facets[PTM_MAX_FACETS][3];

        int ret = get_convex_hull(s->num_nbrs + 1, (const double (*)[3])ch_points, ch, facets);
        ch->ok = ret >= 0;
        if (ret != 0)
                return PTM_NO_ERROR;

        if (ch->num_facets != s->num_facets)
                return PTM_NO_ERROR;                        //incorrect number of facets in convex hull

        int max_degree = graph_degree(s->num_facets, facets, s->num_nbrs, degree);
        if (max_degree > s->max_degree)
                return PTM_NO_ERROR;

        if (s->type == PTM_MATCH_SC)
                for (int i = 0;i<s->num_nbrs;i++)
                        if (degree[i] != 4)
                                return PTM_NO_ERROR;

        double normalized[PTM_MAX_POINTS][3];
        subtract_barycentre(s->num_nbrs + 1, points, normalized);

        int8_t code[2 * PTM_MAX_EDGES];
        int8_t colours[PTM_MAX_POINTS] = {0};
        int8_t canonical_labelling[PTM_MAX_POINTS];
        uint64_t hash = 0;
        ret = canonical_form_coloured(s->num_facets, facets, s->num_nbrs, degree, colours, canonical_labelling, &code[0], &hash);
        if (ret != PTM_NO_ERROR)
                return ret;

        check_graphs(s, hash, canonical_labelling, normalized, res);
        return PTM_NO_ERROR;
}

int match_fcc_hcp_ico(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res)
{
        int num_nbrs = structure_fcc.num_nbrs;
        int num_facets = structure_fcc.num_facets;
        int max_degree = structure_fcc.max_degree;

        int8_t degree[PTM_MAX_NBRS];
        int8_t facets[PTM_MAX_FACETS][3];

        int ret = get_convex_hull(num_nbrs + 1, (const double (*)[3])ch_points, ch, facets);
        ch->ok = ret >= 0;
        if (ret != 0)
                return PTM_NO_ERROR;

        if (ch->num_facets != num_facets)
                return PTM_NO_ERROR;                        //incorrect number of facets in convex hull

        int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
        if (_max_degree > max_degree)
                return PTM_NO_ERROR;

        double normalized[PTM_MAX_POINTS][3];
        subtract_barycentre(num_nbrs + 1, points, normalized);

        int8_t code[2 * PTM_MAX_EDGES];
        int8_t colours[PTM_MAX_POINTS] = {0};
        int8_t canonical_labelling[PTM_MAX_POINTS];
        uint64_t hash = 0;
        ret = canonical_form_coloured(num_facets, facets, num_nbrs, degree, colours, canonical_labelling, &code[0], &hash);
        if (ret != PTM_NO_ERROR)
                return ret;

        if (flags & PTM_CHECK_FCC)        check_graphs(&structure_fcc, hash, canonical_labelling, normalized, res);
        if (flags & PTM_CHECK_HCP)        check_graphs(&structure_hcp, hash, canonical_labelling, normalized, res);
        if (flags & PTM_CHECK_ICO)        check_graphs(&structure_ico, hash, canonical_labelling, normalized, res);
        return PTM_NO_ERROR;
}

int match_dcub_dhex(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res)
{
        int num_nbrs = structure_dcub.num_nbrs;
        int num_facets = structure_fcc.num_facets;
        int max_degree = structure_dcub.max_degree;


        int8_t facets[PTM_MAX_FACETS][3];
        int ret = get_convex_hull(num_nbrs + 1, (const double (*)[3])ch_points, ch, facets);
        ch->ok = ret >= 0;
        if (ret != 0)
                return PTM_NO_ERROR;

        //check for facets with multiple inner atoms
        bool inverted[4] = {false, false, false, false};
        for (int i=0;i<ch->num_facets;i++)
        {
                int n = 0;
                for (int j=0;j<3;j++)
                {
                        if (facets[i][j] <= 3)
                        {
                                inverted[facets[i][j]] = true;
                                n++;
                        }
                }
                if (n > 1)
                        return PTM_NO_ERROR;
        }

        int num_inverted = 0;
        for (int i=0;i<4;i++)
                num_inverted += inverted[i] ? 1 : 0;

        if (ch->num_facets != num_facets + 2 * num_inverted)
                return PTM_NO_ERROR;                        //incorrect number of facets in convex hull

        int8_t degree[PTM_MAX_NBRS];
        int _max_degree = graph_degree(num_facets, facets, num_nbrs, degree);
        if (_max_degree > max_degree)
                return PTM_NO_ERROR;

        int num_found = 0;
        int8_t toadd[4][3];
        for (int i=0;i<ch->num_facets;i++)
        {
                int a = facets[i][0];
                int b = facets[i][1];
                int c = facets[i][2];
                if (a <= 3 || b <= 3 || c <= 3)
                        continue;

                int i0 = (a - 4) / 3;
                int i1 = (b - 4) / 3;
                int i2 = (c - 4) / 3;

                if (i0 == i1 && i0 == i2)
                {
                        if (num_found + num_inverted >= 4)
                                return PTM_NO_ERROR;

                        toadd[num_found][0] = a;
                        toadd[num_found][1] = b;
                        toadd[num_found][2] = c;
                        num_found++;

                        memcpy(&facets[i], &facets[ch->num_facets - 1], 3 * sizeof(int8_t));
                        ch->num_facets--;
                        i--;
                }
        }

        if (num_found + num_inverted != 4)
                return PTM_NO_ERROR;

        for (int i=0;i<num_found;i++)
        {
                int a = toadd[i][0];
                int b = toadd[i][1];
                int c = toadd[i][2];

                int i0 = (a - 4) / 3;

                facets[ch->num_facets][0] = i0;
                facets[ch->num_facets][1] = b;
                facets[ch->num_facets][2] = c;
                ch->num_facets++;

                facets[ch->num_facets][0] = a;
                facets[ch->num_facets][1] = i0;
                facets[ch->num_facets][2] = c;
                ch->num_facets++;

                facets[ch->num_facets][0] = a;
                facets[ch->num_facets][1] = b;
                facets[ch->num_facets][2] = i0;
                ch->num_facets++;
        }

        _max_degree = graph_degree(ch->num_facets, facets, num_nbrs, degree);
        if (_max_degree > max_degree)
                return PTM_NO_ERROR;

        double normalized[PTM_MAX_POINTS][3];
        subtract_barycentre(num_nbrs + 1, points, normalized);

        int8_t code[2 * PTM_MAX_EDGES];
        int8_t colours[PTM_MAX_POINTS] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        int8_t canonical_labelling[PTM_MAX_POINTS];
        uint64_t hash = 0;
        ret = canonical_form_coloured(ch->num_facets, facets, num_nbrs, degree, colours, canonical_labelling, &code[0], &hash);
        if (ret != PTM_NO_ERROR)
                return ret;

        if (flags & PTM_CHECK_DCUB)        check_graphs(&structure_dcub, hash, canonical_labelling, normalized, res);
        if (flags & PTM_CHECK_DHEX)        check_graphs(&structure_dhex, hash, canonical_labelling, normalized, res);

        return PTM_NO_ERROR;
}


static void check_graphs_graphene(        const refdata_t* s,
                                        int num_points,
                                        const double (*ideal_points)[3],
                                        double (*normalized)[3],
                                        int8_t* mapping,
                                        result_t* res)
{
        double G1 = 0, G2 = 0;
        for (int i=0;i<num_points;i++)
        {
                double x1 = ideal_points[i][0];
                double y1 = ideal_points[i][1];
                double z1 = ideal_points[i][2];

                double x2 = normalized[i][0];
                double y2 = normalized[i][1];
                double z2 = normalized[i][2];

                G1 += x1 * x1 + y1 * y1 + z1 * z1;
                G2 += x2 * x2 + y2 * y2 + z2 * z2;
        }
        double E0 = (G1 + G2) / 2;

        double q[4], scale = 0;
        double rmsd = calc_rmsd(num_points, ideal_points, normalized, mapping, G1, G2, E0, q, &scale);
        if (rmsd < res->rmsd)
        {
                res->rmsd = rmsd;
                res->scale = scale;
                res->ref_struct = s;
                memcpy(res->q, q, 4 * sizeof(double));
                memcpy(res->mapping, mapping, sizeof(int8_t) * num_points);
        }
}

int match_graphene(double (*points)[3], result_t* res)
{
        int num_nbrs = structure_graphene.num_nbrs;
        int num_points = num_nbrs + 1;
        const double (*ideal_points)[3] = structure_graphene.points;

        double normalized[PTM_MAX_POINTS][3];
        subtract_barycentre(num_points, points, normalized);

        int8_t mapping[PTM_MAX_POINTS];
        for (int i=0;i<num_points;i++)
                mapping[i] = i;

        for (int i=0;i<2;i++)
        {
                std::swap(mapping[4], mapping[5]);

                for (int j=0;j<2;j++)
                {
                        std::swap(mapping[6], mapping[7]);

                        for (int k=0;k<2;k++)
                        {
                                std::swap(mapping[8], mapping[9]);

                                check_graphs_graphene(&structure_graphene, num_points, ideal_points, normalized, mapping, res);
                        }
                }
        }

        return PTM_NO_ERROR;
}

}

