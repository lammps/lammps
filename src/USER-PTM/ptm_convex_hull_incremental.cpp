/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cmath>
#include <cfloat>
#include <cstring>
#include <cassert>
#include "ptm_convex_hull_incremental.h"
#include "ptm_constants.h"

namespace ptm {

#define VISIBLE 1
#define INVISIBLE 2
#define BOTH 3
#define TOLERANCE 1E-8

static double norm_squared(double* p)
{
        double x = p[0];
        double y = p[1];
        double z = p[2];

        return x*x + y*y + z*z;
}

static double dot_product(const double* a, const double* b)
{
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross_product(double* a, double* b, double* c)
{
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
}

static void calculate_plane_normal(const double (*points)[3], int a, int b, int c, double* plane_normal)
{
        double u[3] = {        points[b][0] - points[a][0],
                        points[b][1] - points[a][1],
                        points[b][2] - points[a][2]        };

        double v[3] = {        points[c][0] - points[a][0],
                        points[c][1] - points[a][1],
                        points[c][2] - points[a][2]        };

        cross_product(u, v, plane_normal);
        double norm = sqrt(norm_squared(plane_normal));
        plane_normal[0] /= norm;
        plane_normal[1] /= norm;
        plane_normal[2] /= norm;
}

static double point_plane_distance(const double* w, const double* plane_point, const double* plane_cross)
{
        return          plane_cross[0] * (plane_point[0] - w[0])
                + plane_cross[1] * (plane_point[1] - w[1])
                + plane_cross[2] * (plane_point[2] - w[2]);
}

static bool calc_max_extent(int num_points, const double (*points)[3], int* min_index, int* max_index)
{
        for (int j=0;j<3;j++)
        {
                double dmin = DBL_MAX, dmax = -DBL_MAX;
                int imin = 0, imax = 0;

                for (int i = 0;i<num_points;i++)
                {
                        double d = points[i][j];
                        if (d < dmin)
                        {
                                dmin = d;
                                imin = i;
                        }
                        if (d > dmax)
                        {
                                dmax = d;
                                imax = i;
                        }
                }

                if (imin == imax)
                        return false;        //degenerate point set

                min_index[j] = imin;
                max_index[j] = imax;
        }

        return true;
}

static bool find_third_point(int num_points, const double (*points)[3], int a, int b, int* p_c)
{
        const double* x1 = points[a];
        const double* x2 = points[b];

        double x2x1[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        double ns_x2x1 = norm_squared(x2x1);

        int bi = -1;
        double max_dist = 0.0;
        for (int i = 0;i<num_points;i++)
        {
                if (i == a || i == b)
                        continue;

                const double* x0 = points[i];

                double x1x0[3] = {x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]};
                double dot = dot_product(x1x0, x2x1);
                double dist = (norm_squared(x1x0) * ns_x2x1 - dot*dot) / ns_x2x1;

                if (dist > max_dist)
                {
                        max_dist = dist;
                        bi = i;
                }
        }

        *p_c = bi;
        return max_dist > TOLERANCE;
}

static bool find_fourth_point(int num_points, const double (*points)[3], int a, int b, int c, int* p_d)
{
        double plane_normal[3];
        calculate_plane_normal(points, a, b, c, plane_normal);


        int bi = -1;
        double max_dist = 0.0;
        for (int i = 0;i<num_points;i++)
        {
                if (i == a || i == b || i == c)
                        continue;

                const double* x0 = points[i];
                double dist = fabs(point_plane_distance(x0, points[a], plane_normal));
                if (dist > max_dist)
                {
                        max_dist = dist;
                        bi = i;
                }
        }

        *p_d = bi;
        return max_dist > TOLERANCE;
}

static int initial_simplex(int num_points, const double (*points)[3], int* initial_vertices)
{
        int min_index[3] = {0};
        int max_index[3] = {0};
        if (!calc_max_extent(num_points, points, min_index, max_index))
                return -1;

        int bi = -1;
        double max_dist = 0.0;
        for (int i = 0;i<3;i++)
        {
                int a = min_index[i], b = max_index[i];
                double delta[3] = {        points[a][0] - points[b][0],
                                        points[a][1] - points[b][1],
                                        points[a][2] - points[b][2]        };
                double dist = norm_squared(delta);
                if (dist > max_dist)
                {
                        bi = i;
                        max_dist = dist;
                }
        }

        //first two points are (a, b)
        int a = min_index[bi], b = max_index[bi], c = -1, d = -1;

        if (!find_third_point(num_points, points, a, b, &c))
                return -2;

        if (!find_fourth_point(num_points, points, a, b, c, &d))
                return -3;

        initial_vertices[0] = a;
        initial_vertices[1] = b;
        initial_vertices[2] = c;
        initial_vertices[3] = d;
        return 0;
}

static bool visible(const double* w, const double* plane_point, const double* plane_normal)
{
        return point_plane_distance(w, plane_point, plane_normal) > 0;
}

void add_facet(const double (*points)[3], int a, int b, int c, int8_t* facet, double* plane_normal, double* barycentre)
{
        calculate_plane_normal(points, a, b, c, plane_normal);
        if (visible(barycentre, points[a], plane_normal))
        {
                plane_normal[0] = -plane_normal[0];
                plane_normal[1] = -plane_normal[1];
                plane_normal[2] = -plane_normal[2];

                facet[0] = b;
                facet[1] = a;
                facet[2] = c;
        }
        else
        {
                facet[0] = a;
                facet[1] = b;
                facet[2] = c;
        }
}

static int initialize_convex_hull(int num_points, const double (*points)[3], int8_t facets[][3], double plane_normal[][3], bool* processed, int* initial_vertices, double* barycentre)
{
        memset(processed, 0, PTM_MAX_POINTS * sizeof(bool));
        memset(barycentre, 0, 3 * sizeof(double));
        int ret = initial_simplex(num_points, points, initial_vertices);
        if (ret != 0)
                return ret;

        for (int i = 0;i<4;i++)
        {
                int a = initial_vertices[i];
                processed[a] = true;

                barycentre[0] += points[a][0];
                barycentre[1] += points[a][1];
                barycentre[2] += points[a][2];
        }
        barycentre[0] /= 4;
        barycentre[1] /= 4;
        barycentre[2] /= 4;

        add_facet(points, initial_vertices[0], initial_vertices[1], initial_vertices[2], facets[0], plane_normal[0], barycentre);
        add_facet(points, initial_vertices[0], initial_vertices[1], initial_vertices[3], facets[1], plane_normal[1], barycentre);
        add_facet(points, initial_vertices[0], initial_vertices[2], initial_vertices[3], facets[2], plane_normal[2], barycentre);
        add_facet(points, initial_vertices[1], initial_vertices[2], initial_vertices[3], facets[3], plane_normal[3], barycentre);
        return 0;
}

int get_convex_hull(int num_points, const double (*points)[3], convexhull_t* ch, int8_t simplex[][3])
{
        assert(        num_points == PTM_NUM_POINTS_FCC
                || num_points == PTM_NUM_POINTS_HCP
                || num_points == PTM_NUM_POINTS_BCC
                || num_points == PTM_NUM_POINTS_ICO
                || num_points == PTM_NUM_POINTS_SC
                || num_points == PTM_NUM_POINTS_DCUB
                || num_points == PTM_NUM_POINTS_DHEX);

        int ret = 0;
        int num_prev = ch->num_prev;
        ch->num_prev = num_points;
        if (!ch->ok || 0)
        {
                ret = initialize_convex_hull(num_points, points, ch->facets, ch->plane_normal, ch->processed, ch->initial_vertices, ch->barycentre);
                if (ret != 0)
                        return ret;

                ch->num_facets = 4;
                num_prev = 0;
        }

        for (int i = num_prev;i<num_points;i++)
        {
                if (ch->processed[i])
                        continue;
                ch->processed[i] = true;

                int num_to_add = 0;
                int8_t to_add[PTM_MAX_FACETS][3];
                int8_t edge_visible[PTM_MAX_POINTS][PTM_MAX_POINTS];
                memset(edge_visible, 0, sizeof(int8_t) * PTM_MAX_POINTS * PTM_MAX_POINTS);
                for (int j = 0;j<ch->num_facets;j++)
                {
                        int a = ch->facets[j][0];
                        int b = ch->facets[j][1];
                        int c = ch->facets[j][2];

                        int u = 0, v = 0, w = 0;

                        double distance = point_plane_distance(points[i], points[a], ch->plane_normal[j]);
                        bool vis = distance > TOLERANCE;
                        if (vis)
                        {
                                u = edge_visible[a][b] |= VISIBLE;
                                edge_visible[b][a] |= VISIBLE;

                                v = edge_visible[b][c] |= VISIBLE;
                                edge_visible[c][b] |= VISIBLE;

                                w = edge_visible[c][a] |= VISIBLE;
                                edge_visible[a][c] |= VISIBLE;

                                memcpy(ch->facets[j], ch->facets[ch->num_facets-1], 3 * sizeof(int8_t));
                                memcpy(ch->plane_normal[j], ch->plane_normal[ch->num_facets-1], 3 * sizeof(double));
                                ch->num_facets--;
                                j--;
                        }
                        else
                        {
                                u = edge_visible[a][b] |= INVISIBLE;
                                edge_visible[b][a] |= INVISIBLE;

                                v = edge_visible[b][c] |= INVISIBLE;
                                edge_visible[c][b] |= INVISIBLE;

                                w = edge_visible[c][a] |= INVISIBLE;
                                edge_visible[a][c] |= INVISIBLE;
                        }

                        if (u == BOTH)
                        {
                                to_add[num_to_add][0] = i;
                                to_add[num_to_add][1] = a;
                                to_add[num_to_add][2] = b;
                                num_to_add++;
                        }

                        if (v == BOTH)
                        {
                                to_add[num_to_add][0] = i;
                                to_add[num_to_add][1] = b;
                                to_add[num_to_add][2] = c;
                                num_to_add++;
                        }

                        if (w == BOTH)
                        {
                                to_add[num_to_add][0] = i;
                                to_add[num_to_add][1] = c;
                                to_add[num_to_add][2] = a;
                                num_to_add++;
                        }
                }

                for (int j = 0;j<num_to_add;j++)
                {
                        if (ch->num_facets >= PTM_MAX_FACETS)
                                return -4;

                        add_facet(points, to_add[j][0], to_add[j][1], to_add[j][2], ch->facets[ch->num_facets], ch->plane_normal[ch->num_facets], ch->barycentre); ch->num_facets++;
                }
        }

        for (int i=0;i<ch->num_facets;i++)
        {
                int a = ch->facets[i][0];
                int b = ch->facets[i][1];
                int c = ch->facets[i][2];
                if (a == 0 || b == 0 || c == 0)
                        return 1;                //central atom contained in convex hull

                simplex[i][0] = a - 1;
                simplex[i][1] = b - 1;
                simplex[i][2] = c - 1;
        }

        return ret;
}

}

