/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef PTM_CONVEX_HULL_INCREMENTAL_H
#define PTM_CONVEX_HULL_INCREMENTAL_H

#include "ptm_constants.h"
#include <cstdbool>
#include <cstdint>

namespace ptm {

typedef struct {
  int8_t facets[PTM_MAX_FACETS][3];
  double plane_normal[PTM_MAX_FACETS][3];
  bool processed[PTM_MAX_POINTS];
  int initial_vertices[4];
  double barycentre[3];
  int num_facets;
  int num_prev;
  bool ok;

} convexhull_t;

void add_facet(const double (*points)[3], int a, int b, int c, int8_t *facet, double *plane_normal,
               double *barycentre);
int get_convex_hull(int num_points, const double (*points)[3], convexhull_t *ch,
                    int8_t simplex[][3]);

}    // namespace ptm

#endif
