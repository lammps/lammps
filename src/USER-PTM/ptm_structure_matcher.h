/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef PTM_STRUCTURE_MATCHER_H
#define PTM_STRUCTURE_MATCHER_H

#include "ptm_initialize_data.h"
#include "ptm_constants.h"


namespace ptm {

typedef struct
{
        double rmsd;
        double scale;
        double q[4];                //rotation in quaternion form (rigid body transformation)
        int8_t mapping[PTM_MAX_POINTS];
        const refdata_t* ref_struct;
} result_t;

int match_general(const refdata_t* s, double (*ch_points)[3], double (*points)[3], convexhull_t* ch, result_t* res);
int match_fcc_hcp_ico(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);
int match_dcub_dhex(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);
int match_graphene(double (*points)[3], result_t* res);

}

#endif

