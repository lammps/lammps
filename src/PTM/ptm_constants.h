/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// clang-format off

#ifndef PTM_CONSTANTS_H
#define PTM_CONSTANTS_H

//------------------------------------
//    definitions
//------------------------------------
#define PTM_NO_ERROR         0


#define PTM_CHECK_FCC        (1 << 0)
#define PTM_CHECK_HCP        (1 << 1)
#define PTM_CHECK_BCC        (1 << 2)
#define PTM_CHECK_ICO        (1 << 3)
#define PTM_CHECK_SC         (1 << 4)
#define PTM_CHECK_DCUB       (1 << 5)
#define PTM_CHECK_DHEX       (1 << 6)
#define PTM_CHECK_GRAPHENE   (1 << 7)
#define PTM_CHECK_DEFAULT    (PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC)
#define PTM_CHECK_ALL        (PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC | PTM_CHECK_DCUB | PTM_CHECK_DHEX | PTM_CHECK_GRAPHENE)

#define PTM_MATCH_NONE       0
#define PTM_MATCH_FCC        1
#define PTM_MATCH_HCP        2
#define PTM_MATCH_BCC        3
#define PTM_MATCH_ICO        4
#define PTM_MATCH_SC         5
#define PTM_MATCH_DCUB       6
#define PTM_MATCH_DHEX       7
#define PTM_MATCH_GRAPHENE   8

#define PTM_ALLOY_NONE       0
#define PTM_ALLOY_PURE       1
#define PTM_ALLOY_L10        2
#define PTM_ALLOY_L12_CU     3
#define PTM_ALLOY_L12_AU     4
#define PTM_ALLOY_B2         5
#define PTM_ALLOY_SIC        6
#define PTM_ALLOY_BN         7


#define PTM_MAX_INPUT_POINTS 19
#define PTM_MAX_NBRS         16
#define PTM_MAX_POINTS       (PTM_MAX_NBRS + 1)
#define PTM_MAX_FACETS       28        //2 * PTM_MAX_NBRS - 4
#define PTM_MAX_EDGES        42        //3 * PTM_MAX_NBRS - 6


//------------------------------------
//    number of neighbours
//------------------------------------
#define PTM_NUM_NBRS_FCC     12
#define PTM_NUM_NBRS_HCP     12
#define PTM_NUM_NBRS_BCC     14
#define PTM_NUM_NBRS_ICO     12
#define PTM_NUM_NBRS_SC       6
#define PTM_NUM_NBRS_DCUB    16
#define PTM_NUM_NBRS_DHEX    16
#define PTM_NUM_NBRS_GRAPHENE 9

#define PTM_NUM_POINTS_FCC  (PTM_NUM_NBRS_FCC + 1)
#define PTM_NUM_POINTS_HCP  (PTM_NUM_NBRS_HCP + 1)
#define PTM_NUM_POINTS_BCC  (PTM_NUM_NBRS_BCC + 1)
#define PTM_NUM_POINTS_ICO  (PTM_NUM_NBRS_ICO + 1)
#define PTM_NUM_POINTS_SC   (PTM_NUM_NBRS_SC  + 1)
#define PTM_NUM_POINTS_DCUB (PTM_NUM_NBRS_DCUB + 1)
#define PTM_NUM_POINTS_DHEX (PTM_NUM_NBRS_DHEX + 1)
#define PTM_NUM_POINTS_GRAPHENE (PTM_NUM_NBRS_GRAPHENE + 1)

extern const int ptm_num_nbrs[9];

//------------------------------------
//    template structures
//------------------------------------

// these point sets have barycentre {0, 0, 0} and are scaled such that the mean neighbour distance is 1

extern const double ptm_template_fcc[PTM_NUM_POINTS_FCC][3];
extern const double ptm_template_hcp[PTM_NUM_POINTS_HCP][3];
extern const double ptm_template_bcc[PTM_NUM_POINTS_BCC][3];
extern const double ptm_template_ico[PTM_NUM_POINTS_ICO][3];
extern const double ptm_template_sc[PTM_NUM_POINTS_SC][3];
extern const double ptm_template_dcub[PTM_NUM_POINTS_DCUB][3];
extern const double ptm_template_dhex[PTM_NUM_POINTS_DHEX][3];
extern const double ptm_template_graphene[PTM_NUM_POINTS_GRAPHENE][3];

#endif

// clang-format on
