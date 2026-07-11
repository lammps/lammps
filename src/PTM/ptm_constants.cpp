// clang-format off
/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "ptm_constants.h"

const int ptm_num_nbrs[9] = {0, PTM_NUM_NBRS_FCC, PTM_NUM_NBRS_HCP, PTM_NUM_NBRS_BCC, PTM_NUM_NBRS_ICO, PTM_NUM_NBRS_SC, PTM_NUM_NBRS_DCUB, PTM_NUM_NBRS_DHEX, PTM_NUM_NBRS_GRAPHENE};

//------------------------------------
//    template structures
//------------------------------------

#define MYSQRT2 1.41421356237309504880
#define MYSQRT3 1.73205080756887729352
#define MYSQRT5 2.23606797749978969640
#define MYSQRT6 2.44948974278317809819
// sqrt((5.0-sqrt(5.0))/10.0):
#define MY5MINUSSQRT5BY10 0.52573111211913360602
// sqrt((5.0+sqrt(5.0))/10.0):
#define MY5PLUSSQRT5BY10 0.85065080835203993218
//these point sets have barycentre {0, 0, 0} and are scaled such that the mean neighbour distance is 1

const double ptm_template_fcc[PTM_NUM_POINTS_FCC][3] = {
                                                                {          0,          0,          0 },
                                                                {  MYSQRT2/2,  MYSQRT2/2,          0 },
                                                                {          0,  MYSQRT2/2,  MYSQRT2/2 },
                                                                {  MYSQRT2/2,          0,  MYSQRT2/2 },
                                                                { -MYSQRT2/2, -MYSQRT2/2,          0 },
                                                                {          0, -MYSQRT2/2, -MYSQRT2/2 },
                                                                { -MYSQRT2/2,          0, -MYSQRT2/2 },
                                                                { -MYSQRT2/2,  MYSQRT2/2,          0 },
                                                                {          0, -MYSQRT2/2,  MYSQRT2/2 },
                                                                { -MYSQRT2/2,          0,  MYSQRT2/2 },
                                                                {  MYSQRT2/2, -MYSQRT2/2,          0 },
                                                                {          0,  MYSQRT2/2, -MYSQRT2/2 },
                                                                {  MYSQRT2/2,          0, -MYSQRT2/2 },
};

const double ptm_template_hcp[PTM_NUM_POINTS_HCP][3] = {
                                                                {          0,          0,          0 },
                                                                {        0.5, -MYSQRT3/2,          0 },
                                                                {         -1,          0,          0 },
                                                                {       -0.5,  MYSQRT3/6, -MYSQRT6/3 },
                                                                {        0.5,  MYSQRT3/6, -MYSQRT6/3 },
                                                                {          0, -MYSQRT3/3, -MYSQRT6/3 },
                                                                {       -0.5,  MYSQRT3/2,          0 },
                                                                {        0.5,  MYSQRT3/2,          0 },
                                                                {          1,          0,          0 },
                                                                {       -0.5, -MYSQRT3/2,          0 },
                                                                {          0, -MYSQRT3/3,  MYSQRT6/3 },
                                                                {        0.5,  MYSQRT3/6,  MYSQRT6/3 },
                                                                {       -0.5,  MYSQRT3/6,  MYSQRT6/3 },
};

const double ptm_template_bcc[PTM_NUM_POINTS_BCC][3] = {
                                                                {                0,                0,                0 },
                                                                { 7*MYSQRT3/3-7./2, 7*MYSQRT3/3-7./2, 7*MYSQRT3/3-7./2 },
                                                                { 7./2-7*MYSQRT3/3, 7*MYSQRT3/3-7./2, 7*MYSQRT3/3-7./2 },
                                                                { 7*MYSQRT3/3-7./2, 7*MYSQRT3/3-7./2, 7./2-7*MYSQRT3/3 },
                                                                { 7./2-7*MYSQRT3/3, 7./2-7*MYSQRT3/3, 7*MYSQRT3/3-7./2 },
                                                                { 7*MYSQRT3/3-7./2, 7./2-7*MYSQRT3/3, 7*MYSQRT3/3-7./2 },
                                                                { 7./2-7*MYSQRT3/3, 7*MYSQRT3/3-7./2, 7./2-7*MYSQRT3/3 },
                                                                { 7./2-7*MYSQRT3/3, 7./2-7*MYSQRT3/3, 7./2-7*MYSQRT3/3 },
                                                                { 7*MYSQRT3/3-7./2, 7./2-7*MYSQRT3/3, 7./2-7*MYSQRT3/3 },
                                                                {   14*MYSQRT3/3-7,                0,                0 },
                                                                {   7-14*MYSQRT3/3,                0,                0 },
                                                                {                0,   14*MYSQRT3/3-7,                0 },
                                                                {                0,   7-14*MYSQRT3/3,                0 },
                                                                {                0,                0,   14*MYSQRT3/3-7 },
                                                                {                0,                0,   7-14*MYSQRT3/3 },
};

const double ptm_template_ico[PTM_NUM_POINTS_ICO][3] = {
                                                                {                     0,                     0,                     0 },
                                                                {                     0,                     0,                     1 },
                                                                {                     0,                     0,                    -1 },
                                                                {    -MY5MINUSSQRT5BY10,        (5+MYSQRT5)/10,            -MYSQRT5/5 },
                                                                {     MY5MINUSSQRT5BY10,       -(5+MYSQRT5)/10,             MYSQRT5/5 },
                                                                {                     0,          -2*MYSQRT5/5,            -MYSQRT5/5 },
                                                                {                     0,           2*MYSQRT5/5,             MYSQRT5/5 },
                                                                {      MY5PLUSSQRT5BY10,       -(5-MYSQRT5)/10,            -MYSQRT5/5 },
                                                                {     -MY5PLUSSQRT5BY10,        (5-MYSQRT5)/10,             MYSQRT5/5 },
                                                                {     -MY5PLUSSQRT5BY10,       -(5-MYSQRT5)/10,            -MYSQRT5/5 },
                                                                {      MY5PLUSSQRT5BY10,        (5-MYSQRT5)/10,             MYSQRT5/5 },
                                                                {     MY5MINUSSQRT5BY10,        (5+MYSQRT5)/10,            -MYSQRT5/5 },
                                                                {    -MY5MINUSSQRT5BY10,       -(5+MYSQRT5)/10,             MYSQRT5/5 },
};

const double ptm_template_sc[PTM_NUM_POINTS_SC][3] = {
                                                                {  0,  0,  0 },
                                                                {  0,  0, -1 },
                                                                {  0,  0,  1 },
                                                                {  0, -1,  0 },
                                                                {  0,  1,  0 },
                                                                { -1,  0,  0 },
                                                                {  1,  0,  0 },
};

const double ptm_template_dcub[PTM_NUM_POINTS_DCUB][3] = {
                                                                {                      0,                      0,                      0 },
                                                                {  4/(MYSQRT3+6*MYSQRT2),  4/(MYSQRT3+6*MYSQRT2),  4/(MYSQRT3+6*MYSQRT2) },
                                                                {  4/(MYSQRT3+6*MYSQRT2), -4/(MYSQRT3+6*MYSQRT2), -4/(MYSQRT3+6*MYSQRT2) },
                                                                { -4/(MYSQRT3+6*MYSQRT2), -4/(MYSQRT3+6*MYSQRT2),  4/(MYSQRT3+6*MYSQRT2) },
                                                                { -4/(MYSQRT3+6*MYSQRT2),  4/(MYSQRT3+6*MYSQRT2), -4/(MYSQRT3+6*MYSQRT2) },
                                                                {  8/(MYSQRT3+6*MYSQRT2),  8/(MYSQRT3+6*MYSQRT2),                      0 },
                                                                {                      0,  8/(MYSQRT3+6*MYSQRT2),  8/(MYSQRT3+6*MYSQRT2) },
                                                                {  8/(MYSQRT3+6*MYSQRT2),                      0,  8/(MYSQRT3+6*MYSQRT2) },
                                                                {                      0, -8/(MYSQRT3+6*MYSQRT2), -8/(MYSQRT3+6*MYSQRT2) },
                                                                {  8/(MYSQRT3+6*MYSQRT2), -8/(MYSQRT3+6*MYSQRT2),                      0 },
                                                                {  8/(MYSQRT3+6*MYSQRT2),                      0, -8/(MYSQRT3+6*MYSQRT2) },
                                                                { -8/(MYSQRT3+6*MYSQRT2), -8/(MYSQRT3+6*MYSQRT2),                      0 },
                                                                {                      0, -8/(MYSQRT3+6*MYSQRT2),  8/(MYSQRT3+6*MYSQRT2) },
                                                                { -8/(MYSQRT3+6*MYSQRT2),                      0,  8/(MYSQRT3+6*MYSQRT2) },
                                                                { -8/(MYSQRT3+6*MYSQRT2),                      0, -8/(MYSQRT3+6*MYSQRT2) },
                                                                { -8/(MYSQRT3+6*MYSQRT2),  8/(MYSQRT3+6*MYSQRT2),                      0 },
                                                                {                      0,  8/(MYSQRT3+6*MYSQRT2), -8/(MYSQRT3+6*MYSQRT2) },
};

const double ptm_template_dhex[PTM_NUM_POINTS_DHEX][3] = {
                                                                {                                   0,                                   0,                                   0 },
                                                                {      -4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)),   -4*MYSQRT3/(3*MYSQRT3+18*MYSQRT2) },
                                                                {                                   0,   -8*MYSQRT6/(3*MYSQRT3+18*MYSQRT2),   -4*MYSQRT3/(3*MYSQRT3+18*MYSQRT2) },
                                                                {       4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)),   -4*MYSQRT3/(3*MYSQRT3+18*MYSQRT2) },
                                                                {                                   0,                                   0,       4*MYSQRT3/(MYSQRT3+6*MYSQRT2) },
                                                                {      -8*MYSQRT2/(MYSQRT3+6*MYSQRT2),                                   0,                                   0 },
                                                                {      -4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)), -16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
                                                                {      -4*MYSQRT2/(MYSQRT3+6*MYSQRT2),       4*MYSQRT6/(MYSQRT3+6*MYSQRT2),                                   0 },
                                                                {       4*MYSQRT2/(MYSQRT3+6*MYSQRT2),      -4*MYSQRT6/(MYSQRT3+6*MYSQRT2),                                   0 },
                                                                {                                   0,   -8*MYSQRT6/(3*MYSQRT3+18*MYSQRT2), -16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
                                                                {      -4*MYSQRT2/(MYSQRT3+6*MYSQRT2),      -4*MYSQRT6/(MYSQRT3+6*MYSQRT2),                                   0 },
                                                                {       4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)), -16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
                                                                {       4*MYSQRT2/(MYSQRT3+6*MYSQRT2),       4*MYSQRT6/(MYSQRT3+6*MYSQRT2),                                   0 },
                                                                {       8*MYSQRT2/(MYSQRT3+6*MYSQRT2),                                   0,                                   0 },
                                                                {                                   0,   -8*MYSQRT6/(3*MYSQRT3+18*MYSQRT2),  16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
                                                                {       4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)),  16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
                                                                {      -4*MYSQRT2/(MYSQRT3+6*MYSQRT2),   4*MYSQRT6/(3*(MYSQRT3+6*MYSQRT2)),  16*MYSQRT3/(3*(MYSQRT3+6*MYSQRT2)) },
};

const double ptm_template_graphene[PTM_NUM_POINTS_GRAPHENE][3] = {
                                                                {                    0,                    0,                    0 },
                                                                {                    0,  -3./11+6*MYSQRT3/11,                    0 },
                                                                {  -3*MYSQRT3/22+9./11,  -3*MYSQRT3/11+3./22,                    0 },
                                                                {  -9./11+3*MYSQRT3/22,  -3*MYSQRT3/11+3./22,                    0 },
                                                                {  -9./11+3*MYSQRT3/22,  -9./22+9*MYSQRT3/11,                    0 },
                                                                {  -3*MYSQRT3/22+9./11,  -9./22+9*MYSQRT3/11,                    0 },
                                                                { -3*MYSQRT3/11+18./11,                    0,                    0 },
                                                                {  -3*MYSQRT3/22+9./11,  -9*MYSQRT3/11+9./22,                    0 },
                                                                {  -9./11+3*MYSQRT3/22,  -9*MYSQRT3/11+9./22,                    0 },
                                                                { -18./11+3*MYSQRT3/11,                    0,                    0 },
};
