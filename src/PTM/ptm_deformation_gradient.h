/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// clang-format off

#ifndef PTM_DEFORMATION_GRADIENT_H
#define PTM_DEFORMATION_GRADIENT_H

#include "ptm_constants.h"

#include <cstdint>
#include <cmath>

namespace ptm {

void calculate_deformation_gradient(int num_points, const double (*ideal_points)[3], int8_t* mapping, double (*normalized)[3], const double (*penrose)[3], double* F, double* res);

const double penrose_sc[PTM_NUM_POINTS_SC][3] = {
        {    0,    0,    0 },
        {    0,    0, -0.5 },
        {    0,    0,  0.5 },
        {    0, -0.5,    0 },
        {    0,  0.5,    0 },
        { -0.5,    0,    0 },
        {  0.5,    0,    0 },
};

const double penrose_fcc[PTM_NUM_POINTS_FCC][3] = {
        {          0,          0,          0 },
        {  sqrt(2)/8,  sqrt(2)/8,          0 },
        {          0,  sqrt(2)/8,  sqrt(2)/8 },
        {  sqrt(2)/8,          0,  sqrt(2)/8 },
        { -sqrt(2)/8, -sqrt(2)/8,          0 },
        {          0, -sqrt(2)/8, -sqrt(2)/8 },
        { -sqrt(2)/8,          0, -sqrt(2)/8 },
        { -sqrt(2)/8,  sqrt(2)/8,          0 },
        {          0, -sqrt(2)/8,  sqrt(2)/8 },
        { -sqrt(2)/8,          0,  sqrt(2)/8 },
        {  sqrt(2)/8, -sqrt(2)/8,          0 },
        {          0,  sqrt(2)/8, -sqrt(2)/8 },
        {  sqrt(2)/8,          0, -sqrt(2)/8 },
};

const double penrose_hcp[PTM_NUM_POINTS_HCP][3] = {
        {           0,           0,           0 },
        {        1./8,  -sqrt(3)/8,           0 },
        {       -1./4,           0,           0 },
        {       -1./8,  sqrt(3)/24, -sqrt(6)/12 },
        {        1./8,  sqrt(3)/24, -sqrt(6)/12 },
        {           0, -sqrt(3)/12, -sqrt(6)/12 },
        {       -1./8,   sqrt(3)/8,           0 },
        {        1./8,   sqrt(3)/8,           0 },
        {        1./4,           0,           0 },
        {       -1./8,  -sqrt(3)/8,           0 },
        {           0, -sqrt(3)/12,  sqrt(6)/12 },
        {        1./8,  sqrt(3)/24,  sqrt(6)/12 },
        {       -1./8,  sqrt(3)/24,  sqrt(6)/12 },
};

const double penrose_hcp_alt1[PTM_NUM_POINTS_HCP][3] = {
        {           0,           0,           0 },
        {        1./4,           0,           0 },
        {       -1./8,  -sqrt(3)/8,           0 },
        {       -1./8, -sqrt(3)/24, -sqrt(6)/12 },
        {           0,  sqrt(3)/12, -sqrt(6)/12 },
        {        1./8, -sqrt(3)/24, -sqrt(6)/12 },
        {       -1./4,           0,           0 },
        {       -1./8,   sqrt(3)/8,           0 },
        {        1./8,   sqrt(3)/8,           0 },
        {        1./8,  -sqrt(3)/8,           0 },
        {        1./8, -sqrt(3)/24,  sqrt(6)/12 },
        {           0,  sqrt(3)/12,  sqrt(6)/12 },
        {       -1./8, -sqrt(3)/24,  sqrt(6)/12 },
};

const double penrose_ico[PTM_NUM_POINTS_ICO][3] = {
        {                          0,                          0,                          0 },
        {                          0,                          0,                       0.25 },
        {                          0,                          0,                      -0.25 },
        { -sqrt(-10*sqrt(5) + 50)/40,          sqrt(5)/40 + 1./8,                -sqrt(5)/20 },
        {  sqrt(-10*sqrt(5) + 50)/40,         -1./8 - sqrt(5)/40,                 sqrt(5)/20 },
        {                          0,                -sqrt(5)/10,                -sqrt(5)/20 },
        {                          0,                 sqrt(5)/10,                 sqrt(5)/20 },
        {   sqrt(10*sqrt(5) + 50)/40,         -1./8 + sqrt(5)/40,                -sqrt(5)/20 },
        {  -sqrt(10*sqrt(5) + 50)/40,         -sqrt(5)/40 + 1./8,                 sqrt(5)/20 },
        {  -sqrt(10*sqrt(5) + 50)/40,         -1./8 + sqrt(5)/40,                -sqrt(5)/20 },
        {   sqrt(10*sqrt(5) + 50)/40,         -sqrt(5)/40 + 1./8,                 sqrt(5)/20 },
        {  sqrt(-10*sqrt(5) + 50)/40,          sqrt(5)/40 + 1./8,                -sqrt(5)/20 },
        { -sqrt(-10*sqrt(5) + 50)/40,         -1./8 - sqrt(5)/40,                 sqrt(5)/20 },
};

const double penrose_bcc[PTM_NUM_POINTS_BCC][3] = {
        {                  0,                  0,                   0 },
        {  3./56 + sqrt(3)/28,  3./56 + sqrt(3)/28,  3./56 + sqrt(3)/28 },
        { -sqrt(3)/28 - 3./56,  3./56 + sqrt(3)/28,  3./56 + sqrt(3)/28 },
        {  3./56 + sqrt(3)/28,  3./56 + sqrt(3)/28, -sqrt(3)/28 - 3./56 },
        { -sqrt(3)/28 - 3./56, -sqrt(3)/28 - 3./56,  3./56 + sqrt(3)/28 },
        {  3./56 + sqrt(3)/28, -sqrt(3)/28 - 3./56,  3./56 + sqrt(3)/28 },
        { -sqrt(3)/28 - 3./56,  3./56 + sqrt(3)/28, -sqrt(3)/28 - 3./56 },
        { -sqrt(3)/28 - 3./56, -sqrt(3)/28 - 3./56, -sqrt(3)/28 - 3./56 },
        {  3./56 + sqrt(3)/28, -sqrt(3)/28 - 3./56, -sqrt(3)/28 - 3./56 },
        {  3./28 + sqrt(3)/14,                   0,                   0 },
        { -sqrt(3)/14 - 3./28,                   0,                   0 },
        {                   0,  3./28 + sqrt(3)/14,                   0 },
        {                   0, -sqrt(3)/14 - 3./28,                   0 },
        {                   0,                   0,  3./28 + sqrt(3)/14 },
        {                   0,                   0, -sqrt(3)/14 - 3./28 },
};

const double penrose_dcub[PTM_NUM_POINTS_DCUB][3] = {
        {                               0,                                0,                                0 },
        {  23./(48*(-sqrt(3) + 6*sqrt(2))),  23./(48*(-sqrt(3) + 6*sqrt(2))),  23./(48*(-sqrt(3) + 6*sqrt(2))) },
        {  23./(48*(-sqrt(3) + 6*sqrt(2))), -23./(-48*sqrt(3) + 288*sqrt(2)), -23./(-48*sqrt(3) + 288*sqrt(2)) },
        { -23./(-48*sqrt(3) + 288*sqrt(2)), -23./(-48*sqrt(3) + 288*sqrt(2)),  23./(48*(-sqrt(3) + 6*sqrt(2))) },
        { -23./(-48*sqrt(3) + 288*sqrt(2)),  23./(48*(-sqrt(3) + 6*sqrt(2))), -23./(-48*sqrt(3) + 288*sqrt(2)) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0 },
        {                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))),  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        {                                0, -23./(-24*sqrt(3) + 144*sqrt(2)), -23./(-24*sqrt(3) + 144*sqrt(2)) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))), -23./(-24*sqrt(3) + 144*sqrt(2)),                                0 },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),                               0, -23./(-24*sqrt(3) + 144*sqrt(2)) },
        { -23./(-24*sqrt(3) + 144*sqrt(2)), -23./(-24*sqrt(3) + 144*sqrt(2)),                                0 },
        {                                0, -23./(-24*sqrt(3) + 144*sqrt(2)),  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),                                0, -23./(-24*sqrt(3) + 144*sqrt(2)) },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0 },
        {                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))), -23./(-24*sqrt(3) + 144*sqrt(2)) },
};

const double penrose_dcub_alt1[PTM_NUM_POINTS_DCUB][3] = {
        {                                0,                                0,                                0 },
        {  23./(48*(-sqrt(3) + 6*sqrt(2))), -23./(-48*sqrt(3) + 288*sqrt(2)),  23./(48*(-sqrt(3) + 6*sqrt(2))) },
        {  23./(48*(-sqrt(3) + 6*sqrt(2))),  23./(48*(-sqrt(3) + 6*sqrt(2))), -23./(-48*sqrt(3) + 288*sqrt(2)) },
        { -23./(-48*sqrt(3) + 288*sqrt(2)), -23./(-48*sqrt(3) + 288*sqrt(2)), -23./(-48*sqrt(3) + 288*sqrt(2)) },
        { -23./(-48*sqrt(3) + 288*sqrt(2)),  23./(48*(-sqrt(3) + 6*sqrt(2))),  23./(48*(-sqrt(3) + 6*sqrt(2))) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        {                                0, -23./(-24*sqrt(3) + 144*sqrt(2)),  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))), -23./(-24*sqrt(3) + 144*sqrt(2)),                                0 },
        {                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))), -23./(-24*sqrt(3) + 144*sqrt(2)) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0, -23./(-24*sqrt(3) + 144*sqrt(2)) },
        {  23./(24*(-sqrt(3) + 6*sqrt(2))),  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0 },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),                                0, -23./(-24*sqrt(3) + 144*sqrt(2)) },
        {                                0, -23./(-24*sqrt(3) + 144*sqrt(2)), -23./(-24*sqrt(3) + 144*sqrt(2)) },
        { -23./(-24*sqrt(3) + 144*sqrt(2)), -23./(-24*sqrt(3) + 144*sqrt(2)),                                0 },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),  23./(24*(-sqrt(3) + 6*sqrt(2))),                                0 },
        { -23./(-24*sqrt(3) + 144*sqrt(2)),                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))) },
        {                                0,  23./(24*(-sqrt(3) + 6*sqrt(2))),  23./(24*(-sqrt(3) + 6*sqrt(2))) },
};


const double penrose_dhex[PTM_NUM_POINTS_DHEX][3] = {
        {                                        0,                                        0,                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {                                        0,  -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {                                        0,                                        0,   23*sqrt(3)/(48*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-24*sqrt(3) + 144*sqrt(2)),                                        0,                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {                                        0,  -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {   23*sqrt(2)/(24*(-sqrt(3) + 6*sqrt(2))),                                        0,                                        0 },
        {                                        0,  -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
};

const double penrose_dhex_alt1[PTM_NUM_POINTS_DHEX][3] = {
        {                                        0,                                        0,                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(3)/(-144*sqrt(3) + 864*sqrt(2)) },
        {                                        0,                                        0,   23*sqrt(3)/(48*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {  -23*sqrt(2)/(-24*sqrt(3) + 144*sqrt(2)),                                        0,                                        0 },
        {   23*sqrt(2)/(24*(-sqrt(3) + 6*sqrt(2))),                                        0,                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
};

const double penrose_dhex_alt2[PTM_NUM_POINTS_DHEX][3] = {
        {                                       0,                                       0,                                       0 },
        {                                       0, -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)), 23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        { -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), 23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        {  23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), 23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        {                                       0,                                       0, -23*sqrt(3)/(-48*sqrt(3) + 288*sqrt(2)) },
        { -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                       0 },
        {                                       0, -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)),  23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {  23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                       0 },
        { -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                       0 },
        { -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        { -23*sqrt(2)/(-24*sqrt(3) + 144*sqrt(2)),                                       0,                                       0 },
        {  23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {  23*sqrt(2)/(24*(-sqrt(3) + 6*sqrt(2))),                                       0,                                       0 },
        {  23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                       0 },
        { -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {  23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), 23*sqrt(6)/(144*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {                                       0, -23*sqrt(6)/(-72*sqrt(3) + 432*sqrt(2)), -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
};

const double penrose_dhex_alt3[PTM_NUM_POINTS_DHEX][3] = {
        {                                        0,                                        0,                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))),  23*sqrt(3)/(144*(-sqrt(3) + 6*sqrt(2))) },
        {                                        0,                                        0,  -23*sqrt(3)/(-48*sqrt(3) + 288*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {   23*sqrt(2)/(24*(-sqrt(3) + 6*sqrt(2))),                                        0,                                        0 },
        {  -23*sqrt(2)/(-24*sqrt(3) + 144*sqrt(2)),                                        0,                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),  -23*sqrt(6)/(-48*sqrt(3) + 288*sqrt(2)),                                        0 },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(3)/(36*(-sqrt(3) + 6*sqrt(2))) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)),   23*sqrt(6)/(48*(-sqrt(3) + 6*sqrt(2))),                                        0 },
        {  -23*sqrt(2)/(-48*sqrt(3) + 288*sqrt(2)), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {                                        0,   23*sqrt(6)/(72*(-sqrt(3) + 6*sqrt(2))),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
        {   23*sqrt(2)/(48*(-sqrt(3) + 6*sqrt(2))), -23*sqrt(6)/(-144*sqrt(3) + 864*sqrt(2)),  -23*sqrt(3)/(-36*sqrt(3) + 216*sqrt(2)) },
};


const double penrose_graphene[PTM_NUM_POINTS_GRAPHENE][3] = {
        {                     0,                     0,                   0 },
        {                     0,  2./63 + 4*sqrt(3)/63,                   0 },
        {    sqrt(3)/63 + 2./21, -2*sqrt(3)/63 - 1./63,                   0 },
        {   -2./21 - sqrt(3)/63, -2*sqrt(3)/63 - 1./63,                   0 },
        {   -2./21 - sqrt(3)/63,  1./21 + 2*sqrt(3)/21,                   0 },
        {    sqrt(3)/63 + 2./21,  1./21 + 2*sqrt(3)/21,                   0 },
        {  2*sqrt(3)/63 + 4./21,                     0,                   0 },
        {    sqrt(3)/63 + 2./21, -2*sqrt(3)/21 - 1./21,                   0 },
        {   -2./21 - sqrt(3)/63, -2*sqrt(3)/21 - 1./21,                   0 },
        { -4./21 - 2*sqrt(3)/63,                     0,                   0 },
};

const double penrose_graphene_alt1[PTM_NUM_POINTS_GRAPHENE][3] = {
        {                     0,                     0,                    0 },
        {   -2./21 - sqrt(3)/63,  1./63 + 2*sqrt(3)/63,                    0 },
        {    sqrt(3)/63 + 2./21,  1./63 + 2*sqrt(3)/63,                    0 },
        {                     0, -4*sqrt(3)/63 - 2./63,                    0 },
        { -4./21 - 2*sqrt(3)/63,                     0,                    0 },
        {   -2./21 - sqrt(3)/63,  1./21 + 2*sqrt(3)/21,                    0 },
        {    sqrt(3)/63 + 2./21,  1./21 + 2*sqrt(3)/21,                    0 },
        {  2*sqrt(3)/63 + 4./21,                     0,                    0 },
        {    sqrt(3)/63 + 2./21, -2*sqrt(3)/21 - 1./21,                    0 },
        {   -2./21 - sqrt(3)/63, -2*sqrt(3)/21 - 1./21,                    0 },
};

}

#endif

// clang-format on
