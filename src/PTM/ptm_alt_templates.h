/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// clang-format off

#ifndef PTM_ALT_TEMPLATES_H
#define PTM_ALT_TEMPLATES_H

#include "ptm_constants.h"

#include <cmath>

const double ptm_template_hcp_alt1[PTM_NUM_POINTS_HCP][3] = {
                                                                {          0,          0,          0 },
                                                                {          1,          0,          0 },
                                                                {       -0.5, -sqrt(3.0)/2,          0 },
                                                                {       -0.5, -sqrt(3.0)/6, -sqrt(6.0)/3 },
                                                                {          0,  sqrt(3.0)/3, -sqrt(6.0)/3 },
                                                                {        0.5, -sqrt(3.0)/6, -sqrt(6.0)/3 },
                                                                {         -1,          0,          0 },
                                                                {       -0.5,  sqrt(3.0)/2,          0 },
                                                                {        0.5,  sqrt(3.0)/2,          0 },
                                                                {        0.5, -sqrt(3.0)/2,          0 },
                                                                {        0.5, -sqrt(3.0)/6,  sqrt(6.0)/3 },
                                                                {          0,  sqrt(3.0)/3,  sqrt(6.0)/3 },
                                                                {       -0.5, -sqrt(3.0)/6,  sqrt(6.0)/3 },
};

const double ptm_template_dcub_alt1[PTM_NUM_POINTS_DCUB][3] = {
                                                                {                      0,                      0,                      0 },
                                                                {  4/(sqrt(3.0)+6*sqrt(2.0)), -4/(sqrt(3.0)+6*sqrt(2.0)),  4/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {  4/(sqrt(3.0)+6*sqrt(2.0)),  4/(sqrt(3.0)+6*sqrt(2.0)), -4/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                { -4/(sqrt(3.0)+6*sqrt(2.0)), -4/(sqrt(3.0)+6*sqrt(2.0)), -4/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                { -4/(sqrt(3.0)+6*sqrt(2.0)),  4/(sqrt(3.0)+6*sqrt(2.0)),  4/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {  8/(sqrt(3.0)+6*sqrt(2.0)),                      0,  8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {                      0, -8/(sqrt(3.0)+6*sqrt(2.0)),  8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {  8/(sqrt(3.0)+6*sqrt(2.0)), -8/(sqrt(3.0)+6*sqrt(2.0)),                      0 },
                                                                {                      0,  8/(sqrt(3.0)+6*sqrt(2.0)), -8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {  8/(sqrt(3.0)+6*sqrt(2.0)),                      0, -8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {  8/(sqrt(3.0)+6*sqrt(2.0)),  8/(sqrt(3.0)+6*sqrt(2.0)),                      0 },
                                                                { -8/(sqrt(3.0)+6*sqrt(2.0)),                      0, -8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {                      0, -8/(sqrt(3.0)+6*sqrt(2.0)), -8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                { -8/(sqrt(3.0)+6*sqrt(2.0)), -8/(sqrt(3.0)+6*sqrt(2.0)),                      0 },
                                                                { -8/(sqrt(3.0)+6*sqrt(2.0)),  8/(sqrt(3.0)+6*sqrt(2.0)),                      0 },
                                                                { -8/(sqrt(3.0)+6*sqrt(2.0)),                      0,  8/(sqrt(3.0)+6*sqrt(2.0)) },
                                                                {                      0,  8/(sqrt(3.0)+6*sqrt(2.0)),  8/(sqrt(3.0)+6*sqrt(2.0)) },
};

const double ptm_template_dhex_alt1[PTM_NUM_POINTS_DHEX][3] = {
        {                                   0,                                   0,                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),   -4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),   -4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),   -4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {                                   0,                                   0,       4*sqrt(3.0)/(sqrt(3.0)+6*sqrt(2.0)) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {      -8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {       8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
};

const double ptm_template_dhex_alt2[PTM_NUM_POINTS_DHEX][3] = {
        {                                   0,                                   0,                                   0 },
        {                                   0,   -8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {                                   0,                                   0,      -4*sqrt(3.0)/(sqrt(3.0)+6*sqrt(2.0)) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {                                   0,   -8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {      -8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),   4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {                                   0,   -8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
};

const double ptm_template_dhex_alt3[PTM_NUM_POINTS_DHEX][3] = {
        {                                   0,                                   0,                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),    4*sqrt(3.0)/(3*sqrt(3.0)+18*sqrt(2.0)) },
        {                                   0,                                   0,      -4*sqrt(3.0)/(sqrt(3.0)+6*sqrt(2.0)) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {      -8*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0,                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),      -4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)),  16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),       4*sqrt(6.0)/(sqrt(3.0)+6*sqrt(2.0)),                                   0 },
        {      -4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {                                   0,    8*sqrt(6.0)/(3*sqrt(3.0)+18*sqrt(2.0)), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
        {       4*sqrt(2.0)/(sqrt(3.0)+6*sqrt(2.0)),  -4*sqrt(6.0)/(3*(sqrt(3.0)+6*sqrt(2.0))), -16*sqrt(3.0)/(3*(sqrt(3.0)+6*sqrt(2.0))) },
};


const double ptm_template_graphene_alt1[PTM_NUM_POINTS_GRAPHENE][3] = {
                                                                        {                    0,                    0,                    0 },
                                                                        {   3*sqrt(3.0)/22-9./11,  -3./22+3*sqrt(3.0)/11,                    0 },
                                                                        {   9./11-3*sqrt(3.0)/22,  -3./22+3*sqrt(3.0)/11,                    0 },
                                                                        {                    0,  -6*sqrt(3.0)/11+3./11,                    0 },
                                                                        { -18./11+3*sqrt(3.0)/11,                    0,                    0 },
                                                                        {   3*sqrt(3.0)/22-9./11,  -9./22+9*sqrt(3.0)/11,                    0 },
                                                                        {   9./11-3*sqrt(3.0)/22,  -9./22+9*sqrt(3.0)/11,                    0 },
                                                                        { -3*sqrt(3.0)/11+18./11,                    0,                    0 },
                                                                        {   9./11-3*sqrt(3.0)/22,  -9*sqrt(3.0)/11+9./22,                    0 },
                                                                        {   3*sqrt(3.0)/22-9./11,  -9*sqrt(3.0)/11+9./22,                    0 },
};

#endif

// clang-format on
