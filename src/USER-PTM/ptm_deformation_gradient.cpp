/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "ptm_deformation_gradient.h"

namespace ptm {

void calculate_deformation_gradient(int num_points, const double (*ideal_points)[3], int8_t* mapping, double (*normalized)[3], const double (*penrose)[3], double* F, double* res)
{
        for (int i = 0;i<3;i++)
        {
                for (int j = 0;j<3;j++)
                {
                        double acc = 0.0;
                        for (int k = 0;k<num_points;k++)
                                acc += penrose[k][j] * normalized[mapping[k]][i];

                        F[i*3 + j] = acc;
                }
        }

        res[0] = 0;
        res[1] = 0;
        res[2] = 0;

        for (int k = 0;k<num_points;k++)
        {
                for (int i = 0;i<3;i++)
                {
                        double acc = 0.0;
                        for (int j = 0;j<3;j++)
                        {
                                acc += F[i*3 + j] * ideal_points[k][j];
                        }

                        double delta = acc - normalized[mapping[k]][i];
                        res[i] += delta * delta;
                }
        }
}

}

