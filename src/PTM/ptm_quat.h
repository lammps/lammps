/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef PTM_QUAT_H
#define PTM_QUAT_H

namespace ptm {

int rotate_quaternion_into_cubic_fundamental_zone(double *q);
int rotate_quaternion_into_diamond_cubic_fundamental_zone(double *q);
int rotate_quaternion_into_icosahedral_fundamental_zone(double *q);
int rotate_quaternion_into_hcp_fundamental_zone(double *q);
int rotate_quaternion_into_hcp_conventional_fundamental_zone(double *q);
int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double *q);

void normalize_quaternion(double *q);
void quaternion_to_rotation_matrix(double *q, double *U);
void rotation_matrix_to_quaternion(double *u, double *q);
double quat_dot(double *a, double *b);
double quat_misorientation(double *q1, double *q2);

}    // namespace ptm

#endif
