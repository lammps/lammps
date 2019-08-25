/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstring>
#include <cmath>
#include <cfloat>
#include <algorithm>


namespace ptm {

#define SIGN(x) (x >= 0 ? 1 : -1)


const double generator_cubic[24][4] = {
        {          1,          0,          0,          0 },
        {  sqrt(2)/2,  sqrt(2)/2,          0,          0 },
        {  sqrt(2)/2,          0,  sqrt(2)/2,          0 },
        {  sqrt(2)/2,          0,          0,  sqrt(2)/2 },
        {  sqrt(2)/2,          0,          0, -sqrt(2)/2 },
        {  sqrt(2)/2,          0, -sqrt(2)/2,          0 },
        {  sqrt(2)/2, -sqrt(2)/2,         -0,         -0 },
        {        0.5,        0.5,        0.5,        0.5 },
        {        0.5,        0.5,        0.5,       -0.5 },
        {        0.5,        0.5,       -0.5,        0.5 },
        {        0.5,        0.5,       -0.5,       -0.5 },
        {        0.5,       -0.5,        0.5,        0.5 },
        {        0.5,       -0.5,        0.5,       -0.5 },
        {        0.5,       -0.5,       -0.5,        0.5 },
        {        0.5,       -0.5,       -0.5,       -0.5 },
        {          0,          1,          0,          0 },
        {          0,  sqrt(2)/2,  sqrt(2)/2,          0 },
        {          0,  sqrt(2)/2,          0,  sqrt(2)/2 },
        {          0,  sqrt(2)/2,          0, -sqrt(2)/2 },
        {          0,  sqrt(2)/2, -sqrt(2)/2,          0 },
        {          0,          0,          1,          0 },
        {          0,          0,  sqrt(2)/2,  sqrt(2)/2 },
        {          0,          0,  sqrt(2)/2, -sqrt(2)/2 },
        {          0,          0,          0,          1 },
};

const double generator_diamond_cubic[12][4] = {
        {    1,    0,    0,    0 },
        {  0.5,  0.5,  0.5,  0.5 },
        {  0.5,  0.5,  0.5, -0.5 },
        {  0.5,  0.5, -0.5,  0.5 },
        {  0.5,  0.5, -0.5, -0.5 },
        {  0.5, -0.5,  0.5,  0.5 },
        {  0.5, -0.5,  0.5, -0.5 },
        {  0.5, -0.5, -0.5,  0.5 },
        {  0.5, -0.5, -0.5, -0.5 },
        {    0,    1,    0,    0 },
        {    0,    0,    1,    0 },
        {    0,    0,    0,    1 },
};

const double generator_hcp[6][4] = {
        {          1,          0,          0,          0 },
        {        0.5,          0,          0,  sqrt(3)/2 },
        {        0.5,          0,          0, -sqrt(3)/2 },
        {          0,  sqrt(3)/2,        0.5,          0 },
        {          0,  sqrt(3)/2,       -0.5,          0 },
        {          0,          0,          1,          0 },
};


const double generator_hcp_conventional[12][4] = {
        {          1,          0,          0,          0 },
        {  sqrt(3)/2,          0,          0,        0.5 },
        {  sqrt(3)/2,          0,          0,       -0.5 },
        {        0.5,          0,          0,  sqrt(3)/2 },
        {        0.5,          0,          0, -sqrt(3)/2 },
        {          0,          1,          0,          0 },
        {          0,  sqrt(3)/2,        0.5,          0 },
        {          0,  sqrt(3)/2,       -0.5,          0 },
        {          0,        0.5,  sqrt(3)/2,          0 },
        {          0,        0.5, -sqrt(3)/2,          0 },
        {          0,          0,          1,          0 },
        {          0,          0,          0,          1 },
};

const double generator_diamond_hexagonal[3][4] = {
        {          1,          0,          0,          0 },
        {        0.5,          0,          0,  sqrt(3)/2 },
        {        0.5,          0,          0, -sqrt(3)/2 },
};

const double generator_icosahedral[60][4] = {
        {                        1,                        0,                        0,                        0 },
        {            (1+sqrt(5))/4,                      0.5,   sqrt(25-10*sqrt(5))/10,   sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,                      0.5,  -sqrt(25-10*sqrt(5))/10,  -sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,            1/(1+sqrt(5)),   sqrt(10*sqrt(5)+50)/20,  -sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,            1/(1+sqrt(5)),  -sqrt(10*sqrt(5)+50)/20,   sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,                        0,   sqrt(50-10*sqrt(5))/10,   sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,                        0,                        0,     sqrt(5./8-sqrt(5)/8) },
        {            (1+sqrt(5))/4,                        0,                        0,    -sqrt(5./8-sqrt(5)/8) },
        {            (1+sqrt(5))/4,                        0,  -sqrt(50-10*sqrt(5))/10,  -sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,           -1/(1+sqrt(5)),   sqrt(10*sqrt(5)+50)/20,  -sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,           -1/(1+sqrt(5)),  -sqrt(10*sqrt(5)+50)/20,   sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,                     -0.5,   sqrt(25-10*sqrt(5))/10,   sqrt(50-10*sqrt(5))/20 },
        {            (1+sqrt(5))/4,                     -0.5,  -sqrt(25-10*sqrt(5))/10,  -sqrt(50-10*sqrt(5))/20 },
        {                      0.5,            (1+sqrt(5))/4,   sqrt(50-10*sqrt(5))/20,  -sqrt(25-10*sqrt(5))/10 },
        {                      0.5,            (1+sqrt(5))/4,  -sqrt(50-10*sqrt(5))/20,   sqrt(25-10*sqrt(5))/10 },
        {                      0.5,                      0.5,   sqrt((5+2*sqrt(5))/20),   sqrt(25-10*sqrt(5))/10 },
        {                      0.5,                      0.5,   sqrt(25-10*sqrt(5))/10,  -sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                      0.5,  -sqrt(25-10*sqrt(5))/10,   sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                      0.5,  -sqrt((5+2*sqrt(5))/20),  -sqrt(25-10*sqrt(5))/10 },
        {                      0.5,            1/(1+sqrt(5)),   sqrt(10*sqrt(5)+50)/20,   sqrt((5+2*sqrt(5))/20) },
        {                      0.5,            1/(1+sqrt(5)),  -sqrt(10*sqrt(5)+50)/20,  -sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                        0,     sqrt((5+sqrt(5))/10),  -sqrt(25-10*sqrt(5))/10 },
        {                      0.5,                        0,   sqrt(50-10*sqrt(5))/10,  -sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                        0,  -sqrt(50-10*sqrt(5))/10,   sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                        0,    -sqrt((5+sqrt(5))/10),   sqrt(25-10*sqrt(5))/10 },
        {                      0.5,           -1/(1+sqrt(5)),   sqrt(10*sqrt(5)+50)/20,   sqrt((5+2*sqrt(5))/20) },
        {                      0.5,           -1/(1+sqrt(5)),  -sqrt(10*sqrt(5)+50)/20,  -sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                     -0.5,   sqrt((5+2*sqrt(5))/20),   sqrt(25-10*sqrt(5))/10 },
        {                      0.5,                     -0.5,   sqrt(25-10*sqrt(5))/10,  -sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                     -0.5,  -sqrt(25-10*sqrt(5))/10,   sqrt((5+2*sqrt(5))/20) },
        {                      0.5,                     -0.5,  -sqrt((5+2*sqrt(5))/20),  -sqrt(25-10*sqrt(5))/10 },
        {                      0.5,           -(1+sqrt(5))/4,   sqrt(50-10*sqrt(5))/20,  -sqrt(25-10*sqrt(5))/10 },
        {                      0.5,           -(1+sqrt(5))/4,  -sqrt(50-10*sqrt(5))/20,   sqrt(25-10*sqrt(5))/10 },
        {            1/(1+sqrt(5)),            (1+sqrt(5))/4,   sqrt(50-10*sqrt(5))/20,   sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),            (1+sqrt(5))/4,  -sqrt(50-10*sqrt(5))/20,  -sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                      0.5,   sqrt((5+2*sqrt(5))/20),  -sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                      0.5,  -sqrt((5+2*sqrt(5))/20),   sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                        0,     sqrt((5+sqrt(5))/10),   sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                        0,                        0,  sqrt(1-1/(2*sqrt(5)+6)) },
        {            1/(1+sqrt(5)),                        0,                        0, -sqrt(1-1/(2*sqrt(5)+6)) },
        {            1/(1+sqrt(5)),                        0,    -sqrt((5+sqrt(5))/10),  -sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                     -0.5,   sqrt((5+2*sqrt(5))/20),  -sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),                     -0.5,  -sqrt((5+2*sqrt(5))/20),   sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),           -(1+sqrt(5))/4,   sqrt(50-10*sqrt(5))/20,   sqrt(10*sqrt(5)+50)/20 },
        {            1/(1+sqrt(5)),           -(1+sqrt(5))/4,  -sqrt(50-10*sqrt(5))/20,  -sqrt(10*sqrt(5)+50)/20 },
        {                        0,                        1,                        0,                        0 },
        {                        0,            (1+sqrt(5))/4,     sqrt(5./8-sqrt(5)/8),                        0 },
        {                        0,            (1+sqrt(5))/4,   sqrt(50-10*sqrt(5))/20,  -sqrt(50-10*sqrt(5))/10 },
        {                        0,            (1+sqrt(5))/4,  -sqrt(50-10*sqrt(5))/20,   sqrt(50-10*sqrt(5))/10 },
        {                        0,            (1+sqrt(5))/4,    -sqrt(5./8-sqrt(5)/8),                        0 },
        {                        0,                      0.5,   sqrt((5+2*sqrt(5))/20),   sqrt(50-10*sqrt(5))/10 },
        {                        0,                      0.5,   sqrt(25-10*sqrt(5))/10,     sqrt((5+sqrt(5))/10) },
        {                        0,                      0.5,  -sqrt(25-10*sqrt(5))/10,    -sqrt((5+sqrt(5))/10) },
        {                        0,                      0.5,  -sqrt((5+2*sqrt(5))/20),  -sqrt(50-10*sqrt(5))/10 },
        {                        0,            1/(1+sqrt(5)),  sqrt(1-1/(2*sqrt(5)+6)),                        0 },
        {                        0,            1/(1+sqrt(5)),   sqrt(10*sqrt(5)+50)/20,    -sqrt((5+sqrt(5))/10) },
        {                        0,            1/(1+sqrt(5)),  -sqrt(10*sqrt(5)+50)/20,     sqrt((5+sqrt(5))/10) },
        {                        0,            1/(1+sqrt(5)), -sqrt(1-1/(2*sqrt(5)+6)),                        0 },
        {                        0,                        0,     sqrt((5+sqrt(5))/10),  -sqrt(50-10*sqrt(5))/10 },
        {                        0,                        0,   sqrt(50-10*sqrt(5))/10,     sqrt((5+sqrt(5))/10) },
};


static void quat_rot(double* r, double* a, double* b)
{
        b[0] = (r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]);
        b[1] = (r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]);
        b[2] = (r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]);
        b[3] = (r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]);
}

static int rotate_quaternion_into_fundamental_zone(int num_generators, const double (*generator)[4], double* q)
{
        double max = 0.0;
        int i = 0, bi = -1;
        for (i=0;i<num_generators;i++)
        {
                const double* g = generator[i];
                double t = fabs(q[0] * g[0] - q[1] * g[1] - q[2] * g[2] - q[3] * g[3]);
                if (t > max)
                {
                        max = t;
                        bi = i;
                }
        }

        double f[4];
        quat_rot(q, (double*)generator[bi], f);
        memcpy(q, &f, 4 * sizeof(double));
        if (q[0] < 0)
        {
                q[0] = -q[0];
                q[1] = -q[1];
                q[2] = -q[2];
                q[3] = -q[3];
        }

        return bi;
}

int rotate_quaternion_into_cubic_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(24, generator_cubic, q);
}

int rotate_quaternion_into_diamond_cubic_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(12, generator_diamond_cubic, q);
}

int rotate_quaternion_into_icosahedral_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(60, generator_icosahedral, q);
}

int rotate_quaternion_into_hcp_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(6, generator_hcp, q);
}

int rotate_quaternion_into_hcp_conventional_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(12, generator_hcp_conventional, q);
}

int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double* q)
{
        return rotate_quaternion_into_fundamental_zone(3, generator_diamond_hexagonal, q);
}

double quat_dot(double* a, double* b)
{
        return          a[0] * b[0]
                + a[1] * b[1]
                + a[2] * b[2]
                + a[3] * b[3];
}

double quat_size(double* q)
{
        return sqrt(quat_dot(q, q));
}

void normalize_quaternion(double* q)
{
        double size = quat_size(q);

        q[0] /= size;
        q[1] /= size;
        q[2] /= size;
        q[3] /= size;
}

void quaternion_to_rotation_matrix(double* q, double* u)
{
        double a = q[0];
        double b = q[1];
        double c = q[2];
        double d = q[3];

        u[0] = a*a + b*b - c*c - d*d;
        u[1] = 2*b*c - 2*a*d;
        u[2] = 2*b*d + 2*a*c;

        u[3] = 2*b*c + 2*a*d;
        u[4] = a*a - b*b + c*c - d*d;
        u[5] = 2*c*d - 2*a*b;

        u[6] = 2*b*d - 2*a*c;
        u[7] = 2*c*d + 2*a*b;
        u[8] = a*a - b*b - c*c + d*d;
}

void rotation_matrix_to_quaternion(double* u, double* q)
{
        double r11 = u[0];
        double r12 = u[1];
        double r13 = u[2];
        double r21 = u[3];
        double r22 = u[4];
        double r23 = u[5];
        double r31 = u[6];
        double r32 = u[7];
        double r33 = u[8];

        q[0] = (1.0 + r11 + r22 + r33) / 4.0;
        q[1] = (1.0 + r11 - r22 - r33) / 4.0;
        q[2] = (1.0 - r11 + r22 - r33) / 4.0;
        q[3] = (1.0 - r11 - r22 + r33) / 4.0;

        q[0] = sqrt(std::max(0., q[0]));
        q[1] = sqrt(std::max(0., q[1]));
        q[2] = sqrt(std::max(0., q[2]));
        q[3] = sqrt(std::max(0., q[3]));

        double m0 = std::max(q[0], q[1]);
        double m1 = std::max(q[2], q[3]);
        double max = std::max(m0, m1);

        int i = 0;
        for (i=0;i<4;i++)
                if (q[i] == max)
                        break;

        if (i == 0)
        {
                q[1] *= SIGN(r32 - r23);
                q[2] *= SIGN(r13 - r31);
                q[3] *= SIGN(r21 - r12);
        }
        else if (i == 1)
        {
                q[0] *= SIGN(r32 - r23);
                q[2] *= SIGN(r21 + r12);
                q[3] *= SIGN(r13 + r31);
        }
        else if (i == 2)
        {
                q[0] *= SIGN(r13 - r31);
                q[1] *= SIGN(r21 + r12);
                q[3] *= SIGN(r32 + r23);
        }
        else if (i == 3)
        {
                q[0] *= SIGN(r21 - r12);
                q[1] *= SIGN(r31 + r13);
                q[2] *= SIGN(r32 + r23);
        }

        normalize_quaternion(q);
}

double quat_quick_misorientation(double* q1, double* q2)
{
        double t = quat_dot(q1, q2);
        t = std::min(1., std::max(-1., t));
        return 2 * t * t - 1;
}

double quat_misorientation(double* q1, double* q2)
{
        return acos(quat_quick_misorientation(q1, q2));
}

}

