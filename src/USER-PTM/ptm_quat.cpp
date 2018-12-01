#include <cstring>
#include <cmath>
#include <cfloat>


namespace ptm {

#define SIGN(x) (x >= 0 ? 1 : -1)
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


#define SQRT_2         1.4142135623730951454746218587388284504414
#define HALF_SQRT_2    0.7071067811865474617150084668537601828575

#define PHI            1.6180339887498949025257388711906969547272
#define HALF_PHI       0.8090169943749474512628694355953484773636

#define INV_PHI        0.6180339887498947915034364086750429123640
#define HALF_INV_PHI   0.3090169943749473957517182043375214561820

#define SQRT_5_        2.23606797749978969640917366873127623544061835961152572427089
#define SQRT_2_3       0.8164965809277260344600790631375275552273
#define SQRT_1_6       0.4082482904638630172300395315687637776136


double generator_cubic[24][4] = {		{1,	0,	0,	0	},
						{0,	1,	0,	0	},
						{0,	0,	1,	0	},
						{0,	0,	0,	1	},
						{0.5,	0.5,	0.5,	0.5	},
						{0.5,	0.5,	-0.5,	0.5	},
						{0.5,	-0.5,	0.5,	0.5	},
						{0.5,	-0.5,	-0.5,	0.5	},
						{-0.5,	0.5,	0.5,	0.5	},
						{-0.5,	0.5,	-0.5,	0.5	},
						{-0.5,	-0.5,	0.5,	0.5	},
						{-0.5,	-0.5,	-0.5,	0.5	},
						{HALF_SQRT_2,	HALF_SQRT_2,	0,	0	},
						{HALF_SQRT_2,	0,	HALF_SQRT_2,	0	},
						{HALF_SQRT_2,	0,	0,	HALF_SQRT_2	},
						{-HALF_SQRT_2,	HALF_SQRT_2,	0,	0	},
						{-HALF_SQRT_2,	0,	HALF_SQRT_2,	0	},
						{-HALF_SQRT_2,	0,	0,	HALF_SQRT_2	},
						{0,	HALF_SQRT_2,	HALF_SQRT_2,	0	},
						{0,	HALF_SQRT_2,	0,	HALF_SQRT_2	},
						{0,	0,	HALF_SQRT_2,	HALF_SQRT_2	},
						{0,	-HALF_SQRT_2,	HALF_SQRT_2,	0	},
						{0,	-HALF_SQRT_2,	0,	HALF_SQRT_2	},
						{0,	0,	-HALF_SQRT_2,	HALF_SQRT_2	}	};

double generator_diamond_cubic[12][4] = {	{1,	0,	0,	0	},
						{0,	1,	0,	0	},
						{0,	0,	1,	0	},
						{0,	0,	0,	1	},
						{0.5,	0.5,	0.5,	0.5	},
						{0.5,	0.5,	-0.5,	0.5	},
						{0.5,	-0.5,	0.5,	0.5	},
						{0.5,	-0.5,	-0.5,	0.5	},
						{-0.5,	0.5,	0.5,	0.5	},
						{-0.5,	0.5,	-0.5,	0.5	},
						{-0.5,	-0.5,	0.5,	0.5	},
						{-0.5,	-0.5,	-0.5,	0.5	}	};

double generator_hcp[6][4] = {			{1, 0, 0, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5},
						{0, SQRT_2_3, -SQRT_1_6, -SQRT_1_6},
						{0, SQRT_1_6, -SQRT_2_3, SQRT_1_6},
						{0, SQRT_1_6, SQRT_1_6, -SQRT_2_3}	};

double generator_diamond_hexagonal[3][4] = {	{1, 0, 0, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5}	};

double generator_icosahedral[60][4] = {		{1, 0, 0, 0},
						{HALF_PHI, -HALF_INV_PHI, -0.5, 0},
						{HALF_PHI, 0, -HALF_INV_PHI, -0.5},
						{HALF_PHI, -0.5, 0, -HALF_INV_PHI},
						{HALF_PHI, HALF_INV_PHI, -0.5, 0},
						{HALF_PHI, 0, HALF_INV_PHI, -0.5},
						{HALF_PHI, -0.5, 0, HALF_INV_PHI},
						{HALF_PHI, 0.5, 0, -HALF_INV_PHI},
						{HALF_PHI, 0, -HALF_INV_PHI, 0.5},
						{HALF_PHI, -HALF_INV_PHI, 0.5, 0},
						{HALF_PHI, 0, HALF_INV_PHI, 0.5},
						{HALF_PHI, HALF_INV_PHI, 0.5, 0},
						{HALF_PHI, 0.5, 0, HALF_INV_PHI},
						{0.5, HALF_PHI, -HALF_INV_PHI, 0},
						{0.5, HALF_PHI, HALF_INV_PHI, 0},
						{0.5, 0.5, 0.5, 0.5},
						{0.5, 0.5, 0.5, -0.5},
						{0.5, 0.5, -0.5, 0.5},
						{0.5, 0.5, -0.5, -0.5},
						{0.5, HALF_INV_PHI, 0, HALF_PHI},
						{0.5, HALF_INV_PHI, 0, -HALF_PHI},
						{0.5, 0, HALF_PHI, -HALF_INV_PHI},
						{0.5, 0, HALF_PHI, HALF_INV_PHI},
						{0.5, 0, -HALF_PHI, -HALF_INV_PHI},
						{0.5, 0, -HALF_PHI, HALF_INV_PHI},
						{0.5, -HALF_INV_PHI, 0, HALF_PHI},
						{0.5, -HALF_INV_PHI, 0, -HALF_PHI},
						{0.5, -0.5, 0.5, 0.5},
						{0.5, -0.5, 0.5, -0.5},
						{0.5, -0.5, -0.5, 0.5},
						{0.5, -0.5, -0.5, -0.5},
						{0.5, -HALF_PHI, -HALF_INV_PHI, 0},
						{0.5, -HALF_PHI, HALF_INV_PHI, 0},
						{HALF_INV_PHI, -HALF_PHI, 0, -0.5},
						{HALF_INV_PHI, 0, -0.5, -HALF_PHI},
						{HALF_INV_PHI, -0.5, -HALF_PHI, 0},
						{HALF_INV_PHI, 0, 0.5, -HALF_PHI},
						{HALF_INV_PHI, -HALF_PHI, 0, 0.5},
						{HALF_INV_PHI, 0.5, -HALF_PHI, 0},
						{HALF_INV_PHI, HALF_PHI, 0, -0.5},
						{HALF_INV_PHI, -0.5, HALF_PHI, 0},
						{HALF_INV_PHI, 0, -0.5, HALF_PHI},
						{HALF_INV_PHI, HALF_PHI, 0, 0.5},
						{HALF_INV_PHI, 0, 0.5, HALF_PHI},
						{HALF_INV_PHI, 0.5, HALF_PHI, 0},
						{0, 1, 0, 0},
						{0, HALF_PHI, -0.5, HALF_INV_PHI},
						{0, HALF_PHI, -0.5, -HALF_INV_PHI},
						{0, HALF_PHI, 0.5, HALF_INV_PHI},
						{0, HALF_PHI, 0.5, -HALF_INV_PHI},
						{0, 0.5, HALF_INV_PHI, -HALF_PHI},
						{0, 0.5, HALF_INV_PHI, HALF_PHI},
						{0, 0.5, -HALF_INV_PHI, -HALF_PHI},
						{0, 0.5, -HALF_INV_PHI, HALF_PHI},
						{0, HALF_INV_PHI, -HALF_PHI, 0.5},
						{0, HALF_INV_PHI, -HALF_PHI, -0.5},
						{0, HALF_INV_PHI, HALF_PHI, 0.5},
						{0, HALF_INV_PHI, HALF_PHI, -0.5},
						{0, 0, 1, 0},
						{0, 0, 0, 1}	};

static void quat_rot(double* r, double* a, double* b)
{
	b[0] = (r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]);
	b[1] = (r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]);
	b[2] = (r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]);
	b[3] = (r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]);
}

static int rotate_quaternion_into_fundamental_zone(int num_generators, double (*generator)[4], double* q)
{
	double max = 0.0;
	int i = 0, bi = -1;
	for (i=0;i<num_generators;i++)
	{
		double* g = generator[i];
		double t = fabs(q[0] * g[0] - q[1] * g[1] - q[2] * g[2] - q[3] * g[3]);
		if (t > max)
		{
			max = t;
			bi = i;
		}
	}

	double f[4];
	quat_rot(q, generator[bi], f);
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

int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double* q)
{
	return rotate_quaternion_into_fundamental_zone(3, generator_diamond_hexagonal, q);
}

double quat_dot(double* a, double* b)
{
	return	  a[0] * b[0]
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

	q[0] = sqrt(MAX(0, q[0]));
	q[1] = sqrt(MAX(0, q[1]));
	q[2] = sqrt(MAX(0, q[2]));
	q[3] = sqrt(MAX(0, q[3]));

	double m0 = MAX(q[0], q[1]);
	double m1 = MAX(q[2], q[3]);
	double max = MAX(m0, m1);

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

double quat_quick_misorientation(double* q1, double* q2)
{
	double t = quat_dot(q1, q2);
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_misorientation(double* q1, double* q2)
{
	return acos(quat_quick_misorientation(q1, q2));
}


double quat_quick_disorientation_cubic(double* q0, double* q1)
{
	double qrot[4];
	double qinv[4] = {q0[0], -q0[1], -q0[2], -q0[3]};
	quat_rot(qinv, q1, qrot);

	rotate_quaternion_into_cubic_fundamental_zone(qrot);
	double t = qrot[0];
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_disorientation_cubic(double* q0, double* q1)
{
	return acos(quat_quick_disorientation_cubic(q0, q1));
}

double quat_quick_disorientation_diamond_cubic(double* q0, double* q1)
{
	double qrot[4];
	double qinv[4] = {q0[0], -q0[1], -q0[2], -q0[3]};
	quat_rot(qinv, q1, qrot);

	rotate_quaternion_into_diamond_cubic_fundamental_zone(qrot);
	double t = qrot[0];
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_disorientation_diamond_cubic(double* q0, double* q1)
{
	return acos(quat_quick_disorientation_diamond_cubic(q0, q1));
}

double quat_quick_disorientation_hcp(double* q0, double* q1)
{
	double qrot[4];
	double qinv[4] = {q0[0], -q0[1], -q0[2], -q0[3]};
	quat_rot(qinv, q1, qrot);

	rotate_quaternion_into_hcp_fundamental_zone(qrot);
	double t = qrot[0];
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_disorientation_hcp(double* q0, double* q1)
{
	return acos(quat_quick_disorientation_hcp(q0, q1));
}

double quat_quick_disorientation_diamond_hexagonal(double* q0, double* q1)
{
	double qrot[4];
	double qinv[4] = {q0[0], -q0[1], -q0[2], -q0[3]};
	quat_rot(qinv, q1, qrot);

	rotate_quaternion_into_diamond_hexagonal_fundamental_zone(qrot);
	double t = qrot[0];
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_disorientation_diamond_hexagonal(double* q0, double* q1)
{
	return acos(quat_quick_disorientation_diamond_hexagonal(q0, q1));
}

double quat_quick_disorientation_icosahedral(double* q0, double* q1)
{
	double qrot[4];
	double qinv[4] = {q0[0], -q0[1], -q0[2], -q0[3]};
	quat_rot(qinv, q1, qrot);

	rotate_quaternion_into_icosahedral_fundamental_zone(qrot);
	double t = qrot[0];
	t = MIN(1, MAX(-1, t));
	return 2 * t * t - 1;
}

double quat_disorientation_icosahedral(double* q0, double* q1)
{
	return acos(quat_quick_disorientation_icosahedral(q0, q1));
}

}

