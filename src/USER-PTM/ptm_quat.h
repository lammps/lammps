#ifndef PTM_QUAT_H
#define PTM_QUAT_H

namespace ptm {

int rotate_quaternion_into_cubic_fundamental_zone(double* q);
int rotate_quaternion_into_diamond_cubic_fundamental_zone(double* q);
int rotate_quaternion_into_icosahedral_fundamental_zone(double* q);
int rotate_quaternion_into_hcp_fundamental_zone(double* q);
int rotate_quaternion_into_diamond_hexagonal_fundamental_zone(double* q);

void normalize_quaternion(double* q);
void quaternion_to_rotation_matrix(double* q, double* U);
void rotation_matrix_to_quaternion(double* u, double* q);
double quat_dot(double* a, double* b);
double quat_quick_misorientation(double* q1, double* q2);
double quat_misorientation(double* q1, double* q2);

double quat_quick_disorientation_cubic(double* q0, double* q1);
double quat_disorientation_cubic(double* q0, double* q1);
double quat_quick_disorientation_diamond_cubic(double* q0, double* q1);
double quat_disorientation_diamond_cubic(double* q0, double* q1);
double quat_quick_disorientation_hcp(double* q0, double* q1);
double quat_disorientation_hcp(double* q0, double* q1);
double quat_quick_disorientation_diamond_hexagonal(double* q0, double* q1);
double quat_disorientation_diamond_hexagonal(double* q0, double* q1);
double quat_quick_disorientation_icosahedral(double* q0, double* q1);
double quat_disorientation_icosahedral(double* q0, double* q1);

}

#endif

