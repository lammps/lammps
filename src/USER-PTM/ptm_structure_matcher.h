#ifndef PTM_STRUCTURE_MATCHER_H
#define PTM_STRUCTURE_MATCHER_H

#include "ptm_initialize_data.h"
#include "ptm_constants.h"


namespace ptm {

typedef struct
{
	double rmsd;
	double scale;
	double q[4];		//rotation in quaternion form (rigid body transformation)
	int8_t mapping[PTM_MAX_POINTS];
	const refdata_t* ref_struct;
} result_t;

int match_general(const refdata_t* s, double (*ch_points)[3], double (*points)[3], convexhull_t* ch, result_t* res);
int match_fcc_hcp_ico(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);
int match_dcub_dhex(double (*ch_points)[3], double (*points)[3], int32_t flags, convexhull_t* ch, result_t* res);

}

#endif

