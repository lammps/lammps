#ifndef PTM_DEFORMATION_GRADIENT_H
#define PTM_DEFORMATION_GRADIENT_H

#include <stdint.h>
#include "ptm_constants.h"

namespace ptm {

void calculate_deformation_gradient(int num_points, const double (*ideal_points)[3], int8_t* mapping, double (*normalized)[3], const double (*penrose)[3], double* F, double* res);

//sc
#define k_sc 0.5
const double penrose_sc[PTM_NUM_POINTS_SC][3] = {	
					{0, 0, 0},
					{0, 0, -k_sc},
					{0, 0, k_sc},
					{0, -k_sc, 0},
					{0, k_sc, 0},
					{-k_sc, 0, 0},
					{k_sc, 0, 0},
				};

//fcc
#define k_fcc 0.17677669529663678216
const double penrose_fcc[PTM_NUM_POINTS_FCC][3] = {
					{0, 0, 0},
					{0, k_fcc, k_fcc},
					{0, -k_fcc, -k_fcc},
					{0, k_fcc, -k_fcc},
					{0, -k_fcc, k_fcc},
					{k_fcc, 0, k_fcc},
					{-k_fcc, 0, -k_fcc},
					{k_fcc, 0, -k_fcc},
					{-k_fcc, 0, k_fcc},
					{k_fcc, k_fcc, -0},
					{-k_fcc, -k_fcc, 0},
					{k_fcc, -k_fcc, 0},
					{-k_fcc, k_fcc, -0},
				};

//hcp
#define k_hcp 0.17677669529663678216
const double penrose_hcp[PTM_NUM_POINTS_HCP][3] = {
					{0, 0, 0},
					{k_hcp, 0, k_hcp},
					{-k_hcp/3, -4*k_hcp/3, -k_hcp/3},
					{k_hcp, k_hcp, 0},
					{-k_hcp/3, -k_hcp/3, -4*k_hcp/3},
					{0, k_hcp, k_hcp},
					{-4*k_hcp/3, -k_hcp/3, -k_hcp/3},
					{-k_hcp, k_hcp, -0},
					{0, k_hcp, -k_hcp},
					{k_hcp, 0, -k_hcp},
					{k_hcp, -k_hcp, 0},
					{-k_hcp, 0, k_hcp},
					{0, -k_hcp, k_hcp},
				};

//ico
#define k_ico 0.13143277802974323576
#define phi 1.61803398874989490253
//((1.0 + sqrt(5)) / 2)
const double penrose_ico[PTM_NUM_POINTS_ICO][3] = {
					{0, 0, 0},
					{0, k_ico, phi*k_ico},
					{0, -k_ico, -phi*k_ico},
					{0, k_ico, -phi*k_ico},
					{0, -k_ico, phi*k_ico},
					{-k_ico, -phi*k_ico, -0},
					{k_ico, phi*k_ico, 0},
					{k_ico, -phi*k_ico, 0},
					{-k_ico, phi*k_ico, -0},
					{-phi*k_ico, 0, -k_ico},
					{phi*k_ico, 0, k_ico},
					{phi*k_ico, 0, -k_ico},
					{-phi*k_ico, 0, k_ico},
				};

//bcc
#define k_bcc 0.11543038598460284017
const double penrose_bcc[PTM_NUM_POINTS_BCC][3] = {
					{0, 0, 0},
					{-k_bcc, -k_bcc, -k_bcc},
					{k_bcc, k_bcc, k_bcc},
					{k_bcc, -k_bcc, -k_bcc},
					{-k_bcc, k_bcc, k_bcc},
					{-k_bcc, k_bcc, -k_bcc},
					{k_bcc, -k_bcc, k_bcc},
					{-k_bcc, -k_bcc, k_bcc},
					{k_bcc, k_bcc, -k_bcc},
					{0, 0, -2*k_bcc},
					{0, 0, 2*k_bcc},
					{0, -2*k_bcc, 0},
					{0, 2*k_bcc, 0},
					{-2*k_bcc, 0, 0},
					{2*k_bcc, 0, -0},
				};

//dcub
#define kdcub 0.07095369570691034689
const double penrose_dcub[PTM_NUM_POINTS_DCUB][3] = {
					{          0,          0,          0 },
					{     -kdcub,      kdcub,      kdcub },
					{     -kdcub,     -kdcub,     -kdcub },
					{      kdcub,     -kdcub,      kdcub },
					{      kdcub,      kdcub,     -kdcub },
					{ -2 * kdcub,          0,  2 * kdcub },
					{ -2 * kdcub,  2 * kdcub,          0 },
					{          0,  2 * kdcub,  2 * kdcub },
					{ -2 * kdcub, -2 * kdcub,          0 },
					{ -2 * kdcub,          0, -2 * kdcub },
					{          0, -2 * kdcub, -2 * kdcub },
					{          0, -2 * kdcub,  2 * kdcub },
					{  2 * kdcub, -2 * kdcub,          0 },
					{  2 * kdcub,          0,  2 * kdcub },
					{          0,  2 * kdcub, -2 * kdcub },
					{  2 * kdcub,          0, -2 * kdcub },
				 	{  2 * kdcub,  2 * kdcub,          0 },
				};


#define kdhex 0.04730246380471011397
const double penrose_dhex[PTM_NUM_POINTS_DHEX][3] = {
					{          0,          0,           0 },
					{     -kdcub,     -kdcub,      -kdcub },
					{      kdcub,     -kdcub,       kdcub },
					{     -kdcub,      kdcub,       kdcub },
					{      kdcub,      kdcub,      -kdcub },
					{     -kdhex, -4 * kdhex,      -kdhex },
					{ -4 * kdhex,     -kdhex,      -kdhex },
					{     -kdhex,     -kdhex,  -4 * kdhex },
					{  2 * kdcub,          0,   2 * kdcub },
					{  2 * kdcub, -2 * kdcub,           0 },
					{          0, -2 * kdcub,   2 * kdcub },
					{          0,  2 * kdcub,   2 * kdcub },
					{ -2 * kdcub,  2 * kdcub,           0 },
					{ -2 * kdcub,          0,   2 * kdcub },
					{  2 * kdcub,  2 * kdcub,           0 },
					{          0,  2 * kdcub,  -2 * kdcub },
					{  2 * kdcub,          0,  -2 * kdcub },
				};
}

#endif


