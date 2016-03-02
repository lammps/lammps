#include "comm.h"
#include "manifold_spine_opt.h"
#include "math_const.h"

#include <math.h>

using namespace LAMMPS_NS;

using namespace MathConst;

manifold_spine_opt::manifold_spine_opt( LAMMPS *lmp, int argc, char **argv ) : Pointers(lmp), lu_table(2048)
{
}


manifold_spine_opt::~manifold_spine_opt()
{
}


/*
 * Equation for spine is:
 * 
 * -x^2 - y^2 + (a^2 - (z/c)^2)*( 1 + A*sin(const(z) *z^2) )^4
 * const(z) = (z < 0) ? B2 : B
 * params[5] = { a, A, B, B2, c }
 *
 * These functions are however optimized quite a bit (but still slow =()
 */


double manifold_spine_opt::g( const double *x )
{
  const double B  = params[2];
  const double B2 = params[3];

  const double x2 = x[0]*x[0];
  const double y2 = x[1]*x[1];
  const double z2 = x[2]*x[2];

  // Delay branching up to here, compute everything that is used by both before.
  if( x[2] > 0 ){
    const double z2_c2 = inv_c2 * z2;
    
    const double Bz2 = B * z2;
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = 1.0 + A4*s4_Bz2;
    const double a2_m_z2_c2 = a2 - z2_c2;

    return -x2 - y2 + a2_m_z2_c2 * one_p_A4s4;
    
  }else{
    const double Bz2 = B2 * z2;
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s3_Bz2 = s_Bz2*s2_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = (1.0 + A4*s4_Bz2);
    const double a2_m_z2 = a2 - z2;

    return -x2 - y2 + a2_m_z2 * one_p_A4s4;
  }
}



void manifold_spine_opt::n( const double *x, double *nn )
{
  const double B  = params[2];
  const double B2 = params[3];

	
  nn[0] = -2*x[0];
  nn[1] = -2*x[1];

  const double z2 = x[2]*x[2];

  // Delay branching up to here, compute everything that is used by both before.
  if( x[2] > 0 ){
    const double z2_c2 = inv_c2 * z2;
    const double z_c2  = inv_c2 * x[2];
    
    const double Bz2 = B * z2;
    const double c_Bz2 = lu_table.cos(Bz2);
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s3_Bz2 = s_Bz2*s2_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = 1.0 + A4*s4_Bz2;
    const double a2_m_z2_c2 = a2 - z2_c2;


    nn[2] = 8*A4*B*x[2] * a2_m_z2_c2 * c_Bz2 * s3_Bz2 - 2.0*z_c2 * one_p_A4s4;
    
  }else{
    const double Bz2 = B2 * z2;
    const double c_Bz2 = lu_table.cos(Bz2);
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s3_Bz2 = s_Bz2*s2_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = (1.0 + A4*s4_Bz2);
    const double a2_m_z2 = a2 - z2;
	    
	    
    nn[2] = 8*A4*B2*x[2] * a2_m_z2 * c_Bz2 * s3_Bz2 - 2.0*x[2] * one_p_A4s4;
  }
}

double manifold_spine_opt::g_and_n( const double *x, double *nn )
{
  const double B  = params[2];
  const double B2 = params[3];
	
  nn[0] = -2*x[0];
  nn[1] = -2*x[1];

  const double x2 = x[0]*x[0];
  const double y2 = x[1]*x[1];
  const double z2 = x[2]*x[2];

  // Delay branching up to here, compute everything that is used by both before.
  if( x[2] > 0 ){
    const double z2_c2 = inv_c2 * z2;
    const double z_c2  = inv_c2 * x[2];
    
    const double Bz2 = B * z2;
    const double c_Bz2 = lu_table.cos(Bz2);
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s3_Bz2 = s_Bz2*s2_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = 1.0 + A4*s4_Bz2;
    const double a2_m_z2_c2 = a2 - z2_c2;


    nn[2] = 8*A4*B*x[2] * a2_m_z2_c2 * c_Bz2 * s3_Bz2 - 2.0*z_c2 * one_p_A4s4;
    return -x2 - y2 + a2_m_z2_c2 * one_p_A4s4;
    
  }else{
    const double Bz2 = B2 * z2;
    const double c_Bz2 = lu_table.cos(Bz2);
    const double s_Bz2 = lu_table.sin(Bz2);
    const double s2_Bz2 = s_Bz2*s_Bz2;
    const double s3_Bz2 = s_Bz2*s2_Bz2;
    const double s4_Bz2 = s2_Bz2*s2_Bz2;
    const double one_p_A4s4 = (1.0 + A4*s4_Bz2);
    const double a2_m_z2 = a2 - z2;
	    
	    
    nn[2] = 8*A4*B2*x[2] * a2_m_z2 * c_Bz2 * s3_Bz2 - 2.0*x[2] * one_p_A4s4;
    return -x2 - y2 + a2_m_z2 * one_p_A4s4;
  }
}


void manifold_spine_opt::post_param_init()
{
  const double a  = params[0];
  const double A  = params[1];
  const double c  = params[4];
  
  c2 = c*c;
  inv_c2 = 1.0 / c2;
  a2 = a*a;
  const double A2 = A*A;
  A4 = A2*A2;
}

