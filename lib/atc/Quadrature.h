/** Quadrature : creat Gaussian quadrature lists  */
#ifndef QUADRATURE_H
#define QUADRATURE_H

using namespace std;

static const int line_ngauss = 10;
static double line_xg[line_ngauss], line_wg[line_ngauss];
/** domain of integration is -1 to 1 */
static void set_line_quadrature(int ng, double* xg, double* wg)
{
  /** integration schemes : 3, 4, 5, 10 point Gauss */
    
  if (ng == 3) {
    xg[0] = 0.0; wg[0] = 8.0/9.0;
    xg[1] = -sqrt(3.0/5.0); wg[1] = 5.0/9.0;
    xg[2] =  sqrt(3.0/5.0); wg[2] = 5.0/9.0;
  }
  else if (ng == 4) {
    xg[0] = -sqrt((3.0/7.0)-(sqrt(120.0)/35.0)); wg[0] = (1.0/2.0)+(5.0/(3.0*sqrt(120.0)));
    xg[1] =  sqrt((3.0/7.0)-(sqrt(120.0)/35.0)); wg[1] = (1.0/2.0)+(5.0/(3.0*sqrt(120.0)));
    xg[2] = -sqrt((3.0/7.0)+(sqrt(120.0)/35.0)); wg[2] = (1.0/2.0)-(5.0/(3.0*sqrt(120.0)));
    xg[3] =  sqrt((3.0/7.0)+(sqrt(120.0)/35.0)); wg[3] = (1.0/2.0)-(5.0/(3.0*sqrt(120.0)));
  }
  else if (ng == 5) {
    xg[0] = 0.0; wg[0] = 0.5688888888888889;
    xg[1] = -sqrt((35.0-sqrt(280.0))/63.0); wg[1] = 0.4786286704993658;
    xg[2] =  sqrt((35.0-sqrt(280.0))/63.0); wg[2] = 0.4786286704993658;
    xg[3] = -sqrt((35.0+sqrt(280.0))/63.0); wg[3] = 0.2369268850561891;
    xg[4] =  sqrt((35.0+sqrt(280.0))/63.0); wg[4] = 0.2369268850561891;
  }
  else if (ng == 10) {
    xg[0] = -0.14887434; wg[0] = 0.29552422;
    xg[1] =  0.14887434; wg[1] = 0.29552422;
    xg[2] = -0.43339539; wg[2] = 0.26926672;
    xg[3] =  0.43339539; wg[3] = 0.26926672;
    xg[4] = -0.67940957; wg[4] = 0.21908636;
    xg[5] =  0.67940957; wg[5] = 0.21908636;
    xg[6] = -0.86506337; wg[6] = 0.14945135;
    xg[7] =  0.86506337; wg[7] = 0.14945135;
    xg[8] = -0.97390653; wg[8] = 0.06667134;
    xg[9] =  0.97390653; wg[9] = 0.06667134;
  }
  else { cout << "Invalid choice of number of Gauss points for Quadrature" << endl; }
  //else { throw ATC_Error(0,"Invalid choice of number of Gauss points for Quadrature"); }

  return;
};

#endif
