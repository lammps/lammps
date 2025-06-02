#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

#define rc  5.50679// A
#define h   0.50037// A
#define E1  2.01458*1e2//eV
#define E2  6.59288*1e-3//eV
#define r01 0.83591//A
#define r02 4.46867//A
#define a1  2.97758//1/A
#define a2  1.54927//1/A
#define d   0.86225*1e-2//A
#define rs1 2.24//A
#define rs2 1.8//A
#define rs3 1.2//A
#define S1  4.//eV/A^4
#define S2  40.//eV/A^4
#define S3  1.15*1e3//eV/A^4
#define a   3.80362
#define r03 -2.19885//A
#define r04 -2.61984*1e2//A
#define b1  0.17394//A^-2
#define b2  5.35661*1e2//1/A
#define F0  -2.28235//eV
#define F2  1.35535//eV
#define q1  -1.27775//eV
#define q2  -0.86074//eV
#define q3  1.78804//eV
#define q4  2.97571//eV
#define Q1  0.4
#define Q2  0.3

#define EOK 0
#define ERROR -1


double M (double r, double r0, double alpha) {
   return exp (-2.*alpha*(r-r0)) - 2.*exp (-alpha*(r-r0));
}

double psi (double x) {
   if (x >= 0.) return 0.;
   else return pow (x, 4.) / (1.+pow (x, 4.) );
}

double H (double x) {
   if (x >= 0) return 1.;
   else return 0.;
}

double V (double r) {
   return ( E1*M(r,r01,a1) + E2*M(r,r02,a2) + d ) * psi ( (r-rc)/h ) -
      H (rs1-r)*S1*pow (rs1-r, 4.) -
      H (rs2-r)*S2*pow (rs2-r, 4.) -
      H (rs3-r)*S3*pow (rs3-r, 4.);
}

double rho (double r) {
   return ( a*exp (-b1*pow (r-r03, 2.)) + exp (-b2*(r-r04)) ) * psi ( (r-rc)/h );
}

double F (double rho_) {
   if (rho_ < 1.)
      return F0 + .5*F2*pow (rho_-1., 2.) 
         + q1*pow (rho_-1.,3.)
         + q2*pow (rho_-1.,4.)
         + q3*pow (rho_-1.,5.)
         + q4*pow (rho_-1.,6.);
   else if (rho_ > 1.)
      return ( F0 + .5*F2*pow (rho_-1., 2.) + q1*pow (rho_-1., 3.) + Q1*pow (rho_-1., 4.)) /
         ( 1. + Q2*pow (rho_-1., 3.) );
   else exit (ERROR);

}

int main (void) {
   int Nr = 10001;
   double rmax = 9.;
   double dr = rmax/(double)Nr;
   int Nrho = 10001;
   double rhomax = rho (0.);
   double drho = rhomax/(double)Nrho;

   int atomic_number = 1;
   double mass = 63.55;
   double lattice_constant = 3.615;
   char lattice_type[] = "FCC";

   int i;

   char LAMMPSFilename[] = "Cu_Mishin1.eam";
   FILE *LAMMPSFile = fopen (LAMMPSFilename, "w");
   if (!LAMMPSFile) exit (ERROR);

   // Header for setfl format
   fprintf (LAMMPSFile, \
         "#-> LAMMPS Potential File in DYNAMO 86 setfl Format <-#\n"\
         "# Mishin Cu EAM1 PRB(2001)63:224106\n"\
         "# Implemented by G. Ziegenhain (2007) gerolf@ziegenhain.com \n"\
         "%d Cu\n"\
         "%d %20.20f %d %20.20f %20.20f\n"\
         "%d %20.20f %20.20f %s\n",
         atomic_number, 
         Nrho, drho, Nr, dr, rc,
         atomic_number, mass, lattice_constant, lattice_type);

   // Embedding function
   for (i = 0; i < Nrho; i++) 
      fprintf (LAMMPSFile, "%20.20f\n", F ((double)i*drho));
   // Density function
   for (i = 0; i < Nr; i++) 
      fprintf (LAMMPSFile, "%20.20f\n", rho ((double)i*dr));
   // Pair potential
   for (i = 0; i < Nr; i++)   
      fprintf (LAMMPSFile, "%20.20f\n", V ((double)i*dr) * (double)i*dr);

   fclose (LAMMPSFile);
   return (EOK);
}
