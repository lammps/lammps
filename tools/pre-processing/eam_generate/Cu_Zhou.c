#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

#define EOK 0
#define ERROR -1

#define re 2.556162
#define fe 1.554485
#define rhoe 22.150141
#define alpha 7.669911
#define beta 4.090619
#define A 0.327584
#define B 0.468735
#define kappa  0.431307
#define lambda 0.86214
#define Fn0 -2.17649
#define Fn1 -0.140035
#define Fn2 0.285621
#define Fn3 -1.750834
#define F0 -2.19
#define F1 0.
#define F2 0.702991
#define F3 0.683705
#define eta 0.921150
#define Fe -2.191675

double V (double r) {
   return ( 
         ( A*exp (-alpha * (r/re-1.) ) )  /  (1. + pow (r/re-kappa, 20.))
         -
         ( B*exp (-beta * (r/re-1.) ) )  /  (1. + pow (r/re-lambda, 20.))
         );
}

double rho (double r) {
   return (
         (fe * exp (-beta * (r/re-1.)))  /  (1. + pow (r/re-lambda, 20.))
         );
}

double F (double rho_) {
   double rhon = .85*rhoe,
          rho0 = 1.15*rhoe;
   if (rho_ < rhon)
      return ( 
            Fn0 * pow (rho_/rhon-1., 0.) +
            Fn1 * pow (rho_/rhon-1., 1.) +
            Fn2 * pow (rho_/rhon-1., 2.) +
            Fn3 * pow (rho_/rhon-1., 3.)
            );
   else if (rhon <= rho_ && rho_ < rho0)
      return (
            F0 * pow (rho_/rhoe-1., 0.) +
            F1 * pow (rho_/rhoe-1., 1.) +
            F2 * pow (rho_/rhoe-1., 2.) +
            F3 * pow (rho_/rhoe-1., 3.) 
            );
   else if (rho0 <= rho_)
      return (
            Fe*(1. - log ( pow (rho_/rhoe, eta) )) * pow (rho_/rhoe, eta)
            );
}

int main (void) {
   int Nr = 10001;
   double rmax = 3.615*2.5;
   double dr = rmax/(double)Nr;
   int Nrho = 10001;
   double rhomax = rho (0.);
   double drho = rhomax/(double)Nrho;

   int atomic_number = 1;
   double mass = 63.55;
   double lattice_constant = 3.615;
   char lattice_type[] = "FCC";

   int i;

   char LAMMPSFilename[] = "Cu_Zhou.eam";
   FILE *LAMMPSFile = fopen (LAMMPSFilename, "w");
   if (!LAMMPSFile) exit (ERROR);

   // Header for setfl format
   fprintf (LAMMPSFile, \
         "#-> LAMMPS Potential File in DYNAMO 86 setfl Format <-#\n"\
         "# Zhou Cu Acta mater(2001)49:4005\n"\
         "# Implemented by G. Ziegenhain (2007) gerolf@ziegenhain.com\n"\
         "%d Cu\n"\
         "%d %20.20f %d %20.20f %20.20f\n"\
         "%d %20.20f %20.20f %s\n",
         atomic_number, 
         Nrho, drho, Nr, dr, rmax,
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
