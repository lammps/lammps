#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

#define EOK 0
#define ERROR -1

#define re 2.74084
#define fe 3.48734
#define rhoe 37.234847
#define alpha 8.900114
#define beta 4.746728
#define A 0.882435
#define B 1.394592
#define kappa 0.139209
#define lambda 0.278417
#define Fn0 -4.946281
#define Fn1 -0.148818
#define Fn2 0.365057
#define Fn3 -4.432406
#define F0 -4.96
#define F1 0.
#define F2 0.661935
#define F3 0.348147
#define eta -0.582714
#define Fe -4.961306

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
            Fe*(1. - log( pow (rho_/rhoe, eta) ) ) * pow (rho_/rhoe, eta)
            );
}

int main (void) {
   int Nr = 10001;
   double rmax = 2.5*3.157;
   double dr = rmax/(double)Nr;
   int Nrho = 10001;
   double rhomax = rho (0.);
   double drho = rhomax/(double)Nrho;

   int atomic_number = 1;
   double mass = 183.84;
   double lattice_constant = 3.157;
   char lattice_type[] = "BCC";

   int i;

   char LAMMPSFilename[] = "W_Zhou.eam";
   FILE *LAMMPSFile = fopen (LAMMPSFilename, "w");
   if (!LAMMPSFile) exit (ERROR);

   // Header for setfl format
   fprintf (LAMMPSFile, \
         "#-> LAMMPS Potential File in DYNAMO 86 setfl Format <-#\n"\
         "# Zhou W Acta mater(2001)49:4005\n"\
         "# Implemented by G. Ziegenhain (2007) gerolf@ziegenhain.com\n"\
         "%d W\n"\
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
