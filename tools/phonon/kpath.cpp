
#include "kpath.h"

#include "global.h"
#include "dynmat.h"
#include "memory.h"
#include "qnodes.h"

#ifdef UseSPG
extern "C"{
#include "spglib.h"
}
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

/* ----------------------------------------------------------------------------
 * Class kPath will help to find the high symmetry k-path for a given lattice.
 * The lattice info is passed via a DynMat object.
 * ---------------------------------------------------------------------------- */
kPath::kPath(DynMat *dm, QNodes *qn)
{
   // create memory 
   memory = new Memory();
 
   // pass the class from main
   q = qn;
   dynmat = dm;
   sysdim = dynmat->sysdim;
   num_atom = dynmat->nucell;
   
   memory->create(atpos, num_atom, 3, "kpath:atpos");
   memory->create(attyp, num_atom,    "kpath:attyp");
   
   // set default, in case system dimension under study is not 3.
   for (int i = 0; i < num_atom; ++i)
   for (int idim = 0; idim < 3; ++idim) atpos[i][idim] = 0.;
 
   for (int i = 0; i < 3; ++i)
   for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;
 
   for (int i = 0; i < 3; ++i) latvec[i][i] = 1.;
 
   // get atomic type info
   for (int i = 0; i < num_atom; ++i) attyp[i] = dynmat->attyp[i];
 
   // get unit cell vector info
   int ndim = 0;
   for (int idim = 0; idim < 3; ++idim)
   for (int jdim = 0; jdim < 3; ++jdim) latvec[jdim][idim] = dynmat->basevec[ndim++];
 
   // get atom position in unit cell; fractional
   for (int i = 0; i < num_atom; ++i)
   for (int idim = 0; idim < sysdim; ++idim) atpos[i][idim] = dynmat->basis[i][idim];
 
   // get the space group number
   double symprec = 1.0e-3;
   double **pos;
   memory->create(pos,num_atom,3,"kpath:pos");
   if (dynmat->symprec > 0.0) symprec = dynmat->symprec;

   for (int i = 0; i < num_atom; ++i)
   for (int j = 0; j < 3; ++j) pos[i][j] = atpos[i][j];
   spgnum  = spg_get_international(symbol, latvec, (double (*)[3])pos, attyp, num_atom, symprec);
   memory->destroy(pos);
}

/* ----------------------------------------------------------------------------
 * Show the lattice info
 * ---------------------------------------------------------------------------- */
void kPath::show_info()
{
   // display the unit cell info read
   puts("--------------------------------------------------------------------------------");
   printf("The basis vectors of the unit cell:\n");
   for (int idim = 0; idim < 3; ++idim){
      printf("  A%d =", idim+1);
      for (int jdim = 0; jdim < 3; ++jdim) printf(" %lg", latvec[jdim][idim]);
      printf("\n");
   }
 
   printf("Atom(s) in the unit cell:\n");
   printf("  No.  type  sx  sy sz\n");
   for (int i = 0; i < MIN(num_atom, NUMATOM); ++i)
      printf("  %d %d %lg %lg %lg\n", i+1, attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
   if (num_atom > NUMATOM) printf("  ... (%d atoms omitted.)\n", num_atom-NUMATOM);
 
   printf("The space group number of your unit cell is: %d => %s\n", spgnum, symbol);
   puts("--------------------------------------------------------------------------------");
}

/* ----------------------------------------------------------------------------
 * Free the memeory used by kPath.
 * ---------------------------------------------------------------------------- */
kPath::~kPath( )
{
   memory->destroy(attyp);
   memory->destroy(atpos);
   delete memory;
   dynmat = NULL;
   q = NULL;
}
  
/* ----------------------------------------------------------------------------
 * Main workhorse
 * ---------------------------------------------------------------------------- */
void kPath::kpath( )
{
  q->ndstr.clear();
  q->qs.clear();
  q->qe.clear();
  q->nqbin.clear();
  
  // angles
  const double la = sqrt(latvec[0][0] * latvec[0][0] + latvec[0][1] * latvec[0][1] + latvec[0][2] * latvec[0][2]);
  const double lb = sqrt(latvec[1][0] * latvec[1][0] + latvec[1][1] * latvec[1][1] + latvec[1][2] * latvec[1][2]);
  const double lc = sqrt(latvec[2][0] * latvec[2][0] + latvec[2][1] * latvec[2][1] + latvec[2][2] * latvec[2][2]);
  const double cosa = sqrt(latvec[1][0] * latvec[2][0] + latvec[1][1] * latvec[2][1] + latvec[1][2] * latvec[2][2])/(lb*lc);
  const double cosg = sqrt(latvec[0][0] * latvec[1][0] + latvec[0][1] * latvec[1][1] + latvec[0][2] * latvec[1][2])/(la*lb);
  double ivec[3][3];
  int ndim = 0;
  for (int idim = 0; idim < 3; ++idim)
  for (int jdim = 0; jdim < 3; ++jdim) ivec[jdim][idim] = dynmat->ibasevec[ndim++];
  const double ka = sqrt(ivec[0][0] * ivec[0][0] + ivec[0][1] * ivec[0][1] + ivec[0][2] * ivec[0][2]);
  const double kb = sqrt(ivec[1][0] * ivec[1][0] + ivec[1][1] * ivec[1][1] + ivec[1][2] * ivec[1][2]);
  const double coskg = sqrt(ivec[0][0] * ivec[1][0] + ivec[0][1] * ivec[1][1] + ivec[0][2] * ivec[1][2])/(ka*kb);

  double *qtmp;
  if (spgnum >= 1 && spgnum <= 2){  // Triclinic system
    if (fabs(coskg) > ZERO){        // A.14, TRI1a and TRI2a

      q->ndstr.push_back("X");
      // X-G
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] =  qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-Y
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Y/L");

      // L-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z/N");

      // N-G
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-M
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("M/R");

      // R-G
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = qtmp[1] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

    } else {        // A.14, TRI1b and TRI2b

      q->ndstr.push_back("X");
      // X-G
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = -0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] =  qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-Y
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[1] = qtmp[2] = 0.; qtmp[0] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Y/L");

      // L-G
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z/N");

      // N-G
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = -0.5; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-M
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("M/R");

      // R-G
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = -0.5; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");
    }

  } else if (spgnum >= 3 && spgnum <= 15){ // Monoclinic

    if (symbol[0] == 'P'){             // MCL-P
      const double eta = (1.-lb*cosa/lc)/(2.*(1.-cosa*cosa));
      const double niu = 0.5 - eta * lc * cosa / lb;

      q->ndstr.push_back("{/Symbol G}");
      // G-Y
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Y");

      // Y-H
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = eta; qtmp[2] = 1.-niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("H");

      // H-C
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = eta; qtmp[2] = 1.-niu;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("C");

      // C-E
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("E");

      // E-M1
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("M_1");

      // M1-A
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = niu;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("A");

      // A-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("X");

      // X-H1
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[2] = 1.-eta; qtmp[1] = niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("H_1/M");

      // M-D
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[2] = eta; qtmp[1] = 1.-niu;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("D");

      // D-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z/Y");

      // Y-D
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("D");

    } else {                  // MCL-C

      if (coskg < 0.){        // MCLC1
        const double xi = (2. - lb * cosa / lc) / (4.*(1.-cosa*cosa));
        const double eta = 0.5 + 2.*xi*lc/lb*cosa;
        const double psi = 0.75 - la * la /(4.*lb*lb*(1.-cosa*cosa));
        const double phi = psi + (0.75 - psi) * lb / lc * cosa;

        q->ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("Y");

        // Y-F
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("F");

        // F-L
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("L");

        // L-I
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = phi; qtmp[1] = 1.-phi; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("I/I_1");

        // I1-Z
        qtmp = new double [3];
        qtmp[0] = 1.-phi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("Z");

        // Z-F1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("F_1/Y");

        // Y-X1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = psi; qtmp[1] = 1.-psi; qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("X_1/X");

        // X-G
        qtmp = new double [3];
        qtmp[0] = 1.-psi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("{/Symbol G}");

        // G-N
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("N/M");

        // M-G
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("{/Symbol G}");

      } else if ( fabs(coskg) < ZERO) {     // MCLC2
        const double xi = (2. - lb * cosa / lc) / (4.*(1.-cosa*cosa));
        const double eta = 0.5 + 2.*xi*lc/lb*cosa;
        const double psi = 0.75 - la * la /(4.*lb*lb*(1.-cosa*cosa));
        const double phi = psi + (0.75 - psi) * lb / lc * cosa;

        q->ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("Y");

        // Y-F
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("F");

        // F-L
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("L");

        // L-I
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = phi; qtmp[1] = 1.-phi; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("I/I_1");

        // I1-Z
        qtmp = new double [3];
        qtmp[0] = 1.-phi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("Z");

        // Z-F1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("F_1/N");

        // N-G
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("{/Symbol G}");

        // G-M
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);

        q->ndstr.push_back("M");

      } else {

        double flag = lb / lc * cosa + lb * lb / (la * la *(1.-cosa*cosa));
        if (fabs(flag) < ZERO){           // MCLC4
          const double miu = 0.25*(1. + lb * lb / (la *la));
          const double del = lb * lc * cosa / (2.*la*la);
          const double xi  = miu - 0.25 + (1. - lb * cosa / lc)/(4.*(1.-cosa*cosa));
          const double eta = 0.5 + 2.*xi*lc/lb*cosa;
          const double phi = 1. + xi - 2.*miu;
          const double psi = eta - 2.*del;
 
          q->ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y");

          // Y-F
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("F");

          // F-H
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("H");

          // H-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Z");

          // Z-I
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("I/H_1");

          // H1-Y1
          qtmp = new double [3];
          qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y_1");

          // Y1-X
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("X");

          // X-G
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");

          // G-N
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("N/M");

          // M-G
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");

        } else if (flag < 1.){            // MCLC3
          const double miu = 0.25*(1. + lb * lb / (la *la));
          const double del = lb * lc * cosa / (2.*la*la);
          const double xi  = miu - 0.25 + (1. - lb * cosa / lc)/(4.*(1.-cosa*cosa));
          const double eta = 0.5 + 2.*xi*lc/lb*cosa;
          const double phi = 1. + xi - 2.*miu;
          const double psi = eta - 2.*del;
 
          q->ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y");

          // Y-F
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("F");

          // F-H
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("H");

          // H-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Z");

          // Z-I
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("I/H_1");

          // I-F1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = phi; qtmp[2] = phi - 1.; qtmp[1] = psi;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("F_1/H_1");

          // H1-Y1
          qtmp = new double [3];
          qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y_1");

          // Y1-X
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("X");

          // X-G
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");

          // G-N
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("N/M");

          // M-G
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");

        } else {                          // MCLC5
          const double xi = (lb*lb/(la*la) + (1.-lb*cosa/lc)/(1.-cosa*cosa))*0.25;
          const double eta = 0.5 + 2.*xi*lc*cosa/lb;
          const double miu = 0.5*eta + lb * lb /(4.*la*la) - lb*lc/(2.*la*la)*cosa;
          const double niu = 2.*miu - xi;
          const double omg = (4.*niu - 1. - lb*lb*(1.-cosa*cosa)/(la*la))*lc/(2.*lb*cosa);
          const double del = xi*lc*cosa/lb + omg*0.5 - 0.25;
          const double rho = 1.-xi*la*la/(lb*lb);
 
          q->ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y");

          // Y-F
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = niu; qtmp[2] = omg;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("F");

          // F-L
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = niu; qtmp[2] = omg;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y");

          // L-I
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = rho; qtmp[1] = 1.-rho; qtmp[2] = 0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("I/I_1");

          // I1-Z
          qtmp = new double [3];
          qtmp[0] = 1.-rho; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Z");

          // Z-H
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("H");

          // H-F1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-omg;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("F_1/H_1");

          // H1-Y1
          qtmp = new double [3];
          qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y_1");

          // Y1-X
          qtmp = new double [3];
          qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("Y_1");

          // X-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");

          // G-N
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("N/M");

          // M-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          q->qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          q->qe.push_back(qtmp);
  
          q->ndstr.push_back("{/Symbol G}");
        }
      }
    }

  } else if (spgnum >= 16 && spgnum <= 74){ // Orthorhombic

    if (symbol[0] == 'P'){             // ORC

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X");

      // X-S
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("S");

      // S-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("{/Symbol G}");

      // G-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Z");

      // Z-U
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("U");

      // U-R
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("R");

      // R-T
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("T");

      // T-Z
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Z/Y");

      // Y-T
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("T/U");

      // U-X
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X/S");

      // S-R
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("R");

    } else if (symbol[0] == 'F'){        // ORCF

      double flag = 1./(la*la) - 1./(lb*lb) - 1./(lc*lc);
      if (fabs(flag) < ZERO){            // ORCF3

        const double xi  = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
        const double eta = 0.25 * (1. + la*la/(lb*lb) + la*la/(lc*lc));

        q->ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y");

        // Y-T
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("T");

        // T-Z
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z");

        // Z-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");

        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("X");

        // X-A1
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("A_1");

        // A1-Y
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y/X");

        // X-A
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("A");

        // A-Z
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z/L");

        // L-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");

      } else if (flag > 0.){             // ORCF1

        const double xi  = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
        const double eta = 0.25 * (1. + la*la/(lb*lb) + la*la/(lc*lc));

        q->ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y");

        // Y-T
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("T");

        // T-Z
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z");

        // Z-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");

        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("X");

        // X-A1
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("A_1");

        // A1-Y
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y/T");

        // T-X1
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 1.; qtmp[1] = qtmp[2] = 1.-eta;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("X_1/X");

        // X-A
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("A");

        // A-Z
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z/L");

        // L-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");

      } else {                           // ORCF2

        const double eta = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
        const double phi = 0.25 * (1. + lc*lc/(lb*lb) - lc*lc/(la*la));
        const double del = 0.25 * (1. + lb*lb/(la*la) - lb*lb/(lc*lc));

        q->ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y");

        // Y-C
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[2] = 0.5-eta; qtmp[1] = 1.-eta;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("C");

        // C-D
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[2] = 0.5-eta; qtmp[1] = 1.-eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5-del; qtmp[2] = 0.5; qtmp[1] = 1.-del;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("D");

        // D-X
        qtmp = new double [3];
        qtmp[0] = 0.5-del; qtmp[2] = 0.5; qtmp[1] = 1.-del;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("X");

        // X-G
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z");

        // Z-D1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5+del; qtmp[1] = 0.5; qtmp[2] = del;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("D_1");

        // D1-H
        qtmp = new double [3];
        qtmp[0] = 0.5+del; qtmp[1] = 0.5; qtmp[2] = del;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("H");

        // H-C
        qtmp = new double [3];
        qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5-eta; qtmp[2] = 1.-eta;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("C/C_1");

        // C1-Z
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.5+eta; qtmp[2] = eta;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Z/X");

        // X-H1
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = phi; qtmp[1] = 0.5+phi; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("H_1/H");

        // H-Y
        qtmp = new double [3];
        qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("Y/L");

        // L-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
    
        q->ndstr.push_back("{/Symbol G}");
      }

    } else if (symbol[0] == 'C'){        // ORCC

      const double xi = 0.25 * (1. + la*la/(lb*lb));

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X");

      // X-S
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("S");

      // S-R
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("R");

      // R-A
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("A");

      // A-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Z");

      // Z-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("{/Symbol G}");

      // G-Y
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Y");

      // Y-X1
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X_1");

      // X1-A1
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("A_1");

      // A1-T
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("T");

      // T-Y
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Y/Z");

      // Z-T
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("T");

    } else {                             // ORCI

      const double xi  = 0.25 * (1. + la*la/(lc*lc));
      const double eta = 0.25 * (1. + lb*lb/(lc*lc));
      const double del = (lb*lb-la*la)/(4.*lc*lc);
      const double miu = (lb*lb+la*la)/(4.*lc*lc);

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = qtmp[2] = xi;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X");

      // X-L
      qtmp = new double [3];
      qtmp[0] = -xi; qtmp[1] = qtmp[2] = xi;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = -miu; qtmp[1] = miu; qtmp[2] = 0.5-del;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("L");

      // L-T
      qtmp = new double [3];
      qtmp[0] = -miu; qtmp[1] = miu; qtmp[2] = 0.5-del;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("T");

      // T-W
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("W");

      // W-R
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("R");

      // R-X1
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = xi; qtmp[1] = 1.-xi; qtmp[2] = -xi;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("X_1");

      // X1-Z
      qtmp = new double [3];
      qtmp[0] = xi; qtmp[1] = 1.-xi; qtmp[2] = -xi;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Z");

      // Z-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("{/Symbol G}");

      // G-Y
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Y");

      // Y-S
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("S");

      // S-W
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("W/L_1");

      // L1-Y
      qtmp = new double [3];
      qtmp[0] = miu; qtmp[1] = -miu; qtmp[2] = 0.5+del;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Y/Y_1");

      // Y1-Z
      qtmp = new double [3];
      qtmp[0] = 1.-eta; qtmp[1] = eta; qtmp[2] = -eta;
      q->qs.push_back(qtmp);
  
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
      q->qe.push_back(qtmp);
  
      q->ndstr.push_back("Z");
    }

  } else if (spgnum >= 75 && spgnum <= 142){ // Tetragonal

    if (symbol[0] == 'P'){             // TET

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("X");

      // X-M
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("M");

      // M-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z");

      // Z-R
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("R");

      // R-A
      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("A");

      // A-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z/X");

      // X-R
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("R/M");

      // M-A
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("A");

    } else {                           // BCT
      if (la > lc){                    // BCT1
        const double eta = 0.25 * (1. + lc*lc/(la*la));

        q->ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("X");

        // X-M
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("M");

        // M-G
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = eta; qtmp[2] = -eta;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Z");

        // Z-P
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = eta; qtmp[2] = -eta;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("P");

        // P-N
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("N");

        // N-Z1
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = -eta; qtmp[1] = 1.-eta; qtmp[2] = eta;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Z_1");

        // Z1-M
        qtmp = new double [3];
        qtmp[0] = -eta; qtmp[1] = 1.-eta; qtmp[2] = eta;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("M");

        // X-P
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("P");

      } else {                         // BCT2

        const double eta = 0.25 * (1. + la*la/(lc*lc));
        const double xi  = la*la/(2.*lc*lc);

        q->ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("X");

        // X-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = xi; qtmp[2] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Y");

        // Y-Sigma
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = xi; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = -eta; qtmp[1] = qtmp[2] = eta;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("{/Symbol S}");

        // Sigma-G
        qtmp = new double [3];
        qtmp[0] = -eta; qtmp[1] = qtmp[2] = eta;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Z");

        // Z-Sigma_1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = eta; qtmp[1] = 1.-eta; qtmp[2] = -eta;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("{/Symbol S}_1");

        // Sigma_1-N
        qtmp = new double [3];
        qtmp[0] = eta; qtmp[1] = 1.-eta; qtmp[2] = -eta;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("N");

        // N-P
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("P");

        // P-Y1
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -xi;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Y_1");

        // Y1-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -xi;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("Z");

        // X-P
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        q->qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
        q->qe.push_back(qtmp);
 
        q->ndstr.push_back("P");
      }
    }
  } else if (spgnum >= 143 && spgnum <= 167){ // Trigonal
    if (cosg > 0.){                           // RHL1

      const double eta = (1.+4.*cosa)/(2.+4.*cosa);
      const double niu = 0.75 - 0.5*eta;

      q->ndstr.push_back("{/Symbol G}");
      // G-L
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("L");

      // L-B1
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = eta - 1.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("B_1/B");

      // B-Z
      qtmp = new double [3];
      qtmp[0] = eta; qtmp[1] = 0.5; qtmp[2] = 1.-eta;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z");

      // Z-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = niu; qtmp[1] = 0.; qtmp[2] = -niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("X/Q");

      // Q-F
      qtmp = new double [3];
      qtmp[0] = 1.-niu; qtmp[1] = niu; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("F");

      // F-P1
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-eta;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("P_1");

      // P1-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-eta;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z/L");

      // L-P
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = eta; qtmp[1] = qtmp[2] = niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("P");

    } else {                                  // RHL2
      const double eta = 0.5*(1.+cosa)/(1.-cosa);
      const double niu = 0.75 - 0.5*eta;

      q->ndstr.push_back("{/Symbol G}");
      // G-P
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 1.-niu; qtmp[1] = -niu;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("P");

      // P-Z
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 1.-niu; qtmp[1] = -niu;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z");

      // Z-Q
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = eta;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Q");

      // Q-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = eta;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("{/Symbol G}");

      // G-F
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("F");

      // F-P1
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = niu; qtmp[1] = qtmp[2] = niu - 1.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("P_1");

      // P1-Q1
      qtmp = new double [3];
      qtmp[0] = niu; qtmp[1] = qtmp[2] = niu - 1.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 1.-eta; qtmp[1] = qtmp[2] = -eta;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Q_1");

      // Q1-L
      qtmp = new double [3];
      qtmp[0] = 1.-eta; qtmp[1] = qtmp[2] = -eta;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("L");

      // L-Z
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);

      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qe.push_back(qtmp);

      q->ndstr.push_back("Z");
    }

  } else if (spgnum >= 168 && spgnum <= 194){ // Hexagonal

    q->ndstr.push_back("{/Symbol G}");
    // G-M
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = qtmp[2] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("M");
   
    // M-K
    qtmp = new double [3];
    qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("K");
   
    // K-G
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = qtmp[2] = 0.;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("{/Symbol G}");
   
    // G-A
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = qtmp[2] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("A");
   
    // A-L
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("L");
   
    // L-H
    qtmp = new double [3];
    qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("H");
   
    // H-A
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("A/L");
   
    // L-M
    qtmp = new double [3];
    qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("M/K");
   
    // K-H
    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
    q->qs.push_back(qtmp);

    qtmp = new double [3];
    qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
    q->qe.push_back(qtmp);

    q->ndstr.push_back("H");
   
  } else if (spgnum >= 195 && spgnum <= 230){ // Cubic

    if (symbol[0] == 'P'){                    // CUB

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("X");
   
      // X-M
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("M");
   
      // M-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("{/Symbol G}");
   
      // G-R
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("R");
   
      // R-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("X/M");
   
      // M-R
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("R");

    } else if (symbol[0] == 'F'){          // FCC

      q->ndstr.push_back("{/Symbol G}");
      // G-X
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("X");
   
      // X-W
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("W");
   
      // W-K
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("K");
   
      // K-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("{/Symbol G}");
   
      // G-L
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("L");
   
      // L-U
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("U");
   
      // U-W
      qtmp = new double [3];
      qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("W");
   
      // W-L
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("L");
   
      // L-K
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("K/U");
   
      // U-X
      qtmp = new double [3];
      qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("X");
   
    } else {                               // BCC

      q->ndstr.push_back("{/Symbol G}");
      // G-H
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("H");
   
      // H-N
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("N");
   
      // N-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("{/Symbol G}");
   
      // G-P
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("P");
   
      // P-H
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("H/P");
   
      // P-N
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
      q->qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      q->qe.push_back(qtmp);
 
      q->ndstr.push_back("N");
   
    }

  } else {
    printf("\nSorry, failed to identify the crystal system, please use the manual mode.\n");
  }
 
   // to determine the number of points along each line, with a step size of 0.05
   const double qs_inv = 1./QSTEP;
   int nbin = q->qs.size();
   for (int is = 0; is < nbin; ++is){
      double *qstr = q->qs[is];
      double *qend = q->qe[is];
      double ql = 0.;
      for (int i = 0; i < 3; ++i) ql += (qend[i] - qstr[i])*(qend[i] - qstr[i]);
      int nqpt = MAX(int(sqrt(ql) * qs_inv + 0.5), 2);
      q->nqbin.push_back(nqpt);
   }

   }

/* ----------------------------------------------------------------------------
 * Show the k-path info
 * ---------------------------------------------------------------------------- */
void kPath::show_path()
{
   if (q == nullptr) return;
   int nbin = q->ndstr.size();
   if (nbin > 0){
      puts("\n--------------------------------------------------------------------------------");
      printf("k-path for the current lattice will be:\n  %s", q->ndstr[0].c_str());
      for (int is = 1; is < nbin; ++is) printf("-%s", q->ndstr[is].c_str());

      printf("\n\nThe fractional coordinates of these paths are:\n");
      for (int is = 0; is < nbin-1; ++is)
         printf("  [%6.4f %6.4f %6.4f] --> [%6.4f %6.4f %6.4f] (%s - %s)\n", q->qs[is][0],
                q->qs[is][1], q->qs[is][2], q->qe[is][0], q->qe[is][1], q->qe[is][2],
                q->ndstr[is].c_str(), q->ndstr[is+1].c_str() );
      puts("--------------------------------------------------------------------------------");
   }

   }
/* ---------------------------------------------------------------------------- */
#endif
