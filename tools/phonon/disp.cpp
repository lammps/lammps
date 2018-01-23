#include <vector>
#include "string.h"
#include "phonon.h"
#include "green.h"
#include "timer.h"
#include "global.h"

#ifdef UseSPG
extern "C"{
#include "spglib.h"
}
#endif

/*------------------------------------------------------------------------------
 * Private method to evaluate the phonon dispersion curves
 *----------------------------------------------------------------------------*/
void Phonon::pdisp()
{
  // ask the output file name and write the header.
  char str[MAXLINE];
  for (int ii = 0; ii < 80; ++ii) printf("=");
  printf("\n");
#ifdef UseSPG
  // ask method to generate q-lines
  int method = 2;
  printf("Please select your method to generate the phonon dispersion:\n");
  printf("  1. Manual, should always work;\n");
  printf("  2. Automatic, works only for 3D crystals (CMS49-299).\nYour choice [2]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) method = atoi(strtok(str," \t\n\r\f"));
  method = 2 - method%2;
  printf("Your  selection: %d\n", method);
#endif
  printf("\nPlease input the filename to output the dispersion data [pdisp.dat]:");
  if (count_words(fgets(str,MAXLINE,stdin)) < 1) strcpy(str, "pdisp.dat");
  char *ptr = strtok(str," \t\n\r\f");
  char *fname = new char[strlen(ptr)+1];
  strcpy(fname,ptr);

  // to store the nodes of the dispersion curve
  std::vector<double> nodes;        nodes.clear();
  std::vector<std::string> ndstr;   ndstr.clear();
  std::vector<double *> qs, qe;     qs.clear(); qe.clear();
  std::vector<int> nqbin;           nqbin.clear();

  // now the calculate the dispersion curve
  double qstr[3], qend[3];
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2+1;
  qend[0] = qend[1] = qend[2] = 0.;

  double *egvs = new double [ndim];
#ifdef UseSPG
  if (method == 1){
#endif
    while (1){
      for (int i = 0; i < 3; ++i) qstr[i] = qend[i];
  
      printf("\nPlease input the start q-point in unit of B1->B3, q to exit [%g %g %g]: ", qstr[0], qstr[1], qstr[2]);
      int n = count_words(fgets(str, MAXLINE, stdin));
      ptr = strtok(str, " \t\n\r\f");
      if ((n == 1) && (strcmp(ptr,"q") == 0)) break;
      else if (n >= 3){
        qstr[0] = atof(ptr);
        qstr[1] = atof(strtok(NULL, " \t\n\r\f"));
        qstr[2] = atof(strtok(NULL, " \t\n\r\f"));
      }
  
      do printf("Please input the end q-point in unit of B1->B3: ");
      while (count_words(fgets(str, MAXLINE, stdin)) < 3);
      qend[0] = atof(strtok(str,  " \t\n\r\f"));
      qend[1] = atof(strtok(NULL, " \t\n\r\f"));
      qend[2] = atof(strtok(NULL, " \t\n\r\f"));
  
      printf("Please input the # of points along the line [%d]: ", nq);
      if (count_words(fgets(str, MAXLINE, stdin)) > 0) nq = atoi(strtok(str," \t\n\r\f"));
      nq = MAX(nq,2);
  
      double *qtmp = new double [3];
      for (int i = 0; i < 3; ++i) qtmp[i] = qstr[i];
      qs.push_back(qtmp);

      qtmp = new double [3];
      for (int i = 0; i < 3; ++i) qtmp[i] = qend[i];
      qe.push_back(qtmp);

      nqbin.push_back(nq);
      ndstr.push_back("");
    }
    ndstr.push_back("");

#ifdef UseSPG
  } else {

    memory->grow(atpos, dynmat->nucell, 3, "pdisp:atpos");
    memory->grow(attyp, dynmat->nucell,    "pdisp:attyp");
  
    num_atom = dynmat->nucell;
    // set default, in case system dimension under study is not 3.
    for (int i = 0; i < dynmat->nucell; ++i)
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

    // display the unit cell info read
    printf("\n");
    for (int ii = 0; ii < 80; ++ii) printf("-"); printf("\n");
    printf("The basis vectors of the unit cell:\n");
    for (int idim = 0; idim < 3; ++idim) printf("  A%d = %lg %lg %lg\n", idim+1, latvec[0][idim], latvec[1][idim], latvec[2][idim]);
    printf("Atom(s) in the unit cell:\n");
    printf("  No.  type  sx  sy sz\n");
    for (int i = 0; i < MIN(num_atom, NUMATOM); ++i) printf("  %d %d %lg %lg %lg\n", i+1, attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
    if (num_atom > NUMATOM) printf("  ... (%d atoms omitted.)\n", num_atom-NUMATOM);

    char symbol[11];
    double symprec = 1.e-4, pos[num_atom][3];
    for (int i = 0; i < num_atom; ++i)
    for (int j = 0; j < 3; ++j) pos[i][j] = atpos[i][j];

    int spgnum  = spg_get_international(symbol, latvec, pos, attyp, num_atom, symprec);
    printf("The space group number of your unit cell is: %d => %s\n", spgnum, symbol);
    for (int ii = 0; ii < 80; ++ii) printf("-"); printf("\n");
  
    // angles
    const double la = sqrt(latvec[0][0] * latvec[0][0] + latvec[0][1] * latvec[0][1] + latvec[0][2] * latvec[0][2]);
    const double lb = sqrt(latvec[1][0] * latvec[1][0] + latvec[1][1] * latvec[1][1] + latvec[1][2] * latvec[1][2]);
    const double lc = sqrt(latvec[2][0] * latvec[2][0] + latvec[2][1] * latvec[2][1] + latvec[2][2] * latvec[2][2]);
    const double cosa = sqrt(latvec[1][0] * latvec[2][0] + latvec[1][1] * latvec[2][1] + latvec[1][2] * latvec[2][2])/(lb*lc);
    const double cosg = sqrt(latvec[0][0] * latvec[1][0] + latvec[0][1] * latvec[1][1] + latvec[0][2] * latvec[1][2])/(la*lb);
    double ivec[3][3];
    ndim = 0;
    for (int idim = 0; idim < 3; ++idim)
    for (int jdim = 0; jdim < 3; ++jdim) ivec[jdim][idim] = dynmat->ibasevec[ndim++];
    const double ka = sqrt(ivec[0][0] * ivec[0][0] + ivec[0][1] * ivec[0][1] + ivec[0][2] * ivec[0][2]);
    const double kb = sqrt(ivec[1][0] * ivec[1][0] + ivec[1][1] * ivec[1][1] + ivec[1][2] * ivec[1][2]);
    const double coskg = sqrt(ivec[0][0] * ivec[1][0] + ivec[0][1] * ivec[1][1] + ivec[0][2] * ivec[1][2])/(ka*kb);
    
    double *qtmp;
    if (spgnum >= 1 && spgnum <= 2){  // Triclinic system
      if (fabs(coskg) > ZERO){        // A.14, TRI1a and TRI2a

        ndstr.push_back("X");
        // X-G
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] =  qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Y/L");

        // L-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Z/N");

        // N-G
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-M
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("M/R");

        // R-G
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = qtmp[1] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

      } else {        // A.14, TRI1b and TRI2b

        ndstr.push_back("X");
        // X-G
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = -0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] =  qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[1] = qtmp[2] = 0.; qtmp[0] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Y/L");

        // L-G
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Z/N");

        // N-G
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = -0.5; qtmp[1] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-M
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("M/R");

        // R-G
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = -0.5; qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");
      }

    } else if (spgnum >= 3 && spgnum <= 15){ // Monoclinic

      if (symbol[0] == 'P'){             // MCL-P
        const double eta = (1.-lb*cosa/lc)/(2.*(1.-cosa*cosa));
        const double niu = 0.5 - eta * lc * cosa / lb;

        ndstr.push_back("{/Symbol G}");
        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Y");

        // Y-H
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = eta; qtmp[2] = 1.-niu;
        qe.push_back(qtmp);

        ndstr.push_back("H");

        // H-C
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = eta; qtmp[2] = 1.-niu;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("C");

        // C-E
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("E");

        // E-M1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = niu;
        qe.push_back(qtmp);

        ndstr.push_back("M_1");

        // M1-A
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = niu;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("A");

        // A-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("X");

        // X-H1
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[2] = 1.-eta; qtmp[1] = niu;
        qe.push_back(qtmp);

        ndstr.push_back("H_1/M");

        // M-D
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[2] = eta; qtmp[1] = 1.-niu;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("D");

        // D-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("Z/Y");

        // Y-D
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("D");

      } else {                  // MCL-C

        if (coskg < 0.){        // MCLC1
          const double xi = (2. - lb * cosa / lc) / (4.*(1.-cosa*cosa));
          const double eta = 0.5 + 2.*xi*lc/lb*cosa;
          const double psi = 0.75 - la * la /(4.*lb*lb*(1.-cosa*cosa));
          const double phi = psi + (0.75 - psi) * lb / lc * cosa;

          ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("Y");
  
          // Y-F
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
          qe.push_back(qtmp);
  
          ndstr.push_back("F");
  
          // F-L
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("L");
  
          // L-I
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = phi; qtmp[1] = 1.-phi; qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("I/I_1");
  
          // I1-Z
          qtmp = new double [3];
          qtmp[0] = 1.-phi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("Z");
  
          // Z-F1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          qe.push_back(qtmp);
  
          ndstr.push_back("F_1/Y");
  
          // Y-X1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = psi; qtmp[1] = 1.-psi; qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("X_1/X");
  
          // X-G
          qtmp = new double [3];
          qtmp[0] = 1.-psi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("{/Symbol G}");
  
          // G-N
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("N/M");
  
          // M-G
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("{/Symbol G}");
  
        } else if ( fabs(coskg) < ZERO) {     // MCLC2
          const double xi = (2. - lb * cosa / lc) / (4.*(1.-cosa*cosa));
          const double eta = 0.5 + 2.*xi*lc/lb*cosa;
          const double psi = 0.75 - la * la /(4.*lb*lb*(1.-cosa*cosa));
          const double phi = psi + (0.75 - psi) * lb / lc * cosa;

          ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("Y");
  
          // Y-F
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
          qe.push_back(qtmp);
  
          ndstr.push_back("F");
  
          // F-L
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 1.-xi; qtmp[2] = 1.-eta;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("L");
  
          // L-I
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = phi; qtmp[1] = 1.-phi; qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("I/I_1");
  
          // I1-Z
          qtmp = new double [3];
          qtmp[0] = 1.-phi; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qe.push_back(qtmp);
  
          ndstr.push_back("Z");
  
          // Z-F1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
          qe.push_back(qtmp);
  
          ndstr.push_back("F_1/N");
  
          // N-G
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("{/Symbol G}");
  
          // G-M
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
  
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
  
          ndstr.push_back("M");
  
        } else {

          double flag = lb / lc * cosa + lb * lb / (la * la *(1.-cosa*cosa));
          if (fabs(flag) < ZERO){           // MCLC4
            const double miu = 0.25*(1. + lb * lb / (la *la));
            const double del = lb * lc * cosa / (2.*la*la);
            const double xi  = miu - 0.25 + (1. - lb * cosa / lc)/(4.*(1.-cosa*cosa));
            const double eta = 0.5 + 2.*xi*lc/lb*cosa;
            const double phi = 1. + xi - 2.*miu;
            const double psi = eta - 2.*del;
   
            ndstr.push_back("{/Symbol G}");
            // G-Y
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y");

            // Y-F
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
            qe.push_back(qtmp);
    
            ndstr.push_back("F");

            // F-H
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qe.push_back(qtmp);
    
            ndstr.push_back("H");

            // H-Z
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("Z");

            // Z-I
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("I/H_1");

            // H1-Y1
            qtmp = new double [3];
            qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y_1");

            // Y1-X
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("X");

            // X-G
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");

            // G-N
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("N/M");

            // M-G
            qtmp = new double [3];
            qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");

          } else if (flag < 1.){            // MCLC3
            const double miu = 0.25*(1. + lb * lb / (la *la));
            const double del = lb * lc * cosa / (2.*la*la);
            const double xi  = miu - 0.25 + (1. - lb * cosa / lc)/(4.*(1.-cosa*cosa));
            const double eta = 0.5 + 2.*xi*lc/lb*cosa;
            const double phi = 1. + xi - 2.*miu;
            const double psi = eta - 2.*del;
   
            ndstr.push_back("{/Symbol G}");
            // G-Y
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y");

            // Y-F
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
            qe.push_back(qtmp);
    
            ndstr.push_back("F");

            // F-H
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 1.-phi; qtmp[2] = 1.-psi;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qe.push_back(qtmp);
    
            ndstr.push_back("H");

            // H-Z
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("Z");

            // Z-I
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("I/H_1");

            // I-F1
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = phi; qtmp[2] = phi - 1.; qtmp[1] = psi;
            qe.push_back(qtmp);
    
            ndstr.push_back("F_1/H_1");

            // H1-Y1
            qtmp = new double [3];
            qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y_1");

            // Y1-X
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("X");

            // X-G
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");

            // G-N
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("N/M");

            // M-G
            qtmp = new double [3];
            qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");

          } else {                          // MCLC5
            const double xi = (lb*lb/(la*la) + (1.-lb*cosa/lc)/(1.-cosa*cosa))*0.25;
            const double eta = 0.5 + 2.*xi*lc*cosa/lb;
            const double miu = 0.5*eta + lb * lb /(4.*la*la) - lb*lc/(2.*la*la)*cosa;
            const double niu = 2.*miu - xi;
            const double omg = (4.*niu - 1. - lb*lb*(1.-cosa*cosa)/(la*la))*lc/(2.*lb*cosa);
            const double del = xi*lc*cosa/lb + omg*0.5 - 0.25;
            const double rho = 1.-xi*la*la/(lb*lb);
   
            ndstr.push_back("{/Symbol G}");
            // G-Y
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y");

            // Y-F
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = miu; qtmp[2] = del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = niu; qtmp[2] = omg;
            qe.push_back(qtmp);
    
            ndstr.push_back("F");

            // F-L
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = niu; qtmp[2] = omg;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y");

            // L-I
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = rho; qtmp[1] = 1.-rho; qtmp[2] = 0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("I/I_1");

            // I1-Z
            qtmp = new double [3];
            qtmp[0] = 1.-rho; qtmp[1] = -qtmp[0]; qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qe.push_back(qtmp);
    
            ndstr.push_back("Z");

            // Z-H
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qe.push_back(qtmp);
    
            ndstr.push_back("H");

            // H-F1
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = xi; qtmp[2] = eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-omg;
            qe.push_back(qtmp);
    
            ndstr.push_back("F_1/H_1");

            // H1-Y1
            qtmp = new double [3];
            qtmp[0] = 1.-xi; qtmp[1] = -xi; qtmp[2] = 1.-eta;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y_1");

            // Y1-X
            qtmp = new double [3];
            qtmp[0] = 1.-miu; qtmp[1] = -miu; qtmp[2] = -del;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("Y_1");

            // X-G
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");

            // G-N
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("N/M");

            // M-G
            qtmp = new double [3];
            qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
            qs.push_back(qtmp);
    
            qtmp = new double [3];
            qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
            qe.push_back(qtmp);
    
            ndstr.push_back("{/Symbol G}");
          }
        }
      }

    } else if (spgnum >= 16 && spgnum <= 74){ // Orthorhombic

      if (symbol[0] == 'P'){             // ORC

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("X");

        // X-S
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("S");

        // S-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("Z");

        // Z-U
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("U");

        // U-R
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("R");

        // R-T
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("T");

        // T-Z
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("Z/Y");

        // Y-T
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("T/U");

        // U-X
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("X/S");

        // S-R
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("R");

      } else if (symbol[0] == 'F'){        // ORCF

        double flag = 1./(la*la) - 1./(lb*lb) - 1./(lc*lc);
        if (fabs(flag) < ZERO){            // ORCF3

          const double xi  = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
          const double eta = 0.25 * (1. + la*la/(lb*lb) + la*la/(lc*lc));

          ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y");

          // Y-T
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
      
          ndstr.push_back("T");

          // T-Z
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z");

          // Z-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");

          // G-X
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qe.push_back(qtmp);
      
          ndstr.push_back("X");

          // X-A1
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
          qe.push_back(qtmp);
      
          ndstr.push_back("A_1");

          // A1-Y
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y/X");

          // X-A
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
          qe.push_back(qtmp);
      
          ndstr.push_back("A");

          // A-Z
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z/L");

          // L-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");

        } else if (flag > 0.){             // ORCF1

          const double xi  = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
          const double eta = 0.25 * (1. + la*la/(lb*lb) + la*la/(lc*lc));

          ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y");

          // Y-T
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
      
          ndstr.push_back("T");

          // T-Z
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z");

          // Z-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");

          // G-X
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qe.push_back(qtmp);
      
          ndstr.push_back("X");

          // X-A1
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
          qe.push_back(qtmp);
      
          ndstr.push_back("A_1");

          // A1-Y
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5-xi; qtmp[2] = 1.-xi;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y/T");

          // T-X1
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 1.; qtmp[1] = qtmp[2] = 1.-eta;
          qe.push_back(qtmp);
      
          ndstr.push_back("X_1/X");

          // X-A
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
          qe.push_back(qtmp);
      
          ndstr.push_back("A");

          // A-Z
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5+xi; qtmp[2] = xi;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z/L");

          // L-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");

        } else {                           // ORCF2

          const double eta = 0.25 * (1. + la*la/(lb*lb) - la*la/(lc*lc));
          const double phi = 0.25 * (1. + lc*lc/(lb*lb) - lc*lc/(la*la));
          const double del = 0.25 * (1. + lb*lb/(la*la) - lb*lb/(lc*lc));

          ndstr.push_back("{/Symbol G}");
          // G-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y");

          // Y-C
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[2] = 0.5-eta; qtmp[1] = 1.-eta;
          qe.push_back(qtmp);
      
          ndstr.push_back("C");

          // C-D
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[2] = 0.5-eta; qtmp[1] = 1.-eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5-del; qtmp[2] = 0.5; qtmp[1] = 1.-del;
          qe.push_back(qtmp);
      
          ndstr.push_back("D");

          // D-X
          qtmp = new double [3];
          qtmp[0] = 0.5-del; qtmp[2] = 0.5; qtmp[1] = 1.-del;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
      
          ndstr.push_back("X");

          // X-G
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");

          // G-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z");

          // Z-D1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5+del; qtmp[1] = 0.5; qtmp[2] = del;
          qe.push_back(qtmp);
      
          ndstr.push_back("D_1");

          // D1-H
          qtmp = new double [3];
          qtmp[0] = 0.5+del; qtmp[1] = 0.5; qtmp[2] = del;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
          qe.push_back(qtmp);
      
          ndstr.push_back("H");

          // H-C
          qtmp = new double [3];
          qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5-eta; qtmp[2] = 1.-eta;
          qe.push_back(qtmp);
      
          ndstr.push_back("C/C_1");

          // C1-Z
          qtmp = new double [3];
          qtmp[0] = 0.5; qtmp[1] = 0.5+eta; qtmp[2] = eta;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Z/X");

          // X-H1
          qtmp = new double [3];
          qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = phi; qtmp[1] = 0.5+phi; qtmp[2] = 0.5;
          qe.push_back(qtmp);
      
          ndstr.push_back("H_1/H");

          // H-Y
          qtmp = new double [3];
          qtmp[0] = 1.-phi; qtmp[1] = 0.5-phi; qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("Y/L");

          // L-G
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
      
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
      
          ndstr.push_back("{/Symbol G}");
        }

      } else if (symbol[0] == 'C'){        // ORCC

        const double xi = 0.25 * (1. + la*la/(lb*lb));

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("X");

        // X-S
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("S");

        // S-R
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("R");

        // R-A
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("A");

        // A-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = xi; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("Z");

        // Z-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("{/Symbol G}");

        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("Y");

        // Y-X1
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("X_1");

        // X1-A1
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("A_1");

        // A1-T
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = 1.-xi; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("T");

        // T-Y
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("Y/Z");

        // Z-T
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("T");

      } else {                             // ORCI

        const double xi  = 0.25 * (1. + la*la/(lc*lc));
        const double eta = 0.25 * (1. + lb*lb/(lc*lc));
        const double del = (lb*lb-la*la)/(4.*lc*lc);
        const double miu = (lb*lb+la*la)/(4.*lc*lc);

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = qtmp[2] = xi;
        qe.push_back(qtmp);
    
        ndstr.push_back("X");

        // X-L
        qtmp = new double [3];
        qtmp[0] = -xi; qtmp[1] = qtmp[2] = xi;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = -miu; qtmp[1] = miu; qtmp[2] = 0.5-del;
        qe.push_back(qtmp);
    
        ndstr.push_back("L");

        // L-T
        qtmp = new double [3];
        qtmp[0] = -miu; qtmp[1] = miu; qtmp[2] = 0.5-del;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("T");

        // T-W
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qe.push_back(qtmp);
    
        ndstr.push_back("W");

        // W-R
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("R");

        // R-X1
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = xi; qtmp[1] = 1.-xi; qtmp[2] = -xi;
        qe.push_back(qtmp);
    
        ndstr.push_back("X_1");

        // X1-Z
        qtmp = new double [3];
        qtmp[0] = xi; qtmp[1] = 1.-xi; qtmp[2] = -xi;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("Z");

        // Z-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("{/Symbol G}");

        // G-Y
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
        qe.push_back(qtmp);
    
        ndstr.push_back("Y");

        // Y-S
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
    
        ndstr.push_back("S");

        // S-W
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qe.push_back(qtmp);
    
        ndstr.push_back("W/L_1");

        // L1-Y
        qtmp = new double [3];
        qtmp[0] = miu; qtmp[1] = -miu; qtmp[2] = 0.5+del;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = eta; qtmp[1] = -eta;
        qe.push_back(qtmp);
    
        ndstr.push_back("Y/Y_1");

        // Y1-Z
        qtmp = new double [3];
        qtmp[0] = 1.-eta; qtmp[1] = eta; qtmp[2] = -eta;
        qs.push_back(qtmp);
    
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
        qe.push_back(qtmp);
    
        ndstr.push_back("Z");
      }

    } else if (spgnum >= 75 && spgnum <= 142){ // Tetragonal

      if (symbol[0] == 'P'){             // TET

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("X");

        // X-M
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("M");

        // M-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);

        ndstr.push_back("{/Symbol G}");

        // G-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Z");

        // Z-R
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("R");

        // R-A
        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("A");

        // A-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("Z/X");

        // X-R
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = 0.; qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("R/M");

        // M-A
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);

        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);

        ndstr.push_back("A");

      } else {                           // BCT
        if (la > lc){                    // BCT1
          const double eta = 0.25 * (1. + lc*lc/(la*la));

          ndstr.push_back("{/Symbol G}");
          // G-X
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("X");

          // X-M
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("M");

          // M-G
          qtmp = new double [3];
          qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
   
          ndstr.push_back("{/Symbol G}");

          // G-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = eta; qtmp[2] = -eta;
          qe.push_back(qtmp);
   
          ndstr.push_back("Z");

          // Z-P
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = eta; qtmp[2] = -eta;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
          qe.push_back(qtmp);
   
          ndstr.push_back("P");

          // P-N
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("N");

          // N-Z1
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = -eta; qtmp[1] = 1.-eta; qtmp[2] = eta;
          qe.push_back(qtmp);
   
          ndstr.push_back("Z_1");

          // Z1-M
          qtmp = new double [3];
          qtmp[0] = -eta; qtmp[1] = 1.-eta; qtmp[2] = eta;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = -0.5; qtmp[1] = qtmp[2] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("M");

          // X-P
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
          qe.push_back(qtmp);
   
          ndstr.push_back("P");

        } else {                         // BCT2

          const double eta = 0.25 * (1. + la*la/(lc*lc));
          const double xi  = la*la/(2.*lc*lc);

          ndstr.push_back("{/Symbol G}");
          // G-X
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("X");

          // X-Y
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = -xi; qtmp[1] = xi; qtmp[2] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("Y");

          // Y-Sigma
          qtmp = new double [3];
          qtmp[0] = -xi; qtmp[1] = xi; qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = -eta; qtmp[1] = qtmp[2] = eta;
          qe.push_back(qtmp);
   
          ndstr.push_back("{/Symbol S}");

          // Sigma-G
          qtmp = new double [3];
          qtmp[0] = -eta; qtmp[1] = qtmp[2] = eta;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qe.push_back(qtmp);
   
          ndstr.push_back("{/Symbol G}");

          // G-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = qtmp[2] = 0.;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("Z");

          // Z-Sigma_1
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.5; qtmp[2] = -0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = eta; qtmp[1] = 1.-eta; qtmp[2] = -eta;
          qe.push_back(qtmp);
   
          ndstr.push_back("{/Symbol S}_1");

          // Sigma_1-N
          qtmp = new double [3];
          qtmp[0] = eta; qtmp[1] = 1.-eta; qtmp[2] = -eta;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("N");

          // N-P
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
          qe.push_back(qtmp);
   
          ndstr.push_back("P");

          // P-Y1
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -xi;
          qe.push_back(qtmp);
   
          ndstr.push_back("Y_1");

          // Y1-Z
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -xi;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
          qe.push_back(qtmp);
   
          ndstr.push_back("Z");

          // X-P
          qtmp = new double [3];
          qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
          qs.push_back(qtmp);
   
          qtmp = new double [3];
          qtmp[0] = qtmp[2] = qtmp[1] = 0.25;
          qe.push_back(qtmp);
   
          ndstr.push_back("P");
        }
      }
    } else if (spgnum >= 143 && spgnum <= 167){ // Trigonal
      if (cosg > 0.){                           // RHL1

        const double eta = (1.+4.*cosa)/(2.+4.*cosa);
        const double niu = 0.75 - 0.5*eta;

        ndstr.push_back("{/Symbol G}");
        // G-L
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("L");

        // L-B1
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 1.-eta; qtmp[2] = eta - 1.;
        qe.push_back(qtmp);
 
        ndstr.push_back("B_1/B");

        // B-Z
        qtmp = new double [3];
        qtmp[0] = eta; qtmp[1] = 0.5; qtmp[2] = 1.-eta;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
 
        ndstr.push_back("Z");

        // Z-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("{/Symbol G}");

        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = niu; qtmp[1] = 0.; qtmp[2] = -niu;
        qe.push_back(qtmp);
 
        ndstr.push_back("X/Q");

        // Q-F
        qtmp = new double [3];
        qtmp[0] = 1.-niu; qtmp[1] = niu; qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("F");

        // F-P1
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-eta;
        qe.push_back(qtmp);
 
        ndstr.push_back("P_1");

        // P1-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 1.-niu; qtmp[2] = 1.-eta;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
 
        ndstr.push_back("Z/L");

        // L-P
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = eta; qtmp[1] = qtmp[2] = niu;
        qe.push_back(qtmp);
 
        ndstr.push_back("P");

      } else {                                  // RHL2
        const double eta = 0.5*(1.+cosa)/(1.-cosa);
        const double niu = 0.75 - 0.5*eta;

        ndstr.push_back("{/Symbol G}");
        // G-P
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 1.-niu; qtmp[1] = -niu;
        qe.push_back(qtmp);
 
        ndstr.push_back("P");

        // P-Z
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 1.-niu; qtmp[1] = -niu;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qe.push_back(qtmp);
 
        ndstr.push_back("Z");

        // Z-Q
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = eta;
        qe.push_back(qtmp);
 
        ndstr.push_back("Q");

        // Q-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = eta;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("{/Symbol G}");

        // G-F
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("F");

        // F-P1
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = -0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = niu; qtmp[1] = qtmp[2] = niu - 1.;
        qe.push_back(qtmp);
 
        ndstr.push_back("P_1");

        // P1-Q1
        qtmp = new double [3];
        qtmp[0] = niu; qtmp[1] = qtmp[2] = niu - 1.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = 1.-eta; qtmp[1] = qtmp[2] = -eta;
        qe.push_back(qtmp);
 
        ndstr.push_back("Q_1");

        // Q1-L
        qtmp = new double [3];
        qtmp[0] = 1.-eta; qtmp[1] = qtmp[2] = -eta;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
 
        ndstr.push_back("L");

        // L-Z
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
 
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qe.push_back(qtmp);
 
        ndstr.push_back("Z");
      }

    } else if (spgnum >= 168 && spgnum <= 194){ // Hexagonal

      ndstr.push_back("{/Symbol G}");
      // G-M
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      qe.push_back(qtmp);
 
      ndstr.push_back("M");
     
      // M-K
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
      qe.push_back(qtmp);
 
      ndstr.push_back("K");
     
      // K-G
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      qe.push_back(qtmp);
 
      ndstr.push_back("{/Symbol G}");
     
      // G-A
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = qtmp[2] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      qe.push_back(qtmp);
 
      ndstr.push_back("A");
     
      // A-L
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      qe.push_back(qtmp);
 
      ndstr.push_back("L");
     
      // L-H
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
      qe.push_back(qtmp);
 
      ndstr.push_back("H");
     
      // H-A
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
      qe.push_back(qtmp);
 
      ndstr.push_back("A/L");
     
      // L-M
      qtmp = new double [3];
      qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = 0.5; qtmp[1] = qtmp[2] = 0.;
      qe.push_back(qtmp);
 
      ndstr.push_back("M/K");
     
      // K-H
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.;
      qs.push_back(qtmp);
 
      qtmp = new double [3];
      qtmp[0] = qtmp[1] = 1./3.; qtmp[2] = 0.5;
      qe.push_back(qtmp);
 
      ndstr.push_back("H");
     
    } else if (spgnum >= 195 && spgnum <= 230){ // Cubic

      if (symbol[0] == 'P'){                    // CUB

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("X");
     
        // X-M
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qe.push_back(qtmp);
   
        ndstr.push_back("M");
     
        // M-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
   
        ndstr.push_back("{/Symbol G}");
     
        // G-R
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("R");
     
        // R-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.; qtmp[1] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("X/M");
     
        // M-R
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.5; qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("R");

      } else if (symbol[0] == 'F'){          // FCC

        ndstr.push_back("{/Symbol G}");
        // G-X
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qe.push_back(qtmp);
   
        ndstr.push_back("X");
     
        // X-W
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
        qe.push_back(qtmp);
   
        ndstr.push_back("W");
     
        // W-K
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
        qe.push_back(qtmp);
   
        ndstr.push_back("K");
     
        // K-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
   
        ndstr.push_back("{/Symbol G}");
     
        // G-L
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("L");
     
        // L-U
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
        qe.push_back(qtmp);
   
        ndstr.push_back("U");
     
        // U-W
        qtmp = new double [3];
        qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
        qe.push_back(qtmp);
   
        ndstr.push_back("W");
     
        // W-L
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.25; qtmp[2] = 0.75;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("L");
     
        // L-K
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.375; qtmp[2] = 0.75;
        qe.push_back(qtmp);
   
        ndstr.push_back("K/U");
     
        // U-X
        qtmp = new double [3];
        qtmp[0] = 0.625; qtmp[1] = 0.25; qtmp[2] = 0.625;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = 0.5; qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("X");
     
      } else {                               // BCC

        ndstr.push_back("{/Symbol G}");
        // G-H
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("H");
     
        // H-N
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("N");
     
        // N-G
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qe.push_back(qtmp);
   
        ndstr.push_back("{/Symbol G}");
     
        // G-P
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qe.push_back(qtmp);
   
        ndstr.push_back("P");
     
        // P-H
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[2] = 0.5; qtmp[1] = -0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("H/P");
     
        // P-N
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = qtmp[2] = 0.25;
        qs.push_back(qtmp);
   
        qtmp = new double [3];
        qtmp[0] = qtmp[1] = 0.; qtmp[2] = 0.5;
        qe.push_back(qtmp);
   
        ndstr.push_back("N");
     
      }

    } else {
      printf("\nSorry, failed to identify the crystal system, please use the manual mode.\n");
    }

    // to determine the number of points along each line, with a step size of 0.05
    const double qs_inv = 1./QSTEP;
    int nbin = qs.size();
    if (nbin > 0) printf("\nPhonon dispersion will be evaluated along lines:\n\t%s", ndstr[0].c_str());
    for (int is = 0; is < nbin; ++is){
      double *qstr = qs[is];
      double *qend = qe[is];
      double ql = 0.;
      for (int i = 0; i < 3; ++i) ql += (qend[i] - qstr[i])*(qend[i] - qstr[i]);
      int nqpt = MAX(int(sqrt(ql) * qs_inv + 0.5), 2);
      nqbin.push_back(nqpt);

      printf("-%s", ndstr[is+1].c_str());
    }
    if (nbin > 0) printf("\n");
  }
#endif

  FILE *fp = fopen(fname, "w");
  fprintf(fp,"# q     qr    freq\n");
  fprintf(fp,"# 2pi/L  2pi/L %s\n", dynmat->funit);

  double qr = 0., dq, q[3], qinc[3];
  int nbin = qs.size();
  for (int is = 0; is < nbin; ++is){
    double *qstr = qs[is];
    double *qend = qe[is];
    int nbin = nqbin[is];
    for (int i = 0; i < 3; ++i) qinc[i] = (qend[i]-qstr[i])/double(nbin-1);
    dq = sqrt(qinc[0]*qinc[0]+qinc[1]*qinc[1]+qinc[2]*qinc[2]);
  
    nodes.push_back(qr);
    for (int i = 0; i < 3; ++i) q[i] = qstr[i];
    for (int ii = 0; ii < nbin; ++ii){
      double wii = 1.;
      dynmat->getDMq(q, &wii);
      if (wii > 0.){
        dynmat->geteigen(egvs, 0);
        fprintf(fp,"%lg %lg %lg %lg ", q[0], q[1], q[2], qr);
        for (int i = 0; i < ndim; ++i) fprintf(fp," %lg", egvs[i]);
      }
      fprintf(fp,"\n");
  
      for (int i = 0; i < 3; ++i) q[i] += qinc[i];
      qr += dq;
    }
    qr -= dq;
    delete []qstr;
    delete []qend;
  }
  qs.clear(); qe.clear();
  if (qr > 0.) nodes.push_back(qr);
  fclose(fp);
  delete []egvs;

  // write the gnuplot script which helps to visualize the result
  int nnd = nodes.size();
  if (nnd > 1){
    const char qmk = char(34); // "
    fp = fopen("pdisp.gnuplot", "w");
    fprintf(fp,"set term post enha colo 20\nset out %cpdisp.eps%c\n\n", qmk, qmk);
    fprintf(fp,"set xlabel %cq%c\n", qmk, qmk);
    fprintf(fp,"set ylabel %cfrequency (THz)%c\n\n", qmk, qmk);
    fprintf(fp,"set xrange [0:%lg]\nset yrange [0:*]\n\n", nodes[nnd-1]);
    fprintf(fp,"set grid xtics\n");
    fprintf(fp,"# {/Symbol G} will give you letter gamma in the label\nset xtics (");
    for (int i = 0; i < nnd-1; ++i) fprintf(fp,"%c%s%c %lg, ", qmk, ndstr[i].c_str(), qmk, nodes[i]);
    fprintf(fp, "%c%s%c %lg)\n\n", qmk, ndstr[nnd-1].c_str(), qmk, nodes[nnd-1]);
    fprintf(fp, "unset key\n\n");
    fprintf(fp, "plot %c%s%c u 4:5 w l lt 1", qmk, fname, qmk);
    for (int i = 1; i < ndim; ++i) fprintf(fp,",\\\n%c%c u 4:%d w l lt 1", qmk, qmk, i+5);
    fclose(fp);

    printf("\nPhonon dispersion data are written to: %s, you can visualize the results\n", fname);
    printf("by invoking: `gnuplot pdisp.gnuplot; gv pdisp.eps`\n");
  }
  for (int ii = 0; ii < 80; ++ii) printf("=");
  printf("\n");

  delete []fname;
  nodes.clear();
  ndstr.clear();

return;
}
