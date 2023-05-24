
#include "phonon.h"

#include "dynmat.h"
#include "global.h"
#include "input.h"
#include "kpath.h"
#include "qnodes.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

/*------------------------------------------------------------------------------
 * Private method to evaluate the phonon dispersion curves
 *----------------------------------------------------------------------------*/
void Phonon::pdisp()
{
  // ask the output file name and write the header.
   char str[MAXLINE];
   puts("================================================================================");

#ifdef UseSPG
  // ask method to generate q-lines
   int method = 2;
   printf("Please select your method to generate the phonon dispersion:\n");
   printf("  1. Manual, should always work;\n");
   printf("  2. Automatic, works only for 3D crystals (CMS49-299).\nYour choice [2]: ");
   input->read_stdin(str);
   if (count_words(str) > 0) method = atoi(strtok(str," \t\n\r\f"));
   method = 2 - method%2;
   printf("Your  selection: %d\n", method);
#endif
   printf("\nPlease input the filename to output the dispersion data [pdisp.dat]:");
   input->read_stdin(str);
   if (count_words(str) < 1) strcpy(str, "pdisp.dat");
   char *ptr = strtok(str," \t\n\r\f");
   char *fname = new char[strlen(ptr)+1];
   strcpy(fname,ptr);
 
   // to store the nodes of the dispersion curve
   QNodes *qnodes = new QNodes();
 
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
         input->read_stdin(str);
         int n = count_words(str);
         ptr = strtok(str, " \t\n\r\f");
         if ((n == 1) && (strcmp(ptr,"q") == 0)) break;
         else if (n >= 3){
           qstr[0] = atof(ptr);
           qstr[1] = atof(strtok(NULL, " \t\n\r\f"));
           qstr[2] = atof(strtok(NULL, " \t\n\r\f"));
         }
     
         while ( 1 ){
            printf("Please input the end q-point in unit of B1->B3: ");
            input->read_stdin(str);
            if (count_words(str) >= 3) break;
         }
         qend[0] = atof(strtok(str,  " \t\n\r\f"));
         qend[1] = atof(strtok(NULL, " \t\n\r\f"));
         qend[2] = atof(strtok(NULL, " \t\n\r\f"));
     
         printf("Please input the # of points along the line [%d]: ", nq);
         input->read_stdin(str);
         if (count_words(str) > 0) nq = atoi(strtok(str," \t\n\r\f"));
         nq = MAX(nq,2);
     
         double *qtmp = new double [3];
         for (int i = 0; i < 3; ++i) qtmp[i] = qstr[i];
         qnodes->qs.push_back(qtmp);
   
         qtmp = new double [3];
         for (int i = 0; i < 3; ++i) qtmp[i] = qend[i];
         qnodes->qe.push_back(qtmp);
   
         qnodes->nqbin.push_back(nq);
         qnodes->ndstr.push_back("");
      }
      qnodes->ndstr.push_back("");
   
#ifdef UseSPG
   } else {
      kPath *kp = new kPath(dynmat, qnodes);
      kp->show_info();
      kp->kpath();
      kp->show_path();
      delete kp;
   }
#endif

   FILE *fp = fopen(fname, "w");
   fprintf(fp,"# q     qr    freq\n");
   fprintf(fp,"# 2pi/L  2pi/L %s\n", dynmat->funit);
 
   double qr = 0., dq, q[3], qinc[3];
   int nbin = qnodes->qs.size();
   qnodes->nodes.clear();
   for (int is = 0; is < nbin; ++is){
     double *qstr = qnodes->qs[is];
     double *qend = qnodes->qe[is];
     int nbin = qnodes->nqbin[is];
     for (int i = 0; i < 3; ++i) qinc[i] = (qend[i]-qstr[i])/double(nbin-1);
     dq = sqrt(qinc[0]*qinc[0]+qinc[1]*qinc[1]+qinc[2]*qinc[2]);
   
     qnodes->nodes.push_back(qr);
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
   if (qr > 0.) qnodes->nodes.push_back(qr);
   fclose(fp);
   delete []egvs;
 
   // write the gnuplot script which helps to visualize the result
   int nnd = qnodes->nodes.size();
   if (nnd > 1){
     const char qmk = char(34); // "
     fp = fopen("pdisp.gnuplot", "w");
     fprintf(fp,"set term post enha colo 20\nset out %cpdisp.eps%c\n\n", qmk, qmk);
     fprintf(fp,"set xlabel %cq%c\n", qmk, qmk);
     fprintf(fp,"set ylabel %cfrequency (THz)%c\n\n", qmk, qmk);
     fprintf(fp,"set xrange [0:%lg]\nset yrange [0:*]\n\n", qnodes->nodes[nnd-1]);
     fprintf(fp,"set grid xtics\n");
     fprintf(fp,"# {/Symbol G} will give you letter gamma in the label\nset xtics (");
     for (int i = 0; i < nnd-1; ++i) fprintf(fp,"%c%s%c %lg, ", qmk, qnodes->ndstr[i].c_str(), qmk, qnodes->nodes[i]);
     fprintf(fp, "%c%s%c %lg)\n\n", qmk, qnodes->ndstr[nnd-1].c_str(), qmk, qnodes->nodes[nnd-1]);
     fprintf(fp, "unset key\n\n");
     fprintf(fp, "plot %c%s%c u 4:5 w l lt 1", qmk, fname, qmk);
     for (int i = 1; i < ndim; ++i) fprintf(fp,",\\\n%c%c u 4:%d w l lt 1", qmk, qmk, i+5);
     fclose(fp);
 
     printf("\nPhonon dispersion data are written to: %s, you can visualize the results\n", fname);
     printf("by invoking: `gnuplot pdisp.gnuplot; gv pdisp.eps`\n");
   }
   puts("================================================================================");
 
   delete []fname;
   delete qnodes;
 
return;
}

/*----------------------------------------------------------------------------*/
