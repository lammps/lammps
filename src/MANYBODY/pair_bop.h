/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   The this work follows the formulation from (a) D.G. Pettifor, et al., Mat.
   Sci. and Eng. A365, 2-13, (2004) and (b) D.A. Murdick, et al., Phys.
   Rev. B 73, 045206 (2006). (c) D.K. Ward, et al., Phys. Rev. B 85, 115206
   (2012)

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bop,PairBOP)

#else

#ifndef LMP_PAIR_BOP_H
#define LMP_PAIR_BOP_H

#include "pair.h"
#include <cstring>

namespace LAMMPS_NS {

  class PairBOP : public Pair {

    public:
      PairBOP(class LAMMPS *);
      virtual ~PairBOP();
      void compute(int, int);
      void settings(int, char **);
      void coeff(int, char **);
      void init_style();
      double init_one(int, int);
      double memory_usage();

    private:

      class tabularFunction {
 
        public:
 
          tabularFunction() {
            size = 0;
            xmin = 0.0;
            xmax = 0.0;
            xmaxsq = 0.0;
            xs = NULL;
            ys = NULL;
            ys1 = NULL;
            ys2 = NULL;
            ys3 = NULL;
            ys4 = NULL;
            ys5 = NULL;
            ys6 = NULL;
          }
          tabularFunction(int n) {
            size = n;
            xmin = 0.0;
            xmax = 0.0;
            xmaxsq = 0.0;
            if (n == 0) {
              xs = NULL;
              ys = NULL;
              ys1 = NULL;
              ys2 = NULL;
              ys3 = NULL;
              ys4 = NULL;
              ys5 = NULL;
              ys6 = NULL;
            } else {
              xs = new double[n];
              ys = new double[n];
              ys1 = new double[n];
              ys2 = new double[n];
              ys3 = new double[n];
              ys4 = new double[n];
              ys5 = new double[n];
              ys6 = new double[n];
            }
          }
          tabularFunction(int n, double x1, double x2) {
            size = n;
            xmin = x1;
            xmax = x2;
            xmaxsq = xmax*xmax;
            if (n == 0) {
              xs = NULL;
              ys = NULL;
              ys1 = NULL;
              ys2 = NULL;
              ys3 = NULL;
              ys4 = NULL;
              ys5 = NULL;
              ys6 = NULL;
            } else {
              xs = new double[n];
              ys = new double[n];
              ys1 = new double[n];
              ys2 = new double[n];
              ys3 = new double[n];
              ys4 = new double[n];
              ys5 = new double[n];
              ys6 = new double[n];
            }
          }

          virtual ~tabularFunction() {
            if (xs) delete [] xs;
            if (ys) delete [] ys;
            if (ys1) delete [] ys1;
            if (ys2) delete [] ys2;
            if (ys3) delete [] ys3;
            if (ys4) delete [] ys4;
            if (ys5) delete [] ys5;
            if (ys6) delete [] ys6;
          }

          void set_xrange(double x1, double x2) {
            xmin = x1;
            xmax = x2;
            xmaxsq = xmax*xmax;
          }

          void set_values(int n, double x1, double x2, double * values)
          {
            reset_size(n);
            xmin = x1;
            xmax = x2;
            xmaxsq = xmax*xmax;
            memcpy(ys,values,n*sizeof(double));
            initialize();
          }

          double get_xmin() {
            return xmin;
          }

          double get_xmax() {
            return xmax;
          }

          double get_xmaxsq() {
            return xmaxsq;
          }

          void value(double x, double &y, int ny, double &y1, int ny1)
          {
            double ps = (x - xmin) * rdx + 1.0;
            int ks = ps;
            if (ks > size-1) ks = size-1;
            ps = ps - ks;
            if (ps > 1.0) ps = 1.0;
            if (ny) y = ((ys3[ks-1]*ps + ys2[ks-1])*ps + ys1[ks-1])*ps + ys[ks-1];
            if (ny1) y1 = (ys6[ks-1]*ps + ys5[ks-1])*ps + ys4[ks-1];
          }

        protected:

          void reset_size(int n) {
            if (n != size) {
              size = n;
              if (xs) delete [] xs;
              xs = new double[n];
              if (ys) delete [] ys;
              ys = new double[n];
              if (ys1) delete [] ys1;
              ys1 = new double[n];
              if (ys2) delete [] ys2;
              ys2 = new double[n];
              if (ys3) delete [] ys3;
              ys3 = new double[n];
              if (ys4) delete [] ys4;
              ys4 = new double[n];
              if (ys5) delete [] ys5;
              ys5 = new double[n];
              if (ys6) delete [] ys6;
              ys6 = new double[n];
            }
          }

          void initialize() {
            rdx = (xmax - xmin) / (size - 1.0);
            for (int i = 0; i < size; i++) {
              xs[i] = xmin + i * rdx;
            }
            rdx = 1.0 / rdx;
            ys1[0] = ys[1] - ys[0];
            ys1[1] = 0.5 * (ys[2] - ys[0]);
            ys1[size-2] = 0.5 * (ys[size-1] - ys[size-3]);
            ys1[size-1] = ys[size-1] - ys[size-2];
            for (int i = 2; i < size-2; i++) {
              ys1[i]=((ys[i-2]-ys[i+2])+ 8.0*(ys[i+1]-ys[i-1]))/12.0;
            }
            for (int i = 0; i < size-1; i++) {
              ys2[i]=3.0*(ys[i+1]-ys[i])-2.0*ys1[i]-ys1[i+1];
              ys3[i]=ys1[i]+ys1[i+1]-2.0*(ys[i+1]-ys[i]);
            }
            ys2[size-1]=0.0;
            ys3[size-1]=0.0;
            for (int i = 0; i < size; i++) {
              ys4[i]=ys1[i]*rdx;
              ys5[i]=2.0*ys2[i]*rdx;
              ys6[i]=3.0*ys3[i]*rdx;
            }
          }

          int size;
          double xmin,xmax,xmaxsq,rdx;
          double *xs, *ys, *ys1, *ys2, *ys3, *ys4, *ys5, *ys6;

      };

      struct PairParameters {
        double cutB, cutBsq, cutL, cutLsq;
        class tabularFunction * betaS;
        class tabularFunction * betaP;
        class tabularFunction * rep;
        class tabularFunction * cphi;
        class tabularFunction * bo;
        PairParameters() {
          cutB = 0.0;
          cutBsq = 0.0;
          cutL = 0.0;
          cutLsq = 0.0;
          betaS = NULL;
          betaP = NULL;
          rep = NULL;
          cphi = NULL;
          bo = NULL;
        };
      };

      struct TripletParameters {
        class tabularFunction * G;
        TripletParameters() {
          G = NULL;
        };
      };

      struct PairList1 {
        double r,dis[3];
        double betaS, dBetaS, betaP, dBetaP, rep, dRep;
        PairList1() {
        };
      };

      struct PairList2 {
        double r,dis[3];
        double rep, dRep;
        PairList2() {
        };
      };

      struct TripleList {
        double G, dG, cosAng, dCosAngi[3], dCosAngj[3], dCosAngk[3];
        TripleList() {
        };
      };

      struct B_SG{
        double dAA[3];
        double dBB[3];
        double dCC[3];
        double dDD[3];
        double dEE1[3];
        double dFF[3];
        double dAAC[3];
        double dSigB1[3];
        double dSigB[3];
        int temp;
        int i;
        int j;
      };

      struct B_PI{
        double dAA[3];
        double dBB[3];
        double dPiB[3];
        int temp;
        int i;
        int j;
      };

      PairParameters    *pairParameters;
      TripletParameters *tripletParameters;

  // Parameters variables

      double small1, small2, small3g, small4, small5, small6, small7, *pi_p;
      double *sigma_c, *sigma_a, *pi_c, *pi_a, *sigma_delta, *pi_delta;
      double *sigma_f,*sigma_k,*small3;
      double *pro_delta, *pro;

      int me;
      int bop_types;                // number of elments in potential
      int npairs;                   // number of element pairs
      int ntriples;                 // number of all triples
      char **elements;              // names of unique elements
      double bytes;

      int otfly;                    // = 1 faster, more memory, = 0 slower, less memory

      PairList1 *pairlist1;
      PairList2 *pairlist2;
      TripleList *triplelist;

      B_SG *bt_sg;
      B_PI *bt_pi;

      int *BOP_index;               // index for neighbor list position
      int *BOP_total;               // index for neighbor list position
      int *BOP_index2;              // index for neighbor list position
      int *BOP_total2;              // index for neighbor list position
      int *neigh_index;             // index for neighbor list position
      int *neigh_index2;            // index for neighbor list position
      int atomlimit;                // current size of atom based list
      int neighlimit;               // current size of neighbor based list
      int neighlimit2;              // current size of neighbor based list
      int neineilimit;              // current size of triple based list
      int sglimit;                  // current size of bt_sg
      int pilimit;                  // current size of bt_pi
      int *cos_index;               // index for neighbor cosine if not using on the fly
      double cutmax;

      void grab(FILE *, int, double *);
      void write_tables(int);

      void gneigh();
      void angle(double, double *, double, double *, double &,
                 double *, double *);
      double SigmaBo(int, int);
      double PiBo(int, int);
      void read_table(char *);
      void allocate();
      void memory_sg(int);
      void memory_pi(int);
      void initial_sg(int);
      void initial_pi(int);

  };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal pair_style command

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory. Check the input script or data file.

E: Pair style BOP requires atom IDs

This is a requirement to use the BOP potential.

E: Pair style BOP requires newton pair on

This is a restriction to use the BOP potential.

E: Pair style bop requires a comm ghost cutoff of at least %lf

Use the comm_modify cutoff to set this. See the pair bop doc page for
more details.

E: All pair coeffs are not set

Self-explanatory.

E: Cannot open BOP potential file %s

Self-explanatory.

E: Incorrect table format check for element types

Self-explanatory.

E: Unsupported BOP potential file format

Self-explanatory.

E: Pair style bop requires system dimension of at least %g

Self-explanatory.

UNDOCUMENTED

*/
