/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   see LLNL copyright notice at bottom of file
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(meam/spline,PairMEAMSpline)

#else

#ifndef LMP_PAIR_MEAM_SPLINE_H
#define LMP_PAIR_MEAM_SPLINE_H

#include "pair.h"

namespace LAMMPS_NS {

/// Set this to 1 if you intend to use MEAM potentials with non-uniform spline knots.
/// Set this to 0 if you intend to use only MEAM potentials with spline knots on a uniform grid.
///
/// With SUPPORT_NON_GRID_SPLINES == 0, the code runs about 50% faster.

#define SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES 0

class PairMEAMSpline : public Pair
{
public:
        PairMEAMSpline(class LAMMPS *);
        virtual ~PairMEAMSpline();
        virtual void compute(int, int);
        void settings(int, char **);
        void coeff(int, char **);
        void init_style();
        void init_list(int, class NeighList *);
        double init_one(int, int);

        int pack_forward_comm(int, int *, double *, int, int *);
        void unpack_forward_comm(int, int, double *);
        int pack_reverse_comm(int, int, double *);
        void unpack_reverse_comm(int, int *, double *);
        double memory_usage();

protected:
  char **elements;              // names of unique elements
  int *map;                     // mapping from atom types to elements
  int nelements;                // # of unique elements

        class SplineFunction {
        public:

                /// Default constructor.
                SplineFunction() : X(NULL), Xs(NULL), Y(NULL), Y2(NULL), Ydelta(NULL), N(0) {}

                /// Destructor.
                ~SplineFunction() {
                        delete[] X;
                        delete[] Xs;
                        delete[] Y;
                        delete[] Y2;
                        delete[] Ydelta;
                }

                /// Initialization of spline function.
                void init(int _N, double _deriv0, double _derivN) {
                        N = _N;
                        deriv0 = _deriv0;
                        derivN = _derivN;
                        delete[] X;
                        delete[] Xs;
                        delete[] Y;
                        delete[] Y2;
                        delete[] Ydelta;
                        X = new double[N];
                        Xs = new double[N];
                        Y = new double[N];
                        Y2 = new double[N];
                        Ydelta = new double[N];
                }

                /// Adds a knot to the spline.
                void setKnot(int n, double x, double y) { X[n] = x; Y[n] = y; }

                /// Returns the number of knots.
                int numKnots() const { return N; }

                /// Parses the spline knots from a text file.
                void parse(FILE* fp, Error* error);

                /// Calculates the second derivatives of the cubic spline.
                void prepareSpline(Error* error);

                /// Evaluates the spline function at position x.
                inline double eval(double x) const
                {
                        x -= xmin;
                        if(x <= 0.0) {  // Left extrapolation.
                                return Y[0] + deriv0 * x;
                        }
                        else if(x >= xmax_shifted) {  // Right extrapolation.
                                return Y[N-1] + derivN * (x - xmax_shifted);
                        }
                        else {
#if SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES
                                // Do interval search.
                                int klo = 0;
                                int khi = N-1;
                                while(khi - klo > 1) {
                                        int k = (khi + klo) / 2;
                                        if(Xs[k] > x) khi = k;
                                        else klo = k;
                                }
                                double h = Xs[khi] - Xs[klo];
                                // Do spline interpolation.
                                double a = (Xs[khi] - x)/h;
                                double b = 1.0 - a; // = (x - X[klo])/h
                                return a * Y[klo] + b * Y[khi] + ((a*a*a - a) * Y2[klo] + (b*b*b - b) * Y2[khi])*(h*h)/6.0;
#else
                                // For a spline with grid points, we can directly calculate the interval X is in.
                                int klo = (int)(x / h);
                                int khi = klo + 1;
                                double a = Xs[khi] - x;
                                double b = h - a;
                                return Y[khi] - a * Ydelta[klo] + ((a*a - hsq) * a * Y2[klo] + (b*b - hsq) * b * Y2[khi]);
#endif
                        }
                }

                /// Evaluates the spline function and its first derivative at position x.
                inline double eval(double x, double& deriv) const
                {
                        x -= xmin;
                        if(x <= 0.0) {  // Left extrapolation.
                                deriv = deriv0;
                                return Y[0] + deriv0 * x;
                        }
                        else if(x >= xmax_shifted) {  // Right extrapolation.
                                deriv = derivN;
                                return Y[N-1] + derivN * (x - xmax_shifted);
                        }
                        else {
#if SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES
                                // Do interval search.
                                int klo = 0;
                                int khi = N-1;
                                while(khi - klo > 1) {
                                        int k = (khi + klo) / 2;
                                        if(Xs[k] > x) khi = k;
                                        else klo = k;
                                }
                                double h = Xs[khi] - Xs[klo];
                                // Do spline interpolation.
                                double a = (Xs[khi] - x)/h;
                                double b = 1.0 - a; // = (x - X[klo])/h
                                deriv = (Y[khi] - Y[klo]) / h + ((3.0*b*b - 1.0) * Y2[khi] - (3.0*a*a - 1.0) * Y2[klo]) * h / 6.0;
                                return a * Y[klo] + b * Y[khi] + ((a*a*a - a) * Y2[klo] + (b*b*b - b) * Y2[khi]) * (h*h) / 6.0;
#else
                                // For a spline with grid points, we can directly calculate the interval X is in.
                                int klo = (int)(x / h);
                                int khi = klo + 1;
                                double a = Xs[khi] - x;
                                double b = h - a;
                                deriv = Ydelta[klo] + ((3.0*b*b - hsq) * Y2[khi] - (3.0*a*a - hsq) * Y2[klo]);
                                return Y[khi] - a * Ydelta[klo] + ((a*a - hsq) * a * Y2[klo] + (b*b - hsq) * b * Y2[khi]);
#endif
                        }
                }

                /// Returns the number of bytes used by this function object.
                double memory_usage() const { return sizeof(*this) + sizeof(X[0]) * N * 3; }

                /// Returns the cutoff radius of this function.
                double cutoff() const { return X[N-1]; }

                /// Writes a Gnuplot script that plots the spline function.
                void writeGnuplot(const char* filename, const char* title = NULL) const;

                /// Broadcasts the spline function parameters to all processors.
                void communicate(MPI_Comm& world, int me);

        private:
                double* X;                                // Positions of spline knots
                double* Xs;                                // Shifted positions of spline knots
                double* Y;                                // Function values at spline knots
                double* Y2;                                // Second derivatives at spline knots
                double* Ydelta;                        // If this is a grid spline, Ydelta[i] = (Y[i+1]-Y[i])/h
                int N;                                        // Number of spline knots
                double deriv0;                        // First derivative at knot 0
                double derivN;                        // First derivative at knot (N-1)
                double xmin;                        // The beginning of the interval on which the spline function is defined.
                double xmax;                        // The end of the interval on which the spline function is defined.
                int isGridSpline;                // Indicates that all spline knots are on a regular grid.
                double h;                                // The distance between knots if this is a grid spline with equidistant knots.
                double hsq;                                // The squared distance between knots if this is a grid spline with equidistant knots.
                double xmax_shifted;        // The end of the spline interval after it has been shifted to begin at X=0.
        };

        /// Helper data structure for potential routine.
        struct MEAM2Body {
                int tag;
                double r;
                double f, fprime;
                double del[3];
        };

        SplineFunction phi;                        // Phi(r_ij)
        SplineFunction rho;                        // Rho(r_ij)
        SplineFunction f;                        // f(r_ij)
        SplineFunction U;                        // U(rho)
        SplineFunction g;                        // g(cos_theta)
        double zero_atom_energy;        // Shift embedding energy by this value to make it zero for a single atom in vacuum.

        double cutoff;              // The cutoff radius

        double* Uprime_values;                // Used for temporary storage of U'(rho) values
        int nmax;                                        // Size of temporary array.
        int maxNeighbors;                        // The last maximum number of neighbors a single atoms has.
        MEAM2Body* twoBodyInfo;                // Temporary array.

        void read_file(const char* filename);
        void allocate();
};

}

#endif
#endif

/* ----------------------------------------------------------------------
 * Spline-based Modified Embedded Atom method (MEAM) potential routine.
 *
 * Copyright (2011) Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Alexander Stukowski (<alex@stukowski.com>).
 * LLNL-CODE-525797 All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2, dated June 1991.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
 * GNU General Public License for more details.
 *
 * Our Preamble Notice
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the
 * Lawrence Livermore National Laboratory under Contract No.
 * DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process,
 * or services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * See file 'pair_spline_meam.cpp' for history of changes.
------------------------------------------------------------------------- */
