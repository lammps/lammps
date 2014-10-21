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

PairStyle(meam/sw/spline,PairMEAMSWSpline)

#else

#ifndef LMP_PAIR_MEAM_SW_SPLINE_H
#define LMP_PAIR_MEAM_SW_SPLINE_H

#include "pair.h"

namespace LAMMPS_NS {

/// Set this to 1 if you intend to use MEAM potentials with non-uniform spline knots.
/// Set this to 0 if you intend to use only MEAM potentials with spline knots on a uniform grid.
///
/// With SUPPORT_NON_GRID_SPLINES == 0, the code runs about 50% faster.

#define SPLINE_MEAMSW_SUPPORT_NON_GRID_SPLINES 0

class PairMEAMSWSpline : public Pair
{
public:
        PairMEAMSWSpline(class LAMMPS *);
        virtual ~PairMEAMSWSpline();
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
#if SPLINE_MEAMSW_SUPPORT_NON_GRID_SPLINES
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
                                //
                                int klo = (int)(x / h);
                                if ( klo > N - 2 ) klo = N - 2;
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
#if SPLINE_MEAMSW_SUPPORT_NON_GRID_SPLINES
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
                                if ( klo > N - 2 ) klo = N - 2;
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
                double F, Fprime;
                double del[3];
        };

        SplineFunction phi;                        // Phi(r_ij)
        SplineFunction rho;                        // Rho(r_ij)
        SplineFunction f;                        // f(r_ij)
        SplineFunction U;                        // U(rho)
        SplineFunction g;                        // g(cos_theta)
        SplineFunction F;                        // F(r_ij)
        SplineFunction G;                        // G(cos_theta)
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
 * Spline-based Modified Embedded Atom Method plus 
 * Stillinger-Weber (MEAM+SW) potential routine.
 *
 * Copyright (2012) Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Robert E. Rudd (<robert.rudd@llnl.gov>).
 * Based on the spline MEAM routine written by Alexander Stukowski 
 * (<alex@stukowski.com>).
 * LLNL-CODE-588032 All rights reserved.
 *
 * The spline-based MEAM+SW format was first devised and used to develop
 * potentials for bcc transition metals by Jeremy Nicklas, Michael Fellinger,
 * and Hyoungki Park at The Ohio State University.
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
 * The precise terms and conditions for copying, distribution and modification
 * follows.
 *
 * GNU Terms and Conditions for Copying, Distribution, and Modification
 *
 * 0.  This License applies to any program or other work which contains a
 * notice placed by the copyright holder saying it may be distributed under
 * the terms of this General Public License.  The "Program," below, refers to
 * any such program or work, and a "work based on the Program" means either
 * the Program or any derivative work under copyright law: that is to say, a
 * work containing the Program or a portion of it, either verbatim or with
 * modifications and/or translated into another language.  (Hereinafter,
 * translation is included without limitation in the term "modification".)
 * Each licensee is addressed as "you."
 *
 * Activities other than copying, distribution and modification are not
 * covered by this License; they are outside its scope.  The act of running
 * the Program is not restricted, and the output from the Program is covered
 * only if its contents constitute a work based on the Program (independent of
 * having been made by running the Program).  Whether that is true depends on
 * what the Program does.  
 *
 * 1.  You may copy and distribute verbatim copies of the Program's source
 * code as you receive it, in any medium, provided that you conspicuously and
 * appropriately publish on each copy an appropriate copyright notice and
 * disclaimer of warranty; keep intact all the notices that refer to this
 * License and to the absence of any warranty; and give any other recipients
 * of the Program a copy of this License along with the Program.
 *
 * You may charge a fee for the physical act of transferring a copy, and you
 * may at your option offer warranty protection in exchange for a fee.
 *
 * 2.  You may modify your copy or copies of the Program or any portion of it,
 * thus forming a work based on the Program, and copy and distribute such
 * modifications or work under the terms of Section 1 above, provided that you
 * also meet all of these conditions:
 *
 *  a)  You must cause the modified files to carry prominent notices stating
 *  that you changed the files and the date of any change.
 *
 *  b)  You must cause any work that you distribute or publish, that in whole
 *  or in part contains or is derived from the Program or any part thereof, to
 *  be licensed as a whole at no charge to all third parties under the terms
 *  of this License.
 *
 *  c)  If the modified program normally reads commands interactively when
 *  run, you must cause it, when started running for such interactive use in
 *  the most ordinary way, to print or display an announcement including an
 *  appropriate copyright notice and a notice that there is no warranty (or
 *  else, saying that you provide a warranty) and that users may redistribute
 *  the program under these conditions, and telling the user how to view a
 *  copy of this License.  (Exception: if the Program itself is interactive
 *  but does not normally print such an announcement, your work based on the
 *  Program is not required to print an announcement.)
 *
 * These requirements apply to the modified work as a whole.  If
 * identifiable sections of that work are not derived from the Program, and
 * can be reasonably considered independent and separate works in
 * themselves, then this License, and its terms, do not apply to those
 * sections when you distribute them as separate work.  But when you
 * distribute the same section as part of a whole which is a work based on
 * the Program, the distribution of the whole must be on the terms of this
 * License, whose permissions for other licensees extend to the entire
 * whole, and thus to each and every part regardless of who wrote it.
 *
 * Thus, it is not the intent of this section to claim rights or contest
 * your rights to work written entirely by you; rather, the intent is to
 * exercise the right to control the distribution of derivative or
 * collective works based on the Program.
 *
 * In addition, mere aggregation of another work not based on the Program
 * with the Program (or with a work based on the Program) on a volume of a
 * storage or distribution medium does not bring the other work under the
 * scope of this License.
 *
 * 3.  You may copy and distribute the Program (or a work based on it, under
 * Section 2) in object code or executable form under the terms of Sections
 * 1 and 2 above provided that you also do one of the following:
 *
 *  a)  Accompany it with the complete corresponding machine-readable source
 *  code, which must be distributed under the terms of Sections 1 and 2 above
 *  on a medium customarily used for software interchange; or,
 *
 *  b)  Accompany it with a written offer, valid for at least three years,
 *  to give any third party, for a charge no more than your cost of
 *  physically performing source distribution, a complete machine-readable
 *  copy of the corresponding source code, to be distributed under the terms
 *  of Sections 1 and 2 above on a medium customarily used for software
 *  interchange; or,
 *
 *  c)  Accompany it with the information you received as to the offer to
 *  distribute corresponding source code.  (This alternative is allowed only
 *  for noncommercial distribution and only if you received the program in
 *  object code or executable form with such an offer, in accord with
 *  Subsection b above.)
 *
 * The source code for a work means the preferred form the work for making
 * modifications to it.  For an executable work, complete source code means
 * all the source code for all modules it contains, plus any associated
 * interface definition files, plus the scripts used to control compilation
 * and installation of the executable.  However, as a special exception, the
 * source code distributed need not include anything that is normally
 * distributed (in either source or binary form) with the major components
 * (compiler, kernel, and so on) of the operating system on which the
 * executable runs, unless that component itself accompanies the executable.
 *
 * If distribution of executable or object code is made by offering access to
 * copy from a designated place, then offering equivalent access to copy the
 * source code from the same place counts as distribution of the source code,
 * even though third parties are not compelled to copy the source along with
 * the object code.
 *
 * 4.  You may not copy, modify, sublicense, or distribute the Program except
 * as expressly provided under this License.  Any attempt otherwise to copy,
 * modify, sublicense or distribute the Program is void, and will
 * automatically terminate your rights under this License.  However, parties
 * who have received copies, or rights, from you under this License will not
 * have their licenses terminated so long as such parties remain in full
 * compliance.
 *
 * 5.  You are not required to accept this License, since you have not signed
 * it.  However, nothing else grants you permission to modify or distribute
 * the Program or its derivative works.  These actions are prohibited by law
 * if you do not accept this License.  Therefore, by modifying or distributing
 * the Program (or any work based on the Program), you indicate your
 * acceptance of this License to do so, and all its terms and conditions for
 * copying, distributing or modifying the Program or works based on it.
 *
 * 6.  Each time you redistribute the Program (or any work based on the
 * Program), the recipient automatically receives a license from the original
 * licensor to copy, distribute or modify the Program subject to these terms
 * and conditions.  You may not impose any further restrictions on the
 * recipients' exercise of the rights granted herein.  You are not responsible
 * for enforcing compliance by third parties to this License.
 *
 * 7.  If, as a consequence of a court judgment or allegation of patent
 * infringement or for any other reason (not limited to patent 
 * issues), conditions are imposed on you (whether by court 
 * order, agreement or otherwise) that contradict the conditions 
 * of this License, they do not excuse you from the conditions 
 * of this License.  If you cannot distribute so as to satisfy
 * simultaneously your obligations under this License and any other pertinent
 * obligations, then as a consequence you may not distribute the Program at
 * all.  For example, if a patent license would not permit royalty-free
 * redistribution of the Program by all those who receive copies directly or
 * indirectly through you, then the only way you could satisfy both it and
 * this License would be to refrain entirely from distribution of the Program.
 *
 * If any portion of this section is held invalid or unenforceable under any
 * particular circumstance, the balance of the section is intended to apply
 * and the section as a whole is intended to apply in other circumstances.
 *
 * It is not the purpose to this section to induce you to infringe any patents
 * or other property right claims or to contest validity of any such claims;
 * this section has the sole purpose of protecting the integrity of the free
 * software distribution system, which is implemented by public license
 * practices.  Many people have made generous contributions to the wide range
 * of software distributed through that system in reliance on consistent
 * application of that system; it is up to the author/donor to decide if he or
 * she is willing to distribute software through any other system and a
 * licensee cannot impose that choice.
 *
 * This section is intended to make thoroughly clear what is believed to be a
 * consequence of the rest of this License.
 *
 * 8.  If the distribution and/or use of the Program is restricted in certain
 * countries either by patents or by copyrighted interfaces, the original
 * copyright holder who places the Program under this License may add an
 * explicit geographical distribution limitation excluding those countries, so
 * that distribution is permitted only in or among countries not thus
 * excluded.  In such case, this License incorporates the limitation as if
 * written in the body of this License.
 *
 * 9.  The Free Software Foundation may publish revised and/or new versions of
 * the General Public License from time to time.  Such new versions will be
 * similar in spirit to the present version, but may differ in detail to
 * address new problems or concerns.
 *
 * Each version is given a distinguishing version number.  If the Program
 * specifies a version number of this License which applies to it and "any
 * later version," you have the option of following the terms and conditions
 * either of that version of any later version published by the Free Software
 * Foundation.  If the Program does not specify a version number of this
 * License, you may choose any version ever published by the Free Software
 * Foundation.
 *
 * 10.  If you wish to incorporate parts of the Program into other free
 * programs whose distribution conditions are different, write to the author
 * to ask for permission.  For software which is copyrighted by the Free
 * Software Foundation, write to the Free Software Foundation; we sometimes
 * make exceptions for this.  Our decision to grant permission will be guided
 * by the two goals of preserving the free status of all derivatives of our
 * free software and or promoting the sharing and reuse of software generally.
 *
 * NO WARRANTY
 *
 * 11.  BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
 * FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
 * OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
 * PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
 * OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
 * TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
 * PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
 * REPAIR OR CORRECTION.
 *
 * 12.  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
 * WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
 * REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
 * INCLUDING ANY GENERAL, SPECIAL INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
 * OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
 * TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
 * YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
 * PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * END OF TERMS AND CONDITIONS 
------------------------------------------------------------------------- */
