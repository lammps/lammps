/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors:
     Ling-Ti Kong

   Contact:
     School of Materials Science and Engineering,
     Shanghai Jiao Tong University,
     800 Dongchuan Road, Minhang,
     Shanghai 200240, CHINA

     konglt@sjtu.edu.cn; konglt@gmail.com
------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(phonon,FixPhonon)

#else

#ifndef FIX_PHONON_H
#define FIX_PHONON_H

#include <complex>
#include "fix.h"
#include <map>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

namespace LAMMPS_NS {

class FixPhonon : public Fix {
 public:
  FixPhonon(class LAMMPS *, int, char **);
  ~FixPhonon();

  int  setmask();
  void init();
  void setup(int);
  void end_of_step();
  void post_run();
  double memory_usage();
  int modify_param(int, char **);

 private:
  int me,nprocs;
  bigint waitsteps;                             // wait these number of timesteps before recording atom positions
  bigint prev_nstep;                            // number of steps from previous run(s); to judge if waitsteps is reached.
  int nfreq, ifreq;                             // after this number of measurement (nfreq), the result will be output once
  int nx,ny,nz,nucell,ntotal;                   // surface dimensions in x- and y-direction, number of atom per unit surface cell
  int neval;                                    // # of evaluations
  int sysdim;                                   // system dimension
  int ngroup, nfind;                            // total number of atoms in group; total number of atoms on this proc
  char *prefix, *logfile;                       // prefix of output file names
  FILE *flog;
  
  double *M_inv_sqrt;

  class FFT3d *fft;                             // to do fft via the fft3d wraper
  int nxlo,nxhi,mysize;                         // size info for local MPI_FFTW
  int mynpt,mynq,fft_nsend;
  int *fft_cnts, *fft_disp;
  int fft_dim, fft_dim2;
  double *fft_data;
  
  int  itag, idx, idq;                          // index variables
  std::map<int,int> tag2surf, surf2tag;         // Mapping info

  double **RIloc;                               // R(r) and index on local proc
  double **RIall;                               // gathered R(r) and index
  double **Rsort;                               // sorted R(r)
  double **Rnow;                                // Current R(r) on local proc
  double **Rsum;                                // Accumulated R(r) on local proc

  int *recvcnts, *displs;                       // MPI related variables

  std::complex<double> **Rqnow;                 // Current R(q) on local proc
  std::complex<double> **Rqsum;                 // Accumulator for conj(R(q)_alpha)*R(q)_beta
  std::complex<double> **Phi_q;                 // Phi's on local proc
  std::complex<double> **Phi_all;               // Phi for all

  void readmap();                               // to read the mapping of gf atoms
  char *mapfile;                                // file name of the map file

  void getmass();                               // to get the mass of each atom in a unit cell

  int nasr;
  void postprocess();                           // to post process the data
  void EnforceASR();                            // to apply acoustic sum rule to gamma point force constant matrix

  char *id_temp;                                // compute id for temperature
  double *TempSum;                              // to get the average temperature vector
  double inv_nTemp;                             // inverse of number of atoms in temperature group
  class Compute *temperature;                   // compute that computes the temperature

  double hsum[6], **basis;
  int *basetype;

  // private methods to do matrix inversion
  void GaussJordan(int, std::complex<double>*);

};
}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix phonon command...

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No atom found for fix phonon!

Self-explanatory. Number of atoms in the group that was passed to
fix-phonon is less than 1.

E: Can not open output file %s"

Self-explanatory.

E: Illegal fix_modify command

Self-explanatory.

E: Could not find fix_modify temp ID

Self-explanatory.

E: Fix_modify temp ID does not compute temperature

Self-explanatory.

E: Cannot open input map file %s

Self-explanatory.

E: Error while reading header of mapping file!

Self-explanatory. The first line of the map file is expected to
contain 4 positive integer numbers.

E: FFT mesh and number of atoms in group mismatch!

Self-explanatory. The product of the 4 numbers should be exactly the
total number of atoms in the group that was passed to fix-phonon.

E: Error while reading comment of mapping file!

Self-explanatory. The second line of the map file should be a comment line.

E: The mapping is incomplete!

Self-explanatory.

E: Error while reading mapping file!

Self-explanatory.

E: The mapping info read is incorrect!

Self-explanatory.

E: Singular matrix in complex GaussJordan!

Self-explanatory.

W: More than one fix phonon defined

Self-explanatory. Just to warn that more than one fix-phonon is defined, but allowed.

*/
