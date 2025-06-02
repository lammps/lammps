#ifndef GLOBAL_H
#define GLOBAL_H

#define ZERO 1.0e-8
#define MAXLINE 512

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

// one can customize the following parameters
#define QSTEP 0.02      // Step size when evaluating phonon dispersion automatically
#define NUMATOM 10      // Maximum # of atoms that will be displayed when printing basis info
#endif
