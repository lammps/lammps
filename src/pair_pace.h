//
// Created by lysogy36 on 27.02.20.
//


#ifdef PAIR_CLASS

PairStyle(pace,PairPACE)

#else

#ifndef LMP_PAIR_TERSOFF_H
#define LMP_PAIR_TERSOFF_H

#include "pair.h"
#include "ace_evaluator.h"
#include "ace_c_basis.h"

namespace LAMMPS_NS {

    class PairPACE : public Pair {
    public:
        PairPACE(class LAMMPS *);

        virtual ~PairPACE();

        virtual void compute(int, int);

        void settings(int, char **);

        void coeff(int, char **);

        virtual void init_style();

        double init_one(int, int);

        // virtual double memory_usage();

    protected:
        ACECTildeBasisSet *basis_set = nullptr;

        ACECTildeEvaluator *ace = nullptr;

        char *potential_file_name;

        virtual void allocate();

        void read_files(char *, char *);

        inline int equal(double *x, double *y);


        double rcutmax;               // max cutoff for all elements
        int nelements;                // # of unique elements
        char **elements;              // names of unique elements

        int *map;                     // mapping from atom types to elements
        int *jlist_local;
        int *type_local;
    };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Tersoff requires atom IDs

This is a requirement to use the Tersoff potential.

E: Pair style Tersoff requires newton pair on

See the newton command.  This is a restriction to use the Tersoff
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/