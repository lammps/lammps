//
// Created by lysogy36 on 27.02.20.
//


#ifdef PAIR_CLASS

PairStyle(pace,PairPACE)

#else

#ifndef LMP_PAIR_PACE_H
#define LMP_PAIR_PACE_H

#include "pair.h"
#include "ace_evaluator.h"
#include "ace_recursive.h"
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

        void *extract(const char *, int &);

        // virtual double memory_usage();

    protected:
        ACECTildeBasisSet *basis_set = nullptr;

        ACERecursiveEvaluator *ace = nullptr;

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
        double **scale;

        bool recursive = false;       // "recursive" option for ACERecursiveEvaluator
    };

}

#endif
#endif