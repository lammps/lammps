/*
Copyright 2021 Yury Lysogorskiy^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
    This FILENAME is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//
// Created by Lysogorskiy Yury on 1.01.22.
//


#ifdef PAIR_CLASS

PairStyle(pace/al,PairPACEActiveLearning)

#else

#ifndef LMP_PAIR_PACE_AL_H
#define LMP_PAIR_PACE_AL_H

#include "pair.h"
#include "dump_custom.h"
#include "ace_b_evaluator.h"
#include "ace_b_basis.h"
#include "ace_recursive.h"

namespace LAMMPS_NS {

    class PairPACEActiveLearning : public Pair {
    public:
        PairPACEActiveLearning(class LAMMPS *);

        virtual ~PairPACEActiveLearning();

        virtual void compute(int, int);

        void settings(int, char **);

        void coeff(int, char **);

        virtual void init_style();

        double init_one(int, int);

        void *extract(const char *, int &);

        // virtual double memory_usage();

    protected:
        ACEBBasisSet *basis_set = nullptr;

        ACECTildeBasisSet *ctilde_basis_set = nullptr;

        ACEBEvaluator* ace = nullptr;

        ACERecursiveEvaluator* rec_ace = nullptr;

        char *potential_file_name;
        int gamma_grade_eval_freq = 1;

        DumpCustom *dump = nullptr;

        char *active_set_inv_filename;
        double gamma_lower_bound = 1.5;
        double gamma_upper_bound = 10;
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
    };

}

#endif
#endif