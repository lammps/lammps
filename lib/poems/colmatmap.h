/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: colmatmap.h                                             *
 *      AUTHORS: See Author List                                           * 
 *      GRANTS: See Grants List                                            *
 *      COPYRIGHT: (C) 2005 by Authors as listed in Author's List          *
 *      LICENSE: Please see License Agreement                              *
 *      DOWNLOAD: Free at www.rpi.edu/~anderk5                             *
 *      ADMINISTRATOR: Prof. Kurt Anderson                                 *
 *                     Computational Dynamics Lab                          *
 *                     Rensselaer Polytechnic Institute                    *
 *                     110 8th St. Troy NY 12180                           * 
 *      CONTACT:        anderk5@rpi.edu                                    *
 *_________________________________________________________________________*/


#ifndef COLMATMAP_H
#define COLMATMAP_H

#include <iostream>

#include "virtualcolmatrix.h"
#include "virtualmatrix.h"

class ColMatrix;

class ColMatMap : public VirtualColMatrix  {
        double** elements;
public:
        ColMatMap();
        ~ColMatMap();
        ColMatMap(const ColMatMap& A);  // copy constructor
        ColMatMap(ColMatrix& A);  // copy constructor
        ColMatMap(int m);  // size constructor

        double& operator_1int (int i); // array access
        double Get_1int(int i) const;
        void Set_1int(int i, double value);
        double BasicGet_1int(int i) const;
        void BasicSet_1int(int i, double value);
        void BasicIncrement_1int(int i, double value);
        void SetElementPointer(int i, double* p);
        double* GetElementPointer(int i);
        void Dim(int m);
        void Const(double value);
        MatrixType GetType() const;
        std::ostream& WriteData(std::ostream& c) const;

        void AssignVM(const VirtualMatrix& A);
        ColMatMap& operator=(const ColMatMap& A); // assignment operator
        ColMatMap& operator=(const ColMatrix& A); // overloaded =
        ColMatMap& operator=(const VirtualColMatrix& A); // overloaded =
        ColMatMap& operator=(const VirtualMatrix& A); // overloaded =
        ColMatMap& operator*=(double b);

        // Fast Matrix Operators
        friend void FastAssign(ColMatMap& A, ColMatMap& C); //C = A
        friend void FastAssign(ColMatMap& A, ColMatrix& C);  // C = A
        friend void FastAssign(ColMatrix& A, ColMatMap& C);  // C = A
        friend void FastForwardEuler(ColMatMap& X,ColMatMap& Xdot, double dt);
        friend void FastForwardEuler(ColMatMap& X,ColMatrix& Xdot, double dt);
        friend void FastCKRK5(ColMatMap& X, ColMatrix& Xi, ColMatrix* f, double* c, double dt);
        friend void FastFRK5(ColMatMap& X, ColMatrix& Xi, ColMatrix* f, double* c, double dt);
};

#endif
