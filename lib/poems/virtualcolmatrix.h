/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: virtualcolmatrix.h                                      *
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


#ifndef VIRTUALCOLMATRIX_H
#define VIRTUALCOLMATRIX_H

#include "virtualmatrix.h"


class VirtualColMatrix : public VirtualMatrix  {
public:
        VirtualColMatrix();
        ~VirtualColMatrix() = default;
        double& operator_2int (int i, int j); // array access
        double Get_2int (int i, int j) const;
        void Set_2int (int i, int j, double value);
        double BasicGet_2int(int i, int j) const;
        void BasicSet_2int(int i, int j, double value);
        void BasicIncrement_2int(int i, int j, double value);

        virtual double& operator_1int (int i) = 0; // array access
        virtual double Get_1int(int i) const = 0;
        virtual void Set_1int(int i, double value) = 0;
        virtual double BasicGet_1int(int i) const = 0;
        virtual void BasicSet_1int(int i, double value) = 0;
        virtual void BasicIncrement_1int(int i, double value) = 0;

};

#endif
