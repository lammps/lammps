/*
 *_________________________________________________________________________*
 *      POEMS: PARALLELIZABLE OPEN SOURCE EFFICIENT MULTIBODY SOFTWARE     *
 *      DESCRIPTION: SEE READ-ME                                           *
 *      FILE NAME: virtualmatrix.h                                              *
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


#ifndef VIRTUALMATRIX_H
#define VIRTUALMATRIX_H
#include <iostream>

namespace POEMS {
enum MatrixType {
	MATRIX = 0,
	COLMATRIX = 1,
	ROWMATRIX = 2,
	MAT3X3 = 3,
	VECT3 = 4,
	MAT6X6 = 5,
	VECT6 = 6,
	COLMATMAP = 7,
	VECT4 = 8,
	MAT4X4 = 9
};

class VirtualMatrix {
protected:
	int numrows, numcols;
public:
	VirtualMatrix();
	virtual ~VirtualMatrix();
	int GetNumRows() const;
	int GetNumCols() const;

	double& operator() (int i, int j); // array access
	double Get(int i, int j) const;
	void Set(int i, int j, double value);
	double BasicGet(int i, int j) const;
	void BasicSet(int i, int j, double value);
	void BasicIncrement(int i, int j, double value);

	double& operator() (int i); // array access
	double Get(int i) const;
	void Set(int i, double value);
	double BasicGet(int i) const;
	void BasicSet(int i, double value);
	void BasicIncrement(int i, double value);

	virtual void Const(double value) = 0;
	virtual MatrixType GetType() const = 0;
	virtual void AssignVM(const VirtualMatrix& A) = 0;
	void Zeros();
	void Ones();
	virtual std::ostream& WriteData(std::ostream& c) const;
	virtual std::istream& ReadData(std::istream& c);

protected:
	virtual double& operator_2int(int i, int j) = 0;
	virtual double& operator_1int(int i);
	virtual double Get_2int(int i, int j) const = 0;
	virtual double Get_1int(int i) const ;
	virtual void Set_2int(int i, int j, double value) = 0;
	virtual void Set_1int(int i, double value);
	virtual double BasicGet_2int(int i, int j) const = 0;
	virtual double BasicGet_1int(int i) const ;
	virtual void BasicSet_2int(int i, int j, double value) = 0;
	virtual void BasicSet_1int(int i, double value);
	virtual void BasicIncrement_2int(int i, int j, double value) = 0;
	virtual void BasicIncrement_1int(int i, double value);

};

// overloaded operators
std::ostream& operator<< (std::ostream& c, const VirtualMatrix& A); // output
std::istream& operator>> (std::istream& c, VirtualMatrix& A); // input
}
#endif
