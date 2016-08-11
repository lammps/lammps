/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */
 
 /*-----------------------Contribution authors------------------------------
  * 			Amin Aramoon, aaramoo1@jhu.edu
  * -------------------------------------------------------------------- */

#ifndef LMP_PRIORITY_LIST_H_
#define LMP_PRIORITY_LIST_H_

#include <vector>

enum CoordinateType {
	bond, angle, dihedral, improper
};

struct Entity {
	int *ID; // * -> -1
	int value;
	int level;
	Entity() {
		ID = NULL;
		value = -1;
		level = 0;
	}
	virtual ~Entity() {
		delete[] ID;
	}
};

class PriorityList {
public:
	char* id;
	CoordinateType type;
	int length;
	PriorityList(char*, CoordinateType);
	virtual ~PriorityList();
	void set(int *, int);
	bool get(int *&, int&);
private:
	int equals(Entity *&, int *&);
	std::vector<Entity*> values;
};

#endif /* PRIORITY_LIST_H_ */
