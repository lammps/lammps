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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "priority_list.h"

using namespace std;

PriorityList::PriorityList(char* tag, CoordinateType t) {

	type = t;

	switch (type) {
	case bond:
		length = 2;
		break;
	case angle:
		length = 3;
		break;
	case dihedral:
		length = 4;
		break;
	case improper:
		length = 4;
		break;
	}

	int n = strlen(tag) + 1;
	id = new char[n];
	strcpy(id, tag);
}

PriorityList::~PriorityList() {
	values.empty();
}

void PriorityList::set(int * t1, int value) {
	int i, p = 0;
	for (i = 0; i < length; ++i)
		p += (t1[i] == 0) ? 1 : 0;
	vector<Entity*>::iterator it = values.begin();
	bool equal = false;
	while (it != values.end()) {
		Entity * e = (*it);
		if (equals(e, t1) && e->level == p) {
			equal = true;
			break;
		} else if (e->level > p) {
			break;
		}
		it++;
	}
	if (equal)
		(*it)->value = value;
	else {
		Entity *ne = new Entity();
		ne->ID = new int[length];
		for (int i = 0; i < length; ++i)
			ne->ID[i] = t1[i];
		ne->value = value;
		ne->level = p;
		values.insert(it, ne);
	}
}

bool PriorityList::get(int *& t1, int & answer) {
	for (vector<Entity*>::iterator it = values.begin(); it != values.end();
			it++) {
		if (equals((*it), t1)) {
			answer = (*it)->value;
			return true;
		}
	}
	return false;
}

int PriorityList::equals(Entity *& e, int *& id) {
	int *ID = e->ID;
	int res = 1, res_rev = 1;
	for (int i = 0; i < length; ++i) {
		res *= (ID[i] == id[i] || ID[i] == 0) ? 1 : 0;
		res_rev *= (ID[i] == id[length - i - 1] || ID[i] == 0) ? 1 : 0;
	}
	if (res == 1 || res_rev == 1)
		return 1;
	return 0;
}
