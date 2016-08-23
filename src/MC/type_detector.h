/* ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Author: Amin Aramoon (Johns Hopkins U) aaramoo1@jhu.edu
------------------------------------------------------------------------- */

#ifndef LMP_TYPE_DETECTOR_H_
#define LMP_TYPE_DETECTOR_H_

#include <stdlib.h>
#include <vector>

namespace LAMMPS_NS {
	
enum{BOND,ANGLE,DIHEDRAL,IMPROPER};


class TypeDetector {

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
			ID = NULL;
		}
	};

public:
	int type; 
	int length;

	TypeDetector() {
	}

	virtual ~TypeDetector() {
		values.empty();
	}

	bool init(const char* input, int style) {
		bool status = true;
		std::vector<char*> args;
		char* s = new char[128];
		strcpy(s, input);
		tokenize(s, args);
		int len = args.size();
		type = style;
		if (type == BOND)
			length = 2;
		else if (type == ANGLE)
			length = 3;
		else if (type == DIHEDRAL)
			length = 4;
		else if (type == IMPROPER)
			length = 4;

		int ntype = len / (length + 1);
		int *t = new int[len];
		int val = 0;
		int iarg = 0;
		for (int num = 0; num < ntype; num++) {
            for (int i = 0; i < length; i++) {
				if (strcmp(args[iarg], "*") == 0)  t[i] = 0;
                else if (isdigit(args[iarg][0])) t[i] = atoi(args[iarg]);
                else status = false;
                iarg++;
            }
            val = atoi(args[iarg++]);			
            if (!status) break;
            this->set(t, val);
        }
        delete[] t;
        delete[] s;
        return status;        
	}

	bool get(int *&t1, int &answer) {
		for (std::vector<Entity *>::iterator it = values.begin();
				it != values.end(); it++) {
			if (equals((*it), t1)) {
				answer = (*it)->value;
				return true;
			}
		}
		return false;
	}

	int get_num_types() {
		int max = 0;
		for (std::vector<Entity *>::iterator it = values.begin();
				it != values.end(); it++) {
			if ((*it)->value > max) {
				max = (*it)->value;
			}
		}
		return max;
	}

private:

	void set(int *&t1, int value) {
		int i, p = 0;
		for (i = 0; i < length; ++i)
			p += (t1[i] == 0) ? 1 : 0;
		std::vector<Entity *>::iterator it = values.begin();
		bool equal = false;
		while (it != values.end()) {
			Entity *e = (*it);
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

	int equals(Entity *&e, int *&id) {
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

	void tokenize(char* &dup, std::vector<char*> & v) {
		const char delim[] = { ' ', ',', ';' };
		char *token = strtok(dup, delim);
		while (token != NULL) {
			v.push_back(token);
			token = strtok(NULL, delim);
		}
	}

	std::vector<Entity*> values;
};

}

#endif /* SRC_TYPE_DETECTOR_H_ */
