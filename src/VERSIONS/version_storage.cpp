#include "version_storage.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;
VStorage::VStorage(LAMMPS *lmp) : Pointers(lmp){
    current_typeset = 0;
    ntype_sets = 0;
    nactive_typesets = 0;
    }

VStorage::~VStorage(){
  for (int i = 0; i < ntype_sets; i++){
      memory->destroy(typeset_map[i]);
  }
}

//Create a new type set. On the first call type set 0 (existing atom types) and type set 1 (new types) are created
int VStorage::add_typeset(const int *types){
	if(!ntype_sets){
		int *initial_replacement = nullptr;
    	initial_replacement = memory->grow(initial_replacement, lmp->atom->nmax, "atom:typeset");
    	for(int i = 0; i < lmp->atom->nmax; i++){
			initial_replacement[i] = lmp->atom->type[i];
    		}
		    typeset_map[ntype_sets++] = initial_replacement;
			nactive_typesets++;
		}

	int *sub_types = nullptr;
 	int m;
	sub_types = memory->grow(sub_types, lmp->atom->nmax, "atom:typeset");
	for(int i = 0; i < lmp->atom->natoms; i++){
        if((m = lmp->atom->map(i+1)) >= 0){
      	    sub_types[m] = types[i];
            }
        }
	typeset_map[ntype_sets] = sub_types;
	nactive_typesets++;
	return ntype_sets++;
	}

//Delete a particular type set
int VStorage::delete_typeset(int typeset_id){
	int res = -1;
	if(typeset_id != current_typeset){
    	memory->destroy(typeset_map[typeset_id]);
    	nactive_typesets--;
    	res = 0;
    	}
	return res;
	}
