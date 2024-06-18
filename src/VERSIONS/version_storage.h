#ifndef LMP_VERSION_STORAGE_H
#define LMP_VERSION_STORAGE_H
#endif

#include "pointers.h"
#include "atom.h"

#include <map>

namespace LAMMPS_NS{

class VStorage : protected Pointers{
    public:
        std::map<int, int*> typeset_map;
        int current_typeset;
        int ntype_sets;
        int nactive_typesets;

        VStorage(class LAMMPS*);
        ~VStorage() override;
        int add_typeset(const int *types);
        int delete_typeset(int typeset_id);
    };
    }
