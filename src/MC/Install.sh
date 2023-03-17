# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# this is default Install.sh for all packages
# if package has an auxiliary library or a file with a dependency,
# then package dir has its own customized Install.sh

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
    if (test $mode = 0) then
        rm -f ../$1
    elif (! cmp -s $1 ../$1) then
        if (test -z "$2" || test -e ../$2) then
            cp $1 ..
            if (test $mode = 2) then
                echo "  updating src/$1"
            fi
        fi
    elif (test -n "$2") then
        if (test ! -e ../$2) then
            rm -f ../$1
        fi
    fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
    test -f ${file} && action $file
done


action fix_atom_swap.cpp
action fix_atom_swap.h
action fix_bond_break.cpp
action fix_bond_break.h
action fix_bond_create_angle.cpp
action fix_bond_create_angle.h
action fix_bond_create.cpp
action fix_bond_create.h
action fix_bond_swap.cpp
action fix_bond_swap.h
action fix_charge_regulation.cpp
action fix_charge_regulation.h
action fix_gcmc.cpp
action fix_gcmc.h
action fix_mol_swap.cpp
action fix_mol_swap.h
action fix_sgcmc.cpp   pair_eam.cpp
action fix_sgcmc.h     pair_eam.h
action fix_tfmc.cpp
action fix_tfmc.h
action fix_widom.cpp
action fix_widom.h
action pair_dsmc.cpp
action pair_dsmc.h
