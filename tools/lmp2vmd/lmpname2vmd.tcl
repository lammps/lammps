# small script to assign atom names to type numbers in LAMMPS .
# (c) 2008 Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>

proc lmptypetoname {mol names} {
    if {"$mol" == "top"} {
        set mol [molinfo top]
    }

    set t 0
    foreach n $names {
        incr t
        set sel [atomselect $mol "type $t"]
        $sel set name $n
        $sel delete
    }
    return 0
}

