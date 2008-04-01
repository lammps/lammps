# small script to assign a radius by type number in VMD
# (c) 2008 Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>

proc lmptypetoradius {mol rlist} {
    if {"$mol" == "top"} {
        set mol [molinfo top]
    }

    set t 0
    foreach r $rlist {
        incr t
        set sel [atomselect $mol "type $t"]
        $sel set radius $r
        $sel delete
    }
    return 0
}

