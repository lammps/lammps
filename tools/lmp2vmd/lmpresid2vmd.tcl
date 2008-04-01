# small script to extract residue info from a lammps data file
# (c) 2008 Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>

proc lmpresidfromdata {mol filename} {

    if {"$mol" == "top"} {
        set mol [molinfo top]
    }

    # create an empty bondlist
    set na [molinfo $mol get numatoms]; # number of atoms of molecule
    set da 0;                           # number of atoms in data file
    set rl {};                          # resid list
    for {set i 0} {$i < $na} {incr i} {
        lappend rl 0
    }

    # open lammps data file
    if {[catch {open $filename r} fp]} {
        puts stderr "could not open file $filename"
        return -1
    }

    # read file line by line until we hit the Bonds keyword
    while {[gets $fp line] >= 0} {
        # pick number of atoms
        regexp {^\s*(\d+)\s+atoms} $line dummy da

        if { [regexp {^\s*Atoms} $line] } {
            # sanity check
            if {$na != $da} {
                puts stderr \
                    "number of atoms in VMD molecule ($na) does not match data file ($da)"
                return -1
            }

            puts "atoms= $da\nnow reading Atoms section"
            break
        }
    }

    # skip one line
    gets $fp line
    # read the Atoms data
    for {set i 0} {$i < $da} {incr i} {
        gets $fp line
        # grep bond numbers from entry and adjust to VMD numbering style
        regexp {^\s*(\d+)\s+(\d+)\s+\d+.*} $line dummy ba bb
        incr ba -1
        lset rl $ba $bb
    }
    close $fp

    set sel [atomselect $mol all]
    $sel set resid $rl
    $sel delete
}

