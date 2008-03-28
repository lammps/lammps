# small script to extract bonding info from a lammps data file
# so that VMD will display the correct bonds for CG-MD.
# (c) 2007 Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>

proc lmpbondsfromdata {mol filename} {

    # create an empty bondlist
    set na [molinfo $mol get numatoms];     # number of atoms
    set nb 0;                           # number of bonds
    set bl {};                          # bond list
    for {set i 0} {$i < $na} {incr i} {
        lappend bl {}
    }

    # open lammps data file
    if {[catch {open $filename r} fp]} {
        puts stderr "could not open file $filename"
        return -1
    }

    # read file line by line until we hit the Bonds keyword
    while {[gets $fp line] >= 0} {
        # pick number of bonds
        regexp {^\s*(\d+)\s+bonds} $line dummy nb

        if { [regexp {^\s*Bonds} $line] } {
            puts "nbonds= $nb\n now reading Bonds section"
            break
        }
    }

    # skip one line
    gets $fp line
    # read the bonds file
    for {set i 0} {$i < $nb} {incr i} {
        gets $fp line
        # grep bond numbers from entry and adjust to VMD numbering style
        regexp {^\s*\d+\s+\d+\s+(\d+)\s+(\d+)} $line dummy ba bb
        incr ba -1
        incr bb -1

        set bn [lindex $bl $ba]
        lappend bn $bb
        set bl [lreplace $bl $ba $ba $bn]

        set bn [lindex $bl $bb]
        lappend bn $ba
        set bl [lreplace $bl $bb $bb $bn]
    }
    close $fp

    set sel [atomselect $mol all]
    $sel setbonds $bl
    $sel delete
}

