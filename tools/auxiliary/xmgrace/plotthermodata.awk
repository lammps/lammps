#!/sw/bin/awk -f

# Usage: ./plotthermodata.awk log.lammps | xmgrace -nxy -

BEGIN {
    if ( stepnum == 0 ) stepnum = 1  # if unset, pick first set of output
    print "@xaxis  label \"Time Step\""
}
/Step/ {
    num += 1
    if ( num != stepnum ) next
    n = 0
    for (i=2; i<=NF; i++) {
        print "@s" n, "legend \"" $i "\""
        n += 1
    }
}
/Step/,/step/ { if ( $0 ~ /[sS]tep/ ) next; else print }
