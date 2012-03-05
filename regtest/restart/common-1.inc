
clear

include ${sfx}.inc

units lj
atom_style ${astyle}
special_bonds lj/coul 1.0 1.0 1.0
newton ${newtn}

read_data ../data/data.${astyle}-${tst}
include ${pstyle}.inc

write_restart ${pstyle}.restart-${sfx}-${tst}

shell ${r2d} ${pstyle}.restart-${sfx}-${tst} &
             ${pstyle}.data-${sfx}-${tst}    &
             ${coeff} > /dev/null            &
      || echo Error processing ${pstyle} ${sfx} ${tst}

