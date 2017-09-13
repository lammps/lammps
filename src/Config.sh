# create lmpconfig.h header file

mode="$1"

# enforce using portable C locale
LC_ALL=C
export LC_ALL

rm -f .tmp.lmpconfig.h.*
tmpfile=.tmp.lmpconfig.h.$$

# header of new config file
cat > ${tmpfile} <<EOF
/* ---------------------------------------------------------------------
   LAMMPS compile time configuration settings
   generated on: $(date)
--------------------------------------------------------------------- */
#ifndef LMPCONFIG_H
#define LMPCONFIG_H

#if defined(LAMMPS_SMALLSMALL) || defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#error "please use 'make config' to set integer size. don't define in makefile"
#endif
EOF

# use default settings

if [ x"${mode}" = x"default" ]
then
    intsize='LAMMPS_SMALLBIG'
else
    # define settings interactively
    cat <<EOF
============================
Choose integer size settings:
============================
1) small/small: 32-bit atom-id and molecule-id, 10-bit imageflags
                32-bit total number of atoms and timesteps.
   recommended setting for 32-bit machines

2) small/big:   32-bit atom-id and molecule-id, 10-bit imageflags
                64-bit total number of atoms and timesteps
   default setting, recommended for 64-bit machines unless big/big is needed

3) big/big:     64-bit atom-id and molecule-id, 21-bit imageflags
                64-bit total number of atoms and timesteps
                (max. number of atoms per MPI rank is always 32-bit)
EOF
    intsize=""
    while [  x"${intsize}" = x"" ]
    do
        echo -n "Your choice (1/[2]/3): "
        read choice
        case ${choice} in
            1) intsize='LAMMPS_SMALLSMALL' ;;
            2|"") intsize='LAMMPS_SMALLBIG'   ;;
            3) intsize='LAMMPS_BIGBIG'     ;; 
            *) echo "Invalid input: ${choice}"     ;;
        esac
    done
    echo "using ${intsize}"
fi

echo "#define ${intsize}" >> ${tmpfile}
echo '#endif' >> ${tmpfile}
mv ${tmpfile} lmpconfig.h
