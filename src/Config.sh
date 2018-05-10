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
#error "Please use 'make config' to set integer size. don't define in makefile"
#endif

EOF

# default settings
type gzip > /dev/null 2>&1 && usegzip='LAMMPS_GZIP' || usegzip=''
type ffmpeg > /dev/null 2>&1 && useffmpeg='LAMMPS_FFMPEG' || useffmpeg=''
singlefft=''
intsize='LAMMPS_SMALLBIG'

if [ x"${mode}" != x"default" ]
then
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
    # gzip
    cat <<EOF
============================
Enable gzip compression via external gzip executable:
============================
With this feature enabled LAMMPS can read and write certain files
(e.g. data files, atom style and custom style dump files) in
in gzip compressed form. Requires an external gzip command.
EOF
    usegzip='x'
    while [  x"${usegzip}" = x"x" ]
    do
        echo -n "Your choice ([y]/n): "
        read choice
        case ${choice} in
            y|"") usegzip='LAMMPS_GZIP'        ;;
            n) usegzip=''                      ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    test -n "${usegzip}" && echo "gzip enabled" || echo "gzip disabled"
    # ffmpeg
    cat <<EOF
============================
Enable ffmpeg movie creation via external ffmpeg executable:
============================
With this feature enabled LAMMPS supports creating movies of trajectories
on-the-fly using the 'movie' dump style.  This is done by piping a stream
of images directly to an external ffmpeg command.
EOF
    useffmpeg='x'
    while [  x"${useffmpeg}" = x"x" ]
    do
        echo -n "Your choice ([y]/n): "
        read choice
        case ${choice} in
            y|"") useffmpeg='LAMMPS_FFMPEG' ;;
            n) useffmpeg=''                      ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    test -n "${useffmpeg}" && echo "ffmpeg enabled" || echo "ffmpeg disabled"
    # single precision FFT
    cat <<EOF
============================
Enable single precision FFTs
============================
With this feature enabled LAMMPS will do 3d-FFTs (e.g. in the PPPM kspace
style) with single precision floating point numbers instead of double
precision. This is a trade-off between accuracy and communication speed.
It is not recommended to enable this feature, unless you are running with
a large number of MPI ranks (100s or 1000s) and kspace performance becomes
a bottleneck. Also using run style verlet/split or using mixed OpenMP/MPI
parallelization via the USER-OMP or USER-INTEL packages should be tried
first for improved (parallel) performance.
EOF
    singlefft='x'
    while [  x"${singlefft}" = x"x" ]
    do
        echo -n "Your choice (y/[n]): "
        read choice
        case ${choice} in
            y) singlefft='FFT_SINGLE' ;;
            n|"") singlefft=''                      ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    test -n "${singlefft}" && echo "Using single precision FFTs" \
         || echo "Using double precision FFTs"
fi

echo "#define ${intsize}" >> ${tmpfile}
test -n "${usegzip}" && echo "#define ${usegzip}" >> ${tmpfile}
test -n "${useffmpeg}" && echo "#define ${useffmpeg}" >> ${tmpfile}
test -n "${singlefft}" && echo "#define ${singlefft}" >> ${tmpfile}
echo '#endif' >> ${tmpfile}
mv ${tmpfile} lmpconfig.h
echo Global compile time configuration done
