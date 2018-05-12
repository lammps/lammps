#!/bin/sh
# create lmpconfig.h header file from lmpconfig.h.in

mode="$1"

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# default settings
if type gzip > /dev/null 2>&1 ;
then \
    usegzip='LAMMPS_GZIP'
else
    usegzip='LAMMPS_GZIP_DISABLED'
fi
if type ffmpeg > /dev/null 2>&1 ;
then \
    useffmpeg='LAMMPS_FFMPEG'
else
    useffmpeg='LAMMPS_FFMPEG_DISABLED'
fi
singlefft='FFT_DOUBLE'
intsizes='LAMMPS_SMALLBIG'
exceptions='LAMMPS_EXCEPTIONS_DISABLED'
packfft='PACK_ARRAY'
memalign='LAMMPS_MEMALIGN_DEFAULT'
alignval=''

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
    intsizes=""
    while [  x"${intsizes}" = x"" ]
    do
        echo -n "Your choice (1/[2]/3): "
        read choice
        case ${choice} in
            1) intsizes='LAMMPS_SMALLSMALL' ;;
            2|"") intsizes='LAMMPS_SMALLBIG' ;;
            3) intsizes='LAMMPS_BIGBIG' ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    echo "using ${intsizes}"

    # gzip
    cat <<EOF
============================
Enable gzip compression via external gzip executable:
============================
With this feature enabled LAMMPS can read and write certain files
(e.g. data files, atom style and custom style dump files) in
in gzip compressed form. Requires an external gzip command.
EOF
    usegzip=''
    while [  x"${usegzip}" = x"" ]
    do
        echo -n 'Your choice ([y]/n): '
        read choice
        case ${choice} in
            y|"") usegzip='LAMMPS_GZIP'        ;;
            n) usegzip='LAMMPS_GZIP_DISABLED'  ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    [ "${usegzip}" = "LAMMPS_GZIP" ] && echo "gzip enabled" || echo "gzip disabled"

    # ffmpeg
    cat <<EOF
============================
Enable ffmpeg movie creation via external ffmpeg executable:
============================
With this feature enabled LAMMPS supports creating movies of trajectories
on-the-fly using the 'movie' dump style.  This is done by piping a stream
of images directly to an external ffmpeg command.
EOF
    useffmpeg=''
    while [  x"${useffmpeg}" = x"" ]
    do
        echo -n "Your choice ([y]/n): "
        read choice
        case ${choice} in
            y|"") useffmpeg='LAMMPS_FFMPEG' ;;
            n)    useffmpeg='LAMMPS_FFMPEG_DISABLED' ;;
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
    singlefft=''
    while [  x"${singlefft}" = x"" ]
    do
        echo -n "Your choice (y/[n]): "
        read choice
        case ${choice} in
            y) singlefft='FFT_SINGLE' ;;
            n|"") singlefft='FFT_DOUBLE' ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    test -n "${singlefft}" && echo "Using single precision FFTs" \
         || echo "Using double precision FFTs"

    # FFT data packing
    cat <<EOF
============================
Choose FFT data packing method:
============================
1) PACK_ARRAY:   use array indices to pack data (default)
2) PACK_POINTER: use pointers to pack data
3) PACK_MEMCPY:  use pointers and memcpy() to pack data

The pack data strategy is usually the most efficient.
EOF
    packfft=""
    while [  x"${packfft}" = x"" ]
    do
        echo -n "Your choice ([1]/2/3): "
        read choice
        case ${choice} in
            1|"") packfft='PACK_ARRAY' ;;
            2) packfft='PACK_POINTER' ;;
            3) packfft='PACK_MEMCPY' ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    echo "using ${packfft}"

    # exceptions
    cat <<EOF
============================
Use C++ exceptions instead of calling exit() or MPI_Abort()
============================
With this feature enabled LAMMPS will throw a C++ exception instead
of terminating. This is most useful when using LAMMPS as a library.
Through exception handling, the calling application will not crash,
but can either try to recover or exit more gracefully.
EOF
    exceptions=''
    while [  x"${exceptions}" = x"" ]
    do
        echo -n "Your choice (y/[n]): "
        read choice
        case ${choice} in
            y) exceptions='LAMMPS_EXCEPTIONS' ;;
            n|"") exceptions='LAMMPS_EXCEPTIONS_DISABLE' ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    if [ "${exceptions}" = "LAMMPS_EXCEPTIONS" ] ;
    then \
        echo "Using exceptions"
    else
        echo "Using exit() or MPI_Abort()"
    fi

    # memory alignment
    cat <<EOF
============================
Select preferred memory alignment for allocated memory
============================
Modern CPUs can benefit from blocks of allocated memory being
allocated on certain address boundaries. Often this preferred
alignment is larger than the default (64 bytes vs. 16 bytes).
You can set this prefernce below. A value of 64 bytes is recommended.
EOF
    memalign=''
    alignval=''
    while [  x"${memalign}" = x"" ]
    do
        echo -n "Your choice (0/32/[64]/128/256): "
        read choice
        case ${choice} in
            32|128|256)
                memalign='LAMMPS_MEMALIGN'
                alignval=${choice}
                ;;
            64|"")
                memalign='LAMMPS_MEMALIGN'
                alignval=64
                ;;
            0) memalign='LAMMPS_MEMALIGN_DISABLE' ;;
            *) echo "Invalid input: ${choice}" ;;
        esac
    done
    test -n "${alignval}" && echo "Using ${alignval} memory alignment" \
         || echo "Using default memory alignment"
fi

# edit template into actual header file

sed -e "s,@INTSIZES@,${intsizes}," \
    -e "s,@WITH_GZIP@,${usegzip}," \
    -e "s,@WITH_FFMPEG@,${useffmpeg}," \
    -e "s,@FFT_PRECISION@,${singlefft}," \
    -e "s,@PACK_FFT@,${packfft}", \
    -e "s,@WITH_EXCEPTIONS@,${exceptions}," \
    -e "s,@MEMALIGN@,${memalign} ${alignval}," \
    lmpconfig.h.in > lmpconfig.h.tmp

if [ "`diff --brief lmpconfig.h lmpconfig.h.tmp`" != "" ] ;
then \
   mv lmpconfig.h.tmp lmpconfig.h
else
    rm lmpconfig.h.tmp
fi
echo Global compile time configuration done
