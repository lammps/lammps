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

# all package files, the GRANULAR sub-models for body/rounded/polyhedron
# and wall/body/polyhedron require the GRANULAR package

action body_nparticle.cpp
action body_nparticle.h
action body_rounded_polygon.cpp
action body_rounded_polygon.h
action body_rounded_polyhedron.cpp
action body_rounded_polyhedron.h
action compute_body_local.cpp
action compute_body_local.h
action compute_temp_body.cpp
action compute_temp_body.h
action fix_nh_body.cpp
action fix_nh_body.h
action fix_nph_body.cpp
action fix_nph_body.h
action fix_npt_body.cpp
action fix_npt_body.h
action fix_nve_body.cpp
action fix_nve_body.h
action fix_nvt_body.cpp
action fix_nvt_body.h
action fix_wall_body_polygon.cpp
action fix_wall_body_polygon.h
action fix_wall_body_polyhedron.cpp
action fix_wall_body_polyhedron.h
action fix_wall_body_polyhedron_granular.cpp granular_model.cpp
action fix_wall_body_polyhedron_granular.h   granular_model.h
action pair_body_nparticle.cpp
action pair_body_nparticle.h
action pair_body_rounded_polygon.cpp
action pair_body_rounded_polygon.h
action pair_body_rounded_polyhedron.cpp
action pair_body_rounded_polyhedron_granular.cpp granular_model.cpp
action pair_body_rounded_polyhedron_granular.h   granular_model.h
action pair_body_rounded_polyhedron.h
