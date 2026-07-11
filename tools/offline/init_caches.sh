#!/bin/bash
echo "##############################################################################"
echo "Initializing LAMMPS offline compilation environment"
echo "##############################################################################"

if [ -z "${LAMMPS_CACHING_DIR}" ]
then
    export LAMMPS_CACHING_DIR="$HOME/.cache/lammps"
    echo "environment variable LAMMPS_CACHING_DIR not set"
    echo "Using default $LAMMPS_CACHING_DIR as cache directory..."
else
    echo "Using $LAMMPS_CACHING_DIR as cache directory..."
fi

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
CACHE_SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

if [ -z "${LAMMPS_DIR}" ]
then
    export LAMMPS_DIR="$(realpath $SCRIPT_DIR/../../)"
    echo "environment variable LAMMPS_DIR not set"
    echo "Using default $LAMMPS_DIR as LAMMPS distribution base directory..."
else
    echo "Using $LAMMPS_DIR as LAMMPS distribution base directory..."
fi

export GITHUB_PROXY_DIR="$LAMMPS_CACHING_DIR/github"
export LOGGING_DIR="$LAMMPS_CACHING_DIR/logs"
export PIP_CACHE_DIR="$LAMMPS_CACHING_DIR/pip"
export HTTP_CACHE_DIR="$LAMMPS_CACHING_DIR/http"

mkdir -p "$GITHUB_PROXY_DIR"
mkdir -p "$LOGGING_DIR"
mkdir -p "$PIP_CACHE_DIR"
mkdir -p "$HTTP_CACHE_DIR"

"${CACHE_SCRIPTS_DIR}/init_pip_cache.sh"
"${CACHE_SCRIPTS_DIR}/init_git_cache.sh"
"${CACHE_SCRIPTS_DIR}/init_http_cache.sh"
echo "##############################################################################"
echo
echo "To activate:"
echo "source \"${SCRIPT_DIR}/use_caches.sh\""
echo
echo "##############################################################################"
