#!/bin/bash

if [ -z "${GITHUB_PROXY_DIR}" ]
then
    echo "Must set GITHUB_PROXY_DIR environment variable"
    exit 1
fi

mkdir -p "$GITHUB_PROXY_DIR"
cd "$GITHUB_PROXY_DIR"

PROJECTS=(
    akohlmey/sphinx-fortran
    mathjax/MathJax
)

for project in ${PROJECTS[@]}
do
    GH_NAMESPACE="$(dirname $project)"
    GH_PROJECT="$(basename $project)"
    mkdir -p "$GH_NAMESPACE"
    git clone --bare "https://github.com/$GH_NAMESPACE/$GH_PROJECT.git" "$GH_NAMESPACE/$GH_PROJECT.git"
done
