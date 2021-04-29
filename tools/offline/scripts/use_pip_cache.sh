#!/bin/bash
# needs to be sourced

if [ -z "${PIP_CACHE_DIR}" ]
then
    echo "Must set PIP_CACHE_DIR environment variable"
    exit 1
fi

export PIP_NO_INDEX=1
export PIP_FIND_LINKS=$PIP_CACHE_DIR

echo "Disabling pip index and use local cache directory..."

function deactivate_pip_cache {
    echo "Removing local pip cache configuration..."
    unset PIP_NO_INDEX
    unset PIP_FIND_LINKS
}
