#!/bin/bash
# needs to be sourced

if [ -z "${PIP_CACHE_DIR}" ]
then
    echo "Must set PIP_CACHE_DIR environment variable"
    exit 1
fi

URL_ENCODED_PATH="$(python3 -c "import urllib.parse; print(urllib.parse.quote(input()))" <<< "$PIP_CACHE_DIR")"

export PIP_NO_INDEX=1
export PIP_FIND_LINKS="file://$URL_ENCODED_PATH"

echo "Disabling pip index and use local cache directory..."

function deactivate_pip_cache {
    echo "Removing local pip cache configuration..."
    unset PIP_NO_INDEX
    unset PIP_FIND_LINKS
}
