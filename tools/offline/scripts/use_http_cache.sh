#!/bin/bash

if [ -z "${HTTP_CACHE_DIR}" ]
then
    echo "Must set HTTP_CACHE_DIR environment variable"
    exit 1
fi

if [ -z "${HTTP_CACHE_PORT}" ]
then
    HTTP_CACHE_PORT=8080
fi

until ! lsof -Pi :$HTTP_CACHE_PORT -sTCP:LISTEN -t >/dev/null
do
    echo "Port ${HTTP_CACHE_PORT} already in use, trying another..."
    ((HTTP_CACHE_PORT=HTTP_CACHE_PORT+1))
    sleep 1
done

if [ -z "${LOGGING_DIR}" ]
then
    echo "Must set LOGGING_DIR environment variable"
    exit 1
fi

python3 -m http.server $HTTP_CACHE_PORT --directory "$HTTP_CACHE_DIR" 2>&1 > "${LOGGING_DIR}/http.log" &
export HTTP_CACHE_PID=$!

export HTTP_CACHE_URL="http://localhost:$HTTP_CACHE_PORT"
export LAMMPS_HTTP_CACHE_CONFIG="$HTTP_CACHE_DIR/proxy.cmake"
echo "Running local HTTP cache server on $HTTP_CACHE_URL (pid: $HTTP_CACHE_PID)"

function deactivate_http_cache {
    echo "Shutting down HTTP cache server..."
    kill $HTTP_CACHE_PID
    unset HTTP_CACHE_PID
    unset LAMMPS_HTTP_CACHE_CONFIG
}
