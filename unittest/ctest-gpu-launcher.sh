#!/bin/sh
# Launcher for running the unit tests of a GPU enabled build concurrently on
# multiple NVIDIA GPUs.  With a resource spec file, CTest assigns a slot on one
# of the listed GPUs to each test and passes it in the environment, e.g. as
# CTEST_RESOURCE_GROUP_0_GPUS=id:2,slots:1 .  This makes the GPU with that id,
# i.e. its nvidia-smi index, the only visible one and then runs the test.
# Without resource allocation, the test is run unchanged.
#
# usage: ctest-gpu-launcher.sh <command> [argument]...

if [ -n "${CTEST_RESOURCE_GROUP_0_GPUS:-}" ]
then
    gpu="${CTEST_RESOURCE_GROUP_0_GPUS#id:}"
    gpu="${gpu%%,*}"
    # nvidia-smi numbers the GPUs by PCI bus ID, but CUDA by default does not
    CUDA_DEVICE_ORDER=PCI_BUS_ID
    CUDA_VISIBLE_DEVICES="${gpu}"
    export CUDA_DEVICE_ORDER CUDA_VISIBLE_DEVICES
fi
exec "$@"
