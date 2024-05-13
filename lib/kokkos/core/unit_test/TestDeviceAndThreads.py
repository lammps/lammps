#@HEADER
# ************************************************************************

#                        Kokkos v. 4.0
#       Copyright (2022) National Technology & Engineering
#               Solutions of Sandia, LLC (NTESS).

# Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.

# Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.

# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

# ************************************************************************
# @HEADER

import unittest
import subprocess
import platform
import os

PREFIX = "$<TARGET_FILE_DIR:Kokkos_CoreUnitTest_DeviceAndThreads>"
EXECUTABLE = "$<TARGET_FILE_NAME:Kokkos_CoreUnitTest_DeviceAndThreads>"
COMMAND = "/".join([PREFIX, EXECUTABLE])


def GetFlag(flag, *extra_args):
    p = subprocess.run([COMMAND, flag, *extra_args], capture_output=True)
    if p.returncode != 0:
        raise Exception(p.stderr.decode("utf-8"))
    return int(p.stdout)

def GetNumThreads(max_threads):
    args = []
    name = platform.system()
    if name == 'Darwin':
        args = ['sysctl', '-n', 'hw.physicalcpu_max']
    elif name == 'Linux':
        args = ['nproc', '--all']
    else:
        args = ['wmic', 'cpu', 'get', 'NumberOfCores']

    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = result.stdout.decode('utf-8')
    phys_cores_count = int(output)
    looplist = [1] + [i*phys_cores_count for i in [1,2,3,4,5,6,7]] \
        if GetFlag("hwloc_enabled") else [1,2,3,4,5]

    for x in looplist:
        if x >= max_threads:
            break
        yield x
    yield max_threads

class KokkosInitializationTestCase(unittest.TestCase):
    def test_num_threads(self):
        max_threads = GetFlag("max_threads")
        if max_threads == 1:
            self.skipTest("no host parallel backend enabled")
        for num_threads in GetNumThreads(max_threads):
            self.assertEqual(
                num_threads,
                GetFlag(
                    "num_threads",
                    "--kokkos-num-threads={}".format(num_threads)))

    def test_num_devices(self):
        if "KOKKOS_VISIBLE_DEVICES" in os.environ:
            self.skipTest("KOKKOS_VISIBLE_DEVICES environment variable is set")
        num_devices = GetFlag("num_devices")
        self.assertNotEqual(num_devices, 0)
        if num_devices == -1:
            self.skipTest("no device backend enabled")
        self.assertGreaterEqual(num_devices, 1)

    def test_device_id(self):
        if "KOKKOS_VISIBLE_DEVICES" in os.environ:
            self.skipTest("KOKKOS_VISIBLE_DEVICES environment variable is set")
        num_devices = GetFlag("num_devices")
        if num_devices == -1:
            self.assertEqual(-1, GetFlag("device_id"))
            self.skipTest("no device backend enabled")
        # by default use the first GPU available for execution
        self.assertEqual(0, GetFlag("device_id"))
        for device_id in range(num_devices):
            self.assertEqual(
                device_id,
                GetFlag(
                    "device_id",
                    "--kokkos-device-id={}".format(device_id)))

    def test_disable_warnings(self):
        self.assertEqual(0, GetFlag("disable_warnings"))
        self.assertEqual(
            0,
            GetFlag(
                "disable_warnings",
                "--kokkos-disable-warnings=0"))
        self.assertEqual(
            1,
            GetFlag(
                "disable_warnings",
                "--kokkos-disable-warnings=1"))

    def test_tune_internals(self):
        self.assertEqual(0, GetFlag("tune_internals"))
        self.assertEqual(
            0,
            GetFlag(
                "tune_internals",
                "--kokkos-tune-internals=0"))
        self.assertEqual(
            1,
            GetFlag(
                "tune_internals",
                "--kokkos-tune-internals=1"))


if __name__ == '__main__':
    unittest.main()
