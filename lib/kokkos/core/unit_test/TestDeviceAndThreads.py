# SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

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
    for x in [1, 2] if GetFlag("hwloc_enabled") else [1, 2, 3, 4, 5]:
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
        self.assertGreaterEqual(GetFlag("device_id"), 0)
        self.assertLessEqual(GetFlag("device_id"), num_devices)
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
