#!/usr/bin/env python3
"""Verify that a Windows binary bundle contains all required DLL files.

Scans all .exe and .dll files in the given directories with objdump,
collects the DLLs they import, and reports every import that is neither
bundled in the application directory nor a well-known Windows system DLL.

The first directory is the application directory (usually "bin"): the
Windows loader resolves DLL names for the whole process from there first,
so it is both scanned and used as the set of provided DLLs.  Additional
directories (Qt plugins, LAMMPS plugins) contain binaries that are loaded
into the same process by absolute path: they are scanned for imports, but
do not provide DLL name resolution.

Typical use on a built installer (NSIS installers unpack with 7z):

    7z x -y -oUNPACKED LAMMPS-64bit-GUI-4Jul2026.exe
    ./check_missing_dlls.py UNPACKED/bin UNPACKED/qt6plugins UNPACKED/plugins

Exit status: 0 = complete, 1 = missing DLLs found, 2 = error.
"""

import argparse
import os
import subprocess
import sys

# Windows system DLLs and API sets that must not be bundled.
# Merged from the skip list in lammps-gui's build_windows_cross_nsis.sh.
SYSTEM_DLLS = {
    'advapi32.dll', 'authz.dll', 'avicap32.dll', 'avrt.dll', 'bcrypt.dll',
    'cabinet.dll', 'cfgmgr32.dll', 'comctl32.dll', 'comdlg32.dll',
    'crypt32.dll', 'd2d1.dll', 'd3d9.dll', 'd3d11.dll', 'd3d12.dll',
    'dbghelp.dll', 'dnsapi.dll', 'dwmapi.dll', 'dwrite.dll', 'dxgi.dll',
    'gdi32.dll', 'gdiplus.dll', 'imm32.dll', 'iphlpapi.dll', 'kernel32.dll',
    'mpr.dll', 'msvcrt.dll', 'msimg32.dll', 'ncrypt.dll', 'netapi32.dll',
    'ntdll.dll', 'ole32.dll', 'oleaut32.dll', 'opengl32.dll', 'psapi.dll',
    'secur32.dll', 'setupapi.dll', 'shcore.dll', 'shell32.dll', 'shlwapi.dll',
    'user32.dll', 'userenv.dll', 'uxtheme.dll', 'version.dll', 'winhttp.dll',
    'winmm.dll', 'winspool.drv', 'wldap32.dll', 'ws2_32.dll', 'wsock32.dll',
    'wtsapi32.dll'
}
SYSTEM_PREFIXES = ('api-ms-win-', 'ext-ms-')

def is_system_dll(name):
    return name in SYSTEM_DLLS or name.startswith(SYSTEM_PREFIXES)


def find_binaries(topdir):
    """Return all .exe and .dll files below topdir."""
    binaries = []
    for root, _, files in os.walk(topdir):
        for f in files:
            if f.lower().endswith(('.exe', '.dll')):
                binaries.append(os.path.join(root, f))
    return sorted(binaries)


def imported_dlls(binary, objdump):
    """Return the DLL names (lower case) imported by a PE binary."""
    try:
        out = subprocess.run([objdump, '-p', binary], capture_output=True,
                             text=True, check=True).stdout
    except (subprocess.CalledProcessError, OSError) as err:
        sys.exit('cannot get imports of %s: %s' % (binary, err))
    imports = set()
    for line in out.splitlines():
        if 'DLL Name:' in line:
            imports.add(line.split('DLL Name:')[1].strip().lower())
    return imports


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('appdir',
                        help='directory with the executables and bundled DLLs')
    parser.add_argument('scandirs', nargs='*',
                        help='additional directories with loaded binaries '
                             '(Qt plugins, LAMMPS plugins)')
    parser.add_argument('--sysroot', metavar='DIR',
                        default='/usr/x86_64-w64-mingw32/sys-root/mingw/bin',
                        help='where to look for copies of missing DLLs '
                             '(default: %(default)s)')
    parser.add_argument('--objdump', default='objdump',
                        help='objdump executable to use (default: %(default)s)')
    parser.add_argument('--skip', action='append', default=[], metavar='NAME',
                        help='treat this DLL name as a system DLL (repeatable)')
    args = parser.parse_args()

    for name in args.skip:
        SYSTEM_DLLS.add(name.lower())

    if not os.path.isdir(args.appdir):
        sys.exit('no such directory: %s' % args.appdir)

    provided = {f.lower() for f in os.listdir(args.appdir)
                if f.lower().endswith('.dll')}
    binaries = find_binaries(args.appdir)
    for d in args.scandirs:
        if not os.path.isdir(d):
            sys.exit('no such directory: %s' % d)
        binaries += find_binaries(d)

    # dll name -> list of binaries that import it
    missing = {}
    for binary in binaries:
        for dll in imported_dlls(binary, args.objdump):
            if dll not in provided and not is_system_dll(dll):
                missing.setdefault(dll, []).append(os.path.basename(binary))

    if not missing:
        print('OK: all imports of %d binaries satisfied by %d bundled DLLs'
              % (len(binaries), len(provided)))
        return 0

    for dll in sorted(missing):
        print('MISSING: %s' % dll)
        print('  needed by: %s' % ', '.join(sorted(set(missing[dll]))))
        for f in sorted(os.listdir(args.sysroot)) if os.path.isdir(args.sysroot) else []:
            if f.lower() == dll:
                print('  candidate: %s' % os.path.join(args.sysroot, f))
    print('ERROR: %d missing DLL file(s)' % len(missing))
    return 1


if __name__ == '__main__':
    sys.exit(main())
