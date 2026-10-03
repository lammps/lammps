#!/usr/bin/env python

# Script to build windows installer packages for LAMMPS
# (c) 2017,2018,2019,2020,2021,2022,2023,2024,2025,2026 Axel Kohlmeyer <akohlmey@gmail.com>

from __future__ import print_function
import sys,os,shutil,glob,re,subprocess,tarfile,gzip,time,inspect
try: from urllib.request import urlretrieve as geturl
except: from urllib import urlretrieve as geturl

try:
    import multiprocessing
    numcpus = multiprocessing.cpu_count()
except:
    numcpus = 1

# helper functions

def error(str=None):
    if not str: print(helpmsg)
    else: print(sys.argv[0],"ERROR:",str)
    sys.exit()

def getbool(arg,keyword):
    if arg in ['yes','Yes','Y','y','on','1','True','true']:
        return True
    elif arg in ['no','No','N','n','off','0','False','false']:
        return False
    else:
        error("Unknown %s option: %s" % (keyword,arg))

def fullpath(path):
    return os.path.abspath(os.path.expanduser(path))

def getexe(url,name):
    gzname = name + ".gz"
    geturl(url,gzname)
    with gzip.open(gzname,'rb') as gz_in:
      with open(name,'wb') as f_out:
        shutil.copyfileobj(gz_in,f_out)
    gz_in.close()
    f_out.close()
    os.remove(gzname)

def system(cmd):
    print('System command: ', cmd)
    try:
        txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    except subprocess.CalledProcessError as e:
        print("Command '%s' returned non-zero exit status" % e.cmd)
        print(e.output.decode('UTF-8', 'replace'))
        error("System command failed")
    return txt.decode('UTF-8')

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# record location and name of python script
homedir, exename = os.path.split(os.path.abspath(inspect.getsourcefile(lambda:0)))

# default settings help message and default settings

parflag = 'no'
pythonflag  = False
guiflag = False
revflag = system('git rev-parse --abbrev-ref HEAD').strip()
verbose = True
gitdir  = os.path.abspath(os.path.join(homedir,'..','..'))

helpmsg = """
Usage: python %s -p <mpi> -y <yes|no> -a <yes|no> -u <yes|no>

Flags (all flags are optional, defaults listed below):
  -p : select message passing parallel build (default value: %s)
    -p ms       : build an MPI parallel version with MS-MPI SDK 10.1
    -p no       : build a serial version using MPI STUBS library
  -y : select python support (default value: %s)
    -y yes      : build with python included
    -y no       : build without python
  -u : select whether to include the LAMMPS GUI (default value: %s)
    -u yes      : build includes LAMMPS GUI
    -u no       : build does not include LAMMPS GUI

Example:
  python %s -p ms -y no -u no
""" % (exename,parflag,pythonflag,guiflag,exename)

# parse arguments

argv = sys.argv
argc = len(argv)
i = 1

while i < argc:
    if i+1 >= argc:
        print("\nMissing argument to flag:",argv[i])
        error()
    elif argv[i] == '-p':
        parflag = argv[i+1]
    elif argv[i] == '-y':
        pythonflag = getbool(argv[i+1],"python")
    elif argv[i] == '-u':
        guiflag = getbool(argv[i+1],"gui")
    else:
        print("\nUnknown flag:",argv[i])
        error()
    i+=2

# checks
if parflag != 'no' and parflag != 'ms':
    error("Unsupported parallel flag %s" % parflag)
if pythonflag and guiflag:
    error("May only include either Python or LAMMPS GUI")

# test for valid revision name format: branch names, release tags, or commit hashes
rev1 = re.compile("^(stable|release|develop|maintenance)$")
rev2 = re.compile(r"^(patch|stable)_\d+(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\d{4}$")
use_gitrev=False
if not rev1.match(revflag) and not rev2.match(revflag):
    newflag=system('git rev-parse HEAD').strip()
    print("Using revision flag %s for branch %s" % (newflag,revflag))
    revflag=newflag
    use_gitrev=True

# create working directory
if pythonflag:
    builddir = os.path.join(fullpath('.'),"tmp-%s-%s-python" % (parflag,revflag))
elif guiflag:
    builddir = os.path.join(fullpath('.'),"tmp-%s-%s-gui" % (parflag,revflag))
else:
    builddir = os.path.join(fullpath('.'),"tmp-%s-%s-noadmin" % (parflag,revflag))
shutil.rmtree(builddir,True)
try:
    os.mkdir(builddir)
except:
    error("Cannot create temporary build folder: %s" % builddir)

# check for prerequisites and set up build environment
cc_cmd = which('x86_64-w64-mingw32-gcc')
cxx_cmd = which('x86_64-w64-mingw32-g++')
fc_cmd = which('x86_64-w64-mingw32-gfortran')
ar_cmd = which('x86_64-w64-mingw32-ar')
size_cmd = which('x86_64-w64-mingw32-size')
nsis_cmd = which('makensis')
lmp_size = 'smallbig'

print("""
Settings: building LAMMPS revision %s for 64-bit Windows 10+ with %d CPUs
Message passing  : %s
Home folder      : %s
Source folder    : %s
Build folder     : %s
C compiler       : %s
C++ compiler     : %s
Fortran compiler : %s
Library archiver : %s
""" % (revflag,numcpus,parflag,homedir,gitdir,builddir,cc_cmd,cxx_cmd,fc_cmd,ar_cmd))

# switch to build folder
os.chdir(builddir)

# download what is not automatically downloaded by CMake
print("Downloading third party tools")
url='http://download.lammps.org/thirdparty'
print("FFmpeg")
getexe("%s/ffmpeg-win64.exe.gz" % url,"ffmpeg.exe")
print("gzip")
getexe("%s/gzip.exe.gz" % url,"gzip.exe")

if parflag == "mpi" or parflag == "ms":
    mpiflag = "on"
else:
    mpiflag = "off"

print("Configuring build with CMake")
cmd = "mingw64-cmake -D CMAKE_BUILD_TYPE=Release"
cmd += " -C %s/cmake/presets/mingw-cross.cmake -S %s/cmake" % (gitdir,gitdir)
cmd += " -C %s/cmake/presets/kokkos-openmp.cmake" % gitdir
cmd += " -DBUILD_SHARED_LIBS=on -DBUILD_MPI=%s -DBUILD_OMP=ON" % mpiflag
if parflag == 'ms':
  cmd += " -DUSE_MSMPI=on"
if guiflag:
  cmd += " -DBUILD_LAMMPS_GUI=on -DDOWNLOAD_POTENTIALS=off -DQt6_DIR=/usr/x86_64-w64-mingw32/sys-root/mingw/lib/cmake/Qt6"
cmd += " -DPKG_INTEL=no"
cmd += " -DCMAKE_CXX_COMPILER_LAUNCHER=ccache"
cmd += " -DPKG_PLUGIN=yes"
cmd += " -DCMAKE_CXX_STANDARD=20"
if pythonflag: cmd += " -DPKG_PYTHON=yes"

print("Running: ",cmd)
txt = system(cmd)
if verbose: print(txt)

# create qt.conf file
if guiflag:
  with open("qt.conf", "w") as qtconf:
    qtconf.write("[Paths]\nPlugins = ../qt6plugins\n")
    qtconf.close()

print("Compiling")
system("cmake --build . --parallel %d" % numcpus)
print("Done")


# define version flag of the installers:
# - use current timestamp, when pulling from develop (for daily builds)
# - parse version from src/version.h when pulling from stable, release, or specific tag
# - otherwise use revflag, i.e. the commit hash
version = revflag
if revflag == 'stable' or revflag == 'release' or rev2.match(revflag):
  with open(os.path.join(gitdir,"src","version.h"),'r') as v_file:
    verexp = re.compile(r'^.*"(\w+) (\w+) (\w+)".*$')
    vertxt = v_file.readline()
    verseq = verexp.match(vertxt).groups()
    version = "".join(verseq)
elif revflag == 'develop' or revflag == 'maintenance':
    version = time.strftime('%Y-%m-%d')

if not pythonflag and not guiflag:
  print("Configuring pace plugin build with CMake")
  cmd = "mingw64-cmake -D CMAKE_BUILD_TYPE=Release"
  cmd += " -S %s/examples/PACKAGES/pace/plugin -B paceplugin" % gitdir
  cmd += " -DBUILD_SHARED_LIBS=on -DBUILD_MPI=%s -DBUILD_OMP=ON" % mpiflag
  cmd += " -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DLAMMPS_SOURCE_DIR=%s/src" % gitdir
  cmd += " -DLAMMPS_VERSION=%s" % version
  if parflag == 'ms': cmd += " -DUSE_MSMPI=on"
  cmd += " -DCMAKE_CXX_STANDARD=20"

  print("Running: ",cmd)
  txt = system(cmd)
  if verbose: print(txt)
  print("Done")

  print("Compiling and building installer")
  txt = system("cmake --build paceplugin --target package --parallel %d" % numcpus)
  if verbose: print(txt)
  for exe in glob.glob('paceplugin/LAMMPS*plugin*.exe'):
    shutil.move(exe,os.path.join('..',os.path.basename(exe)))
  print("Done")

  print("Configuring plumed plugin build with CMake")
  cmd = "mingw64-cmake -D CMAKE_BUILD_TYPE=Release"
  cmd += " -S %s/examples/PACKAGES/plumed/plugin -B plumedplugin" % gitdir
  cmd += " -DBUILD_SHARED_LIBS=on -DBUILD_MPI=%s -DBUILD_OMP=ON" % mpiflag
  cmd += " -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DLAMMPS_SOURCE_DIR=%s/src" % gitdir
  cmd += " -DLAMMPS_VERSION=%s" % version
  if parflag == 'ms': cmd += " -DUSE_MSMPI=on"
  cmd += " -DCMAKE_CXX_STANDARD=20"

  print("Running: ",cmd)
  txt = system(cmd)
  if verbose: print(txt)
  print("Done")

  print("Compiling and building installer")
  txt = system("cmake --build plumedplugin --target package --parallel %d" % numcpus)
  if verbose: print(txt)
  for exe in glob.glob('plumedplugin/LAMMPS*plugin*.exe'):
    shutil.move(exe,os.path.join('..',os.path.basename(exe)))
  print("Done")

  print("Cloning lammps-plugin package")
  if revflag == 'stable' or revflag == 'release' or rev2.match(revflag):
      txt = system("git clone -b %s --depth 1 git@github.com:lammps/lammps-plugins.git" % revflag)
  else:
      txt = system("git clone -b develop --depth 1 git@github.com:lammps/lammps-plugins.git")
  if verbose: print(txt)
  print("Configuring LAMMPS plugin collection build with CMake")
  cmd = "mingw64-cmake -D CMAKE_BUILD_TYPE=Release"
  cmd += " -S lammps-plugins -B build_plugins"
  cmd += " -DBUILD_SHARED_LIBS=on -DBUILD_MPI=%s -DBUILD_OMP=ON" % mpiflag
  cmd += " -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DLAMMPS_SOURCE_DIR=%s/src" % gitdir
  cmd += " -DLAMMPS_VERSION=%s" % version
  if parflag == 'ms': cmd += " -DUSE_MSMPI=on"
  cmd += " -DCMAKE_CXX_STANDARD=20"

  print("Running: ",cmd)
  txt = system(cmd)
  if verbose: print(txt)
  print("Done")

  print("Compiling and building installer")
  txt = system("cmake --build build_plugins --target package")
  if verbose: print(txt)
  for exe in glob.glob('build_plugins/LAMMPS*plugin*.exe'):
    shutil.move(exe,os.path.join('..',os.path.basename(exe)))
  print("Done")

print("Building PDF manual")
os.chdir(os.path.join(gitdir,"doc"))
txt = system("make upgrade")
if verbose: print(txt)
txt = system("make pdf")
if verbose: print(txt)
shutil.move("Manual.pdf",os.path.join(builddir,"LAMMPS-Manual.pdf"))
print("Done")

# switch back to build folder and copy/process files for inclusion in installer
print("Collect and convert files for the Installer package")
os.chdir(builddir)
shutil.copytree(os.path.join(gitdir,"examples"),os.path.join(builddir,"examples"),symlinks=False, ignore_dangling_symlinks=True)
shutil.copytree(os.path.join(gitdir,"bench"),os.path.join(builddir,"bench"),symlinks=False, ignore_dangling_symlinks=False)
shutil.copytree(os.path.join(gitdir,"tools"),os.path.join(builddir,"tools"),symlinks=False, ignore_dangling_symlinks=False)
shutil.copytree(os.path.join(gitdir,"python","lammps"),os.path.join(builddir,"python","lammps"),symlinks=False,ignore_dangling_symlinks=False)
shutil.copytree(os.path.join(gitdir,"potentials"),os.path.join(builddir,"potentials"),symlinks=False,ignore_dangling_symlinks=False)
shutil.copy(os.path.join(gitdir,"README"),os.path.join(builddir,"README.txt"))
shutil.copy(os.path.join(gitdir,"LICENSE"),os.path.join(builddir,"LICENSE.txt"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","colvars-refman-lammps.pdf"),os.path.join(builddir,"Colvars-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"tools","createatoms","Manual.pdf"),os.path.join(builddir,"CreateAtoms-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","kspace.pdf"),os.path.join(builddir,"Kspace-Extra-Info.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","pair_gayberne_extra.pdf"),os.path.join(builddir,"PairGayBerne-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","pair_resquared_extra.pdf"),os.path.join(builddir,"PairReSquared-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","PDLammps_overview.pdf"),os.path.join(builddir,"PDLAMMPS-Overview.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","PDLammps_EPS.pdf"),os.path.join(builddir,"PDLAMMPS-EPS.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","PDLammps_VES.pdf"),os.path.join(builddir,"PDLAMMPS-VES.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","SPH_LAMMPS_userguide.pdf"),os.path.join(builddir,"SPH-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","MACHDYN_LAMMPS_userguide.pdf"),os.path.join(builddir,"MACHDYN-Manual.pdf"))
shutil.copy(os.path.join(gitdir,"doc","src","PDF","CG-DNA.pdf"),os.path.join(builddir,"CG-DNA-Manual.pdf"))

# prune outdated inputs, too large files, or examples of packages we don't bundle
for d in ['accelerate','kim','mscg','PACKAGES/quip','PACKAGES/vtk']:
    shutil.rmtree(os.path.join("examples",d),True)
for d in ['FERMI','KEPLER']:
    shutil.rmtree(os.path.join("bench",d),True)
shutil.rmtree("tools/msi2lmp/test",True)
if os.path.exists("potentials/C_10_10.mesocnt"):
    os.remove("potentials/C_10_10.mesocnt")
if os.path.exists("potentials/TABTP_10_10.mesont"):
    os.remove("potentials/TABTP_10_10.mesont")
if os.path.exists("examples/PACKAGES/mesont/C_10_10.mesocnt"):
    os.remove("examples/PACKAGES/mesont/C_10_10.mesocnt")
if os.path.exists("examples/PACKAGES/mesont/TABTP_10_10.mesont"):
    os.remove("examples/PACKAGES/mesont/TABTP_10_10.mesont")

# convert text files to CR-LF conventions
txt = system("unix2dos LICENSE.txt README.txt tools/msi2lmp/README")
if verbose: print(txt)
txt = system("find bench examples potentials python tools/msi2lmp/frc_files -type f -print | xargs unix2dos || :")
if verbose: print(txt)
# mass rename README to README.txt
txt = system('for f in $(find tools bench examples potentials python -name README -print); do  mv -v $f $f.txt; done')
if verbose: print(txt)
# mass rename in.<name> to in.<name>.lmp
txt = system('for f in $(find bench examples -name in.\\* -print); do  mv -v $f $f.lmp; done')
if verbose: print(txt)
print("Done")

print("Configuring and building installer")
os.chdir(builddir)
if pythonflag:
    nsisfile = os.path.join(homedir,"lammps-python.nsis")
elif guiflag:
    nsisfile = os.path.join(homedir,"lammps-gui.nsis")
elif parflag == 'ms':
    nsisfile = os.path.join(homedir,"lammps-msmpi.nsis")
else:
    nsisfile = os.path.join(homedir,"lammps-noadmin.nsis")

shutil.copy(nsisfile,os.path.join(builddir,"lammps.nsis"))
shutil.copy(os.path.join(homedir,"FileAssociation.nsh"),os.path.join(builddir,"FileAssociation.nsh"))
shutil.copy(os.path.join(homedir,"lammps.ico"),os.path.join(builddir,"lammps.ico"))
shutil.copy(os.path.join(homedir,"lammps-text-logo-wide.bmp"),os.path.join(builddir,"lammps-text-logo-wide.bmp"))

mingwdir = '/usr/x86_64-w64-mingw32/sys-root/mingw/bin/'

if parflag == 'ms':
    txt = system("makensis -DMINGW=%s -DVERSION=%s-MSMPI -DLMPREV=%s lammps.nsis" % (mingwdir,version,revflag))
    if verbose: print(txt)
else:
    txt = system("makensis -DMINGW=%s -DVERSION=%s -DLMPREV=%s lammps.nsis" % (mingwdir,version,revflag))
    if verbose: print(txt)

# clean up after successful build
os.chdir('..')

# check GUI package for missing dll files
if guiflag:
    print("Checking GUI installer for missing DLL files")
    unpack = os.path.join(builddir, "dllcheck")
    for exe in glob.glob('LAMMPS-*GUI-%s*.exe' % version):
        print("Checking: ", exe);
        shutil.rmtree(unpack, True)
        system("7z x -y -o%s %s" % (unpack, exe))
        dirs = [os.path.join(unpack, d) for d in ('bin', 'qt6plugins', 'plugins')
               if os.path.isdir(os.path.join(unpack, d))]
        print(system("%s --sysroot %s %s" % (os.path.join(homedir, 'check_missing_dlls.py'), mingwdir, ' '.join(dirs))))
    shutil.rmtree(unpack, True)

print("Cleaning up...")
shutil.rmtree(builddir,True)

print("Done.")

