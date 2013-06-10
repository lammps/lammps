#!/usr/bin/env python
# if necessary, edit preceding line to point to your Python
# or launch as "python Make.py ..."

# Purpose: manage LAMMPS packages, external libraries, and builds
#          create a version of LAMMPS specific to an input script(s)
# Syntax: Make.py switch args ...
# Help: type Make.py or Make.py -h

# -----------------------------------------------------------------------
# LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
# http://lammps.sandia.gov, Sandia National Laboratories
# Steve Plimpton, sjplimp@sandia.gov
#
# Copyright (2003) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.
#
# See the README file in the top-level LAMMPS directory.
# -----------------------------------------------------------------------

# hi-level tasks that can be performed:

# install/uninstall packages and build the associated external libs
#   use -p and -u and -e
# install packages needed for input script(s)
#   use -i and -p
# create new dir with only the source code needed for input script(s)
#   use -i and -n
# build LAMMPS, either in src or new dir
#   use -b

# data structures:

# packstd list = ["COLLOID", "GPU", ...]
# packuser list = ["USER-SPH", "USER-EWALDN", ...]
# cfiles,hfiles dicts of LAMMPS *.cpp,*.h files =
#   {"fix_nve.cpp": "", "angle_diplole.cpp": "USER-MISC", ...}
# cdep,hdep dicts of files with list of LAMMPS *.h files they include = 
#   {"fix_nve.cpp": ['fix_nve.h', 'atom.h', 'force.h', ...'], ...}
# creq,hreq lists of src files needed to build LAMMPS for input scripts =
#   ["main.cpp", "atom.cpp", ...]

# ISSUES:

# see KLUDGE comments
# add more general input script parsing, see NOTE comments
#   including another file
#   variable substitution could be used to create new command (hard)
#   single/double quotes could affect args and hide nested commands
#   hybrid styles can have required sub-styles
#   if-then-else and run every could nest other commands as args
# special case: linalg sometimes needed with ATC external lib
#               could scan atc/Makefile.lammps for it
# can long makes print incrementally to screen as they progress?
# add a verbose output option?

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------

import sys,glob,re,os,commands

# -----------------------------------------------------------------------
# variables

# support files in lammps/src neeeded for a build

support = ["Makefile","Make.sh","Makefile.package.empty",
           "Makefile.package.settings.empty"]

# packages that have external libs with their external lib dir

extlibs = {"USER-ATC": "atc", "USER-AWPMD": "awpmd", "USER-COLVARS": "colvars",
           "USER-CUDA": "cuda","GPU": "gpu","MEAM": "meam", "POEMS": "poems",
           "REAX": "reax"}

# help messages

syntax_message = """
Make.py switch args ...
Switches: -i file1 file2 ...
          -p package1 package2 ...
          -u package1 package2 ...
          -e package1 arg1 arg2 package2 ...
          -o dir
          -b machine
          -s suffix1 suffix2 ...
          -l dir
          -j N
          -h switch1 switch2 ...
Help: Make.py or Make.py -h
"""

ihelp = """
-i file1 file2 ...
scan input scripts, make list of needed files and packages
"""

phelp = """
-p package1 package2 ...
specify zero or more packages
can use "all" or "standard" or "user"
install packages required by -i input scripts
also install listed packages
"""

uhelp = """
-u package1 package2 ...
specify zero or more packages
can specify "all" or "standard" or "user", or package w/out user prefix
uninstall packages not required by -i input scripts
also uninstall listed packages
"""

ehelp = """
-e package1 arg1 arg2 package2 ...
build external libs for zero or more packages
can specify package w/out user prefix
each package can have 1 or 2 optional args
  arg1 = suffix for lib Makefile
  arg2 = suffix for lib Makefile.lammps
example packages with external libs: "gpu", "meam", "user-cuda"
build external libs for packages required by -i input scripts
also build external libs for installed packages if none listed
also build external libs for listed packages
if arg1 specified, use Makefile.arg1 to build lib, else use Makefile
if arg2 specified, copy Makefile.lammps.arg2 to Makefile.lammps,
  else copy Makefile.lammps.arg1 to Makefile.lammps, else do nothing
"""

ohelp = """
-o newdir
create newdir with exactly the needed src files for -i input scripts
newdir will be at same level as src so LAMMPS build will work
requires use of -i, cannot use with -p or -u
newdir will include STUBS and MAKE, but no package dirs
"""

bhelp = """
-b machine
build LAMMPS in src or newdir using Makefile.machine
assumes external package libs are already built
if no -o, build in src
if -o, build in newdir, including STUBS build if needed
"""

shelp = """
-s suffix1 suffix2 ...
specify optimization suffixes to add to style names in -i input scripts
requires use of -i
same suffixes as LAMMPS -s command-line switch allows
example suffixes = "opt", "gpu", "cuda", "omp"
"""

lhelp = """
-l dir
specify LAMMPS home dir
required if not running Make.py from within lammps/src
"""

jhelp = """
-j N
perform LAMMPS and external lib builds in parallel
adds "-j N" switch to "make" commands
N = 0 performs serial make
default N = 16
"""

hhelp = """
Make.py lists all switches
Make.py -h -i -b ... gives help on each switch
Make.py -h gives this message
"""

# -----------------------------------------------------------------------
# functions

# print syntax message

def syntax():
  print syntax_message
  sys.exit()

# help messages

def helptxt(switches):
  if not switches: switches = ["-h"]
  for sw in switches:
    if sw == "-i": print ihelp
    elif sw == "-p": print phelp
    elif sw == "-u": print uhelp
    elif sw == "-e": print ehelp
    elif sw == "-o": print ohelp
    elif sw == "-b": print bhelp
    elif sw == "-s": print shelp
    elif sw == "-l": print lhelp
    elif sw == "-j": print jhelp
    elif sw == "-h": print hhelp
    else: error("Invalid help switch")

# print error message

def error(str):
  print "ERROR:",str
  sys.exit()
  
# convert *.cpp file to *.h file

def c2h(file):
  return file[:-4] + ".h"
  
# convert *.h file to *.cpp file

def h2c(file):
  return file[:-2] + ".cpp"
  
# append file to list
# but only if file is not already in list
# but only if file is in cfiles or hfiles

def addfile(file,list):
  if file in list: return
  if file not in cfiles and file not in hfiles: return
  list.append(file)

# grab list of LAMMPS package dirs from src/Makefile
  
def lammps_packages():
  txt = open(srcdir + "/Makefile",'r').read()
  pattern = "PACKAGE = (.*?)\n\n"
  match = re.search(pattern,txt,re.DOTALL)
  packstd = match.group(1).replace("\\","").upper().split()
  pattern = "PACKUSER = (.*?)\n\n"
  match = re.search(pattern,txt,re.DOTALL)
  packuser = match.group(1).replace("\\","").upper().split()
  return packstd,packuser

# uninstall package in LAMMPS src dir
# package can be upper-case or lower-case

def uninstall_package(package):
  package = package.lower()
  cmd = "cd %s; make no-%s" % (srcdir,package)
  commands.getoutput(cmd)
  print "  uninstalled",package
  
# install package in LAMMPS src dir
# package can be upper-case or lower-case

def install_package(package):
  package = package.lower()
  cmd = "cd %s; make yes-%s" % (srcdir,package)
  commands.getoutput(cmd)
  print "  installed",package

# re-parse list of packages input via -e switch
# one or two flags may be appended to a package
# convert into dict with value = list of flags

def reparse_packages(packages):
  dict = {}
  last = ""
  for package in packages:
    upper = package.upper()
    flag = 0
    if upper in packstd: flag = 1
    elif upper in packuser: flag = 1
    elif "USER-" + upper in packuser: flag = 1
    if flag:
      dict[package] = []
      last = package
    else:
      if not last: error("Invalid package in -e list")
      buildflags = dict[last]
      if len(buildflags) == 2:
        error("Too many machine flags appended to package %s" % last)
      buildflags.append(package)
  return dict

# detect which packages are installed via "make package-status" command
# return as list

def detect_installed_packages():
  cmd = "cd %s; make package-status" % srcdir
  output = commands.getoutput(cmd)
  pattern = "Installed\s*(.*): package (.*)"
  status = re.findall(pattern,output)
  installed = []
  for entry in status:
    if entry[0] == "YES": installed.append(entry[1])
  return installed
  
# create dict of all LAMMPS *.cpp and *.h files with no duplicates
# key = file
# value = package dir if file in a package dir
# value = empty if file is only in lammps/src
# discard style_*.h files

def lammps_files():
  files = []
  for dir in packstd:
    files += glob.glob(dir + "/*.cpp")
    files += glob.glob(dir + "/*.h")
  for dir in packuser:
    files += glob.glob(dir + "/*.cpp")
    files += glob.glob(dir + "/*.h")
  files += glob.glob("*.cpp")
  files += glob.glob("*.h")

  cfiles = {}
  hfiles = {}
  for file in files:
    dir = os.path.dirname(file)
    file = os.path.basename(file)
    if file.find("style_") == 0: continue
    if file[-2:] == ".h":
      if file not in hfiles: hfiles[file] = dir
    else:
      if file not in cfiles: cfiles[file] = dir
      
  return cfiles,hfiles

# parse each cfile,hfile to find list of LAMMPS *.h files it includes
# create cdep,hdep = dict with key = filename, value = list of *.h files
# KLUDGE for accelerator_cuda.h to ignore the *.h files it includes
#        when "cuda" is not explicitly listed as a suffix,
#        this is b/c it has an ifdef that includes several *cuda*.h files

def dependencies():
  pattern = re.compile('^#include "(.+?)"',re.MULTILINE)

  cdep = {}
  for file in cfiles:
    if cfiles[file]: file = cfiles[file] + "/" + file
    txt = open(file,'r').read()
    possible = re.findall(pattern,txt)
    depends = []
    for dfile in possible:
      if dfile in hfiles: depends.append(dfile)
    cdep[os.path.basename(file)] = depends

  hdep = {}
  for file in hfiles:
    if hfiles[file]: file = hfiles[file] + "/" + file
    txt = open(file,'r').read()
    possible = re.findall(pattern,txt)
    depends = []
    for dfile in possible:
      if dfile in hfiles: depends.append(dfile)
    if file == "accelerator_cuda.h" and "cuda" not in suffixes: depends = []
    hdep[os.path.basename(file)] = depends

  return cdep,hdep

# add style files referenced in input script to lookup list
# parse input script lines to look for them
# add suffix for certain styles

def need_from_input_script(file,lookup):
  lines = open(file,'r').readlines()

  # skip comment and blank lines
  # NOTE: need to treat concatenated lines
  # NOTE: need to handle quoted args when split into words
  # NOTE: need to treat if and run every and print, which nest commands

  cmds = {}
  for line in lines:
    if line.find("#") >= 0: line = line[:line.find("#")]
    line = line.strip()
    if len(line) == 0: continue
    words = line.split()
    cmds[words[0]] = words[1:]

  # for commands with styles, generate file associated with style name
  # suffixes induce multiple file variants
  # final else handles command styles if "CommandStyle" found in *.h file
    
  for cmd in cmds:
    args = cmds[cmd]
    files = []
    if cmd == "atom_style":
      files.append("atom_vec_" + args[0])
      for suffix in suffixes: files.append("atom_vec_" + args[0] + "/" + suffix)
    elif cmd == "pair_style":
      files.append("pair_" + args[0])
      for suffix in suffixes: files.append("pair_" + args[0] + "/" + suffix)
    elif cmd == "bond_style":
      files.append("bond_" + args[0])
      for suffix in suffixes: files.append("bond_" + args[0] + "/" + suffix)
    elif cmd == "angle_style":
      files.append("angle_" + args[0])
      for suffix in suffixes: files.append("angle_" + args[0] + "/" + suffix)
    elif cmd == "dihedral_style":
      files.append("dihedral_" + args[0])
      for suffix in suffixes: files.append("dihedral_" + args[0] + "/" + suffix)
    elif cmd == "improper_style":
      files.append("improper_" + args[0])
      for suffix in suffixes: files.append("improper_" + args[0] + "/" + suffix)
    elif cmd == "kspace_style":
      files.append(args[0])
      for suffix in suffixes: files.append(args[0] + "/" + suffix)
    elif cmd == "run_style":
      files.append(args[0])
      for suffix in suffixes: files.append(args[0] + "/" + suffix)
    elif cmd == "min_style":
      files.append("min_" + args[0])
    elif cmd == "compute":
      files.append("compute_" + args[2])
      for suffix in suffixes: files.append("compute_" + args[2] + "/" + suffix)
    elif cmd == "dump":
      files.append("dump_" + args[2])
    elif cmd == "fix":
      files.append("fix_" + args[2])
      for suffix in suffixes: files.append("fix_" + args[2] + "/" + suffix)
    elif cmd == "region":
      files.append("region_" + args[1])
    else:
      tmpfile = cmd + ".cpp"
      if tmpfile not in cfiles: continue
      if cfiles[tmpfile]: tmpfile = cfiles[tmpfile] + "/" + tmpfile
      txt = open(c2h(tmpfile)).read()
      if not "CommandStyle" in txt: continue
      files.append(cmd)

    for file in files:
      file = file.replace("/","_") + ".cpp"
      addfile(file,lookup)

# add additional CPP files required to build LAMMPS to lookup list
# always add main.cpp
# always add atom_vec_atomic,
#   since Atom class creates it directly
# always add compute_temp, compute_pe, compute_pressure,
#   since Output class creates these computes directly
# always add verlet and min_cg,
#   since Update class creates them directly
# always add src/neigh_*.cpp since are part of Neighbor class,
#   but only neighbor.cpp is inferred by neighbor.h
# KLUDGE:
#   likewise add USER-OMP/neigh_*.cpp if any USER-OMP files are in lookup
#   likewise add USER-CUDA/neigh_*.cpp if any USER-CUDA files are in lookup
      
def need_additional(lookup):
  addfile("main.cpp",lookup)
  addfile("atom_vec_atomic.cpp",lookup)
  addfile("compute_temp.cpp",lookup)
  addfile("compute_pe.cpp",lookup)
  addfile("compute_pressure.cpp",lookup)
  addfile("verlet.cpp",lookup)
  addfile("min_cg.cpp",lookup)

  user_cuda = user_omp = 0
  for file in lookup:
    if file[-2:] == ".h":
      if hfiles[file] == "USER-CUDA": user_cuda = 1
      if hfiles[file] == "USER-OMP": user_omp = 1
    else:
      if cfiles[file] == "USER-CUDA": user_cuda = 1
      if cfiles[file] == "USER-OMP": user_omp = 1
    
  for file in cfiles:
    if file.find("neigh_") == 0:
      if cfiles[file] == "": addfile(file,lookup)
      if cfiles[file] == "USER-CUDA" and user_cuda: addfile(file,lookup)
      if cfiles[file] == "USER-OMP" and user_omp: addfile(file,lookup)

# creq,hreq = list of all CPP and H files needed to build LAMMPS
# inferred from lookup list and each file's dependencies in cdep,hdep
# when an H file is added, also add corresponding CPP file to lookup
# addfile insures no duplicates in creq,hreq

def resolve_dependencies(lookup):
  creq = []
  hreq = []

  i = 0
  while i < len(lookup):
    file = lookup[i]
    if file[-2:] == ".h":
      hreq.append(file)
      dfiles = hdep[file]
      addfile(h2c(file),lookup)
    else:
      creq.append(file)
      dfiles = cdep[file]
    for dfile in dfiles: addfile(dfile,lookup)
    i += 1
    
  return creq,hreq

# -----------------------------------------------------------------------
# main program

# defaults

infiles = []
packflag = 0
packages = []
unpackflag = 0
unpackages = []
extflag = 0
extpackages = []
newdir = ""
build = ""
suffixes = []
lammpsdir = ""
jmake = 16
helpflag = 0
help = []

# parse command-line args

args = sys.argv
narg = len(args)

if narg == 1: syntax()

iarg = 1
while iarg < narg:
  if args[iarg] == "-i":
    if iarg+2 > narg: syntax()
    iarg += 1
    while iarg < narg and args[iarg][0] != '-':
      infiles.append(args[iarg])
      iarg += 1
  elif args[iarg] == "-p":
    packflag = 1
    iarg += 1
    while iarg < narg and args[iarg][0] != '-':
      packages.append(args[iarg])
      iarg += 1
  elif args[iarg] == "-u":
    unpackflag = 1
    iarg += 1
    while iarg < narg and args[iarg][0] != '-':
      unpackages.append(args[iarg])
      iarg += 1
  elif args[iarg] == "-e":
    extflag = 1
    iarg += 1
    while iarg < narg and args[iarg][0] != '-':
      extpackages.append(args[iarg])
      iarg += 1
  elif args[iarg] == "-o":
    if iarg+2 > narg: syntax()
    newdir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    if iarg+2 > narg: syntax()
    build = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-s":
    if iarg+2 > narg: syntax()
    iarg += 1
    while iarg < narg and args[iarg][0] != '-':
      suffixes.append(args[iarg])
      iarg += 1
  elif args[iarg] == "-l":
    if iarg+2 > narg: syntax()
    lammpsdir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-j":
    if iarg+2 > narg: syntax()
    jmake = int(args[iarg+1])
    iarg += 1
  elif args[iarg] == "-h":
    helpflag = 1
    iarg += 1
    while iarg < narg and args[iarg][0] == '-':
      help.append(args[iarg])
      iarg += 1
  else: syntax()

# help

if helpflag: helptxt(help)
  
# error check

if newdir:
  if not infiles: error("Cannot use -o without -i")
  if packflag: error("Cannot use -o with -p")
  if unpackages: error("Cannot use -o with -u")

if suffixes:
  if not infiles: error("Cannot use -s without -i")

# setup strings
# srcdir = LAMMPS src dir
# check that srcdir is valid via 2 "empty" files and version.h
# libdir = LAMMPS lib dir
# newdir = new src dir
# makestr = -j N for parallel make
  
if lammpsdir: srcdir = lammpsdir + "/src"
else: srcdir = "."

if not os.path.isfile(srcdir + "/Makefile.package.empty") \
      or not os.path.isfile(srcdir + "/Makefile.package.settings.empty") \
      or not os.path.isfile(srcdir + "/version.h"):
  error("Not running in LAMMPS src dir, use -l switch")

if srcdir == ".": libdir = "../lib"
else: libdir = srcdir[:-3] + "lib"

if newdir:
  if srcdir == ".": newdir = "../" + newdir
  else: newdir = srcdir[:-3] + newdir

if jmake == 0: makestr = ""
else: makestr = "-j %d" % jmake
  
# packstd,packuser = list of package dirs in lammps/src

packstd,packuser = lammps_packages()

# -i switch
# cfiles,hfiles = all LAMMPS src files
# cdep,hdep = src files with list of *.h files each depends on
# lookup = needed *.cpp files = default + input script styles
# creq,hreq = files needed to build LAMMPS = lookup + dependencies
# packneed = needed packages

if infiles:
  print "Scanning input scripts ..."
  cfiles,hfiles = lammps_files()
  cdep,hdep = dependencies()
  lookup = []
  for file in infiles: need_from_input_script(file,lookup)
  need_additional(lookup)
  creq,hreq = resolve_dependencies(lookup)
  
  packneed = []
  for file in creq:
    if cfiles[file] and cfiles[file] not in packneed:
      packneed.append(cfiles[file])
  for file in hreq:
    if hfiles[file] and hfiles[file] not in packneed:
      packneed.append(hfiles[file])

  print "  scripts require %d LAMMPS *.cpp files" % len(creq)
  print "  scripts require %d LAMMPS *.h files" % len(hreq)
  print "  scripts require %d packages:" % len(packneed),
  print packneed

# -u switch
# if -i, add unneeded packages to list
# add explicitly named packages to list
# uninstall packages in list only if installed

if unpackflag:
  print "Uninstalling packages ..."
  
  list = []
  if infiles:
    for package in packstd:
      if package not in packneed: list.append(package)
    for package in packuser:
      if package not in packneed: list.append(package)
  for entry in unpackages:
    entry = entry.upper()
    if entry == "ALL":
      for package in packstd:
        if package not in list: list.append(package)
      for package in packuser:
        if package not in list: list.append(package)
    elif entry == "STANDARD":
      for package in packstd:
        if package not in list: list.append(package)
    elif entry == "USER":
      for package in packuser:
        if package not in list: list.append(package)
    elif entry in packstd:
      if entry not in list: list.append(entry)
    elif entry in packuser:
      if entry not in list: list.append(entry)
    elif "USER-" + entry in packuser:
      if "USER-" + entry not in list: list.append("USER-" + entry)
    else: error("Invalid package %s to uninstall" % entry.lower())

  installed_packages = detect_installed_packages()
  for package in list:
    if package in installed_packages: uninstall_package(package)
    
# -p switch
# if -i, add needed packages to list
# add explicitly named packages to list
# install packages in list only if uninstalled

if packflag:
  print "Installing packages ..."
  
  list = []
  if infiles:
    for package in packneed: list.append(package)
  for entry in packages:
    entry = entry.upper()
    if entry == "ALL":
      for package in packstd:
        if package not in list: list.append(package)
      for package in packuser:
        if package not in list: list.append(package)
    elif entry == "STANDARD":
      for package in packstd:
        if package not in list: list.append(package)
    elif entry == "USER":
      for package in packuser:
        if package not in list: list.append(package)
    elif entry in packstd:
      if entry not in list: list.append(entry)
    elif entry in packuser:
      if entry not in list: list.append(entry)
    elif "USER-" + entry in packuser:
      if "USER-" + entry not in list: list.append("USER-" + entry)
    else: error("Invalid package %s to install" % entry.lower())

  installed_packages = detect_installed_packages()
  for package in list:
    if package not in installed_packages: install_package(package)

# -e switch
# add explicitly named libs to list first, so can include arg1/arg2
# if -i, add needed libs to list
# if none named, add libs for installed packages to list
# use Makefile.arg1 and Makefile.lammps.arg1/arg2 if specified
# force rebuild by doing clean first
    
if extflag:
  print "Building external libraries ..."

  list = {}
  packages_with_flags = reparse_packages(extpackages)
  for package in packages_with_flags:
    upper = package.upper()
    if upper in extlibs and upper not in list:
      list[upper] = packages_with_flags[package]
    upper = "USER-" + upper
    if upper in extlibs and upper not in list:
      list[upper] = packages_with_flags[package]

  if infiles:
    for package in packneed:
      if package in extlibs and package not in list: list[package] = []
  if not extpackages:
    installed_packages = detect_installed_packages()
    for package in installed_packages:
      if package in extlibs and package not in list: list[package] = []
    
  for package in list:
    lib = extlibs[package]
    print "  lib/%s/lib%s.a for package %s ..." % (lib,lib,package)
    args = list[package]

    # makefile = Makefile by default unless arg1 specified
    
    makefile = "Makefile"
    if args: makefile = "Makefile.%s" % args[0]
    if not os.path.isfile("%s/%s/%s" % (libdir,lib,makefile)):
      error("%s for external lib %s does not exist" % (makefile,lib))

    # no makefile_lammps by default unless arg2 or arg1 specified
      
    makefile_lammps = ""
    if args and len(args) == 2:
      makefile_lammps = "Makefile.lammps.%s" % args[1]
      if not os.path.isfile("%s/%s/%s" % (libdir,lib,makefile_lammps)):
        error("%s for external lib %s does not exist" % (makefile,lib))
    elif args:
      makefile_lammps = "Makefile.lammps.%s" % args[0]

    # target lib file
    # KLUDGE for cuda, since it is liblammpscuda.a instead of libcuda.a
      
    libfile = "%s/%s/lib%s.a" % (libdir,lib,lib)
    if lib == "cuda": libfile = "%s/%s/liblammpscuda.a" % (libdir,lib)
    
    # remove lib file

    if os.path.exists(libfile): os.remove(libfile)
    
    # make clean

    print "    cleaning ..."
    cmd = "cd %s/%s; make -f %s clean" % (libdir,lib,makefile)
    output = commands.getoutput(cmd)

    # build the external lib and check if successful
      
    print "    building ..."
    cmd = "cd %s/%s; make %s -f %s" % (libdir,lib,makestr,makefile)
    output = commands.getoutput(cmd)
    if os.path.exists(libfile):
      print "  successfully built lib/%s/lib%s.a" % (lib,lib)
    else:
      print output
      error("Build of lib/%s/lib%s.a was unsuccessful" % (lib,lib))

    # if it exists, copy makefile_lammps to Makefile.lammps

    if makefile_lammps:
      if os.path.isfile("%s/%s/%s" % (libdir,lib,makefile_lammps)):
        cmd = "cd %s/%s; cp %s Makefile.lammps" % (libdir,lib,makefile_lammps) 
        commands.getoutput(cmd)

# -n switch
# create newdir with exactly the needed src files for scripts
# newdir has support files and STUBS and MAKE, but no package dirs
# insure Makefile.package and Makefile.package.settings are correct by
#   installing all needed packages in newdir, then delete package dirs

if newdir:
  print "Creating %s ..." % newdir
  if os.path.exists(newdir): error("Directory %s already exists" % newdir)
  os.mkdir(newdir)
  os.mkdir(newdir + "/MAKE")
  os.mkdir(newdir + "/STUBS")
    
  for file in support:
    commands.getoutput("cp %s/%s %s" % (srcdir,file,newdir))

  files = glob.glob("%s/MAKE/Makefile.*" % srcdir)
  fstr = " ".join(files)
  commands.getoutput("cp %s %s/MAKE" % (fstr,newdir))

  commands.getoutput("cp %s/STUBS/Makefile %s/STUBS" % (srcdir,newdir))
  commands.getoutput("cp %s/STUBS/mpi.c %s/STUBS" % (srcdir,newdir))
  commands.getoutput("cp %s/STUBS/mpi.h %s/STUBS" % (srcdir,newdir))

  for package in packneed:
    os.mkdir("%s/%s" % (newdir,package))
    commands.getoutput("cp %s/%s/* %s/%s" % (srcdir,package,newdir,package))
    commands.getoutput("cd %s; make yes-%s" % (newdir,package.lower()))
    commands.getoutput("rm -r %s/%s" % (newdir,package))

  commands.getoutput("rm -r %s/*.cpp" % (newdir))
  commands.getoutput("rm -r %s/*.h" % (newdir))
    
  for file in creq:
    if cfiles[file]: file = cfiles[file] + "/" + file
    commands.getoutput("cp %s/%s %s" % (srcdir,file,newdir))
  for file in hreq:
    if hfiles[file]: file = hfiles[file] + "/" + file
    commands.getoutput("cp %s/%s %s" % (srcdir,file,newdir))

# -b switch
# build LAMMPS in src or newdir
# assumes external package libs are already built
# build STUBS lib if needed and not already built
    
if build:
  print "Building LAMMPS for %s ..." % build
  if newdir: builddir = newdir
  else: builddir = srcdir

  makefile = "%s/MAKE/Makefile.%s" % (builddir,build)
  if not os.path.exists(makefile):
    error("%s does not exist" % makefile)

  txt = open(makefile,"r").read()
  if "../STUBS" in txt:
    if not os.path.exists("%s/STUBS/libmpi.a" % builddir):
      cmd = "cd %s/STUBS; make" % builddir
      output = commands.getoutput(cmd)
      if os.path.exists("%s/STUBS/libmpi.a" % builddir):
        print "Successfully built STUBS/libmpi.a"
      else:
        print output
        error("Build of STUBS/libmpi.a was unsuccessful")
    
  if os.path.exists("%s/lmp_%s" % (builddir,build)):
    os.remove("%s/lmp_%s" % (builddir,build))
  cmd = "cd %s; make %s %s" % (builddir,makestr,build)
  output = commands.getoutput(cmd)
  if os.path.exists("%s/lmp_%s" % (builddir,build)):
    print "Successfully built %s/lmp_%s" % (builddir,build)
  else:
    print output
    error("Build of %s/lmp_%s was unsuccessful" % (builddir,build))
