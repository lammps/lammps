#!/usr/bin/python

# Make.py tool for building LAMMPS and its package libs

import sys,os,commands,re,copy

# switch abbrevs
# switch classes = created class for each switch
# lib classes = all auxiliary libs in LAMMPS plus "all"
# extra classes = packages that need extra build options
# setargs = allowed settings
# actionargs = allowed actions (also lib-dir and machine)

abbrevs = "adhjmoprsv"

switchclasses = ("actions","dir","help","jmake","makefile",
                 "output","packages","redo","settings","verbose")
libclasses = ("atc","awpmd","colvars","cuda","gpu",
              "meam","poems","qmmm","reax")
extraclasses = ("intel","kokkos")

setargs = ("gzip","#gzip","ffmpeg","#ffmpeg","smallbig","bigbig","smallsmall")
actionargs = ("lib-all","file","clean","exe")

# ----------------------------------------------------------------
# functions
# ----------------------------------------------------------------

# print error str and exit

def error(str):
  print "ERROR:",str
  sys.exit()

# store command-line args as sw = dict of key/value
# key = switch letter, value = list of following args
# order = list of switches in order specified
# enforce no switch more than once
# once any arg is an action, store remaining args as -a switch
# can specify explicit -a
  
def parse_args(args):
  narg = len(args)
  sw = {}
  order = []
  iarg = 0
  while iarg < narg:
    if args[iarg][0] != '-':
      switch = 'a'
      first = iarg
    else:
      switch = args[iarg][1:]
      first = iarg+1
    if switch in sw: error("Duplicate switch %s" % args[iarg])
    order.append(switch)
    if switch == 'a':
      sw[switch] = args[first:]
      break
    last = first+1
    while last < narg and args[last][0] != '-' and \
          args[last] not in actionargs and \
          not args[last].startswith("lib-"):
      last += 1
    sw[switch] = args[first:last]
    iarg = last
    
  return sw,order

# convert switches in sw back to a string, in switch_order
# just append action args

def switch2str(switches,switch_order):
  txt = ""
  for switch in switch_order:
    if txt: txt += ' '
    if switch == 'a': txt += ' '.join(switches[switch])
    else:
      txt += "-%s" % switch
      if switches[switch]: txt += ' ' + ' '.join(switches[switch])
  return txt

# ----------------------------------------------------------------
# classes, one per switch
# ----------------------------------------------------------------

# actions

class Actions:
  def __init__(self,list):
    self.inlist = list[:]
    
  def help(self):
    return """
Actions:
  possible actions: lib-all, lib-dir, file, clean, exe or machine
  can specify zero or more actions in any order
    except machine must be last (if used)
    each action can appear no more than once
    if switches used and machine is only action, prefix by "-a" switch
  actions happen in order below, indpendent of specified order
  some actions depend on installed packages
    installed packages = currently installed + result of -p switch
  lib-all or lib-dir = build auxiliary libraries
    lib-all builds all auxiliary libs needed by installed packages
    lib-dir builds a specific lib whether package installed or not
      dir is any dir in lib directory (atc, cuda, meam, etc) except linalg
      can be specified multiple times for different dirs
  file = create src/MAKE/MINE/Makefile.auto
    use -m switch for Makefile.machine to start from
      else use existing Makefile.auto
    adds settings needed for installed accelerator packages
  clean = invoke "make clean-auto", insures full build
    useful if compiler flags have changed
  exe or machine = build LAMMPS
    machine can be any existing Makefile.machine suffix
      machine is simply converted to exe, as well as:
        "-m machine" added if -m switch not specified
        "-o machine" added if -o switch not specified
        if either "-m"  or "-o" are specified, they are not overridden
    exe builds using Makefile.auto
    if no file action, first generates a src/MAKE/MINE/Makefile.auto
      use -m switch to make copy of existing Makefile.machine
        or Makefile.auto must already exist
      unlike file action, this does not change Makefile.auto
    does not invoke and lib actions, since libs could be previously built
    produces src/lmp_auto or error message if unsuccessful
"""
  
  def check(self):
    alist = []
    nlib = 0
    for i,one in enumerate(self.inlist):
      if one in alist: error("An action is duplicated")
      if one.startswith("lib-"):
        lib = one[4:]
        if lib != "all" and lib not in libclasses: error("Actions are invalid")
        alist.insert(nlib,one)
        nlib += 1
      elif one == "file":
        if nlib == 0: alist.insert(0,"file")
        else: alist.insert(1,"file")
      elif one == "clean":
        if nlib == 0: alist.insert(0,"clean")
        elif "file" not in alist: alist.insert(1,"clean")
        else: alist.insert(2,"clean")
      elif one == "exe": alist.append("exe")
      # allow last action to be unknown in case is a machine (checked in setup)
      elif i == len(self.inlist)-1: alist.append(one)
      else: error("Actions are invalid")
    self.alist = alist

  # if last action is unknown, assume machine and convert to exe
  # only done if action is a suffix for an existing Makefile.machine
  # return machine if conversion done, else None
    
  def setup(self):
    machine = self.alist[-1]
    if machine in actionargs or machine.startswith("lib-"): return None
    make = MakeReader(machine,2)
    self.alist[-1] = "exe"
    return machine

  # build libraries needed in installed packages
  
  def lib(self,suffix):
    if suffix != "all":
      print "  building library",suffix
      exec("%s.build()" % suffix)
    else:
      final = packages.final
      for one in packages.lib:
        if final[one]:
          if "user" in one: pkg = one[5:]
          else: pkg = one
          print "  building library",pkg
          str = "%s.build()" % pkg
          exec(str)

  # read existing Makefile.machine
  # if tweak = 1: tweak file if accelerator packages will be installed
  # if tweak = 0: no change to file
  # write out new Makefile.auto
  
  def file(self,tweak):
    if makefile: machine = makefile.machine
    else: machine = "auto"
    make = MakeReader(machine,1)
    
    if tweak:
      final = packages.final
      if final["user-omp"]:
        make.addvar("CCFLAGS","-fopenmp")
        make.addvar("CCFLAGS","-restrict")
        make.addvar("LINKFLAGS","-fopenmp")

      if final["user-intel"]:
        if intel.mode == "cpu":
          make.addvar("CCFLAGS","-fopenmp")
          make.addvar("CCFLAGS","-DLAMMPS_MEMALIGN=64")
          make.addvar("CCFLAGS","-restrict")
          make.addvar("CCFLAGS","-xHost")
          make.addvar("CCFLAGS","-fno-alias")
          make.addvar("CCFLAGS","-ansi-alias")
          make.addvar("CCFLAGS","-override-limits")
          make.addvar("LINKFLAGS","-fopenmp")
          make.delvar("CCFLAGS","-DLMP_INTEL_OFFLOAD")
          make.delvar("LINKFLAGS","-offload")
        elif intel.mode == "phi":
          make.addvar("CCFLAGS","-fopenmp")
          make.addvar("CCFLAGS","-DLAMMPS_MEMALIGN=64")
          make.addvar("CCFLAGS","-restrict")
          make.addvar("CCFLAGS","-xHost")
          make.addvar("CCFLAGS","-DLMP_INTEL_OFFLOAD")
          make.addvar("CCFLAGS","-fno-alias")
          make.addvar("CCFLAGS","-ansi-alias")
          make.addvar("CCFLAGS","-override-limits")
          make.addvar("CCFLAGS",'-offload-option,mic,compiler,' +
                      '"-fp-model fast=2 -mGLOB_default_function_attrs=' +
                      '\\"gather_scatter_loop_unroll=4\\""')
          make.addvar("LINKFLAGS","-fopenmp")
          make.addvar("LINKFLAGS","-offload")
        else: error("Specified user-intel package w/out cpu or phi suffix")

      if final["kokkos"]:
        if kokkos.mode == "omp":
          make.addvar("OMP","yes")
          make.delvar("CUDA")
          make.delvar("MIC")
        elif kokkos.mode == "cuda":
          make.addvar("OMP","yes")
          make.addvar("CUDA","yes")
          make.delvar("MIC")
          if kokkos.archflag:
            make.delvar("CCFLAGS","-arch=sm_*")
            make.addvar("CCFLAGS","-arch=sm_%s" % kokkos.arch)
        elif kokkos.mode == "phi":
          make.addvar("OMP","yes")
          make.addvar("MIC","yes")
          make.delvar("CUDA")
        else: error("Specified kokkos package w/out omp or cuda or phi suffix")

      if settings:
        list = settings.inlist
        for one in list:
          if one == "gzip": make.addvar("LMP_INC","-DLAMMPS_GZIP")
          elif one == "#gzip": make.delvar("LMP_INC","-DLAMMPS_GZIP")
          elif one == "ffmpeg": make.addvar("LMP_INC","-DLAMMPS_FFMPEG")
          elif one == "#ffmpeg": make.delvar("LMP_INC","-DLAMMPS_FFMPEG")
          elif one == "smallbig":
            make.delvar("LMP_INC","-DLAMMPS_BIGBIG")
            make.delvar("LMP_INC","-DLAMMPS_SMALLSMALL")
          elif one == "bigbig":
            make.delvar("LMP_INC","-DLAMMPS_SMALLBIG")
            make.delvar("LMP_INC","-DLAMMPS_SMALLSMALL")
            make.addvar("LMP_INC","-DLAMMPS_BIGBIG")
          elif one == "smallsmall":
            make.delvar("LMP_INC","-DLAMMPS_SMALLBIG")
            make.delvar("LMP_INC","-DLAMMPS_BIGBIG")
            make.addvar("LMP_INC","-DLAMMPS_SMALLSMALL")
          
    make.write("%s/MAKE/MINE/Makefile.auto" % dir.src,1)
    print "Created src/MAKE/MINE/Makefile.auto"

  # invoke "make clean-auto" to force clean before build
    
  def clean(self):
    str = "cd %s; make clean-auto" % dir.src
    commands.getoutput(str)
    if verbose: print "Performed make clean-auto"

  # build LAMMPS using Makefile.auto and -j setting
  # delete lmp_auto first, so can detect if build fails
    
  def exe(self):
    if "file" not in self.alist: self.file(0)
    commands.getoutput("cd %s; rm -f lmp_auto" % dir.src)
    if jmake: str = "cd %s; make -j %d auto" % (dir.src,jmake.n)
    else: str = "cd %s; make auto" % dir.src
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/lmp_auto" % dir.src):
      if not verbose: print txt
      error('Unsuccessful "make auto"')
    else: print "Created src/lmp_auto"
    
# dir switch

class Dir:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
      
  def help(self):
    return """
-d dir
  dir = LAMMPS home dir
  if -d not specified, working dir must be lammps/src
"""

  def check(self):
    if self.inlist != None and len(self.inlist) != 1:
      error("-d args are invalid")
      
  # if inlist = None, check that cwd = lammps/src
  # store cwd and lammps dir
  # derive src,make,lib dirs from lammps dir
  # check that they all exist
  
  def setup(self):
    self.cwd = os.getcwd()
    if self.inlist == None: self.lammps = ".."
    else: self.lammps = self.inlist[0]
    self.lammps = os.path.realpath(self.lammps)
    self.src = self.lammps + "/src"
    self.make = self.lammps + "/src/MAKE"
    self.lib = self.lammps + "/lib"
    if not os.path.isdir(self.lammps): error("LAMMPS home dir is invalid")
    if not os.path.isdir(self.src): error("LAMMPS src dir is invalid")
    if not os.path.isdir(self.lib): error("LAMMPS lib dir is invalid")

# help switch

class Help:
  def __init__(self,list): pass

  def help(self):
    return """
Syntax: Make.py switch args ... {action1} {action2} ...
  actions:
    lib-all, lib-dir, clean, file, exe or machine
    zero or more actions, in any order (machine must be last)
  switches:
    -d (dir), -j (jmake), -m (makefile), -o (output),
    -p (packages), -r (redo), -s (settings), -v (verbose)
  switches for libs:
    -atc, -awpmd, -colvars, -cuda
    -gpu, -meam, -poems, -qmmm, -reax
  switches for extra package options:
    -intel, -kokkos

add -h switch to command line to print this message
  and help on other specified switches or actions
"""

# jmake switch
  
class Jmake:
  def __init__(self,list):
    self.inlist = list[:]
  
  def help(self):
    return """
-j N
  use N procs for performing parallel make commands
  used when building a lib or LAMMPS itself
  if -j not specified, serial make commands run on single core
"""

  def check(self):
    if len(self.inlist) != 1: error("-j args are invalid")
    if not self.inlist[0].isdigit(): error("-j args are invalid")
    n = int(self.inlist[0])
    if n <= 0: error("-j args are invalid")
    self.n = n
        
# makefile switch

class Makefile:
  def __init__(self,list):
    self.inlist = list[:]
  
  def help(self):
    return """
-m machine
  use src/MAKE/Makefile.machine as starting point to create Makefile.auto
  if -m not specified, file/exe actions alter existing Makefile.auto  
"""
  
  def check(self):
    if len(self.inlist) != 1: error("-m args are invalid")
    self.machine = self.inlist[0]
    
# output switch

class Output:
  def __init__(self,list):
    self.inlist = list[:]

  def help(self):
    return """
-o machine
  copy final src/lmp_auto to lmp_machine in working dir
  if -o not specified, exe action only produces src/lmp_auto
"""

  def check(self):
    if len(self.inlist) != 1: error("-o args are invalid")
    self.machine = self.inlist[0]

# packages switch
  
class Packages:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]

  def help(self):
    return """
-p = package1 package2 ...
  list of packages to install or uninstall in order specified
  operates on set of packages currently installed
  valid package names:
    and LAMMPS standard or user package (type "make package" to see list)
    prefix by yes/no to install/uninstall (see abbrevs)
      yes-molecule, yes-user-atc, no-molecule, no-user-atc
  can use LAMMPS categories (type "make package" to see list)
    all = all standard and user packages
    std (or standard) = all standard packages
    user = all user packages
    lib = all standard and user packages with auxiliary libs
  can abbreviate package names and yes/no
    omp = user-omp = yes-user-omp
    #omp = #user-omp = no-user-omp
    user = yes-user, #user = no-user
    all = yes-all, #all = no-all
  when action performed, list is processed in order,
    as if typed "make yes/no" for each
  if "orig" or "original" is last package in list,
    set of installed packages will be restored to original (current) list
    after "build" action is performed
  if -p not specified, currently installed packages are not changed
"""

  def check(self):
    if self.inlist != None and not self.inlist: error("-p args are invalid")

  def setup(self):
      
    # extract package lists from src/Makefile

    make = MakeReader("%s/Makefile" % dir.src)
    std = make.getvar("PACKAGE")
    user = make.getvar("PACKUSER")
    lib = make.getvar("PACKLIB")
    all = std + user
    
    # plist = command line args expanded to yes-package or no-package
        
    plist = []
    if self.inlist:
      for one in self.inlist:
        if one in std:
          plist.append("yes-%s" % one)
        elif one in user:
          plist.append("yes-%s" % one)
        elif "user-"+one in user:
          plist.append("yes-user-%s" % one)
        elif one == '#' and one[1:] in std:
          plist.append("no-%s" % one[1:])
        elif one[0] == '#' and one[1:] in user:
          plist.append("no-%s" % one[1:])
        elif one[0] == '#' and "user-"+one[1:] in user:
          plist.append("no-user-%s" % one[1:])
        elif one == "std" or one == "standard" or one == "user" or \
              one == "lib" or one == "all": 
          plist.append("yes-%s" % one)
        elif one == "#std" or one == "#standard" or one == "#user" or \
              one == "#lib" or one == "#all": 
          plist.append("no-%s" % one[1:])
        elif one == "orig": plist.append(one)
        else: error("Invalid package name %s" % one)
      if "orig" in plist and plist.index("orig") != len(plist)-1:
        error('-p orig arg must be last')
      if plist.count("orig") > 1: error('-p orig arg must be last')

    # original = dict of all packages
    # key = package name, value = 1 if currently installed, else 0

    original = {}
    str = "cd %s; make ps" % dir.src
    output = commands.getoutput(str).split('\n')
    pattern = "Installed\s+(\w+): package (\S+)"
    for line in output:
      m = re.search(pattern,line)
      if not m: continue
      pkg = m.group(2).lower()
      if pkg not in all: error('Package list does not math "make ps" results')
      if m.group(1) == "NO": original[pkg] = 0      
      elif m.group(1) == "YES": original[pkg] = 1

    # final = dict of all packages after plist applied to original
    # key = package name, value = 1 if installed, else 0
        
    final = copy.deepcopy(original)
    for i,one in enumerate(plist):
      if "yes" in one:
        pkg = one[4:]
        yes = 1
      else:
        pkg = one[3:]
        yes = 0
      if pkg in all:
        final[pkg] = yes
      elif pkg == "std":
        for pkg in std: final[pkg] = yes
      elif pkg == "user":
        for pkg in user: final[pkg] = yes
      elif pkg == "lib":
        for pkg in lib: final[pkg] = yes
      elif pkg == "all":
        for pkg in all: final[pkg] = yes

    self.std = std
    self.user = user
    self.lib = lib
    self.all = all
    self.plist = plist
    self.original = original
    self.final = final
    
  # install packages in plist
    
  def install(self):
    if self.plist: print "Installing packages ..."
    for one in self.plist:
      if one == "orig": continue
      commands.getoutput("cd %s; make %s" % (dir.src,one))
    if self.plist and verbose:
      txt = commands.getoutput("cd %s; make ps" % dir.src)
      print "Package status after installation:"
      print txt
      
  # restore packages to original list if requested
  # order of re-install should not matter matter b/c of Depend.sh
  
  def uninstall(self):
    if not self.plist or self.plist[-1] != "orig": return
    print "Restoring packages to original state ..."
    commands.getoutput("cd %s; make no-all" % dir.src)
    for one in self.all:
      if self.original[one]:
        commands.getoutput("cd %s; make yes-%s" % (dir.src,one))
    if verbose:
      txt = commands.getoutput("cd %s; make ps" % dir.src)
      print "Restored package status:"
      print txt
      
# redo switch
    
class Redo:
  def __init__(self,list):
    self.inlist = list[:]
  
  def help(self):
    return """
-r file label1 label2 ...
  all args are optional
  redo file format:
    blank lines and lines starting with "#" are skipped
    other lines are treated as commands
    each command is a list of Make.py args, as if typed at command-line
    commands can have leading label, followed by ":"
    commands cannot contain a "-r" switch
  if no args, execute previous command from src/Make.py.last
  if one arg, execute all commands from specified file
    unlabeled or labeled commands are all executed
  if multiple args, execute only matching labeled commands from file
  if other switches are specified,
    they replace matching switches in file command(s)
    or they are added to file command(s)
  if other actions are specified,
    they are added to any actions in file command(s), without de-duplication
"""
  
  def check(self):
    if len(self.inlist) == 0:
      self.dir = 1
      self.file = "Make.py.last"
      self.labels = []
    else:
      self.dir = 0
      self.file = self.inlist[0]
      self.labels = self.inlist[1:]

  # read redo file
  # self.commands = list of commands to execute
      
  def setup(self):
    if self.dir: file = "%s/%s" % (dir.src,self.file)
    else: file = self.file
    if not os.path.isfile(file): error("Redo file %s does not exist" % file)
    lines = open(file,'r').readlines()
    
    cmdlines = []
    for line in lines:
      line = line.strip()
      if not line or line[0] == '#' : continue
      cmdlines.append(line)

    # if no labels, add all file commands to command list
    # if labels, make a dict with key = label, value = command
    #   and discard unlabeled commands
      
    dict = {}
    commands = []
    for line in cmdlines:
      words = line.split()
      if "-r" in words: error("Redo command cannot contain -r switch")
      if words[0][-1] == ':': label = words[0][:-1]
      else: label = None
      if not self.labels:
        if label: commands.append(' '.join(words[1:]))
        else: commands.append(line)
      else:
        if not label: continue
        dict[label] = ' '.join(words[1:])

    # extract labeled commands from dict and add to command list
        
    for label in self.labels:
      if label not in dict: error("Redo label not in redo file")
      commands.append(dict[label])

    self.commands = commands

# settings switch

class Settings:
  def __init__(self,list):
    self.inlist = list[:]
  
  def help(self):
    return """
-s set1 set2 ...
  possible settings = gzip smallbig bigbig smallsmall
  add each setting as LAMMPS setting to created Makefile.auto
  if -s not specified, no settings are changed in Makefile.auto
"""
  
  def check(self):
    if not self.inlist: error("-s args are invalid")
    for one in self.inlist:
      if one not in setargs: error("-s args are invalid")
  
# verbose switch

class Verbose:
  def __init__(self,list):
    self.inlist = list[:]
  
  def help(self):
    return """
-v (no arguments)
  produce verbose output as Make.py executes
  if -v not specified, minimal output is produced
"""
  
  def check(self):
    if len(self.inlist): error("-v args are invalid")

# ----------------------------------------------------------------
# classes, one per LAMMPS lib
# ----------------------------------------------------------------

# ATC lib

class ATC:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "g++"
    self.lammpsflag = 0

  def help(self):
    return """
-atc make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = g++)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-atc args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-atc args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-atc args are invalid")

  def build(self):
    libdir = dir.lib + "/atc"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)
    
    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libatc.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/atc library")
    else: print "Created lib/atc library"
    
# AWPMD lib

class AWPMD:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "mpicc"
    self.lammpsflag = 0

  def help(self):
    return """
-awpmd make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = mpicc)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-awpmd args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-awpmd args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-awpmd args are invalid")

  def build(self):
    libdir = dir.lib + "/awpmd"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libawpmd.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/awpmd library")
    else: print "Created lib/awpmd library"

# COLVARS lib

class COLVARS:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "g++"
    self.lammpsflag = 0

  def help(self):
    return """
-colvars make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = g++)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-colvars args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-colvars args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-colvars args are invalid")

  def build(self):
    libdir = dir.lib + "/colvars"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libcolvars.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/colvars library")
    else: print "Created lib/colvars library"

# CUDA lib

class CUDA:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.mode = "double"
    self.arch = "31"

  def help(self):
    return """
-cuda mode=double arch=31
  all args are optional and can be in any order
  mode = double or mixed or single (def = double)
  arch = M (def = 31)
    M = 31 for Kepler
    M = 20 for CC2.0 (GF100/110, e.g. C2050,GTX580,GTX470)
    M = 21 for CC2.1 (GF104/114,  e.g. GTX560, GTX460, GTX450)
    M = 13 for CC1.3 (GF200, e.g. C1060, GTX285)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-cuda args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-cuda args are invalid")
      if words[0] == "mode": self.mode = words[1]
      elif words[0] == "arch": self.arch = words[1]
      else: error("-cuda args are invalid")
    if self.mode != "double" and self.mode != "mixed" and \
          self.mode != "single":
      error("-cuda args are invalid")
    if not self.arch.isdigit(): error("-cuda args are invalid")
          
  def build(self): 
    libdir = dir.lib + "/cuda"
    commands.getoutput("cd %s; make clean" % libdir)
    if self.mode == "double": n = 2
    elif self.mode == "mixed": n = 3
    elif self.mode == "single": n = 1
    str = "cd %s; make prec=%d arch=%s" % (libdir,n,self.arch)
    commands.getoutput(str)
    if jmake: str = "cd %s; make -j %d" % (libdir,jmake.n)
    else: str = "cd %s; make" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/liblammpscuda.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/cuda library")
    else: print "Created lib/cuda library"

# GPU lib

class GPU:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "linux.double"
    self.lammpsflag = self.modeflag = self.archflag = 0

  def help(self):
    return """
-gpu make=suffix lammps=suffix2 mode=double arch=N
  all args are optional and can be in any order
  make = use Makefile.suffix (def = linux.double)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
  mode = double or mixed or single (def = CUDA_PREC in makefile)
  arch = 31 (Kepler) or 21 (Fermi) (def = CUDA_ARCH in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-gpu args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-gpu args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      elif words[0] == "mode":
        self.mode = words[1]
        self.modeflag = 1
      elif words[0] == "arch":
        self.arch = words[1]
        self.archflag = 1
      else: error("-gpu args are invalid")
      if self.modeflag and (self.mode != "double" and
                            self.mode != "mixed" and
                            self.mode != "single"):
        error("-gpu args are invalid")
      if self.archflag and not self.arch.isdigit():
        error("-gpu args are invalid")

  def build(self):
    libdir = dir.lib + "/gpu"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.modeflag:
      if self.mode == "double":
        make.setvar("CUDA_PRECISION","-D_DOUBLE_DOUBLE")
      elif self.mode == "mixed":
        make.setvar("CUDA_PRECISION","-D_SINGLE_DOUBLE")
      elif self.mode == "single":
        make.setvar("CUDA_PRECISION","-D_SINGLE_SINGLE")
    if self.archflag:
      make.setvar("CUDA_ARCH","-arch=sm_%s" % self.arch)
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libgpu.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/gpu library")
    else: print "Created lib/gpu library"

# MEAM lib

class MEAM:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "gfortran"
    self.lammpsflag = 0

  def help(self):
    return """
-meam make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = gfortran)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-meam args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-meam args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-meam args are invalid")

  def build(self):
    libdir = dir.lib + "/meam"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    # do not use -j for MEAM build, parallel build does not work
    str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libmeam.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/meam library")
    else: print "Created lib/meam library"

# POEMS lib

class POEMS:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "g++"
    self.lammpsflag = 0

  def help(self):
    return """
-poems make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = g++)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-poems args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-poems args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-poems args are invalid")

  def build(self):
    libdir = dir.lib + "/poems"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libpoems.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/poems library")
    else: print "Created lib/poems library"

# QMMM lib

class QMMM:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "gfortran"
    self.lammpsflag = 0

  def help(self):
    return """
-qmmm make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = gfortran)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-qmmm args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-qmmm args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-qmmm args are invalid")

  def build(self):
    libdir = dir.lib + "/qmmm"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libqmmm.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/qmmm library")
    else: print "Created lib/qmmm library"

# REAX lib

class REAX:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.make = "gfortran"
    self.lammpsflag = 0

  def help(self):
    return """
-reax make=suffix lammps=suffix2
  all args are optional and can be in any order
  make = use Makefile.suffix (def = gfortran)
  lammps = use Makefile.lammps.suffix2 (def = EXTRAMAKE in makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-reax args are invalid")
    for one in self.inlist:
      words = one.split('=')
      if len(words) != 2: error("-reax args are invalid")
      if words[0] == "make": self.make = words[1]
      elif words[0] == "lammps": 
        self.lammps = words[1]
        self.lammpsflag = 1
      else: error("-reax args are invalid")

  def build(self):
    libdir = dir.lib + "/reax"
    make = MakeReader("%s/Makefile.%s" % (libdir,self.make))
    if self.lammpsflag:
      make.setvar("EXTRAMAKE","Makefile.lammps.%s" % self.lammps)
    make.write("%s/Makefile.auto" % libdir)

    commands.getoutput("cd %s; make -f Makefile.auto clean" % libdir)
    if jmake: str = "cd %s; make -j %d -f Makefile.auto" % (libdir,jmake.n)
    else: str = "cd %s; make -f Makefile.auto" % libdir
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/libreax.a" % libdir) or \
          not os.path.isfile("%s/Makefile.lammps" % libdir):
      if not verbose: print txt
      error("Unsuccessful build of lib/reax library")
    else: print "Created lib/reax library"

# ----------------------------------------------------------------
# extra classes for intel and kokkos package options
# ----------------------------------------------------------------

# Intel class

class Intel:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.mode = "cpu"

  def help(self):
    return """
-intel mode
  mode = cpu or phi (def = cpu)
    build Intel package for CPU or Xeon Phi
"""

  def check(self):
    if self.inlist == None: return
    if len(self.inlist) != 1: error("-intel args are invalid")
    self.mode = self.inlist[0]
    if self.mode != "cpu" and self.mode != "phi":
      error("-intel args are invalid")

# Kokkos class

class Kokkos:
  def __init__(self,list):
    if list == None: self.inlist = None
    else: self.inlist = list[:]
    self.mode = "omp"
    self.archflag = 0
    
  def help(self):
    return """
-kokkos mode arch=N
  mode is not optional, arch is optional
  mode = omp or cuda or phi (def = omp if -kokkos is not used)
    build Kokkos package for omp or cuda or phi
  arch = 31 (Kepler) or 21 (Fermi) (def = -arch setting in Makefile)
"""

  def check(self):
    if self.inlist != None and len(self.inlist) == 0:
      error("-kokkos args are invalid")

    if self.inlist == None: return
    if len(self.inlist) < 1: error("-kokkos args are invalid")
    self.mode = self.inlist[0]
    if self.mode != "omp" and self.mode != "cuda" and self.mode != "phi":
      error("-kokkos args are invalid")
    for one in self.inlist[1:]:
      words = one.split('=')
      if len(words) != 2: error("-kokkos args are invalid")
      if words[0] == "arch":
        self.arch = words[1]
        self.archflag = 1
      else: error("-kokkos args are invalid")

# ----------------------------------------------------------------
# auxiliary classes
# ----------------------------------------------------------------

# read, tweak, and write a Makefile

class MakeReader:

  # read a makefile
  # flag = 0 if file is full path name
  # flag = 1,2 if file is suffix for any Makefile.machine under src/MAKE
  #   look for this file in same order that src/Makefile does
  #   if flag = 1, read the file
  #   if flag = 2, just check if file exists

  def __init__(self,file,flag=0):
    if flag == 0:
      if not os.path.isfile(file): error("Makefile %s does not exist" % file)
      lines = open(file,'r').readlines()
    else:
      mfile = "%s/MAKE/MINE/Makefile.%s" % (dir.src,file)
      if not os.path.isfile(mfile):
        mfile = "%s/MAKE/Makefile.%s" % (dir.src,file)
        if not os.path.isfile(mfile):
          mfile = "%s/MAKE/OPTIONS/Makefile.%s" % (dir.src,file)
          if not os.path.isfile(mfile):
            mfile = "%s/MAKE/MACHINES/Makefile.%s" % (dir.src,file)
            if not os.path.isfile(mfile):
              error("Makefile.%s does not exist" % file)
      if flag == 1: lines = open(mfile,'r').readlines()
      else: return

    # scan lines of makefile
    # if not a variable line, just copy to newlines
    # if a variable line, concatenate any continuation lines
    # convert variable to var dict entry: key = name, value = list of words
    #   discard any portion of value string with a comment char
    # varinfo = list of variable info: (name, name with whitespace for print)
    # add index into varinfo to newlines
    # addindex = index of LAMMPS-specific line to add KOKKOS vars before it
    
    var = {}
    varinfo = []
    newlines = []
    pattern = "(\S+\s+=\s+)(.*)"
    multiline = 0
    self.addindex = 0
    
    for line in lines:
      line = line[:-1]
      if "LAMMPS-specific settings" in line: self.addindex = len(newlines)
      if multiline:
        if '#' in line: line = line[:line.find('#')]
        morevalues = line.split()
        values = values[:-1] + morevalues
        if values[-1] != '\\':
          var[name] = values
          multiline = 0
          newlines.append(str(len(varinfo)))
          varinfo.append((name,namewhite))
        continue
      varflag = 1
      if len(line.strip()) == 0: varflag = 0
      elif line.lstrip()[0] == '#': varflag = 0
      else:
        m = re.match(pattern,line)
        if not m: varflag = 0
      if varflag:
        namewhite = m.group(1)
        name = namewhite.split()[0]
        if name in var:
          error("Makefile variable %s appears more than once" % name)
        remainder = m.group(2)
        if '#' in remainder: remainder = remainder[:remainder.find('#')]
        values = remainder.split()
        if values and values[-1] == '\\': multiline = 1
        else:
          var[name] = values
          newlines.append(str(len(varinfo)))
          varinfo.append((name,namewhite))
      else:
        newlines.append(line)

    self.var = var
    self.varinfo = varinfo
    self.lines = newlines
             
  # return list of values associated with var
  # return None if var not defined
    
  def getvar(self,var):
    if var in self.var: return self.var[var]
    else: return None

  # set var to single value
  # if var not defined, error
  
  def setvar(self,var,value):
    if var not in self.var: error("Variable %s not in makefile" % var)
    self.var[var] = [value]
    
  # add value to var
  # do not add if value already defined by var
  # if var not defined,
  #   create new variable 2 lines before LAMMPS-specific settings section
  
  def addvar(self,var,value):
    if var in self.var:
      if value not in self.var[var]: self.var[var].append(value)
    else:
      if self.addindex == 0: error("No LAMMPS-specific settings line " +
                                   "in makefile, needed to add a variable")
      self.var[var] = [value]
      varwhite = "%s =\t\t" % var
      self.lines.insert(self.addindex-2,str(len(self.varinfo)))
      self.varinfo.append((var,varwhite))
      
  # if value = None, remove entire var
  #   no need to update lines or varinfo, write() will ignore deleted vars
  # else remove value from var
  # value can have trailing '*' to remove wildcard match
  # if var or value not defined, ignore it
      
  def delvar(self,var,value=None):
    if var not in self.var: return
    if not value: del self.var[var]
    elif value and value[-1] != '*':
      if value not in self.var[var]: return
      self.var[var].remove(value)
    else:
      value = value[:-1]
      values = self.var[var]
      dellist = []
      for i,one in enumerate(values):
        if one.startswith(value): dellist.append(i)
      while dellist: values.pop(dellist.pop())
      self.var[var] = values
      
  # write stored makefile lines to file, using vars that may have been updated
  # do not write var if not in dict, since has been deleted
  # wrap var values into multiple lines if needed
  # file = 1 if this is Makefile.auto, change 1st line to use "auto"
    
  def write(self,file,flag=0):
    fp = open(file,'w')
    for i,line in enumerate(self.lines):
      if not line.isdigit():
        if flag and i == 0:
          line = "# auto = makefile auto-generated by Make.py"
        print >>fp,line
      else:
        index = int(line)
        name = self.varinfo[index][0]
        txt = self.varinfo[index][1]
        if name not in self.var: continue
        values = self.var[name]
        print >>fp,"%s%s" % (txt,' '.join(values))
  
# ----------------------------------------------------------------
# main program
# ----------------------------------------------------------------

# parse command-line args
# switches dict: key = switch letter, value = list of args
# switch_order = list of switches in order
# will possibly be merged with redo file args below
        
cmd_switches,cmd_switch_order = parse_args(sys.argv[1:])

# check for redo switch, process redo file
# redolist = list of commands to execute

redoflag = 0
redolist = []

if 'r' in cmd_switches and 'h' not in cmd_switches:
  redoflag = 1
  redo = Redo(cmd_switches['r'])
  redo.check()
  if 'd' in cmd_switches:
    dir = Dir(cmd_switches['d'])
    dir.check()
  else: dir = Dir(None)
  dir.setup()
  redo.setup()
  redolist = redo.commands
  redoindex = 0
  del dir
  del redo
  if not redolist: error("No commands to execute from redo file")

# loop over Make.py commands
# if no redo switch, loop once for command-line command
# if redo, loop over one or more commands from redo file

while 1:
      
  # if redo:
  #   parse next command from redo file
  #   use command-line switches to add/replace file command switches
  #     if actions in both, are just concatenated
  #     do not add -r, since already processed
  #       and don't want -r swtich to appear in Make.py.last file
  #   print resulting new command
  # else just use command-line switches

  if redoflag:
    if redoindex == len(redolist): break
    args = redolist[redoindex].split()
    switches,switch_order = parse_args(args)
    redoindex += 1
    
    for switch in cmd_switches:
      if switch == 'r': continue
      if switch == 'a':
        if switch in switches:
          switches[switch] = switches[switch] + cmd_switches[switch]
        else:
          switches[switch] = cmd_switches[switch]
          switch_order.append('a')
      else:
        if switch not in switches:
          if 'a' in switches: switch_order.insert(-1,switch)
          else: switch_order.append(switch)
        switches[switch] = cmd_switches[switch]

    argstr = switch2str(switches,switch_order)
    print "Redo command: Make.py",argstr
  else:
    switches = cmd_switches
    switch_order = cmd_switch_order

  # initialize all class variables to None

  for one in switchclasses: exec("%s = None" % one)
  for one in libclasses: exec("%s = None" % one)
  for one in extraclasses: exec("%s = None" % one)

  # classes = dictionary of created classes
  # key = switch, value = class instance

  classes = {}
  for switch in switches:
    if len(switch) == 1 and switch in abbrevs:
      i = abbrevs.index(switch)
      capitalized = switchclasses[i][0].upper() + switchclasses[i][1:]
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (switchclasses[i],switch,capitalized,switch)
      exec(txt)
    elif len(switch) > 1 and switch in libclasses:
      i = libclasses.index(switch)
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (libclasses[i],switch,libclasses[i].upper(),switch)
      exec(txt)
    elif len(switch) > 1 and switch in extraclasses:
      i = extraclasses.index(switch)
      capitalized = extraclasses[i][0].upper() + extraclasses[i][1:]
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (extraclasses[i],switch,capitalized,switch)
      exec(txt)
    else: error("Unknown command-line switch -%s" % switch)

  # print help messages and exit

  if help or (actions and "-h" in actions.inlist) or not switches:
    if not help: help = Help(None)
    print help.help()
    for switch in switch_order:
      if switch == "h": continue
      print classes[switch].help()[1:]
    sys.exit()

  # create needed default classes if not specified with switch
  # dir and packages, all lib classes so defaults are set

  if not dir: dir = Dir(None)
  if not packages: packages = Packages(None)

  for one in libclasses:
    txt = "if not %s: %s = %s(None)" % (one,one,one.upper())
    exec(txt)

  for one in extraclasses:
    capitalized = one[0].upper() + one[1:]
    txt = "if not %s: %s = %s(None)" % (one,one,capitalized)
    exec(txt)

  # error check on args for all classes

  for switch in classes: classes[switch].check()

  # prep for action
  # actions.setup() detects if last action = machine
  # if yes, induces addition of "-m" and "-o" switches
  
  dir.setup()
  packages.setup()

  if actions:
    machine = actions.setup()
    if machine:
      switches['a'][-1] = "exe"
      if 'm' not in switches:
        switches['m'] = [machine]
        switch_order.insert(-1,'m')
        classes['m'] = Makefile(switches['m'])
        classes['m'].check()
      if 'o' not in switches:
        switches['o'] = [machine]
        switch_order.insert(-1,'o')
        classes['o'] = Makefile(switches['o'])
        classes['o'].check()

  # perform actions

  packages.install()

  if actions:
    for action in actions.alist:
      print "Action %s ..." % action
      if action.startswith("lib-"): actions.lib(action[4:])
      elif action == "file": actions.file(1)
      elif action == "clean": actions.clean()
      elif action == "exe": actions.exe()

  packages.uninstall()
  
  # create output file if requested and exe action performed

  if output and actions and "exe" in actions.alist:
    txt = "cp %s/lmp_auto %s/lmp_%s" % (dir.src,dir.cwd,output.machine)
    commands.getoutput(txt)
      
  # write current Make.py command to src/Make.py.last

  argstr = switch2str(switches,switch_order)
  open("%s/Make.py.last" % dir.src,'w').write(argstr + '\n')

  # if not redoflag, done

  if not redoflag: break
