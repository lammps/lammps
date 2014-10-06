#!/usr/bin/python

# Make.py tool for managing packages and their auxiliary libs,
#   auto-editing machine Makefiles, and building LAMMPS
# Sytax: Make.py -h (for help)

import sys,os,commands,re,copy

# switch abbrevs
# switch classes = created class for each switch
# lib classes = auxiliary package libs
# build classes = build options with defaults
# make classes = makefile options with no defaults
# setargs = makefile settings
# actionargs = allowed actions (also lib-dir and machine)

abbrevs = "adhjmoprsv"

switchclasses = ("actions","dir","help","jmake","makefile",
                 "output","packages","redo","settings","verbose")
libclasses = ("atc","awpmd","colvars","cuda","gpu",
              "meam","poems","qmmm","reax")
buildclasses = ("intel","kokkos")
makeclasses = ("cc","mpi","fft","jpg","png")

setargs = ("gzip","#gzip","ffmpeg","#ffmpeg","smallbig","bigbig","smallsmall")
actionargs = ("lib-all","file","clean","exe")

# ----------------------------------------------------------------
# functions
# ----------------------------------------------------------------

# if flag = 1, print str and exit
# if flag = 0, print str as warning and do not exit

def error(str,flag=1):
  if flag:
    print "ERROR:",str
    sys.exit()
  else:
    print "WARNING:",str

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
    last = first
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

# check if compiler works with ccflags on dummy one-line tmpauto.cpp file
# return 1 if successful, else 0
# warn = 1 = print warning if not successful, warn = 0 = no warning

def compile_check(compiler,ccflags,warn):
  open("tmpauto.cpp",'w').write("int main(int, char **) {}")
  str = "%s %s -c tmpauto.cpp" % (compiler,ccflags)
  txt = commands.getoutput(str)
  if txt or not os.path.isfile("tmpauto.o"):
    flag = 0
    if warn:
      print str
      if txt: print txt
      else: print "compile produced no output"
  else: flag = 1
  os.remove("tmpauto.cpp")
  if os.path.isfile("tmpauto.o"): os.remove("tmpauto.o")
  return flag

# check if linker works with linkflags on tmpauto.o file
# return 1 if successful, else 0
# warn = 1 = print warning if not successful, warn = 0 = no warning

def link_check(linker,linkflags,warn):
  open("tmpauto.cpp",'w').write("int main(int, char **) {}")
  str = "%s %s -o tmpauto tmpauto.cpp" % (linker,linkflags)
  txt = commands.getoutput(str)
  if txt or not os.path.isfile("tmpauto"):
    flag = 0
    if warn:
      print str
      if txt: print txt
      else: print "link produced no output"
  else: flag = 1
  os.remove("tmpauto.cpp")
  if os.path.isfile("tmpauto"): os.remove("tmpauto")
  return flag

# ----------------------------------------------------------------
# switch classes, one per switch
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
    if switches used and machine is only action, must prefix with "-a" switch
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

  # build one or more auxiliary package libraries
  
  def lib(self,suffix):
    if suffix != "all":
      print "building library",suffix,"..."
      str = "%s.build()" % suffix
      exec(str)
    else:
      final = packages.final
      for one in packages.lib:
        if final[one]:
          if "user" in one: pkg = one[5:]
          else: pkg = one
          print "building library",pkg,"..."
          str = "%s.build()" % pkg
          exec(str)

  # read Makefile.machine
  # if caller = "file", edit via switches
  # if caller = "exe", just read
  # write out new Makefile.auto
  
  def file(self,caller):

    # if caller = "file", create from mpi or read from makefile.machine or auto
    # if caller = "exe" and "file" action already invoked, read from auto
    # if caller = "exe" and no "file" action, read from makefile.machine or auto
    
    if caller == "file":
      if makefile and makefile.machine == "none":
        if cc and mpi: machine = "mpi"
        else: error("Cannot create makefile unless -cc and -mpi are used")
      elif makefile: machine = makefile.machine
      else: machine = "auto"
    elif caller == "exe" and "file" in self.alist:
      machine = "auto"
    elif caller == "exe" and "file" not in self.alist:
      if makefile and makefile.machine == "none":
        error("Cannot build with makefile = none")
      elif makefile: machine = makefile.machine
      else: machine = "auto"

    make = MakeReader(machine,1)

    # change makefile settings to user specifications
        
    precompiler = ""
    if caller == "file":
            
      # add compiler/linker and default CCFLAGS,LINKFLAGS
      # if cc.wrap, add wrapper setting for nvcc or mpi = ompi/mpich
      # if cc.wwrap, add 2nd wrapper setting for mpi = ompi/mpich
      # precompiler = env variable setting for OpenMPI wrapper compiler
      
      if cc:
        make.setvar("CC",cc.compiler)
        make.setvar("LINK",cc.compiler)
        if cc.wrap:
          abbrev = cc.abbrev
          if abbrev == "nvcc":
            make.addvar("CC","-ccbin=%s" % cc.wrap)
            make.addvar("LINK","-ccbin=%s" % cc.wrap)
          elif abbrev == "mpi":
            txt = commands.getoutput("mpicxx -show")
            if "-lmpich" in txt:
              make.addvar("CC","-cxx=%s" % cc.wrap)
              make.addvar("LINK","-cxx=%s" % cc.wrap)
            elif "-lmpi" in txt:
              make.addvar("OMPI_CXX",cc.wrap,"cc")
              precompiler = "env OMPI_CXX=%s " % cc.wrap
            else: error("Could not add MPI wrapper compiler, " +
                        "did not recognize OpenMPI or MPICH")
        if cc.wwrap:
          txt = commands.getoutput("mpicxx -show")
          if "-lmpich" in txt:
            make.addvar("CC","-Xcompiler -cxx=%s" % cc.wwrap)
            make.addvar("LINK","-Xcompiler -cxx=%s" % cc.wwrap)
          elif "-lmpi" in txt:
            make.addvar("OMPI_CXX",cc.wwrap,"cc")
            precompiler = "env OMPI_CXX=%s " % cc.wwrap
          else: error("Could not add MPI wrapper compiler, " +
                      "did not recognize OpenMPI or MPICH")
        make.setvar("CCFLAGS","-g")
        make.addvar("CCFLAGS","-O3")
        make.setvar("LINKFLAGS","-g")
        make.addvar("LINKFLAGS","-O")

# add MPI settings

      if mpi:
        make.delvar("MPI_INC","*")
        make.delvar("MPI_PATH","*")
        make.delvar("MPI_LIB","*")
        if mpi.style == "mpi":
          make.addvar("MPI_INC","-DMPICH_SKIP_MPICXX")
          make.addvar("MPI_INC","-DOMPI_SKIP_MPICXX=1")
        elif mpi.style == "mpich":
          make.addvar("MPI_INC","-DMPICH_SKIP_MPICXX")
          make.addvar("MPI_INC","-DOMPI_SKIP_MPICXX=1")
          if mpi.dir: make.addvar("MPI_INC","-I%s/include" % mpi.dir)
          if mpi.dir: make.addvar("MPI_PATH","-L%s/lib" % mpi.dir)
          make.addvar("MPI_LIB","-lmpich")
          make.addvar("MPI_LIB","-lmpl")
          make.addvar("MPI_LIB","-lpthread")
        elif mpi.style == "ompi":
          make.addvar("MPI_INC","-DMPICH_SKIP_MPICXX")
          make.addvar("MPI_INC","-DOMPI_SKIP_MPICXX=1")
          if mpi.dir: make.addvar("MPI_INC","-I%s/include" % mpi.dir)
          if mpi.dir: make.addvar("MPI_PATH","-L%s/lib" % mpi.dir)
          make.addvar("MPI_LIB","-lmpi")
          make.addvar("MPI_LIB","-lmpi_cxx")
        elif mpi.style == "serial":
          make.addvar("MPI_INC","-I../STUBS")
          make.addvar("MPI_PATH","-L../STUBS")
          make.addvar("MPI_LIB","-lmpi_stubs")

      # add accelerator package CCFLAGS and LINKFLAGS and variables
      # pre = "" if compiler not nvcc,
      #   else "-Xcompiler " to pass flag thru to wrapper compiler
          
      compiler = precompiler + ' '.join(make.getvar("CC"))
      linker = precompiler + ' '.join(make.getvar("LINK"))
      if "nvcc" in compiler: pre = "-Xcompiler "
      else: pre = ""
            
      final = packages.final
      if final["opt"]:
        if compile_check(compiler,pre + "-restrict",0):
          make.addvar("CCFLAGS",pre + "-restrict")
          
      if final["user-omp"]:
        if compile_check(compiler,pre + "-restrict",0):
          make.addvar("CCFLAGS",pre + "-restrict")
        #if "nvcc" not in compiler:
          if compile_check(compiler,pre + "-fopenmp",1):
            make.addvar("CCFLAGS",pre + "-fopenmp")
            make.addvar("LINKFLAGS",pre + "-fopenmp")

      if final["user-intel"]:
        if intel.mode == "cpu":
          if compile_check(compiler,pre + "-fopenmp",1):
            make.addvar("CCFLAGS",pre + "-fopenmp")
            make.addvar("LINKFLAGS",pre + "-fopenmp")
          make.addvar("CCFLAGS",pre + "-DLAMMPS_MEMALIGN=64")
          if compile_check(compiler,pre + "-restrict",1):
            make.addvar("CCFLAGS",pre + "-restrict")
          if compile_check(compiler,pre + "-xHost",1):
            make.addvar("CCFLAGS",pre + "-xHost")
            make.addvar("LINKFLAGS",pre + "-xHost")
          if compile_check(compiler,pre + "-fno-alias",1):
            make.addvar("CCFLAGS",pre + "-fno-alias")
          if compile_check(compiler,pre + "-ansi-alias",1):
            make.addvar("CCFLAGS",pre + "-ansi-alias")
          if compile_check(compiler,pre + "-override-limits",1):
            make.addvar("CCFLAGS",pre + "-override-limits")
          make.delvar("CCFLAGS",pre + "-DLMP_INTEL_OFFLOAD")
          make.delvar("LINKFLAGS",pre + "-offload")
        elif intel.mode == "phi":
          if compile_check(compiler,pre + "-fopenmp",1):
            make.addvar("CCFLAGS",pre + "-fopenmp")
            make.addvar("LINKFLAGS",pre + "-fopenmp")
          make.addvar("CCFLAGS",pre + "-DLAMMPS_MEMALIGN=64")
          if compile_check(compiler,pre + "-restrict",1):
            make.addvar("CCFLAGS",pre + "-restrict")
          if compile_check(compiler,pre + "-xHost",1):
            make.addvar("CCFLAGS",pre + "-xHost")
          make.addvar("CCFLAGS",pre + "-DLMP_INTEL_OFFLOAD")
          if compile_check(compiler,pre + "-fno-alias",1):
            make.addvar("CCFLAGS",pre + "-fno-alias")
          if compile_check(compiler,pre + "-ansi-alias",1):
            make.addvar("CCFLAGS",pre + "-ansi-alias")
          if compile_check(compiler,pre + "-override-limits",1):
            make.addvar("CCFLAGS",pre + "-override-limits")
          if compile_check(compiler,pre + '-offload-option,mic,compiler,' +
                          '"-fp-model fast=2 -mGLOB_default_function_attrs=' +
                          '\\"gather_scatter_loop_unroll=4\\""',1):
            make.addvar("CCFLAGS",pre + '-offload-option,mic,compiler,' +
                        '"-fp-model fast=2 -mGLOB_default_function_attrs=' +
                        '\\"gather_scatter_loop_unroll=4\\""')
          if link_check(linker,"-offload",1):
            make.addvar("LINKFLAGS",pre + "-offload")

      if final["kokkos"]:
        if kokkos.mode == "omp":
          make.addvar("OMP","yes","lmp")
          make.delvar("CUDA")
          make.delvar("MIC")
        elif kokkos.mode == "cuda":
          if "nvcc" not in compiler:
            error("Kokkos/cuda build appears to not be " +
                  "using NVIDIA nvcc compiler",0)
          make.addvar("OMP","yes","lmp")
          make.addvar("CUDA","yes","lmp")
          make.delvar("MIC")
          if kokkos.archflag:
            make.delvar("CCFLAGS","-arch=sm_*")
            make.addvar("CCFLAGS","-arch=sm_%s" % kokkos.arch)
        elif kokkos.mode == "phi":
          make.addvar("OMP","yes","lmp")
          make.addvar("MIC","yes","lmp")
          make.delvar("CUDA")

      # add LMP settings
      
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
          
      # add FFT, JPG, PNG settings

      if fft:
        make.delvar("FFT_INC","*")
        make.delvar("FFT_PATH","*")
        make.delvar("FFT_LIB","*")
        if fft.mode == "none": make.addvar("FFT_INC","-DFFT_NONE")
        else:
          make.addvar("FFT_INC","-DFFT_%s" % fft.mode.upper())
          make.addvar("FFT_LIB",fft.lib)
          if fft.dir:
            make.addvar("FFT_INC","-I%s/include" % fft.dir)
            make.addvar("FFT_PATH","-L%s/lib" % fft.dir)
          else:
            if fft.incdir: make.addvar("FFT_INC","-I%s" % fft.incdir)
            if fft.libdir: make.addvar("FFT_PATH","-L%s" % fft.libdir)

      if jpg:
        if jpg.on == 0:
          make.delvar("LMP_INC","-DLAMMPS_JPEG")
          make.delvar("JPG_LIB","-ljpeg")
        else:
          make.addvar("LMP_INC","-DLAMMPS_JPEG")
          make.addvar("JPG_LIB","-ljpeg")
          if jpg.dir:
            make.addvar("JPG_INC","-I%s/include" % jpg.dir)
            make.addvar("JPG_PATH","-L%s/lib" % jpg.dir)
          else:
            if jpg.incdir: make.addvar("JPG_INC","-I%s" % jpg.incdir)
            if jpg.libdir: make.addvar("JPG_PATH","-L%s" % jpg.libdir)

      if png:
        if png.on == 0:
          make.delvar("LMP_INC","-DLAMMPS_PNG")
          make.delvar("JPG_LIB","-lpng")
        else:
          make.addvar("LMP_INC","-DLAMMPS_PNG")
          make.addvar("JPG_LIB","-lpng")
          if png.dir:
            make.addvar("JPG_INC","-I%s/include" % png.dir)
            make.addvar("JPG_PATH","-L%s/lib" % png.dir)
          else:
            if png.incdir: make.addvar("JPG_INC","-I%s" % png.incdir)
            if png.libdir: make.addvar("JPG_PATH","-L%s" % png.libdir)

    # write out Makefile.auto
    # unless caller = "exe" and "file" action already invoked

    if caller == "file" or "file" not in self.alist:
      make.write("%s/MAKE/MINE/Makefile.auto" % dir.src,1)
      print "Created src/MAKE/MINE/Makefile.auto"

    # test full compile and link
    # unless caller = "file" and "exe" action will be invoked later

    if caller == "file" and "exe" in self.alist: return
    compiler = precompiler + ' '.join(make.getvar("CC"))
    ccflags = ' '.join(make.getvar("CCFLAGS"))
    linker = precompiler + ' '.join(make.getvar("LINK"))
    linkflags = ' '.join(make.getvar("LINKFLAGS"))
    if not compile_check(compiler,ccflags,1):
      error("Test of compilation failed")
    if not link_check(linker,linkflags,1): error("Test of link failed")

  # invoke "make clean-auto" to force clean before build
    
  def clean(self):
    str = "cd %s; make clean-auto" % dir.src
    commands.getoutput(str)
    if verbose: print "Performed make clean-auto"

  # build LAMMPS using Makefile.auto and -j setting
  # invoke self.file() first, to test makefile compile/link
  # delete existing lmp_auto, so can detect if build fails

  def exe(self):
    self.file("exe")
    commands.getoutput("cd %s; rm -f lmp_auto" % dir.src)
    if jmake: str = "cd %s; make -j %d auto" % (dir.src,jmake.n)
    else: str = "cd %s; make auto" % dir.src
    txt = commands.getoutput(str)
    if verbose: print txt
    if not os.path.isfile("%s/lmp_auto" % dir.src):
      if not verbose: print txt
      error('Unsuccessful "make auto"')
    elif not output: print "Created src/lmp_auto"
    
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
  switches for build and makefile options:
    -intel, -kokkos, -cc, -mpi, -fft, -jpg, -png

add -h switch to command line to print this message
  and help on other specified switches or actions
  add -a switch if not seeing action help
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
  use Makefile.machine under src/MAKE as starting point to create Makefile.auto
  if machine = "none", file action will create Makefile.auto from scratch
    must use -cc and -mpi switches to specify compiler and MPI
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
    all = all standard and user packages (also none = no-all)
    std (or standard) = all standard packages
    user = all user packages
    lib = all standard and user packages with auxiliary libs
  can abbreviate package names and yes/no
    omp = user-omp = yes-user-omp
    ^omp = ^user-omp = no-user-omp
    user = yes-user, ^user = no-user
    all = yes-all, ^all = none = no-all
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
        elif one == '^' and one[1:] in std:
          plist.append("no-%s" % one[1:])
        elif one[0] == '^' and one[1:] in user:
          plist.append("no-%s" % one[1:])
        elif one[0] == '^' and "user-"+one[1:] in user:
          plist.append("no-user-%s" % one[1:])
        elif one == "std" or one == "standard" or one == "user" or \
              one == "lib" or one == "all": 
          plist.append("yes-%s" % one)
        elif one == "^std" or one == "^standard" or one == "^user" or \
              one == "^lib" or one == "^all": 
          plist.append("no-%s" % one[1:])
        elif one == "none": plist.append("no-all")
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
    file = self.file
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
# lib classes, one per LAMMPS auxiliary lib
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
# build classes for intel, kokkos build options
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
# makefile classes for CC, MPI, JPG, PNG, FFT settings
# ----------------------------------------------------------------

# Cc class

class Cc:
  def __init__(self,list):
    self.inlist = list[:]
    self.compiler = self.abbrev = ""
    self.wrap = self.wabbrev = ""
    self.wwrap = ""

  def help(self):
    return """
-cc compiler wrap=wcompiler wwrap=wwcompiler
  change CC setting in makefile
  compiler is required, all other args are optional
  compiler = any string with g++ or icc or icpc or nvcc
             or mpi (or mpicxx, mpiCC, mpiicpc, etc)
    can be compiler name or full path to compiler
    mpi by itself is changed to mpicxx
  wcompiler = g++ or icc or icpc or mpi (or mpicxx, mpiCC, mpiicpc, etc)
    can only be used when compiler is a wrapper (mpi or nvcc)
    mpi and variants can only be used with compiler = nvcc
    mpi by itself is changed to mpicxx
    specify compiler for wrapper to use in CC setting
  wwcompiler = g++ or icc or icpc
    can only be used when wcompiler is itself a wrapper (mpi)
    specify compiler for wrapper of wrapper to use in CC setting
"""

  def check(self):
    if len(self.inlist) < 1: error("-cc args are invalid")
    self.compiler = self.inlist[0]
    if self.compiler == "mpi":
      self.compiler = "mpicxx"
      self.abbrev = "mpi"
    elif self.compiler.startswith("mpi"):
      self.abbrev = "mpi"
    elif self.compiler == "g++" or self.compiler == "icc" or \
          self.compiler == "icpc" or self.compiler == "nvcc":
      self.abbrev = self.compiler
    elif "mpi" in self.compiler: self.abbrev = "mpi"
    elif "g++" in self.compiler: self.abbrev = "g++"
    elif "icc" in self.compiler: self.abbrev = "icc"
    elif "icpc" in self.compiler: self.abbrev = "icpc"
    elif "nvcc" in self.compiler: self.abbrev = "nvcc"
    else: error("-cc args are invalid")
    for one in self.inlist[1:]:
      words = one.split('=')
      if len(words) != 2: error("-cc args are invalid")
      if words[0] == "wrap":
        if self.abbrev != "mpi" and self.abbrev != "nvcc":
          error("-cc compiler is not a wrapper")
        self.wrap = words[1]
        if self.wrap == "mpi":
          self.wrap = "mpicxx"
          self.wabbrev = "mpi"
        elif self.wrap.startswith("mpi"):
          self.wabbrev = "mpi"
        elif self.compiler == "g++" or self.compiler == "icc" or \
              self.compiler == "icpc":
          self.wabbrev = self.wrap
      elif words[0] == "wwrap":
        self.wwrap = words[1]
        if self.wabbrev != "mpi": error("-cc wrap is not a wrapper")
        if self.wwrap != "g++" and self.wwrap != "icc" and self.wwrap != "icpc":
          error("-cc args are invalid")
      else: error("-cc args are invalid")

# Mpi class

class Mpi:
  def __init__(self,list):
    self.inlist = list[:]
    self.style = self.dir = ""
                
  def help(self):
    return """
-mpi style dir=path
  change MPI settings in makefile
  style is required, all other args are optional
  style = mpi or mpich or ompi or serial
    mpi = no MPI settings (assume compiler is MPI wrapper)
    mpich = use explicit settings for MPICH
    ompi = use explicit settings for OpenMPI
    serial = use settings for src/STUBS library
  dir = path for MPICH or OpenMPI directory
    add -I and -L settings for include and lib sub-dirs
"""

  def check(self):
    if len(self.inlist) < 1: error("-mpi args are invalid")
    self.style = self.inlist[0]
    if self.style != "mpi" and self.style != "mpich" and \
      self.style != "ompi" and self.style != "serial":
      error("-mpi args are invalid")
    for one in self.inlist[1:]:
      words = one.split('=')
      if len(words) != 2: error("-mpi args are invalid")
      if words[0] == "dir": self.dir = words[1]
      else: error("-mpi args are invalid")

# Fft class

class Fft:
  def __init__(self,list):
    self.inlist = list[:]
    self.dir = self.incdir = self.libdir = ""

  def help(self):
    return """
-fft mode lib=libname dir=homedir idir=incdir ldir=libdir
  change FFT settings in makefile
  mode is required, all other args are optional
  removes all current FFT variable settings
  mode = none or fftw or fftw3 of ...
    adds -DFFT_MODE setting
  lib = name of FFT library to link with (def is libname = mode)
    adds -lliblib setting, e.g. -llibfftw3
  dir = home dir for include and library files (def = none)
    adds -Idir/include and -Ldir/lib settings
    if set, overrides idir and ldir args
  idir = dir for include file (def = none)
    adds -Iidir setting
  ldir = dir for library file (def = none)
    adds -Lldir setting
"""

  def check(self):
    if not len(self.inlist): error("-fft args are invalid")
    self.mode = self.inlist[0]
    self.lib = "-l%s" % self.mode
    for one in self.inlist[1:]:
      words = one.split('=')
      if len(words) != 2: error("-fft args are invalid")
      if words[0] == "lib": self.lib = "-l%s" % words[1]
      elif words[0] == "dir": self.dir = words[1]
      elif words[0] == "idir": self.incdir = words[1]
      elif words[0] == "ldir": self.libdir = words[1]
      else: error("-fft args are invalid")

# Jpg class

class Jpg:
  def __init__(self,list):
    self.inlist = list[:]
    self.on = 1
    self.dir = self.incdir = self.libdir = ""

  def help(self):
    return """
-jpg flag dir=homedir idir=incdir ldir=libdir
  change JPG settings in makefile
  all args are optional, flag must come first if specified
  flag = yes or no (def = yes)
    include or exclude JPEG support
    adds/removes -DLAMMPS_JPEG and -ljpeg settings
  dir = home dir for include and library files (def = none)
    adds -Idir/include and -Ldir/lib settings
    if set, overrides idir and ldir args
  idir = dir for include file (def = none)
    adds -Iidir setting
  ldir = dir for library file (def = none)
    adds -Lldir setting
"""

  def check(self):
    for i,one in enumerate(self.inlist):
      if one == "no" and i == 0: self.on = 0
      elif one == "yes" and i == 0: self.on = 1
      else:
        words = one.split('=')
        if len(words) != 2: error("-jpeg args are invalid")
        if words[0] == "dir": self.dir = words[1]
        elif words[0] == "idir": self.incdir = words[1]
        elif words[0] == "ldir": self.libdir = words[1]
        else: error("-jpeg args are invalid")

# Png class

class Png:
  def __init__(self,list):
    self.inlist = list[:]
    self.on = 1
    self.dir = self.incdir = self.libdir = ""

  def help(self):
    return """
-png flag dir=homedir idir=incdir ldir=libdir
  change PNG settings in makefile
  all args are optional, flag must come first if specified
  flag = yes or no (def = yes)
    include or exclude PNG support
    adds/removes -DLAMMPS_PNG and -lpng settings
  dir = home dir for include and library files (def = none)
    adds -Idir/include and -Ldir/lib settings
    if set, overrides idir and ldir args
  idir = dir for include file (def = none)
    adds -Iidir setting
  ldir = dir for library file (def = none)
    adds -Lldir setting
"""

  def check(self):
    for i,one in enumerate(self.inlist):
      if one == "no" and i == 0: self.on = 0
      elif one == "yes" and i == 0: self.on = 1
      else:
        words = one.split('=')
        if len(words) != 2: error("-png args are invalid")
        if words[0] == "dir": self.dir = words[1]
        elif words[0] == "idir": self.incdir = words[1]
        elif words[0] == "ldir": self.libdir = words[1]
        else: error("-png args are invalid")

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
    # ccindex = index of "CC =" line, to add OMPI var before it
    # lmpindex = index of "LAMMPS-specific settings" line to add KOKKOS vars before it
    
    var = {}
    varinfo = []
    newlines = []
    pattern = "(\S+\s+=\s+)(.*)"
    multiline = 0
    self.ccindex = self.lmpindex = 0
    
    for line in lines:
      line = line[:-1]
      if "CC =" in line: self.ccindex = len(newlines)
      if "LAMMPS-specific settings" in line: self.lmpindex = len(newlines)
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
  #   create new variable using where
  #   where="cc", line before "CC =" line, use ":="
  #   where="lmp", 2 lines before "LAMMPS-specific settings" line, use "="
  
  def addvar(self,var,value,where=""):
    if var in self.var:
      if value not in self.var[var]: self.var[var].append(value)
    else:
      if not where:
        error("Variable %s with value %s is not in makefile" % (var,value))
      if where == "cc":
        if not self.ccindex: error("No 'CC =' line in makefile to add variable")
        index = self.ccindex
        varwhite = "%s :=\t\t" % var
      elif where == "lmp":
        if not self.lmpindex: error("No 'LAMMPS-specific settings line' " +
                                    "in makefile to add variable")
        index = self.lmpindex - 2
        varwhite = "%s =\t\t" % var
      self.var[var] = [value]
      varwhite = "%s =\t\t" % var
      self.lines.insert(index,str(len(self.varinfo)))
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

if "v" in cmd_switches:
  print "Command-line parsing:"
  for switch in cmd_switch_order:
    print "  %s: %s" % (switch,' '.join(cmd_switches[switch]))

# check for redo switch, process redo file
# redolist = list of commands to execute

redoflag = 0
redolist = []

if 'r' in cmd_switches and 'h' not in cmd_switches:
  redoflag = 1
  redo = Redo(cmd_switches['r'])
  redo.check()
  redo.setup()
  redolist = redo.commands
  redoindex = 0
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
  for one in buildclasses: exec("%s = None" % one)
  for one in makeclasses: exec("%s = None" % one)
  
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
    elif switch in libclasses:
      i = libclasses.index(switch)
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (libclasses[i],switch,libclasses[i].upper(),switch)
      exec(txt)
    elif switch in buildclasses:
      i = buildclasses.index(switch)
      capitalized = buildclasses[i][0].upper() + buildclasses[i][1:]
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (buildclasses[i],switch,capitalized,switch)
      exec(txt)
    elif switch in makeclasses:
      i = makeclasses.index(switch)
      capitalized = makeclasses[i][0].upper() + makeclasses[i][1:]
      txt = '%s = classes["%s"] = %s(switches["%s"])' % \
          (makeclasses[i],switch,capitalized,switch)
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
  # dir and packages plus lib and build classes so defaults are set

  if not dir: dir = Dir(None)
  if not packages: packages = Packages(None)

  for one in libclasses:
    txt = "if not %s: %s = %s(None)" % (one,one,one.upper())
    exec(txt)

  for one in buildclasses:
    capitalized = one[0].upper() + one[1:]
    txt = "if not %s: %s = %s(None)" % (one,one,capitalized)
    exec(txt)

  # error check on args for all classes

  for switch in classes: classes[switch].check()

  # prep for action
  # actions.setup() detects if last action = machine
  # if yes, induce addition of "-m" and "-o" switches
  
  dir.setup()
  packages.setup()

  if actions:
    machine = actions.setup()
    if machine:
      switches['a'][-1] = "exe"
      if 'm' not in switches:
        switches['m'] = [machine]
        switch_order.insert(-1,'m')
        makefile = classes['m'] = Makefile(switches['m'])
        makefile.check()
      if 'o' not in switches:
        switches['o'] = [machine]
        switch_order.insert(-1,'o')
        output = classes['o'] = Makefile(switches['o'])
        output.check()

  # perform actions

  packages.install()

  if actions:
    for action in actions.alist:
      print "Action %s ..." % action
      if action.startswith("lib-"): actions.lib(action[4:])
      elif action == "file": actions.file("file")
      elif action == "clean": actions.clean()
      elif action == "exe": actions.exe()

  packages.uninstall()
  
  # create output file if requested and exe action performed

  if output and actions and "exe" in actions.alist:
    txt = "cp %s/lmp_auto %s/lmp_%s" % (dir.src,dir.cwd,output.machine)
    commands.getoutput(txt)
    print "Created lmp_%s in %s" % (output.machine,dir.cwd)

  # write current Make.py command to src/Make.py.last

  argstr = switch2str(switches,switch_order)
  open("%s/Make.py.last" % dir.src,'w').write(argstr + '\n')

  # if not redoflag, done

  if not redoflag: break
