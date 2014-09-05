#!/usr/local/bin/python

# Syntax: build.py target1 target2 ...
#         targets:
#         cpu, opt, omp,
#         gpu/double, gpu/mixed, gpu/single,
#         cuda/double, cuda/mixed, cuda/single,
#         intel/cpu, intel/phi,
#         kokkos/omp, kokkos/phi, kokkos/cuda
#         gpu = gpu/double + gpu/mixed + gpu/single
#         cuda = cuda/double + cuda/mixed + cuda/single
#         intel = intel/cpu + intel/phi
#         kokkos = kokkos/omp + kokkos/phi + kokkos/cuda
#         all = cpu + opt + omp + gpu + cuda + intel + kokkos

# create exectuables for different packages
# MUST set lmpdir to path of LAMMPS home directory

import sys,commands,os

lmpdir = "~/lammps"

# build LAMMPS
# copy makefile into src/MAKE as Makefile.foo, then remove it

def build_lammps(makefile,pkg):
  print "Building LAMMPS with %s and %s packages ..." % (makefile,pkg)
  commands.getoutput("cp %s %s/src/MAKE/Makefile.foo" % (makefile,lmpdir))
  cwd = os.getcwd()
  os.chdir(os.path.expanduser(lmpdir + "/src"))
  str = "make clean-foo"
  txt = commands.getoutput(str)
  str = "make no-all"
  txt = commands.getoutput(str)
  for package in pkg:
    str = "make yes-%s" % package
    txt = commands.getoutput(str)
    print txt
  str = "make -j 16 foo"
  txt = commands.getoutput(str)
  os.remove("MAKE/Makefile.foo")
  os.chdir(cwd)

# build GPU library in LAMMPS
# copy makefile into lib/gpu as Makefile.foo, then remove it
  
def build_gpu(makefile):
  print "Building GPU lib with %s ..." % makefile
  commands.getoutput("cp %s %s/lib/gpu/Makefile.foo" % (makefile,lmpdir))
  cwd = os.getcwd()
  os.chdir(os.path.expanduser(lmpdir + "/lib/gpu"))
  str = "make -f Makefile.foo clean"
  txt = commands.getoutput(str)
  str = "make -j 16 -f Makefile.foo"
  txt = commands.getoutput(str)
  os.remove("Makefile.foo")
  os.chdir(cwd)

# build CUDA library in LAMMPS
# set precision and arch explicitly as options to make in lib/cuda
  
def build_cuda(precision,arch):
  print "Building USER-CUDA lib with %s and arch sm_%d ..." % (precision,arch)
  cwd = os.getcwd()
  os.chdir(os.path.expanduser(lmpdir + "/lib/cuda"))
  str = "make clean"
  txt = commands.getoutput(str)
  if precision == "double": pflag = 2
  elif precision == "mixed": pflag = 4
  elif precision == "single": pflag = 1
  str = "make -j 16 precision=%d arch=%s" % (pflag,arch)
  txt = commands.getoutput(str)

  os.chdir(cwd)

# main program
# convert target keywords into target flags

cpu = opt = omp = 0
gpu_double = gpu_mixed = gpu_single = 0
cuda_double = cuda_mixed = cuda_single = 0
intel_cpu = intel_phi = 0
kokkos_omp = kokkos_phi = kokkos_cuda = 0
  
targets = sys.argv[1:]
for target in targets:
  if target == "cpu": cpu = 1
  elif target == "opt": opt = 1
  elif target == "omp": omp = 1
  elif target == "gpu/double": gpu_double = 1
  elif target == "gpu/mixed": gpu_mixed = 1
  elif target == "gpu/single": gpu_single = 1
  elif target == "gpu": gpu_double = gpu_mixed = gpu_single = 1
  elif target == "cuda/double": cuda_double = 1
  elif target == "cuda/mixed": cuda_mixed = 1
  elif target == "cuda/single": cuda_single = 1
  elif target == "cuda": cuda_double = cuda_mixed = cuda_single = 1
  elif target == "intel/cpu": intel_cpu = 1
  elif target == "intel/phi": intel_phi = 1
  elif target == "intel": intel_cpu = intel_phi = 1
  elif target == "kokkos/omp": kokkos_omp = 1
  elif target == "kokkos/phi": kokkos_phi = 1
  elif target == "kokkos/cuda": kokkos_cuda = 1
  elif target == "kokkos": kokkos_omp = kokkos_phi = kokkos_cuda = 1
  else: print "Target",target,"is unknown"

# CPU

if cpu:
  build_lammps(makefile = "Makefile.cpu", pkg = [])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_cpu" % lmpdir)

# OPT

if opt:
  build_lammps(makefile = "Makefile.opt", pkg = ["opt"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_opt" % lmpdir)

# OMP

if omp:
  build_lammps(makefile = "Makefile.omp", pkg = ["user-omp"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_omp" % lmpdir)

# GPU, 3 precisions

if gpu_double:
  build_gpu(makefile = "Makefile.gpu.double")
  build_lammps(makefile = "Makefile.gpu", pkg = ["gpu"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_gpu_double" % lmpdir)

if gpu_mixed:
  build_gpu(makefile = "Makefile.gpu.mixed")
  build_lammps(makefile = "Makefile.gpu", pkg = ["gpu"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_gpu_mixed" % lmpdir)

if gpu_single:
  build_gpu(makefile = "Makefile.gpu.single")
  build_lammps(makefile = "Makefile.gpu", pkg = ["gpu"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_gpu_single" % lmpdir)

# CUDA, 3 precisions

if cuda_double:
  build_cuda(precision = "double", arch = 35)
  build_lammps(makefile = "Makefile.cuda", pkg = ["kspace","user-cuda"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_cuda_double" % lmpdir)

if cuda_mixed:
  build_cuda(precision = "mixed", arch = 35)
  build_lammps(makefile = "Makefile.cuda", pkg = ["kspace","user-cuda"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_cuda_mixed" % lmpdir)

if cuda_single:
  build_cuda(precision = "single", arch = 35)
  build_lammps(makefile = "Makefile.cuda", pkg = ["kspace","user-cuda"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_cuda_single" % lmpdir)

# INTEL, CPU and Phi

if intel_cpu:
  build_lammps(makefile = "Makefile.intel.cpu", pkg = ["user-intel"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_intel_cpu" % lmpdir)

if intel_phi:
  build_lammps(makefile = "Makefile.intel.phi", pkg = ["user-intel","user-omp"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_intel_phi" % lmpdir)

# KOKKOS, all variants

if kokkos_omp:
  build_lammps(makefile = "Makefile.kokkos.omp", pkg = ["kokkos"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_kokkos_omp" % lmpdir)

if kokkos_phi:
  build_lammps(makefile = "Makefile.kokkos.phi", pkg = ["kokkos"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_kokkos_phi" % lmpdir)

if kokkos_cuda:
  build_lammps(makefile = "Makefile.kokkos.cuda", pkg = ["kokkos"])
  print commands.getoutput("mv %s/src/lmp_foo ./lmp_kokkos_cuda" % lmpdir)
