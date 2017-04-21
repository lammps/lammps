#!/usr/bin/env python

#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov

#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.

#   See the README file in the top-level LAMMPS directory.


# regression.py tool for numerical comparisions of benchmark log files
# Created by Stan Moore (SNL), stamoor at sandia.gov
# based on benchmark.py created by Reese Jones (SNL)
# requires log.py from Pizza.py toolkit, http://pizza.sandia.gov/

usage = """
  regression.py: numerical comparisions of logs and corresponding benchmarks
  Syntax: regression.py <descriptor> <LAMMPS_args> <test_dirs> <options>
    descriptor = any string without spaces, appended to log files
    LAMMPS_args = string to launch the benchmark calculation
      the path to the executable must be an absolute path
      e.g. ~/lammps/src/lmp_g++ or "mpirun -np 4 ~/lammps/src/lmp_g++ -v x 10"
    test_dirs = list of one or more dirs to recursively search for scripts
      scripts = any in.* file
    options = one or more keyword/value pairs
      wildcards are expanded by Python, not the shell, so it may be necessary to escape * as \*
      -exclude <subdir1 subdir2* ...>
        do not run tests from these sub-dirs or their children
        default = none
      -only <subdir1 subdir2* ...>
        only run tests from these sub-dirs or their children
        default = none
      -customonly <file1 file2* ...>
        only run tests from sub-dirs that contain these files
        default = none
      -custom <file_prefix>
        read options from this file_prefix plus test name in each sub-dir, if it exists
        valid options are: launch, descriptors, tolerance, error_norm, relative_error
        the number of launches and descriptors must match
        lines in file have syntax like:
          descriptors = ["1","4","8"]
          error_norm = "L1"
        default = "options"
      -error_norm <"L1" or "L2" or "max">
        metric for comparing a column of output to gold standard answer
        these are vector norms, treating the column as a vector
        default = "max"
      -relative_error <"True" or "False">
        treat tolerance as a relative error or not
        default = "False"
      -tolerance <float>
        column difference > tolerance = fail, else success
        default = 1.0e-7
      -logread <dir module>
        path for where to find the log-file reading Python module
        default = . log (log.py in this dir or somewhere Python can find it)
"""

import sys
import os
import math
import re
from operator import itemgetter
from glob import glob
import time

import shutil
import platform

#====================================================
### global variables
#====================================================
nrows = 0
not_same_rows = set()
auto_rebless_flag = False
min_same_rows = -1
custom_file = "options"
default_error_norm = "max"
default_relative_error = False
default_tolerance = 1.e-7
logread = [".","log"]

#====================================================
### constants
#====================================================
fail_pattern = re.compile("FAIL");
warn_pattern = re.compile("WARNING");

#====================================================
### date
#====================================================
def date():
 return time.asctime()

#====================================================
### timer
#====================================================
def start():
  global dt
  dt = -(time.time())
def stop():
  global dt
  dt += (time.time())
  return dt

#====================================================
### run a regression test
#====================================================
def run_test(test,lmps,descriptor):
  global not_same_rows
  global nrows
  msg = ""
  input = "in."+test
  log = "log."+descriptor+"."+test
  stdout = "stdout."+descriptor+"."+test
  new_flag = False

  # print test header
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "dir =",os.getcwd()
  print "test =",test
  sys.stdout.flush()
  if (custom_flag):
    print "descriptor =",descriptor
    print "norm =",error_norm
    print "tolerance =",tolerance
    print "relative error =",relative_error
    sys.stdout.flush()

  # check if gold standard exists, if not create it
  system_name = platform.system()
  gold_standard = glob("log.*"+system_name.lower()+"."+descriptor+"."+test)
  if (len(gold_standard) > 0):
    ref = (gold_standard)[0];
    print "gold standard =",ref
    sys.stdout.flush()
  else:
    new_flag = True
    msg += add_test(test,lmps,descriptor)
    gold_standard = glob("log.*"+system_name.lower()+"."+descriptor+"."+test)
    if not (len(gold_standard) > 0):
      raise Exception("No logfile found")
    ref = (gold_standard)[0];

  # compare current run to gold standard
  msg += "==== comparing "+log+" with "+ref+" ====\n"
  # remove old log and stdout files
  if (os.path.isfile(log)): os.remove(log)
  if (os.path.isfile(stdout)): os.remove(stdout)
  # run the test
  os.system(lmps+" -in "+input+" -log "+log+" >& "+stdout);
  # check if a log file was generated
  if (not os.path.isfile(log)) :
    msg += "!!! no "+log+"\n";
    msg += "!!! test "+test+" FAILED\n"
    return msg
  # extract data
  [passing,cdict,cdata,emsg] = extract_data(log,stdout);
  [passing,bdict,bdata,emsg] = extract_data(ref,stdout);
  msg += emsg
  fail = False
  if (not passing) : fail = True
  if (fail) :
    msg += "!!! test "+test+" FAILED\n"
    if (new_flag):
      os.remove(ref)
    return msg
  # compare columns
  not_same_rows.clear()
  cols = range(len(bdict))
  if (len(cdata) != len(bdata)):
    msg += "!!! data size "+str(len(cdata))+" does not match data "+str(len(bdata))+" in "+ref+"\n";
    msg += "!!! test "+test+" FAILED\n"
    return msg 
  i = 0
  for name in bdict:
    [passing,cmsg] = compare(name,cdata[cols[i]],bdata[cols[i]]);
    i += 1
    msg += cmsg
    if (not passing) : fail = True

  # print out results
  nsame_rows = nrows - len(not_same_rows)
  if (not fail):
    if (nrows == nsame_rows):
      msg += "\nAll rows are identical\n"
    else:
      msg += "\nWARNING: Only "+str(nsame_rows)+" out of "+str(nrows)+" rows are identical\n"
      if (auto_rebless_flag):
        dmy = time.strftime("%d%b%y")
        hms = time.strftime("%H:%M:%S")
        shutil.copyfile(ref,"old_"+ref+"_"+dmy+"_"+hms)
        shutil.copyfile(log,ref)
        msg += "WARNING: Gold standard for test "+test+" has been auto-reblessed\n"
  if (fail) :
    msg += "!!! test "+test+" FAILED\n"
  else : 
    msg += "*** test "+test+" passed\n"
  return msg

#====================================================
### add a regression test
#====================================================
def add_test(test,lmps,descriptor):
  input = "in."+test;
  log = "log."+descriptor+"."+test
  stdout = "stdout."+descriptor+"."+test
  msg = "==== generating gold standard for test "+test+" ====\n"
  if (os.path.isfile(log)): os.remove(log)
  if (os.path.isfile(stdout)): os.remove(stdout)
  os.system(lmps+" -in "+input+" -log "+log+" >& "+stdout);
  if (not os.path.isfile(log)) :
    msg += "!!! no "+log+"\n";
    msg += "!!! test "+test+" FAILED\n"
    return msg
  dmy = time.strftime("%d%b%y")
  system_name = platform.system()
  shutil.copyfile(log,"log."+dmy+"."+system_name.lower()+"."+descriptor+"."+test)
  return msg

#====================================================
### extract data from log file
#====================================================
def extract_data(file,stdout):
  msg = ""
  dictionary = [];
  data = []
  read = False

  msg = error_check(file,stdout)
  if (msg != ""):
    return [False,dictionary,data,msg]

  try:
    if logreader.__name__ == "log": lg= logreader(file)
    elif logreader.__name__ == "olog": lg= logreader(file,"Step")
    else: raise Exception("Unknown log reader")
  except:
    msg += "Invalid logfile found\n"
    return [False,dictionary,data,msg]

  if (len(lg.names) <= 0):
    msg += "Invalid logfile found\n"
    return [False,dictionary,data,msg]

  for name in lg.names:
    if (name == "CPU"): continue
    dictionary.append(name)
    data.append(lg.get(name))
  return [True,dictionary,data,msg]

#====================================================
### check log and stdout file for errors
#====================================================
def error_check(file,stdout):
  msg = ""
  text_file = open(file, "r")
  lines = text_file.readlines()
  text_file.close()
  num_lines = int(len(lines))
  # check for errors
  for i in xrange(num_lines):
    if "ERROR" in lines[i] or "exited on signal" in lines[i]:
      msg += lines[i]

  text_file = open(stdout, "r")
  lines = text_file.readlines()
  text_file.close()
  num_lines = int(len(lines))
  # check for errors
  for i in xrange(num_lines):
    if "ERROR" in lines[i] or "exited on signal" in lines[i]:
      msg += lines[i]

  return msg

#====================================================
### compare columns of current run and gold standard
#====================================================
def compare(name,colA,colB):
  msg = ""
  err1 = 0.
  err2 = 0.
  errmax = 0.
  norm1 = 0.
  norm2 = 0.
  normmax = 0.
  nsame_rows = 0
  global nrows
  global not_same_rows
  n = len(colB)
  nrows = n
  if (len(colA) != len(colB)):
    msg = "Cannot compare columns\n"
    return [False,msg];
  for i in range(n):
    vA = float(colA[i])
    vB = float(colB[i])
    norm1 += abs(vB)
    norm2 += vB*vB
    normmax = max(normmax,abs(vB))
    dv = vA-vB
    if (abs(dv) > tolerance):
      not_same_rows.add(i)
    else:
      nsame_rows += 1
    err1 += abs(dv)
    err2 += dv*dv
    errmax = max(errmax,abs(dv))
  norm1 /= n
  norm2 = math.sqrt(norm2/n)
  err1 /= n
  err2 = math.sqrt(err2/n)

  if error_norm == "L1":
    err = err1
    norm = norm1
  elif error_norm == "L2":
    err = err2
    norm = norm2
  elif error_norm == "max":
    err = errmax
    norm = normmax
  else:
    raise Exception("Invalid error norm")

  if (relative_error and norm > tolerance):
    err /= norm
    norm = 1.0

  if (norm > tolerance) :  
    msg = "{0:7s}  error {1:4} wrt norm {2:7}\n".format(name,err,norm)
  else :
    msg = "{0:7s}           error {1:4}\n"               .format(name,err)
  pass_flag = False
  if min_same_rows >= 0:
    pass_flag = err < tolerance or nsame_rows >= min_same_rows or nsame_rows == nrows;
  else:
    pass_flag = err < tolerance
  return [pass_flag,msg];

#====================================================
### run and time tests
#====================================================
def execute(test):
  global tolerance,error_norm,relative_error,custom_flag
  global launch,descriptors
  msg = ""
  os.chdir(test[0])
  tolerance = default_tolerance
  error_norm = default_error_norm
  relative_error = default_relative_error
  launch = [default_lmps]
  descriptors = [default_descriptor]
  options = custom_file+"."+test[1]
  custom_flag = False
  if os.path.isfile(os.path.join(test[0],options)):
    custom_flag = True
    exec(open(options).read(),globals())
  for i,lmp_args in enumerate(launch):
    start()
    msg += run_test(test[1],lmp_args,descriptors[i])
    elapsed_time = stop()
    msg += "elapsed time = "+str(elapsed_time)+" s for test "+test[1]+"\n"
    msg += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
  os.chdir(home)
  return msg

#====================================================
### parse inputs
#====================================================
def init() :
  global default_descriptor, ntests, default_lmps, home
  global auto_rebless_flag, min_same_rows, custom_file
  global default_error_norm, default_relative_error, default_tolerance
  global logread

  # parse input arguments
  if (len(sys.argv) < 4):
    print usage
    sys.exit(1)
  default_descriptor = sys.argv[1]
  default_lmps = sys.argv[2]
  top_dir = os.path.abspath(sys.argv[3])
  dirs = [name for name in os.listdir(top_dir) if os.path.isdir(os.path.join(top_dir, name))]
  tests = []
  exclude_dirs = []
  only_dirs = []
  only_files = []
  os.chdir(top_dir)
  home = os.getcwd()
  cnt = 4
  keywords = ["-auto-rebless","-min-same-rows","-exclude","-only",
              "-customonly","-custom","-error_norm",
              "-relative_error","-tolerance","-logread"]
  while (cnt < len(sys.argv)): 
    option = sys.argv[cnt]
    if ("-auto-rebless" == option):
      flag = sys.argv[cnt+1]
      if (flag == "True"):
        auto_rebless_flag = True
      elif (flag == "False"):
        auto_rebless_flag = False
      else:
        raise Exception("Invalid optional arguements for regression.py")
      cnt += 2
    elif ("-min-same-rows" == option):
      min_same_rows = int(sys.argv[cnt+1])
      cnt += 2
    elif ("-exclude" == option):
      cnt += 1
      while sys.argv[cnt] not in keywords:
        names = [name for name in glob(os.path.join(top_dir, sys.argv[cnt])) if os.path.isdir(name)]
        if len(names) == 0:
          raise Exception("Directory " + sys.argv[cnt] + " not found")
        for name in names:
          exclude_dirs.append(name)
        cnt += 1
        if cnt == len(sys.argv): break
    elif ("-only" == option):
      cnt += 1
      while sys.argv[cnt] not in keywords:
        names = [name for name in glob(os.path.join(top_dir, sys.argv[cnt])) if os.path.isdir(name)]
        if len(names) == 0:
          raise Exception("Directory " + sys.argv[cnt] + " not found")
        for name in names:
          only_dirs.append(name)
        cnt += 1
        if cnt == len(sys.argv): break
    elif ("-customonly" == option):
      cnt += 1
      while sys.argv[cnt] not in keywords:
        only_files.append(sys.argv[cnt])
        cnt += 1
        if cnt == len(sys.argv): break
    elif ("-custom" == option):
      custom_file = sys.argv[cnt+1]
      cnt += 2
    elif ("-error_norm" == option):
      default_error_norm = sys.argv[cnt+1]
      cnt += 2
    elif ("-relative_error" == option):
      flag = sys.argv[cnt+1]
      if (flag == "True"):
        default_relative_error = True
      elif (flag == "False"):
        default_relative_error = False
      else:
        raise Exception("Invalid optional arguements for regression.py")
      cnt += 2
    elif ("-tolerance" == option):
      default_tolerance = float(sys.argv[cnt+1])
      cnt += 2
    elif ("-logread" == option):
      logread = [os.path.abspath(sys.argv[cnt+1]),sys.argv[cnt+2]]
      cnt += 3
    else:
      raise Exception("Invalid optional arguements for regression.py")

  # recursively loop through directories to get tests
  for x in os.walk(top_dir):
    dir = x[0];
    if ".svn" in dir: continue

    # exclude these directories
    if len(exclude_dirs) > 0:
      skip = False
      for i in exclude_dirs:
        for y in os.walk(i):
          if y[0] == dir: skip = True
      if skip: continue

    # only include these directories
    if len(only_dirs) > 0:
      skip = True
      for i in only_dirs:
        for y in os.walk(i):
          if y[0] == dir: skip = False
      if skip: continue

    # only include directories with these files
    if len(only_files) > 0:
      skip = True
      for name in only_files:
        for path in glob(os.path.join(dir,name)):
          if os.path.isfile(path):
            skip = False
      if skip: continue

    os.chdir(dir);
    for path in glob("./in.*"):
      test = path[5:];
      tests.append([dir,test])
    os.chdir(home)
  ntests = len(tests)
  #tests.sort()

  # print header
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "start: ",date()
  print "ntests:",ntests
  print "default descriptor:",default_descriptor
  print "default norm:",default_error_norm
  print "default tolerance:",default_tolerance
  print "default relative error:",default_relative_error
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print
  print "subdirs =",dirs
  print
  sys.stdout.flush()

  return tests

#====================================================
### main
#====================================================
if __name__ == '__main__':
  tests = init()

  if logread[0] != ".": sys.path.append(logread[0])
  strcmd = "from %s import %s as logreader" % (logread[1],logread[1])
  exec strcmd

  nfails = 0
  fail_list = []
  nwarnings = 0
  warn_list = []

  # run the tests

  for test in tests:
    msg = execute(test)
    if (fail_pattern.search(msg)) : 
      nfails += 1
      fail_list.append(test)
    elif (warn_pattern.search(msg)) : 
      nwarnings += 1
      warn_list.append(test)
    print msg
    sys.stdout.flush()

  # print out results

  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "end:",date()
  if (nfails == 0):
    print ntests,"tests passed"
    print "*** no failures ***"
  else:
    print "!!!",nfails,"of",ntests,"tests failed"
    for test in fail_list:
      print test
  if (nwarnings > 0):
    print "\n!!! Warnings were generated in the following tests"
    for test in warn_list:
      print test
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  sys.stdout.flush()
