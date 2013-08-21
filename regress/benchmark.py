#!/usr/bin/env python 
"""
  function: numerical comparisions of logs and corresponding benchmarks
  usage: benchmark.py <nprocs> <njobs> <dirs>
"""
import sys
import os
import math
import re
from operator import itemgetter
from glob import glob
import time
import multiprocessing as mp
try:
  import Queue as queue # 2.6
except ImportError:
  import queue # 3.0

#====================================================
### constants
#====================================================
thermo_pattern  = re.compile("^Step "); # fragile
data_pattern    = re.compile("\s*\d"); # fragile
fail_pattern    = re.compile("FAIL"); 
tol = 1.e-6 # 1.e-10
arch = "openmpi"
src_path = "../src/" #relative to home
exe_path = "../"+src_path

#====================================================
### date
#====================================================
def date():
 return time.asctime()

#====================================================
### timer
#====================================================
## NOTE these don't seem to work how I expect them to
def start():
  global dt
  dt = -(time.clock())
def stop():
  global dt
  dt += (time.clock())
  return dt

#====================================================
### run a benchmark
#====================================================
def run_test(test):
  input = "in."+test;
  log = "log."+test
  stdout = "stdout."+test
  ref = (glob(log+"*."+str(np)))[0];
  msg = "==== comparing "+log+" with "+ref+" ====\n"
  if (os.path.isfile(log)): os.remove(log)
  if (os.path.isfile(stdout)): os.remove(stdout)
  os.system(lmps+input+" >& "+stdout);
  if (not os.path.isfile(log)) :
    msg += "!!! no "+log+"\n";
    msg += "!!! test "+test+" FAILED\n"
    return msg
  [cdict,cdata] = extract_data(log);
  [bdict,bdata] = extract_data(ref);
  cols = range(len(bdict)) 
  if (len(cdata) != len(bdata)):
    msg += "!!! data size "+str(len(cdata))+" does not match data "+str(len(bdata))+" in "+ref+"\n";
    msg += "!!! test "+test+" FAILED\n"
    return msg
  fail = False 
  i = 0
  for name in bdict:
    [passing,cmsg] = compare(name,cdata[cols[i]],bdata[cols[i]]);
    i += 1
    msg += cmsg
    if (not passing) : fail = True
  if (fail) :
    msg += "!!! test "+test+" FAILED\n"
  else : 
    msg += "*** test "+test+" passed\n"
  return msg

#====================================================
### extract data from log file
#====================================================
def extract_data(file):
  dictionary = [];
  data = []
  read = False
  for line in open(file):
    if (read and data_pattern.match(line)) :
      cols = line.split();
      data.append(cols)
    else :
      read = False
    if (thermo_pattern.match(line)):
      dictionary = line.split();
      read = True
  return [dictionary,data]

#====================================================
### compare columns of current and benchmark
#====================================================
def compare(name,col1,col2):
  err = 0.
  norm1 = 0.
  norm2 = 0.
  n = len(col2)
  for i in range(n):
    v1 = float(col1[i])
    v2 = float(col2[i])
    norm1 += v1*v1
    norm2 += v2*v2
    dv = v1-v2
    err   += dv*dv
  norm1 /= n
  norm2 /= n
  err   /= n
  if (norm2 > tol) :  
    msg = "{0:7s}  relative error {1:4} wrt norm {2:7}\n".format(name,err,norm2)
  else :
    msg = "{0:7s}           error {1:4}\n"               .format(name,err)
  return [(err < tol),msg];

#################################################################
class Worker(mp.Process):
  def __init__(self, work_queue, result_queue):
    mp.Process.__init__(self)
    self.work_queue = work_queue
    self.result_queue = result_queue
  def run(self):
    while True:
      try:
        job = self.work_queue.get_nowait()
      except queue.Empty:
        break
      #print(">>> starting " + str(job[1]) + " ...")
      os.chdir(job[0])
      start()
      msg = run_test(job[1])
      elapsed_time = stop()
      msg += "elapsed time "+str(elapsed_time)+"\n"
      os.chdir(home)
      self.result_queue.put([job[1],msg])

#====================================================
### parse
#====================================================
def init() :
  global np, njobs, ntests, lmps, arch, home
  home = os.getcwd()
  if (len(sys.argv) < 4) :
    print "usage: benchmark.py <nprocs> <njobs> <test_dirs>"
    sys.exit(1)
  np = int(sys.argv[1])
  njobs = int(sys.argv[2])
  lmps = "../"+src_path+"lmp_"+arch+" -in "
  if (np > 1): 
    lmps = "mpirun -np "+str(np)+" "+lmps
  else:
    arch = "serial"
    lmps = exe_path+"lmp_"+arch+" -in "
  pool = mp.Pool(njobs)
  dirs = sys.argv[3:]
  tests = []
  for dir in dirs:
    os.chdir(dir);
    for path in glob("./in.*"):
      test = path[5:];
      tests.append([dir,test])
    os.chdir(home)
  ntests = len(tests)
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "start: ",date()
  print "arch:",arch,
  print "nprocs:",np
  print "ntests:",ntests,
  print "njobs:",njobs
  print "relative tolerance:",tol
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print
  return tests

#====================================================
### build executable
#====================================================
def build(arch):
  os.system("cd ..; svn update >& svn_update.log")
  os.system("cd ../src; make no-atc >& /dev/null")
  os.system("cd ../src; make clean-all >& /dev/null")
  #os.system("cd ../src; make yes-all >& /dev/null")
  os.system("cd ../src; make yes-dipole >& /dev/null")
  sys.stdout.flush()
  print "** building ",arch,"...",
  os.system("cd "+src_path+"; make -j "+str(np)+" "+arch+" >& build_"+arch+".log")
  if (not os.path.isfile(src_path+"lmp_"+arch)) :
    print "!!! build ",arch," FAILED"
    sys.exit()
  else:
    print "done"
  print

#====================================================
### main
#====================================================
if __name__ == '__main__':
  tests = init()
  build(arch)
  work_queue = mp.Queue()
  for test in tests:
    work_queue.put(test)
  result_queue = mp.Queue()
  nfails = 0
  fail_list = []
  for i in range(njobs):
    w = Worker(work_queue, result_queue)
    w.start()
  for i in range(ntests):
    [test,msg] = result_queue.get()
    if (fail_pattern.search(msg)) : 
      nfails += 1
      fail_list.append(test)
      #print msg # can print only if failed
    print msg # can print only if failed
    #print test, " passed"
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  print "end:",date()
  if (nfails == 0):
    print "*** no failures ***"
  else :
    print "!!!",nfails,"of",ntests,"tests failed"
    for test in fail_list:
      print test
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
