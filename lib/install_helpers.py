import hashlib,os,subprocess,sys

# default help message

defhelp = """
Syntax from src dir: make lib-libname args="-m machine -e suffix"
Syntax from lib dir: python Install.py -m machine -e suffix

libname = name of lib dir (e.g. atc, h5md, meam, poems, etc)
specify -m and optionally -e, order does not matter

  -m = peform a clean followed by "make -f Makefile.machine"
       machine = suffix of a lib/Makefile.* file
  -e = set EXTRAMAKE variable in Makefile.machine to Makefile.lammps.suffix
       does not alter existing Makefile.machine

Examples:

make lib-poems args="-m serial" # build POEMS lib with same settings as in the serial Makefile in src
make lib-colvars args="-m mpi"  # build USER-COLVARS lib with same settings as in the mpi Makefile in src
make lib-meam args="-m ifort"   # build MEAM lib with custom Makefile.ifort (using Intel Fortran)
"""

# print error message or help
def error(str=None,help=None):
  if not str:
    if not help:
        print(defhelp)
    else:
        print(help)
  else: print("ERROR",str)
  sys.exit()

# try to auto-detect the maximum number of available CPUs
def get_cpus():
  try:
    import multiprocessing
    n_cpus = multiprocessing.cpu_count()
  except:
    n_cpus = 1
  return n_cpus

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

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

def geturl(url,fname):
  success = False

  if which('curl') != None:
    cmd = 'curl -L -o "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling curl failed with: %s" % e.output.decode('UTF-8'))

  if not success and which('wget') != None:
    cmd = 'wget -O "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    error("Failed to download source code with 'curl' or 'wget'")
  return

def checkmd5sum(md5sum,fname):
    with open(fname,'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(81920)
            if not data:
                break
            m.update(data)
    fh.close()
    return m.hexdigest() == md5sum

