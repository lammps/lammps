import hashlib
import os
import re
import subprocess

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

def geturl(url, fname):
  success = False

  if which('curl') != None:
    cmd = 'curl -L -o "%s" %s' % (fname, url)
    try:
      subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling curl failed with: %s" % e.output.decode('UTF-8'))

  if not success and which('wget') != None:
    cmd = 'wget -O "%s" %s' % (fname, url)
    try:
      subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    raise Exception("Failed to download source code with 'curl' or 'wget' from " + url)
  return

def checkmd5sum(md5sum, fname):
    with open(fname, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(81920)
            if not data:
                break
            m.update(data)
    fh.close()
    return m.hexdigest() == md5sum

def getfallback(lib, url):
  archive = re.sub(r'^https://.*/([^/]+gz)', r'-\1', url)
  return 'https://download.lammps.org/thirdparty/' + lib + archive
