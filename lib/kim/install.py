#!usr/local/python

# install.py tool to setup the kim-api library
# used to automate the steps described in the README file in this dir

import sys,os,re,urllib,commands

help = """
Syntax: install.py -v version
        specify one or more options, order does not matter
        -v = version of kim-api to download and work with
             default = kim-api-v1.8.2 (current as of June 2017)
"""

def error():
  print help
  sys.exit()

# parse args

args = sys.argv

dir = "."
version = "kim-api-v1.8.2"

iarg = 1
while iarg < len(args):
  if args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2
  else: error()

dir = os.path.abspath(dir)
url = "https://s3.openkim.org/kim-api/%s.tgz" % version

# download and unpack tarball

print "Downloading kim-api tarball ..."
urllib.urlretrieve(url,"%s/%s.tgz" % (dir,version))
print "Unpacking kim-api tarball ..."
cmd = "cd %s; rm -rf %s; tar zxvf %s.tgz" % (dir,version,version)
txt = commands.getstatusoutput(cmd)
if txt[0] != 0: error()

# configure kim-api
print "Configuring kim-api ..."
cmd = "cd %s/%s; ./configure --prefix='%s'" % (dir,version,dir)
txt = commands.getstatusoutput(cmd)
print txt[1]
if txt[0] != 0: error()

# build kim-api

print "Building kim-api ..."
cmd = "cd %s/%s; make" % (dir,version)
txt = commands.getstatusoutput(cmd)
print txt[1]
if txt[0] != 0: error()

# install kim-api

print "Installing kim-api ..."
cmd = "cd %s/%s; make install" % (dir,version)
txt = commands.getstatusoutput(cmd)
print txt[1]
if txt[0] != 0: error()

cmd = "cd %s/%s; make install-set-default-to-v1" %(dir,version)
txt = commands.getstatusoutput(cmd)
print txt[1]
if txt[0] != 0: error()
