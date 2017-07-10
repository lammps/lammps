#!usr/local/python

# install.py tool to setup the kim-api library
# used to automate the steps described in the README file in this dir

import sys,os,re,urllib,commands

help = """
Syntax: install.py -v version -c kim-dir -b kim-model-name -a kim-name
        specify one or more options, order does not matter
        -v = version of kim-api to download and work with
             default = kim-api-v1.8.2 (current as of June 2017)
        -c = create Makefile.KIM_DIR within lammps lib/kim to configure lammps
             for use with the kim-api library installed at "kim-dir" (absolute
             path).  default = this dir
        -b = build kim-api and kim model where kim-model-name can be a specific
             openkim.org model name (such as
             "EAM_Dynamo_Ackland_W__MO_141627196590_002") or the keyword
             "OpenKIM" to install all compatible models from the openkim.org
             site.
        -a = add kim-name openkim.org item (model driver or model) to existing
             kim-api instalation.
"""

def error():
  print help
  sys.exit()

# parse args

args = sys.argv

thisdir = os.environ['PWD']
dir = thisdir
version = "kim-api-v1.8.2"

dirflag = 0
buildflag = 0
addflag = 0

iarg = 1
while iarg < len(args):
  if args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-c":
    dirflag = 1
    if iarg+2 > len(args): error()
    dir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = 1
    if iarg+2 > len(args): error()
    modelname = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-a":
    addflag = 1
    if iarg+2 > len(args): error()
    addmodelname = args[iarg+1]
    iarg += 2
  else: error()

thisdir = os.path.abspath(thisdir)
dir = os.path.abspath(dir)
url = "https://s3.openkim.org/kim-api/%s.tgz" % version

# download and unpack tarball


if not os.path.isfile("%s/Makefile.KIM_DIR" % thisdir):
  open("%s/Makefile.KIM_DIR" % thisdir, 'w').write("KIM_INSTALL_DIR=%s" % dir)
  open("%s/Makefile.KIM_Config" % thisdir, 'w').write("include %s/lib/kim-api/Makefile.KIM_Config" % dir)
  print "Created %s/Makefile.KIM_DIR : using %s" % (thisdir,dir)
else:
  if dirflag == 1:
    open("%s/Makefile.KIM_DIR" % thisdir, 'w').write("KIM_INSTALL_DIR=%s" % dir)
    open("%s/Makefile.KIM_Config" % thisdir, 'w').write("include %s/lib/kim-api/Makefile.KIM_Config" % dir)
    print "Updated %s/Makefile.KIM_DIR : using %s" % (thisdir,dir)


if buildflag == 1:
  # download kim-api
  print "Downloading kim-api tarball ..."
  urllib.urlretrieve(url,"%s/%s.tgz" % (thisdir,version))
  print "Unpacking kim-api tarball ..."
  cmd = "cd %s; rm -rf %s; tar zxvf %s.tgz" % (thisdir,version,version)
  txt = commands.getstatusoutput(cmd)
  if txt[0] != 0: error()

  # configure kim-api
  print "Configuring kim-api ..."
  cmd = "cd %s/%s; ./configure --prefix='%s'" % (thisdir,version,dir)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()

  # build kim-api
  print "Configuring model : %s" % modelname
  cmd = "cd %s/%s; make add-%s" % (thisdir,version,modelname)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()
  #
  print "Building kim-api ..."
  cmd = "cd %s/%s; make" % (thisdir,version)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()

  # install kim-api
  print "Installing kim-api ..."
  cmd = "cd %s/%s; make install" % (thisdir,version)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()
  #
  cmd = "cd %s/%s; make install-set-default-to-v1" %(thisdir,version)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()

  # remove source files
  print "Removing kim-api source and build files ..."
  cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" % (thisdir,version,version)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()

if addflag == 1:
  # download model
  url = "https://openkim.org/download/%s.tgz" % addmodelname
  print "Downloading item tarball ..."
  urllib.urlretrieve(url,"%s/%s.tgz" % (thisdir,addmodelname))
  print "Unpacking item tarball ..."
  cmd = "cd %s; tar zxvf %s.tgz" % (thisdir,addmodelname)
  txt = commands.getstatusoutput(cmd)
  if txt[0] != 0: error()
  #
  print "Building item ..."
  cmd = "cd %s/%s; make; make install" %(thisdir,addmodelname)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()
  #
  print "Removing kim item source and build files ..."
  cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" %(thisdir,addmodelname,addmodelname)
  txt = commands.getstatusoutput(cmd)
  print txt[1]
  if txt[0] != 0: error()

