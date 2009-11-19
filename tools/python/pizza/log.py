# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# log tool

oneline = "Read LAMMPS log files and extract thermodynamic data"

docstr = """
l = log("file1")                     read in one or more log files
l = log("log1 log2.gz")              can be gzipped
l = log("file*")                     wildcard expands to multiple files
l = log("log.lammps",0)              two args = store filename, but don't read

  incomplete and duplicate thermo entries are deleted

time = l.next()                      read new thermo info from file

  used with 2-argument constructor to allow reading thermo incrementally
  return time stamp of last thermo read
  return -1 if no new thermo since last read

nvec = l.nvec                        # of vectors of thermo info
nlen = l.nlen                        length of each vectors
names = l.names                      list of vector names
t,pe,... = l.get("Time","KE",...)    return one or more vectors of values
l.write("file.txt")	 	     write all vectors to a file
l.write("file.txt","Time","PE",...)  write listed vectors to a file

  get and write allow abbreviated (uniquely) vector names
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   nvec = # of vectors
#   nlen = length of each vector
#   names = list of vector names
#   ptr = dictionary, key = name, value = index into data for which column
#   data[i][j] = 2d array of floats, i = 0 to # of entries, j = 0 to nvecs-1
#   style = style of LAMMPS log file, 1 = multi, 2 = one, 3 = gran
#   firststr = string that begins a thermo section in log file
#   increment = 1 if log file being read incrementally
#   eof = ptr into incremental file for where to start next read

# Imports and external programs

import sys, re, glob
from os import popen

try: tmp = PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class log:

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.nvec = 0
    self.names = []
    self.ptr = {}
    self.data = []

    # flist = list of all log file names

    words = list[0].split()
    self.flist = []
    for word in words: self.flist += glob.glob(word)
    if len(self.flist) == 0 and len(list) == 1:
      raise StandardError,"no log file specified"

    if len(list) == 1:
      self.increment = 0
      self.read_all()
    else:
      if len(self.flist) > 1:
        raise StandardError,"can only incrementally read one log file"
      self.increment = 1
      self.eof = 0

  # --------------------------------------------------------------------
  # read all thermo from all files
  
  def read_all(self):
    self.read_header(self.flist[0])
    if self.nvec == 0: raise StandardError,"log file has no values"

    # read all files

    for file in self.flist: self.read_one(file)
    print

    # sort entries by timestep, cull duplicates
    
    self.data.sort(self.compare)
    self.cull()
    self.nlen = len(self.data)
    print "read %d log entries" % self.nlen

  # --------------------------------------------------------------------

  def next(self):
    if not self.increment: raise StandardError,"cannot read incrementally"

    if self.nvec == 0:
      try: open(self.flist[0],'r')
      except: return -1
      self.read_header(self.flist[0])
      if self.nvec == 0: return -1

    self.eof = self.read_one(self.flist[0],self.eof)
    return int(self.data[-1][0])

  # --------------------------------------------------------------------

  def get(self,*keys):
    if len(keys) == 0:
      raise StandardError, "no log vectors specified"

    map = []
    for key in keys:
      if self.ptr.has_key(key):
        map.append(self.ptr[key])
      else:
        count = 0
        for i in range(self.nvec):
	  if self.names[i].find(key) == 0:
	    count += 1
	    index = i
        if count == 1:
          map.append(index)
        else:
          raise StandardError, "unique log vector %s not found" % key

    vecs = []
    for i in range(len(keys)):
      vecs.append(self.nlen * [0])
      for j in xrange(self.nlen):
        vecs[i][j] = self.data[j][map[i]]

    if len(keys) == 1: return vecs[0]
    else: return vecs

  # --------------------------------------------------------------------

  def write(self,filename,*keys):
    if len(keys):
      map = []
      for key in keys:
        if self.ptr.has_key(key):
          map.append(self.ptr[key])
        else:
          count = 0
          for i in range(self.nvec):
	    if self.names[i].find(key) == 0:
	      count += 1
	      index = i
          if count == 1:
            map.append(index)
          else:
            raise StandardError, "unique log vector %s not found" % key
    else:
      map = range(self.nvec)

    f = open(filename,"w")
    for i in xrange(self.nlen):
      for j in xrange(len(map)):
        print >>f,self.data[i][map[j]],
      print >>f
    f.close()

  # --------------------------------------------------------------------

  def compare(self,a,b):
    if a[0] < b[0]:
      return -1
    elif a[0] > b[0]:
      return 1
    else:
      return 0

  # --------------------------------------------------------------------

  def cull(self):
    i = 1
    while i < len(self.data):
      if self.data[i][0] == self.data[i-1][0]: del self.data[i]
      else: i += 1

  # --------------------------------------------------------------------

  def read_header(self,file):
    str_multi = "----- Step"
    str_one = "Step "

    if file[-3:] == ".gz":
      txt = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r').read()
    else:
      txt = open(file).read()

    if txt.find(str_multi) >= 0:
      self.firststr = str_multi
      self.style = 1
    elif txt.find(str_one) >= 0:
      self.firststr = str_one
      self.style = 2
    else:
      return

    if self.style == 1:
      s1 = txt.find(self.firststr)
      s2 = txt.find("\n--",s1)
      pattern = "\s(\S*)\s*="
      keywords = re.findall(pattern,txt[s1:s2])
      keywords.insert(0,"Step")
      i = 0
      for keyword in keywords:
	self.names.append(keyword)
        self.ptr[keyword] = i
        i += 1

    else:
      s1 = txt.find(self.firststr)
      s2 = txt.find("\n",s1)
      line = txt[s1:s2]
      words = line.split()
      for i in range(len(words)):
	self.names.append(words[i])
        self.ptr[words[i]] = i

    self.nvec = len(self.names)

  # --------------------------------------------------------------------

  def read_one(self,*list):

    # if 2nd arg exists set file ptr to that value
    # read entire (rest of) file into txt

    file = list[0]
    if file[-3:] == ".gz":
      f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'rb')
    else:
      f = open(file,'rb')

    if len(list) == 2: f.seek(list[1])
    txt = f.read()
    if file[-3:] == ".gz": eof = 0
    else: eof = f.tell()
    f.close()

    start = last = 0
    while not last:

      # chunk = contiguous set of thermo entries (line or multi-line)
      # s1 = 1st char on 1st line of chunk
      # s2 = 1st char on line after chunk
      # set last = 1 if this is last chunk in file, leave 0 otherwise
      # set start = position in file to start looking for next chunk
      # rewind eof if final entry is incomplete

      s1 = txt.find(self.firststr,start)
      s2 = txt.find("Loop time of",start+1)

      if s1 >= 0 and s2 >= 0 and s1 < s2:    # found s1,s2 with s1 before s2
        if self.style == 2:
	  s1 = txt.find("\n",s1) + 1
      elif s1 >= 0 and s2 >= 0 and s2 < s1:  # found s1,s2 with s2 before s1
        s1 = 0
      elif s1 == -1 and s2 >= 0:             # found s2, but no s1
	last = 1
        s1 = 0
      elif s1 >= 0 and s2 == -1:             # found s1, but no s2
        last = 1
        if self.style == 1:
          s2 = txt.rfind("\n--",s1) + 1
        else:
	  s1 = txt.find("\n",s1) + 1
          s2 = txt.rfind("\n",s1) + 1
	eof -= len(txt) - s2
      elif s1 == -1 and s2 == -1:            # found neither
                                             # could be end-of-file section
					     # or entire read was one chunk

        if txt.find("Loop time of",start) == start:   # end of file, so exit
	  eof -= len(txt) - start                     # reset eof to "Loop"
	  break

	last = 1                                      # entire read is a chunk
        s1 = 0
        if self.style == 1:
          s2 = txt.rfind("\n--",s1) + 1
        else:
          s2 = txt.rfind("\n",s1) + 1
	eof -= len(txt) - s2
	if s1 == s2: break

      chunk = txt[s1:s2-1]
      start = s2

      # split chunk into entries
      # parse each entry for numeric fields, append to data
  
      if self.style == 1:
        sections = chunk.split("\n--")
        pat1 = re.compile("Step\s*(\S*)\s")
        pat2 = re.compile("=\s*(\S*)")
        for section in sections:
          word1 = [re.search(pat1,section).group(1)]
          word2 = re.findall(pat2,section)
          words = word1 + word2
          self.data.append(map(float,words))

      else:
        lines = chunk.split("\n")
        for line in lines:
          words = line.split()
          self.data.append(map(float,words))

      # print last timestep of chunk

      print int(self.data[len(self.data)-1][0]),
      sys.stdout.flush()

    return eof
