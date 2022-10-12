# Pizza.py toolkit, https://lammps.github.io/pizza
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# for python3 compatibility
from __future__ import print_function

# gnu tool

oneline = "Create plots via GnuPlot plotting program"

docstr = """
g = gnu()              start up GnuPlot
g.stop()               shut down GnuPlot process

g.plot(a)                      plot vector A against linear index
g.plot(a,b)            plot B against A
g.plot(a,b,c,d,...)        plot B against A, D against C, etc
g.mplot(M,N,S,"file",a,b,...)  multiple plots saved to file0000.eps, etc

  each plot argument can be a tuple, list, or Numeric/NumPy vector
  mplot loops over range(M,N,S) and create one plot per iteration
    last args are same as list of vectors for plot(), e.g. 1, 2, 4 vectors
    each plot is made from a portion of the vectors, depending on loop index i
      Ith plot is of b[0:i] vs a[0:i], etc
    series of plots saved as file0000.eps, file0001.eps, etc
    if use xrange(),yrange() then plot axes will be same for all plots

g("plot 'file.dat' using 2:3 with lines")      execute string in GnuPlot

g.enter()                       enter GnuPlot shell
gnuplot> plot sin(x) with lines         type commands directly to GnuPlot
gnuplot> exit, quit                 exit GnuPlot shell

g.export("data",range(100),a,...)       create file with columns of numbers

  all vectors must be of equal length
  could plot from file with GnuPlot command: plot 'data' using 1:2 with lines

g.select(N)                    figure N becomes the current plot

  subsequent commands apply to this plot

g.hide(N)                  delete window for figure N
g.save("file")                 save current plot as file.eps

Set attributes for current plot:

g.erase()                      reset all attributes to default values
g.aspect(1.3)                  aspect ratio
g.xtitle("Time")               x axis text
g.ytitle("Energy")             y axis text
g.title("My Plot")             title text
g.title("title","x","y")       title, x axis, y axis text
g.xrange(xmin,xmax)            x axis range
g.xrange()                     default x axis range
g.yrange(ymin,ymax)            y axis range
g.yrange()                     default y axis range
g.xlog()                       toggle x axis between linear and log
g.ylog()                       toggle y axis between linear and log
g.label(x,y,"text")            place label at x,y coords
g.curve(N,'r')                 set color of curve N

  colors: 'k' = black, 'r' = red, 'g' = green, 'b' = blue
          'm' = magenta, 'c' = cyan, 'y' = yellow
"""

# History
#   8/05, Matt Jones (BYU): original version
#   9/05, Steve Plimpton: added mplot() method

# ToDo list
#   allow choice of JPG or PNG or GIF when saving ?
#     can this be done from GnuPlot or have to do via ImageMagick convert ?
#     way to trim EPS plot that is created ?
#   hide does not work on Mac aqua
#   select does not pop window to front on Mac aqua

# Variables
#   current = index of current figure (1-N)
#   figures = list of figure objects with each plot's attributes
#             so they aren't lost between replots

# Imports and external programs

import os
import sys

try: from DEFAULTS import PIZZA_GNUPLOT
except ImportError: PIZZA_GNUPLOT = "gnuplot -p"
try: from DEFAULTS import PIZZA_GNUTERM
except ImportError: PIZZA_GNUTERM = "x11"

# Class definition

class gnu:

  # --------------------------------------------------------------------

  def __init__(self):
    self.GNUPLOT = os.popen(PIZZA_GNUPLOT,'w')
    self.file = "tmp.gnu"
    self.figures = []
    self.select(1)

  # --------------------------------------------------------------------

  def stop(self):
    self.__call__("quit")
    del self.GNUPLOT

  # --------------------------------------------------------------------

  def __call__(self,command):
    self.GNUPLOT.write(command + '\n')
    self.GNUPLOT.flush()

  # --------------------------------------------------------------------

  def enter(self):
    while 1:
      if sys.version_info[0] == 3:
        command = input("gnuplot> ")
      else:
        command = raw_input("gnuplot> ")
      if command == "quit" or command == "exit": return
      self.__call__(command)

  # --------------------------------------------------------------------
  # write plot vectors to files and plot them

  def plot(self,*vectors):
    if len(vectors) == 1:
      file = self.file + ".%d.1" % self.current
      linear = range(len(vectors[0]))
      self.export(file,linear,vectors[0])
      self.figures[self.current-1].ncurves = 1
    else:
      if len(vectors) % 2: raise Exception("vectors must come in pairs")
      for i in range(0,len(vectors),2):
        file = self.file + ".%d.%d" % (self.current,i/2+1)
        self.export(file,vectors[i],vectors[i+1])
      self.figures[self.current-1].ncurves = len(vectors)/2
    self.draw()

  # --------------------------------------------------------------------
  # create multiple plots from growing vectors, save to numbered files
  # don't plot empty vector, create a [0] instead

  def mplot(self,start,stop,skip,file,*vectors):
    n = 0
    for i in range(start,stop,skip):
      partial_vecs = []
      for vec in vectors:
        if i: partial_vecs.append(vec[:i])
        else: partial_vecs.append([0])
      self.plot(*partial_vecs)

      if n < 10:     newfile = file + "000" + str(n)
      elif n < 100:  newfile = file + "00" + str(n)
      elif n < 1000: newfile = file + "0" + str(n)
      else:          newfile = file + str(n)

      self.save(newfile)
      n += 1

  # --------------------------------------------------------------------
  # write list of equal-length vectors to filename

  def export(self,filename,*vectors):
    n = len(vectors[0])
    for vector in vectors:
      if len(vector) != n: raise Exception("vectors must be same length")
    f = open(filename,'w')
    nvec = len(vectors)
    for i in range(n):
      for j in range(nvec):
        print(str(vectors[j][i])+" ",file=f,end='')
      print ("",file=f)
    f.close()

  # --------------------------------------------------------------------
  # select plot N as current plot

  def select(self,n):
    self.current = n
    if len(self.figures) < n:
      for i in range(n - len(self.figures)):
        self.figures.append(figure())
    cmd = "set term " + PIZZA_GNUTERM + ' ' + str(n)
    self.__call__(cmd)
    if self.figures[n-1].ncurves: self.draw()

  # --------------------------------------------------------------------
  # delete window for plot N

  def hide(self,n):
    cmd = "set term %s close %d" % (PIZZA_GNUTERM,n)
    self.__call__(cmd)

  # --------------------------------------------------------------------
  # save plot to file.eps
  # final re-select will reset terminal
  # do not continue until plot file is written out
  #   else script could go forward and change data file
  #   use tmp.done as semaphore to indicate plot is finished

  def save(self,file):
    self.__call__("set terminal postscript enhanced solid lw 2 color portrait")
    cmd = "set output '%s.eps'" % file
    self.__call__(cmd)
    if os.path.exists("tmp.done"): os.remove("tmp.done")
    self.draw()
    self.__call__("!touch tmp.done")
    while not os.path.exists("tmp.done"): continue
    self.__call__("set output")
    self.select(self.current)

  # --------------------------------------------------------------------
  # restore default attributes by creating a new fig object

  def erase(self):
    fig = figure()
    fig.ncurves = self.figures[self.current-1].ncurves
    self.figures[self.current-1] = fig
    self.draw()

  # --------------------------------------------------------------------

  def aspect(self,value):
    self.figures[self.current-1].aspect = value
    self.draw()

  # --------------------------------------------------------------------

  def xrange(self,*values):
    if len(values) == 0:
      self.figures[self.current-1].xlimit = 0
    else:
      self.figures[self.current-1].xlimit = (values[0],values[1])
    self.draw()

  # --------------------------------------------------------------------

  def yrange(self,*values):
    if len(values) == 0:
      self.figures[self.current-1].ylimit = 0
    else:
      self.figures[self.current-1].ylimit = (values[0],values[1])
    self.draw()

  # --------------------------------------------------------------------

  def label(self,x,y,text):
    self.figures[self.current-1].labels.append((x,y,text))
    self.figures[self.current-1].nlabels += 1
    self.draw()

  # --------------------------------------------------------------------

  def nolabels(self):
    self.figures[self.current-1].nlabel = 0
    self.figures[self.current-1].labels = []
    self.draw()

  # --------------------------------------------------------------------

  def title(self,*strings):
    if len(strings) == 1:
      self.figures[self.current-1].title = strings[0]
    else:
      self.figures[self.current-1].title = strings[0]
      self.figures[self.current-1].xtitle = strings[1]
      self.figures[self.current-1].ytitle = strings[2]
    self.draw()

  # --------------------------------------------------------------------

  def xtitle(self,label):
    self.figures[self.current-1].xtitle = label
    self.draw()

  # --------------------------------------------------------------------

  def ytitle(self,label):
    self.figures[self.current-1].ytitle = label
    self.draw()

  # --------------------------------------------------------------------

  def xlog(self):
    if self.figures[self.current-1].xlog:
      self.figures[self.current-1].xlog = 0
    else:
      self.figures[self.current-1].xlog = 1
    self.draw()

  # --------------------------------------------------------------------

  def ylog(self):
    if self.figures[self.current-1].ylog:
      self.figures[self.current-1].ylog = 0
    else:
      self.figures[self.current-1].ylog = 1
    self.draw()

  # --------------------------------------------------------------------

  def curve(self,num,color):
    fig = self.figures[self.current-1]
    while len(fig.colors) < num: fig.colors.append(0)
    fig.colors[num-1] = colormap[color]
    self.draw()

  # --------------------------------------------------------------------
  # draw a plot with all its settings
  # just return if no files of vectors defined yet

  def draw(self):
    fig = self.figures[self.current-1]
    if not fig.ncurves: return

    cmd = 'set size ratio ' + str(1.0/float(fig.aspect))
    self.__call__(cmd)

    cmd = 'set title ' + '"' + fig.title + '"'
    self.__call__(cmd)
    cmd = 'set xlabel ' + '"' + fig.xtitle + '"'
    self.__call__(cmd)
    cmd = 'set ylabel ' + '"' + fig.ytitle + '"'
    self.__call__(cmd)

    if fig.xlog: self.__call__("set logscale x")
    else: self.__call__("unset logscale x")
    if fig.ylog: self.__call__("set logscale y")
    else: self.__call__("unset logscale y")
    if fig.xlimit:
      cmd = 'set xr [' + str(fig.xlimit[0]) + ':' + str(fig.xlimit[1]) + ']'
      self.__call__(cmd)
    else: self.__call__("set xr [*:*]")
    if fig.ylimit:
      cmd = 'set yr [' + str(fig.ylimit[0]) + ':' + str(fig.ylimit[1]) + ']'
      self.__call__(cmd)
    else: self.__call__("set yr [*:*]")

    self.__call__("set nolabel")
    for i in range(fig.nlabels):
      x = fig.labels[i][0]
      y = fig.labels[i][1]
      text = fig.labels[i][2]
      cmd = 'set label ' + '\"' + text + '\" at ' + str(x) + ',' + str(y)
      self.__call__(cmd)

    self.__call__("set key off")
    cmd = 'plot '
    for i in range(int(fig.ncurves)):
      file = self.file + ".%d.%d" % (self.current,i+1)
      if len(fig.colors) > i and fig.colors[i]:
        cmd += "'" + file + "' using 1:2 with line %d, " % fig.colors[i]
      else:
        cmd += "'" + file + "' using 1:2 with lines, "
    self.__call__(cmd[:-2])

# --------------------------------------------------------------------
# class to store settings for a single plot

class figure:

  def __init__(self):
    self.ncurves = 0
    self.colors  = []
    self.title   = ""
    self.xtitle  = ""
    self.ytitle  = ""
    self.aspect  = 1.3
    self.xlimit  = 0
    self.ylimit  = 0
    self.xlog    = 0
    self.ylog    = 0
    self.nlabels = 0
    self.labels  = []

# --------------------------------------------------------------------
# line color settings

colormap = {'k':-1, 'r':1, 'g':2, 'b':3, 'm':4, 'c':5, 'y':7}
