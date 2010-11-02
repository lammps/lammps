# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# gl tool

oneline = "3d interactive visualization via OpenGL"

docstr = """
g = gl(d)                   create OpenGL display for data in d

  d = atom snapshot object (dump, data)

g.bg("black")               set background color (def = "black")
g.size(N)		    set image size to NxN
g.size(N,M)		    set image size to NxM
g.rotate(60,135)            view from z theta and azimuthal phi (def = 60,30)
g.shift(x,y)                translate by x,y pixels in view window (def = 0,0)
g.zoom(0.5)                 scale image by factor (def = 1)
g.box(0/1/2)                0/1/2 = none/variable/fixed box
g.box(0/1/2,"green")        set box color
g.box(0/1/2,"red",4)        set box edge thickness
g.file = "image"            file prefix for created images (def = "image")

g.show(N)                   show image of snapshot at timestep N
  
g.all()                     make images of all selected snapshots
g.all(P)                    images of all, start file label at P
g.all(N,M,P)                make M images of snapshot N, start label at P

g.pan(60,135,1.0,40,135,1.5)    pan during all() operation
g.pan()                         no pan during all() (default)

  args = z theta, azimuthal phi, zoom factor at beginning and end
  values at each step are interpolated between beginning and end values

g.select = "$x > %g*3.0"    string to pass to d.aselect.test() during all()
g.select = ""               no extra aselect (default)
				
  %g varies from 0.0 to 1.0 from beginning to end of all()
  
g.acol(2,"green")		   set atom colors by atom type (1-N)
g.acol([2,4],["red","blue"])	   1st arg = one type or list of types
g.acol(0,"blue")	           2nd arg = one color or list of colors
g.acol(range(20),["red","blue"])   if list lengths unequal, interpolate
g.acol(range(10),"loop")           assign colors in loop, randomly ordered

  if 1st arg is 0, set all types to 2nd arg
  if list of types has a 0 (e.g. range(10)), +1 is added to each value
  interpolate means colors blend smoothly from one value to the next

g.arad([1,2],[0.5,0.3])            set atom radii, same rules as acol()

g.bcol()			   set bond color, same args as acol()
g.brad()			   set bond thickness, same args as arad()

g.tcol()			   set triangle color, same args as acol()
g.tfill()			   set triangle fill, 0 fill, 1 line, 2 both

g.lcol()                           set line color, same args as acol()
g.lrad()                           set line thickness, same args as arad()

g.adef()                           set atom/bond/tri/line properties to default
g.bdef()			   default = "loop" for colors, 0.45 for radii
g.tdef()  			   default = 0.25 for bond/line thickness
g.ldef()  			   default = 0 fill

  by default 100 types are assigned
  if atom/bond/tri/line has type > # defined properties, is an error
  
from vizinfo import colors         access color list
print colors                       list defined color names and RGB values
colors["nickname"] = [R,G,B]       set new RGB values from 0 to 255

  140 pre-defined colors: red, green, blue, purple, yellow, black, white, etc

Settings specific to gl tool:

g.q(10)                     set quality of image (def = 5)
g.axis(0/1)                 turn xyz axes off/on
g.ortho(0/1)                perspective (0) vs orthographic (1) view
g.clip('xlo',0.25)          clip in xyz from lo/hi at box fraction (0-1)
g.reload()                  force all data to be reloaded
g.cache = 0/1               turn off/on GL cache lists (def = on)
theta,phi,x,y,scale,up = g.gview()   grab all current view parameters
g.sview(theta,phi,x,y,scale,up)      set all view parameters

  data reload is necessary if dump selection is used to change the data
  cache lists usually improve graphics performance
  gview returns values to use in other commands:
    theta,phi are args to rotate()
    x,y are args to shift()
    scale is arg to zoom()
    up is a 3-vector arg to sview()
"""

# History
#   9/05, Steve Plimpton (SNL): original version

# ToDo list
#   when do aselect with select str while looping N times on same timestep
#     would not let you grow # of atoms selected

# Variables
#   ztheta = vertical angle from z-azis of viewpoint
#   azphi = azimuthal angle of viewpoint
#   xshift,yshift = xy translation of scene (in pixels)
#   distance = size of simulation box (largest dim)
#   eye = viewpoint distance from center of scene
#   file = filename prefix to use for images produced
#   boxflag = 0/1/2 for drawing simulation box: none/variable/fixed
#   bxcol = color of box
#   bxthick = thickness of box lines
#   bgcol = color of background
#   vizinfo = scene attributes
#   center[3] = center point of simulation box
#   view[3] = direction towards eye in simulation box (unit vector)
#   up[3] = screen up direction in simulation box (unit vector)
#   right[3] = screen right direction in simulation box (unit vector)

# Imports and external programs

from math import sin,cos,sqrt,pi,acos
from OpenGL.Tk import *
from OpenGL.GLUT import *
import Image
from vizinfo import vizinfo

# Class definition

class gl:

# --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data
    self.root = None
    self.xpixels = 512
    self.ypixels = 512
    self.ztheta = 60
    self.azphi = 30
    self.scale = 1.0
    self.xshift = self.yshift = 0
    
    self.file = "image"
    self.boxflag = 0
    self.bxcol = [1,1,0]
    self.bxthick = 0.3
    self.bgcol = [0,0,0]
    self.labels = []
    self.panflag = 0
    self.select = ""

    self.axisflag = 0
    self.orthoflag = 1
    self.nslices = 5
    self.nstacks = 5
    self.nsides = 10
    self.theta_amplify = 2
    self.shiny = 2
    
    self.clipflag = 0
    self.clipxlo = self.clipylo = self.clipzlo = 0.0
    self.clipxhi = self.clipyhi = self.clipzhi = 1.0

    self.nclist = 0
    self.calllist = [0]         # indexed by 1-Ntype, so start with 0 index
    self.cache = 1
    self.cachelist = 0

    self.boxdraw = []
    self.atomdraw = []
    self.bonddraw = []
    self.tridraw = []
    self.linedraw = []

    self.ready = 0
    self.create_window()

    self.vizinfo = vizinfo()
    self.adef()
    self.bdef()
    self.tdef()
    self.ldef()
    
    self.center = 3*[0]
    self.view = 3*[0]
    self.up = 3*[0]
    self.right = 3*[0]
    self.viewupright()

  # --------------------------------------------------------------------

  def bg(self,color):
    from vizinfo import colors
    self.bgcol = [colors[color][0]/255.0,colors[color][1]/255.0,
                  colors[color][2]/255.0]
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def size(self,xnew,ynew=None):
    self.xpixels = xnew
    if not ynew: self.ypixels = self.xpixels
    else: self.ypixels = ynew
    self.create_window()
    
  # --------------------------------------------------------------------

  def axis(self,value):
    self.axisflag = value
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def create_window(self):
    if self.root: self.root.destroy()
    
    from __main__ import tkroot
    self.root = Toplevel(tkroot)
    self.root.title('Pizza.py gl tool')

    self.w = MyOpengl(self.root,width=self.xpixels,height=self.ypixels,
                      double=1,depth=1)
    self.w.pack(expand=YES)
#    self.w.pack(expand=YES,fill=BOTH)
    
    glViewport(0,0,self.xpixels,self.ypixels)
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)

    self.rtrack = self.xpixels
    if self.ypixels > self.xpixels: self.rtrack = self.ypixels

    self.w.redraw = self.redraw
    self.w.parent = self
    self.w.tkRedraw()
    tkroot.update_idletasks()              # force window to appear
    
  # --------------------------------------------------------------------

  def clip(self,which,value):
    if which == "xlo":
      self.clipxlo = value
      if value > self.clipxhi: self.clipxlo = self.clipxhi
    elif which == "xhi":
      self.clipxhi = value
      if value < self.clipxlo: self.clipxhi = self.clipxlo
    elif which == "ylo":
      self.clipylo = value
      if value > self.clipyhi: self.clipylo = self.clipyhi
    elif which == "yhi":
      self.clipyhi = value
      if value < self.clipylo: self.clipyhi = self.clipylo
    elif which == "zlo":
      self.clipzlo = value
      if value > self.clipzhi: self.clipzlo = self.clipzhi
    elif which == "zhi":
      self.clipzhi = value
      if value < self.clipzlo: self.clipzhi = self.clipzlo

    oldflag = self.clipflag
    if self.clipxlo > 0 or self.clipylo > 0 or self.clipzlo > 0 or \
       self.clipxhi < 1 or self.clipyhi < 1 or self.clipzhi < 1:
      self.clipflag = 1
    else: self.clipflag = 0

    if oldflag == 0 and self.clipflag == 0: return
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def q(self,value):
    self.nslices = value
    self.nstacks = value
    self.make_atom_calllist()
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def ortho(self,value):
    self.orthoflag = value
    self.w.tkRedraw()

  # --------------------------------------------------------------------
  # set unit vectors for view,up,right from ztheta,azphi
  # assume +z in scene should be up on screen (unless looking down z-axis)
  # right = up x view

  def viewupright(self):
    self.view[0] = cos(pi*self.azphi/180) * sin(pi*self.ztheta/180)
    self.view[1] = sin(pi*self.azphi/180) * sin(pi*self.ztheta/180)
    self.view[2] = cos(pi*self.ztheta/180)

    if self.ztheta == 0.0:
      self.up[0] = cos(pi*self.azphi/180)
      self.up[1] = -sin(pi*self.azphi/180)
      self.up[2] = 0.0
    elif self.ztheta == 180.0:
      self.up[0] = cos(pi*self.azphi/180)
      self.up[1] = sin(pi*self.azphi/180)
      self.up[2] = 0.0
    else:
      dot = self.view[2]		   # dot = (0,0,1) . view
      self.up[0] = -dot*self.view[0]       # up projected onto v = dot * v
      self.up[1] = -dot*self.view[1]       # up perp to v = up - dot * v
      self.up[2] = 1.0 - dot*self.view[2]

    self.up = vecnorm(self.up)
    self.right = veccross(self.up,self.view)

  # --------------------------------------------------------------------
  # reset ztheta,azphi and thus view,up.right
  # called as function from Pizza.py
  
  def rotate(self,ztheta,azphi):
    self.ztheta = ztheta
    self.azphi = azphi
    self.viewupright()
    self.setview()
    self.w.tkRedraw()

  # --------------------------------------------------------------------
  # return all view params to reproduce current display via sview()

  def gview(self):
    return self.ztheta,self.azphi,self.xshift,self.yshift,self.scale,self.up

  # --------------------------------------------------------------------
  # set current view, called by user with full set of view params
  # up is not settable via any other call, all other params are

  def sview(self,ztheta,azphi,xshift,yshift,scale,up):
    self.ztheta = ztheta
    self.azphi = azphi
    self.xshift = xshift
    self.yshift = yshift
    self.scale = scale
    self.up[0] = up[0]
    self.up[1] = up[1]
    self.up[2] = up[2]
    self.up = vecnorm(self.up)
    self.view[0] = cos(pi*self.azphi/180) * sin(pi*self.ztheta/180)
    self.view[1] = sin(pi*self.azphi/180) * sin(pi*self.ztheta/180)
    self.view[2] = cos(pi*self.ztheta/180)
    self.right = veccross(self.up,self.view)
    self.setview()
    self.w.tkRedraw()

  # --------------------------------------------------------------------
  # rotation triggered by mouse trackball
  # project old,new onto unit trackball surf
  # rotate view,up around axis of rotation = old x new
  # right = up x view
  # reset ztheta,azphi from view
  
  def mouse_rotate(self,xnew,ynew,xold,yold):

    # change y pixels to measure from bottom of window instead of top
    
    yold = self.ypixels - yold
    ynew = self.ypixels - ynew

    # vold = unit vector to (xold,yold) projected onto trackball
    # vnew = unit vector to (xnew,ynew) projected onto trackball
    # return (no rotation) if either projection point is outside rtrack

    vold = [0,0,0]
    vold[0] = xold - (0.5*self.xpixels + self.xshift)
    vold[1] = yold - (0.5*self.ypixels + self.yshift)
    vold[2] = self.rtrack*self.rtrack - vold[0]*vold[0] - vold[1]*vold[1]
    if vold[2] < 0: return
    vold[2] = sqrt(vold[2])
    vold = vecnorm(vold)

    vnew = [0,0,0]
    vnew[0] = xnew - (0.5*self.xpixels + self.xshift)
    vnew[1] = ynew - (0.5*self.ypixels + self.yshift)
    vnew[2] = self.rtrack*self.rtrack - vnew[0]*vnew[0] - vnew[1]*vnew[1]
    if vnew[2] < 0: return
    vnew[2] = sqrt(vnew[2])
    vnew = vecnorm(vnew)

    # rot = trackball rotation axis in screen ref frame = vold x vnew
    # theta = angle of rotation = sin(theta) for small theta
    # axis = rotation axis in body ref frame described by right,up,view

    rot = veccross(vold,vnew)
    theta = sqrt(rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2])
    theta *= self.theta_amplify

    axis = [0,0,0]
    axis[0] = rot[0]*self.right[0] + rot[1]*self.up[0] + rot[2]*self.view[0]
    axis[1] = rot[0]*self.right[1] + rot[1]*self.up[1] + rot[2]*self.view[1]
    axis[2] = rot[0]*self.right[2] + rot[1]*self.up[2] + rot[2]*self.view[2]
    axis = vecnorm(axis)
    
    # view is changed by (axis x view) scaled by theta
    # up is changed by (axis x up) scaled by theta
    # force up to be perp to view via up_perp = up - (up . view) view
    # right = up x view

    delta = veccross(axis,self.view)
    self.view[0] -= theta*delta[0]
    self.view[1] -= theta*delta[1]
    self.view[2] -= theta*delta[2]
    self.view = vecnorm(self.view)

    delta = veccross(axis,self.up)
    self.up[0] -= theta*delta[0]
    self.up[1] -= theta*delta[1]
    self.up[2] -= theta*delta[2]

    dot = vecdot(self.up,self.view)
    self.up[0] -= dot*self.view[0]
    self.up[1] -= dot*self.view[1]
    self.up[2] -= dot*self.view[2]
    self.up = vecnorm(self.up)

    self.right = veccross(self.up,self.view)

    # convert new view to ztheta,azphi

    self.ztheta = acos(self.view[2])/pi * 180.0
    if (self.ztheta == 0.0): self.azphi = 0.0
    else: self.azphi = acos(self.view[0]/sin(pi*self.ztheta/180.0))/pi * 180.0
    if self.view[1] < 0: self.azphi = 360.0 - self.azphi
    self.setview()
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def shift(self,x,y):
    self.xshift = x;
    self.yshift = y;
    self.setview()
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def zoom(self,scale):
    self.scale = scale
    self.setview()
    self.w.tkRedraw()

  # --------------------------------------------------------------------
  # set view params needed by redraw
  # input:  center = center of box
  #         distance = size of scene (longest box length)
  #         scale = zoom factor (1.0 = no zoom)
  #         xshift,yshift = translation factor in pixels
  #         view = unit vector from center to viewpoint
  #         up = unit vector in up direction in scene
  #         right = unit vector in right direction in scene
  # output: eye = distance to view scene from
  #         xto,yto,zto = point to look to
  #         xfrom,yfrom,zfrom = point to look from
  
  def setview(self):
    if not self.ready: return                  # no distance since no scene yet
    
    self.eye = 3 * self.distance / self.scale
    xfactor = 0.5*self.eye*self.xshift/self.xpixels
    yfactor = 0.5*self.eye*self.yshift/self.ypixels
    
    self.xto = self.center[0] - xfactor*self.right[0] - yfactor*self.up[0]
    self.yto = self.center[1] - xfactor*self.right[1] - yfactor*self.up[1]
    self.zto = self.center[2] - xfactor*self.right[2] - yfactor*self.up[2]

    self.xfrom = self.xto + self.eye*self.view[0]
    self.yfrom = self.yto + self.eye*self.view[1]
    self.zfrom = self.zto + self.eye*self.view[2]

  # --------------------------------------------------------------------
  # box attributes, also used for triangle lines
  
  def box(self,*args):
    self.boxflag = args[0]
    if len(args) > 1:
      from vizinfo import colors
      self.bxcol = [colors[args[1]][0]/255.0,colors[args[1]][1]/255.0,
                    colors[args[1]][2]/255.0]
    if len(args) > 2: self.bxthick = args[2]
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------
  # grab all selected snapshots from data object
  # add GL-specific info to each bond
  
  def reload(self):
    print "Loading data into gl tool ..."
    data = self.data

    self.timeframes = []
    self.boxframes = []
    self.atomframes = []
    self.bondframes = []
    self.triframes = []
    self.lineframes = []

    box = []
    if self.boxflag == 2: box = data.maxbox()

    flag = 0
    while 1:
      which,time,flag = data.iterator(flag)
      if flag == -1: break
      time,boxone,atoms,bonds,tris,lines = data.viz(which)
      if self.boxflag < 2: box = boxone
      if bonds: self.bonds_augment(bonds)

      self.timeframes.append(time)
      self.boxframes.append(box)
      self.atomframes.append(atoms)
      self.bondframes.append(bonds)
      self.triframes.append(tris)
      self.lineframes.append(lines)
      
      print time,
      sys.stdout.flush()
    print

    self.nframes = len(self.timeframes)
    self.distance = compute_distance(self.boxframes[0])
    self.center = compute_center(self.boxframes[0])
    self.ready = 1
    self.setview()

  # --------------------------------------------------------------------

  def nolabel(self):
    self.cachelist = -self.cachelist
    self.labels = []
  
  # --------------------------------------------------------------------
  # show a single snapshot
  # distance from snapshot box or max box for all selected steps
  
  def show(self,ntime):
    data = self.data
    which = data.findtime(ntime)
    time,box,atoms,bonds,tris,lines = data.viz(which)
    if self.boxflag == 2: box = data.maxbox()
    self.distance = compute_distance(box)
    self.center = compute_center(box)

    if bonds: self.bonds_augment(bonds)

    self.boxdraw = box
    self.atomdraw = atoms
    self.bonddraw = bonds
    self.tridraw = tris
    self.linedraw = lines

    self.ready = 1
    self.setview()
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
    self.save()
    
  # --------------------------------------------------------------------

  def pan(self,*list):
    if len(list) == 0: self.panflag = 0
    else:
      self.panflag = 1
      self.ztheta_start = list[0]
      self.azphi_start = list[1]
      self.scale_start = list[2]
      self.ztheta_stop = list[3]
      self.azphi_stop = list[4]
      self.scale_stop = list[5]
      
  # --------------------------------------------------------------------

  def all(self,*list):
    data = self.data
    if len(list) == 0:
      nstart = 0
      ncount = data.nselect
    elif len(list) == 1:
      nstart = list[0]
      ncount = data.nselect
    else:
      ntime = list[0]
      nstart = list[2]
      ncount = list[1]

    if self.boxflag == 2: box = data.maxbox()

    # loop over all selected steps
    # distance from 1st snapshot box or max box for all selected steps
    # recompute box center on 1st step or if panning

    if len(list) <= 1:

      n = nstart
      i = flag = 0
      while 1:
        which,time,flag = data.iterator(flag)
        if flag == -1: break

        fraction = float(i) / (ncount-1)
        
        if self.select != "":
          newstr = self.select % fraction
          data.aselect.test(newstr,time)
        time,boxone,atoms,bonds,tris,lines = data.viz(which)

        if self.boxflag < 2: box = boxone
        if n == nstart: self.distance = compute_distance(box)

        if n < 10:     file = self.file + "000" + str(n)
        elif n < 100:  file = self.file + "00" + str(n)
        elif n < 1000: file = self.file + "0" + str(n)
        else:          file = self.file + str(n)

        if self.panflag:
          self.ztheta = self.ztheta_start + \
                        fraction*(self.ztheta_stop - self.ztheta_start)
          self.azphi = self.azphi_start + \
                       fraction*(self.azphi_stop - self.azphi_start)
          self.scale = self.scale_start + \
                          fraction*(self.scale_stop - self.scale_start)
          self.viewupright()

	if n == nstart or self.panflag: self.center = compute_center(box)

        if bonds: self.bonds_augment(bonds)

        self.boxdraw = box
        self.atomdraw = atoms
        self.bonddraw = bonds
        self.tridraw = tris
        self.linedraw = lines

        self.ready = 1
        self.setview()
        self.cachelist = -self.cachelist
        self.w.tkRedraw()
        self.save(file)
        
        print time,
        sys.stdout.flush()
        i += 1
        n += 1

    # loop ncount times on same step
    # distance from 1st snapshot box or max box for all selected steps
    # recompute box center on 1st step or if panning

    else:

      which = data.findtime(ntime)

      n = nstart
      for i in range(ncount):
        fraction = float(i) / (ncount-1)

        if self.select != "":
          newstr = self.select % fraction
          data.aselect.test(newstr,ntime)
        time,boxone,atoms,bonds,tris,lines = data.viz(which)

        if self.boxflag < 2: box = boxone
        if n == nstart: self.distance = compute_distance(box)

        if n < 10:     file = self.file + "000" + str(n)
        elif n < 100:  file = self.file + "00" + str(n)
        elif n < 1000: file = self.file + "0" + str(n)
        else:          file = self.file + str(n)

        if self.panflag:
          self.ztheta = self.ztheta_start + \
                        fraction*(self.ztheta_stop - self.ztheta_start)
          self.azphi = self.azphi_start + \
                       fraction*(self.azphi_stop - self.azphi_start)
          self.scale = self.scale_start + \
                          fraction*(self.scale_stop - self.scale_start)
          self.viewupright()

	if n == nstart or self.panflag: self.center = compute_center(box)

        if bonds: self.bonds_augment(bonds)

        self.boxdraw = box
        self.atomdraw = atoms
        self.bonddraw = bonds
        self.tridraw = tris
        self.linedraw = lines

        self.ready = 1
        self.setview()
        self.cachelist = -self.cachelist
        self.w.tkRedraw()
        self.save(file)

        print n,
        sys.stdout.flush()
        n += 1

    print "\n%d images" % ncount

  # --------------------------------------------------------------------

  def display(self,index):
    self.boxdraw = self.boxframes[index]
    self.atomdraw = self.atomframes[index]
    self.bonddraw = self.bondframes[index]
    self.tridraw = self.triframes[index]
    self.linedraw = self.lineframes[index]

    self.ready = 1
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
    return (self.timeframes[index],len(self.atomdraw))

  # --------------------------------------------------------------------
  # draw the GL scene
  
  def redraw(self,o):
    # clear window to background color
    
    glClearColor(self.bgcol[0],self.bgcol[1],self.bgcol[2],0)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # not ready if no scene yet
    
    if not self.ready: return

    # set view from eye, distance, 3 lookat vectors (from,to,up)
    
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    if self.orthoflag:
      glOrtho(-0.25*self.eye,0.25*self.eye,-0.25*self.eye,0.25*self.eye,
              self.eye-2*self.distance,self.eye+2*self.distance)
    else:
      gluPerspective(30.0,1.0,0.01,10000.0)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(self.xfrom,self.yfrom,self.zfrom,self.xto,self.yto,self.zto,
              self.up[0],self.up[1],self.up[2])

    # draw scene from display list if caching allowed and list hasn't changed
    # else redraw and store as new display list if caching allowed
    
    if self.cache and self.cachelist > 0: glCallList(self.cachelist);
    else:
      if self.cache:
        if self.cachelist < 0: glDeleteLists(-self.cachelist,1)
        self.cachelist = glGenLists(1)
        glNewList(self.cachelist,GL_COMPILE_AND_EXECUTE)
        
      # draw box, clip-box, xyz axes, lines

      glDisable(GL_LIGHTING)

      if self.boxflag:
        self.draw_box(0)
        if self.clipflag: self.draw_box(1)
      if self.axisflag: self.draw_axes()

      ncolor = self.vizinfo.nlcolor
      for line in self.linedraw:
        itype = int(line[1])
        if itype > ncolor: raise StandardError,"line type too big"
        red,green,blue = self.vizinfo.lcolor[itype]
        glColor3f(red,green,blue)
        thick = self.vizinfo.lrad[itype]
	glLineWidth(thick)
        glBegin(GL_LINES)
        glVertex3f(line[2],line[3],line[4])
        glVertex3f(line[5],line[6],line[7])
        glEnd()

      glEnable(GL_LIGHTING)

      # draw non-clipped scene = atoms, bonds, triangles

# draw atoms as collection of points
# cannot put PointSize inside glBegin
#   so probably need to group atoms by type for best performance
#   or just allow one radius
# need to scale radius appropriately with box size
#   or could leave it at absolute value
# use POINT_SMOOTH to enable anti-aliasing and round points
# multiple timesteps via vcr::play() is still not fast
#  caching makes it fast for single frame, but multiple frames is slow
# need to enable clipping

#      if not self.clipflag:
#        glDisable(GL_LIGHTING)
#        glEnable(GL_POINT_SMOOTH)
#        glPointSize(self.vizinfo.arad[int(self.atomdraw[0][1])])
#        glBegin(GL_POINTS)
#        for atom in self.atomdraw:
#          red,green,blue = self.vizinfo.acolor[int(atom[1])]
#          glColor(red,green,blue)
#          glVertex3d(atom[2],atom[3],atom[4])
#        glEnd()
#        glEnable(GL_LIGHTING)

      if not self.clipflag:
        for atom in self.atomdraw:
          glTranslatef(atom[2],atom[3],atom[4]);
          glCallList(self.calllist[int(atom[1])]);
          glTranslatef(-atom[2],-atom[3],-atom[4]);

        if self.bonddraw:
          bound = 0.25 * self.distance
          ncolor = self.vizinfo.nbcolor
          for bond in self.bonddraw:
            if bond[10] > bound: continue
            itype = int(bond[1])
            if itype > ncolor: raise StandardError,"bond type too big"
            red,green,blue = self.vizinfo.bcolor[itype]
            rad = self.vizinfo.brad[itype]
            glPushMatrix()
            glTranslatef(bond[2],bond[3],bond[4])
            glRotatef(bond[11],bond[12],bond[13],0.0)
            glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,[red,green,blue,1.0]);
            glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,self.shiny);
            obj = gluNewQuadric()
            gluCylinder(obj,rad,rad,bond[10],self.nsides,self.nsides)
            glPopMatrix()

        if self.tridraw:
          fillflag = self.vizinfo.tfill[int(self.tridraw[0][1])]
    
          if fillflag != 1:
            if fillflag:
              glEnable(GL_POLYGON_OFFSET_FILL)
              glPolygonOffset(1.0,1.0)
            glBegin(GL_TRIANGLES)
            ncolor = self.vizinfo.ntcolor
            for tri in self.tridraw:
              itype = int(tri[1])
              if itype > ncolor: raise StandardError,"tri type too big"
              red,green,blue = self.vizinfo.tcolor[itype]
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,[red,green,blue,1.0]);
              glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,self.shiny);
              glNormal3f(tri[11],tri[12],tri[13])
              glVertex3f(tri[2],tri[3],tri[4])
              glVertex3f(tri[5],tri[6],tri[7])
              glVertex3f(tri[8],tri[9],tri[10])
            glEnd()
            if fillflag: glDisable(GL_POLYGON_OFFSET_FILL)

          if fillflag:
            glDisable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
            glLineWidth(self.bxthick)
            glColor3f(self.bxcol[0],self.bxcol[1],self.bxcol[2])
            glBegin(GL_TRIANGLES)
            for tri in self.tridraw:
              glVertex3f(tri[2],tri[3],tri[4])
              glVertex3f(tri[5],tri[6],tri[7])
              glVertex3f(tri[8],tri[9],tri[10])
            glEnd()
            glEnable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)

      # draw clipped scene = atoms, bonds, triangles

      else:
        box = self.boxdraw
        xlo = box[0] + self.clipxlo*(box[3] - box[0])
        xhi = box[0] + self.clipxhi*(box[3] - box[0])
        ylo = box[1] + self.clipylo*(box[4] - box[1])
        yhi = box[1] + self.clipyhi*(box[4] - box[1])
        zlo = box[2] + self.clipzlo*(box[5] - box[2])
        zhi = box[2] + self.clipzhi*(box[5] - box[2])

        for atom in self.atomdraw:
          x,y,z = atom[2],atom[3],atom[4]
          if x >= xlo and x <= xhi and y >= ylo and y <= yhi and \
                 z >= zlo and z <= zhi:
            glTranslatef(x,y,z);
            glCallList(self.calllist[int(atom[1])]);
            glTranslatef(-x,-y,-z);

        if self.bonddraw:
          bound = 0.25 * self.distance
          ncolor = self.vizinfo.nbcolor
          for bond in self.bonddraw:
            xmin = min2(bond[2],bond[5])
            xmax = max2(bond[2],bond[5])
            ymin = min2(bond[3],bond[6])
            ymax = max2(bond[3],bond[6])
            zmin = min2(bond[4],bond[7])
            zmax = max2(bond[4],bond[7])
            if xmin >= xlo and xmax <= xhi and \
                   ymin >= ylo and ymax <= yhi and zmin >= zlo and zmax <= zhi:
              if bond[10] > bound: continue
              itype = int(bond[1])
              if itype > ncolor: raise StandardError,"bond type too big"
              red,green,blue = self.vizinfo.bcolor[itype]
              rad = self.vizinfo.brad[itype]
              glPushMatrix()
              glTranslatef(bond[2],bond[3],bond[4])
              glRotatef(bond[11],bond[12],bond[13],0.0)
              glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,[red,green,blue,1.0]);
              glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,self.shiny);
              obj = gluNewQuadric()
              gluCylinder(obj,rad,rad,bond[10],self.nsides,self.nsides)
              glPopMatrix()

        if self.tridraw:  
          fillflag = self.vizinfo.tfill[int(self.tridraw[0][1])]

          if fillflag != 1:
            if fillflag:
              glEnable(GL_POLYGON_OFFSET_FILL)
              glPolygonOffset(1.0,1.0)
            glBegin(GL_TRIANGLES)
            ncolor = self.vizinfo.ntcolor
            for tri in self.tridraw:
              xmin = min3(tri[2],tri[5],tri[8])
              xmax = max3(tri[2],tri[5],tri[8])
              ymin = min3(tri[3],tri[6],tri[9])
              ymax = max3(tri[3],tri[6],tri[9])
              zmin = min3(tri[4],tri[7],tri[10])
              zmax = max3(tri[4],tri[7],tri[10])
              if xmin >= xlo and xmax <= xhi and \
                     ymin >= ylo and ymax <= yhi and \
                     zmin >= zlo and zmax <= zhi:
                itype = int(tri[1])
                if itype > ncolor: raise StandardError,"tri type too big"
                red,green,blue = self.vizinfo.tcolor[itype]
                glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,
                             [red,green,blue,1.0]);
                glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,self.shiny);
                glNormal3f(tri[11],tri[12],tri[13])
                glVertex3f(tri[2],tri[3],tri[4])
                glVertex3f(tri[5],tri[6],tri[7])
                glVertex3f(tri[8],tri[9],tri[10])
            glEnd()
            if fillflag: glDisable(GL_POLYGON_OFFSET_FILL)

          if fillflag:
            glDisable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
            glLineWidth(self.bxthick)
            glColor3f(self.bxcol[0],self.bxcol[1],self.bxcol[2])
            glBegin(GL_TRIANGLES)
            for tri in self.tridraw:
              xmin = min3(tri[2],tri[5],tri[8])
              xmax = max3(tri[2],tri[5],tri[8])
              ymin = min3(tri[3],tri[6],tri[9])
              ymax = max3(tri[3],tri[6],tri[9])
              zmin = min3(tri[4],tri[7],tri[10])
              zmax = max3(tri[4],tri[7],tri[10])
              if xmin >= xlo and xmax <= xhi and \
                     ymin >= ylo and ymax <= yhi and \
                     zmin >= zlo and zmax <= zhi:
                glVertex3f(tri[2],tri[3],tri[4])
                glVertex3f(tri[5],tri[6],tri[7])
                glVertex3f(tri[8],tri[9],tri[10])
            glEnd()
            glEnable(GL_LIGHTING)
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
      
      if self.cache: glEndList()

    glFlush()

  # --------------------------------------------------------------------
  # make new call list for each atom type
  # called when atom color/rad/quality is changed
  
  def make_atom_calllist(self):
    # extend calllist array if necessary
    
    if self.vizinfo.nacolor > self.nclist:
      for i in range(self.vizinfo.nacolor-self.nclist): self.calllist.append(0)
      self.nclist = self.vizinfo.nacolor

    # create new calllist for each atom type
    
    for itype in xrange(1,self.vizinfo.nacolor+1):
      if self.calllist[itype]: glDeleteLists(self.calllist[itype],1)
      ilist = glGenLists(1)
      self.calllist[itype] = ilist
      glNewList(ilist,GL_COMPILE)
      red,green,blue = self.vizinfo.acolor[itype]
      rad = self.vizinfo.arad[itype]
      glColor3f(red,green,blue);
      
#      glPointSize(10.0*rad)
#      glBegin(GL_POINTS)
#      glVertex3f(0.0,0.0,0.0)
#      glEnd()
      
      glMaterialfv(GL_FRONT,GL_EMISSION,[red,green,blue,1.0]);
      glMaterialf(GL_FRONT,GL_SHININESS,self.shiny);
      glutSolidSphere(rad,self.nslices,self.nstacks)
      glEndList()

  # --------------------------------------------------------------------
  # augment bond info returned by viz() with info needed for GL draw
  # info = length, theta, -dy, dx for bond orientation
  
  def bonds_augment(self,bonds):
    for bond in bonds:
      dx = bond[5] - bond[2]
      dy = bond[6] - bond[3]
      dz = bond[7] - bond[4]
      length = sqrt(dx*dx + dy*dy + dz*dz)
      dx /= length
      dy /= length
      dz /= length
      theta = acos(dz)*180.0/pi
      bond += [length,theta,-dy,dx]

  # --------------------------------------------------------------------

  def draw_box(self,flag):
    xlo,ylo,zlo,xhi,yhi,zhi = self.boxdraw

    if flag:
      tmp = xlo + self.clipxlo*(xhi - xlo)
      xhi = xlo + self.clipxhi*(xhi - xlo)
      xlo = tmp
      tmp = ylo + self.clipylo*(yhi - ylo)
      yhi = ylo + self.clipyhi*(yhi - ylo)
      ylo = tmp
      tmp = zlo + self.clipzlo*(zhi - zlo)
      zhi = zlo + self.clipzhi*(zhi - zlo)
      zlo = tmp

    glLineWidth(self.bxthick)
    glColor3f(self.bxcol[0],self.bxcol[1],self.bxcol[2])
    
    glBegin(GL_LINE_LOOP)
    glVertex3f(xlo,ylo,zlo)
    glVertex3f(xhi,ylo,zlo)
    glVertex3f(xhi,yhi,zlo)
    glVertex3f(xlo,yhi,zlo)
    glEnd()

    glBegin(GL_LINE_LOOP)
    glVertex3f(xlo,ylo,zhi)
    glVertex3f(xhi,ylo,zhi)
    glVertex3f(xhi,yhi,zhi)
    glVertex3f(xlo,yhi,zhi)
    glEnd()

    glBegin(GL_LINES)
    glVertex3f(xlo,ylo,zlo)
    glVertex3f(xlo,ylo,zhi)
    glVertex3f(xhi,ylo,zlo)
    glVertex3f(xhi,ylo,zhi)
    glVertex3f(xhi,yhi,zlo)
    glVertex3f(xhi,yhi,zhi)
    glVertex3f(xlo,yhi,zlo)
    glVertex3f(xlo,yhi,zhi)
    glEnd()

  # --------------------------------------------------------------------

  def draw_axes(self):
    xlo,ylo,zlo,xhi,yhi,zhi = self.boxdraw

    delta = xhi-xlo
    if yhi-ylo > delta: delta = yhi-ylo
    if zhi-zlo > delta: delta = zhi-zlo
    delta *= 0.1
      
    glLineWidth(self.bxthick)

    glBegin(GL_LINES)
    glColor3f(1,0,0)
    glVertex3f(xlo-delta,ylo-delta,zlo-delta)
    glVertex3f(xhi-delta,ylo-delta,zlo-delta)
    glColor3f(0,1,0)
    glVertex3f(xlo-delta,ylo-delta,zlo-delta)
    glVertex3f(xlo-delta,yhi-delta,zlo-delta)
    glColor3f(0,0,1)
    glVertex3f(xlo-delta,ylo-delta,zlo-delta)
    glVertex3f(xlo-delta,ylo-delta,zhi-delta)
    glEnd()

  # --------------------------------------------------------------------

  def save(self,file=None):
    self.w.update()      # force image on screen to be current before saving it
        
    pstring = glReadPixels(0,0,self.xpixels,self.ypixels,
                           GL_RGBA,GL_UNSIGNED_BYTE)
    snapshot = Image.fromstring("RGBA",(self.xpixels,self.ypixels),pstring)
    snapshot = snapshot.transpose(Image.FLIP_TOP_BOTTOM)

    if not file: file = self.file
    snapshot.save(file + ".png")

  # --------------------------------------------------------------------
  
  def adef(self):
    self.vizinfo.setcolors("atom",range(100),"loop")
    self.vizinfo.setradii("atom",range(100),0.45)
    self.make_atom_calllist()
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
    
  # --------------------------------------------------------------------

  def bdef(self):
    self.vizinfo.setcolors("bond",range(100),"loop")
    self.vizinfo.setradii("bond",range(100),0.25)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def tdef(self):
    self.vizinfo.setcolors("tri",range(100),"loop")
    self.vizinfo.setfills("tri",range(100),0)  
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def ldef(self):
    self.vizinfo.setcolors("line",range(100),"loop")  
    self.vizinfo.setradii("line",range(100),0.25)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def acol(self,atypes,colors):
    self.vizinfo.setcolors("atom",atypes,colors)
    self.make_atom_calllist()
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
    
  # --------------------------------------------------------------------

  def arad(self,atypes,radii):
    self.vizinfo.setradii("atom",atypes,radii)  
    self.make_atom_calllist()
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
    
  # --------------------------------------------------------------------

  def bcol(self,btypes,colors):
    self.vizinfo.setcolors("bond",btypes,colors)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
  
  # --------------------------------------------------------------------

  def brad(self,btypes,radii):
    self.vizinfo.setradii("bond",btypes,radii)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()
  
  # --------------------------------------------------------------------

  def tcol(self,ttypes,colors):
    self.vizinfo.setcolors("tri",ttypes,colors)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def tfill(self,ttypes,flags):
    self.vizinfo.setfills("tri",ttypes,flags)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def lcol(self,ltypes,colors):
    self.vizinfo.setcolors("line",ltypes,colors)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

  # --------------------------------------------------------------------

  def lrad(self,ltypes,radii):
    self.vizinfo.setradii("line",ltypes,radii)
    self.cachelist = -self.cachelist
    self.w.tkRedraw()

# --------------------------------------------------------------------
# derived class from Togl's Opengl
# overwrite redraw, translate, rotate, scale methods
# latter 3 are mouse-motion methods

class MyOpengl(Opengl):
  def __init__(self, master, cnf={}, **kw):
    args = (self,master,cnf)
    Opengl.__init__(*args,**kw)
    Opengl.autospin_allowed = 0
    
  # redraw Opengl scene
  # call parent redraw() method
  
  def tkRedraw(self,*dummy):
    if not self.initialised: return
    self.tk.call(self._w,'makecurrent')
    self.redraw(self)
    self.tk.call(self._w,'swapbuffers')

  # left button translate
  # access parent xshift/yshift and call parent trans() method
  
  def tkTranslate(self,event):
    dx = event.x - self.xmouse
    dy = event.y - self.ymouse
    x = self.parent.xshift + dx
    y = self.parent.yshift - dy
    self.parent.shift(x,y)
    self.tkRedraw()
    self.tkRecordMouse(event)

  # middle button trackball
  # call parent mouse_rotate() method

  def tkRotate(self,event):
    self.parent.mouse_rotate(event.x,event.y,self.xmouse,self.ymouse)
    self.tkRedraw()
    self.tkRecordMouse(event)

  # right button zoom
  # access parent scale and call parent zoom() method
  
  def tkScale(self,event):
    scale = 1 - 0.01 * (event.y - self.ymouse)
    if scale < 0.001: scale = 0.001
    elif scale > 1000: scale = 1000
    scale *= self.parent.scale
    self.parent.zoom(scale)
    self.tkRedraw()
    self.tkRecordMouse(event)

# --------------------------------------------------------------------
# draw a line segment

def segment(p1,p2):
  glVertex3f(p1[0],p1[1],p1[2])
  glVertex3f(p2[0],p2[1],p2[2])

# --------------------------------------------------------------------
# normalize a 3-vector to unit length

def vecnorm(v):
  length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return [v[0]/length,v[1]/length,v[2]/length]

# --------------------------------------------------------------------
# dot product of two 3-vectors

def vecdot(v1,v2):
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

# --------------------------------------------------------------------
# cross product of two 3-vectors

def veccross(v1,v2):
  v = [0,0,0]
  v[0] = v1[1]*v2[2] - v1[2]*v2[1]
  v[1] = v1[2]*v2[0] - v1[0]*v2[2]
  v[2] = v1[0]*v2[1] - v1[1]*v2[0]
  return v

# --------------------------------------------------------------------
# return characteristic distance of simulation domain = max dimension

def compute_distance(box):
  distance = box[3]-box[0]
  if box[4]-box[1] > distance: distance = box[4]-box[1]
  if box[5]-box[2] > distance: distance = box[5]-box[2]
  return distance

# --------------------------------------------------------------------
# return center of box as 3 vector

def compute_center(box):
  c = [0,0,0]
  c[0] = 0.5 * (box[0] + box[3])
  c[1] = 0.5 * (box[1] + box[4])
  c[2] = 0.5 * (box[2] + box[5])
  return c

# --------------------------------------------------------------------
# return min of 2 values

def min2(a,b):
  if b < a: a = b
  return a

# --------------------------------------------------------------------
# return max of 2 values

def max2(a,b):
  if b > a: a = b
  return a

# --------------------------------------------------------------------
# return min of 3 values

def min3(a,b,c):
  if b < a: a = b
  if c < a: a = c
  return a

# --------------------------------------------------------------------
# return max of 3 values

def max3(a,b,c):
  if b > a: a = b
  if c > a: a = c
  return a
