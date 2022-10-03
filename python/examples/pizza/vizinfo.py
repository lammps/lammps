# Pizza.py toolkit, https://lammps.github.io/pizza
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# vizinfo class, not a top-level Pizza.py tool

# History
#   8/05, Matt Jones (BYU): original version
#   9/05, Steve Plimpton: added 140-color table

# ToDo list

# Variables

# Imports and external programs

import types

# Class definition

class vizinfo:
  """
  Information holder for Pizza.py visualization tools
  
  acolor,bcolor,tcolor,lcolor = RGB values for each atom/bond/tri/line type
  arad = radius of each atom type
  brad,lrad = thickness of each bond/line type
  tfill = fill flag for each triangle type

  all of these arrays are indexed by object type which runs 1-Ntype
  nacolor,nbcolor,ntcolor,nlcolor,narad,nbrad,nlrad,ntfill
    are # of types each array holds
  actual length is nacolor+1 so that array can be indexed by 1-Ntype

  setcolors() = set atom/bond/tri/line colors
  setradii() = set atom/bond/line radii/thickness
  setfill() = set triangle fill factor
  extend() = grow an array
  """
  
  # --------------------------------------------------------------------

  def __init__(self):
     self.acolor = []
     self.arad = []
     self.bcolor = []
     self.brad = []
     self.tcolor = []
     self.tfill = []
     self.lcolor = []
     self.lrad = []
     self.nacolor = self.narad = 0
     self.nbcolor = self.nbrad = 0
     self.ntcolor = self.ntfill = 0
     self.nlcolor = self.nlrad = 0
     
  # --------------------------------------------------------------------
  # set color RGB for which = atoms, bonds, triangles
  
  def setcolors(self,which,ids,rgbs):

    # convert args into lists if single values
    # if arg = 0, convert to full-range list
    
    if type(ids) is types.IntType and ids == 0:
      if which == "atom": ids = range(self.nacolor)
      if which == "bond": ids = range(self.nbcolor)
      if which == "tri": ids = range(self.ntcolor)
      if which == "line": ids = range(self.nlcolor)
    if type(ids) is not types.ListType and type(ids) is not types.TupleType:
      ids = [ids]
    if type(rgbs) is not types.ListType and type(rgbs) is not types.TupleType:
      rgbs = [rgbs]

    # if list of types has a 0, increment each type value

    if 0 in ids:
      for i in range(len(ids)): ids[i] += 1

    # extend storage list if necessary
    # extend other arrays for same "which" so that gl::make_atom_calllist
    #   has valid arrays to work with

    if which == "atom":
      if max(ids) > self.nacolor:
        self.nacolor = self.extend(self.acolor,max(ids))
        self.nacolor = self.extend(self.arad,max(ids))
    if which == "bond":
      if max(ids) > self.nbcolor:
        self.nbcolor = self.extend(self.bcolor,max(ids))
        self.nbcolor = self.extend(self.brad,max(ids))
    if which == "tri":
      if max(ids) > self.ntcolor:
        self.ntcolor = self.extend(self.tcolor,max(ids))
        self.ntcolor = self.extend(self.tfill,max(ids))
    if which == "line":
      if max(ids) > self.nlcolor:
        self.nlcolor = self.extend(self.lcolor,max(ids))
        self.nlcolor = self.extend(self.lrad,max(ids))
    
    # set color for each type
    # if list lengths match, set directly, else interpolate
    # convert final color from 0-255 to 0.0-1.0
    
    ntypes = len(ids)
    nrgbs = len(rgbs)

    for i in range(ntypes):
      id = ids[i]

      if rgbs[0] == "loop":
        list = colors.keys()
        red,green,blue = colors[list[i % len(colors)]]
      elif ntypes == nrgbs:
        red,green,blue = colors[rgbs[i]]
      else:
        r = i/float(ntypes-1) * float(nrgbs-1)
        jlo = int(r)
        jhi = jlo + 1
        if jhi == nrgbs: jhi = nrgbs - 1
        clo = colors[rgbs[jlo]]
        chi = colors[rgbs[jhi]]
        delta = r - jlo
        red = clo[0] + delta*(chi[0]-clo[0])
        green = clo[1] + delta*(chi[1]-clo[1])
        blue = clo[2] + delta*(chi[2]-clo[2])

      color = [red/255.0,green/255.0,blue/255.0]

      if which == "atom": self.acolor[id] = color
      if which == "bond": self.bcolor[id] = color
      if which == "tri":  self.tcolor[id] = color
      if which == "line": self.lcolor[id] = color
  
  # --------------------------------------------------------------------
  # set radii for which = atoms, bonds, lines

  def setradii(self,which,ids,radii):

    # convert args into lists if single values
    # if arg = 0, convert to full-range list
    
    if type(ids) is types.IntType and ids == 0:
      if which == "atom": ids = range(self.narad)
      if which == "bond": ids = range(self.nbrad)
      if which == "line": ids = range(self.nlrad)
    if type(ids) is not types.ListType and type(ids) is not types.TupleType:
      ids = [ids]
    if type(radii) is not types.ListType and \
           type(radii) is not types.TupleType:
      radii = [radii]

    # if list of types has a 0, increment each type value

    if 0 in ids:
      for i in range(len(ids)): ids[i] += 1

    # extend storage list if necessary
    # extend other arrays for same "which" so that gl::make_atom_calllist
    #   has valid arrays to work with

    if which == "atom":
      if max(ids) > self.narad:
        self.narad = self.extend(self.arad,max(ids))
        self.narad = self.extend(self.acolor,max(ids))
    if which == "bond":
      if max(ids) > self.nbrad:
        self.nbrad = self.extend(self.brad,max(ids))
        self.nbrad = self.extend(self.bcolor,max(ids))
    if which == "line":
      if max(ids) > self.nlrad:
        self.nlrad = self.extend(self.lrad,max(ids))
        self.nlrad = self.extend(self.lcolor,max(ids))

    # set radius for each type
    # if list lengths match, set directly, else interpolate

    ntypes = len(ids)
    nradii = len(radii)

    for i in range(ntypes):
      id = ids[i]

      if ntypes == nradii: rad = radii[i]
      else:
        r = i/float(ntypes-1) * float(nradii-1)
        jlo = int(r)
        jhi = jlo + 1
        if jhi == nradii: jhi = nradii - 1
        rlo = radii[jlo]
        rhi = radii[jhi]
        delta = r - jlo
        rad = rlo + delta*(rhi-rlo)

      if which == "atom": self.arad[id] = rad
      if which == "bond": self.brad[id] = rad
      if which == "line": self.lrad[id] = rad
    
  # --------------------------------------------------------------------
  # set triangle fill style
  # 0 = fill only, 1 = line only, 2 = fill and line
  
  def setfills(self,which,ids,fills):

    # convert args into lists if single values
    # if arg = 0, convert to full-range list
    
    if type(ids) is types.IntType and ids == 0:
      ids = range(self.ntfill)
    if type(ids) is not types.ListType and type(ids) is not types.TupleType:
      ids = [ids]
    if type(fills) is not types.ListType and \
           type(fills) is not types.TupleType:
      fills = [fills]

    # if list of types has a 0, increment each type value

    if 0 in ids:
      for i in range(len(ids)): ids[i] += 1

    # extend storage list if necessary
    # extend other arrays for same "which" so that gl::make_atom_calllist
    #   has valid arrays to work with

    if max(ids) > self.ntfill:
      self.ntfill = self.extend(self.tfill,max(ids))
      self.ntfill = self.extend(self.tcolor,max(ids))

    # set fill flag for each type
    # if list lengths match, set directly, else set types to 1st fill value

    if len(fills) == len(ids):
      for i in range(len(ids)): self.tfill[ids[i]] = int(fills[i])
    else:
      for id in ids: self.tfill[id] = int(fills[0])
    
  # --------------------------------------------------------------------

  def extend(self,array,n):
    for i in range(n-len(array)+1): array.append(0)
    return n

# --------------------------------------------------------------------
# dictionary of 140 color names and associated RGB values

colors = {}

colors["aliceblue"] = [240, 248, 255]
colors["antiquewhite"] = [250, 235, 215]
colors["aqua"] = [0, 255, 255]
colors["aquamarine"] = [127, 255, 212]
colors["azure"] = [240, 255, 255]
colors["beige"] = [245, 245, 220]
colors["bisque"] = [255, 228, 196]
colors["black"] = [0, 0, 0]
colors["blanchedalmond"] = [255, 255, 205]
colors["blue"] = [0, 0, 255]
colors["blueviolet"] = [138, 43, 226]
colors["brown"] = [165, 42, 42]
colors["burlywood"] = [222, 184, 135]
colors["cadetblue"] = [95, 158, 160]
colors["chartreuse"] = [127, 255, 0]
colors["chocolate"] = [210, 105, 30]
colors["coral"] = [255, 127, 80]
colors["cornflowerblue"] = [100, 149, 237]
colors["cornsilk"] = [255, 248, 220]
colors["crimson"] = [220, 20, 60]
colors["cyan"] = [0, 255, 255]
colors["darkblue"] = [0, 0, 139]
colors["darkcyan"] = [0, 139, 139]
colors["darkgoldenrod"] = [184, 134, 11]
colors["darkgray"] = [169, 169, 169]
colors["darkgreen"] = [0, 100, 0]
colors["darkkhaki"] = [189, 183, 107]
colors["darkmagenta"] = [139, 0, 139]
colors["darkolivegreen"] = [85, 107, 47]
colors["darkorange"] = [255, 140, 0]
colors["darkorchid"] = [153, 50, 204]
colors["darkred"] = [139, 0, 0]
colors["darksalmon"] = [233, 150, 122]
colors["darkseagreen"] = [143, 188, 143]
colors["darkslateblue"] = [72, 61, 139]
colors["darkslategray"] = [47, 79, 79]
colors["darkturquoise"] = [0, 206, 209]
colors["darkviolet"] = [148, 0, 211]
colors["deeppink"] = [255, 20, 147]
colors["deepskyblue"] = [0, 191, 255]
colors["dimgray"] = [105, 105, 105]
colors["dodgerblue"] = [30, 144, 255]
colors["firebrick"] = [178, 34, 34]
colors["floralwhite"] = [255, 250, 240]
colors["forestgreen"] = [34, 139, 34]
colors["fuchsia"] = [255, 0, 255]
colors["gainsboro"] = [220, 220, 220]
colors["ghostwhite"] = [248, 248, 255]
colors["gold"] = [255, 215, 0]
colors["goldenrod"] = [218, 165, 32]
colors["gray"] = [128, 128, 128]
colors["green"] = [0, 128, 0]
colors["greenyellow"] = [173, 255, 47]
colors["honeydew"] = [240, 255, 240]
colors["hotpink"] = [255, 105, 180]
colors["indianred"] = [205, 92, 92]
colors["indigo"] = [75, 0, 130]
colors["ivory"] = [255, 240, 240]
colors["khaki"] = [240, 230, 140]
colors["lavender"] = [230, 230, 250]
colors["lavenderblush"] = [255, 240, 245]
colors["lawngreen"] = [124, 252, 0]
colors["lemonchiffon"] = [255, 250, 205]
colors["lightblue"] = [173, 216, 230]
colors["lightcoral"] = [240, 128, 128]
colors["lightcyan"] = [224, 255, 255]
colors["lightgoldenrodyellow"] = [250, 250, 210]
colors["lightgreen"] = [144, 238, 144]
colors["lightgrey"] = [211, 211, 211]
colors["lightpink"] = [255, 182, 193]
colors["lightsalmon"] = [255, 160, 122]
colors["lightseagreen"] = [32, 178, 170]
colors["lightskyblue"] = [135, 206, 250]
colors["lightslategray"] = [119, 136, 153]
colors["lightsteelblue"] = [176, 196, 222]
colors["lightyellow"] = [255, 255, 224]
colors["lime"] = [0, 255, 0]
colors["limegreen"] = [50, 205, 50]
colors["linen"] = [250, 240, 230]
colors["magenta"] = [255, 0, 255]
colors["maroon"] = [128, 0, 0]
colors["mediumaquamarine"] = [102, 205, 170]
colors["mediumblue"] = [0, 0, 205]
colors["mediumorchid"] = [186, 85, 211]
colors["mediumpurple"] = [147, 112, 219]
colors["mediumseagreen"] = [60, 179, 113]
colors["mediumslateblue"] = [123, 104, 238]
colors["mediumspringgreen"] = [0, 250, 154]
colors["mediumturquoise"] = [72, 209, 204]
colors["mediumvioletred"] = [199, 21, 133]
colors["midnightblue"] = [25, 25, 112]
colors["mintcream"] = [245, 255, 250]
colors["mistyrose"] = [255, 228, 225]
colors["moccasin"] = [255, 228, 181]
colors["navajowhite"] = [255, 222, 173]
colors["navy"] = [0, 0, 128]
colors["oldlace"] = [253, 245, 230]
colors["olive"] = [128, 128, 0]
colors["olivedrab"] = [107, 142, 35]
colors["orange"] = [255, 165, 0]
colors["orangered"] = [255, 69, 0]
colors["orchid"] = [218, 112, 214]
colors["palegoldenrod"] = [238, 232, 170]
colors["palegreen"] = [152, 251, 152]
colors["paleturquoise"] = [175, 238, 238]
colors["palevioletred"] = [219, 112, 147]
colors["papayawhip"] = [255, 239, 213]
colors["peachpuff"] = [255, 239, 213]
colors["peru"] = [205, 133, 63]
colors["pink"] = [255, 192, 203]
colors["plum"] = [221, 160, 221]
colors["powderblue"] = [176, 224, 230]
colors["purple"] = [128, 0, 128]
colors["red"] = [255, 0, 0]
colors["rosybrown"] = [188, 143, 143]
colors["royalblue"] = [65, 105, 225]
colors["saddlebrown"] = [139, 69, 19]
colors["salmon"] = [250, 128, 114]
colors["sandybrown"] = [244, 164, 96]
colors["seagreen"] = [46, 139, 87]
colors["seashell"] = [255, 245, 238]
colors["sienna"] = [160, 82, 45]
colors["silver"] = [192, 192, 192]
colors["skyblue"] = [135, 206, 235]
colors["slateblue"] = [106, 90, 205]
colors["slategray"] = [112, 128, 144]
colors["snow"] = [255, 250, 250]
colors["springgreen"] = [0, 255, 127]
colors["steelblue"] = [70, 130, 180]
colors["tan"] = [210, 180, 140]
colors["teal"] = [0, 128, 128]
colors["thistle"] = [216, 191, 216]
colors["tomato"] = [253, 99, 71]
colors["turquoise"] = [64, 224, 208]
colors["violet"] = [238, 130, 238]
colors["wheat"] = [245, 222, 179]
colors["white"] = [255, 255, 255]
colors["whitesmoke"] = [245, 245, 245]
colors["yellow"] = [255, 255, 0]
colors["yellowgreen"] = [154, 205, 50]
