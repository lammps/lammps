# state file generated using paraview version 5.8.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

import glob;
import os;
import sys;
import pdb;

# @base_dir
#base_dir = './output/simulation_polymer4/fene_test003/'
#base_dir = './output/simulation_polymer4/harmonic_test008/'
base_dir = './'

print("base_dir = " + str(base_dir));

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1365, 747]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView1.OrientationAxesOutlineColor = [0.49019607843137253, 0.49019607843137253, 0.49019607843137253]
renderView1.CenterOfRotation = [0.0, 1e-20, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.0, -1355.156412055249, 0.0]
renderView1.CameraFocalPoint = [0.0, 1e-20, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 350.74028853269766
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitle = 'x'
renderView1.AxesGrid.YTitle = 'y'
renderView1.AxesGrid.ZTitle = 'z'
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.GridColor = [0.49019607843137253, 0.49019607843137253, 0.49019607843137253]
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.UseCustomBounds = 1
renderView1.AxesGrid.CustomBounds = [-0.15, 1.15, -0.15, 1.15, -0.6, 0.6]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Unstructured Grid Reader'
path_cwd = os.getcwd();
#os.chdir(base_dir + '/vtk');
#file_list = sorted(glob.glob('Particles_????????.vtp'));
#pdb.set_trace();
file_list = sorted(glob.glob(base_dir + '/vtk/Particles_????????.vtp'));
#os.chdir(path_cwd);
print("Particles:");
print("file_list[0] = " + str(file_list[0]));
print("file_list[-1] = " + str(file_list[-1]));
print("len(file_list) = " + str(len(file_list)));

polymerBeads = XMLPolyDataReader(FileName=file_list);
RenameSource("polymer",polymerBeads);
polymerBeads.PointArrayStatus = ['id', 'type', 'vx', 'fx']


file_list = sorted(glob.glob(base_dir + '/vtk/Particles_????????_boundingBox.vtu'));
#os.chdir(path_cwd);
print("Bounding Box:");
print("file_list[0] = " + str(file_list[0]));
print("file_list[-1] = " + str(file_list[-1]));
print("len(file_list) = " + str(len(file_list)));
boundingBox1 = XMLUnstructuredGridReader(FileName=file_list);
RenameSource("bounding_box1",boundingBox1);

# create a new 'Glyph'
glyph1 = Glyph(Input=polymerBeads,GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 33.459614562988286
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'
RenameSource("all_beads",glyph1);

# create a new 'Programmable Source'
programmableSource1 = ProgrammableSource()
RenameSource("programmableSource1",programmableSource1);
programmableSource1.Script = """#This script generates a helix curve.
#This is intended as the script of a 'Programmable Source'
import math

numPts = 80 # Points along Helix
length = 8.0 # Length of Helix
rounds = 3.0 # Number of times around

#Get a vtk.PolyData object for the output
pdo = self.GetPolyDataOutput()

#This will store the points for the Helix
newPts = vtk.vtkPoints()
for i in range(0, numPts):
   #Generate the Points along the Helix
   x = i*length/numPts
   y = math.sin(i*rounds*2*math.pi/numPts)
   z = math.cos(i*rounds*2*math.pi/numPts)
   #Insert the Points into the vtkPoints object
   #The first parameter indicates the reference.
   #value for the point. Here we add them sequentially.
   #Note that the first point is at index 0 (not 1).
   newPts.InsertPoint(i, x,y,z)

#Add the points to the vtkPolyData object
#Right now the points are not associated with a line - 
#it is just a set of unconnected points. We need to
#create a 'cell' object that ties points together
#to make a curve (in this case). This is done below.
#A 'cell' is just an object that tells how points are
#connected to make a 1D, 2D, or 3D object.
pdo.SetPoints(newPts)

#Make a vtkPolyLine which holds the info necessary
#to create a curve composed of line segments. This
#really just hold constructor data that will be passed
#to vtkPolyData to add a new line.
aPolyLine = vtk.vtkPolyLine()

#Indicate the number of points along the line
aPolyLine.GetPointIds().SetNumberOfIds(numPts)
for i in range(0,numPts):
   #Add the points to the line. The first value indicates
   #the order of the point on the line. The second value
   #is a reference to a point in a vtkPoints object. Depends
   #on the order that Points were added to vtkPoints object.
   #Note that this will not be associated with actual points
   #until it is added to a vtkPolyData object which holds a
   #vtkPoints object.
   aPolyLine.GetPointIds().SetId(i, i)

#Allocate the number of 'cells' that will be added. We are just
#adding one vtkPolyLine 'cell' to the vtkPolyData object.
pdo.Allocate(1, 1)

#Add the poly line 'cell' to the vtkPolyData object.
pdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())

#The Helix is ready to plot! Click 'Apply'."""
programmableSource1.ScriptRequestInformation = ''
programmableSource1.PythonPath = ''


# create a new 'Programmable Filter'
bondsFilter1 = ProgrammableFilter(Input=polymerBeads)
RenameSource("bonds_filter1",bondsFilter1);
bondsFilter1.OutputDataSetType = 'vtkMolecule'
bondsFilter1.Script = """import numpy as np;

pdi = self.GetPolyDataInput();
pdo = self.GetMoleculeOutput()

#num_points = pdi.GetNumberOfPoints()
num_points = 100; num_dim = 3;
bond_list = []; X_list = np.zeros((num_points,num_dim));
for i in range(0, num_points):
  coord = pdi.GetPoint(i);
  #x, y, z = coord[:3]
  #X_list[i,:] = np.array([x,y,z]);
  X_list[i,:] = coord[:3];
  pdo.AppendAtom(1,X_list[i,0],X_list[i,1],X_list[i,2]);

flag = True;
if flag:
  ell = 170; # threshold distance
  ell_sq = ell*ell;
  for i in range(0, num_points - 1):
    X1 = X_list[i,:]; X2 = X_list[i + 1,:];
    norm_sq = np.sum(np.power(X1 - X2,2));    
    if norm_sq < ell_sq:
      pdo.AppendBond(i,i+1,1);

"""
bondsFilter1.RequestInformationScript = ''
bondsFilter1.RequestUpdateExtentScript = ''
bondsFilter1.PythonPath = ''

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from polymerBeads
polymerBeadsDisplay = Show(polymerBeads, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
polymerBeadsDisplay.Representation = 'Wireframe'
polymerBeadsDisplay.AmbientColor = [0.0, 0.3333333333333333, 0.4980392156862745]
polymerBeadsDisplay.ColorArrayName = [None, '']
polymerBeadsDisplay.DiffuseColor = [0.0, 0.3333333333333333, 0.4980392156862745]
polymerBeadsDisplay.Opacity = 0.26
polymerBeadsDisplay.LineWidth = 5.0
polymerBeadsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
polymerBeadsDisplay.SelectOrientationVectors = 'None'
polymerBeadsDisplay.ScaleFactor = 40.400000000000006
polymerBeadsDisplay.SelectScaleArray = 'None'
polymerBeadsDisplay.GlyphType = 'Arrow'
polymerBeadsDisplay.GlyphTableIndexArray = 'None'
polymerBeadsDisplay.GaussianRadius = 2.02
polymerBeadsDisplay.SetScaleArray = [None, '']
polymerBeadsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
polymerBeadsDisplay.OpacityArray = [None, '']
polymerBeadsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
polymerBeadsDisplay.DataAxesGrid = 'GridAxesRepresentation'
polymerBeadsDisplay.PolarAxes = 'PolarAxesRepresentation'
polymerBeadsDisplay.ScalarOpacityUnitDistance = 699.7485262578264

# show data from polymerBeads
polymerBeadsDisplay = Show(polymerBeads, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
polymerBeadsDisplay.Representation = 'Surface'
polymerBeadsDisplay.ColorArrayName = [None, '']
polymerBeadsDisplay.OSPRayScaleArray = 'fx'
polymerBeadsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
polymerBeadsDisplay.SelectOrientationVectors = 'None'
polymerBeadsDisplay.ScaleFactor = 10.023999786376955
polymerBeadsDisplay.SelectScaleArray = 'None'
polymerBeadsDisplay.GlyphType = 'Arrow'
polymerBeadsDisplay.GlyphTableIndexArray = 'None'
polymerBeadsDisplay.GaussianRadius = 0.5011999893188477
polymerBeadsDisplay.SetScaleArray = ['POINTS', 'fx']
polymerBeadsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
polymerBeadsDisplay.OpacityArray = ['POINTS', 'fx']
polymerBeadsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
polymerBeadsDisplay.DataAxesGrid = 'GridAxesRepresentation'
polymerBeadsDisplay.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
polymerBeadsDisplay.ScaleTransferFunction.Points = [-388037.97892, 0.0, 0.5, 0.0, 388037.97892, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
polymerBeadsDisplay.OpacityTransferFunction.Points = [-388037.97892, 0.0, 0.5, 0.0, 388037.97892, 1.0, 0.5, 0.0]
polymerBeadsDisplay.Visibility = 0;

# create a new 'Calculator'
calculator2 = Calculator(Input=polymerBeads)
calculator2.ResultArrayName = 'indicator_2'
calculator2.Function = 'exp(100*(type - 2))'
RenameSource("indicator_2",calculator2);

# create a new 'Glyph'
glyph2 = Glyph(Input=calculator2,
    GlyphType='Sphere')
glyph2.OrientationArray = ['POINTS', 'No orientation array']
glyph2.ScaleArray = ['POINTS', 'indicator_2']
glyph2.ScaleFactor = 10.10439338684082
glyph2.GlyphTransform = 'Transform2'
glyph2.GlyphMode = 'All Points'
RenameSource("tracer_beads",glyph2);

# create a new 'Calculator'
calculator1 = Calculator(Input=polymerBeads)
calculator1.ResultArrayName = 'indicator_1'
calculator1.Function = 'exp(-100*(type - 1))'
RenameSource("indicator_1",calculator1);

# create a new 'Glyph'
glyph1 = Glyph(Input=calculator1,
    GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'indicator_1']
glyph1.ScaleFactor = 18.175
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'
RenameSource("polymer_beads",glyph1);

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from glyph1
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [1.0, 0.6666666666666666, 1.0]
glyph1Display.ColorArrayName = [None, '']
glyph1Display.DiffuseColor = [1.0, 0.6666666666666666, 1.0]
glyph1Display.OSPRayScaleArray = 'fx'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 39.98500061035156
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 1.9992500305175782
glyph1Display.SetScaleArray = ['POINTS', 'fx']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'fx']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-100697.49604, 0.0, 0.5, 0.0, 100697.49604, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-100697.49604, 0.0, 0.5, 0.0, 100697.49604, 1.0, 0.5, 0.0]

glyph1Display.Visibility = 0;

# show data from glyph2
glyph2Display = Show(glyph2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph2Display.Representation = 'Surface'
glyph2Display.AmbientColor = [0.0, 0.6666666666666666, 1.0]
glyph2Display.ColorArrayName = ['POINTS', '']
glyph2Display.DiffuseColor = [0.0, 0.6666666666666666, 1.0]
glyph2Display.Opacity = 0.2
glyph2Display.OSPRayScaleArray = 'indicator_2'
glyph2Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph2Display.SelectOrientationVectors = 'None'
glyph2Display.ScaleFactor = 40.41757354736328
glyph2Display.SelectScaleArray = 'indicator_2'
glyph2Display.GlyphType = 'Arrow'
glyph2Display.GlyphTableIndexArray = 'indicator_2'
glyph2Display.GaussianRadius = 2.020878677368164
glyph2Display.SetScaleArray = ['POINTS', 'indicator_2']
glyph2Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph2Display.OpacityArray = ['POINTS', 'indicator_2']
glyph2Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph2Display.DataAxesGrid = 'GridAxesRepresentation'
glyph2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [1.3838965267367376e-87, 0.0, 0.5, 0.0, 3.720075976020836e-44, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [1.3838965267367376e-87, 0.0, 0.5, 0.0, 3.720075976020836e-44, 1.0, 0.5, 0.0]


# show data from bondsFilter1
bondsFilter1Display = Show(bondsFilter1, renderView1, 'PVMoleculeRepresentation')

# trace defaults for the display properties.
bondsFilter1Display.Representation = 'Molecule'
bondsFilter1Display.AtomicRadiusFactor = 5.0
bondsFilter1Display.AtomicRadiusArrayName = 'Atomic Numbers'
bondsFilter1Display.BondRadius = 2.5
bondsFilter1Display.BondColorMode = 0
bondsFilter1Display.BondColor = [0.0, 0.6666666666666666, 1.0]


# setup the color legend parameters for each legend in this view

# get color transfer function/color map for 'AtomicNumbers'
atomicNumbersLUT = GetColorTransferFunction('AtomicNumbers')
atomicNumbersLUT.InterpretValuesAsCategories = 1
atomicNumbersLUT.AnnotationsInitialized = 1
atomicNumbersLUT.ShowCategoricalColorsinDataRangeOnly = 1
atomicNumbersLUT.Annotations = ['0', 'Xx', '1', 'H', '2', 'He', '3', 'Li', '4', 'Be', '5', 'B', '6', 'C', '7', 'N', '8', 'O', '9', 'F', '10', 'Ne', '11', 'Na', '12', 'Mg', '13', 'Al', '14', 'Si', '15', 'P', '16', 'S', '17', 'Cl', '18', 'Ar', '19', 'K', '20', 'Ca', '21', 'Sc', '22', 'Ti', '23', 'V', '24', 'Cr', '25', 'Mn', '26', 'Fe', '27', 'Co', '28', 'Ni', '29', 'Cu', '30', 'Zn', '31', 'Ga', '32', 'Ge', '33', 'As', '34', 'Se', '35', 'Br', '36', 'Kr', '37', 'Rb', '38', 'Sr', '39', 'Y', '40', 'Zr', '41', 'Nb', '42', 'Mo', '43', 'Tc', '44', 'Ru', '45', 'Rh', '46', 'Pd', '47', 'Ag', '48', 'Cd', '49', 'In', '50', 'Sn', '51', 'Sb', '52', 'Te', '53', 'I', '54', 'Xe', '55', 'Cs', '56', 'Ba', '57', 'La', '58', 'Ce', '59', 'Pr', '60', 'Nd', '61', 'Pm', '62', 'Sm', '63', 'Eu', '64', 'Gd', '65', 'Tb', '66', 'Dy', '67', 'Ho', '68', 'Er', '69', 'Tm', '70', 'Yb', '71', 'Lu', '72', 'Hf', '73', 'Ta', '74', 'W', '75', 'Re', '76', 'Os', '77', 'Ir', '78', 'Pt', '79', 'Au', '80', 'Hg', '81', 'Tl', '82', 'Pb', '83', 'Bi', '84', 'Po', '85', 'At', '86', 'Rn', '87', 'Fr', '88', 'Ra', '89', 'Ac', '90', 'Th', '91', 'Pa', '92', 'U', '93', 'Np', '94', 'Pu', '95', 'Am', '96', 'Cm', '97', 'Bk', '98', 'Cf', '99', 'Es', '100', 'Fm', '101', 'Md', '102', 'No', '103', 'Lr', '104', 'Rf', '105', 'Db', '106', 'Sg', '107', 'Bh', '108', 'Hs', '109', 'Mt', '110', 'Ds', '111', 'Rg', '112', 'Cn', '113', 'Uut', '114', 'Uuq', '115', 'Uup', '116', 'Uuh', '117', 'Uus', '118', 'Uuo']
atomicNumbersLUT.ActiveAnnotatedValues = ['1']
atomicNumbersLUT.IndexedColors = [0.06999313344014649, 0.5000076295109483, 0.7000076295109483, 1.0, 0.6666666666666666, 1.0, 0.8500038147554742, 1.0, 1.0, 0.8, 0.5000076295109483, 1.0, 0.7600061036087586, 1.0, 0.0, 1.0, 0.7100022888532845, 0.7100022888532845, 0.5000076295109483, 0.5000076295109483, 0.5000076295109483, 0.05000381475547418, 0.05000381475547418, 1.0, 1.0, 0.05000381475547418, 0.05000381475547418, 0.7000076295109483, 1.0, 1.0, 0.7000076295109483, 0.8899977111467154, 0.9600061036087587, 0.6699931334401464, 0.3600061036087587, 0.9499961852445258, 0.5400015259021896, 1.0, 0.0, 0.7499961852445258, 0.6500038147554742, 0.6500038147554742, 0.5000076295109483, 0.6, 0.6, 1.0, 0.5000076295109483, 0.0, 1.0, 1.0, 0.19000534065766383, 0.11999694819562066, 0.9400015259021897, 0.11999694819562066, 0.5000076295109483, 0.820004577706569, 0.8899977111467154, 0.5600061036087587, 0.2500038147554742, 0.8299992370489052, 0.23999389639124133, 1.0, 0.0, 0.9000076295109484, 0.9000076295109484, 0.9000076295109484, 0.7499961852445258, 0.7600061036087586, 0.779995422293431, 0.6500038147554742, 0.6500038147554742, 0.6699931334401464, 0.5400015259021896, 0.6, 0.779995422293431, 0.6099946593423362, 0.48000305180437935, 0.779995422293431, 0.5000076295109483, 0.48000305180437935, 0.779995422293431, 0.4399938963912413, 0.48000305180437935, 0.779995422293431, 0.3600061036087587, 0.48000305180437935, 0.7600061036087586, 1.0, 0.48000305180437935, 0.379995422293431, 0.4899977111467155, 0.5000076295109483, 0.6899977111467155, 0.7600061036087586, 0.5600061036087587, 0.5600061036087587, 0.4, 0.5600061036087587, 0.5600061036087587, 0.7400015259021897, 0.5000076295109483, 0.8899977111467154, 1.0, 0.6299992370489051, 0.0, 0.6500038147554742, 0.16000610360875867, 0.16000610360875867, 0.3600061036087587, 0.7199969481956207, 0.820004577706569, 0.4399938963912413, 0.179995422293431, 0.6899977111467155, 0.0, 1.0, 0.0, 0.579995422293431, 1.0, 1.0, 0.579995422293431, 0.8800030518043793, 0.8800030518043793, 0.4500038147554742, 0.7600061036087586, 0.7900053406576638, 0.3300068665598535, 0.7100022888532845, 0.7100022888532845, 0.22999923704890515, 0.620004577706569, 0.620004577706569, 0.14000152590218967, 0.5600061036087587, 0.5600061036087587, 0.03999389639124132, 0.4899977111467155, 0.5499961852445259, 0.0, 0.40999465934233614, 0.5199969481956207, 0.8800030518043793, 0.8800030518043793, 1.0, 1.0, 0.8500038147554742, 0.5600061036087587, 0.6500038147554742, 0.4599984740978103, 0.4500038147554742, 0.4, 0.5000076295109483, 0.5000076295109483, 0.620004577706569, 0.39000534065766385, 0.7100022888532845, 0.8299992370489052, 0.48000305180437935, 0.0, 0.579995422293431, 0.0, 0.579995422293431, 0.25999847409781035, 0.620004577706569, 0.6899977111467155, 0.3400015259021897, 0.0899977111467155, 0.5600061036087587, 0.0, 0.7900053406576638, 0.0, 0.4399938963912413, 0.8299992370489052, 1.0, 1.0, 1.0, 0.779995422293431, 0.8500038147554742, 1.0, 0.779995422293431, 0.779995422293431, 1.0, 0.779995422293431, 0.6399938963912413, 1.0, 0.779995422293431, 0.5600061036087587, 1.0, 0.779995422293431, 0.379995422293431, 1.0, 0.779995422293431, 0.26999313344014647, 1.0, 0.779995422293431, 0.19000534065766383, 1.0, 0.779995422293431, 0.11999694819562066, 1.0, 0.779995422293431, 0.0, 1.0, 0.6099946593423362, 0.0, 0.9000076295109484, 0.4599984740978103, 0.0, 0.8299992370489052, 0.31999694819562063, 0.0, 0.7499961852445258, 0.220004577706569, 0.0, 0.6699931334401464, 0.14000152590218967, 0.30000762951094834, 0.7600061036087586, 1.0, 0.30000762951094834, 0.6500038147554742, 1.0, 0.13000686655985352, 0.579995422293431, 0.8399938963912413, 0.14999618524452582, 0.4899977111467155, 0.6699931334401464, 0.14999618524452582, 0.4, 0.5900053406576639, 0.0899977111467155, 0.3300068665598535, 0.5300068665598535, 0.9600061036087587, 0.9300068665598535, 0.820004577706569, 0.8, 0.820004577706569, 0.11999694819562066, 0.7100022888532845, 0.7100022888532845, 0.7600061036087586, 0.6500038147554742, 0.3300068665598535, 0.30000762951094834, 0.3400015259021897, 0.3499961852445258, 0.379995422293431, 0.620004577706569, 0.3100022888532845, 0.7100022888532845, 0.6699931334401464, 0.3600061036087587, 0.0, 0.4599984740978103, 0.3100022888532845, 0.26999313344014647, 0.25999847409781035, 0.5100022888532845, 0.5900053406576639, 0.25999847409781035, 0.0, 0.4, 0.0, 0.4899977111467155, 0.0, 0.4399938963912413, 0.6699931334401464, 0.979995422293431, 0.0, 0.7300068665598535, 1.0, 0.0, 0.6299992370489051, 1.0, 0.0, 0.5600061036087587, 1.0, 0.0, 0.5000076295109483, 1.0, 0.0, 0.420004577706569, 1.0, 0.3300068665598535, 0.3600061036087587, 0.9499961852445258, 0.4699931334401465, 0.3600061036087587, 0.8899977111467154, 0.5400015259021896, 0.3100022888532845, 0.8899977111467154, 0.6299992370489051, 0.20999465934233616, 0.8299992370489052, 0.7000076295109483, 0.11999694819562066, 0.8299992370489052, 0.7000076295109483, 0.11999694819562066, 0.7300068665598535, 0.7000076295109483, 0.05000381475547418, 0.6500038147554742, 0.7400015259021897, 0.05000381475547418, 0.5300068665598535, 0.779995422293431, 0.0, 0.4, 0.8, 0.0, 0.3499961852445258, 0.820004577706569, 0.0, 0.3100022888532845, 0.8500038147554742, 0.0, 0.26999313344014647, 0.8800030518043793, 0.0, 0.220004577706569, 0.9000076295109484, 0.0, 0.179995422293431, 0.9100022888532845, 0.0, 0.14999618524452582, 0.9199969481956206, 0.0, 0.14000152590218967, 0.9300068665598535, 0.0, 0.13000686655985352, 0.9400015259021897, 0.0, 0.11999694819562066, 0.9499961852445258, 0.0, 0.1100022888532845, 0.9600061036087587, 0.0, 0.10000762951094835, 0.9700007629510948, 0.0, 0.0899977111467155, 0.979995422293431, 0.0, 0.08000305180437933, 0.9900053406576639, 0.0, 0.06999313344014649, 1.0, 0.0, 0.05999847409781033]
atomicNumbersLUT.IndexedOpacities = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

# get color legend/bar for atomicNumbersLUT in view renderView1
atomicNumbersLUTColorBar = GetScalarBar(atomicNumbersLUT, renderView1)
atomicNumbersLUTColorBar.Title = 'Atomic Numbers'
atomicNumbersLUTColorBar.ComponentTitle = ''

# set color bar visibility
atomicNumbersLUTColorBar.Visibility = 1

# show data from boundingBox1
boundingBox1Display = Show(boundingBox1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
boundingBox1Display.Representation = 'Wireframe'
boundingBox1Display.AmbientColor = [0.0, 0.0, 0.0]
boundingBox1Display.ColorArrayName = [None, '']
boundingBox1Display.DiffuseColor = [0.0, 0.0, 0.0]
boundingBox1Display.Opacity = 0.27
boundingBox1Display.LineWidth = 3.0
boundingBox1Display.OSPRayScaleFunction = 'PiecewiseFunction'
boundingBox1Display.SelectOrientationVectors = 'None'
boundingBox1Display.ScaleFactor = 40.5
boundingBox1Display.SelectScaleArray = 'None'
boundingBox1Display.GlyphType = 'Arrow'
boundingBox1Display.GlyphTableIndexArray = 'None'
boundingBox1Display.GaussianRadius = 2.025
boundingBox1Display.SetScaleArray = [None, '']
boundingBox1Display.ScaleTransferFunction = 'PiecewiseFunction'
boundingBox1Display.OpacityArray = [None, '']
boundingBox1Display.OpacityTransferFunction = 'PiecewiseFunction'
boundingBox1Display.DataAxesGrid = 'GridAxesRepresentation'
boundingBox1Display.PolarAxes = 'PolarAxesRepresentation'
boundingBox1Display.ScalarOpacityUnitDistance = 701.4805770653953


# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(programmableSource1)
# ----------------------------------------------------------------
