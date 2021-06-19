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

import glob,os,sys; #pdb;

# @base_dir
#base_dir = './output/simulation_particle4/fene_test003/'
#base_dir = './output/simulation_particle4/harmonic_test008/'
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

particleBeads = XMLPolyDataReader(FileName=file_list);
RenameSource("particles",particleBeads);
particleBeads.PointArrayStatus = ['id', 'type', 'vx', 'fx']


file_list = sorted(glob.glob(base_dir + '/vtk/Particles_????????_boundingBox.vtu'));
#os.chdir(path_cwd);
print("Bounding Box:");
print("file_list[0] = " + str(file_list[0]));
print("file_list[-1] = " + str(file_list[-1]));
print("len(file_list) = " + str(len(file_list)));
boundingBox1 = XMLUnstructuredGridReader(FileName=file_list);
RenameSource("bounding_box1",boundingBox1);

# create a new 'Glyph'
glyph1 = Glyph(Input=particleBeads,GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 15
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'
RenameSource("all_beads",glyph1);

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from particleBeads
particleBeadsDisplay = Show(particleBeads, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
particleBeadsDisplay.Representation = 'Wireframe'
particleBeadsDisplay.AmbientColor = [0.0, 0.3333333333333333, 0.4980392156862745]
particleBeadsDisplay.ColorArrayName = [None, '']
particleBeadsDisplay.DiffuseColor = [0.0, 0.3333333333333333, 0.4980392156862745]
particleBeadsDisplay.Opacity = 0.26
particleBeadsDisplay.LineWidth = 5.0
particleBeadsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
particleBeadsDisplay.SelectOrientationVectors = 'None'
particleBeadsDisplay.ScaleFactor = 40.400000000000006
particleBeadsDisplay.SelectScaleArray = 'None'
particleBeadsDisplay.GlyphType = 'Arrow'
particleBeadsDisplay.GlyphTableIndexArray = 'None'
particleBeadsDisplay.GaussianRadius = 2.02
particleBeadsDisplay.SetScaleArray = [None, '']
particleBeadsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
particleBeadsDisplay.OpacityArray = [None, '']
particleBeadsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
particleBeadsDisplay.DataAxesGrid = 'GridAxesRepresentation'
particleBeadsDisplay.PolarAxes = 'PolarAxesRepresentation'
particleBeadsDisplay.ScalarOpacityUnitDistance = 699.7485262578264

# show data from particleBeads
particleBeadsDisplay = Show(particleBeads, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
particleBeadsDisplay.Representation = 'Surface'
particleBeadsDisplay.ColorArrayName = [None, '']
particleBeadsDisplay.OSPRayScaleArray = 'fx'
particleBeadsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
particleBeadsDisplay.SelectOrientationVectors = 'None'
particleBeadsDisplay.ScaleFactor = 10.023999786376955
particleBeadsDisplay.SelectScaleArray = 'None'
particleBeadsDisplay.GlyphType = 'Arrow'
particleBeadsDisplay.GlyphTableIndexArray = 'None'
particleBeadsDisplay.GaussianRadius = 0.5011999893188477
particleBeadsDisplay.SetScaleArray = ['POINTS', 'fx']
particleBeadsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
particleBeadsDisplay.OpacityArray = ['POINTS', 'fx']
particleBeadsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
particleBeadsDisplay.DataAxesGrid = 'GridAxesRepresentation'
particleBeadsDisplay.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
particleBeadsDisplay.ScaleTransferFunction.Points = [-388037.97892, 0.0, 0.5, 0.0, 388037.97892, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
particleBeadsDisplay.OpacityTransferFunction.Points = [-388037.97892, 0.0, 0.5, 0.0, 388037.97892, 1.0, 0.5, 0.0]
particleBeadsDisplay.Visibility = 0;

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

glyph1Display.Visibility = 1;

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
SetActiveSource(glyph1)

# ----------------------------------------------------------------


