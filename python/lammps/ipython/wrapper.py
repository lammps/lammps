# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

################################################################################
# IPython/Jupyter Notebook additions
# Written by Richard Berger <richard.berger@outlook.com>
################################################################################

class wrapper(object):
  """ lammps API IPython Wrapper

  This is a wrapper class that provides additional methods on top of an
  existing :py:class:`lammps` instance. It provides additional methods
  that allow create and/or embed visualizations created by native LAMMPS
  commands.

  There is no need to explicitly instantiate this class. Each instance
  of :py:class:`lammps` has a :py:attr:`ipython <lammps.ipython>` property
  that returns an instance.

  :param lmp: instance of the :py:class:`lammps` class
  :type  lmp: lammps
  """
  def __init__(self, lmp):
      self.lmp = lmp

  def image(self, filename="snapshot.png", group="all", color="type", diameter="type",
            size=None, view=None, center=None, up=None, zoom=1.0, background_color="white"):
    """ Generate image using write_dump command and display it

    See :doc:`dump image <dump_image>` for more information.

    :param filename: Name of the image file that should be generated. The extension determines whether it is PNG or JPEG
    :type filename: string
    :param group: the group of atoms write_image should use
    :type group: string
    :param color: name of property used to determine color
    :type color: string
    :param diameter: name of property used to determine atom diameter
    :type diameter: string
    :param size: dimensions of image
    :type size: tuple (width, height)
    :param view: view parameters
    :type view: tuple (theta, phi)
    :param center: center parameters
    :type center: tuple (flag, center_x, center_y, center_z)
    :param up: vector pointing to up direction
    :type up: tuple (up_x, up_y, up_z)
    :param zoom: zoom factor
    :type zoom: float
    :param background_color: background color of scene
    :type background_color: string

    :return: Image instance used to display image in notebook
    :rtype: :py:class:`IPython.core.display.Image`
    """
    cmd_args = [group, "image", filename, color, diameter]

    if size is not None:
      width = size[0]
      height = size[1]
      cmd_args += ["size", width, height]

    if view is not None:
      theta = view[0]
      phi = view[1]
      cmd_args += ["view", theta, phi]

    if center is not None:
      flag = center[0]
      Cx = center[1]
      Cy = center[2]
      Cz = center[3]
      cmd_args += ["center", flag, Cx, Cy, Cz]

    if up is not None:
      Ux = up[0]
      Uy = up[1]
      Uz = up[2]
      cmd_args += ["up", Ux, Uy, Uz]

    if zoom is not None:
      cmd_args += ["zoom", zoom]

    cmd_args.append("modify backcolor " + background_color)

    self.lmp.cmd.write_dump(*cmd_args)
    from IPython.core.display import Image
    return Image(filename)

  def video(self, filename):
    """
    Load video from file

    Can be used to visualize videos from :doc:`dump movie <dump_image>`.

    :param filename: Path to video file
    :type filename: string
    :return: HTML Video Tag used by notebook to embed a video
    :rtype: :py:class:`IPython.display.HTML`
    """
    from IPython.display import HTML
    return HTML("<video controls><source src=\"" + filename + "\"></video>")
