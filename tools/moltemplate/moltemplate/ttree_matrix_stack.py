# Author: Andrew Jewett (jewett.aij@gmail.com)
#         http://www.moltemplate.org
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2012, Regents of the University of California
# All rights reserved.

import random, math
from collections import deque
from array import array

#try:
#    from StringIO import StringIO
#except ImportError:
#    from io import StringIO

try:
    from .ttree_lex import InputError, ErrorLeader, OSrcLoc
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from ttree_lex import InputError, ErrorLeader, OSrcLoc



def MultMat(dest, A, B):
    """ Multiply two matrices together. Store result in "dest".

    """
    I = len(A)
    J = len(B[0])
    K = len(B)  # or len(A[0])
    for i in range(0, I):
        for j in range(0, J):
            dest[i][j] = 0.0
            for k in range(0, K):
                dest[i][j] += A[i][k] * B[k][j]

def Transpose(M):
    return [ [M[j][i] for j in range(0, len(M))]
             for i in range(0, len(M[0])) ]


def TransposeInPlace(M):
    N = len(M)
    for i in range(0, N):
        for j in range(0, i):
            M_ij = M[i][j]
            M_ji = M[j][i]
            M[i][j] = M_ji
            M[j][i] = M_ij
            


def MatToStr(M):
    strs = []
    for i in range(0, len(M)):
        for j in range(0, len(M[i])):
            strs.append(str(M[i][j]) + ' ')
        strs.append('\n')
    return(''.join(strs))


def LinTransform(dest, M, x):
    """ Multiply matrix M by 1-dimensioal array (vector) "x" (from the right).
        Store result in 1-dimensional array "dest".
        In this function, wetreat "x" and "dest" as a column vectors.
        (Not row vectors.)

    """
    I = len(M)
    J = len(x)
    for i in range(0, I):
        dest[i] = 0.0
        for j in range(0, J):
            dest[i] += M[i][j] * x[j]


def AffineTransform(dest, M, x):
    """ This function performs an affine transformation on vector "x".
        Multiply 3-dimensional vector "x" by first three columns of 3x4
        matrix M.  Add to this the final column of M.  Store result in "dest":
    dest[0] = M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2]  +  M[0][3]
    dest[1] = M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2]  +  M[1][3]
    dest[2] = M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2]  +  M[2][3]

    """
    D = len(M)
    #assert(len(M[0]) == D+1)
    for i in range(0, D):
        dest[i] = 0.0
        for j in range(0, D):
            dest[i] += M[i][j] * x[j]
        dest[i] += M[i][D]  # (translation offset stored in final column)


def AffineCompose(dest, M2, M1):
    """
    Multiplication for pairs of 3x4 matrices is technically undefined.
    However what we want to do is compose two affine transformations: M1 and M2

    3x4 matrices are used to define rotations/translations
    x' = M[0][0]*x + M[0][1]*y + M[0][2]*z + M[0][3]
    y' = M[1][0]*x + M[1][1]*y + M[1][2]*z + M[1][3]
    z' = M[2][0]*x + M[2][1]*y + M[2][2]*z + M[2][3]

    We want to create a new 3x4 matrix representing an affine transformation
    (M2 M1), defined so that when (M2 M1) is applied to vector x, the result is
    M2 (M1 x).  In other words:
        first, affine transformation M1 is applied to to x
         then, affine transformation M2 is applied to (M1 x)

    """

    D = len(M1)
    #assert(len(M1[0]) == D+1)
    #assert(len(M2[0]) == D+1)
    for i in range(0, D):
        dest[i][D] = 0.0
        for j in range(0, D + 1):
            dest[i][j] = 0.0
            for k in range(0, D):
                dest[i][j] += M2[i][k] * M1[k][j]
        dest[i][D] += M2[i][D]


def CopyMat(dest, source):
    for i in range(0, len(source)):
        for j in range(0, len(source[i])):
            dest[i][j] = source[i][j]


class AffineStack(object):
    """
    This class defines a matrix stack used to define compositions of affine
    transformations of 3 dimensional coordinates (rotation and translation).
    Affine transformations are represented using 3x4 matrices.
    (Coordinates of atoms are thought of as column vectors: [[x],[y],[z]],
     although they are represented internally in the more ordinary way [x,y,z].
     To aplly an affine transformation to a vector, multiply the vector
     by the matrix, from the left-hand side, as explained in the comments for:
     AffineTransform(dest, M, x)
    Note: The last column of the 3x4 matrix stores a translational offset.
          This bears similarity with the original OpenGL matrix stack
          http://content.gpwiki.org/index.php/OpenGL:Tutorials:Theory
    (OpenGL uses 4x4 matrices.  We don't need the final row of these matrices,
     because in OpenGL, these rows are used for perspective transformations.)
   http://en.wikipedia.org/wiki/Homogeneous_coordinates#Use_in_computer_graphics

    """

    def __init__(self):
        self.stack = None
        self.M = None
        self._tmp = None
        self.Clear()

    def Clear(self):
        self.stack = deque([])
        self.M = [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0]]  # (identity, initially)
        self._tmp = [[1.0, 0.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0]]

    def PushRight(self, M):
        # Push a copy of matrix self.M onto the stack
        #  We make no distinction between "right" and "left" here.
        #  All transformations are pushed onto the stack in the same way.
        #  (The "right" and "left" refer to whether the matrix is multiplied
        #   on the right of left hand side.  Because not all matrices need be
        #   invertible, we require that matrices be popped from the stack
        #   in the reverse order they were pushed.  This prevents the ability
        #   to push and pop matrices to either end of the stack in an arbitrary
        #   order (like append(), appendleft(), pop(), popleft()).)
        self.stack.append([[self.M[i][j] for j in range(0, len(self.M[i]))]
                           for i in range(0, len(self.M))])
        #  The "Right" and "Left" refer to whether the new matrix is multiplied
        #  on the right or left side of the culmulatie matrix product.
        # Afterwards, self._tmp = self.M * M
        AffineCompose(self._tmp, self.M, M)

        # sys.stderr.write('DEBUG: PushLeft()\n' +
        #                 MatToStr(self._tmp) + '\n  =  \n' +
        #                 MatToStr(M) +  '\n  *  \n' +
        #                 MatToStr(self.M) + '\n')

        CopyMat(self.M, self._tmp)          # Copy self._tmp into self.M

    def PushLeft(self, M):
        # Push a copy of matrix self.M onto the stack
        #  We make no distinction between right and left here.
        #  All transformations are pushed onto the stack in the same way.
        #  (The "right" and "left" refer to whether the matrix is multiplied
        #   on the right of left hand side.  Because not all matrices need be
        #   invertible, we require that matrices be popped from the stack
        #   in the reverse order they were pushed.  This prevents the ability
        #   to push and pop matrices to either end of the stack in an arbitrary
        #   order (like append(), appendleft(), pop(), popleft()).)
        self.stack.append([[self.M[i][j] for j in range(0, len(self.M[i]))]
                           for i in range(0, len(self.M))])
        #  The "Right" and "Left" refer to whether the new matrix is multiplied
        #  on the right or left side of the culmulatie matrix product.
        # Afterwards, self._tmp = M * self.M
        AffineCompose(self._tmp, M, self.M)

        # sys.stderr.write('DEBUG: PushLeft()\n' +
        #                 MatToStr(self._tmp) + '\n  =  \n' +
        #                 MatToStr(M) +  '\n  *  \n' +
        #                 MatToStr(self.M) + '\n')

        CopyMat(self.M, self._tmp)          # Copy self.tmp into self.M

    def Pop(self):
        CopyMat(self.M, self.stack.pop())
        #  (No need to return a matrix,"self.M",after popping.
        #   The caller can directly access self.M later.)
        # return self.M

    def PopRight(self):
        self.Pop()

    def PopLeft(self):
        self.Pop()

    def PushCommandsRight(self,
                          text,  # text containing affine transformation commands
                          # The next two arguments are optional:
                          src_loc=OSrcLoc(),   # for debugging
                          xcm=None):  # position of center of object
        """Generate affine transformation matrices from simple text commands
           (such as \"rotcm(90,0,0,1)\" and \"move(0,5.0,0)".
            Chains of "rotcm", "movecm", "rot", and "move" commands
            can also be strung together:
               \"rotcm(90,0,0,1).move(0,5.0,0)\"
           Commands ending in \"cm\" are carried out relative to center-of-mass
           (average position) of the object, and consequently require
           an additional argument (\"xcm\").

           """
        self.PushRight(AffineStack.CommandsToMatrix(text, src_loc, xcm))

    def PushCommandsLeft(self,
                         text,  # text containing affine transformation commands
                         # The next two arguments are optional:
                         src_loc=OSrcLoc(),   # for debugging
                         xcm=None):  # position of center of object
        self.PushLeft(AffineStack.CommandsToMatrix(text, src_loc, xcm))

    def __len__(self):
        return 1 + len(self.stack)

    @staticmethod
    def CommandsToMatrix(text,  # text containing affine transformation commands
                         src_loc=OSrcLoc(),   # for debugging
                         xcm=None):  # position of center of object
        Mdest = [[1.0, 0.0, 0.0, 0.0], [
            0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        M = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        Mtmp = [[1.0, 0.0, 0.0, 0.0], [
            0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]

        transform_commands = text.split(').')
        for transform_str in transform_commands:
            if transform_str.find('move(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[5:i_paren_close].split(',')
                if (len(args) != 3):
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires 3 numerical arguments.')
                M = [[1.0, 0.0, 0.0, float(args[0])],
                     [0.0, 1.0, 0.0, float(args[1])],
                     [0.0, 0.0, 1.0, float(args[2])]]
                AffineCompose(Mtmp, M, Mdest)
                CopyMat(Mdest, Mtmp)

            # if transform_str.find('movecm(') == 0:
            # #     assert(xcm != None)
            #    i_paren_close = transform_str.find(')')
            #    if i_paren_close == -1:
            #        i_paren_close = len(transform_str)
            #    args = transform_str[8:i_paren_close].split(',')
            #    if (len(args) != 3):
            #        raise InputError('Error near '+ErrorLeader(src_loc.infile, src_loc.lineno)+':\n'
            #                         '       Invalid command: \"'+transform_str+'\"\n'
            #                         '       This command requires 3 numerical arguments.')
            #    M = [[1.0, 0.0, 0.0, float(args[0])-(xcm[0])],
            #         [0.0, 1.0, 0.0, float(args[1])-(xcm[1])],
            #         [0.0, 0.0, 1.0, float(args[2])-(xcm[2])]]
            #    AffineCompose(Mtmp, M, Mdest)
            #    CopyMat(Mdest, Mtmp)

            elif transform_str.find('move_rand(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[10:i_paren_close].split(',')

                seed = 1
                if len(args) in (2,4,7):
                    seed = int(args[0])
                random.seed(seed)
                if len(args) == 1:
                    sigma = float(args[1])
                    x = random.gauss(0.0, sigma)
                    y = random.gauss(0.0, sigma)
                    z = random.gauss(0.0, sigma)
                elif len(args) == 2:
                    # seed = int(args[0])  this was already handled above
                    sigma = float(args[1])
                    x = random.gauss(0.0, sigma)
                    y = random.gauss(0.0, sigma)
                    z = random.gauss(0.0, sigma)
                elif len(args) == 3:
                    x = random.gauss(0.0, float(args[0]))
                    y = random.gauss(0.0, float(args[1]))
                    z = random.gauss(0.0, float(args[2]))
                elif len(args) == 4:
                    # seed = int(args[0])   this was already handled above
                    x = random.gauss(0.0, float(args[1]))
                    y = random.gauss(0.0, float(args[2]))
                    z = random.gauss(0.0, float(args[3]))
                elif len(args) == 6:
                    x_min = float(args[0])
                    x_max = float(args[1])
                    y_min = float(args[2])
                    y_max = float(args[3])
                    z_min = float(args[4])
                    z_max = float(args[5])
                    x = x_min + (x_max - x_min)*(random.random()-0.5)
                    y = y_min + (y_max - y_min)*(random.random()-0.5)
                    z = z_min + (z_max - z_min)*(random.random()-0.5)
                        
                elif len(args) == 7:
                    # seed = int(args[0])  this was already handled above
                    x_min = float(args[1])
                    x_max = float(args[2])
                    y_min = float(args[3])
                    y_max = float(args[4])
                    z_min = float(args[5])
                    z_max = float(args[6])
                    x = x_min + (x_max - x_min)*(random.random()-0.5)
                    y = y_min + (y_max - y_min)*(random.random()-0.5)
                    z = z_min + (z_max - z_min)*(random.random()-0.5)
                else:
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 1, 2, 3, 4, 6 or 7 numerical arguments.  Either:\n'
                                     '           move_rand(gauss_sigma) or\n'
                                     '           move_rand(seed, gauss_sigma) or\n'
                                     '           move_rand(gauss_sigma_x, gauss_sigma_y, gauss_sigma_z) or\n'
                                     '           move_rand(seed, gauss_sigma_x, gauss_sigma_y, gauss_sigma_z) or\n'
                                     '           move_rand(x_min, x_max, y_min, y_max, z_min, z_max) or\n'
                                     '           move_rand(seed, x_min, x_max, y_min, y_max, z_min, z_max)\n')

                M = [[1.0, 0.0, 0.0, x],
                     [0.0, 1.0, 0.0, y],
                     [0.0, 0.0, 1.0, z]]
                AffineCompose(Mtmp, M, Mdest)
                CopyMat(Mdest, Mtmp)

            elif transform_str.find('rot(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[4:i_paren_close].split(',')
                center_v = None
                if (len(args) == 7):
                    center_v = [float(args[4]), float(args[5]), float(args[6])]
                elif (len(args) != 4):
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 4 or 7 numerical arguments.  Either:\n'
                                     '           rot(angle, axisX, axisY, axiZ)  or \n'
                                     '           rot(angle, axisX, axisY, axiZ, centerX, centerY, centerZ)')
                M[0][3] = 0.0  # RotMatAXYZ() only modifies 3x3 submatrix of M
                M[1][3] = 0.0  # The remaining final column must be zeroed by hand
                M[2][3] = 0.0
                RotMatAXYZ(M,
                           float(args[0]) * math.pi / 180.0,
                           float(args[1]),
                           float(args[2]),
                           float(args[3]))
                if (center_v == None):
                    AffineCompose(Mtmp, M, Mdest)
                    CopyMat(Mdest, Mtmp)
                else:
                    # Move "center_v" to the origin
                    moveCentToOrig = [[1.0, 0.0, 0.0, -center_v[0]],
                                      [0.0, 1.0, 0.0, -center_v[1]],
                                      [0.0, 0.0, 1.0, -center_v[2]]]
                    AffineCompose(Mtmp, moveCentToOrig, Mdest)
                    CopyMat(Mdest, Mtmp)
                    # Rotate the coordinates (relative to the origin)
                    AffineCompose(Mtmp, M, Mdest)  # M is the rotation matrix
                    CopyMat(Mdest, Mtmp)
                    # Move the origin back to center_v
                    moveCentBack = [[1.0, 0.0, 0.0, center_v[0]],
                                    [0.0, 1.0, 0.0, center_v[1]],
                                    [0.0, 0.0, 1.0, center_v[2]]]
                    AffineCompose(Mtmp, moveCentBack, Mdest)
                    CopyMat(Mdest, Mtmp)

            elif ((transform_str.find('quat(') == 0) or
                  (transform_str.find('quatT(') == 0)):
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[5:i_paren_close].split(',')
                center_v = None
                if (len(args) == 7):
                    center_v = [float(args[4]), float(args[5]), float(args[6])]
                elif (len(args) != 4):
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 4 or 7 numerical arguments.  Either:\n'
                                     '           rot(angle, axisX, axisY, axiZ)  or \n'
                                     '           rot(angle, axisX, axisY, axiZ, centerX, centerY, centerZ)')
                M[0][3] = 0.0  # RotMatAXYZ() only modifies 3x3 submatrix of M
                M[1][3] = 0.0  # The remaining final column must be zeroed by hand
                M[2][3] = 0.0
                q = (float(args[0]),
                     float(args[1]),
                     float(args[2]),
                     float(args[3]))
                Quaternion2Matrix(q, M)
                if (transform_str.find('quatT(') == 0):
                    TransposeInPlace(M)
                if (center_v == None):
                    AffineCompose(Mtmp, M, Mdest)
                    CopyMat(Mdest, Mtmp)
                else:
                    # Move "center_v" to the origin
                    moveCentToOrig = [[1.0, 0.0, 0.0, -center_v[0]],
                                      [0.0, 1.0, 0.0, -center_v[1]],
                                      [0.0, 0.0, 1.0, -center_v[2]]]
                    AffineCompose(Mtmp, moveCentToOrig, Mdest)
                    CopyMat(Mdest, Mtmp)
                    # Rotate the coordinates (relative to the origin)
                    AffineCompose(Mtmp, M, Mdest)  # M is the rotation matrix
                    CopyMat(Mdest, Mtmp)
                    # Move the origin back to center_v
                    moveCentBack = [[1.0, 0.0, 0.0, center_v[0]],
                                    [0.0, 1.0, 0.0, center_v[1]],
                                    [0.0, 0.0, 1.0, center_v[2]]]
                    AffineCompose(Mtmp, moveCentBack, Mdest)
                    CopyMat(Mdest, Mtmp)

            # # elif transform_str.find('rotcm(') == 0:
            # #     assert(xcm != None)
            # #     i_paren_close = transform_str.find(')')
            # #     if i_paren_close == -1:
            # #         i_paren_close = len(transform_str)
            # #     args = transform_str[6:i_paren_close].split(',')
            # #     if (len(args) != 4):
            # #         raise InputError('Error near '+ErrorLeader(src_loc.infile, src_loc.lineno)+':\n'
            # #                          '       Invalid command: \"'+transform_str+'\"\n'
            # #                          '       This command requires 4 numerical arguments.')
            # #
            # #     moveCMtoOrig = [[1.0, 0.0, 0.0, -xcm[0]],
            # #                     [0.0, 1.0, 0.0, -xcm[1]],
            # #                     [0.0, 0.0, 1.0, -xcm[2]]]
            # #     AffineCompose(Mtmp, moveCMtoOrig, Mdest)
            # #     CopyMat(Mdest, Mtmp)
            # #     M[0][3] = 0.0#RotMatAXYZ() only modifies 3x3 submatrix of M
            # #     M[1][3] = 0.0#The remaining final column must be zeroed by hand
            # #     M[2][3] = 0.0
            # #     RotMatAXYZ(M,
            # #                float(args[0])*math.pi/180.0,
            # #                float(args[1]),
            # #                float(args[2]),
            # #                float(args[3]))
            # #     AffineCompose(Mtmp, M, Mdest)
            # #     CopyMat(Mdest, Mtmp)
            # #     moveCmBack = [[1.0, 0.0, 0.0, xcm[0]],
            # #                   [0.0, 1.0, 0.0, xcm[1]],
            # #                   [0.0, 0.0, 1.0, xcm[2]]]
            # #     AffineCompose(Mtmp, moveCmBack, Mdest)
            # #     CopyMat(Mdest, Mtmp)

            elif transform_str.find('rot_rand(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[9:i_paren_close].split(',')

                seed = 1
                if len(args) in (2,6):
                    seed = int(args[0])
                random.seed(seed)
                raxis = [0.0, 0.0, 0.0]
                if len(args) < 5:
                    # choose a random rotation axis
                    raxis_len = 0.0
                    while (not ((0.01<raxis_len) and (raxis_len <= 1.0))):
                        raxis = [-1+2.0*(random.random()-0.5) for d in range(0,3)]
                        raxis_len = math.sqrt(raxis[0]**2 + raxis[1]**2 + raxis[2]**2)
                    for d in range(0,3):
                        raxis[d] /= raxis_len
                if len(args) == 0:
                    angle_min = angle_max = 2*math.pi
                elif len(args) == 1:
                    angle_min = 0.0
                    angle_max = float(args[0]) * math.pi / 180.0,
                elif len(args) == 5:
                    angle_min = float(args[0])
                    angle_max = float(args[1])
                    raxis[0]  = float(args[2])
                    raxis[1]  = float(args[3])
                    raxis[2]  = float(args[4])
                elif len(args) == 6:
                    seed      = int(args[0])
                    angle_min = float(args[1])
                    angle_max = float(args[2])
                    raxis[0]  = float(args[3])
                    raxis[1]  = float(args[4])
                    raxis[2]  = float(args[5])
                else:
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 0, 1, 2, 5 or 6 numerical arguments. Either:\n'
                                     '           rot_rand()  or \n'
                                     '           rot_rand(delta_angle)  or \n'
                                     '           rot_rand(seed, delta_angle)  or \n'
                                     '           rot_rand(angle_min, angle_max, axisX, axisY, axiZ) or\n'
                                     '           rot_rand(seed, angle_min, angle_max, axisX, axisY, axiZ)')

                angle = angle_min + (angle_max - angle_min)*(random.random() - 0.5)

                M[0][3] = 0.0  # RotMatAXYZ() only modifies 3x3 submatrix of M
                M[1][3] = 0.0  # The remaining final column must be zeroed by hand
                M[2][3] = 0.0
                RotMatAXYZ(M,
                           angle,
                           raxis[0], raxis[1], raxis[2])
                AffineCompose(Mtmp, M, Mdest)
                CopyMat(Mdest, Mtmp)


            elif transform_str.find('rotvv(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[6:i_paren_close].split(',')
                center_v = None
                if (len(args) == 9):
                    center_v = [float(args[6]), float(args[7]), float(args[8])]
                elif (len(args) != 6):
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 6 or 9 numerical arguments.  Either:\n'
                                     '           rotvv(Xold,Yold,Zold,Xnew,Ynew,Znew)  or \n'
                                     '           rotvv(Xold,Yold,Zold,Xnew,Ynew,Znew,centerX,centerY,centerZ)')
                M[0][3] = 0.0  # RotMatXYZXYZ() only modifies 3x3 submatrix of M
                M[1][3] = 0.0  # The remaining final column must be zeroed by hand
                M[2][3] = 0.0
                RotMatXYZXYZ(M,
                             float(args[0]),
                             float(args[1]),
                             float(args[2]),
                             float(args[3]),
                             float(args[4]),
                             float(args[5]))
                if (center_v == None):
                    AffineCompose(Mtmp, M, Mdest)
                    CopyMat(Mdest, Mtmp)
                else:
                    # Move "center_v" to the origin
                    moveCentToOrig = [[1.0, 0.0, 0.0, -center_v[0]],
                                      [0.0, 1.0, 0.0, -center_v[1]],
                                      [0.0, 0.0, 1.0, -center_v[2]]]
                    AffineCompose(Mtmp, moveCentToOrig, Mdest)
                    CopyMat(Mdest, Mtmp)
                    # Rotate the coordinates (relative to the origin)
                    AffineCompose(Mtmp, M, Mdest)  # M is the rotation matrix
                    CopyMat(Mdest, Mtmp)
                    # Move the origin back to center_v
                    moveCentBack = [[1.0, 0.0, 0.0, center_v[0]],
                                    [0.0, 1.0, 0.0, center_v[1]],
                                    [0.0, 0.0, 1.0, center_v[2]]]
                    AffineCompose(Mtmp, moveCentBack, Mdest)
                    CopyMat(Mdest, Mtmp)

            elif transform_str.find('scale(') == 0:
                i_paren_close = transform_str.find(')')
                if i_paren_close == -1:
                    i_paren_close = len(transform_str)
                args = transform_str[6:i_paren_close].split(',')

                if (len(args) == 1):
                    scale_v = [float(args[0]), float(args[0]), float(args[0])]
                    center_v = [0.0, 0.0, 0.0]
                elif (len(args) == 3):
                    scale_v = [float(args[0]), float(args[1]), float(args[2])]
                    center_v = [0.0, 0.0, 0.0]
                elif (len(args) == 4):
                    scale_v = [float(args[0]), float(args[0]), float(args[0])]
                    center_v = [float(args[1]), float(args[2]), float(args[3])]
                elif (len(args) == 6):
                    scale_v = [float(args[0]), float(args[1]), float(args[2])]
                    center_v = [float(args[3]), float(args[4]), float(args[5])]
                else:
                    raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                     '       Invalid command: \"' + transform_str + '\"\n'
                                     '       This command requires either 1, 3, 4, or 6 numerical arguments. Either:\n'
                                     '           scale(ratio), or \n'
                                     '           scale(ratioX, ratioY, ratioZ),\n'
                                     '           scale(ratio, centerX, centerY, centerZ), or\n'
                                     '           scale(ratioX, ratioY, ratioZ, centerX, centerY, centerZ)')

                ScaleMat(M, scale_v)

                # Now worry about translation:
                for d in range(0, 3):
                    M[d][3] = center_v[d] * (1.0 - scale_v[d])

                AffineCompose(Mtmp, M, Mdest)
                CopyMat(Mdest, Mtmp)

            # # elif transform_str.find('scalecm(') == 0:
            # #     assert(xcm != None)
            # #     i_paren_close = transform_str.find(')')
            # #     if i_paren_close == -1:
            # #         i_paren_close = len(transform_str)
            # #     args = transform_str[8:i_paren_close].split(',')
            # #
            # #     moveCMtoOrig = [[1.0, 0.0, 0.0, -xcm[0]],
            # #                     [0.0, 1.0, 0.0, -xcm[1]],
            # #                     [0.0, 0.0, 1.0, -xcm[2]]]
            # #     AffineCompose(Mtmp, moveCMtoOrig, Mdest)
            # #     CopyMat(Mdest, Mtmp)
            # #
            # #     M[0][3] = 0.0 #ScaleMat() only modifies 3x3 submatrix of M
            # #     M[1][3] = 0.0 #The remaining final column must be zeroed by hand
            # #     M[2][3] = 0.0
            # #     if (len(args) == 1):
            # #         ScaleMat(M, args[0])
            # #     elif (len(args) == 3):
            # #         ScaleMat(M, args)
            # #     else:
            # #         raise InputError('Error near '+ErrorLeader(src_loc.infile, src_loc.lineno)+':\n'
            # #                          '       Invalid command: \"'+transform_str+'\"\n'
            # #                          '       This command requires either 1 or 3 numerical arguments.')
            # #
            # #     AffineCompose(Mtmp, M, Mdest)
            # #     CopyMat(Mdest, Mtmp)
            # #     moveCmBack = [[1.0, 0.0, 0.0, xcm[0]],
            # #                   [0.0, 1.0, 0.0, xcm[1]],
            # #                   [0.0, 0.0, 1.0, xcm[2]]]
            # #     AffineCompose(Mtmp, moveCmBack, Mdest)
            # #     CopyMat(Mdest, Mtmp)

            #elif transform_str.find('read_xyz(') == 0:
            #    i_paren_close = transform_str.find(')')
            #    if i_paren_close == -1:
            #        i_paren_close = len(transform_str)
            #    args = transform_str[4:i_paren_close].split(',')
            #    if (len(args) != 1):
            #        raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
            #                         '       Invalid command: \"' + transform_str + '\"\n'
            #                         '       This command expects the name of a file in XYZ format.\n')
            #    file_name = args[0]
            #    if (not file_name in coord_files):
            #        f = open(file_name, 'r')
            #        f.close()
            #        f.readline() # skip the first 2 lines
            #        f.readline() # of an .xyz file (header)
            #        for line in f:
            #            tokens = line.split()
            #            if (len(tokens) != 3) and (len(tokens) != 0):
            #                raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
            #                                 '       Invalid command: \"' + transform_str + '\"\n'
            #                                 '       This command expects the name of a file in XYZ format.\n')
            #            crds.append((float(tokens[0]),
            #                         float(tokens[1]),
            #                         float(tokens[2]))
            #        self.coord_files[file_name] = crds
            #    else:
            #        crds = self.coord_files[file_name]

            else:
                raise InputError('Error near ' + ErrorLeader(src_loc.infile, src_loc.lineno) + ':\n'
                                 '       Unknown transformation command: \"' + transform_str + '\"\n')

        return Mdest


class MultiAffineStack(object):

    def __init__(self, which_stack=None):
        self.tot_stack = None
        self.stack_lookup = None
        self.stack_keys = None
        self.stacks = None
        self.M = None
        self.error_if_substack_empty = False
        self.coord_files = {}
        self.Clear()

    def Clear(self):
        self.tot_stack = AffineStack()
        self.stack_lookup = {}
        self.stack_keys = deque([])
        self.stacks = deque([])
        self.M = self.tot_stack.M
        self.error_if_substack_empty = False
        self.coord_files = {}

    def _Update(self):
        self.tot_stack.Clear()
        self.M = self.tot_stack.M
        assert(len(self.stacks) > 0)
        for stack in self.stacks:
            self.tot_stack.PushRight(stack.M)

    def PushStack(self, which_stack):
        stack = AffineStack()
        self.stack_keys.append(which_stack)
        self.stack_lookup[which_stack] = stack
        self.stacks.append(stack)
        self.tot_stack.PushRight(stack.M)

    def PopStack(self):
        assert(len(self.stacks) > 0)
        self.tot_stack.PopRight()
        which_stack = self.stack_keys.pop()
        del self.stack_lookup[which_stack]
        self.stacks.pop()

    def Push(self, M, which_stack=None, right_not_left=True):
        if len(self.stacks) == 0:
            self.PushStack(which_stack)
        if which_stack == None:
            stack = self.stacks[-1]
            if right_not_left:
                # This should copy the matrix M into stack.M
                stack.PushRight(M)
            else:
                stack.PushLeft(M)
        else:
            stack = self.stack_lookup[which_stack]
            if right_not_left:
                stack.PushRight(M)
            else:
                stack.PushLeft(M)
        if stack == self.stacks[-1]:
            self.tot_stack.PopRight()  # Replace the last matrix on self.tot_stack
            # Note: Always use tot_stack.PopRight (even if
            # right_not_left=False)
            self.tot_stack.PushRight(stack.M)  # with the the updated version.
            # Note: We could call self._Update(M) here, but that is slower.
        else:
            self._Update()

    def PushRight(self, M, which_stack=None):
        self.Push(M, which_stack, right_not_left=True)

    def PushLeft(self, M, which_stack=None):
        self.Push(M, which_stack, right_not_left=False)

    def PushCommandsRight(self,
                          text,  # text containing affine transformation commands
                          # The next two arguments are optional:
                          src_loc=OSrcLoc(),   # for debugging
                          xcm=None,
                          which_stack=None):  # position of center of object
        """Generate affine transformation matrices from simple text commands
           (such as \"rotcm(90,0,0,1)\" and \"move(0,5.0,0)".
            Chains of "rotcm", "movecm", "rot", and "move" commands
            can also be strung together:
               \"rotcm(90,0,0,1).move(0,5.0,0)\"
           Commands ending in \"cm\" are carried out relative to center-of-mass
           (average position) of the object, and consequently require
           an additional argument (\"xcm\").

           """
        self.PushRight(AffineStack.CommandsToMatrix(text, src_loc, xcm),
                       which_stack)

    def PushCommandsLeft(self,
                         text,  # text containing affine transformation commands
                         # The next two arguments are optional:
                         src_loc=OSrcLoc(),   # for debugging
                         xcm=None,            # position of center of object
                         which_stack=None):
        self.PushLeft(AffineStack.CommandsToMatrix(text, src_loc, xcm),
                      which_stack)

    def Pop(self, which_stack=None, right_not_left=True):
        #empty_stack_error = False

        if which_stack == None:
            stack = self.stacks[-1]
            if len(stack) >= 1:
                if right_not_left:
                    stack.PopRight()
                else:
                    stack.PopLeft()
                # Note: We could call self._Update(M) here, but that is slower
                self.tot_stack.PopRight()  # Replace the last matrix on self.tot_stack
                # Note: Always use tot_stack.PopRight (even if
                # right_not_left=False)
                # with the the updated version.
                self.tot_stack.PushRight(stack.M)
            else:
                assert(False)
            # OPTIONAL CODE BELOW AUTOMATICALLY INVOKES self.PopStack() WHEN
            # THE stacks[-1].stack IS EMPTY.  PROBABLY DOES NOT WORK.  IGNORE
            #    if (not self.error_if_substack_empty):
            #        if right_not_left:
            #            assert(len(self.stacks) > 0)
            #            self.PopStack()
            #        else:
            #            assert(False)
            #    else:
            #        empty_stack_error = True

        else:
            stack = self.stack_lookup[which_stack]
            if len(stack) > 1:
                if right_not_left:
                    stack.PopRight()
                else:
                    stack.PopLeft()
                self._Update()
            else:
                assert(False)
                #empty_stack_error = True

    def PopRight(self, which_stack=None):
        self.Pop(which_stack, right_not_left=True)

    def PopLeft(self, which_stack=None):
        self.Pop(which_stack, right_not_left=True)


def ScaleMat(dest, scale):
    for i in range(0, len(dest)):
        for j in range(0, len(dest[i])):
            dest[i][j] = 0.0
    if ((type(scale) is float) or (type(scale) is int)):
        for i in range(0, len(dest)):
            dest[i][i] = scale
    else:
        for i in range(0, len(dest)):
            dest[i][i] = scale[i]


def RotMatAXYZ(dest, angle, axis_x, axis_y, axis_z):

    r = math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z)

    X = 1.0
    Y = 0.0
    Z = 0.0
    if r > 0.0:  # check for non-sensical input
        X = axis_x / r
        Y = axis_y / r
        Z = axis_z / r
    else:
        angle = 0.0

    # angle *= math.pi/180.0 # "angle" is assumed to be in degrees
    #    on second thought, let the caller worry about angle units.
    c = math.cos(angle)
    s = math.sin(angle)

    dest[0][0] = X * X * (1 - c) + c
    dest[1][1] = Y * Y * (1 - c) + c
    dest[2][2] = Z * Z * (1 - c) + c

    dest[0][1] = X * Y * (1 - c) - Z * s
    dest[0][2] = X * Z * (1 - c) + Y * s

    dest[1][0] = Y * X * (1 - c) + Z * s
    dest[2][0] = Z * X * (1 - c) - Y * s

    dest[1][2] = Y * Z * (1 - c) - X * s
    dest[2][1] = Z * Y * (1 - c) + X * s

    #   formula from these sources:
    # http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
    # also check
    # http://www.manpagez.com/man/3/glRotate/
    #   some pdb test commands:
    # from lttree_matrixstack import *
    # r = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
    # RotMatAXYZ(r, 90.0, 0.0, 0.0, 1.0)


def Quaternion2Matrix(q, M):
    "convert a quaternion (q) to a 3x3 rotation matrix (M)"""

    M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3])
    M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3])
    M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3])
    M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
    M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
    M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
    M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
    M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
    M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);



def Matrix2Quaternion(M, q):
    """convert a 3x3 rotation matrix (M) to a quaternion (q)
       http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    """
    tr = M[0][0] + M[1][1] + M[2][2]
    if tr > 0:
        S = math.sqrt(tr+1.0) * 2                        # S=4*qw 
        qw = 0.25 * S
        qx = (M[2][1] - M[1][2]) / S
        qy = (M[0][2] - M[2][0]) / S
        qz = (M[1][0] - M[0][1]) / S
    elif (M[0][0] > M[1][1]) and (M[0][0] > M[2][2]):
        S = math.sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2   # S=4*qx 
        qw = (M[2][1] - M[1][2]) / S
        qx = 0.25 * S
        qy = (M[0][1] + M[1][0]) / S
        qz = (M[0][2] + M[2][0]) / S
    elif (M[1][1] > M[2][2]):
        S = math.sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2   # S=4*qy
        qw = (M[0][2] - M[2][0]) / S
        qx = (M[0][1] + M[1][0]) / S
        qy = 0.25 * S
        qz = (M[1][2] + M[2][1]) / S
    else:
        S = math.sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2   # S=4*qz
        qw = (M[1][0] - M[0][1]) / S
        qx = (M[0][2] + M[2][0]) / S
        qy = (M[1][2] + M[2][1]) / S
        qz = 0.25 * S
    q[0] = qw
    q[1] = qx
    q[2] = qy
    q[3] = qz


def MultQuat(dest, q1, q2):
    """ multiply 2 quaternions and store the result in "qdest"
      (q1[0] + i*q1[1] + j*q1[2] + k*q1[3])
       *
      (q2[0] + i*q2[1] + j*q3[2] + k*q3[3])
    https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
    """
    dest[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3]
    dest[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2]
    dest[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1]
    dest[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]



def CrossProd(dest, A, B):
    dest[0] = (A[1] * B[2] - B[1] * A[2])
    dest[1] = (A[2] * B[0] - B[2] * A[0])
    dest[2] = (A[0] * B[1] - B[0] * A[1])


def DotProd(A, B):
    c = 0.0
    for d in range(0, len(A)):
        c += A[d] * B[d]
    return c


def Length(A):
    L = 0.0
    for x in A:
        L += x * x
    return math.sqrt(L)


def Normalize(dest, source):
    assert(len(dest) == len(source))
    L = Length(source)
    for d in range(0, len(source)):
        dest[d] = source[d] / L


def RotMatXYZXYZ(dest,
                 xold, yold, zold,
                 xnew, ynew, znew):
    A = [xold, yold, zold]
    B = [xnew, ynew, znew]
    axis = [0.0, 0.0, 0.0]
    CrossProd(axis, A, B)
    La = Length(A)
    Lb = Length(B)
    Lc = Length(axis)
    sinAng = Lc / (La * Lb)
    cosAng = DotProd(A, B) / (La * Lb)
    if Lc > 0.0:
        Normalize(axis, axis)
        angle = math.atan2(sinAng, cosAng)
    else:
        axis = [1.0, 0.0, 0.0]
        angle = 0.0
    RotMatAXYZ(dest, angle, axis[0], axis[1], axis[2])
