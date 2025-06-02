#! /usr/freeware/bin/python

#
# This is dump2trj, a program written by Keir E. Novik to convert
# Lammps position dump files to Amber trajectory files.
#
# Copyright 2000, 2001 Keir E. Novik; all rights reserved.
# 
# Modified by Vikas Varshney, U Akron, 5 July 2005, as described in README
#

#============================================================

def Convert_files():
    'Handle the whole conversion process'

    print
    print 'Welcome to dump2trj, a program to convert Lammps position dump files to\nAmber trajectory format!'
    print

    Basename_list = Find_dump_files()

    for Basename in Basename_list:
        t = Trajectory()
        if t.Read_dump(Basename):
            t.Write_trj(Basename)
        del t
        print

#============================================================

def Find_dump_files():
    'Look for sets of Lammps position dump files to process'

    '''If passed something on the command line, treat it as a list of
    files to process.  Otherwise, look for *.dump in the current
    directory.
    '''

    import os, sys

    Basename_list = []

    # Extract basenames from command line
    for Name in sys.argv[1:]:
        if Name[-5:] == '.dump':
            Basename_list.append(Name[:-5])
        else:
            Basename_list.append(Name)

    if Basename_list == []:
        print 'Looking for Lammps dump files...',
        Dir_list = os.listdir('.')
        for Filename in Dir_list:
            if Filename[-5:] == '.dump':
                Basename_list.append(Filename[:-5])
        Basename_list.sort()
        if Basename_list != []:
            print 'found',
            for i in range(len(Basename_list)-1):
                print Basename_list[i] + ',',
            print Basename_list[-1] + '\n'
    
    if Basename_list == []:
        print 'none.\n'

    return Basename_list
    
#============================================================

class Snapshot:
    def __init__(self, The_trajectory):
        'Initialize the Snapshot class'
        
        self.timestep = The_trajectory.timestep
        self.atoms = The_trajectory.atoms
        self.xlo = The_trajectory.xlo
        self.xhi = The_trajectory.xhi
        self.ylo = The_trajectory.ylo
        self.yhi = The_trajectory.yhi
        self.zlo = The_trajectory.zlo
        self.zhi = The_trajectory.zhi

    #--------------------------------------------------------

    def Read_dump(self, Lines):
        'Read a snapshot (timestep) from a Lammps position dump file'

        '''Trajectory.Read_dump() will pass us only the lines we need.
        '''

        self.Atom_list = Lines
        
    #--------------------------------------------------------

    def Write_trj(self, F):
        'Write a snapshot (timestep) to an Amber trajectory file'

        '''The Atom_list must be sorted, as it may not be in order
        (for example, in a parallel Lammps simulation).
        '''

        import string

        xBOX = (self.xhi - self.xlo)
        yBOX = (self.yhi - self.ylo)
        zBOX = (self.zhi - self.zlo)

        Min = min(self.xlo, self.ylo, self.zlo)
        Max = max(self.xhi, self.yhi, self.zhi, xBOX, yBOX, zBOX)
        if Min <= -1000 or Max >= 10000:
            print '(error: coordinates too large!)'
            return

        Print_list = []

        for Line in NumericalSort(self.Atom_list):
            Item_list = string.split(Line)
            x = xBOX * (Float(Item_list[2])+Float(Item_list[5])) # Modified main box x-coordinate to actual x-coordinate
            y = yBOX * (Float(Item_list[3])+Float(Item_list[6])) # Modified main box y-coordinate to actual y-coordinate 
            z = zBOX * (Float(Item_list[4])+Float(Item_list[7])) # Modified main box z-coordinate to actual z-coordinate
            Print_list.append('%(x)8.3f' % vars())
            Print_list.append('%(y)8.3f' % vars())
            Print_list.append('%(z)8.3f' % vars())

            if len(Print_list) > 9:
                Line = ''
                for j in range(10):
                    Line = Line + Print_list[j]
                Line = Line + '\n'
                Print_list = Print_list[10:]

                try:
                    F.write(Line)
                except IOError, Detail:
                    print '(error:', Detail[1] + '!)'
                    F.close()
                    return

        if len(Print_list) > 0:
            Line = ''
            for j in range(len(Print_list)):
                Line = Line + Print_list[j]
            Line = Line + '\n'

            try:
                F.write(Line)
            except IOError, Detail:
                print '(error:', Detail[1] + '!)'
                F.close()
                return
                
        Line = '%(xBOX)8.3f%(yBOX)8.3f%(zBOX)8.3f\n' % vars()
        try:
            F.write(Line)
        except IOError, Detail:
            print '(error:', Detail[1] + '!)'
            F.close()
            return

#============================================================

class Trajectory:

    def Read_dump(self, Basename):
        'Read a Lammps position dump file'

        import string, sys

        Filename = Basename + '.dump'

        print 'Reading', Filename + '...',
        sys.stdout.flush()

        try:
            F = open(Filename)
        except IOError, Detail:
            print '(error:', Detail[1] + '!)'
            return 0

        try:
            Lines = F.readlines()
        except IOError, Detail:
            print '(error:', Detail[1] + '!)'
            F.close()
            return 0

        F.close()

        i = 0
        self.Snapshot_list = []

        # Parse the dump
        while i < len(Lines):

            if string.find(Lines[i], 'ITEM: TIMESTEP') != -1:
                # Read the timestep
                self.timestep = int(Lines[i+1])
                i = i + 2

            elif string.find(Lines[i], 'ITEM: NUMBER OF ATOMS') != -1:
                # Read the number of atoms
                self.atoms = int(Lines[i+1])
                i = i + 2

            elif string.find(Lines[i], 'ITEM: BOX BOUNDS') != -1:
                # Read the periodic box boundaries
                Item_list = string.split(Lines[i+1])
                self.xlo = Float(Item_list[0])
                self.xhi = Float(Item_list[1])
                Item_list = string.split(Lines[i+2])
                self.ylo = Float(Item_list[0])
                self.yhi = Float(Item_list[1])
                Item_list = string.split(Lines[i+3])
                self.zlo = Float(Item_list[0])
                self.zhi = Float(Item_list[1])
                i = i + 4

            elif string.find(Lines[i], 'ITEM: ATOMS') != -1:
                # Read atom positions
                self.Snapshot_list.append(Snapshot(self))
                Start = i + 1
                End = Start + self.atoms
                self.Snapshot_list[-1].Read_dump(Lines[Start:End])
                i = i + self.atoms + 1
                
            else:
                print '(error: unknown line in file!)'
                return

        print 'done.'
        return 1

    #--------------------------------------------------------

    def Write_trj(self, Basename):
        'Write an Amber trajectory file'

        import os, sys

        Filename = Basename + '.mdcrd'

        Dir_list = os.listdir('.')
        i = 1
        while Filename in Dir_list:
            Filename = Basename + `i` + '.mdcrd'
            i = i + 1
        del i

        print 'Writing', Filename + '...',
        sys.stdout.flush()

        try:
            F = open(Filename, 'w')
        except IOError, Detail:
            print '(error:', Detail[1] + '!)'
            return

        try:
            F.write(Basename + '\n')
        except IOError, Detail:
            print '(error:', Detail[1] + '!)'
            F.close()
            return

        for S in self.Snapshot_list:
            S.Write_trj(F)

        F.close()
        print 'done.'

#============================================================

def Float(s):
    'Return the string s as a float, if possible'

    try:
        x = float(s)
    except ValueError:
        if s[-1] == ',':
            s = s[:-1]
        x = float(s)

    return x

#============================================================

def NumericalSort(String_list):
    'Sort a list of strings by the integer value of the first element'

    import string

    Working_list = []
    for s in String_list:
        Working_list.append((int(string.split(s)[0]), s))

    Working_list.sort()

    Return_list = []
    for Tuple in Working_list:
        Return_list.append(Tuple[1])
    return Return_list

#============================================================

Convert_files()
