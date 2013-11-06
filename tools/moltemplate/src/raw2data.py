#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dump2data import *

g_program_name = 'raw2data.py'
g_date_str     = '2013-9-19'
g_version_str  = 'v0.42'

#######  Main Code Below: #######
sys.stderr.write(g_program_name+' '+g_version_str+' '+g_date_str)
#if sys.version < '3':
#    sys.stderr.write('  (python version < 3)\n')
#else:
sys.stderr.write('\n')

try:
    data_settings = DataSettings()
    misc_settings = MiscSettings()
    misc_settings.multi = False

    warning_strings = []
    ParseArgs(sys.argv, 
              misc_settings, 
              data_settings, 
              warning_strings)

    frame_coords = defaultdict(list)
    frame_coords_ixiyiz = defaultdict(list)
    frame_vects = defaultdict(list)
    frame_velocities = defaultdict(list)
    frame_xlo_str = frame_xhi_str = None
    frame_ylo_str = frame_yhi_str = None
    frame_zlo_str = frame_zhi_str = None
    frame_xy_str  = frame_xz_str  = frame_yz_str = None
    frame_natoms  = -1
    frame_timestep_str = ''
    i_atomid = i_atomtype = i_molid = -1
    i_x  = i_y  = i_z  = i_xu  = i_yu  = i_zu  = -1
    i_xs = i_ys = i_zs = i_xsu = i_ysu = i_zsu = -1

    dump_column_names = []

    #num_frames_in = -1
    num_frames_out = 0
    finished_reading_frame = False
    read_last_frame = False

    #in_coord_file = open('traj.raw','r')
    #in_coord_file = open('tmp_atom_coords.dat','r')
    in_coord_file = sys.stdin
    read_last_frame = False
    while True:

        line = '\n'
        while (line != '') and (line.strip() == ''):
            line = in_coord_file.readline()

        if line == '': # if EOF
            break

        frame_coords = defaultdict(list)
        while line.strip() != '':
            n_crds = len(frame_coords)
            #sys.stdout.write("n_crds="+str(n_crds)+": \""+line.strip()+"\"\n")
            frame_coords[str(n_crds+1)] = line.split()
            line = in_coord_file.readline()

        while (line != '') and (line.strip() == ''):
            line = in_coord_file.readline()

        # Parse the DATA file specified by the user
        # and replace appropriate lines or fields with
        # the corresponding text from the input file.
        if misc_settings.multi:
            out_file_name = data_settings.file_name + '.'\
                + str(num_frames_out)
            sys.stderr.write('  (creating file \"'+out_file_name+'\")\n')
            out_file = open(out_file_name, 'w')
        else:
            out_file = sys.stdout

        WriteFrameToData(out_file,
                         None,
                         misc_settings,
                         data_settings,
                         frame_natoms,
                         frame_coords,
                         frame_coords_ixiyiz,
                         frame_vects,
                         frame_velocities,
                         frame_xlo_str, frame_xhi_str, 
                         frame_ylo_str, frame_yhi_str, 
                         frame_zlo_str, frame_zhi_str,
                         frame_xy_str, frame_xz_str, frame_yz_str)

except (ValueError, InputError) as err:
    sys.stderr.write('\n'+str(err)+'\n')
    sys.exit(-1)
