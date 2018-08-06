#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, io

try:
    from .dump2data import *
    from .extract_lammps_data import lammps_data_sections
except (ImportError, SystemError, ValueError):
    # not installed as a package
    from dump2data import *
    from extract_lammps_data import lammps_data_sections

g_program_name = 'raw2data.py'
g_date_str = '2016-12-21'
g_version_str = 'v0.44.0'

    #######  Main Code Below: #######
def main():
    sys.stderr.write(g_program_name + ' ' + g_version_str + ' ' + g_date_str)
    sys.stderr.write('\n')

    try:
        data_settings = DataSettings()
        misc_settings = MiscSettings()
        misc_settings.multi = False

        # First process any arguments which are specific to "raw2data.py"
        # (and remove them before passing them to dump2data.ParseArgs())
        sort_data_file_by_atom_id = False
        argv = [arg for arg in sys.argv]
        i = 1
        while i < len(argv):
            if argv[i].lower() == '-ignore-atom-id':
                sort_data_file_by_atom_id = True
                del argv[i:i+1]
            if argv[i].lower() == '-sort':
                sort_data_file_by_atom_id = False
                del argv[i:i+1]
            else:
                i += 1
            #(perhaps later I'll add some additional argumets)

        warning_strings = []
        ParseArgs(argv,
                  misc_settings,
                  data_settings,
                  warning_strings)

        frame_atom_order = []

        # The atoms in the "Atoms" section of the data file might be out 
        # of order. Extract the text from that section, and figure out the 
        # atom-ID assigned to each atom (number in the first column).
        # Then assign the these atom-ID numbers to entries in the 
        # coordinates list (frame_coords) in the same order.  We want the
        # lines of text in the coordinate file to pasted to the end of
        # the lines of text in the "Atoms" section of the data file
        # in the same order, regardless of the atom-ID numbers in that file.
        # Counterintuitively, the only way to do that is to assign the same
        # atom-ID numbers to the coordinate list in frame_coords".  So we
        # have to extract those numbers from the data file.
        in_atoms_section = False
        num_blank_lines = 0
        for line in data_settings.contents:
            ic = line.find('#')
            line = line[:ic]    #(this also removes the newline character)
            tokens = line.strip().split()
            if line.strip() == 'Atoms':
                in_atoms_section = True
            elif line.strip() in lammps_data_sections:
                in_atoms_section = False
            elif in_atoms_section:
                if len(tokens) > 0:
                    if sort_data_file_by_atom_id:
                        atomid = tokens[0]
                        frame_atom_order.append(atomid)
                        #sys.stderr.write('atomid=\"'+str(atomid)+'\"\n')
                    else:
                        frame_atom_order.append(str(len(frame_atom_order)+1))
                else:
                    num_blank_lines += 1
                    if num_blank_lines > 1:
                        in_atoms_section = False

        frame_coords = defaultdict(list)
        frame_coords_ixiyiz = defaultdict(list)
        frame_vects = defaultdict(list)
        frame_velocities = defaultdict(list)
        frame_xlo_str = frame_xhi_str = None
        frame_ylo_str = frame_yhi_str = None
        frame_zlo_str = frame_zhi_str = None
        frame_xy_str = frame_xz_str = frame_yz_str = None
        frame_natoms = -1
        frame_timestep_str = ''
        i_atomid = i_atomtype = i_molid = -1
        i_x = i_y = i_z = i_xu = i_yu = i_zu = -1
        i_xs = i_ys = i_zs = i_xsu = i_ysu = i_zsu = -1

        dump_column_names = []

        #num_frames_in = -1
        num_frames_out = 0
        finished_reading_frame = False
        read_last_frame = False

        in_coord_file = sys.stdin
        #in_coord_file = open('tmp_atom_coords.dat','r')

        read_last_frame = False
        while True:

            line = '\n'
            while (line != '') and (line.strip() == ''):
                line = in_coord_file.readline()

            if line == '':  # if EOF
                break

            #frame_vects = defaultdict(list)
            frame_coords = defaultdict(list)
            while line.strip() != '':
                n_crds = len(frame_coords)
                #sys.stdout.write("n_crds="+str(n_crds)+": \""+line.strip()+"\"\n")
                frame_coords[frame_atom_order[n_crds]] = line.split()
                #frame_vects[frame_atom_order[n_crds]] = ['0.0','0.0','0.0']
                line = in_coord_file.readline()

            # Check to see if there are any blank lines at this location in the file
            # If there are, it means we are at a new "frame" (in the trajectory).
            # Skip over any blank line(s) separating this frame from the next frame
            # so that the next time we enter the loop, we are at the beginning
            # of a new frame.
            while (line != '') and (line.strip() == ''):
                line = in_coord_file.readline()

            # Parse the DATA file specified by the user
            # and replace appropriate lines or fields with
            # the corresponding text from the input file.
            if misc_settings.multi:
                out_file_name = data_settings.file_name + '.'\
                    + str(num_frames_out)
                sys.stderr.write('  (creating file \"' + out_file_name + '\")\n')
                out_file = open(out_file_name, 'w')
            else:
                out_file = sys.stdout

            WriteFrameToData(out_file,
                             None,
                             misc_settings,
                             data_settings,
                             None,
                             frame_natoms,
                             frame_coords,
                             frame_coords_ixiyiz,
                             frame_vects,
                             frame_velocities,
                             None,
                             None,
                             frame_xlo_str, frame_xhi_str,
                             frame_ylo_str, frame_yhi_str,
                             frame_zlo_str, frame_zhi_str,
                             frame_xy_str, frame_xz_str, frame_yz_str)

    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return


if __name__ == '__main__':
    main()
