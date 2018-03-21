function [varargout] = readdump_one(varargin)
% Read LAMMPS dump file one timestep at a time
% Input 
%       Dump file name with path 
%       Starting file pointer position in dump file
%       Number of columns in the dump file
% Output is in the form of a structure with following variables
% .timestep     --> Vector containing all time steps
% .Natoms       --> Vector containing number of atoms at each time step
% .x_bound      --> [t,2] array with xlo,xhi at each time step
% .y_bound      --> [t,2] array with ylo,yhi at each time step
% .z_bound      --> [t,2] array with zlo,zhi at each time step
% .atom_data    --> 2 dimensional array with data 
% .position     --> file pointer for reading next time 
%               
% Example
%       data = readdump_one('dump.LAMMPS',0,5); 
%           Reads the first timestep in the file dump.LAMMPS
%           Specify position = -1 if only the last dump step in the
%               file is needed
%
% See also readdump, scandump
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.

try
    dump = fopen(varargin{1},'r');
catch
    error('Dumpfile not found!');
end
position = varargin{2}; % from beg of file
ncol = varargin{3}; %number of columns

i=1;
t=0;
done = 0;
last_status = 0;
if position ~= -1
    fseek(dump,position,'bof');
else
    last_status = 1;
end
while done == 0 & last_status == 0
    id = fgetl(dump);
    if (strncmpi(id,'ITEM: TIMESTEP',numel('ITEM: TIMESTEP')))
            if t == 0
                timestep(i) = str2num(fgetl(dump));
                t=1;
            end
    else
     if (strcmpi(id,'ITEM: NUMBER OF ATOMS',numel('ITEM: NUMBER OF ATOMS')))
            Natoms = str2num(fgetl(dump));
     else
      if (strcmpi(id,'ITEM: BOX BOUNDS',numel('ITEM: BOX BOUNDS')))
            x_bound(1,:) = str2num(fgetl(dump));
            y_bound(1,:) = str2num(fgetl(dump));
            z_bound(1,:) = str2num(fgetl(dump));
      else
       if (strncmpi('ITEM: ATOMS',numel('ITEM: ATOMS')))
            atom_data = zeros(Natoms,ncol);%Allocate memory for atom data
            for j = 1 : 1: Natoms
                atom_data(j,:) = str2num(fgetl(dump));
            end
            done = 1;
            p = ftell(dump);
       end
      end 
     end
    end
end

% Getting only the last step
if last_status == 1
    % First get the position of the beginning of the last step in the
    % dumpfile
    while ~feof(dump)
        temp = fgetl(dump);
        if length(temp) == 14
            if strcmpi(temp,'ITEM: TIMESTEP')
                p = ftell(dump); % starting position of line next to the header
            end
        end
    end
    fclose(dump);
    dump = fopen(varargin{1},'r');
    fseek(dump,p,'bof');
    % getting Timestep
    timestep = str2num(fgetl(dump));
    
    while ~feof(dump)
        id = fgetl(dump);
     if (strcmpi(id,'ITEM: NUMBER OF ATOMS'))
            Natoms = str2num(fgetl(dump));
     else
      if (strcmpi(id,'ITEM: BOX BOUNDS'))
            x_bound(1,:) = str2num(fgetl(dump));
            y_bound(1,:) = str2num(fgetl(dump));
            z_bound(1,:) = str2num(fgetl(dump));
      else
       if (strcmpi(id(1:11),'ITEM: ATOMS'))
            atom_data = zeros(Natoms,ncol);%Allocate memory for atom data
            for j = 1 : 1: Natoms
                atom_data(j,:) = str2num(fgetl(dump));
            end
       end
      end 
     end
    end
  

        
end

%----------Outputs-------------

%OUTPUTS IN SAME VARIABLE STRUCTURE
varargout{1}.timestep = timestep;
varargout{1}.Natoms = Natoms;
varargout{1}.x_bound = x_bound;
varargout{1}.y_bound = y_bound;
varargout{1}.z_bound = z_bound;
varargout{1}.atom_data = atom_data;
varargout{1}.position = p; %gives postion of ITEM: TIMESTEP line
%------------------------------

fclose(dump);

