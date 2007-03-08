function [varargout] = scandump(varargin)
% Function to scan LAMMPS dump file
% Input is dumpfile name
% Output is a structure with the following members
% .timestep --> column vector with all the timesteps
% .natoms   --> column vector with number of atoms in each timestep
% .position --> column vector with position (for input into readdump)
% .ncol     --> column vector with number of columns
% .boxbound --> 3 by 3 by N (number of time steps) 3D array with 
%               [xlo xhi
%                ylo yhi
%                zlo zhi];
% Example
%       dump = scandump('dump.LAMMPS'); 
%
% See also readdump_one, scandump
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

i = 1;
out.position(i,1) = 0;

while ~feof(dump)
    id = fgetl(dump);
    switch id
        case 'ITEM: TIMESTEP'
            out.timestep(i,1) = str2num(fgetl(dump));
        case 'ITEM: NUMBER OF ATOMS'
            out.Natoms(i,1) = str2num(fgetl(dump));
        case 'ITEM: BOX BOUNDS'
            out.boxbound(1,:,i) = str2num(fgetl(dump));
            out.boxbound(2,:,i) = str2num(fgetl(dump));
            out.boxbound(3,:,i) = str2num(fgetl(dump));
        case 'ITEM: ATOMS'
            t = str2num(fgetl(dump));
            out.ncol(i,1) = length(t);
            for j = 2 : 1 : out.Natoms(i,1)
                fgetl(dump);
            end
            out.position(i+1,1) = ftell(dump);
            i = i+1;
    end
end

% Output
varargout{1} = out;
            
