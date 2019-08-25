function lmp2cfg(varargin)
% Converts LAMMPS dump file to Extended CFG Format (No velocity) to be used
% with AtomEye (http://164.107.79.177/Archive/Graphics/A/)
% Input : 
%   Necessary (in order)
%       timestep,Natoms,x_bound,y_bound,z_bound,H,atom_data,mass,
%       cfg file name, dumpfile name
%   Optional
%       'vel'        --> 'yes' (if velocity data needs to be written)
%                    --> Default is 'no'
%       'atomtype'   --> {Atom symbol}
%       'autonumber' --> 'yes' (default) numbers cfgfiles consecutively based on timestep
%       'autolog'    --> 'on' (default) writes details of conversion to a log file 
%       'rotation'   --> To specify rotation (Transform in cfg file) say 'rotation' followed by
%                        number of rotations,the axes of rotation ('X' or 'Y' or 'Z')and the 
%                        angle in degrees. CAUTION : Do remember that rotations are
%                        non-commutative and hence the order of rotations must be the same as
%                        intended
%       'aux'        --> Auxiliary data must be a two column array with Quantity and Type in each
%                        column
%
%
% THE DATA MUST BE SCALED (0 to 1)IN ORDER FOR ATOMEYE TO READ PROPERLY
% Real coordinates x = s * H,  x, s are 1x3 row vectors
%
% Example
%       lmp2cfg_all(data.timestep,data.Natoms,data.xbound,data.ybound,...
%                   data.zbound,H,data.atom_data,mass,cfgfile,dumpfile,...
%                   'vel','no','aux',{'sxx' 'Lammps O/P'
%                                  'syy' 'Lammps O/P'
%                                  'szz' 'Lammps O/P'
%                                  'sxy' 'Lammps O/P'
%                                  'sxz' 'Lammps O/P'
%                                  'syz' 'Lammps O/P'},...
%                   'atomtype',{'Ni'
%                               'Al'},...
%                   'autonumber','on',...
%                   'autolog','on',...
%                   'rotation',2,'X',45,'Y',30);
%
% See also readdump_all, readdump_one, scandump
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.

%------------ Defaults
vel = 'no'; 
auxstatus = 0; 
atom_typestatus = 0;
rotstatus = 0;
autonumber_status = 1; % default is ON
autolog_status = 1;
%----------------------------------

%-----Required Input --------------
timestep    = varargin{1}
Natoms      = varargin{2};
x_bound     = varargin{3};
y_bound     = varargin{4};
z_bound     = varargin{5};
H           = varargin{6};
atom_data   = varargin{7};
mass        = varargin{8};  %arrange it in order of atom type , ie 1 - mass1 etc
filename    = varargin{9}; %Just give one file name

if length(varargin) < 10
    dumpfilename = ['dump.' filename];
else
    dumpfilename = varargin{10};
end



if length(varargin) > 10
    i=11;
    while i<length(varargin)
        id = varargin{i};
        switch id
            case 'vel' 
                vel = varargin{i+1};
                i=i+2;
            case 'aux' 
                auxiliary = varargin{i+1};
                auxstatus = 1;
                i=i+2;
            case 'atomtype'
                atom_type = varargin{i+1};
                atom_typestatus = 1;
                i=i+2;
            case 'rotation'
                nrot = varargin{i+1}; % number of rotations
                i = i+2;
                if nrot <=0 
                    error('Number of rotations must be a positive value');
                end
                rotstatus = 1;
                for j = 1 : 1 : nrot
                    rotaxes{j} = varargin{i};
                    i = i + 1;
                    rotangle(j) = varargin{i};
                    i = i + 1;
                end
            case 'autonumber'
                numberstate = varargin{i+1};
                i = i +2;
                if strcmpi(numberstate,'ON')
                    autonumber_status = 1;
                elseif strcmpi(numberstate,'OFF')
                    autonumber_status = 0;
                end
            case 'autolog'
                logstate = varargin{i+1};
                i = i + 2;
                if strcmpi(logstate,'ON')
                    autolog_status = 1;
                elseif strcmpi(logstate,'OFF')
                    autolog_status = 0;
                end
        end
    end
end

% Calculating Transformation matrix [T]
if rotstatus == 0
    T = eye(3); % Identity matrix
elseif rotstatus == 1
    T = eye(3);
    for j = 1 : 1 : nrot
        B = beta(rotangle(j)*pi/180,rotaxes{j});
        T = B*T;
    end
end
    
T = T'; % because Transform = [beta]transpose
        
        

%----------------------------------
nfiles = length(timestep);
%----Default Atom type names------------
if atom_typestatus == 0
    atom_type = {'Au'
                 'Ni'
                 'Zn'
                 'H'
                 'O'
                 'Cu'
                 'Al'
                 'Ag'
                 'C'
                 'Si'};
end
             
%--------Sorting atom types ---------
[s,id] = sort(atom_data(:,2,:));
%---------Writing CFG files---------------
for i = 1 : 1 : nfiles
    if autonumber_status == 1
        n = [filename '_' num2str(i) '.cfg'];
    elseif autonumber_status == 0
        n = [filename '.cfg'];
    end
    name{i}=n;
    fid = fopen(name{i},'w+');
    fprintf(fid,'Number of particles = %d\n',Natoms(i));
    fprintf(fid,'A = 1.0000000000 Angstrom \n');
    % Writing [H]
    fprintf(fid,'H0(1,1) = %f A \n',H(1,1));
    fprintf(fid,'H0(1,2) = %f A \n',H(1,2));
    fprintf(fid,'H0(1,3) = %f A \n',H(1,3));
    fprintf(fid,'H0(2,1) = %f A \n',H(2,1));
    fprintf(fid,'H0(2,2) = %f A \n',H(2,2));
    fprintf(fid,'H0(2,3) = %f A \n',H(2,3));
    fprintf(fid,'H0(3,1) = %f A \n',H(3,1));
    fprintf(fid,'H0(3,2) = %f A \n',H(3,2));
    fprintf(fid,'H0(3,3) = %f A \n',H(3,3));
    % Writing [T]
    fprintf(fid,'Transform(1,1) = %f \n',T(1,1));
    fprintf(fid,'Transform(1,2) = %f \n',T(1,2));
    fprintf(fid,'Transform(1,3) = %f \n',T(1,3));
    fprintf(fid,'Transform(2,1) = %f \n',T(2,1));
    fprintf(fid,'Transform(2,2) = %f \n',T(2,2));
    fprintf(fid,'Transform(2,3) = %f \n',T(2,3));
    fprintf(fid,'Transform(3,1) = %f \n',T(3,1));
    fprintf(fid,'Transform(3,2) = %f \n',T(3,2));
    fprintf(fid,'Transform(3,3) = %f \n',T(3,3));
    if strcmpi(vel,'no')
        fprintf(fid,'.NO_VELOCITY. \n');
    end
    fprintf(fid,'entry_count = %d \n',length(atom_data(1,:,i))-2);
    if auxstatus == 1
        for k = 1 : 1 : length(auxiliary(:,1))
            fprintf(fid,'auxiliary[%d] = %s [%s]\n',k-1,auxiliary{k,1},...
                                                     auxiliary{k,2});
        end
    end
    aid = atom_data(id(1,1,i),2,i);
%    aid = 1;
    atom_change = 1;
    for j = 1 : 1 : Natoms(i)
        if atom_change == 1
            fprintf(fid,'%f\n',mass(aid));
            fprintf(fid,'%s\n',atom_type{aid});
        end
        atom_change = 0;
        fprintf(fid,'%f\t',atom_data(id(j,1,i),3:length(atom_data(1,:,i)),i));
        fprintf(fid,'\n');
        if j ~= Natoms
            if atom_data(id(j,1,i),2,i) ~=  atom_data(id(j+1,1,i),2,i)
                atom_change = 1;
                aid = atom_data(id(j+1,1,i),2,i);
%                aid = aid+1;
            end
        end
    end
fclose(fid);
end 
if autolog_status == 1
    flog = fopen([filename '_lmp2cfg.log'],'w+');
    fprintf(flog,'----------------------------------------------------------\n');
    fprintf(flog,['LAMMPS DUMP to CFG file conversion :\t' datestr(now) '\n']);
    fprintf(flog,'----------------------------------------------------------\n');
    fprintf(flog,'LAMMPS Dump file : \t %s \n\n',dumpfilename);
    for i = 1 : 1 : nfiles
        fprintf(flog,'Timestep : %d --> \t\t %s \n',timestep(i),name{i});
    end
    fclose(flog);
end

%---------- Function to calculate beta
function b = beta(angle,axes)
        switch axes
            case 'X' % X axes
                b = [1  0           0
                     0  cos(angle)  sin(angle)
                     0  -sin(angle) cos(angle)];
            case 'Y' % Y axes
                b = [cos(angle)     0   sin(angle)
                     0              1   0
                     -sin(angle)    0   cos(angle)];
            case 'Z' % Z axes
                b = [cos(angle)  sin(angle) 0
                     -sin(angle) cos(angle) 0
                     0           0          1];
        end
end
% --------------------------------------------------

end % For main function
    

