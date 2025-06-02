function varargout = readEAM(varargin)
% Function to read EAM potential files
% Input 
%       1: EAM potential file name
%       2: file type --> 'FUNCFL' for single element file
%                    --> 'SETFL'  for multiple element file
% Output : Structure with members       
% .embed    : Embedding function
% .pair     : Pair Potential
% .elecden  : Electron Density
% .nrho     : Number of points for electron density
% .drho     : increment of r for electron density
% .nr       : Number of points for pair potential
% .dr       : increment of r for pair potential
% .rcut     : cut-off distance
% All output is in exactly the same units as in the EAM file. For
% multielement SETFL files the members are multidimensional arrays with
% each element data stored in  (:,:,ELEM)
%
% Example
%       eam = readEAM('cuu3.eam','funcfl');
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.

if length(varargin) < 2
    error('Too few input parameters : Need Filename and filetype');
elseif length(varargin) > 2
    error('Too many input parameters');
end

filename = varargin{1};
type = varargin{2};

try
    fid = fopen(filename,'r');
catch
    error('EAM file not found!');
end

if strcmpi(type,'FUNCFL')
    t = 1;
elseif strcmpi(type,'SETFL')
    t = 2;
else
    error(['Unknown file type : "' type '"']);
end

switch t
    case 1 %FUNCFL
        % Read line 1 : comment line
        a = fgetl(fid); 
        
        % read Line 2
        rem = fgetl(fid);           % ielem amass blat lat
        [ielem,rem] = strtok(rem);
        [amass,rem] = strtok(rem);
        [blat,rem] = strtok(rem);
        [lat,rem] = strtok(rem);
        ielem = str2num(ielem);     % Atomic Number
        amass = str2num(amass);     % Atomic Mass
        blat = str2num(blat);       % Lattice Constant
        
        % Read Line 3
        a = str2num(fgetl(fid));
        nrho = a(1);
        drho = a(2);
        nr = a(3);
        dr = a(4);
        rcut = a(5);
        
        % Reading embedding function
        for i = 1 : 1 : nrho/5
            embed(i,:) = str2num(fgetl(fid));
        end
        
        % Reading pair potential
        for i = 1 : 1 : nr/5
            pair(i,:) = str2num(fgetl(fid));
        end
        
        % Reading electron density function
        for i = 1 : 1 : nr/5
            elecden(i,:) = str2num(fgetl(fid));
        end
        
        % Output
        out.embed       = embed;
        out.pair        = pair;
        out.elecden     = elecden;
        out.nrho        = nrho;
        out.drho        = drho;
        out.nr          = nr;
        out.dr          = dr;
        out.rcut        = rcut;
        
        varargout{1} = out;
% --------------------------------------------------------------------

    case 2      % SETFL 
        
        % Read lines 1 - 3 : comment lines
        a = fgetl(fid);
        a = fgetl(fid);
        a = fgetl(fid);
        
        % Atom types
        ntypes = str2num(fgetl(fid));
        
        % Read Global information 
        a = str2num(fgetl(fid));
        nrho      = a(1);
        drho      = a(2);
        nr        = a(3);
        dr        = a(4);
        rcut   = a(5);
        
        % Read element specific Data
        % Embedding function and Electron Density
        for elem = 1 : 1 : ntypes
            rem = fgetl(fid);           % ielem amass blat lat
            [ielem1,rem]    = strtok(rem);
            [amass1,rem]    = strtok(rem);
            [blat1,rem]     = strtok(rem);
            [lat1,rem]      = strtok(rem);
            ielem(elem)     = str2num(ielem1);     % Atomic Number
            amass(elem)     = str2num(amass1);     % Atomic Mass
            blat(elem)      = str2num(blat1);      % Lattice Constant
            lat(elem,:)       = lat1;                % Lattice type
            
            % Reading embedding function
            for i = 1 : 1 : nrho/5
                embed(i,:,elem) = str2num(fgetl(fid));
            end
            
            % Reading electron density function
            for i = 1 : 1 : nr/5
                elecden(i,:,elem) = str2num(fgetl(fid));
            end
        end
        
        % Pair Potentials
        n_pair = ntypes * (ntypes + 1) / 2;
        for np = 1 : 1 : n_pair
            for i = 1 : 1 : nr/5
                pair(i,:,np) = str2num(fgetl(fid));
            end
        end
        
        % Output
        out.embed       = embed;
        out.elecden     = elecden;
        out.pair        = pair;
        out.nrho        = nrho;
        out.drho        = drho;
        out.nr          = nr;
        out.dr          = dr;
        out.rcut        = rcut;
        
        varargout{1} = out;
        
end

