function varargout = readrdf(varargin)
% Function to read Radial Distribution Funtion output from LAMMPS
% Input
% 'bin'       --> number of bins in rdf histogram
% 'runtime'   --> Run length of each of the run commands
% 'step'      --> rdf Dump step for each of the run commands
% 'ncol'      --> number of columns in the file
% 'final'     --> 'YES' indicates that only the average values will be extracted
% If only the averages are needed, you don't need to specify 'runtime',
% 'step' and 'ncol'
%
% Output is in the form of a structure with following variables
% 
% if 'final' was given as 'NO'
% .rdf_step_data      --> Matrix with RDF for each of the dumped steps
% .rdf_ave_data       --> Matrix with average RDF for each RUN
% .rdf_ave_data       --> Matrix with average RDF for each RUN
% 
% if 'final' was given as 'YES'
% .rdf_ave_data       --> Matrix with average RDF for each RUN
% 
% Example
%       rdf = readrdf('inputfile','bin',100,'runtime',[30000;20000],...
%                     'step',[100;100],'ncol',3,'final','no')
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.


rdf_name = varargin{1}; % RDF File name

% Setting all status checks to zero
bin_status = 0;
runtime_status = 0;
step_status = 0;
ncol_status = 0;
final_status = 'no';



for id = 2 : 1 : length(varargin)
    if strcmpi(varargin{id},'bin')
        rdf_bin = varargin{id+1}; % RDF bins
        bin_status = 1;
    elseif strcmpi(varargin{id},'runtime')
        rdf_runtime = varargin{id+1}; % Runtimes
        runtime_status = 1;
    elseif strcmpi(varargin{id},'step')
        rdf_step = varargin{id+1}; % Runtimes
        step_status = 1;
    elseif strcmpi(varargin{id},'ncol')
        ncol = varargin{id+1}; % Runtimes
        ncol_status = 1;
    elseif strcmpi(varargin{id},'final')
        final_status = varargin{id+1};
    end
end

if ncol_status == 0
    ncol = 3;
end

% Check for errors in input arguments
if bin_status == 0
    error('No bin specified');
elseif step_status == 1 && runtime_status == 0
    error('Step size specified without Runtime');
elseif step_status == 0 && runtime_status == 1
    error('Runtime specified without Step size');
end
if step_status == 1 && runtime_status == 1
    if length(rdf_runtime) ~= length(rdf_step)
        error('Runtime and Step size do not match');
    end
end

% Preallocating memory if runtime and step size are provided
if step_status == 1 && runtime_status == 1 && strcmpi(final_status,'no')
    total_steps = round(sum(rdf_runtime./rdf_step));
    rdf_step_data = zeros(rdf_bin,ncol,total_steps);
    rdf_ave_data = zeros(rdf_bin,ncol,length(rdf_runtime));
elseif strcmpi(final_status,'yes') && step_status == 1 && runtime_status == 1
    rdf_ave_data = zeros(rdf_bin,ncol,length(rdf_runtime));
end


try
    rdf_file = fopen(rdf_name,'r');
catch
    error('RDF file not found!');
end



if strcmpi(final_status,'yes')
    run_id = 1; % Run id..
    while feof(rdf_file) ~= 1
        rdf_data = fgetl(rdf_file);
        if strcmpi(rdf_data(1:11),'RUN AVERAGE')
            fgetl(rdf_file); % to skip the title
            for id = 1 : 1 : rdf_bin
                rdf_ave_data(id,:,run_id) = str2num(fgetl(rdf_file));
            end
            run_id = run_id + 1;
        end
        
    end
else
    run_id = 1;
    id = 1;
     while feof(rdf_file) ~= 1
        rdf_data = fgetl(rdf_file);
        if strcmpi(rdf_data(1:8),'TIMESTEP')
            timestep(id,1) = str2num(rdf_data(10:length(rdf_data)));
            fgetl(rdf_file); % to skip the title
            for j = 1 : 1 : rdf_bin
                rdf_step_data(j,:,id) = str2num(fgetl(rdf_file));
            end
            id = id+1;
        elseif strcmpi(rdf_data(1:11),'RUN AVERAGE')
            fgetl(rdf_file); % to skip the title
            for j = 1 : 1 : rdf_bin
                rdf_ave_data(j,:,run_id) = str2num(fgetl(rdf_file));
            end
            run_id = run_id + 1;
        end
     end
end
 
fclose(rdf_file);

if strcmpi(final_status,'no')
    out_data.timestep = timestep;
    out_data.rdf_step_data = rdf_step_data;
    out_data.rdf_ave_data = rdf_ave_data;
else
    out_data.rdf_ave_data = rdf_ave_data;
end
varargout{1} = out_data;
