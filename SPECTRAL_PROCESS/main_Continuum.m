% main script to Continuum processing

% get_TPXO_barotropic_tide;
% check_extract;
% create_wkb;
% create_spectra;
% compute_peak;
% compute_GM;

% written by aleboyer@ucsd.edu
%% FOR LETICIA if you want you can generate your own metadta data file using the GMACMD library 
%  I created the folloing one mooring_metadata_3000
Meta= load('mooring_metadata_3000.mat');

%% Set paths 

Environment.gmacmdroot='/Volumes/GoogleDrive/My Drive/DATA/';
Environment.N0=5.2e-3;
Environment.dsp.flag=0;
% The climato is what I used to do the WKB stretching you might want to
% change to Jeffrey's vertival structures.
load(fullfile(Environment.gmacmdroot, 'GMACMD/Levitus/levitus_climatology.mat'),'climatology');
Environment.climatology=climatology;
% Environment.Source=dir(fullfile(Environment.gmacmdroot,'GMACMDlocal/*.mat'));
Environment.Source=Meta.Source_mat;

%%
addpath(genpath('/Volumes/DataDrive/GMACMD/TPXO'));
addpath('/Volumes/GoogleDrive/My Drive/TOOLBOXES/GM/GarrettMunkMatlab/');

%% Environment parameters
Environment.dof=5;
Environment.pds="MySineSpec";  % options: MySineSpec or pwelch
Environment.Lwin=30;           % 30 days FFT windows "length" in days. TODO find a better name 
Environment.Lgap=6;            % 6 hours. Minimal length gap inside the timeseries aloowed for linear interpolation
Environment.add_TPXO=0;        % flag to add TPXO
Environment.add_levitus=1;     % flag to add levitus yearly climatology to the mooring structure
Environment.add_block=1;       % flag to correct the time series abd add block to the mooring structure
Environment.add_spectra=1;     % flag to compute the sliding FFT over the timeseries with a window of Lwin length 
Environment.add_KEpeak=1;      % flag to identifie the main KE peak.
Environment.add_GM=1;          % flag to add the GM spectrum
Environment.add_NIshape=1;

%LETICIA: I created a folder inside my perosnal GMACMD repo. 
% you do not have to do this here but 
% you want to have a folder of "processed data" somewhere 
Environment.savepath=fullfile(Environment.gmacmdroot,'GMACMD/testleticia');

%%
Environment.noblock_ind=0;
% for m=1:length(Environment.Source)
%  10/21/2021 M=1200
% for m=1:length(Environment.Source)
for m=1:1

    
    Environment.file=fullfile(Environment.gmacmdroot,Environment.Source{1}(1:end-4));
    mooring=load(Environment.file);
    % ALB I am not working near the equator
    if abs(mooring.latitude)>1
        namesplit=strsplit(Environment.Source{1},'/');
        Environment.name=namesplit{end}(1:end-4);
        Environment.moviename=sprintf('%sGMACMD/MOVIE/%s_f%i.avi',Environment.gmacmdroot,Environment.name,m);
        Environment.mooringnumber=m;
        
        if ~isfolder(Environment.savepath)
            mkdir(Environment.savepath)
        end
        Epsilon_KE_parmetrization(mooring,Environment)
        fprintf('mooring %i \n',Environment.mooringnumber);
        close all
    end
    
end

