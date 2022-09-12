function gmacmd_sliding_fft(mooring,Environment)

% mooring is a structure similar to the one find in GMACMD
% Environement

%% make a time axis
mooring.timeaxis = mooring.begintime+mooring.increment/1440*(0:length(mooring.u)-1);
%% get Coriolis frequency 
speed_earth=2*pi/86400;

mooring.f=2*speed_earth*sind(-mooring.latitude); % -latitude because we are looking at rotary spectra. f is clockwise/omega negative in the northern hemisphere.
mooring.f=mooring.f*86400/2/pi;

%% get TPXO tidal output to remove the barotropic tides from the recorded velocities
%   TODO: function for TPXO
if(Environment.add_TPXO)    
    [mooring.utrope,~]=tmd_tide_pred([Environment.gmacmdroot,'GMACMD/TPXO/TMD/DATA/Model_tpxo7.2'],...
        mooring.timeaxis,mooring.latitude,mooring.longitude,'u');
    [mooring.vtrope,~]=tmd_tide_pred([Environment.gmacmdroot,'GMACMD/TPXO/TMD/DATA/Model_tpxo7.2'],...
        mooring.timeaxis,mooring.latitude,mooring.longitude,'v');
end

%% stratification

if(Environment.add_levitus)    
    mooring  = add_n2levitus(mooring,Environment);
end

%% correct for small gaps (<6 hours) and split the time series in blocks is needed
%  add the percentage of data filled
if(Environment.add_block)
mooring  = correct_timeserie(mooring,Environment.Lgap,Environment.Lwin); % check for plateau spikes and holes.Build blocks of more than 30 days
end

if isfield(mooring,'block')
    if(Environment.add_spectra)
        mooring  = sliding_spectra(mooring);
    end
    if(Environment.add_KEpeak)
        mooring  = get_KE_peak(mooring,Environment);
    end
    if(Environment.add_GM)
        mooring  = get_GM(mooring,Environment);
    end
    
    save(fullfile(Environment.savepath,[Environment.name '.mat']),'mooring')
    
    if Environment.dsp.flag==1
        Environment.nbblock=1;
        func_GMACMDmovie(mooring,Environment);
    end
else
    Environment.noblock_ind=Environment.noblock_ind+1;
    Environment.noblock_files{Environment.noblock_ind}=Environment.file;
end


end