function mooring = get_NIshape(mooring)

% loop over the blocks of individual time serie
for ii=1:length(mooring.block)

    % increment block counter
    T=length(mooring.block{ii}.timespec);
    mooring.block{ii}.NIslope=nan.*(1:T);
    tempo_NIslope=zeros(1,T);
%     count=count+1;
    % find the f peak (if any)
    idx_ni=cellfun(@(x) strcmp(x,'f'),{mooring.block{ii}.KESource.name}); % position of f in the selected KE peak
    idx_peak=[mooring.block{ii}.KESource(~idx_ni).idx_range]; % position (idx) of the NI peak in the spectrum.
    
    % create data no peak
    ip=find(mooring.block{ii}.freq>0); % positive freq
    im=find(mooring.block{ii}.freq<0); % neg freq
    data=mooring.block{ii}.spec; % time avg of the time evolving spectrum
    % temporary spectrum
    tempo=data;
    % nan all peaks that are not NI
    tempo(:,idx_peak)=nan;
    % fold the spectra
    data1d_nopeak=tempo(:,ip)+fliplr(tempo(:,im));
    freq=mooring.block{ii}.freq_pos;
    fi=abs(mooring.f);
    
    loc_freq_stretch=(freq-fi)./fi;
    data1d_nopeak(:,loc_freq_stretch<0)=NaN;
    idx_nonnan=find(~isnan(data1d_nopeak(1,:)));
    
    log10data=log10(data1d_nopeak(:,idx_nonnan(1):idx_nonnan(end)));
    
    log10data= ...
        fillmissing(log10data,'linear',2);
    data1d_nopeak(:,idx_nonnan(1):idx_nonnan(end))=10.^log10data;
    
    loc_freq_stretch(loc_freq_stretch<=0)=nan;
    idxNI=find(log10(loc_freq_stretch)>-1.5 & log10(loc_freq_stretch)<-0.5);
    % only use blocks where I have non nan idxNi > length(idxNI)
    fprintf('length idxNI %i\r\n',length(idxNI))
    fprintf('good idxNI %i\r\n',sum(~isnan(data1d_nopeak(1,idxNI))))
    
    if length(idxNI)>5
        if sum(~isnan(data1d_nopeak(1,idxNI)))>=.5*length(idxNI)
            %enough data near NI, do the job ...
            for t=1:T
                NIslope_poly=polyfit(log10(freq(idxNI)), ...
                    log10(data1d_nopeak(t,idxNI)),...
                    1);
                tempo_NIslope(t)=NIslope_poly(1);
                % store in mooring
                mooring.block{ii}.NIslope(t)=NIslope_poly(1);
            end
        end
    end
end
end % end for idx looppc
