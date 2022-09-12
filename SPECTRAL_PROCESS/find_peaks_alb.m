function [pks,locs,frange,foco_peak,idx_range]=find_peaks_alb(spectrum,freq)

% finds peak in a spectrum
% I smooth the spectrum a LOT, usually using a window with a size of 1/3 of
% the total length. The size of the window is defined by nb_segment (nb_segment=3)
% The ratio spectrum./smooth_spectrum hight lights the peak. I am selecting
% the peak that above the mean (in log10 base) of the this ratio.


 spectrum(spectrum==0)=NaN;
% Fc=50/100;% 
% [b,a]=butter(3,Fc);
% sm_spectrum=filtfilt(b,a,spectrum);
% sm_spectrum=smoothdata(spectrum,'movmean',5);
% ddspec = interp1(freq(2:end-1),...
%     filtfilt(b,a,diff(spectrum,2)),freq);
% ddspec = filtfilt(b,a,diff(spectrum,2));
% dspec = interp1(freq(1:end-1)+nanmean(diff(freq)),...
%                  diff(sm_spectrum./max(sm_spectrum),1),freq);

% get the select peak that are higher than 1/10 of the std FOCO
freqp=freq(freq>0);
freqm=freq(freq<0);
spectp=log10(spectrum(freq>0));
spectm=log10(spectrum(freq<0));
pp=polyfit(freqp,spectp,5);
pm=polyfit(freqm,spectm,5);
trend_spectp=10.^polyval(pp,freqp);
trend_spectm=10.^polyval(pm,freqm);
% trend_spectp=smoothdata(spectp,'movmean',length(spectp)/10);
% trend_spectm=smoothdata(spectm,'movmean',length(spectm)/10);
trend_spectrum=spectrum.*nan;

trend_spectrum(freq>0)=trend_spectp;
trend_spectrum(freq<0)=trend_spectm;

sm_spectrum=smoothdata(spectrum,'movmean',5);
% trend_spectrum=smoothdata(spectrum,'movmean',length(spectrum)/30);
detrend_spectrum=log10(sm_spectrum./trend_spectrum);
detrend_spectrum(abs(freq)>=6)=nan;
peak_height=.5*nanstd(detrend_spectrum);

% peak_height=nanstd(sm_spectrum(:))./10;
% [~,locs,~,~]=findpeaks(sm_spectrum,freq,'MinPeakHeight',peak_height);
% [~,locs,~,~]=findpeaks(detrend_spectrum,freq,'MinPeakProminence',peak_height);
[~,locs,~,~]=findpeaks(detrend_spectrum,freq,'MinPeakProminence',peak_height);
locs=unique(locs);

% center the freq vecteur on the peak
center_freq=@(x) (freq-x);
fmfp=arrayfun(center_freq,locs,'un',0);
ind_fp=cellfun(@(x) find(x==0),fmfp,'un',0);
indp=cellfun(@(x) (x>0),fmfp,'un',0);
indm=cellfun(@(x) (x<0),fmfp,'un',0);

% take both side of the first derivative of the spectra
% to get the edges of the peak (~ 0 of the 1st derivative)
% WATCH OUT: we suppose that a peak is alway surrounded be smaller
% peak and thus the change of sign in the derivative can help us to
% define its shape (width)
% frequency of a peak must fit into edge 1 and edge 2
clear dspecp; clear dspecm; clear ind_edge1;
clear ind_edge2; clear edge1; clear edge2;
% dspecp=cellfun(@(x) abs(dspec(x)),indp,'un',0);
% dspecm=cellfun(@(x) abs(dspec(x)),indm,'un',0);
detspecp=cellfun(@(x) detrend_spectrum(x),indp,'un',0);
detspecm=cellfun(@(x) detrend_spectrum(x),indm,'un',0);

nb_peaks=numel(locs);
%%
for i=1:nb_peaks
%     idx1=find(freq==locs(i));
%     idx2=find(freq==locs(i+1));
%     slot_spectrum=spectrum(idx1:idx2);
%     idx_trough=idx1+find(slot_spectrum==min(slot_spectrum));
%     if i==nb_peaks-1
%         idxpeak=ind_fp{i+1};
%         ddspecp=diff(dspecp{i+1});
%         edge2{i+1}=freq(idxpeak+2*find(ddspecp<0,1,'first')+1);
%     end
%     if i==1
%         idxpeak=ind_fp{i};
%         ddspecm=diff(fliplr(dspecm{i}));
%         edge1{i}=freq(idxpeak-2*find(ddspecm<0,1,'first'));
%     end
%     edge2{i}=freq(idx_trough-1);
%     edge1{i+1}=freq(idx_trough+1);
% 
%         idxpeak=ind_fp{i};
%         ddspecp=diff(dspecp{i});
%         edge2{i}=freq(idxpeak+2*find(ddspecp<0,1,'first')+1);
%         edge1{i}=freq(idxpeak-2*find(ddspecm<0,1,'first'));
        idxpeak=ind_fp{i};
        if(find(detspecp{i}<=0,1,'first')>10)
            edge2{i}=freq(idxpeak+find(detspecp{i}==min(detspecp{i}(1:10)))+1);
        else
            edge2{i}=freq(idxpeak+find(detspecp{i}<=0,1,'first')+1);
        end
        
        if(find(fliplr(detspecm{i})<=0,1,'first')>10)
            tempo=fliplr(detspecm{i});
            edge1{i}=freq(idxpeak-find(tempo==min(tempo(1:10)))+1);
        else
            edge1{i}=freq(idxpeak-find(fliplr(detspecm{i})<=0,1,'first'));
        end
        if isempty(edge1{i})
            edge1{i}=freq(max(1,idxpeak-floor(5*nanmean(diff(freq)))));
        end
        if isempty(edge2{i})
            edge2{i}=freq(min(length(freq),floor(idxpeak+5*nanmean(diff(freq)))));
        end

end

% ind_edge2 = cellfun(@(x,y) find(x<abs(dspec(y)),1,'first'),dspecp,ind_fp,'un',0); % always 1 but with this I am sure to start at the right place
% edge2     = cellfun(@(x,y) freq(x+y),ind_fp,ind_edge2,'un',0);
% ind_edge1 = cellfun(@(x,y) find(fliplr(x)<abs(dspec(y)),1,'first'),dspecm,ind_fp,'un',0);% always 1 but with this I am sure to start at the right place
% edge1     = cellfun(@(x,y) freq(x-y),ind_fp,ind_edge1,'un',0);

frange=cell(nb_peaks,1);
idx_range=cell(nb_peaks,1);
foco_peak=cell(nb_peaks,1);
pks=cell(nb_peaks,1);
clear locs
locs=cell(nb_peaks,1);
for i=1:nb_peaks
    x= edge1{i};    
    y= edge2{i};
    if x==0
        x=1e-12;
    end
    if y==0
        y=-1e-12;
    end
    frange{i}    = freq(freq>x & freq<y & sign(freq)==sign(x));
    idx_range{i} = find(freq>x & freq<y & sign(freq)==sign(x));
    foco_peak{i} = spectrum(idx_range{i});
    pks{i}=max(foco_peak{i});
    locs{i}=frange{i}(foco_peak{i}==pks{i});
end

emptypeak = cellfun(@(x)length(x)<5,idx_range);
frange    = frange(~emptypeak);
idx_range = idx_range(~emptypeak);
foco_peak = foco_peak(~emptypeak);
pks       = cell2mat(pks(~emptypeak));
locs      = cell2mat(locs(~emptypeak));
%%
% close all
% semilogy(freq,spectrum,'b')
% hold on
% semilogy(freq,sm_spectrum,'r')
% semilogy(freq,trend_spectrum,'k')
% scatter(locs,pks,40,'k','filled')
% 
% cellfun(@(x,y) scatter(x,y,20,'d','filled'),frange,foco_peak,'un',0)

