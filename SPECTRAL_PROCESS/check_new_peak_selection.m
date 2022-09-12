close all
peak_enhance_spec(peak_enhance_spec==0)=NaN;
% Fc=.1;% 1/4 hour filter
% [b,a]=butter(5,Fc,'low');
smspec = peak_enhance_spec;
smspec(smspec==0) = nan;
smspec=fillmissing(smspec,'linear');

 smspec =10.^smoothdata(log10(smspec),'movmean',length(freq)/3);
% smspec =medfilt1(peak_enhance_spec,30);

ax(1)=subplot(211);
semilogy(freq,peak_enhance_spec,'b')
hold on
semilogy(freq,smspec,'r')

% threshold=nanmean(log10(peak_enhance_spec./smspec));
threshold=nanmean(peak_enhance_spec./smspec);
grid on
ax(2)=subplot(212);
% plot(freq,log10(peak_enhance_spec./smspec))
semilogy(freq,peak_enhance_spec./smspec)
hold on
semilogy(freq,freq.*0+threshold,'r')
linkaxes(ax,'x')

%%
idx=find(peak_enhance_spec./smspec>threshold);

splitidx=[find(diff(idx)>1) numel(idx)];
nb_peaks=length(splitidx);

%find the max FOCO of the peak
peak_freq=zeros(nb_peaks,1);
local_peaks=zeros(nb_peaks,1);
start_local_idx=1;
for i=1:nb_peaks
    local_idx=idx(start_local_idx:splitidx(i));
    local_foco=peak_enhance_spec(local_idx);
    peak_freq(i)=freq(local_idx(local_foco==max(local_foco)));
    local_peaks(i)=max(local_foco);
    start_local_idx=splitidx(i)+1;
end

hold on
loglog(peak_freq,local_peaks);

semilogy(ax(1),freq(idx),peak_enhance_spec(idx),'m')


