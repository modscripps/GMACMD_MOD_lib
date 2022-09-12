function mooring = get_KE_peak(mooring,Environment)

%Environment.dsp is a structure 
% dsp.flag ==1 plot and save a plot showing the different KE peaks defined
% for this data set. The peak can change between blocks and moorings

% trick for conditional cellfun
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    
clear TC
% define tidal component TC
TC{1}.freq  =1;TC{1}.name='D1';TC{2}.freq=2;TC{2}.name='D2';
TC{3}.freq  =3;TC{3}.name='D3';TC{4}.freq=4;TC{4}.name='D4';
TC{5}.freq  =5;TC{5}.name='D5';TC{6}.freq=6;TC{6}.name='D6';

TC{7}.freq  =-1;TC{7}.name  ='-D1';TC{8}.freq  =-2;TC{8}.name='-D2';
TC{9}.freq  =-3;TC{9}.name  ='-D3';TC{10}.freq =-4;TC{10}.name='-D4';
TC{11}.freq =-5;TC{11}.name ='-D5';TC{12}.freq =-6;TC{12}.name='-D6';


% number of tidal component
NB_TC=length(TC);

% get inertial frequency
it=NB_TC;
%02/12/2019 ALB change abs(mooring.f+TC{t}.freq) to mooring.f+TC{t}.freq
for t=1:NB_TC
    it=it+1;TC{it}.freq=mooring.f+TC{t}.freq;TC{it}.name=['f' TC{t}.name];
end
TC{it+1}.freq=mooring.f;TC{it+1}.name='f';TC{it+2}.freq=2*mooring.f;
TC{it+2}.name='2f';NB_TC=it+2;

for ii=1:length(mooring.block)
    freq= mooring.block{ii}.freq;
    data= mooring.block{ii}.spec;
    N_freq=numel(freq);
    
    % get mean and std to determine real peaks for spectra
    peak_enhance_spec=mooring.block{ii}.meanspec; 

    %get the peaks
    [pks,locs,frange,foco_peak,idx_range]=find_peaks_alb(peak_enhance_spec,freq);
    nb_peaks=numel(pks);

    % if any of the theoretical tidal componant are included in the
    % peak definition
    edge1=cellfun(@(x) x(1),frange);
    edge2=cellfun(@(x) x(end),frange);

    peak_width=edge2-edge1;
    idx_TC=zeros(nb_peaks,1);
    TC_freq=cellfun(@(x) x.freq,TC);
    for i=1:nb_peaks
        local_mask=sign(TC_freq)==sign(locs(i));
        wh_freq=abs(TC_freq-locs(i));
        wh_freq(~local_mask)=NaN;
        idx_freq=find(wh_freq==min(wh_freq) & ...
                      peak_width(i)>min(wh_freq));
        switch length(idx_freq)
            case 0
                %no find
            case 1    
                idx_TC(i)=idx_freq;
            case 2 
                idx_TC(i)=idx_freq(1);
        end
    end
    
    good_peak=find(idx_TC>0);
    nb_peaks=length(good_peak);
    pks=pks(good_peak);
    locs=locs(good_peak);
    frange=frange(good_peak);
    foco_peak=foco_peak(good_peak);
    idx_range=idx_range(good_peak);
    
    edge1=cellfun(@(x) x(1),frange);
    edge2=cellfun(@(x) x(end),frange);

    clear KESource;
    KESource=[TC{idx_TC(good_peak)}];
    for i=1:nb_peaks
        KESource(i).frange=frange{i};
        KESource(i).foco=foco_peak{i};
        KESource(i).idx_range=idx_range{i};
        KESource(i).locs=locs(i);
    end
    
        
        % force the closest peak to f to be the f peak
    if  any(edge2>=mooring.f & edge1<=mooring.f)
        bool1=edge2>mooring.f & edge1<mooring.f;
        if (KESource(bool1).freq-locs(bool1))>(mooring.f-locs(bool1))
            KESource(bool1).name='f';
            lowf_continuum=edge1(find(edge1<mooring.f,1,'last'));
        else
            lowf_continuum=edge2(pks==max(pks));
        end
    else
        lowf_continuum=edge2(pks==max(pks));
    end
    lowf_continuum=lowf_continuum(1);
    
    % force the closest peak to 2f to be the 2f peak
    if  any(edge1<=2.*mooring.f & edge2>=2.*mooring.f)
        bool1=edge1<2.*mooring.f & edge2>2.*mooring.f;
        if (sign(edge1(bool1)) == sign(edge2(bool1)))
            KESource(edge1<2.*mooring.f & edge2>2.*mooring.f).name='2f';
        end
    end

    % compute the KE as function of "time" one point every 24 hours
    % (check sliding spectra)
    freq_res=nanmean(diff(freq));
    for s=1:nb_peaks
        if size(data,1)==1 % case where there the total time serie is 30 day (~ one spectrum)
            KESource(s).KE=sum(data(:,KESource(s).idx_range)).* freq_res;
        else % general case
            KESource(s).KE=sum(data(:,KESource(s).idx_range),2).* freq_res;
        end
        KESource(s).fp=locs(s);
        KESource(s).edge1=edge1(s);
        KESource(s).edge2=edge2(s);        
    end
    if nb_peaks==0
        s=0;
    end
    KESource(s+1).name='Mesoscale-CW';
    KESource(s+1).KE=sum(data(:,ceil(N_freq/2)+[-1 -2]),2).* freq_res;
    KESource(s+1).fp=nanmean(freq(ceil(N_freq/2)+[-1 -2]));
    KESource(s+1).edge1=freq(ceil(N_freq/2)-1);
    KESource(s+1).edge2=freq(ceil(N_freq/2)-2);
    KESource(s+1).frange=freq(ceil(N_freq/2)+[-1 -2]);
    KESource(s+1).foco=data(:,ceil(N_freq/2)+[-1 -2]);
    KESource(s+1).idx_range=ceil(N_freq/2)+[-1 -2];
    
    KESource(s+2).name='Mesoscale-CCW';
    KESource(s+2).KE=sum(data(:,ceil(N_freq/2)+[1 2]),2).* freq_res;
    KESource(s+2).fp=nanmean(freq(ceil(N_freq/2)+[1 2]));
    KESource(s+2).edge2=freq(ceil(N_freq/2)+2);
    KESource(s+2).edge1=freq(ceil(N_freq/2)+1);
    KESource(s+2).frange=freq(ceil(N_freq/2)+[1 2]);
    KESource(s+2).foco=data(:,ceil(N_freq/2)+[1 2]);
    KESource(s+2).idx_range=ceil(N_freq/2)+[1 2];

    
    mooring.block{ii}.KESource=KESource;
    
    
    if Environment.dsp.flag==0
        % define colormap for plots and earth rotation speed
        figure
        peak_enhance_spec(peak_enhance_spec==0)=nan;
        cmap=colormap(jet(100));
        ax=axes;
        
        plot(ax,freq,peak_enhance_spec,'k','linewidth',2);
        hold(ax,'on')
        for t=1:NB_TC
            plot(ax,[TC{t}.freq,TC{t}.freq],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',.5)
        end
        fill(ax,[-lowf_continuum freq(freq>=-lowf_continuum & freq<=lowf_continuum) lowf_continuum],...
            [min(peak_enhance_spec) peak_enhance_spec(freq>=-lowf_continuum & freq<=lowf_continuum) ...
            min(peak_enhance_spec)],[1 1 1])
        fill(ax,[freq(1) freq(freq<=-lowf_continuum) -lowf_continuum],...
            [min(peak_enhance_spec) peak_enhance_spec(freq<=-lowf_continuum) ...
            min(peak_enhance_spec)],.2*[1 1 1])
        fill(ax,[lowf_continuum freq(freq>=lowf_continuum) freq(end)],...
            [min(peak_enhance_spec) peak_enhance_spec(freq>=lowf_continuum) ...
            min(peak_enhance_spec)],.2*[1 1 1])
        colorscale=linspace(0,max(abs(locs)),100);
        for s=1:length(locs)
            fill(ax,[edge1(s) freq(freq<=edge2(s) & freq>=edge1(s)) edge2(s)],...
                [min(peak_enhance_spec) peak_enhance_spec(freq<=edge2(s) & freq>=edge1(s)) ...
                min(peak_enhance_spec)],cmap(find(abs(locs(s))<=colorscale,1,'first'),:))
            text(ax,locs(s),pks(s), KESource(s).name,'backgroundcolor','y','fontsize',25)
        end
        plot(ax,mooring.f.*[1 1],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',2)
        hold(ax,'off')
        grid(ax,'on')
        set(ax,'xlim',[-10,10])
        set(ax,'fontsize',20)
        set(ax,'ylim',[0,1.1*max(pks)])
        xlabel('cpd','fontsize',20)
        ylabel('m^2 s^{-2} /cpd','fontsize',20)
        ax.YScale='log';
        fig=gcf;fig.PaperPosition=[0 0 50 15];
        title(sprintf('%s,%i-block%i',Environment.name,ii))
        print(fullfile(Environment.savepath,['Get_peak' Environment.name '_block' num2str(ii) '.png']),'-dpng2')
    end
    
    
end
       

