function mooring = sliding_spectra(mooring,Environment)

Lwin=Environment.Lwin;
K=3; % multi taper param 
frsm=1/Lwin;  %smooth 30 jours ?

minute_per_day=1440;

for ii=1:length(mooring.block)
    % get rotary velocities
    data=complex(mooring.block{ii}.uL0,mooring.block{ii}.vL0);
    % number of sample in a day because we averaging over a day.
    samplefreq=round(minute_per_day/mooring.increment); % mooring increment needs to be in min. It is the default units for GMACMD
    % number of sample in 30 days
    Tmax=ceil(Lwin./abs(mooring.f))*abs(mooring.f);
    if(isnan(Tmax))
        Tmax=Lwin;
    end
    point=floor(Tmax*samplefreq); %Length in sample of a Lwin day period
    % length in sample of the time serie minus 1 days
    L1=length(data)-point;
    data=data(:);
    % compute the first Lwin block spectrum
    [f,Pblock,Rblock]=MySineSpec(data(1:1+point),samplefreq,K);
    freq_res=nanmean(diff(f));%(point+1)/samplefreq;
    
    %initialize the spectra array in struct mooring
    mooring.block{ii}.spec=zeros(floor(L1/samplefreq),length(Pblock));
    mooring.block{ii}.timespec=zeros(floor(L1/samplefreq),1);
    mooring.block{ii}.freq=f;
    i1=find(f>0);
    i2=find(f<0);
    mooring.block{ii}.freq_pos=f(i1);

    count=0;
    %tspec=nanmean(mooring.block{ii}.time(1:1+point));
    tspec=0;
    if L1/point>1
        for iii=1:L1
            [~,P,Rot]=MySineSpec(data(iii:iii+point),samplefreq,K);
            localt=nanmean(mooring.block{ii}.time(iii:iii+point));
            % sum for daily average
            if mod(iii,round(samplefreq))>0 
                Rblock=Rblock+Rot;
                Pblock=Pblock+P;
                tspec=tspec+localt;
            else
                count=count+1;
                mooring.block{ii}.spec(count,:)=(Pblock+P)./samplefreq;  % clockwise frequency negative.
                mooring.block{ii}.timespec(count)=(tspec+localt)./samplefreq;
                mooring.block{ii}.Rot(count,:)=(Rblock+Rot)./samplefreq;
                Pblock=0;tspec=0;Rblock=0;
            end
            if count==floor(L1/samplefreq)
                break;
            end
        end
        % ALB smooth spectrum : average on a moving window I could use
        % smoothdata
        spec1=SmoothSpec(f(i1),mooring.block{ii}.spec(count,i1),frsm);
        spec2=SmoothSpec(-f(i2),mooring.block{ii}.spec(count,i2),frsm);
        mooring.block{ii}.Rot(count,:)= ...
                    SmoothSpec(f(i1),mooring.block{ii}.Rot(count,:),frsm);
        mooring.block{ii}.spec(count,i1)=spec1;
        mooring.block{ii}.spec(count,i2)=spec2;
    else % case if the time serie is exactly 30 days
        Pblock=Pblock./L1;
        Rblock=Rblock./L1;
        tspec=nanmean(mooring.block{ii}.time(1:1+point));
        tspec=tspec./L1;
        for iii=1:L1
            [~,P,Rot]=MySineSpec(data(iii:iii+point),samplefreq,K);
            localt=nanmean(mooring.block{ii}.time(iii:1+point));
            mooring.block{ii}.spec(1,:)=Pblock+P./L1;
            mooring.block{ii}.Rot(1,:)=Rblock+Rot./L1;
            tspec=tspec+localt./L1;
        end
        mooring.block{ii}.timespec=tspec;
        spec1=SmoothSpec(f(i1),mooring.block{ii}.spec(1,i1),frsm);
        spec2=SmoothSpec(-f(i2),mooring.block{ii}.spec(1,i2),frsm);
        mooring.block{ii}.spec(1,i1)=spec1;
        mooring.block{ii}.spec(1,i2)=spec2;
        mooring.block{ii}.Rot(1,:)= ...
                    SmoothSpec(f(i1),mooring.block{ii}.Rot(1,:),frsm);
    end
    [spl,spu]=SmoothSpecCLim(f,K,frsm);
    mean_spec=mean(mooring.block{ii}.spec,1);
    std_spec=std(mooring.block{ii}.spec,[],1);
    mean_Rot=mean(mooring.block{ii}.Rot,1);
    std_Rot=std(mooring.block{ii}.Rot,[],1);
    mooring.block{ii}.meanspec=mean_spec+std_spec;
    mooring.block{ii}.spl=spl;
    mooring.block{ii}.spu=spu;
    mooring.block{ii}.meanRot=mean_Rot+std_Rot;
    
    % KE inside spectra for 3<omega<10 cpd on the cyclonic side;
    mooring.block{ii}.KE_f3f10CW=nansum(mooring.block{ii}.spec(:,mooring.block{ii}.freq<-3 & mooring.block{ii}.freq>-10),2)*freq_res;    
    mooring.block{ii}.KE_f3f10CCW=nansum(mooring.block{ii}.spec(:,mooring.block{ii}.freq>3 & mooring.block{ii}.freq<10),2)*freq_res;    
end
