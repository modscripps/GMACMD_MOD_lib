function [E0,fom]=find_Snorm(Conti,GM)    
    % compute the GM spectrum for all the E0 "guess" and compute the least square 
    
    %dof = nb taper
    dof=3;
    sig_lnS=5/4*dof^(-7/9);
     
    Conti=Conti(:);GM=GM(:);
    GM(GM==Inf)=nan;
    inan=(~isnan(GM) & ~isnan(Conti));
    
    E0=nanmedian(GM(inan).\Conti(inan));
    
    fom=log10(Conti(inan)./(E0.*GM(inan)));
    fom=nanvar(fom)./sig_lnS;
    
        
        
    
    

