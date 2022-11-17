function [Signal_glued, cfun_glue, Gmask, midpoint] = ...
        glue(AN_withBG, PC_withBG, ANlim, PClim, AN, PC, range, Rlim)
    
% find region where PC and AN are valid
ANmask = (AN_withBG > ANlim(1)) & (AN_withBG < ANlim(2)); 
PCmask = (PC_withBG > PClim(1)) & (PC_withBG < PClim(2)); 
Gmask = ANmask & PCmask & range > Rlim;
tmp = (1:length(AN)); 
if sum(Gmask) > 3
    tmp = tmp(Gmask); 
    midpoint = floor((tmp(1) + tmp(end))/2.);

    X_glue = AN(Gmask);
    Y_glue = PC(Gmask);
    cfun_glue = fit(X_glue,Y_glue,'poly1')
    
    Signal_glued = cfun_glue(AN);
    Signal_glued(midpoint:end) = PC(midpoint:end);
else
    disp('WARNING!!! NOT ENOUGH DATA TO GLUE')
    Signal_glued = PC*nan;
    cfun_glue = nan;
    midpoint = nan;
end

end

