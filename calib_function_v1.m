% calib help
function [fitAB, fitA] = ... 
        calib_function_v1(sonde_alt, sonde_w,...
                          altitude, Pn2, Ph2o,...
                          jdi, jdf, snr_n2, snr_h2o)
    % bug: using bg_std is incorrect
    %                          jdi, jdf, bg_n2_std, bg_h2o_std)



%% calculate h2o as ratio of Ph2o and Pn2

% smooth parameter
sm = 11;

sm_h2o = smooth(Ph2o,sm);
sm_n2 = smooth(Pn2,sm);
sm_h2o = sm_h2o./sm_n2  ; % PC_H2o / Glue_N2

% without smooth, just for plotting
h2o=Ph2o./Pn2;

%% Region for calibration

% initial range
calib_H_inf = 1600;
calib_H_sup = 10000;

%mask_sm = ((sm_h2o./bg_h2o_std) < 10) | ((sm_n2./bg_n2_std) < 10) ;
mask_sm = (snr_n2 < 3); % | (snr_h2o < 10);
mask_sm = mask_sm | altitude>calib_H_sup | altitude<calib_H_inf;

% Plot SNR and region used for calibration
figure (123);
clf; hold on; grid on; box on
plot(snr_h2o,altitude/1000,'b'); 
plot(snr_n2 ,altitude/1000,'r'); 
plot(snr_h2o(mask_sm),altitude(mask_sm)/1000,'xb')
plot(snr_n2 (mask_sm),altitude(mask_sm)/1000,'xr')
legend('h2o pc','n2 glue')
xlabel('Signal / Noise')
set(gca,'xscale','log')
yl=ylabel('Altitude a.s.l. (km)');
set(yl,'units','normalized','position',[-0.05 0.5 0])

% Plot the uncalibrated lidar signal
figure (730)
clf; hold on; grid on; box on
plot(h2o,altitude/1000)
plot(sm_h2o,altitude/1000)
xlabel('Uncalibrated lidar mixing ratio (a.u.) ')
yl=ylabel('Altitude a.s.l. (km)');
set(yl,'units','normalized','position',[-0.05 0.5 0])
xlim([0 0.020])
ylim([0 7])
t=text(.02, .02, '  b  ', 'units', 'normalize',...
       'Horiz','left','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);


%% Apply mask to crop out the calibration region

%sm_h2o(mask_sm) = NaN;
%h2o(mask_sm) = NaN; 

% update range to crop sonde data
altlim = altitude(~mask_sm);
calib_H_inf = altlim(1);
calib_H_sup = altlim(end);

disp(sprintf('calib_H = [%f, %f]',calib_H_inf, calib_H_sup))

% interpolate to sonde altitudes
h2o_lidar = interp1q(altitude,sm_h2o, sonde_alt);

% mask sonde data (to remove the NaN above)
mask = (sonde_alt > calib_H_inf) & (sonde_alt < calib_H_sup) ;

X = h2o_lidar(mask);
Y = sonde_w(mask);

% fit with A*X
myfunc = fittype( @(p1,x) p1*x );
[constante_cfit,gof,output] = fit(X,Y,myfunc, 'startpoint', [750.]);
fitA.cfit = constante_cfit;
fitA.gof  = gof;
fitA.out  = output; 
fitA.cfit

% fit with A*X + B
[constante_cfit,gof,output] = fit(X,Y,'poly1');
fitAB.cfit = constante_cfit;
fitAB.gof  = gof;
fitAB.out  = output; 
fitAB.cfit

clear constante_cfit gof output

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

figure(1538)
clf; hold on; grid on; box on
plot(Pn2,altitude/1000,'g')
plot(Ph2o,altitude/1000,'b')
legend('Nitrogen','Water vapor')
set(gca,'xscale','log','xminorgrid','off')
xlabel('Signal (a.u.)')
yl=ylabel('Altitude a.s.l. (km)');
set(yl,'units','normalized','position',[-0.05 0.5 0])
ylim([0 7])
t=text(.02, .02, '  a  ', 'units', 'normalize',...
       'Horiz','left','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);

figure(1539)
clf
hold on
plot(sonde_w,sonde_alt./1000,'*-')
plot(sonde_w(mask),sonde_alt(mask)./1000,'or')
hold off
title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
legend('Razao de Mistura de H2O sonda','Pontos utilizados na calibração')
xlabel('Razao de Mistura [g/kg]')
ylabel('Altitude asl (km)')
grid on


figure(1540)
clf; hold on; grid on; box on
plot(X,Y,'*b')

plot(fitAB.cfit,'predfunc')
fitresultAB = evalc('fitAB.cfit');
text(0.45,0.2, fitresultAB(21:end-1),'units','normalized', ...
     'Interpreter','none','backgroundcolor','white','edgecolor','black')

plot(fitA.cfit,'g') % 'predfunc')
fitresultA = evalc('fitA.cfit');
text(0.11,0.8, fitresultA(23:end-1),'units','normalized', ...
     'Interpreter','none','backgroundcolor','white','edgecolor','black')

legend('Data','Fit A*X+B','Fit A*X','Location','northeast');
%legend1 = legend(gca,'show');
%set(legend1,'Location','southeast');
xlabel('Uncalibrated lidar mixing ratio (a.u.)')
yl=ylabel('Sounding mixing ratio (g/kg)');
set(yl,'units','normalized','position',[-0.07 0.5 0])
xticks(0:0.004:0.02)
yticks(0:4:20)
xlim([0 0.016])
ylim([0 16])
t=text(.98, .02, '  a  ', 'units', 'normalize',...
       'Horiz','right','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);


figure(1541)
clf
hold on
plot(X,Y-fitAB.cfit(X),'ob')
hold off
% legend('Razão de Mistura de H2O sonda')
title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
xlabel('Razão de Mistura Lidar Não Calibrado')
ylabel('Residuo [g/kg]')
grid on


figure (1542)
clf; hold on; grid on; box on
plot(h2o.*fitAB.cfit.p1 + fitAB.cfit.p2, altitude/1000,'b')
%plot(sm_h2o.*fitAB.cfit.p1 + fitAB.cfit.p2, altitude/1000,'b')
plot(sonde_w,sonde_alt./1000,'o-r')
%plot(h2o.*829, altitude/1000,'-k')
%plot(h2o.*(829+82), altitude/1000,':','color','#909090')
%plot(h2o.*(829-82), altitude/1000,':','color','#909090')

plot(h2o_lidar.*fitAB.cfit.p1 + fitAB.cfit.p2, sonde_alt./1000,'sm')
plot(h2o_lidar(mask).*fitAB.cfit.p1 + fitAB.cfit.p2, sonde_alt(mask)./1000,'sm','MarkerFaceColor','m')
plot(sonde_w(mask),sonde_alt(mask)./1000,'or','MarkerFaceColor','r')

legend('Calibrated lidar profile','Reference sounding','Best estimate')
xlabel('Mixing ratio (g/kg)')
yl=ylabel('Altitude a.s.l. (km)');
set(yl,'units','normalized','position',[-0.07 0.5 0])
ylim([0 7])
xlim([0 20])
xticks(0:4:20)
t=text(.98, .02, '  b  ', 'units', 'normalize',...
       'Horiz','right','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);

end

