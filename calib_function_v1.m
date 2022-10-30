% calib help
function [constante_cfit,gof,output] = ...
        calib_function_v1(sonde_alt, sonde_w,...
                          altitude, Pn2, Ph2o,...
                          jdi, jdf, bg_n2_std, bg_h2o_std)


sm = 11;

sm_h2o = smooth(Ph2o,sm);
sm_n2 = smooth(Pn2,sm);
% sm_n2_stdm = Pn2_stdm./sm.^0.5;
%sm_n2_stdm = Pn2_stdm;
% sm_n2_std = bg_n2_std
% num_std = bg_h2o_std
%sm_n2_stdm = bg_stdm_n2;
h2o = sm_h2o./sm_n2  ; % PC_H2o / N2

%mask_sm_n2 = (sm_n2./sm_n2_std) < 10;
mask_sm = ((sm_h2o./bg_h2o_std) < 10) | ((sm_n2./bg_n2_std) < 10) ;
mask_sm = mask_sm | altitude>10000 | altitude<200;
h2o(mask_sm) = NaN;

figure (333);
clf
hold on
plot(altitude,sm_n2)
plot(altitude(mask_sm),sm_n2(mask_sm),'o')
xlim([0,15000])
hold off

figure(334)
clf
hold on 
plot(altitude,sm_h2o)
plot(altitude(mask_sm),sm_h2o(mask_sm),'o')
xlim([0,15000])
hold off

figure (730)
clf; hold on; grid on; box on
plot(h2o,altitude/1000)
xlabel('Uncalibrated lidar mixing ratio (a.u.) ')
yl=ylabel('Altitude a.s.l. (km)');
set(yl,'units','normalized','position',[-0.05 0.5 0])
ylim([0 7])
t=text(.02, .02, '  b  ', 'units', 'normalize',...
       'Horiz','left','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);

% [sm_n2(1000:1500),altitude(1000:1500)./1000,sm_h2o(1000:1500),h2o(1000:1500)]
% bg_n2_std
% h2o = smooth(Ph2o,sm) ./ smooth(Pn2,sm)  ; % PC_H2o / N2
% % h2o = Ph2o ./ Pn2  ; % PC_H2o / N2


% mask_aux = ~isnan(h2o);
% h2o_lidar = interp1q(altitude(mask_aux),h2o(mask_aux), sonde_alt);

h2o_lidar = interp1q(altitude,h2o, sonde_alt);
% [h2o_lidar.*1000,sonde_alt./1000,sonde_w]
% h2o_lidar = interp1q(smooth(altitude,sm),smooth(h2o,sm), sonde_alt);

calib_H_inf = 200;        
% mask_sinal = h2o_lidar<0 
% sup = sonde_alt(mask_sinal)
sup = altitude(mask_sm & altitude>500);
sup(1)
sup(end)
if numel(sup)>0
    calib_H_sup = sup(1) - 50 ;
else
    calib_H_sup = 10000;
end
%calib_H_sup = 10000
% calib_H_sup = 3500;
disp(sprintf('calib_H_sup = %f',calib_H_sup))
mask = sonde_alt > calib_H_inf & sonde_alt < calib_H_sup ;

X = h2o_lidar(mask);
Y = sonde_w(mask);

% size(sonde_alt)
% size(mask)
% [sonde_alt/1000 mask h2o_lidar*1000 sonde_w ]

% mask = ~isnan(X) & ~isnan(Y);
%[X Y mask]
%constante_cfit = fit(X(mask),Y(mask),'poly1')
%myfunc = fittype( @(p1,x) p1*x )
%[constante_cfit,gof,output] = fit(X,Y,myfunc)
[constante_cfit,gof,output] = fit(X,Y,'poly1')

% constante_calib.p2 = cfun_calib.p2;
% constante_calib.p1 = cfun_calib;



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
% plot(constante_cfit)
%plot(a)
plot(constante_cfit,'predfunc')
legend('Data','Fit','Location','northeast');
%legend1 = legend(gca,'show');
%set(legend1,'Location','southeast');
xlabel('Uncalibrated lidar mixing ratio (a.u.)')
yl=ylabel('Sounding mixing ratio (g/kg)');
set(yl,'units','normalized','position',[-0.07 0.5 0])
xticks(0:0.004:0.02)
yticks(0:4:20)
ylim([0 20])
t=text(.98, .02, '  a  ', 'units', 'normalize',...
       'Horiz','right','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);


figure(1541)
clf
hold on
plot(X,Y-constante_cfit(X),'ob')
hold off
% legend('Razão de Mistura de H2O sonda')
title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
xlabel('Razão de Mistura Lidar Não Calibrado')
ylabel('Residuo [g/kg]')
grid on
 %%
 
disp(['p1 = ' num2str(constante_cfit.p1) ])
%disp(['p2 = ' num2str(constante_cfit.p2) ])
figure (1542)
clf; hold on; grid on; box on
%plot(h2o.*constante_cfit.p1, altitude/1000,'b')
plot(h2o.*constante_cfit.p1 + constante_cfit.p2, altitude/1000,'b')
plot(sonde_w,sonde_alt./1000,'o-r')
plot(h2o.*829, altitude/1000,'-k')
plot(h2o.*(829+82), altitude/1000,':','color','#909090')
plot(h2o.*(829-82), altitude/1000,':','color','#909090')
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

