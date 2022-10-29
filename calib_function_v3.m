% calib help
function [constante_cfit] = calib_function_v3(sonde_alt,sonde_ro,altitude,Pn2,Ph2o,jdi,jdf,Pn2_stdm,bg_n2_std,bg_h2o_std)


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
mask_sm = ((sm_h2o./bg_h2o_std) < 10) | ((sm_n2./bg_n2_std) < 10)  ;

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

%[sm_n2(1000:1500),altitude(1000:1500)./1000,sm_h2o(1000:1500),h2o(1000:1500)]
%bg_n2_std
% h2o = smooth(Ph2o,sm) ./ smooth(Pn2,sm)  ; % PC_H2o / N2
% % h2o = Ph2o ./ Pn2  ; % PC_H2o / N2


% mask_aux = ~isnan(h2o);
% h2o_lidar = interp1q(altitude(mask_aux),h2o(mask_aux), sonde_alt);

%h2o_lidar = interp1q(altitude,h2o, sonde_alt);
ro = interplq(sonde_alt,sonde_ro,altitude)
%[h2o_lidar.*1000,sonde_alt./1000,sonde_w]
% h2o_lidar = interp1q(smooth(altitude,sm),smooth(h2o,sm), sonde_alt);

calib_H_inf = 200;        
% mask_sinal = h2o_lidar<0 
% sup = sonde_alt(mask_sinal)
sup = altitude(mask_sm & altitude>500);
sup(1)
sup(end)
if numel(sup)>0
    calib_H_sup = sup(1) - 50 
else
    calib_H_sup = 10000;
end
%calib_H_sup = 10000
% calib_H_sup = 3500;

mask = sonde_alt > calib_H_inf & sonde_alt < calib_H_sup ;

X = h2o_lidar(mask);
Y = sonde_w(mask);

size(sonde_alt)
size(mask)
[sonde_alt/1000 mask h2o_lidar*1000 sonde_w ]

% mask = ~isnan(X) & ~isnan(Y);
%[X Y mask]
%constante_cfit = fit(X(mask),Y(mask),'poly1')
constante_cfit = fit(X,Y,'poly1')
% constante_calib.p2 = cfun_calib.p2;
% constante_calib.p1 = cfun_calib;



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%
figure(1538)
clf
range=(1:numel(altitude))*7.5;
plot(Pn2.*range.*range,altitude/1000,'g')
hold on
plot(Ph2o,altitude/1000,'b')
hold off
legend('Nitrogen','Water vapor')
set(gca,'xscale','log')
grid on 
title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
xlabel('Signal (a.u.)')
ylabel('Altitude a.s.l. (km)')
ylim([0 7])

%%
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
clf
hold on
plot(X,Y,'*b')
% plot(constante_cfit)
plot(constante_cfit,'predfunc')
hold off
% legend('Razão de Mistura de H2O sonda')
%title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
title('Ajuste linear que fornece a constante de calibração')
xlabel('Razão de Mistura Lidar Não Calibrado')
ylabel('Razão de Mistura Sonda [g/kg]')
grid on

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

figure (1542)
clf
hold on
plot(h2o.*constante_cfit.p1, altitude/1000,'b')
plot(h2o.*constante_cfit.p1 + constante_cfit.p2, altitude/1000,'r')
plot(sonde_w,sonde_alt./1000,'*-g')
legend('Razão de mistura do lidar sem somar o coef lin (g/kg)','Razão de mistra do lidar calibrado (g/kg)','Razão de mistura da radiossonda (g/kg)')
hold off
%ylim([0 15])
%xlim([-2 20])
%title(['data from ' datestr(jdi) ' to ' datestr(jdf)])
title('Sinais do lidar e radiossonda sobrepostos')
xlabel('Razão de mistura g/kg')
ylabel('Altitude')


grid on
























end

