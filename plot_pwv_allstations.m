close all
fclose all;
clear all

%% Duke Reserve 

% read files for 2011 and 2012
disp('Reading reserva 2011')
T_reserva_2011 = readtable('reserva_2011.txt','Delimiter',',');
disp('Reading reserva 2012')
T_reserva_2012 = readtable('reserva_2012.txt','Delimiter',',');

% concatenate 
T_reserva = vertcat(T_reserva_2011, T_reserva_2012);

% convert year + dayofyear to Matlab's julian date
% old bug: we were not subtracting -1
jd_gps_utc_reserva = T_reserva.dayofyear - 1 + datenum(T_reserva.Year,1,1); 

% change to local time (Lidar is in LT)
jd_gps_reserva = jd_gps_utc_reserva - datenum(00,00,00,04,00,00); 

% convert julian date to datetime
dt_gps_reserva = datetime(datestr(jd_gps_reserva));

pwv_gps_reserva = T_reserva.pwv; 
pwv_gps_reserva(pwv_gps_reserva<0)=NaN;

%% Embrapa
disp('Reading embrapa 2011')
T_embrapa = readtable('2011_EMBP_PWV_met_data_int.txt','Delimiter',' ','MultipleDelimsAsOne',true);

% convert year + dayofyear to Matlab's julian date
% old bug: we were not subtracting -1
jd_gps_utc_embrapa = T_embrapa.dayofyear - 1 + datenum(T_embrapa.Year,1,1); 

% change to local time (Lidar is in LT)
jd_gps_embrapa = jd_gps_utc_embrapa - datenum(00,00,00,04,00,00);   

% convert julian date to datetime
dt_gps_embrapa = datetime(datestr(jd_gps_embrapa));

pwv_gps_embrapa = T_embrapa.pwv; 
pwv_gps_embrapa(pwv_gps_embrapa<0)=NaN;

%%
figure(1)
clf; hold on; grid on; box on
plot(dt_gps_embrapa, pwv_gps_embrapa*10,'.')
plot(dt_gps_reserva, pwv_gps_reserva *10,'.')
% BUG: labels were reserved
legend('Embrapa', 'Duke Reserve')
title('PWV All Stations')
ylabel('PWV (mm)')

%fim