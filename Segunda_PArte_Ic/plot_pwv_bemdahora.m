close all
fclose all
clear all
%%
T_reserva_2012 = readtable('reserva_2012.txt','Delimiter',',');

jd_gps_utc_reserva_2012 = T_reserva_2012.dayofyear - 1 + datenum(2012,1,1); %encontramos erro pois estavamos somando um dia a mais
jd_gps_reserva_2012 = jd_gps_utc_reserva_2012 - datenum(00,00,00,04,00,00);   % hora local
a = datetime(datestr(jd_gps_reserva_2012));
pwv_gps_reserva_2012 = T_reserva_2012.pwv; pwv_gps_reserva_2012(pwv_gps_reserva_2012<0)=NaN;
%%
T = readtable('teste.txt','Delimiter',',');

jd_gps_utc = T.dayofyear - 1 + datenum(2011,1,1); %encontramos erro pois estavamos somando um dia a mais
jd_gps = jd_gps_utc - datenum(00,00,00,04,00,00);   % hora local
b = datetime(datestr(jd_gps));
pwv_gps = T.pwv; pwv_gps(pwv_gps<0)=NaN;
%%
T_reserva_2011 = readtable('reserva_2011.txt','Delimiter',',');

jd_gps_utc_reserva_2011 = T_reserva_2011.dayofyear - 1 + datenum(2011,1,1); %encontramos erro pois estavamos somando um dia a mais
jd_gps_reserva_2011 = jd_gps_utc_reserva_2011 - datenum(00,00,00,04,00,00);   % hora local
c = datetime(datestr(jd_gps_reserva_2011));
pwv_gps_reserva_2011 = T_reserva_2011.pwv; pwv_gps_reserva_2011(pwv_gps_reserva_2011<0)=NaN;
%%
d = [c;a];
pwv_gps_reserva = [pwv_gps_reserva_2011;pwv_gps_reserva_2012];
%%
figure(1)
clf
hold on
%plot(a,pwv_gps_reserva_2012*10,'.')
plot(b,pwv_gps*10,'.')
plot(d,pwv_gps_reserva *10,'.')
legend('Reserva','Embrapa')
title('Pwv todas estações')
xlabel('Data')
ylabel('PWV (mm)')
hold off
grid on