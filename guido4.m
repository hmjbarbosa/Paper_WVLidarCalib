%% Example script for doing the WV calibration. Reads 1h of
% nighttime data and performs the calibration of 
%
%             h2o_lidar = PC_H2O / AN_N2 
clear all

addpath('./libs/')
addpath('./sc/')

% datadir
basedir='/Users/hbarbosa/DATA/lidar/';

%	change	display	for	all	plots
set(0,	'DefaultAxesFontSize',	12,	'DefaultLineLineWidth',	2)

% exemplo 4 ------------------------------------

% open various files 
flist=dirpath([basedir 'data/2014/9/14/'],'RM1491420*');
[head, phy]=profile_read_many(flist,10, 0.004);

% print header from file #1
head(1)

% print size of data
size(phy(1).data)

% altudes in m
alt(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);

% altitude above sea level
altitude = alt + head(1).alt;

% times in Matlab's julian date
for i=1:numel(head)
	jd(i)=head(i).jdi;
end

% print first time 
datestr(jd(1))

% remove the background from all channels

for i=1:head(1).nch
    % average of all datafiles
	PcomBG(:,i) = mean(phy(i).data(:,:),2);

    % comput BG from last 5000 bins (37.5 km)
	bg(i) = mean(PcomBG(end-5000:end,i));
	desvio(i) = std(PcomBG(end-5000:end,i), 1);

    % remove BG and correct range
	P(:,i) = PcomBG(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	

% read sounding 
% sonda=csvread('sonda.csv',2)
sonda=csvread('sonda2.csv',2);

%%
% Calibration

h2o = P(:,5) ./ P(:,3) ; % PC_H2o / AN_N2

% interpolate lidar data to sonde levels
sm = 11;
h2o_lidar = interp1q(smooth(altitude,sm),smooth(h2o,sm), sonda(:,2));

% use data above 150 m and up to 10km
mask = sonda(:,2) > 150 & sonda(:,2) < 10000;

% linear fit
X = h2o_lidar(mask);
Y = sonda(mask,6);
cfun_calib = fit(X,Y,'poly1')

% apply the calibration to the lidar profile
h2o_calib = cfun_calib(h2o);

%% plots

figure(9)
clf; hold on; grid on
plot(sonda(:,6), sonda(:,2)/1e3,'om')
xlim([0 20])
ylim([0 7])
xlabel('h2o sonde (g/kg)')
legend('h2o sonde')

figure(10)
clf; hold on; grid on
plot(h2o, altitude/1e3,'-')
plot(smooth(h2o,11), altitude/1e3,'-')
plot(h2o_lidar, sonda(:,2)./1e3,'*k')
xlim([0 2])
ylim([0 7])
xlabel('h2o lidar uncalibrated [a.u.]')
ylabel('km')
legend('h2o lidar', 'h2o lidar smoothed', 'h2o lidar interp')

figure(11)
clf;
plot(X,Y,'*')
hold on
plot(X,cfun_calib(X))
hold off
grid on
xlabel('h2o lidar uncalibrated [a.u.]')
ylabel('h2o sonde [g/kg]')
title('linear calibration')
text(0.04,10.3,['y = ' num2str(cfun_calib.p1,3) '*x + ' num2str(cfun_calib.p2,2)],...
     'fontsize',14)

%%

figure(12)
clf; hold on; grid on
plot(h2o_calib, altitude/1e3,'-')
plot(smooth(h2o_calib,sm), altitude/1e3,'-')
plot(sonda(:,6), sonda(:,2)/1e3,'o')
xlim([-1 20])
ylim([0 7])
title('Calibrated Sinal')
xlabel('H2O Mixing ratio [g/Kg]')
ylabel('Altitude a.s.l. [km]')
legend('lidar','sinal smooth','sonde')

%fim