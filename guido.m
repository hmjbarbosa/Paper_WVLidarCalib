%% Example script to plot the RCS for all channels in a single file
close all
fclose all
clear all

% datadir
basedir='/Users/hbarbosa/DATA/lidar/';

% change display for all plots
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2)

% exemple 1 ------------------------------------

% open only 1 file
fname=[basedir 'data/2014/9/14/RM1491400.000']; % night
%fname=[basedir 'data/2014/9/14/RM1491414.054']; % day

[head, phy] = profile_read(fname);

% print header from file
head

% print first channel 
head.ch(1)

% altudes in m
alt(:,1)=head.ch(1).binw*(1:head.ch(1).ndata);

% remove the background from all channels

for i=1:head.nch
	bg(i) = mean(phy(end-1000:end,i));
	P(:,i) = phy(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	

% grafico
figure(1)
clf; hold on; grid on
plot((Pr2),alt/1e3)
legend('an355','pc355','an387','pc387','pc408')
set(gca,'xscale','log')
xlabel('RCS [a.u.]')
ylabel('Altitude [km]')
ylim([0,20])
title(fname)


%fim