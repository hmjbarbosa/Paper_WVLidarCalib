clear all

addpath('./libs/')
addpath('./sc/')

% datadir
basedir='/Users/hbarbosa/DATA/lidar/';

% exemplo 2 ------------------------------------

% open various files 
flist=dirpath([basedir 'data/2014/9/14/'],'RM1491420*');
[head, phy]=profile_read_many(flist);

% print header from file #1
head(10)

% print size of data
size(phy(1).data)

% altudes in m
alt(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);

% times in Matlab's julian date
for i=1:numel(head)
	jd(i)=head(i).jdi;
end

% print first time 
datestr(jd(1))

% fazer figura nova

figure(1)
clf
gplot2(phy(1).data,[], jd,alt/1.e3)
title('channel 1, uncorrected, mV')
ylabel('km')
xlabel('time')
datetick('x','keeplimits')

% remove the background from all channels and all profiles

for i=1:head(1).nch
	for j=1:numel(head)
		bg(i,j) = mean(phy(i).data(end-1000:end,j));
		P(:,j,i) = phy(i).data(:,j) - bg(i,j);
		Pr2(:,j,i) = P(:,j,i).*alt.*alt;
	end
end	

figure(2)
clf
gplot2(P(:,:,1),[], jd,alt/1.e3)
title('channel 1, BG corrected, mV')
ylabel('km')
xlabel('time')
datetick('x','keeplimits')

% agora testando o efeito de fazer uma media no tempo

figure(3)
clf; hold on
plot(phy(1).data(:,1),alt/1.e3)
plot(mean(phy(1).data(:,1:10),2),alt/1.e3)
plot(mean(phy(1).data(:,:),2),alt/1.e3)
ylim([100 alt(end)/1e3])
ylabel('km')
xlabel('channel 1, mV')
legend('1 profile','10 profile','all')
title('averaging')

% olhando o perfil inteiro, media no tempo

figure(4)
clf
plot(mean(phy(1).data(:,:),2),alt/1.e3)
hold on
plot(mean(phy(2).data(:,:),2),alt/1.e3)
plot(mean(phy(3).data(:,:),2),alt/1.e3)
plot(mean(phy(4).data(:,:),2),alt/1.e3)
plot(mean(phy(5).data(:,:),2),alt/1.e3)
ylabel('km')
legend('an355','pc355','an387','pc387','pc408')
ylim([0 20])
set(gca,'xscale','log')

%%


figure(10)
clf
%set(gcf,'units','normalized','outerposition',[0 0 0.8 1])
gplot2(Pr2(:,:,1),[0:0.1:1.5].*1e7, jd,alt/1.e3)
title('channel 1, RCS')
ylabel('Altitude [km]')
xlabel('Time [LT]')
ylim([0 15])
datetick('x','keeplimits')
hc = colorbar;


%fim