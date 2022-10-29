clear all

addpath('./libs/')
addpath('./sc/')

% datadir
basedir='/Users/hbarbosa/DATA/lidar/';

% exemple 3 ------------------------------------

% open various files 
flist=dirpath([basedir 'data/2014/9/14/'],'RM1491420*');
[head, phy]=profile_read_many(flist,10, 0.004);

% print header from file #1
head(1)

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

% Plot with time mean profiles

figure(6)
clf; hold on; grid on
plot(alt/1.e3, Pr2)
set(gca,'yscale','log')
ylabel('RCS [a.u.]')
xlabel('Altitude [km]')
xlim([0,20])
legend('an355','pc355','an387','pc387','pc408')


%% 
% Elastic channel

ANmask = PcomBG(:,1) > 1.55 & PcomBG(:,1) < 250; % mV
disp(['Points in ANmask = ' num2str(sum(ANmask))])
PCmask = PcomBG(:,2) < 12.; % MHz
disp(['Points in PCmask = ' num2str(sum(PCmask))])

X_glue = P(ANmask & PCmask,1);
Y_glue = P(ANmask & PCmask,2);
cfun_glue = fit(X_glue,Y_glue,'poly1')

Signal_glued = cfun_glue(P(:,1));
Signal_glued(PCmask) = P(PCmask,2); 

figure(8)
clf; hold on; grid on
plot(X_glue,Y_glue,'o')
plot(cfun_glue)
text(0.08,10.3,['y = ' num2str(cfun_glue.p1,3) '*x + ' num2str(cfun_glue.p2,2)], ...
     'fontsize',14)

legend('Data','Fit')
xlabel('AN [mV]')
ylabel('PC [MHz]')
title('Channels 1 + 2')


figure(9)
clf; hold on; grid on
plot(alt/1e3,Pr2(:,1),'-b')
plot(alt/1e3,Pr2(:,2),'-r')
plot(alt/1e3,Signal_glued.*alt.*alt,'-')
plot(alt(ANmask)/1e3,Pr2(ANmask,1),'o')
plot(alt(PCmask)/1e3,Pr2(PCmask,2),'o')
xlim([0 10])
ylim([1e5 2e9])
set(gca,'yscale','log')
legend('an355','pc355','Glue')
ylabel('RCS [a.u.]')
xlabel('Altitude [km]')
title('Channels 1 + 2, 355 nm')

%%
% Nitrogen channel

ANmask = PcomBG(:,3) > 2.06 & PcomBG(:,3) < 50; % mV 
disp(['Points in ANmask = ' num2str(sum(ANmask))])
PCmask = PcomBG(:,4) < 12.; % MHz
disp(['Points in PCmask = ' num2str(sum(PCmask))])

X_glue = P(ANmask & PCmask,3);
Y_glue = P(ANmask & PCmask,4);
cfun_glue = fit(X_glue,Y_glue,'poly1')

Signal_glued = cfun_glue(P(:,3));
Signal_glued(PCmask) = P(PCmask,4); 

figure(18)
clf; hold on; grid on
plot(X_glue,Y_glue,'o')
plot(cfun_glue)
text(0.04,10.3,['y = ' num2str(cfun_glue.p1,3) '*x + ' num2str(cfun_glue.p2,2)], ...
     'fontsize',14)

legend('Data','Fit')
xlabel('AN [mV]')
ylabel('PC [MHz]')
title('Channels 3 + 4')

figure(19)
clf; hold on; grid on
plot(alt/1e3,Pr2(:,3),'-b')
plot(alt/1e3,Pr2(:,4),'-r')
plot(alt/1e3,Signal_glued.*alt.*alt,'-')
plot(alt(ANmask)/1e3,Pr2(ANmask,3),'o')
plot(alt(PCmask)/1e3,Pr2(PCmask,4),'o')
xlim([0 10])
ylim([1e5 2e8])
set(gca,'yscale','log')
legend('an387','pc387','Glue')
ylabel('RCS [a.u.]')
xlabel('Altitude [km]')
title('Channels 3 + 4, 387 nm')

%fim