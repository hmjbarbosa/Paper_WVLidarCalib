%% It seems this script does not do what it's name imply. 
% Only some random plots, probably not usefull and safe to delete. 

%close all
clear all

addpath('./libs/')
addpath('./sc/')

% datadir
basedir='/Users/hbarbosa/DATA/lidar/';

%[nfile, head, chphy]=profile_read_dir('2015/6/10',[],[],[],4000);

[nfile, head, chphy]=profile_read_dir([basedir 'data/2012_11_29_telescope/NocLight'],[],[],[],4000);

% background for each profile 
bg=mean(chphy(1).data(end-500:end,:))*0;

% range (m)
ndata=head(1).ch(1).ndata;
r=7.5*(1:ndata);
r2=r.*r;

clear jd phy
j=0;
%for i=1:nfile
for i=1308:1368
    j=j+1;
    jd(j)=head(i).jdi;
    phy(:,j)=chphy(1).data(:,i)-bg(i);
end

clear mr2
for i=1:j
  mr2(:,i)=r2;
end

figure(1); clf; 
gplot2(phy, [], jd, r/1e3)
caxis(quantile(reshape(phy,1,numel(phy)),[0.01 0.99]))
colorbar; grid
ylabel('km')
xlabel('time')
datetick('x','keeplimits')

figure(2); clf; 
gplot2(phy.*mr2, [], jd, r/1e3)
caxis(quantile(reshape(phy.*mr2,1,numel(phy)),[0.01 0.99]))
colorbar; grid
ylabel('km')
xlabel('time')
datetick('x','keeplimits')

figure(3); clf; 
plot(mod(jd,1)*24,'*')
ylabel('hour of day')
xlabel('profile #')
grid

figure(4); clf; 
plot(r/1e3, mean((phy.*mr2)'))
ylabel('P*R2')
xlabel('km')
grid

figure(5); clf; 
plot(r/1e3, mean(phy'))
ylabel('P')
xlabel('km')
grid
%