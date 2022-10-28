clear all

%%
addpath('./sc/')
%%
%	change	display	for	all	plots
set(0,	'DefaultAxesFontSize',	12,	'DefaultLineLineWidth',	2)
%%
% exemplo 4 ------------------------------------

% abrir varios arquivos
% flist=dirpath('data/2014/9/15/','RM1491500*');
flist=dirpath('data/2014/9/14/','RM1491420*');

[head, phy]=profile_read_many(flist,10, 0.004);

% conteudo do cabecalho do arquivo #10

head(10)

% tamanho dos dados

size(phy(1).data)

% alturas
alt(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);

% altitude  
altitude = alt + head(1).alt;

% horarios dos arquivos
for i=1:numel(head)
	jd(i)=head(i).jdi;
end
% mostra o primeiro horario na tela
datevec(jd(1))
datestr(jd(1))

% subtraindo o ruido dos 118 perfis

for i=1:head(1).nch
	PcomBG(:,i) = mean(phy(i).data(:,:),2);

	bg(i) = mean(PcomBG(end-5000:end,i));
	desvio(i) = std(PcomBG(end-5000:end,i), 1);

	P(:,i) = PcomBG(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	

%%
% calibrando vapor de agua

h2o = P(:,5) ./ P(:,3) ; % PC_H2o / AN_N2

figure(9)
clf; hold on
ylim([0 7])
plot(h2o*12, altitude/1e3,'-')

plot(smooth(h2o,11)*12, altitude/1e3,'-')

xlabel('h2o')
ylabel('km')
grid on

% lendo a sondagem
% sonda=csvread('sonda.csv',2)
sonda=csvread('sonda2.csv',2)

plot(sonda(:,6), sonda(:,2)/1e3,'o')


%%

% interpolar dados lidar nos níveis da sonda
sm = 11;
h2o_lidar = interp1q(smooth(altitude,sm),smooth(h2o,sm), sonda(:,2));

figure(15)
clf
plot(h2o,altitude./1000,'-',h2o_lidar,sonda(:,2)./1000,'*k')
grid on
ylim([0 7])


mask = sonda(:,2) > 150 & sonda(:,2) < 10000 ;

X = h2o_lidar(mask);
Y = sonda(mask,6);

cfun_calib = fit(X,Y,'poly1')

figure(16)
clf;
plot(X,Y,'*')
hold on
plot(X,cfun_calib(X))
hold off
grid on
xlabel('Sinal do lidar [u.a.]')
ylabel('Razão de mistura de H_2O (Sonda) [g/kg]')
title('texto')

saveas(figure(16),'figuras/calib.png')

%%

h2o_calibado = cfun_calib(h2o);

figure(90)
clf; hold on
ylim([0 7])
plot(h2o_calibado, altitude/1e3,'-')

plot(smooth(h2o_calibado,sm), altitude/1e3,'-')

xlabel('h2o')
ylabel('km')
grid on

% lendo a sondagem
% sonda=csvread('sonda.csv',2)
sonda=csvread('sonda2.csv',2)

plot(sonda(:,6), sonda(:,2)/1e3,'o')

%%
%%%%%%%%%%%%%%%%%%%


% 
% % interpolar sonda nos níveis do lidar
% h2o_sonde = interp1q(sonda(:,2),sonda(:,6),altitude)
% 
% figure(10)
% clf
% plot(sonda(:,6),sonda(:,2)./1000,'or',h2o_sonde,altitude./1000,'-')
% grid on
% ylim([0 7])
% 
% 
% mask = altitude > 150 & altitude < 10000 ;
% 
% X = h2o(mask);
% Y = h2o_sonde(mask);
% 
% cfun = fit(X,Y,'poly1')
% 
% figure(10)
% clf;
% plot(X,Y,'.')
% hold on
% plot(X,cfun(X))
% hold off
% grid on
% xlabel('Sinal do lidar [u.a.]')
% ylabel('Razão de mistura [g/kg]')
% title('texto')
% 
% % % saveas(figure(10),'figuras/calib.png')
% 




%fim