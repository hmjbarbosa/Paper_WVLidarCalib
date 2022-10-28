%%
close all
fclose all
clear all

% change display for all plots
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2)
%%
addpath('./IC Henrique/')
addpath('./data/')
addpath('./sc/')
%%
% abrir varios arquivos
% flist=dirpath('data/2014/9/15/','RM1491500*');
flist=dirpath('data/2014/9/14/','RM1491420*');

[head, phy]=profile_read_many(flist);

%%
% alturas
alt(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);

%%
% altitude  
altitude = alt + head(1).alt;

%%
% horarios dos arquivos
for i=1:numel(head)
	jd(i)=head(i).jdi;
end

%%
%tirando o bg
for i=1:head(1).nch
	PcomBG(:,i) = mean(phy(i).data(:,:),2);

	bg(i) = mean(PcomBG(end-5000:end,i));
	desvio(i) = std(PcomBG(end-5000:end,i), 1);

	P(:,i) = PcomBG(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	


%% achar a regiao onde vale o photocount e o analogic 
ANmask = PcomBG(:,3) > 2.06 & PcomBG(:,3) < 50; % em mV do AN
sum(ANmask)
PCmask = PcomBG(:,4) < 12.; % em MHz do PC
sum(PCmask)

X_glue = P(ANmask & PCmask,3);
Y_glue = P(ANmask & PCmask,4);
cfun_glue = fit(X_glue,Y_glue,'poly1')

Signal_glued = cfun_glue(P(:,3));Signal_glued(PCmask) = P(PCmask,4); 

%% razão de mistura do lidar
% % % % % h2o = P(:,5) ./ Signal_glued  ; % PC_H2o / N2

%%
% lendo a sondagem
% sonda=csvread('sonda.csv',2)
sonda=csvread('sonda2.csv',2)

%% 
% % % % % % interpolar dados lidar nos níveis da sonda
% % % % % sm = 11;
% % % % % h2o_lidar = interp1q(smooth(altitude,sm),smooth(h2o,sm), sonda(:,2));
%%

[constante_cfit] = calib_function_v0(sonda(:,2),sonda(:,6),altitude,Signal_glued,P(:,5));


%%


% % % % % mask = sonda(:,2) > 150 & sonda(:,2) < 10000 ;
% % % % % 
% % % % % X = h2o_lidar(mask);
% % % % % Y = sonda(mask,6);
% % % % % 
% % % % % cfun_calib = fit(X,Y,'poly1')
% % % % % constante_de_calibracao = cfun_calib.p1

%%



