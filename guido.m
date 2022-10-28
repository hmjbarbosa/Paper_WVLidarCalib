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
% exemplo 1 ------------------------------------
% abre so um arquivo
fname='data/2014/9/15/RM1491500.001'; % noite
%fname='data/2014/9/15/RM1491514.205'; % dia

[head, phy] = profile_read(fname);

% conteudo do header na tela
head

% conteudo do primeiro canal
head.ch(1)

% removendo o ruido
alt(:,1)=head.ch(1).binw*(1:head.ch(1).ndata);

for i=1:head.nch
	bg(i) = mean(phy(end-1000:end,i));
	P(:,i) = phy(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	

% grafico
figure(1)
clf
plot((Pr2),alt/1e3)
legend('an355','pc355','an387','pc387','pc408')
set(gca,'xscale','log')
grid on
ylim([0,10])
title(fname)


%fim