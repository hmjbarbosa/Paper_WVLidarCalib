% clear all

addpath('./sc/')

% exemplo 2 ------------------------------------

% abrir varios arquivos
% flist=dirpath('data/2014/9/15/','RM1491500*');
flist=dirpath('data/2014/9/14/','RM1491420*');
flist=dirpath('D:/Dados/lidar/data/2014/9/14/','RM14914*');

[head, phy]=profile_read_many(flist);

% conteudo do cabecalho do arquivo #10

head(10)

% tamanho dos dados

size(phy(1).data)

% alturas
alt(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);
% horarios dos arquivos
for i=1:numel(head)
	jd(i)=head(i).jdi;
end
% mostra o primeiro horario na tela
datevec(jd(1))
datestr(jd(1))

% fazer figura nova

figure(1)
clf
gplot2(phy(1).data,[], jd,alt/1.e3)
title('canal 1')
ylabel('km')
xlabel('tempo')

% subtraindo o ruido dos 118 perfis

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
title('canal 1 sem BG')
ylabel('km')
xlabel('tempo')

% agora testando o efeito de fazer uma media no tempo

figure(3)
clf
plot(phy(1).data(:,1),alt/1.e3)
hold on
plot(mean(phy(1).data(:,1:10),2),alt/1.e3)
plot(mean(phy(1).data(:,:),2),alt/1.e3)
ylim([100 alt(end)/1e3])
ylabel('km')
legend('1 perfil','10 perfis','todos')

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

%%


figure(10)
clf
set(gcf,'units','normalized','outerposition',[0 0 0.8 1])
gplot2(Pr2(:,:,1),[0:0.1:1.5].*1e7, jd,alt/1.e3)
title('canal 1 sem BG')
ylabel('Altura [km]')
xlabel('Horário [LT]')
ylim([0 15])
datetick('x','keeplimits')

hc = colorbar;
ylabel(hc,'Sinal corrigido pela distância [U.A.]','fontsize',14)

saveas(gca,'figuras/RCS_2D.png')

% figure(11)
% clf
% set(gcf,'units','normalized','outerposition',[0 0 0.8 1])
% gplot2(log10(Pr2(:,:,1)),[6:0.1:7.5], jd,alt/1.e3)
% title('canal 1 sem BG')
% ylabel('km')
% xlabel('tempo')
% ylim([0 15])


%fim