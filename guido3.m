clear all

addpath('./sc/')

% exemplo 3 ------------------------------------

% abrir varios arquivos
% flist=dirpath('data/2014/9/15/','RM1491500*');
flist=dirpath('data/2014/9/14/','RM1491420*');

[head, phy]=profile_read_many(flist,10, 0.004);

% conteudo do cabecalho do arquivo #10

head(10)
%%
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
%%
% subtraindo o ruido dos 118 perfis

for i=1:head(1).nch
	PcomBG(:,i) = mean(phy(i).data(:,:),2);

	bg(i) = mean(PcomBG(end-5000:end,i));
	desvio(i) = std(PcomBG(end-5000:end,i), 1);

	P(:,i) = PcomBG(:,i) - bg(i);
	Pr2(:,i) = P(:,i).*alt.*alt;
end	

% calcular a regiao onde confiamos em AN e PC
% por enquanto, para o elastico (canais 1 e 2) apenas
%% 
ANmask = PcomBG(:,1) > 1.55 & PcomBG(:,1) < 250; % em mV do AN
sum(ANmask)
PCmask = PcomBG(:,2) < 12.; % em MHz do PC
sum(PCmask)

% olhando o perfil inteiro, media no tempo

figure(6)
clf
plot(Pr2,alt/1.e3)
ylabel('km')
ylim([0,10])
legend('an355','pc355','an387','pc387','pc408')
grid on


% tentando colar na mao o AN e PC do 355
figure(7)
clf; hold on
%xlim([0 4000])
plot(alt/1e3,Pr2(:,1),'-')
plot(alt/1e3,Pr2(:,2),'-')

plot(alt(ANmask)/1e3,(P(ANmask,1)*61.0-0.5413).*alt(ANmask).*alt(ANmask),'-')

plot(alt(ANmask)/1e3,Pr2(ANmask,1),'o')
plot(alt(PCmask)/1e3,Pr2(PCmask,2),'o')

legend('an355','pc355')
grid on

figure(8)
clf;
plot(P(ANmask & PCmask,1), P(ANmask & PCmask,2),'o')
xlabel('AN')
ylabel('PC')
grid on


%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
ANmask = PcomBG(:,3) > 2.06 & PcomBG(:,3) < 50; % em mV do AN
sum(ANmask)
PCmask = PcomBG(:,4) < 12.; % em MHz do PC
sum(PCmask)


% tentando colar na mao o AN e PC do 355
figure(17)
clf; hold on
%xlim([0 4000])
plot(alt/1e3,Pr2(:,3),'-')
plot(alt/1e3,Pr2(:,4),'-')

% plot(alt(ANmask)/1e3,(P(ANmask,3)*67.747-0.10173).*alt(ANmask).*alt(ANmask),'-')

plot(alt(ANmask)/1e3,Pr2(ANmask,3),'.')
plot(alt(PCmask)/1e3,Pr2(PCmask,4),'.')

legend('an387','pc387')
grid on

%%

X_glue = P(ANmask & PCmask,3);
Y_glue = P(ANmask & PCmask,4);
cfun_glue = fit(X_glue,Y_glue,'poly1')

Signal_glued = cfun_glue(P(:,3));Signal_glued(PCmask) = P(PCmask,4); 

%%
figure(18)
clf;
plot(X_glue,Y_glue,'o')
hold on
% plot(X_glue,cfun_glue(X_glue))
plot(cfun_glue)
hold off
text(0.04,10.3,['y = ' num2str(cfun_glue.p1,3) '*x + ' num2str(cfun_glue.p2,2)])

legend('asfas','afs')
xlabel('AN')
ylabel('PC')
grid on

figure(19)
clf; hold on
%xlim([0 4000])
plot(alt/1e3,Pr2(:,3),'-')
plot(alt/1e3,Pr2(:,4),'-')

% plot(alt(ANmask)/1e3,(P(ANmask,3)*cfun_glue.p1 + cfun_glue.p2).*alt(ANmask).*alt(ANmask),'-')
plot(alt/1e3,Signal_glued.*alt.*alt,'-')

plot(alt(ANmask)/1e3,Pr2(ANmask,3),'.')
plot(alt(PCmask)/1e3,Pr2(PCmask,4),'.')

xlim([0 20])

% ylim([1e5 1e9])
% set(gca,'yscale','log')


legend('an387','pc387','Colado')
grid on











%fim