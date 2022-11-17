%close all
fclose all
clear all
%%

% change display for all plots
set(0, 'DefaultAxesFontSize', 16, 'DefaultLineLineWidth', 2)


%prepro_dir = 'C:\Users\guido\Google Drive\IC Henrique\prepro\';
%prepro_dir = 'G:\.shortcut-targets-by-id\10F2V1DwX9Abl02tg9-a2rlR8Rp5yM4tm\prepro\'
prepro_dir = './prepro/'
figure_dir = './prepro/figs/'

flist = dir([prepro_dir 'Calib*.mat']);

%%
jdi = NaN(length(flist),1);
p1 = jdi;
p2 = jdi;

for i = 1:length(flist)
    
    i
   
    clear Data_temp
    Data_temp = load([prepro_dir flist(i).name]);
    jdi(i,1) = Data_temp.jd(1);
    p1(i,1) = Data_temp.constante_cfit.p1;
    p2(i,1) = Data_temp.constante_cfit.p2;
    sse(i,1) = Data_temp.constante_gof.sse;
    rsquare(i,1) = Data_temp.constante_gof.rsquare;
    dfe(i,1) = Data_temp.constante_gof.dfe;
    adjrsquare(i,1) = Data_temp.constante_gof.adjrsquare;
    rmse(i,1) = Data_temp.constante_gof.rmse;
    X = confint(Data_temp.constante_cfit);
    erro_p1(i,1) = (X(2,1)-X(1,1))/4;
    erro_p2(i,1) = (X(2,2)-X(1,2))/4;
    s_max250(i,1) = Data_temp.s_max250(1);
    s_max3000(i,1) = Data_temp.s_max3000(1);
end

dti = datetime(datestr(jdi));


%% new statistical tests

Flag_T = array2table(zeros(length(jdi), 1),'VariableNames',{'FitGood'});
Flag_T.jdi = jdi(:);

Flag_T.eval_p1 = abs(erro_p1./p1)<0.2;
disp(['good p1 = ' num2str(sum(Flag_T.eval_p1)) '  bad p1 = ' num2str(sum(~Flag_T.eval_p1))])

Flag_T.eval_p2 = abs(p2./erro_p2)<1 ; 
disp(['good p2 = ' num2str(sum(Flag_T.eval_p2)) '  bad p2 = ' num2str(sum(~Flag_T.eval_p2))])

Flag_T.eval_r2 = rsquare>0.8;
disp(['good r2 = ' num2str(sum(Flag_T.eval_r2)) '  bad r2 = ' num2str(sum(~Flag_T.eval_r2))])

Flag_T.neb = s_max3000./s_max250<1;
disp(['sem neblina = ' num2str(sum(~Flag_T.neb)) '  com neb = ' num2str(sum(Flag_T.neb))])

Flag_T.FitGood = Flag_T.eval_p1  &  Flag_T.eval_p2 & Flag_T.eval_r2 & (~Flag_T.neb);
disp(['good geral = ' num2str(sum(Flag_T.FitGood)) '  bad geral = ' num2str(sum(~(Flag_T.FitGood)))])


%% qualidade do ajuste                            

figure (10)
clf 
subplot(2,2,1); hold on; grid on; box on
edges = [0:0.05:1];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(rsquare(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(rsquare(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('R2')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,2); hold on; grid on; box on
edges = [0:0.15:4.5];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(rmse(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(rmse(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('RMSE')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,3); hold on; grid on; box on
edges = [0:2:50];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(dfe(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(dfe(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('NGL')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,4); hold on; grid on; box on
edges = [0:10:150];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(sse(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(sse(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('SSE')
ylabel('freq relativa (0-1)')
legend('ruim','bom')


%% P1 contra qualidade do ajuste                            

figure (20)
clf
subplot(2,2,1); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p1(tmp),rsquare(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p1(tmp),rsquare(tmp),'bo')
xlabel('p1')
ylabel('R2')
legend('ruim','bom')
%%
subplot(2,2,2); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p1(tmp),rmse(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p1(tmp),rmse(tmp),'bo')
xlabel('p1')
ylabel('RMSE')
legend('ruim','bom')
%%
subplot(2,2,3); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p1(tmp),dfe(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p1(tmp),dfe(tmp),'bo')
xlabel('p1')
ylabel('NGL')
legend('ruim','bom')
%%
subplot(2,2,4); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p1(tmp),sse(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p1(tmp),sse(tmp),'bo')
xlabel('p1')
ylabel('SSE')
legend('ruim','bom')

%% P2 contra qualidade do ajuste                            

figure (30)
clf
subplot(2,2,1); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p2(tmp),rsquare(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p2(tmp),rsquare(tmp),'bo')
xlabel('p2')
ylabel('R2')
legend('ruim','bom')
%%
subplot(2,2,2); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p2(tmp),rmse(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p2(tmp),rmse(tmp),'bo')
xlabel('p2')
ylabel('RMSE')
legend('ruim','bom')
%%
subplot(2,2,3); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p2(tmp),dfe(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p2(tmp),dfe(tmp),'bo')
xlabel('p2')
ylabel('NGL')
legend('ruim','bom')
%%
subplot(2,2,4); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p2(tmp),sse(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p2(tmp),sse(tmp),'bo')
xlabel('p2')
ylabel('SSE')
legend('ruim','bom')


%% P1 e P2 e os erros

%essa figura acho que mostra bem que os dias de fit bom não etm coef lin alto
figure (40)
clf
subplot(2,2,1); hold on; grid on; box on
edges = [-15:1:15];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(p2(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(p2(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('p2')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,2); hold on; grid on; box on
edges = [0:1:15];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(erro_p2(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(erro_p2(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('erro p2')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,3); hold on; grid on; box on
edges = [0:25:1500];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(p1(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(p1(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('p1')
ylabel('freq relativa (0-1)')
legend('ruim','bom')
%%
subplot(2,2,4); hold on; grid on; box on
edges = [0:25:600];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(erro_p1(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(erro_p1(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
xlabel('erro p1')
ylabel('freq relativa (0-1)')
legend('ruim','bom')



%% correlacao entre p1 e p2 e os erros

figure(300)
clf
subplot(2,2,1); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(erro_p1(tmp)./p1(tmp),p1(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(erro_p1(tmp)./p1(tmp),p1(tmp),'bo')
xlabel('Erro relativo em p1 (0-1)')
ylabel('p1 (g/kg)')
legend('ruim','bom')
%%
subplot(2,2,2); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p2(tmp)./erro_p2(tmp),p2(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p2(tmp)./erro_p2(tmp),p2(tmp),'bo')
xlabel('p2 / Erro em p2')
ylabel('p2 (g/kg)')
legend('ruim','bom')
%%
subplot(2,2,3); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(p1(tmp),p2(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(p1(tmp),p2(tmp),'bo')
xlabel('p1')
ylabel('p2')
legend('ruim','bom')
%%
subplot(2,2,4); hold on; grid on; box on
tmp = Flag_T.FitGood == 0; 
plot(erro_p1(tmp),erro_p2(tmp),'ro')
tmp = Flag_T.FitGood == 1; 
plot(erro_p1(tmp),erro_p2(tmp),'bo')
xlabel('Erro p1')
ylabel('Erro p2')
legend('ruim','bom')
%%

%% Figura 4 do paper - histograma das contantes

figure(100)
clf
%subplot(2,2,1); 
hold on; grid on; box on
tmp = Flag_T.FitGood == 1; 
histogram(p1,linspace(-500, 1500, 50),'facecolor', 'red' , 'Normalization', 'probability')
histogram(p1(tmp),linspace(-500, 1500, 50), 'facecolor', 'blue' , 'Normalization', 'probability')
%histogram(p1(~tmp),linspace(-500, 1500, 50), 'Normalization', 'probability')
xlabel('Calibration Constant (g/kg)')
ylabel('Frequency (0-1)')
legend(['All = '      num2str(mean(p1),'%.0f')       ' \pm ' num2str(std(p1),'%.0f') ' g/kg'],...
       ['Good = '     num2str(mean(p1(tmp )),'%.0f') ' \pm ' num2str(std(p1(tmp )),'%.0f') ' g/kg'],...
       'Location','northwest')
%       ['Rejected = ' num2str(mean(p1(~tmp)),'%.0f') ' \pm ' num2str(std(p1(~tmp)),'%.0f') ' g/kg'],...
t=text(.98, .02, '  a  ', 'units', 'normalize',...
       'Horiz','right','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);

saveas(figure(100),[figure_dir '_hist.png' ])
set(figure(100),'paperunits','points','papersize',[560 420])
print(figure(100),[figure_dir '_hist.pdf' ],'-dpdf','-bestfit')

%%

figure(101)
clf
%subplot(2,2,1); 
hold on; grid on; box on
tmp = Flag_T.FitGood == 1; 
histogram(p2,linspace(-10, 10, 50),'facecolor', 'red' , 'Normalization', 'probability')
histogram(p2(tmp),linspace(-10, 10, 50), 'facecolor', 'blue' , 'Normalization', 'probability')
%histogram(p2(~tmp),linspace(-500, 1500, 50), 'Normalization', 'probability')
xlabel('Calibration Bias (g/kg)')
ylabel('Frequency (0-1)')
legend(['All = '      num2str(mean(p2),'%.2f')       ' \pm ' num2str(std(p2),'%.2f') ' g/kg'],...
       ['Good = '     num2str(mean(p2(tmp )),'%.2f') ' \pm ' num2str(std(p2(tmp )),'%.2f') ' g/kg'],...
       'Location','northwest')
%       ['Rejected = ' num2str(mean(p2(~tmp)),'%.2f') ' \pm ' num2str(std(p2(~tmp)),'%.2f') ' g/kg'],...
t=text(.98, .02, '  b  ', 'units', 'normalize',...
       'Horiz','right','Vert', 'bot',...
       'backgroundcolor','white','edgecolor','black','fontsize',12);

saveas(figure(101),[figure_dir '_hist_p2.png' ])
set(figure(101),'paperunits','points','papersize',[560 420])
print(figure(101),[figure_dir '_hist_p2.pdf' ],'-dpdf','-bestfit')




%% Figure 5 do paper - serie temporal das constantes

figure (200)
clf; hold on; grid on; box on
mask_temp1 = jdi<=jdi(137);
mask_temp2 = jdi>=jdi(138); %192
tmp = Flag_T.FitGood == 1; 

jdi1 = mean(jdi(mask_temp1));
p1t1 = p1(mask_temp1 & tmp);
t1 = jdi(mask_temp1 & tmp)-jdi1;
[c,d] = fit(t1,p1t1,'poly1')
dtifine1=linspace(min(dti(mask_temp1 & tmp)), max(dti(mask_temp1 & tmp)), 50)';
t1fine=linspace(min(t1),max(t1), 50)';
conf1 = tinv(0.975, d.dfe)*sqrt((d.sse/d.dfe/(d.dfe+2))*((d.dfe+2)*0+1+(t1fine - mean(t1)).^2 / var(t1)));

jdi2 = mean(jdi(mask_temp2));
p1t2 = p1(mask_temp2 & tmp);
t2 = jdi(mask_temp2 & tmp)-jdi2;
[e,f] = fit(t2,p1t2,'poly1')

dtifine2=linspace(min(dti(mask_temp2 & tmp)), max(dti(mask_temp2 & tmp)), 50)';
t2fine=linspace(min(t2),max(t2), 50)';
conf2 = tinv(0.975, f.dfe)*sqrt((f.sse/f.dfe/(f.dfe+2))*((f.dfe+2)*0+1+(t2fine - mean(t2)).^2 / var(t2)));

plot(dti(tmp ),p1(tmp ),'o')  % bons
plot(dti((~tmp)),p1((~tmp)),'.r') % ruins
plot(dti(mask_temp1 & tmp),c(t1), '-c','linewidth',2) % fit < 2013
plot(dti(mask_temp2 & tmp),e(t2), '-m','linewidth',2) % fit > 2014
                                                     
% confidence intervals
plot(dtifine1,c(t1fine)+conf1, '--c','linewidth',1.5) % fit < 2013
plot(dtifine1,c(t1fine)-conf1, '--c','linewidth',1.5) 
plot(dtifine2,e(t2fine)+conf2, '--m','linewidth',1.5) % fit < 2013
plot(dtifine2,e(t2fine)-conf2, '--m','linewidth',1.5) 

xlabel('Date')
ylabel('Calibration constant (g/kg)')
legend('Good calibrations','Rejected calibrations','2011 - 2012','2014 - 2016')

% get standard error instead of 95% conf interval
level=2*tcdf(-1,d.dfe);
cc1sig=confint(c,1-level);
stdp1=cc1sig(:,1)-c.p1;
stdp2=cc1sig(:,2)-c.p2;
lineA = sprintf('period 1 = (%.2f \\pm %.2f)*(JD-JDmean) + (%.0f \\pm %.0f)' , ...
                c.p1, stdp1(2), c.p2, stdp2(2) );

level=2*tcdf(-1,f.dfe);
cc1sig=confint(e,1-level);
stdp1=cc1sig(:,1)-e.p1;
stdp2=cc1sig(:,2)-e.p2;
lineB = sprintf('period 2 = (%.2f \\pm %.2f)*(JD-JDmean) + (%.0f \\pm %.0f)' , ...
                e.p1, stdp1(2), e.p2, stdp2(2) );

text(0.1, 0.1, {lineA, lineB}, 'units', 'normalize',...
     'backgroundcolor','white','edgecolor','black','fontsize',12)
print(figure(200),[figure_dir 'fitboas_calibracoes_.png' ],'-dpng')
set(figure(200),'paperunits','points','papersize',[560 420])
print(figure(200),[figure_dir 'fitboas_calibracoes_.pdf' ],'-dpdf','-bestfit')

figure(123) % sim
clf; hold on; grid on; box on
plot(t1,p1t1,'o')
plot(c,'predfunc')
plot(c,'predobs')
xlabel('dias relativo ao meio do periodo #1')
ylabel('calibracao (g/kg)')

%% Igual a Figure 5 do paper, mas para o termo linear

figure (500)
clf; hold on; grid on; box on
mask_temp1 = jdi<=jdi(137);
mask_temp2 = jdi>=jdi(192);
tmp = Flag_T.FitGood == 1; 

jdi1 = mean(jdi(mask_temp1));
p2t1 = p2(mask_temp1 & tmp);
t1 = jdi(mask_temp1 & tmp)-jdi1;
[c,d] = fit(t1,p2t1,'poly1');

jdi2 = mean(jdi(mask_temp2));
p2t2 = p2(mask_temp2 & tmp);
t2 = jdi(mask_temp2 & tmp)-jdi2;
[e,f] = fit(t2,p2t2,'poly1');

%plot(dti(tmp ),p2(tmp ),'o')  % bons
errorbar(dti(tmp ),p2(tmp ),erro_p2(tmp),'o')  % bons
plot(dti((~tmp)),p2((~tmp)),'.r') % ruins
plot(dti(mask_temp1 & tmp),c(t1), '-','linewidth',2) % fit < 2013
plot(dti(mask_temp2 & tmp),e(t2), '-','linewidth',2) % fit > 2014

xlabel('Date')
ylabel('Calibration bias (g/kg)')
legend('Good calibrations','Rejected calibrations','2011 - 2012','2014 - 2016')

% get standard error instead of 95% conf interval
level=2*tcdf(-1,d.dfe);
cc1sig=confint(c,1-level);
stdp1=cc1sig(:,1)-c.p1;
stdp2=cc1sig(:,2)-c.p2;
lineA = sprintf('period 1 = (%.2g \\pm %.2g)*(JD-JDmean) + (%.2g \\pm %.2g)' , ...
                c.p1, stdp1(2), c.p2, stdp2(2) );

level=2*tcdf(-1,f.dfe);
cc1sig=confint(e,1-level);
stdp1=cc1sig(:,1)-e.p1;
stdp2=cc1sig(:,2)-e.p2;
lineB = sprintf('period 2 = (%.2g \\pm %.2g)*(JD-JDmean) + (%.2g \\pm %.2g)' , ...
                e.p1, stdp1(2), e.p2, stdp2(2) );

text(0.1, 0.1, {lineA, lineB}, 'units', 'normalize',...
     'backgroundcolor','white','edgecolor','black','fontsize',12)
%print(figure(500),[figure_dir 'fitboas_calibracoes_.png' ],'-dpng')
%set(figure(500),'paperunits','points','papersize',[560 420])
%print(figure(500),[figure_dir 'fitboas_calibracoes_.pdf' ],'-dpdf','-bestfit')

return
%%
%aqui elimina os dados desse plot menores do que 1.5
figure (327)
clf 
hold on 
plot(dti,s_max3000'./s_max250')
hold off
grid on


%     jdi(i) = Data_temp.jd(1);
%     p1(i) = Data_temp.constante_cfit.p1;
%     p2(i) = Data_temp.constante_cfit.p2;
%     sse(i) = Data_temp.constante_gof.sse;
%     rsquare(i) = Data_temp.constante_gof.rsquare;
%     dfe(i) = Data_temp.constante_gof.dfe;
%     adjrsquare(i) = Data_temp.constante_gof.adjrsquare;
%     rmse(i) = Data_temp.constante_gof.rmse;


%%
mask_temp = Flag_T.year == 2011 & Flag_T.month >= 7
mask_temp0 = mask_temp  |Flag_T.year == 2012
mask_temp10 = mask_temp0|Flag_T.year == 2014 & Flag_T.month == 3|Flag_T.year == 2014 & Flag_T.month == 7|Flag_T.year == 2015 & Flag_T.month == 2|Flag_T.year == 2015 & Flag_T.month == 7 & Flag_T.day == 24|Flag_T.year == 2015 & Flag_T.month == 8 & Flag_T.day == 7|Flag_T.year == 2015 & Flag_T.month == 9 & Flag_T.day == 30|Flag_T.year == 2016 & Flag_T.month == 2|Flag_T.year == 2016 & Flag_T.month == 7|Flag_T.year == 2016 & Flag_T.month == 8 
mask_temp1 =mask_temp10 & Flag_T.FitGood == 0;
mask_temp2 = mask_temp1 & Flag_T.cloud_2km == 1;
mask_temp3 = mask_temp2 & Flag_T.cloud_5km == 1;
%%
%mask_season = Flag_T.FitGood == 0
mask_season1 = Flag_T.month == 6;
mask_season2 = mask_season1|Flag_T.month ==7;
mask_season3 = mask_season2|Flag_T.month ==8;
mask_season4 = mask_season3|Flag_T.month ==9;
mask_season5 = mask_season4 & Flag_T.FitGood == 0;
mask_season6 = ~mask_season4 & Flag_T.FitGood == 0;

%% 




%%
figure(1)
clf
hold on
plot(dti(mask_temp10),p1(mask_temp10),'.')
plot(dti(mask_temp1),p1(mask_temp1),'or')
plot(dti(mask_temp2),p1(mask_temp2),'xg','markersize',14)
plot(dti(mask_temp3),p1(mask_temp3),'s','markersize',14)
legend('todas as calibrações','apenas as com bons ajustes','ajuste bom & sem nuvem<2km','ajuste bom & sem nuvem< 5 km', 'ajuste bom & sem nuvem<10km')
title('Calibrações')
xlabel('data da calibração')
ylabel('constante de calibração')
hold off
grid on


figure(2)
clf
hold on
plot(dti(mask_temp10),p2(mask_temp10),'.')
plot(dti(mask_temp1),p2(mask_temp1),'or')
plot(dti(mask_temp2),p2(mask_temp2),'xg','markersize',14)
plot(dti(mask_temp3),p2(mask_temp3),'s','markersize',14)


hold off
grid on

%%

figure(3)
clf
histogram(p1(mask_temp1),20)
xlabel('Constante de calibração')
ylabel('ocorrências')
legend('Média = 808 Desvpad = 85')
title('Histograma com os dados bons')
%%
figure (4)
clf
histogram(p1(mask_temp10),20)
xlabel('Constante de calibração')
ylabel('ocorrências')
legend('Média = 611 Desvpad = 310')
title('Histograma com todos os dados')
%%
x0 = datenum(2011,07,1,0,0,0);
x = jdi(mask_temp3) - x0;
y = p1(mask_temp3)
%%
g = []
for i = 1:length(jdi(mask_temp1))
    g(i) = i - 1
end
%%
calib_cfit = fit(x',y','poly1')

figure (5)
clf; hold on; grid on; box on
plot(dti(mask_temp10),p1(mask_temp10)/1e3,'.')
plot(dti(mask_temp1),p1(mask_temp1)/1e3,'or')
plot(dti(mask_temp2),p1(mask_temp2)/1e3,'xg','markersize',14)
%plot(dti(mask_temp3),p1(mask_temp3)/1e3,'s','markersize',14)
%plot(x+x0,calib_cfit(x))
plot(dti(mask_temp3),calib_cfit(x)/1e3)
legend('todas as calibrações (n=230)','apenas as com bons ajustes(n=168)','ajuste bom & sem nuvem<2km (n=118)','ajuste bom & sem nuvem< 5 km (n = 108)','calib_cfit(x) = 0.05*x + 790')
xlabel('Date')
ylabel('Calibration Constant (kg/kg)')
%%
j_a = Flag_T.year == 2011 & Flag_T.month >= 7 & Flag_T.month <= 8 & Flag_T.FitGood == 0;
x = std(p1(j_a))
sigma_porcento = (x/mean(p1(j_a)))*100
%%
saveas(figure(1),[figure_dir '_boas_calibracoes_.png' ])
saveas(figure(2),[figure_dir '_bons_coeflin_.png' ])
saveas(figure(3),[figure_dir '_hist_bom_.png' ])
saveas(figure(4),[figure_dir '_hist_tudo_.png' ])
saveas(figure(5),[figure_dir 'ajuste_lin_todascalib_.png' ])


%%