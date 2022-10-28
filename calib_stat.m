%close all
fclose all
clear all
%%


% [filenames_aux, pathname, filterindex] = uigetfile('.mat','MAT-files with gluing data (.mat)',...
%     ['C:\Users\guido\Google Drive\prepro\Calib*'], 'MultiSelect', 'on');
% filenames = cellstr(filenames_aux);


% Data_temp = load('D:\prepro\Calib_2011_02_12_20_00_00.mat');

%%
figure_dir = './figs/'
%%

%prepro_dir = 'C:\Users\guido\Google Drive\IC Henrique\prepro\';
%prepro_dir = 'G:\.shortcut-targets-by-id\10F2V1DwX9Abl02tg9-a2rlR8Rp5yM4tm\prepro\'
prepro_dir = './prepro/'
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


%%
dti = datetime(datestr(jdi));



%%

%%%%%%%%%%%
Flag_T_inp = readtable(['flags_calib_3.txt']);
Flag_T_inp.jd_dia = datenum(Flag_T_inp.year,Flag_T_inp.month,Flag_T_inp.day);


%%%

Mflag_temp = NaN(length(jdi),size(Flag_T_inp,2));
for  i = 1:length(jdi)
    Ind = find(Flag_T_inp.jd_dia == floor(jdi(i)));
    
    if length(Ind) == 1
        
        Mflag_temp(i,:) = table2array(Flag_T_inp(Ind,:));
        
    elseif length(Ind) > 2
        'Tabela com algum erro!!!'
        i
        break
    end
    
end

%% new statistical tests

%Flag_T = array2table(Mflag_temp,'VariableNames',Flag_T_inp.Properties.VariableNames);
Flag_T = array2table(zeros(length(jdi), 1),'VariableNames',{'FitGood'});
Flag_T.jdi = jdi(:);

Flag_T.eval_p1 = abs(erro_p1./p1)<0.2;
disp(['good p1 = ' num2str(sum(Flag_T.eval_p1)) '  bad p1 = ' num2str(sum(~Flag_T.eval_p1))])

Flag_T.eval_p2 = abs(p2./erro_p2)<3 ; 
disp(['good p2 = ' num2str(sum(Flag_T.eval_p2)) '  bad p2 = ' num2str(sum(~Flag_T.eval_p2))])

Flag_T.eval_r2 = rsquare>0.8;
disp(['good r2 = ' num2str(sum(Flag_T.eval_r2)) '  bad r2 = ' num2str(sum(~Flag_T.eval_r2))])

Flag_T.neb = s_max3000./s_max250<1;
disp(['sem neblina = ' num2str(sum(~Flag_T.neb)) '  com neb = ' num2str(sum(Flag_T.neb))])

Flag_T.FitGood = Flag_T.eval_p1  &  Flag_T.eval_p2 & Flag_T.eval_r2 & (~Flag_T.neb);
disp(['good geral = ' num2str(sum(Flag_T.FitGood)) '  bad geral = ' num2str(sum(~(Flag_T.FitGood)))])

return
clear Mflag_temp
%%

%% qualidade do ajuste                            

figure (10)
clf 
subplot(2,2,1); hold on; grid on; box on
edges = [0:0.05:1];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(rsquare(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1),'r-o')
[n1,x1] = histc(rsquare(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1),'b-o')
[n2,x2] = histc(rsquare,edges);
plot(bins,n2(1:end-1),'k-o')
xlabel('R2')
ylabel('contagens')
legend('bom','ruim','todos')
%%
subplot(2,2,2); hold on; grid on; box on
edges = [0:0.15:4.5];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(rmse(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1),'r-o')
[n1,x1] = histc(rmse(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1),'b-o')
[n2,x2] = histc(rmse,edges);
plot(bins,n2(1:end-1),'k-o')
xlabel('RMSE')
ylabel('contagens')
legend('bom','ruim','todos')
%%
subplot(2,2,3); hold on; grid on; box on
edges = [0:2:50];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(dfe(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1),'r-o')
[n1,x1] = histc(dfe(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1),'b-o')
[n2,x2] = histc(dfe,edges);
plot(bins,n2(1:end-1),'k-o')
xlabel('NGL')
ylabel('contagens')
legend('bom','ruim','todos')
%%
subplot(2,2,4); hold on; grid on; box on
edges = [0:10:150];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(sse(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1),'r-o')
[n1,x1] = histc(sse(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1),'b-o')
[n2,x2] = histc(sse,edges);
plot(bins,n2(1:end-1),'k-o')
xlabel('SSE')
ylabel('contagens')
legend('bom','ruim','todos')


%% P1 contra qualidade do ajuste                            

figure (20)
clf
subplot(2,2,1); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),rsquare(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p1(tmp),rsquare(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),rsquare(tmp),'bo')
plot(p1,rsquare,'k.')
xlabel('p1')
ylabel('R2')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
subplot(2,2,2); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),rmse(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p1(tmp),rmse(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),rmse(tmp),'bo')
plot(p1,rmse,'k.')
xlabel('p1')
ylabel('RMSE')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
subplot(2,2,3); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),dfe(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p1(tmp),dfe(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),dfe(tmp),'bo')
plot(p1,dfe,'k.')
xlabel('p1')
ylabel('NGL')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
subplot(2,2,4); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),sse(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p1(tmp),sse(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p1(tmp),sse(tmp),'bo')
plot(p1,sse,'k.')
xlabel('p1')
ylabel('SSE')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')

%% P2 contra qualidade do ajuste                            

figure (30)
clf
subplot(2,2,1); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),rsquare(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p2(tmp),rsquare(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),rsquare(tmp),'bo')
plot(p2,rsquare,'k.')
xlabel('p2')
ylabel('R2')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
subplot(2,2,2); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),rmse(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p2(tmp),rmse(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),rmse(tmp),'bo')
plot(p2,rmse,'k.')
xlabel('p2')
ylabel('RMSE')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
%aqui deu pra ver bem que as sem nuvem tem os menores coef lin
subplot(2,2,3); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),dfe(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p2(tmp),dfe(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),dfe(tmp),'bo')
plot(p2,dfe,'k.')
xlabel('p2')
ylabel('NGL')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
subplot(2,2,4); hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),sse(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(p2(tmp),sse(tmp),'go')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(p2(tmp),sse(tmp),'bo')
plot(p2,sse,'k.')
xlabel('p2')
ylabel('SSE')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')

%%
%subplot(2,2,4); hold on; grid on; box on
%tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
%plot(p2(tmp)./erro_p2(tmp),p1(tmp),'ro')
%tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
%plot(p2(tmp)./erro_p2(tmp),p1(tmp),'go')
%tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
%plot(p2(tmp)./erro_p2(tmp),p1(tmp),'bo')
%xlabel('p2 / Erro em p2')
%ylabel('valor de p1 (g/kg)')
%legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%

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
[n2,x2] = histc(p2,edges);
plot(bins,n2(1:end-1)/sum(n2),'k-o')
xlabel('p2')
ylabel('freq relativa (0-1)')
legend('bom','ruim','todos')
%%
subplot(2,2,2); hold on; grid on; box on
edges = [0:1:15];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(erro_p2(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(erro_p2(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
[n2,x2] = histc(erro_p2,edges);
plot(bins,n2(1:end-1)/sum(n2),'k-o')
xlabel('erro p2')
ylabel('freq relativa (0-1)')
legend('bom','ruim','todos')
%%
subplot(2,2,3); hold on; grid on; box on
edges = [0:25:1500];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(p1(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(p1(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
[n2,x2] = histc(p1,edges);
plot(bins,n2(1:end-1)/sum(n2),'k-o')
xlabel('p1')
ylabel('freq relativa (0-1)')
legend('bom','ruim','todos')
%%
subplot(2,2,4); hold on; grid on; box on
edges = [0:25:600];
bins = edges(1:end-1)+(edges(2)-edges(1))/2;
[n0,x0] = histc(erro_p1(Flag_T.FitGood == 0),edges);
plot(bins,n0(1:end-1)/sum(n0),'r-o')
[n1,x1] = histc(erro_p1(Flag_T.FitGood == 1),edges);
plot(bins,n1(1:end-1)/sum(n1),'b-o')
[n2,x2] = histc(erro_p1,edges);
plot(bins,n2(1:end-1)/sum(n2),'k-o')
xlabel('erro p1')
ylabel('freq relativa (0-1)')
legend('bom','ruim','todos')


%%
figure (325)
clf; hold on; grid on; box on
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 1;
plot(erro_p1(tmp)./p1(tmp),p1(tmp),'ro')
tmp =  Flag_T.cloud_5km == 0 &  Flag_T.cloud_2km == 0;
plot(erro_p1(tmp)./p1(tmp),p1(tmp),'ko')
tmp =  Flag_T.cloud_5km == 1 &  Flag_T.cloud_2km == 1;
plot(erro_p1(tmp)./p1(tmp),p1(tmp),'bo')
xlabel('Erro em p1 (0-1)')
ylabel('valor de p1 (g/kg)')
legend('entre 2 e 5 ','abaixo de 2','sem nuvem')
%%
figure (326)
clf; hold on; grid on; box on
mask_temp = jdi<=jdi(137);
mask_temp0 = jdi>=jdi(192);
eval_p1 = abs(erro_p1./p1)<0.2;
disp(['good p1 = ' num2str(sum(eval_p1)) '  bad p1 = ' num2str(sum(~eval_p1))])
eval_p2 = abs(p2./erro_p2)<3 ; 
disp(['good p2 = ' num2str(sum(eval_p2)) '  bad p2 = ' num2str(sum(~eval_p2))])
eval_r2 = rsquare>0.8;
disp(['good r2 = ' num2str(sum(eval_r2)) '  bad r2 = ' num2str(sum(~eval_r2))])
tmp = eval_p1  &  eval_p2 & eval_r2 ;
neb = s_max3000./s_max250<1;
disp(['sem neblina = ' num2str(sum(~neb)) '  com neb = ' num2str(sum(neb))])
disp(['good geral = ' num2str(sum(tmp & ~neb)) '  bad geral = ' num2str(sum(~(tmp & ~neb)))])
x0 = mean(jdi(mask_temp));
a = p1(mask_temp & tmp & ~neb);
b = jdi(mask_temp & tmp & ~neb)-x0;
[c,d] = fit(b',a','poly1')
x0 = mean(jdi(mask_temp0));
a = p1(mask_temp0 & tmp & ~neb);
x = jdi(mask_temp0 & tmp & ~neb)-x0;
[e,f] = fit(x',a','poly1')
plot(dti(tmp ),p1(tmp ),'o')
plot(dti((~tmp)),p1((~tmp)),'.r')
%neb = s_max3000./s_max250<1;
%plot(dti(neb ),p1(neb ),'sc')          
plot(dti(mask_temp & tmp & ~neb),c(b), '-')
plot(dti(mask_temp0 & tmp & ~neb),e(x), '-')
xlabel('Date')
ylabel('Calibration constants (g/kg)')
%legend('Good calibrations','Not good calibrations','Foggy days','2011 - 2012','2014 - 2016')
legend('Good calibrations','Rejected calibrations','2011 - 2012','2014 - 2016')
% constant = (-0.3 +-0.07)*(JD-JDmean) + (0.78+-0.01)
% constant = (0.13 +-0.04)*(JD-JDmean) + (0.86+-0.01)
saveas(figure(326),[figure_dir 'fitboas_calibracoes_.png' ])

%%
figure(3)
clf; hold on; box on; grid on
histogram(p1/1000,linspace(0, 1.5, 50))
histogram(p1(tmp)/1000,linspace(0, 1.5, 50))
xlabel('Calibration Constants kg/kg')
ylabel('Ocurrences')
legend(['Média =' num2str(mean(p1(tmp ))/1000) 'kg/kg' 'Desvpad =' num2str(std(p1(tmp)))])
title('Histogram only with the "good" calibrations')
%%
figure (4)
clf
histogram(p1/1000,linspace(0, 1.5, 50))
xlabel('Calibration Constants kg/kg')
ylabel('Ocurrences')
legend(['Mean =' num2str(mean(p1)/1000) 'kg/kg' 'Std =' num2str(std(p1(mask_temp)/1000))])
title('Histogram with all the calibrations')
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