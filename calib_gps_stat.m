close all
fclose all
clear all
%%
figure_dir = './prepro3/figures/statistic'
%%

prepro_dir = '../prepro2/';
flist = dir([prepro_dir 'Calib*.mat']);

jdi = NaN(1,length(flist));
p1 = jdi;
p2 = jdi;

for i = 1:length(flist)
    
    i
    
    clear Data_temp
    Data_temp = load([prepro_dir flist(i).name]);
    
    jdi(i) = Data_temp.jd(1);
    p1(i) = Data_temp.constante_cfit.p1;
    p2(i) = Data_temp.constante_cfit.p2;
    calib_pwv_sonde(i) = Data_temp.calib_pwv_sonde(1);
    calib_pwv_gps(i) = Data_temp.calib_pwv_gps(1);
    calib_temp_gps(i) = Data_temp.calib_temp_gps(1);
    calib_press_gps(i) = Data_temp.calib_press_gps(1);
    calib_gps_constante(i) = Data_temp.constante_gps(1);
    dfe(i) = Data_temp.constante_gof.dfe;
    rsquare(i) = Data_temp.constante_gof.rsquare;
    rmse(i) = Data_temp.constante_gof.rmse;
    sse(i) = Data_temp.constante_gof.sse;
    X = confint(Data_temp.constante_cfit);
    erro_p1(i) = (X(2,1)-X(1,1))/4;
    erro_p2(i) = (X(2,2)-X(1,2))/4;
    s_max250(i) = Data_temp.s_max250(1);
    s_max3000(i) = Data_temp.s_max3000(1);
    
end

%%
dti = datetime(datestr(jdi));


%% 
figure (1)
clf 
hold on 
plot(dti,p1,'o')
plot(dti,calib_gps_constante*1000,'o')
legend('Cte calib sonda', 'Cte gps')
xlabel('Data')
ylabel('constante de calibração (g/kg)')
title('Embrapa')
hold off
grid on

%%

dif_constante = -1*(calib_gps_constante*1000 - p1);
dif_pwv = -1*(calib_pwv_gps*10 - calib_pwv_sonde);

figure (2)
clf 
hold on 
plot(dif_pwv,dif_constante,'o')
ylabel('Constante sonda - Constante gps [g/kg]')
xlabel('PWV sonda - PWV gps [mm]')
title('Embrapa')
hold off
grid on

%%

mean(p1)
std(p1)
mean(calib_gps_constante)
std(calib_gps_constante)

%%

figure (3)
clf 
hold on 
plot(dti,calib_pwv_gps*10,'o')
plot(dti,calib_pwv_sonde,'o')
legend('pwv gps', 'pwv sonda')
xlabel('Data')
ylabel('PWV (mm)')
title('Embrapa')
hold off
grid on
%%

x0 = datenum(2011,07,17,0,0,0);
x = jdi - x0;
y = p1;
%%
z = calib_gps_constante*1000;
%%
calib_cfit = fit(x',y','poly1')
calib_cfit_gps = fit(x',z','poly1')
figure (8)
clf
hold on
plot(dti,p1,'o')
plot(dti,calib_gps_constante*1000,'o')
plot(x+x0,calib_cfit(x),'b')
plot(x+x0,calib_cfit_gps(x))
hold off
legend('cte sonda','cte gps','sonda(x) = -0,9*x + 890','gps(x) = 0,3*x+856')
title('Calibrações Embrapa')
xlabel('data da calibração')
ylabel('constante de calibração [g/kg]')
grid on

return

%%
prepro_dir = 'C:\Users\guido\Google Drive\prepro3\';
flist = dir([prepro_dir 'Calib*.mat']);

jdi_reserva = NaN(1,length(flist));
p1_reserva = jdi_reserva;
p2_reserva = jdi_reserva;

for i = 1:length(flist)
    
    i
    
    clear Data_temp
    Data_temp = load([prepro_dir flist(i).name]);
    
    jdi_reserva(i) = Data_temp.jd(1);
    p1_reserva(i) = Data_temp.constante_cfit.p1;
    p2_reserva(i) = Data_temp.constante_cfit.p2;
    calib_pwv_sonde_reserva(i) = Data_temp.calib_pwv_sonde(1);
    calib_pwv_gps_reserva(i) = Data_temp.calib_pwv_gps(1);
    calib_temp_gps_reserva(i) = Data_temp.calib_temp_gps(1);
    calib_press_gps_reserva(i) = Data_temp.calib_press_gps(1);
    calib_gps_constante_reserva(i) = Data_temp.constante_gps(1);
     dfe_reserva(i) = Data_temp.constante_gof.dfe;
    rsquare_reserva(i) = Data_temp.constante_gof.rsquare;
    rmse_reserva(i) = Data_temp.constante_gof.rmse;
    sse_reserva(i) = Data_temp.constante_gof.sse;
    X = confint(Data_temp.constante_cfit);
    erro_p1_reserva(i) = (X(2,1)-X(1,1))/4;
    erro_p2_reserva(i) = (X(2,2)-X(1,2))/4;
    s_max250_reserva(i) = Data_temp.s_max250(1);
    s_max3000_reserva(i) = Data_temp.s_max3000(1);
    
end

dti_reserva = datetime(datestr(jdi_reserva));

%%
%%
%aqui elimina os dados desse plot menores do que 1.5
plot([dti;dti_reserva],[s_max3000'./s_max250' ; s_max3000_reserva'./s_max250_reserva'])
 %%
mask = abs(erro_p1./p1)<0.2 &  abs(p2./erro_p2)<2 & rsquare>0.8;
neb = s_max3000'./s_max250'<1.5;
mask1 = abs(erro_p1_reserva./p1_reserva)<0.2 &  abs(p2_reserva./erro_p2_reserva)<2 & rsquare_reserva>0.8;
neb1 = s_max3000_reserva'./s_max250_reserva'<1.5;
b = horzcat(p1,p1_reserva);
a = horzcat(dti',dti_reserva');
c = horzcat(mask,mask1);
d = horzcat(neb',neb1');
figure (4)
clf
hold on 
plot(a(c & ~d),b(c & ~d),'ok')
plot(dti_reserva(mask1 & ~neb1'),calib_gps_constante_reserva(mask1 & ~neb1')*1000,'or')
plot(dti(mask & ~neb'),calib_gps_constante(mask & ~neb')*1000,'ob')

plot(a(d),b(d),'sg')
plot(a(~c),b(~c),'.k')

plot(dti_reserva(neb1),calib_gps_constante_reserva(neb1)*1000,'sg')
plot(dti(neb),calib_gps_constante(neb)*1000,'sg')

plot(dti_reserva(~mask1),calib_gps_constante_reserva(~mask1)*1000,'.r')
plot(dti(~mask),calib_gps_constante(~mask)*1000,'.b')

legend('Cte sonda', 'Cte gps reserva','Cte gps Embrapa','Dias com neblina','Calibrações ruins')
xlabel('Data')
title('Reserva')
ylabel('constante de calibração [g/kg]')
hold off
grid on

%%

mean(p1_reserva)
std(p1_reserva)
mean(calib_gps_constante_reserva)
std(calib_gps_constante_reserva)

%%

%%

mask = abs(erro_p1./p1)<0.2 &  abs(p2./erro_p2)<2 & rsquare>0.8;
neb = s_max3000'./s_max250'<1.5;
dif_constante = -1*(calib_gps_constante(mask & ~neb')*1000 - p1(mask & ~neb'));
dif_pwv = -1*(calib_pwv_gps(mask & ~neb')*10 - calib_pwv_sonde(mask & ~neb'));

mask = abs(erro_p1_reserva./p1_reserva)<0.2 &  abs(p2_reserva./erro_p2_reserva)<2 & rsquare_reserva>0.8;
neb = s_max3000_reserva'./s_max250_reserva'<1.5;
dif_constante_reserva = -1*(calib_gps_constante_reserva(mask & ~neb')*1000 - p1_reserva(mask & ~neb'));
dif_pwv_reserva = -1*(calib_pwv_gps_reserva(mask & ~neb')*10 - calib_pwv_sonde_reserva(mask & ~neb'));

a = horzcat(dif_constante,dif_constante_reserva);
b = horzcat(dif_pwv,dif_pwv_reserva);
c = fit(b',a','poly1')
figure (13)
clf 
hold on 
plot(b,a,'o')
plot(b,c(b))
ylabel('Constante sonda - Constante gps [g/kg]')
xlabel('PWV sonda - PWV gps [mm]')
title('Relação diferença pwv e cte calibração')
hold off
grid on
%%
figure (12)
clf 
hold on 
plot(x,y*10,'o')
plot(x,z,'o')
legend('pwv gps', 'pwv sonda')
xlabel('Data')
ylabel('PWV (mm)')
title('COmparação de PWV GNSS e Sonda')
hold off
grid on
