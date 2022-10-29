%%
% close all
% fclose all
clear all

% change display for all plots
set(0, 'DefaultAxesFontSize', 16, 'DefaultLineLineWidth', 2)
%%
%addpath('../IC Henrique/')
%addpath('../data/')
addpath('../../matlab/sc/')
%%
load('sonde.mat')

jd_sonde_utc = [sonde(:).jd];
jd_sonde = jd_sonde_utc - datenum(00,00,00,04,00,00);   % hora local
pwv_sonde = [sonde(:).pwat];

%%
T = readtable('reserva_2011.txt','Delimiter',',');

jd_gps_utc = T.dayofyear - 1 + datenum(2011,1,1); %encontramos erro pois estavamos somando um dia a mais
jd_gps = jd_gps_utc - datenum(00,00,00,04,00,00);   % hora local
%%
% abrir varios arquivos
%     datadir = 'data';
%datadir = 'F:\Dados lidars\lidar\data' 
% datadir = 'D:\Dados\lidar\data';
%datadir = 'D:\Dados\lidar\data'
datadir = '/Users/hbarbosa/DATA/lidar/data'
%%

figure_dir = './hmjb/figures/';
prepro_dir = '../prepro3/';

%%


%% ler planilha em format CSV
%% deve ter pelo menos colunas para altura maxima e minima

%% Y M D H Hmin Hmax FitBom? NumvenBaixa(1/0) Media Alta 
%%
%%

%%
% iok = 0
clear calib_pwv_sonde calib_pwv_gps
for yu = datenum(2011,08,27,20,00,00):datenum(2011,08,27,20,00,00)
    
    clearvars -except yu sonde jd_sonde jd_sonde_utc datadir figure_dir prepro_dir calib_pwv_sonde calib_pwv_gps jd_gps T calib_temp_gps calib_press_gps calib_sonde_cfit calib_gps_constante data_das_calib
    
    %
    
    
    [sonde_dist_dias, I_sonde]= min(abs(jd_sonde - yu));
    
    
    if sonde_dist_dias > 0
        
        continue
        
    end
    %
    
    
%     jdi=datenum(2014,09,14,20,00,00);       % data inicial
%     jdf= jdi + 1/24;                        % data final
     
   
    jdi = jd_sonde(I_sonde)- 30/(60*24);
    jdf = jd_sonde(I_sonde)+ 30/(60*24);
    
    
    % dbin_MAO = [10];
    dbin_MAO = [9 -0 9 -1 -1];
    dtime_MAO =  0.004;
    [nfile, head, phy] = profile_read_dates(datadir, jdi, jdf, dbin_MAO, dtime_MAO);
    
    
    
    if nfile == 0
        continue
    end
    
    pwv_gps = T.pwv(jd_gps>=jdi & jd_gps<=jdf); pwv_gps(pwv_gps<0)=NaN;
    press_gps = T.press(jd_gps>=jdi & jd_gps<=jdf); press_gps(press_gps<0)=NaN;
    temp_gps = T.temp(jd_gps>=jdi & jd_gps<=jdf); temp_gps(temp_gps<0)=NaN;
    if numel(pwv_gps)== 0 || isnan(nanmean(pwv_gps))  %aqui não é isnan, o que é?
        continue
    end
%     iok = iok+1
    calib_pwv_sonde = sonde(I_sonde).pwat;
    calib_pwv_gps = nanmean(pwv_gps);
    calib_temp_gps = nanmean(temp_gps);
    calib_press_gps = nanmean(press_gps);
    
    
    %disp(['================================ valor de ' num2str(iok)])
    
    
   
    %%
    % alturas
    range(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);
    
    %%
    % altitude
    altitude = range.*cosd(5) + head(1).alt;
    
    %%
    % horarios dos arquivos
    for i=1:numel(head)
        jd(i)=head(i).jdi;
    end
    
    %% Subtrair BG do sinal bruto
    
    % mask_bg = size(phy(ch).data,1)-1000:size(phy(ch).data,1);
    mask_bg = range > range(end-500) & range < range(end);
    
    for ch=1:head(1).nch
        phy(ch).bgcorr = phy(ch).data.*NaN;
        for jj = 1:length(jd)
            phy(ch).bgcorr(:,jj) = phy(ch).data(:,jj) - mean(phy(ch).data(mask_bg,jj));
        end
    end
    
    
    %%
    aux = phy(1).bgcorr.*(range*ones(1,length(jd))).^2;
    
    r_aux = range<20000;
    figure(11)
    clf
    %set(gcf,'units','normalized','outerposition',[0 0 0.8 1])
    
    [MAT_aux,jd_aux] = crop_tempo_0(aux,jd,1,1); %%% Essa função esta errada (bem porca). Arrumar depois
    gplot2(MAT_aux(r_aux,:),[0:1:10].*1e6, jd_aux,range(r_aux)/1.e3)
    % gplot2(aux,[0:0.1:1.5].*1e7, jd,alt/1.e3)
    %title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    title(['Sinal do lidar ao longo de 1 hora' datestr(jd(1)) ' to ' datestr(jd(end))] )
    ylabel('Altura [km]')
    xlabel('Horário [LT]')
    ylim([0 15])
    datetick('x','keeplimits')
    
    hc = colorbar;
    ylabel(hc,'Sinal corrigido pela distância [U.A.]','fontsize',14)
    
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_RCS_2D_.png' ])
    
    % pause
    %%
    
    %%
    %tirando o bg
    for i=1:head(1).nch
        PcomBG(:,i) = mean(phy(i).data(:,:),2);
        PcomBG_std(:,i) = std(phy(i).data(:,:),1,2);
        PcomBG_stdm(:,i) = std(phy(i).data(:,:),1,2)./size(phy(i).data(:,:),2).^0.5;
        
        bg(i) = mean(PcomBG(mask_bg,i));
        desvio(i) = std(PcomBG(mask_bg,i), 1);
        bg_stdm(i) = std(PcomBG(mask_bg,i), 1)./size(PcomBG(mask_bg,i),1).^0.5;
        %hmjb
        % So faz a correcao de BG se for AN
        % PC nao precisa porque eh de noite, nao tem ruido
        if (head(1).ch(i).photons==0)
            P(:,i) = PcomBG(:,i) - bg(i);
        else
            P(:,i) = PcomBG(:,i);
        end
        
        Pr2(:,i) = P(:,i).*range.*range;
    end
    
    
    %% achar a regiao onde vale o photocount e o analogic
    ANmask = PcomBG(:,3) > 2.06 & PcomBG(:,3) < 50 & range > 100; % em mV do AN
    sum(ANmask)
    PCmask = PcomBG(:,4) < 12.; % em MHz do PC
    sum(PCmask)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_glue = P(ANmask & PCmask,3);
    Y_glue = P(ANmask & PCmask,4);
    cfun_glue = fit(X_glue,Y_glue,'poly1')
    
    
    Signal_glued = cfun_glue(P(:,3));
    Signal_glued(PCmask) = P(PCmask,4);
    %sinal máximo perto do chão 
    a = floor(250/7.5);
    b = floor(3000/7.5);
    s_max250 = max(Signal_glued(1:a));
    s_max3000 = max(Signal_glued(a:b));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(100)
    clf
    hold on
    plot(X_glue,Y_glue,'.r')
    plot(cfun_glue,'b')
    hold off
    title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    xlabel('AN Signal [mV]')
    ylabel('PC Signal [MHz]')
    grid on
    
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_fit_glue_.png' ])
    
    
    %%
    
    figure(110)
    clf
    hold on
    plot(range/1000,P(:,3),range/1000,P(:,4),range/1000,Signal_glued,'k')
    plot(range(ANmask & PCmask)/1000,P(ANmask & PCmask,3),'o',range(ANmask & PCmask)/1000,P(ANmask & PCmask,4),'o')
    plot(range/1000,P(:,5),'c')
   %plot([0,15],[1,1].*bg(4),'-r')
    plot([0,15],[1,1].*(3*desvio(4)),'--r')
    plot([0,15],[1,1].*(10*desvio(4)),'--r')
    %plot([0,15],[1,1].*bg(5),'-b')
    plot([0,15],[1,1].*(3*desvio(5)),'--b')
    plot([0,15],[1,1].*(10*desvio(5)),'--b')
    hold off
    legend('AN','PC','Glued','Fit region AN','Fit region PC')
    grid on
    set(gca,'yscale','log')
    xlim([0 15])
    title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    xlabel('Altura [m]')
    ylabel(['Sinal PC ou AN [MHz ou mv]'])
    
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_Signal_glue_.png' ])
    
    %%
    
% %     figure(110000)
% %     clf
% %     hold on
% %     plot(range/1000,P(:,3),range/1000,P(:,4),range/1000,Signal_glued,'k')
% %     plot(range(ANmask & PCmask)/1000,P(ANmask & PCmask,3),'o',range(ANmask & PCmask)/1000,P(ANmask & PCmask,4),'o')
% %     plot(range/1000,smooth(P(:,4),11) + PcomBG_stdm(:,4),'r')
% %     hold off
% %     legend('AN','PC','Glued','Fit region AN','Fit region PC')
% %     grid on
% %     set(gca,'yscale','log')
% %     xlim([0 15])
% %     title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
% %     xlabel('Altura [m]')
% %     ylabel(['Sinal PC ou AN [MHz ou mv]'])
% %     
% %     
% %     figure(1100000)
% %     clf
% %     plot(range/1000,smooth(P(:,4),11)./PcomBG_stdm(:,4),'r')
% %     grid on
    %%
    % lendo a sondagem
    
    altitude_snd = sonde(I_sonde).alt;  % [m]
    w_snd = sonde(I_sonde).MIXR;        % [g/kg];
    
    mask = altitude_snd >= 0 & w_snd >= 0;
    
    altitude_snd = altitude_snd(mask);
    w_snd = w_snd(mask);    
    
    %%
    
    [constante_cfit,constante_gof,constante_output] = ...
        calib_function_v0_gps(altitude_snd, w_snd, ...
                              altitude, Signal_glued,P(:,5),...
                              jd(1),jd(end),desvio(4),desvio(5));
    %[constante_gps] = calib_function_gps(press_gps,temp_gps,pwv_gps,altitude,Signal_glued,P(:,5),jd(1),jd(end),desvio(4),desvio(5),head(1).alt);

    % [constante_cfit] = calib_function_v0(sonda(:,2),sonda(:,6),altitude,P(:,4),P(:,5));
    
%     %%
%     calib_sonde_cfit{iok} = constante_cfit;
%     calib_gps_constante(iok) = constante_gps;
%     data_das_calib(iok) = jdi
    %%
    fname=[datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_lidar_uncalib_h2o_' ]
    saveas(figure(730),[figure_dir fname '.png' ])
    set(figure(730),'paperunits','points','papersize',[560 420])
    print(figure(730),[figure_dir fname '.pdf' ],'-dpdf','-bestfit')

    fname=[datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_lidar_signal_'];
    saveas(figure(1538),[figure_dir fname '.png' ])
    set(figure(1538),'paperunits','points','papersize',[560 420])
    print(figure(1538),[figure_dir fname '.pdf' ],'-dpdf','-bestfit')

    saveas(figure(1539),[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_MXR_sonde_mask_.png' ])

    fname=[datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_fit_calib_'];
    saveas(figure(1540),[figure_dir fname '.png' ])
    set(figure(1540),'paperunits','points','papersize',[560 420])
    print(figure(1540),[figure_dir fname '.pdf' ],'-dpdf','-bestfit')
        
    saveas(figure(1541),[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_residuos_fit_calib.png' ])
   
    fname=[datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_Signal_calibrado' ]
    saveas(figure(1542),[figure_dir fname '.png' ])
    set(figure(1542),'paperunits','points','papersize',[560 420])
    print(figure(1542),[figure_dir fname '.pdf' ],'-dpdf','-bestfit')

    %%
    
    %fname_save = [prepro_dir 'Calib_' datestr(yu,'yyyy_mm_dd_HH_MM_SS') '.mat' ];
    %save(fname_save,'constante_cfit','constante_gof','constante_output','altitude_snd','w_snd',...
    %    'range','altitude','Signal_glued','P',...
    %    'I_sonde','jd','calib_pwv_sonde','calib_pwv_gps','calib_temp_gps','calib_press_gps',...
    %    'constante_gps','jdi','s_max250','s_max3000')
    
%     pause
    
    
end
%     save([prepro_dir 'compara.mat'],'calib_pwv_sonde',...
%     'calib_pwv_gps','calib_temp_gps','calib_press_gps',...
%     'calib_sonde_cfit','calib_gps_constante','data_das_calib')
