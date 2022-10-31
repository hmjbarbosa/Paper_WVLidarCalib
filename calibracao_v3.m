%%
% close all
% fclose all
clear all

% change display for all plots
set(0, 'DefaultAxesFontSize', 12, 'DefaultLineLineWidth', 2)
%%

addpath('./libs/')
addpath('./sc/')

%% Read sonde data 
load('sonde.mat')

jd_sonde_utc = [sonde(:).jd];
jd_sonde = jd_sonde_utc - datenum(00,00,00,04,00,00);   % local time

%% GPS data from Reserva Duke

% porque os dados da Duke ao inves da Embrapa? ta errado isso...

T = readtable('reserva_2011.txt','Delimiter',',');
jd_gps_utc = T.dayofyear - 1 + datenum(2011,1,1);   % old bug: we were not subtracting -1
jd_gps = jd_gps_utc - datenum(00,00,00,04,00,00);   % local time

%% Paths for data processing

% datadir
datadir='/Users/hbarbosa/DATA/lidar/data/';

figure_dir = './hmjb/figures/';
prepro_dir = './prepro4/';

%% Loop over all dates to be processed

iok = 0;
clear calib_pwv_sonde calib_pwv_gps calib_temp_gps calib_press_gps ...
    calib_sonde_cfit calib_gps_constante data_das_calib;
for yu = datenum(2011,08,23,20,00,00):datenum(2011,08,23,20,00,00)
    
    clearvars -except yu sonde jd_sonde jd_sonde_utc datadir figure_dir prepro_dir ...
        calib_pwv_sonde calib_pwv_gps jd_gps T calib_temp_gps calib_press_gps ...
        calib_sonde_cfit calib_gps_constante data_das_calib iok
    
    % --------------------------------------------------------------------------
    % RADIO SONDE
    % --------------------------------------------------------------------------

    % find closest sounding (in time)
    [sonde_dist_dias, I_sonde]= min(abs(jd_sonde - yu));
    
    % skip if more than 1 day
    if sonde_dist_dias > 0
        continue
    end
    
    % select lidar data around time of sounding
    jdi = jd_sonde(I_sonde)- 30/(60*24);
    jdf = jd_sonde(I_sonde)+ 30/(60*24);

    altitude_sonde = sonde(I_sonde).alt;  % [m]
    wvmix_sonde = sonde(I_sonde).MIXR;    % [g/kg];
    
    mask = altitude_sonde >= 0 & wvmix_sonde >= 0;
    
    altitude_sonde = altitude_sonde(mask);
    wvmix_sonde = wvmix_sonde(mask);    
    
    % --------------------------------------------------------------------------
    % GPS
    % --------------------------------------------------------------------------

    % crop GPS data around time of sounding
    pwv_gps = T.pwv(jd_gps>=jdi & jd_gps<=jdf); 
    pwv_gps(pwv_gps<0)=NaN;
    
    press_gps = T.press(jd_gps>=jdi & jd_gps<=jdf); 
    press_gps(press_gps<0)=NaN;
    
    temp_gps = T.temp(jd_gps>=jdi & jd_gps<=jdf); 
    temp_gps(temp_gps<0)=NaN;
    
    % skip if there is no GPS data
    if numel(pwv_gps)== 0 || isnan(nanmean(pwv_gps))  %aqui não é isnan, o que é?
        continue
    end

    iok=iok+1;
    calib_pwv_sonde(iok) = sonde(I_sonde).pwat;
    calib_pwv_gps(iok) = nanmean(pwv_gps);
    calib_temp_gps(iok) = nanmean(temp_gps);
    calib_press_gps(iok) = nanmean(press_gps);

    % --------------------------------------------------------------------------
    % LIDAR
    % --------------------------------------------------------------------------
    
    % old: fixed dbin=10 for all channels
    % dbin_MAO = [10];
    % new: 
    dbin_MAO = [9 -0 9 -1 -1];
    dtime_MAO =  0.004;
    
    % read lidar data around time of sounding
    [nfile, head, phy] = profile_read_dates(datadir, jdi, jdf, dbin_MAO, dtime_MAO);
    
    % skip if there are no lidar files
    if nfile == 0
        continue
    end
        
    % range (m)
    range(:,1)=head(1).ch(1).binw*(1:head(1).ch(1).ndata);
    
    % altitude a.s.l. (m) 
    % note: lidar is tiled 5deg to avoid specular reflection on ice crystals
    altitude = range.*cosd(5) + head(1).alt;
    
    % time (julian date)
    for i=1:numel(head)
        jd(i)=head(i).jdi;
    end
    
    % Subtrair BG do sinal bruto
    
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
    
    %% BG removal 

    for i=1:head(1).nch
        % mean P, std P, std of the mean P
        P_withBG(:,i) = mean(phy(i).data(:,:),2);
        P_withBG_std(:,i) = std(phy(i).data(:,:),1,2);
        P_withBG_stdm(:,i) = std(phy(i).data(:,:),1,2)./size(phy(i).data(:,:),2).^0.5;
        
        % mean BG, std BG, std of the mean BG
        bg(i) = mean(P_withBG(mask_bg,i));
        bg_std(i) = std(P_withBG(mask_bg,i), 1);
        bg_stdm(i) = std(P_withBG(mask_bg,i), 1)./size(P_withBG(mask_bg,i),1).^0.5;

        % correct for BG if channel is AN
        % note: WV calib at night only, and PC has no BG
        if (head(1).ch(i).photons==0)
            P(:,i) = P_withBG(:,i) - bg(i);
        else
            P(:,i) = P_withBG(:,i);
        end
        
        % RCS
        Pr2(:,i) = P(:,i).*range.*range;
    end
    
    %% Gluing
    
    % find region where PC and AN are valid
    ANmask = P_withBG(:,3) > 2.06 & P_withBG(:,3) < 50; 
    PCmask = P_withBG(:,4) < 12.; 
    %ANmask = (P_withBG(:,3) > bg(3)*1.01) & (P_withBG(:,3) < head(1).ch(3).discr*1e3/2.); 
    %PCmask = P_withBG(:,4) < 45.; 
    Gmask = ANmask & PCmask & range > 100.; 
    
    X_glue = P(Gmask,3);
    Y_glue = P(Gmask,4);
    cfun_glue = fit(X_glue,Y_glue,'poly1')
    
    Signal_glued = cfun_glue(P(:,3));
    Signal_glued(PCmask) = P(PCmask,4);

    % Get signal near ground x PBL
    % to exclude fog near ground
    a = floor(250/7.5);
    b = floor(3000/7.5);
    s_max250 = max(Signal_glued(1:a));
    s_max3000 = max(Signal_glued(a:b));
    
    figure(100)
    clf; hold on; grid on; box on
    plot(X_glue,Y_glue,'.r')
    plot(cfun_glue,'b')
    fitresult = evalc('cfun_glue');
    text(0.45,0.2, fitresult(21:end-1),'units','normalized', ...
         'Interpreter','none','backgroundcolor','white','edgecolor','black')
    title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    xlabel('AN Signal [mV]')
    ylabel('PC Signal [MHz]')
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_fit_glue_.png' ])
    
    figure(101)
    clf; hold on; grid on; box on
    plot(X_glue,Y_glue - cfun_glue(X_glue),'.r')
    xlabel('AN Signal [mV]')
    ylabel('Residue [MHz]')
    set(gca,'xscale','log')
    
    figure(110)
    clf; hold on; grid on; box on
    % nitrogen
    plot(range/1e3, P(:,3), 'b')      % AN 387
    plot(range/1e3, P(:,4), 'r')      % PC 387 
    plot(range/1e3, Signal_glued,'k') % Glue 387
    plot(range(Gmask)/1e3, P(Gmask,3),'o') % region used in the gluing
    plot(range(Gmask)/1e3, P(Gmask,4),'o') % region used in the gluing
    % vapor
    plot(range/1e3, P(:,5),'c') % PC 408
    % bg for channels 387 and 408
    %plot([0,30],[1,1].*bg(4),'-r')
    plot([0,30],[1,1].*(3*bg_std(4)),'--r')
    plot([0,30],[1,1].*(10*bg_std(4)),'--r')
    %plot([0,30],[1,1].*bg(5),'-b')
    plot([0,30],[1,1].*(3*bg_std(5)),'--c')
    plot([0,30],[1,1].*(10*bg_std(5)),'--c')
    legend('AN387','PC387','Glued387','Fit region AN','Fit region PC', 'PC408')

    set(gca,'yscale','log')
    xlim([0 30])
    title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    xlabel('Range [m]')
    ylabel(['Signal PC or AN [MHz or mv]'])
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_Signal_glue_.png' ])
    
    % --------------------------------------------------------------------------
    % CALIBRATION
    % --------------------------------------------------------------------------
    
    [constante_cfit,constante_gof,constante_output] = ...
        calib_function_v1(altitude_sonde, wvmix_sonde, ...
                          altitude, Signal_glued, P(:,5), ...
                          jd(1), jd(end), bg_std(4), bg_std(5));
    
    [constante_gps] = calib_function_gps(press_gps, temp_gps, pwv_gps, ...
                                         altitude, Signal_glued, P(:,5), ...
                                         jd(1), jd(end), bg_std(4), bg_std(5), head(1).alt);

    %%
    calib_sonde_cfit{iok} = constante_cfit;
    calib_gps_constante(iok) = constante_gps;
    data_das_calib(iok) = jdi
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
    
    fname_save = [prepro_dir 'Calib_' datestr(yu,'yyyy_mm_dd_HH_MM_SS') '.mat' ];
    save(fname_save,'constante_cfit','constante_gof','constante_output',...
         'jd','jdi','altitude','range','Signal_glued','P',...
         'altitude_sonde','wvmix_sonde','I_sonde',...
         'calib_pwv_sonde','calib_pwv_gps','calib_temp_gps','calib_press_gps','constante_gps',...
         's_max250','s_max3000')
    
%     pause
    
    
end
save([prepro_dir 'compara.mat'],'calib_pwv_sonde',...
     'calib_pwv_gps','calib_temp_gps','calib_press_gps',...
     'calib_sonde_cfit','calib_gps_constante','data_das_calib')

%fim