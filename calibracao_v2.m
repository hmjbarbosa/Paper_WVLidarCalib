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

%% Paths for data processing

% datadir
datadir='/Users/hbarbosa/DATA/lidar/data/';

figure_dir = './hmjb/figures/';
prepro_dir = './prepro4/';

%% read DC 

% old: fixed dbin=10 for all channels
% dbin_MAO = [10];
% new: 
dbin_MAO = [9 -0 9 -1 -1];
dtime_MAO =  0.004;

DCdir='2012_11_29_telescope/Dark_Current_19h34_night/';
[nfileDC, headDC, phyDC] = profile_read_dir([datadir DCdir], dbin_MAO, dtime_MAO);
    
%% Loop over all dates to be processed

for yu = datenum(2011,08,27,20,00,00):datenum(2011,08,27,20,00,00)
    
    clearvars -except yu sonde jd_sonde jd_sonde_utc datadir figure_dir prepro_dir ...
        dbin_MAO dtime_MAO nfileDC headDC phyDC
    
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
    jdi = jd_sonde(I_sonde) - 0.5/24;
    jdf = jdi + 1./24;
    
    altitude_sonde = sonde(I_sonde).alt;  % [m]
    wvmix_sonde = sonde(I_sonde).MIXR;    % [g/kg];
    
    mask = altitude_sonde >= 0 & wvmix_sonde >= 0;
    
    altitude_sonde = altitude_sonde(mask);
    wvmix_sonde = wvmix_sonde(mask);    
   
    figure(50)
    clf; hold on; grid on; box on;
    plot(sonde(I_sonde).RH, sonde(I_sonde).alt/1e3)
    xlabel('Relative Humidity [%]')
    ylabel('Altitude asl [km]')

    
    % --------------------------------------------------------------------------
    % LIDAR
    % --------------------------------------------------------------------------
        
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
    
    %% BG removal 
    
    % mask_bg = size(phy(ch).data,1)-1000:size(phy(ch).data,1);
    mask_bg = range > range(end-1000) & range < range(end);
        
    for ch=1:head(1).nch
        phy(ch).bgcorr = phy(ch).data.*NaN;
        % crop DC to current length
        DC(:,ch) = mean(phyDC(ch).data(1:head(1).ch(1).ndata, :), 2);
        % make sure DC has zero average
        DC(:,ch) = DC(:,ch) - mean(DC(mask_bg,ch));
        for jj = 1:length(jd)
            phy(ch).bg(jj) = mean(phy(ch).data(mask_bg,jj));
            if (head(1).ch(ch).photons==0)
                phy(ch).bgcorr(:,jj) = phy(ch).data(:,jj) - phy(ch).bg(jj) - DC(:,ch);
            else
                phy(ch).bgcorr(:,jj) = phy(ch).data(:,jj);
            end
        end
    end
    
    R2mat = (range*ones(1,length(jd))).^2;
    
    % PLOT Signal x Time
    for ch=1:head(1).nch
        figure(10+ch); clf; 
        hold on; grid on; box on
        tmp = log10(phy(ch).bgcorr.*R2mat);
        tmp(phy(ch).bgcorr<=0) = nan;
        [h, bar]=gplot2(tmp,[], jd, range/1.e3);
        %set(gca,'colorscale','log') 
        if (head(1).ch(ch).photons==1) 
            tag=[num2str(head(1).ch(ch).wlen) '_PC'];
        else
            tag=[num2str(head(1).ch(ch).wlen) '_AN']; 
        end
        title([datestr(jd(1),'yyyy/mm/dd') ' ' tag],'interpreter','none')
        ylabel('Altitude [km]')
        xlabel('Time [LT]')
        ylim([0 15])
        datetick('x','keeplimits')    
        ylabel(bar,'Log10 RCS [A.U.]','fontsize',14)
        %saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_RCS_2D_' tag '.png' ])
    end
    
    for ch=1:head(1).nch
        % mean P, std P, std of the mean P
        P_withBG(:,ch) = mean(phy(ch).data(:,:),2);
        P_withBG_std(:,ch) = std(phy(ch).data(:,:),1,2);
        
        % mean BG, std BG, std of the mean BG
        bg(ch) = mean(P_withBG(mask_bg,ch));
        bg_std(ch) = std(P_withBG(mask_bg,ch), 1);

        % correct for BG if channel is AN
        % note: WV calib at night only, and PC has no BG
        if (head(1).ch(ch).photons==0)
            P(:,ch) = P_withBG(:,ch) - bg(ch) - DC(:,ch);
            SNR(:,ch) = P(:,ch)./P_withBG_std(:,ch);
        else
            P(:,ch) = P_withBG(:,ch);
            % get total number of photons
            % poisson uncertainty is sqrt(N), hence SNR is N/sqrt(N) = sqrt(N)
            SNR(:,ch) = sqrt(P(:,ch)* head(1).nshoots / 20);
        end
        
        % RCS
        Pr2(:,ch) = P(:,ch).*range.*range;
    end
    
    % PLOT MEAN SIGNAL / Signal-to-noise ratio
    for ch=1:head(1).nch
        figure(20+ch); clf; 
        subplot(1,2,1)
        hold on; grid on; box on
        plot(range.*range.*(P(:,ch)), range/1.e3, 'b-')
        if (head(1).ch(ch).photons==1) 
            tag=[num2str(head(1).ch(ch).wlen) '_PC'];
        else
            tag=[num2str(head(1).ch(ch).wlen) '_AN']; 
        end
        title([datestr(jd(1),'yyyy/mm/dd') ' ' tag],'interpreter','none')
        ylabel('Altitude [km]')
        xlabel('RCS')
        set(gca,'xscale','log')
        ylim([0 30])
        
        subplot(1,2,2)
        plot(SNR(:,ch), range/1.e3, 'r-')
        hold on; grid on; box on
        xlabel('S/N ratio')
        title(['SNR(5km)=' num2str(mean(SNR(650:680,ch)))])
        %saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_RCS_1D_' tag '.png' ])
    end
        
    %% Gluing
    
    % find region where PC and AN are valid
    %ANmask = P_withBG(:,3) > 2.005 & P_withBG(:,3) < 50; 
    %PCmask = P_withBG(:,4) < 12.; 
    maxAN = head(1).ch(3).discr*1e3; % mV
    %resolAN = maxAN/2^head(1).ch(3).bits; % mV
    ANmask = (P_withBG(:,3) > bg(3) + bg_std(3)*10 ) & (P_withBG(:,3) < maxAN/2.); 
    PCmask = P_withBG(:,4) < 45.; 
    Gmask = ANmask & PCmask & range > 100.; 
    
    X_glue = P(Gmask,3);
    Y_glue = P(Gmask,4);
    cfun_glue = fit(X_glue,Y_glue,'poly1')
    
    Signal_glued = cfun_glue(P(:,3));
    Signal_glued(PCmask) = P(PCmask,4);

    % snr from glue
    SNR_glue = sqrt(Signal_glued* head(1).nshoots / 20);

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
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_fit_glue.png' ])
    
    figure(101)
    clf; hold on; grid on; box on
    plot(X_glue,Y_glue - cfun_glue(X_glue),'or')
    xlabel('AN Signal [mV]')
    ylabel('Residue [MHz]')
    set(gca,'xscale','log')
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_fit_glue_resid.png' ])
    
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
    % nao faz sentido o BG dos canais PC. so devemos olhar do AN.
    %plot([0,30],[1,1].*bg(4),'-r')
    plot([0,30],[1,1].*(3*bg_std(3)),'--b')
    plot([0,30],[1,1].*(10*bg_std(3)),'--b')
    %plot([0,30],[1,1].*bg(5),'-b')
    %plot([0,30],[1,1].*(3*bg_std(5)),'--c')
    %plot([0,30],[1,1].*(10*bg_std(5)),'--c')
    legend('AN387','PC387','Glued387','Fit region AN','Fit region PC', 'PC408')

    set(gca,'yscale','log')
    xlim([0 30])
    title(['data from ' datestr(jd(1)) ' to ' datestr(jd(end))])
    xlabel('Range [m]')
    ylabel(['Signal PC or AN [MHz or mv]'])
    saveas(gca,[figure_dir datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_Signal_glue.png' ])
    
    % --------------------------------------------------------------------------
    % CALIBRATION
    % --------------------------------------------------------------------------
    
    [fitAB,fitA] = ...
        calib_function_v1(altitude_sonde, wvmix_sonde, ...
                          altitude, Signal_glued, P(:,5), ...
                          jd(1), jd(end), SNR_glue, SNR(:,5));
    % bug: we should not use the STD of the BG, but instead the SNR
    % of the signal at each point. 
    %                      jd(1), jd(end), bg_std(4), bg_std(5));

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
   
    fname=[datestr(yu,'yyyy_mm_dd_HH_MM_SS') '_Signal_calibrado' ];
    saveas(figure(1542),[figure_dir fname '.png' ])
    set(figure(1542),'paperunits','points','papersize',[560 420])
    print(figure(1542),[figure_dir fname '.pdf' ],'-dpdf','-bestfit')

    %%
    
    fname_save = [prepro_dir 'Calib_' datestr(yu,'yyyy_mm_dd_HH_MM_SS') '.mat' ];
    save(fname_save,'fitAB','fitA',...
         'jd','jdi','altitude','range','Signal_glued','P',...
         'altitude_sonde','wvmix_sonde','I_sonde', 'cfun_glue',...
         's_max250','s_max3000')
    
end

%fim
