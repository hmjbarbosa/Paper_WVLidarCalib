clear all
close al

% abrir varios arquivos
datadir = '/Users/hbarbosa/DATA/lidar/data'

for yu = datenum(2011,08,27,20,00,00):datenum(2011,08,27,20,00,00)
    
    clearvars -except yu sonde jd_sonde jd_sonde_utc datadir figure_dir prepro_dir
    
    %%
    
    
    [sonde_dist_dias, I_sonde]= min(abs(jd_sonde - yu))
    
    
    if sonde_dist_dias > 0
        
        continue
        
    end
    %%
    
    
    % jdi=datenum(2014,09,14,20,00,00);       % data inicial
    % jdf= jdi + 1/24;                        % data final
    
   
    jdi = jd_sonde(I_sonde) - 0.5/24;
    jdf = jdi + 1/24;
    
    
    % dbin_MAO = [10];
    dbin_MAO = [9 -0 9 -1 -1];
    dtime_MAO =  0.004;
    [nfile, head, phy] = profile_read_dates(datadir, jdi, jdf, dbin_MAO, dtime_MAO)
    
    
end

