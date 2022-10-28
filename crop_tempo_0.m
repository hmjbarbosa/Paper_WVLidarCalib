function [Pr2_crop,jdz_crop,jdz_cropf] = crop_tempo_0(Pr2,jdz,Nmin,opt);


jdi = jdz(1) - (jdz(2)-jdz(1))/2; % assume que jdz é midpoint e dt constante
jdf = jdz(end) + (jdz(2)-jdz(1))/2;


% % Nmin = 30; % x minutos de mï¿½dia
% Nmin = 5; % x minutos de mï¿½dia
rtemp = datenum(0,0,0,0,Nmin,0);

ri = jdi:rtemp:(jdf);

TMP = NaN(size(Pr2,1),length(ri));
jd_TMP = NaN(size(ri));
jd_TMPf = NaN(size(ri));

tic
% nch = heads(1).nch;

for ch = 1:1
    for i = 1:length(ri);
        
        ind=find((jdz >= ri(i)) & (jdz < ri(i) + rtemp)); % assume que jdz é midpoint e dt constante
        
        if ~isempty(ind)
%             jd_TMP(i) = jd(ind(1));
            jd_TMP(i) = ri(i) + rtemp/2;
            jd_TMPf(i) = ri(i) + rtemp;
%             nshoots_TMP(i) = sum(nshoots(ind));
            TMP(:,i) = nanmean(Pr2(:,ind),2);
            
        end
        
%             pause
    end
    
    if opt == 0;
        cond = (jd_TMP > 0);
        Pr2_crop = TMP(:,cond);
        jdz_crop = jd_TMP(cond);
        jdz_cropf = jd_TMPf(cond);
        %     nshoots_z = nshoots_TMP(cond);
    elseif opt == 1
        
%         Pr2_crop = TMP;
%         jdz_crop = jd_TMP;
%         jdz_cropf = jd_TMPf;
        
        Pr2_crop = TMP;
        jdz_crop = ri+rtemp./2;
        jdz_cropf = ri+rtemp;
        
    end
end
toc
