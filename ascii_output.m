clear all
close all

type={ 'AN' 'PC' };

basedir{1}='2012_11_29_telescope';
basedir{2}='2014_03_16_telescope';
basedir{3}='2014_03_27_telescope';
basedir{4}='2014_06_26_telescope';
basedir{5}='2014_06_27_telescope';
basedir{6}='2014_09_09_beforealign';
basedir{7}='2014_09_09_afteralign';
basedir{8}='2015_06_09_afteralign';
basedir{9}='2015_06_10_aftertrigger';

for ibase=1:9
%% DARK_CURRENT
if (1)

dark=[];
list=dirpath(basedir{ibase},'Dark_Current*');
for j=1:length(list)
%for j=length(list):length(list)
    [tmp tag]=fileparts(list{j});

    [nfile, head, chphy]=profile_read_dir(list{j});
    if (nfile==0) continue; end

    for i=1:head(1).nch
        % save dark current for rayleigh calculation
        dark(:,i)=mean(chphy(i).data,2);
        
        % output
        modo=sprintf('%04d_%s',head(1).ch(i).wlen, type{head(1).ch(i).photons+1});
        fname=[basedir{ibase} '/ascii/' tag '_' modo '_Manaus.txt'];
        fid=fopen(fname,'w+t');
        fprintf(fid,'Manaus\n'); %1
        fprintf(fid,'Manaus\n'); %2
        fprintf(fid,[modo '\n']); %3
        fprintf(fid,[datestr(head(1).jdi,'yyyy/mm/dd HH:MM') '\n']); %4
        fprintf(fid,'%4.2f min\n',(head(1).jdf-head(1).jdi)*1440); %5
        fprintf(fid,'Altitude\tDC\n'); %6
        for j=1:head(1).ch(1).ndata
            fprintf(fid,'%f\t%e\n',j*head(1).ch(1).binw/1e3,...
                    mean(chphy(i).data(j,:),2));
        end
        fclose(fid);
    end
end
% save dark current for rayleigh calculation
%for i=1:head(1).nch
%    dark(:,i)=mean(chphy(i).data,2);
%end
clear nfile head chphy tmp tag list modo fname fid i j
end

%% RAYLEIGH

list=dirpath(basedir{ibase},'Rayleigh*');
for j=1:length(list)
    [tmp tag]=fileparts(list{j});
    
    [nfile, head, chphy]=profile_read_dir(list{j});
    if (nfile==0) continue; end

    % for the rayleigh fit
    cte=constants(390.);
    tmp=dirpath([basedir{ibase} '/' tag '/'], '82332_*');
    snd=read_sonde_Wyoming(tmp{1});
    lambda=[355 355 387 387]*1e-9;
    mol=molecular(lambda, snd, cte);
    lidar_altitude=head(1).alt;
    r_bin=head(1).ch(1).binw;
    clear P Pr2;
    for i=1:head(1).nch
        P(:,i)=mean(chphy(i).data,2);
        % correct for non-linear dark current
        if ~isempty(dark)
            P(:,i)=P(:,i)-dark(:,i);
        end
        % correct for background
        P(:,i)=P(:,i)-mean(P(end-4000:end,i),1);
    end
    % correct for analog displacement
    P=P(10:end,:);
    % cut just the first 30km
    %P=P(1:4000,:);
    alt=(1:size(P,1))'.*r_bin;
    altsq=alt.*alt;
    %bottomlayer=7e3;
    %toplayer=10e3;
    debug=-1;
    rayleigh_fit
    %return
    
    %figure(3); clf; hold on
    %plot(mean(chphy(1).data,2))  
    %plot(dark(:,1),'r') 
    %plot(mean(chphy(1).data(:,:),2)-dark(:,1)+1.5,'k')
    %return

    % print output
    for i=1:4
        modo=sprintf('%04d_%s',head(1).ch(i).wlen, type{head(1).ch(i).photons+1});
        fname=[basedir{ibase} '/ascii/' tag '_' modo '_Manaus.txt'];
        fid=fopen(fname,'w+t');
        fprintf(fid,'Manaus\n'); %1
        fprintf(fid,'Manaus\n'); %2
        fprintf(fid,[modo '\n']); %3
        fprintf(fid,[datestr(head(1).jdi,'yyyy/mm/dd HH:MM') '\n']); %4
        fprintf(fid,'%4.2f min\n',(head(1).jdf-head(1).jdi)*1440); %5
        fprintf(fid,'radiosounding'); %6
        fprintf(fid,'7-10 km'); %7
        fprintf(fid,'Altitude\tRCS\tattBackMol\n'); %8
        for j=1:length(alt)
            fprintf(fid,'%f\t%e\t%e\n',alt(j)/1e3,...
                    P(j,i).*altsq(j), Pr2_mol(j,i));
        end
        fclose(fid);
    end
end
clear nfile head chphy tmp tag list modo fname fid i j

%return

%% ZERO_BIN

list=dirpath(basedir{ibase},'Zero_Bin*');
for j=1:length(list)
    [tmp tag]=fileparts(list{j});

    [nfile, head, chphy]=profile_read_dir(list{j});
    if (nfile==0) continue; end

    for i=1:head(1).nch
        modo=sprintf('%04d_%s',head(1).ch(i).wlen, type{head(1).ch(i).photons+1});
        fname=[basedir{ibase} '/ascii/' tag '_' modo '_Manaus.txt'];
        fid=fopen(fname,'w+t');
        fprintf(fid,'Manaus\n'); %1
        fprintf(fid,'Manaus\n'); %2
        fprintf(fid,[modo '\n']); %3
        fprintf(fid,[datestr(head(1).jdi,'yyyy/mm/dd HH:MM') '\n']); %4
        fprintf(fid,'%4.2f min\n',(head(1).jdf-head(1).jdi)*1440); %5
        fprintf(fid,'Altitude\tSignal\n'); %6
        for j=1:head(1).ch(1).ndata
            fprintf(fid,'%f\t%e\n',j*head(1).ch(1).binw/1e3,mean(chphy(i).data(j,:),2));
         end
        fclose(fid);
    end
end
clear nfile head chphy tmp tag list modo fname fid i j

%% BIN_SHIFT

list=dirpath(basedir{ibase},'Bin_Shift*');
for j=1:length(list)
    [tmp tag]=fileparts(list{j});

    [nfile, head, chphy]=profile_read_dir(list{j});
    if (nfile==0) continue; end

    for i=2:2:head(1).nch
        modo=sprintf('%04d',head(1).ch(i).wlen);
        fname=[basedir{ibase} '/ascii/' tag '_' modo '_Manaus.txt'];
        fid=fopen(fname,'w+t');
        fprintf(fid,'Manaus\n'); %1
        fprintf(fid,'Manaus\n'); %2
        fprintf(fid,[modo '_AN_PC\n']); %3
        fprintf(fid,[datestr(head(1).jdi,'yyyy/mm/dd HH:MM') '\n']); %4
        fprintf(fid,'%4.2f min\n',(head(1).jdf-head(1).jdi)*1440); %5
        fprintf(fid,'Altitude\tSignal_AN\tSignal_PC\n'); %6
        for j=1:head(1).ch(1).ndata
            fprintf(fid,'%f\t%e\t%e\n',j*head(1).ch(1).binw/1e3,...
                    mean(chphy(i-1).data(j,:),2),mean(chphy(i).data(j,:),2));
        end
        fclose(fid);
    end
end
clear nfile head chphy tmp tag list modo fname fid i j

%% TELECOVER

list1=dirpath(basedir{ibase},'north*');
list2=dirpath(basedir{ibase},'east*');
list3=dirpath(basedir{ibase},'south*');
list4=dirpath(basedir{ibase},'west*');

for j=1:length(list4)
    [nfile1, head1, chphy1]=profile_read_dir(list1{j});
    [nfile2, head2, chphy2]=profile_read_dir(list2{j});
    [nfile3, head3, chphy3]=profile_read_dir(list3{j});
    [nfile4, head4, chphy4]=profile_read_dir(list4{j});
    if (j+1<=length(list1))
        [nfile5, head5, chphy5]=profile_read_dir(list1{j+1});
    else
        nfile5=0;
        clear chphy5
        for i=1:head1(1).nch
            for k=1:head1(1).ch(1).ndata
                chphy5(i).data(k,1)=nan;
            end
        end
    end

    if (nfile1+nfile2+nfile3+nfile4+nfile5==0) continue; end
    
    for i=1:head1(1).nch
        modo=sprintf('%04d_%s',head1(1).ch(i).wlen, type{head1(1).ch(i).photons+1});
        fname=[basedir{ibase} '/ascii/telecover' num2str(j) '_' modo '_Manaus.txt'];
        fid=fopen(fname,'w+t');
        fprintf(fid,'Manaus\n'); %1
        fprintf(fid,'Manaus\n'); %2
        fprintf(fid,[modo '\n']); %3
        fprintf(fid,[datestr(head1(1).jdi,'yyyy/mm/dd HH:MM') '\n']); %4
        fprintf(fid,'Altitude\tN\tE\tS\tW\tN2\n'); %5
        for j=1:head1(1).ch(1).ndata
            fprintf(fid,'%f\t%e\t%e\t%e\t%e\t%e\n',j*head1(1).ch(1).binw/1e3,...
                    mean(chphy1(i).data(j,:),2),mean(chphy2(i).data(j,:),2),...
                    mean(chphy3(i).data(j,:),2),mean(chphy4(i).data(j,:),2),...
                    mean(chphy5(i).data(j,:),2));
        end
        fclose(fid);
    end
end


end