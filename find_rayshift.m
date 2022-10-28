%close all

%clear all
%[nfile, head, chphy]=profile_read_dir('2015/6/10',[],[],[],4000);

%[nfile, head, chphy]=profile_read_dir('2012_11_29_telescope/NocLight',[],[],[],4000);


bg=mean(chphy(1).data(end-500:end,:))*0;

ndata=head(1).ch(1).ndata;
r=7.5*(1:ndata);
r2=r.*r;

clear jd phy
j=0;
%for i=1:nfile
for i=1308:1368
    j=j+1;
    jd(j)=head(i).jdi;
    phy(:,j)=chphy(1).data(:,i)-bg(i);
end

clear mr2
for i=1:j
  mr2(:,i)=r2;
end

figure(1); clf; 
imagesc(phy)
caxis(quantile(reshape(phy,1,numel(phy)),[0.01 0.99]))
colorbar; grid

figure(2); clf; 
imagesc(phy.*mr2)
caxis(quantile(reshape(phy.*mr2,1,numel(phy)),[0.01 0.99]))
colorbar; grid

figure(3); clf; 
plot(mod(jd,1)*24,'*')
grid

figure(4); clf; 
plot(mean((phy.*mr2)'))
grid

figure(5); clf; 
plot(mean(phy'))
grid
%