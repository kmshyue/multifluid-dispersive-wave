function ploteta_airy(j,tj)  
km = 1e3;
h0 = 4e3;
%
% linearized SGN solution
fname = ['../RC300_small_airy/Airy_RC300_t',num2str(tj),'.mat']
load(fname,'r','eta')
%
airy_r = r;
airy_eta = eta;
%
clf
n1 = j+10000;
fname = ['fort.',num2str(n1)];
fname(6) = 't';
fid  = fopen(fname);
t1   = fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
fclose(fid);
%
fname(6) = 'c';
fid    = fopen(fname);
mfluid_data  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
mfluid_data = mfluid_data';
%
plot(mfluid_data(:,1)/km,mfluid_data(:,3)-h0,'r-',...
     'LineWidth',1)
hold on
plot(airy_r/km,airy_eta,'k-',...
     'LineWidth',2)
%
title(['time ', num2str(t1),' sec ($RC=300$m)'],...
       'fontsize',20,'interpreter','latex')
%
legend('mfluid',...
       'Airy',...
       'fontsize',20,'interpreter','latex',...
       'Location','best',...
       'box','off')
%
ylim([-50 50])
%
xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
ylabel('surface displacement (m)','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
grid on
%
%framest = [num2str(j)];
%pname = ['crater_mfluid_airy_eta' framest];
%printpdf(pname)      
