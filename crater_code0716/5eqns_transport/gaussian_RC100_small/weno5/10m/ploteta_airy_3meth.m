function ploteta_airy(j,tj)  
km = 1e3;
h0 = 4e3;
%
% linearized SGN solution
fname = ['../../airy/Airy_RC100_t',num2str(tj),'.mat'];
load(fname,'r','eta')
%
airy_r = r;
airy_eta = eta;
%
clf
% weno 5
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
weno5_data  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
weno5_data = weno5_data';
%
% weno 3
fname = ['../../weno3/mesh2/fort.',num2str(n1)];
fname(24) = 'c';
fid    = fopen(fname);
weno3_data  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
weno3_data = weno3_data';
%
% muscl 
fname = ['../../muscl/mesh2/fort.',num2str(n1)];
fname(24) = 'c';
fid    = fopen(fname);
muscl_data  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
muscl_data = muscl_data';
%
plot(muscl_data(:,1)/km,muscl_data(:,3)-h0,'-',...
     weno3_data(:,1)/km,weno3_data(:,3)-h0,'-',...
     weno5_data(:,1)/km,weno5_data(:,3)-h0,'-',...
     'LineWidth',1)
hold on
plot(airy_r/km,airy_eta,'k-',...
     'LineWidth',2)
%
title(['time ', num2str(t1),' sec (Gaussian $RC=100$m)'],...
       'fontsize',20,'interpreter','latex')
%
legend('muscl  $10$m',...
       'weno 3 $10$m',...
       'weno 5 $10$m',...
       'Airy',...
       'fontsize',20,'interpreter','latex',...
       'Location','NorthEast',...
       'box','off')
%
ylim([-2 2])
%
xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
ylabel('surface displacement (m)','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
grid on
%
%framest = [num2str(j)];
%pname = ['crater_mfluid_gaussian_RC100_airy_eta' framest];
%printpdf(pname)      
