function ploteta(j)  
km = 1e3;
h0 = 4e3;
%
% linearized SGN solution
%fname = ['../../airy/Airy_RC100_t',num2str(tj),'.mat'];
%load(fname,'r','eta')
%
%airy_r = r;
%airy_eta = eta;
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
wave_data  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
wave_data = wave_data';
%
plot(wave_data(:,1)/km,wave_data(:,3)-h0,'r-',...
     'LineWidth',1)
hold on
%plot(airy_r/km,airy_eta,'k-',...
%     'LineWidth',2)
%
title(['time ', num2str(t1),' sec (Gaussian $RC=100$m)'],...
       'fontsize',20,'interpreter','latex')
%
%legend('wave $10$m',...
%       'Airy',...
%       'fontsize',20,'interpreter','latex',...
%       'Location','best',...
%       'box','off')
%
%ylim([-2 2])
%
xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
ylabel('surface displacement (m)','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
grid on
%
%framest = [num2str(j)];
%pname = ['crater_wave_gaussian_RC100_airy_eta' framest];
%printpdf(pname)      
