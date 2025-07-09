function ploteta(j)  
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
data1  = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
data1 = data1';
%
km = 1e3;
h0 = 4*km;
plot(data1(:,1)/km,data1(:,3)-h0,'b-',...
     'LineWidth',1)
%
[min(min(data1(:,3)-h0)) max(max(data1(:,3)-h0))]
%
title(['surface displacement ($RC=300$m)'],...
       'fontsize',20,'interpreter','latex')
%
legend(['$t=$', num2str(t1),'s'],...
       'fontsize',20,'interpreter','latex',...
       'Location','NorthWest',...
       'Box','off')
%
%ylim([-1200 200])
%xlim([0 1/2])
%
xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
ylabel('m','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
grid on
%pname = ['crater_5Eqns_transport_RC300_eta_t0'];
%printpdf(pname)      
