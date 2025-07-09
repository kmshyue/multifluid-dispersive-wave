function plotur(j)           
% plot the average radial velocity
clf
%
km = 1e3;
%
nj = j+10000;
fname = ['fort.',num2str(nj)];
fname(6) = 't';
%
fid = fopen(fname);
t1 = fscanf(fid,'%g',1);        fscanf(fid,'%s',1);
meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
fclose(fid);
%
% data set
fname(6) = 'c';
fid    = fopen(fname);
data   = fscanf(fid,'%g',[3 inf]);
status = fclose(fid);
data   = data';
fid    = fopen(fname);
%
plot(data(:,1)/km,data(:,2),'b-',...
     -data(:,1)/km,data(:,2),'b-',...
     'LineWidth',1);

title(['depth-averaged radial velocity (RC=$300$m)'],...
       'fontsize',20,'interpreter','latex')%
%
grid on

%xlim([0 50])
ylim([-1/10 1/10])
xlabel('radial distance(km)','fontsize',20,'interpreter','latex')
ylabel('m/s','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
%pname = 'crater_5eqns_transport_sgEOS_RC300_ur';
%printpdf(pname)
end
