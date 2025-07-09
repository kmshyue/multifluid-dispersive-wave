function plotetaN(N1,N2)           
% read time and number of mgrids:
clf
%
k  = 0;
km = 1e3;
%
for j=N1:N2
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
    plot(data(:,1)/km,data(:,3)-4000,'-',...
         'LineWidth',1)
%
    if j==N1 
       hold on
    end
%
    tt(k+1) = t1;
    k = k+1;
end

title(['surface displacement after impact ($RC=300$m)'],...
       'fontsize',20,'interpreter','latex')%
%
legend(['$t=$', num2str(tt(1)),'s'],...
       ['$t=$', num2str(tt(2)),'s'],...
       ['$t=$', num2str(tt(3)),'s'],...
       ['$t=$', num2str(tt(4)),'s'],...
       ['$t=$', num2str(tt(5)),'s'],...
       ['$t=$', num2str(tt(6)),'s'],...
       ['$t=$', num2str(tt(7)),'s'],...
       ['$t=$', num2str(tt(6)),'s'],...
       ['$t=$', num2str(tt(9)),'s'],...
       ['$t=$', num2str(tt(10)),'s'],...
       ['$t=$', num2str(tt(11)),'s'],...
       'fontsize',20,'interpreter','latex',...
       'Location','NorthWest',...
       'Box','off')
%
grid on

%xlim([0 50])
%ylim([-10 10])
xlabel('radial distance(km)','fontsize',20,'interpreter','latex')
ylabel('m','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
pname = 'crater_5eqns_transport_sgEOS_RC300_eta';
printpdf(pname)
end
