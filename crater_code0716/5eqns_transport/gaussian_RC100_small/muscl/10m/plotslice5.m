function plotslice5(m,N,k)
%
clf
%
km = 1e3;
h0 = 4*km;
%
kk=0;
for j=0:N
    % read time and number of grids:
    n1 = j+10000;
    fname = ['fort.',num2str(n1)];
    fname(6) = 't';
%
    fid = fopen(fname);
    t1 = fscanf(fid,'%g',1);        fscanf(fid,'%s',1);
    meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
    ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
    fclose(fid);
%
    tt(kk+1) = t1;
%
    qmin = 1e6;
    qmax = -1e6;
%
    gname = 'fort.g0000';
    fid = fopen(gname);
    mx = fscanf(fid,'%d',1);  
    my = fscanf(fid,'%d',1);     
    mgrid = fscanf(fid,'%g %g',[2 inf]);
    status = fclose(fid);
    mgrid = mgrid';
%
    x = reshape(mgrid(:,1),mx,my);
    y = reshape(mgrid(:,2),mx,my);
    x = x/km;
    y = y/km;
%
    % data set 1
    fname(6) = 'q';
    fid      = fopen(fname);
%
    gridno = fscanf(fid,'%d',1);     fscanf(fid,'%s',1);
    level = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
    mx = fscanf(fid,'%d',1);         fscanf(fid,'%s',1);
    my = fscanf(fid,'%d',1);         fscanf(fid,'%s',1);   

    xlow = fscanf(fid,'%g',1);       fscanf(fid,'%s',1);
    ylow = fscanf(fid,'%g',1);       fscanf(fid,'%s',1);
    dx = fscanf(fid,'%g',1);         fscanf(fid,'%s',1);
    dy = fscanf(fid,'%g',1);         fscanf(fid,'%s',1);
%
    data = fscanf(fid,'%g',[meqn,mx*my]);
    data = data';
%
    if m==1
       q1 = reshape(data(:,1),mx,my)+...
            reshape(data(:,2),mx,my);
    elseif m==2
       q1 = reshape(data(:,3),mx,my);  
    elseif m==3
       q1 = reshape(data(:,4),mx,my);  
    elseif m==4
       q1 = reshape(data(:,5),mx,my);  
    elseif m==5
       q1 = 1-reshape(data(:,6),mx,my);  
    end
%
    status = fclose(fid);
%
    xslice = x(k,1);
    yslice = y(k,:);
    qslice = q1(k,:);

    if m==5
       eta(kk+1) = sum(qslice)*dx-h0
    end
%
    kk = kk+1;
%
    if j==0 
       hold on
    end

    plot(yslice,qslice,'-',...
         'LineWidth',2)

%
    title(['cross-sectional plot at r=', num2str(xslice),'km',...
           '(MUSCL, 10m)'],...
           'fontsize',20,'interpreter','latex')
    if j==0 
       hold on
    end
end
%
if m ==1
    ylabel('$\rho$','fontsize',20,'interpreter','latex')
%
%   fname = ['crater_mfluid_transport_rho' framest];
elseif m == 2
    ylabel('$u$','fontsize',20,'interpreter','latex')
%
%   fname = ['crater_mfluid_transport_u' framest];
elseif m == 3
    ylabel('$v$','fontsize',20,'interpreter','latex')
%
   fname = ['crater_muscl_gaussian_RC100_v_slice'];
elseif m == 4
    ylabel('$p$','fontsize',20,'interpreter','latex')
%
%   fname = ['crater_mfluid_transport_p' framest];
elseif m == 5
    ylabel('$\alpha_{w}$','fontsize',20,'interpreter','latex')

%   fname = ['crater_muscl_gaussian_RC100_vofw_slice'];
end
%
if m==5
   legend(['t =',num2str(tt(1)), ', $\eta=$',num2str(eta(1)),'m'],...
       ['t =',num2str(tt(2)),'s', ', $\eta=$',num2str(eta(2)),'m'],...
       ['t =',num2str(tt(3)),'s', ', $\eta=$',num2str(eta(3)),'m'],...
       ['t =',num2str(tt(4)),'s', ', $\eta=$',num2str(eta(4)),'m'],...
       ['t =',num2str(tt(5)),'s', ', $\eta=$',num2str(eta(5)),'m'],...
       ['t =',num2str(tt(6)),'s', ', $\eta=$',num2str(eta(6)),'m'],...
      'fontsize',20,'interpreter','latex',...
      'Location','NorthEast',...
      'Box','off')
elseif m==3
   legend(['t =',num2str(tt(1))],...
       ['t =',num2str(tt(2)),'s'],...
       ['t =',num2str(tt(3)),'s'],...
       ['t =',num2str(tt(4)),'s'],...
       ['t =',num2str(tt(5)),'s'],...
       ['t =',num2str(tt(6)),'s'],...
      'fontsize',20,'interpreter','latex',...
      'Location','SouthWest',...
      'Box','off')
end
grid on
%
xlabel('vetrical distance (km)','fontsize',20,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%axis off
%
%printpdf(fname)
