function crater_5eqns_eta_movie(N)           
clf
km = 1e3;
h0 = 4e3;
%
for j=0:N
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
    data_5eqns  = fscanf(fid,'%g',[3 inf]);
    status = fclose(fid);
    data_5eqns = data_5eqns';
%
    plot(data_5eqns(:,1)/km,data_5eqns(:,3)-h0,'b-',...
         'LineWidth',1)
%
    title(['time t=', num2str(t1),' seconds after impact ($RC=300$m)'],...
       'fontsize',20,'interpreter','latex')
%
    legend('$2$-phase flow solution',...
       'fontsize',20,'interpreter','latex',...
       'Location','NorthWest',...
       'box','off')
%
%    ylim([-10 10])
%
    xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
    ylabel('water height (m)','fontsize',20,'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
    grid on
    F(j+1) = getframe(gcf);
end
%
v = VideoWriter('crater_5eqns_RC300_eta.mp4','MPEG-4');
open(v);
writeVideo(v,F);
close(v);
end
