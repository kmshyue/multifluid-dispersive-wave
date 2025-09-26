c
c ---------------------------------------------------------------
c
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &           auxl,auxr,wave,s,amdq,apdq)
      implicit double precision (a-h,o-z)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension amdq(1-mbc:maxm+mbc,meqn)
      dimension apdq(1-mbc:maxm+mbc,meqn)
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension ql_state(20),qr_state(20)
      dimension aux1(30),aux2(30)
      dimension delta(20)
      dimension wave_local(20,20)
      dimension s_local(20)
      dimension speeds(20,2)
      dimension rot(4)
c
c     # Riemann solver in the normal direction for the Euler
c     # equations on a curvilinear grid
c
      call get_aux_locations_n(ixy,mcapa,locrot,locarea)
c
      do i=2-mbc,mx+mbc-1
         rot(1) = auxl(i,locrot)
         rot(2) = auxl(i,locrot+1)
         call compute_tangent(rot)

         do m=1,meqn
            ql_state(m) = qr(i-1,m)
            qr_state(m) = ql(i,m)
            enddo
c
         call state2_rotate(rot,ql_state)
         call state2_rotate(rot,qr_state)
c
         do ma=1,maux
            aux1(ma) = auxr(i-1,ma)
            aux2(ma) = auxl(i,ma)
            enddo
c
         do m=1,meqn
            delta(m) = qr_state(m)-ql_state(m)
            enddo
c
         call rp_solver(ql_state,qr_state,aux1,aux2,
     &        delta,s_local,wave_local,meqn)
c
         area = auxl(i,locarea)
         do mw=1,mwaves
            call state2_rotate_tr(rot,wave_local(1,mw))
            speeds(mw,1) = area*dmin1(s_local(mw),0.d0)
            speeds(mw,2) = area*dmax1(s_local(mw),0.d0)
            s(i,mw) = speeds(mw,1)+speeds(mw,2)
            enddo
c
         do m=1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw=1,mwaves
               wave(i,m,mw) = wave_local(m,mw)
               amdq(i,m) = amdq(i,m)+speeds(mw,1)*wave(i,m,mw)
               apdq(i,m) = apdq(i,m)+speeds(mw,2)*wave(i,m,mw)
               enddo
            enddo
         enddo
      return
c
      end
