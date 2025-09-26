c
c ---------------------------------------------------------------
c           
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &           aux1,aux2,aux3,ilr,wave,s,asdq,bmasdq,bpasdq)
      implicit double precision (a-h,o-z)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension asdq(1-mbc:maxm+mbc,meqn)
      dimension bmasdq(1-mbc:maxm+mbc,meqn)
      dimension bpasdq(1-mbc:maxm+mbc,meqn)
      dimension aux1(1-mbc:maxm+mbc,maux)
      dimension aux2(1-mbc:maxm+mbc,maux)
      dimension aux3(1-mbc:maxm+mbc,maux)
      dimension speeds(20,2),s_local(20),wave_local(20,20)
      dimension qlm(20),qrm(20),qlp(20),qrp(20)
      dimension deltam(20),deltap(20)
      dimension rotm(4),rotp(4)
c
c     # Riemann solver in the transverse direction for the Euler
c     # equations on a curvilinear grid
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Use the same idea as in rpn2 but now rotate into the direction
c     # normal to the cell edge above or below this cell.
c
      call get_aux_locations_t(ixy,mcapa,locrot,locarea)
c
      do i=2-mbc,mx+mbc
         i1 = i+ilr-2
c
         do m=1,meqn
            qlm(m) = qr(i-1,m)   
            qlp(m) = qr(i-1,m)
            qrm(m) = ql(i,m)    
            qrp(m) = ql(i,m)     
c
            deltap(m) = asdq(i,m)
            deltam(m) = asdq(i,m)
            enddo
c
c        # minus & plus transverse directions 
         do j=1,2
            rotm(j) = aux2(i1,locrot+j-1)
            rotp(j) = aux3(i1,locrot+j-1)
            enddo
         call compute_tangent(rotm)
         call compute_tangent(rotp)
c
         call state2_rotate(rotm,qlm)
         call state2_rotate(rotm,qrm)
         call state2_rotate(rotp,qlp)
         call state2_rotate(rotp,qrp)
c
         call state2_rotate(rotm,deltam)
         call state2_rotate(rotp,deltap)
c
c        # solve for minus side
         call rpt_solver(qlm,qrm,deltam,s_local,wave_local,meqn)
c
         area = aux2(i1,locarea)
         do mw=1,mwaves
            call state2_rotate_tr(rotm,wave_local(1,mw))
            speeds(mw,1) = area*dmin1(s_local(mw),0.d0)
            enddo
c
         do m=1,meqn
            bmasdq(i,m) = 0.d0
            do mw=1,mwaves
               bmasdq(i,m) = bmasdq(i,m)+
     &                       speeds(mw,1)*wave_local(m,mw)
               enddo
            enddo
c
c        # solve for plus side
         call rpt_solver(qlp,qrp,deltap,s_local,wave_local,meqn)
c
         area = aux3(i1,locarea)
         do mw = 1,mwaves
            call state2_rotate_tr(rotp,wave_local(1,mw))
            speeds(mw,2) = area*dmax1(s_local(mw),0.d0)
            enddo
c
         do m=1,meqn
            bpasdq(i,m) = 0.d0
            do mw=1,mwaves
               bpasdq(i,m) = bpasdq(i,m)+ 
     &                       speeds(mw,2)*wave_local(m,mw)
               enddo
            enddo
         enddo                         
      return
c      
      end
