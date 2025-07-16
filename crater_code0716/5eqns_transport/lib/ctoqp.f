c
c ----------------------------------------------------------
c
      subroutine c2prm_1d(maxm,meqn,mbc,maux,mx,q,qp,aux)
      implicit double precision (a-h,o-z)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension qp(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
c
c     # change conservative to primitive
c     # 5-equation transport model 
c     # hybrid ideal gas and stiffened gas EOS
c     # 2-phase case
c
      do i=1-mbc,mx+mbc
         qp(i,1) = q(i,1)/q(i,6)
         qp(i,2) = q(i,2)/(1.0d0-q(i,6))
         qp(i,3) = q(i,3)/(q(i,1)+q(i,2))
         qp(i,4) = q(i,4)/(q(i,1)+q(i,2))
c
         engK  = 0.5d0*(q(i,3)**2+q(i,4)**2)/
     &                 (q(i,1)+q(i,2))
         rhoe  = q(i,5)-engK
         grue0 = 1.d0/(q(i,6)/(geos(1)-1.d0)+
     &           (1.0d0-q(i,6))/(geos(2)-1.d0))
         pref0 = (q(i,6)*geos(1)*bn(1)/(geos(1)-1.d0))+
     &      ((1.0d0-q(i,6))*geos(2)*bn(2)/(geos(2)-1.d0))
c
         qp(i,5) = grue0*(rhoe-pref0)
         qp(i,6) = q(i,6)
         enddo
      return
c
      end
c
c --------------------------------------------------------------
c
      subroutine prm2c_edge(maxm,meqn,mbc,maux,mx,
     &           ql,qr,uu,q,aux)
      implicit double precision (a-h,o-z)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      dimension q(1-mbc:maxm+mbc,meqn)
      dimension aux(1-mbc:maxm+mbc,maux)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension uu(1-mbc:maxm+mbc,meqn,2)
c
c     # change primitive variables to conservative variables
c     # 5-equation transport model 
c     # hybrid ideal gas and stiffened gas EOS
c     # 2-phase case
c
      do i=0,mx+1
c        # left cell-edge 
         ql(i,1) = uu(i,1,1)*uu(i,6,1)
         ql(i,2) = uu(i,2,1)*(1.0d0-uu(i,6,1))
         ql(i,3) = (ql(i,1)+ql(i,2))*uu(i,3,1)
         ql(i,4) = (ql(i,1)+ql(i,2))*uu(i,4,1)
c
         engK  = 0.5d0*(ql(i,3)**2+ql(i,4)**2)/
     &           (ql(i,1)+ql(i,2))
         rhoe1 = (uu(i,5,1)+geos(1)*bn(1))/(geos(1)-1.d0)
         rhoe2 = (uu(i,5,1)+geos(2)*bn(2))/(geos(2)-1.d0)
         rhoe  = uu(i,6,1)*rhoe1+(1.0d0-uu(i,6,1))*rhoe2
c
         ql(i,5) = rhoe+engK
         ql(i,6) = uu(i,6,1)
c
c        # right cell-edge
         qr(i,1) = uu(i,1,2)*uu(i,6,2)
         qr(i,2) = uu(i,2,2)*(1.0d0-uu(i,6,2))
         qr(i,3) = (qr(i,1)+qr(i,2))*uu(i,3,2)
         qr(i,4) = (qr(i,1)+qr(i,2))*uu(i,4,2)
c
         engK  = 0.5d0*(qr(i,3)**2+qr(i,4)**2)/
     &           (qr(i,1)+qr(i,2))
         rhoe1 = (uu(i,5,2)+geos(1)*bn(1))/(geos(1)-1.d0)
         rhoe2 = (uu(i,5,2)+geos(2)*bn(2))/(geos(2)-1.d0)
         rhoe  = uu(i,6,2)*rhoe1+(1.0d0-uu(i,6,2))*rhoe2
c
         qr(i,5) = rhoe+engK
         qr(i,6) = uu(i,6,2)
         enddo
c
c     # positivity-preserving check
      do i=0,mx+1
         if ((ql(i,1) .lt. 0.d0) .or.
     &       (ql(i,2) .lt. 0.d0) .or.
     &       (ql(i,6) .lt. 0.d0) .or.
     &       (1.0d0-ql(i,6) .lt. 0.d0) .or.
     &       (qr(i,1) .lt. 0.d0) .or.
     &       (qr(i,2) .lt. 0.d0) .or.
     &       (qr(i,6) .lt. 0.d0) .or.
     &       (1.0d0-qr(i,6) .lt. 0.d0)) then
c             write(66,*) 'i=',i
c             write(66,*) 'left-state'
c             write(66,*) (qr(i-1,m),m=1,meqn)
c             write(66,*) 'right-state'
c             write(66,*) (ql(i,m),m=1,meqn)
c             stop
c
              do m=1,meqn
                 ql(i,m) = q(i,m)
                 qr(i,m) = q(i,m)
                 enddo
         endif
         enddo
      return
c
      end
