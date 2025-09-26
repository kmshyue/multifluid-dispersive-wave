c
c----------------------------------------------------------------
c
      subroutine src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &           xlower,ylower,dx,dy,q,qold,aux,auxold,
     &           work,mwork,t,dt)
      implicit double precision(a-h,o-z)
      common /cgrav/ grav
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension work(mwork)
      dimension qloc(20)
      dimension psi(20)
c
c     # source terms 

c     # (1) cylindrical symmetry about y-axis (so x=radius)
c     # (2) gravity
c
c     # forward Euler method
c
c     # cylindrical symmetry about y-axis (so x=radius)
c
      ndim = 2

      do i=1,mx
         xc = xlower+dble(i-1)*dx+0.5d0*dx
         do j=1,my
            do m=1,meqn
               qloc(m) = q(i,j,m)
               enddo
c
            call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,rhoh0,p0,
     &           grue0,zeta0,zfa0,zfb0,c20)
c
            gaxis  = -dble(ndim-1)/xc
c
            psi(1) = gaxis*qloc(1)*vx0
            psi(2) = gaxis*qloc(2)*vx0
            psi(3) = gaxis*qloc(3)*vx0
            psi(4) = gaxis*qloc(4)*vx0
            psi(5) = gaxis*rhoh0*vx0
            psi(6) = 0.0d0     
c
            do m=1,meqn
               q(i,j,m) = q(i,j,m)+dt*psi(m)
               enddo
            enddo   
         enddo   
c
c     # gravity
      do i=1,mx
         do j=1,my
            psi(4) = -grav*(q(i,j,1)+q(i,j,2))
            psi(5) = -grav*q(i,j,4)
c
            q(i,j,4) = q(i,j,4)-dt*psi(4)
            q(i,j,5) = q(i,j,5)+dt*psi(5)
            enddo   
         enddo   
      return
c
      end
