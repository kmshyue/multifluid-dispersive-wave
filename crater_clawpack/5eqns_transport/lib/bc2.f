c
c ----------------------------------------------------------
c
      subroutine bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &           dx,dy,q,maux,aux,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      parameter(rhoeps = 1.0d-8)
      common /grav_steady_info/ q0,ubc,H0,gS0,phibc,
     &        geos0,bn0
      common /cgrav/ grav
      integer mthbc(4)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qbc(20)
      dimension rot(4),uv(2)
      external grav_steady_sg,dgrav_steady_sg
c
c     # left boundary:
      go to (100,110,120,130) mthbc(1)+1
      goto 199
c
 100  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
      stop
      go to 199
c
 110  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(1,j,m)
               enddo
            enddo
         enddo
      go to 199

 120  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
               enddo
            enddo
         enddo
      go to 199
c
 130  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         do j=1-mbc,my+mbc
            q(1-ibc,j,3) = -q(1-ibc,j,3)
            enddo
         enddo
      go to 199
c
 199  continue
c     # right boundary:
      go to (200,210,220,230) mthbc(2)+1
      goto 299
c
 200  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

 210  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
               enddo
            enddo
         enddo
      go to 299

 220  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
               enddo
            enddo
         enddo
      go to 299

 230  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            do j=1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         do j=1-mbc,my+mbc
            q(mx+ibc,j,3) = -q(mx+ibc,j,3)
            enddo
         enddo
      go to 299

 299  continue
c     # bottom boundary:
      go to (300,310,320,330) mthbc(3)+1
      goto 399
c
 300  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
 310  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
               enddo
            enddo
         enddo
      go to 399

 320  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
               enddo
            enddo
         enddo
      go to 399

 330  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
c     # (for a general quadrilateral grid)
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            rot(1) = aux(i,1,4)
            rot(2) = aux(i,1,5)
            call compute_tangent(rot)
            uv(1) = q(i,jbc,3)
            uv(2) = q(i,jbc,4)
            call rotate2(rot,uv)
            uv(1) = -uv(1)
            call rotate2_tr(rot,uv)
            q(i,1-jbc,3) = uv(1)
            q(i,1-jbc,4) = uv(2)
            enddo
         enddo
      go to 399

 399  continue
c     # top boundary:
      go to (400,410,420,430) mthbc(4)+1
      goto 499
c
 400  continue
c     # gravitational source terms
c     # steady-state equilibrium conditions
c
      do jbc=1,mbc
         do i=1-mbc,mx+mbc
            xc = xlower+dble(i-1)*dx+0.5d0*dx
            yc = ylower+dble(my+1-jbc-1)*dy+0.5d0*dy
c
            call mapc2p(xc,yc,xp,yp)
            phi0 = grav*yp
c
            ybc = ylower+dble(my+1-jbc+1-1)*dy+0.5d0*dy
            call mapc2p(xc,ybc,xp,yp)
            phibc = grav*yp
c
            do m=1,meqn
               qbc(m) = q(i,my+1-jbc,m)
               enddo
c
            call ctomfd(qbc,rhoa0,rhob0,u0,v0,rhoh0,p0,
     &           grue0,pref0,zfa0,zfb0,c20)
c
            rho0 = rhoa0+rhob0
            Y10  = rhoa0/rho0
            Y20  = rhob0/rho0
c
c           # constant momentum
            q0 = rho0*v0
            w0 = rho0*u0*v0
c
c           # constant entropy
            geos0 = grue0+1.0d0
            bn0   = -pref0*grue0/(grue0+1.0d0)
            S0    = (p0+bn0)/rho0**geos0
c
c           # constant modified enthalpy
            H0 = (qbc(5)+p0)/rho0+phi0
c
c           # transverse velocity
            ubc = u0
c
            gS0 = geos0/(geos0-1.0d0)*S0
c
c           # solve for density
            rho00 = rho0
            rhobc = zero1d_aux(rho00,grav_steady_sg,
     &              dgrav_steady_sg,rhoeps,iflag)
c
            if (iflag .ne. 1) then
c                write(66,*) 'error in top bc: time=',t,
c     &                       ', iflag=',iflag
c                write(66,*) 'qbc=',(qbc(m),m=1,meqn)
c                write(66,*) 'rho0=',rho0,u0,v0,p0
c
c                stop
            else 
                vbc = q0/rhobc
                pbc = S0*rhobc**geos0-bn0
c
                rho1_bc = Y10*rhobc/zfa0
                rho2_bc = Y20*rhobc/zfb0
c
                call prmtoc(qbc,rho1_bc,rho2_bc,ubc,vbc,
     &               pbc,pbc,zfa0,zfb0)
            endif
c
            do m=1,meqn
               q(i,my+jbc,m) = qbc(m)
               enddo
            enddo
         enddo
      go to 499
c
 410  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
               enddo
            enddo
         enddo
      go to 499
c
 420  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
               enddo
            enddo
         enddo
      go to 499
c
 430  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            do i=1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
               enddo
            enddo
         enddo
c
c     # negate the normal velocity:
c     # (for a general quadrilateral grid)
c
      do jbc=1,mbc 
         do i=1-mbc,mx+mbc
            rot(1) = aux(i,my+1,4)
            rot(2) = aux(i,my+1,5)
            call compute_tangent(rot)
c
            uv(1) = q(i,my+1-jbc,3)
            uv(2) = q(i,my+1-jbc,4)
            call rotate2(rot,uv)
c
            uv(1) = -uv(1)
            call rotate2_tr(rot,uv)
c
            q(i,my+jbc,3) = uv(1)
            q(i,my+jbc,4) = uv(2)
            enddo
         enddo
      go to 499
c
 499  continue
      return
c
      end
c
c ------------------------------------------------------------
c
      double precision function grav_steady_sg(rho)
      implicit double precision (a-h,o-z)
      common /grav_steady_info/ q0,u0,H0,gS0,phi,geos0,bn0
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      rhoK = 0.5d0*(q0**2/rho+rho*u0**2)
      grav_steady_sg = rho*(H0-phi)-
     &             (gS0*rho**geos0+rhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', grav_steady_sg=',grav_steady_sg
      return
c
      end
c
c ------------------------------------------------------------
c
      double precision function dgrav_steady_sg(rho)
      implicit double precision (a-h,o-z)
      common /grav_steady_info/ q0,u0,H0,gS0,phi,geos0,bn0
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      drhoK = 0.5d0*(-q0**2/rho**2+u0**2)
      dgrav_steady_sg = (H0-phi)-
     &       (geos0*gS0*rho**(geos0-1.0d0)+drhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', dgrav_steady_sg=',dgrav_steady_sg
      return
      end
