c
c ------------------------------------------------------------
c
      subroutine claw2(maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,work,mwork,
     &           info,bc2,rpn2,rpt2,src2,b4step2,a4step2)
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2,a4step2
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension work(mwork)
      dimension mthlim(mwaves),method(10),nv(2),mthbc(4)
      dimension dtv(5),cflv(4)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c     # Solves a hyperbolic system of conservation laws in two space 
c     # dimensions of the general form
c  
c     capa * q_t + A q_x + B q_y = psi
c
      maxm   = max0(maxmx,maxmy)
      info   = 0
      t      = tstart
      maxn   = nv(1)
      dt     = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin  = dt
      dtmax  = dt
      nv(2)  = 0
c
c     # check for errors in data:
c     ---------------------------
c
      if (mx.gt.maxmx .or. my.gt.maxmy .or. mbc.lt.2) then
          info = 1
          write(6,*) 'CLAW2 ERROR...  check mx,maxmx,my,maxmy,mbc'
          go to 900
      endif
c
      if (method(1) .eq. 0) then
c         # fixed size time steps.  Compute the number of steps:
          if (tend .lt. tstart) then
c             # single step mode 
              maxn = 1
          else
              maxn = (tend-tstart+1d-10)/dt
              if (dabs(maxn*dt-(tend-tstart)) .gt.
     &                 1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 write(6,*) 
     &               'CLAW2 ERROR... dt does not divide (tend-tstart)'
                 go to 900
              endif
           endif
      endif
c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
          info = 3
          write(6,*) 'CLAW2 ERROR...  cflv(2) > cflv(1)'
          go to 900
      endif
c
      if (method(6).gt.method(7)) then
          info = 5
          write(6,*) 'CLAW2 ERROR...  method(6) > method(7)'
          go to 900
      endif
c
      if (method(2) .lt. method(3)) then
          info = 6
          write(6,*) 'CLAW2 ERROR...  method(3) > method(2)'
          go to 900
      endif
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
      i0qadd  = 1
      i0fadd  = i0qadd+(maxm+2*mbc)*meqn
      i0gadd  = i0fadd+(maxm+2*mbc)*meqn
      i0q1d   = i0gadd+2*(maxm+2*mbc)*meqn 
      i0dtdx1 = i0q1d+(maxm+2*mbc)*meqn  
      i0dtdy1 = i0dtdx1+(maxm+2*mbc)
      i0qwork = i0dtdy1+(maxm+2*mbc)
      i0awork = i0qwork+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
      i0aux1  = i0awork+(maxmx+2*mbc)*(maxmy+2*mbc)*maux
      i0aux2  = i0aux1+(maxm+2*mbc)*maux
      i0aux3  = i0aux2+(maxm+2*mbc)*maux
      i0next  = i0aux3+(maxm+2*mbc)*maux  !# next free space
      mused   = i0next-1                  !# space already used
      mwork1  = mwork-mused               !# remaining space 
c
c     # main loop
c
      if (maxn .eq. 0) go to 900
      do 100 n=1,maxn
         told = t   
c
c        # adjust dt to hit tend exactly if we're near end of computation
c        # (unless tend < tstart, which is a flag to take only a single step)
         if ((told+dt .gt. tend) .and. 
     &       (tstart  .lt. tend)) dt = tend-told
c
c        # store dt and t in the common block comxyt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
         dycom = dy
c
c        # extend data from grid to bordering boundary cells:
         call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &        dx,dy,q,maux,aux,told,dt,mthbc)
c
         call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwork))
c
         call copyaux2(maxmx,maxmy,maux,mbc,mx,my,aux,work(i0awork))
c
         call b4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &        xlower,ylower,dx,dy,told,dt,
     &        q,work(i0qwork),aux,work(i0awork),
     &        work(i0next),mwork1,mthlim,mwaves,mthbc,
     &        cflv(2))
c
         if (method(5) .eq. 2) then
c            # source terms over a half time step:
             dt2 = dt/2.d0
c
             call src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &            xlower,ylower,dx,dy,
     &            q,work(i0qwork),aux,work(i0awork),
     &            work(i0next),mwork1,told,dt2)
c
             call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &            dx,dy,q,maux,aux,told,dt,mthbc)
         endif
c
c        # take one step on the conservation law
         if (method(3) .ge. 0) then
c            # unsplit version
c             
             call step2(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,
     &            mx,my,q,work(i0qwork),aux,work(i0awork),
     &            xlower,ylower,dx,dy,told,dt,
     &            method,mthlim,mthbc,cfl,
     &            work(i0qadd),work(i0fadd),work(i0gadd),
     &            work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &            work(i0aux1),work(i0aux2),work(i0aux3),
     &            work(i0next),mwork1,rpn2,rpt2)
         else
c            # dimensional splitting
c            # Godunov splitting
c
             call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,
     &            mx,my,q,work(i0qwork),aux,work(i0awork),
     &            xlower,ylower,dx,dy,told,dt,
     &            method,mthlim,mthbc,cfl,
     &            work(i0qadd),work(i0fadd),work(i0gadd),
     &            work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &            work(i0aux1),work(i0aux2),work(i0aux3),
     &            work(i0next),mwork1,rpn2,rpt2)
         endif
c
c        # after step2
         call a4step2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &        xlower,ylower,dx,dy,told,dt,
     &        q,work(i0qwork),aux,work(i0awork),
     &        work(i0next),mwork1,mthlim,mwaves,mthbc)
c
         if (method(5) .eq. 1) then
c            # source terms over a full time step:
c
             call src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &            xlower,ylower,dx,dy,
     &            q,work(i0qwork),aux,work(i0awork),
     &            work(i0next),mwork1,told,dt)
         elseif (method(5) .eq. 2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!

             call src2(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &            xlower,ylower,dx,dy,
     &            q,work(i0qwork),aux,work(i0awork),
     &            work(i0next),mwork1,told,dt2)
         endif
c
c        # end of step2
         call an2step(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &        xlower,ylower,dx,dy,told,dt,q,work(i0qwork),
     &        aux,work(i0awork),work(i0next),mwork1,mthbc)
c
         t = told+dt
c
         if (method(4).eq.1) then
c            # verbose mode
             write(6,601) n,cfl,dt,t
  601        format('CLAW2... Step',i6,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
         endif
c
c        # choose new time step if variable time step
         if (method(1) .eq. 1) then
             if (cfl .gt. cflv(1)) then
                 write(6,*) 'cfl .gt. cflv(1): cfl=',cfl,', dt=',dt
c                 stop
c
                 call estdt(maxmx,maxmy,mbc,mx,my,meqn,maux,
     &                 dx,dy,dt0,cflv(2),q,aux)

                 dt = dt0
             else
                 if (cfl .eq. 0.d0) then
                     dt = dtv(2)
                 else
                     dt = dmin1(dtv(2),dt*cflv(2)/cfl)
                 endif
             endif
c
             dtmin = dmin1(dt,dtmin)
             dtmax = dmax1(dt,dtmax)
         endif
c
c        # see if we are done:
c
         nv(2) = nv(2)+1
         if (t .ge. tend) go to 900
  100    continue
c
  900  continue
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          write(6,*) 'CLAW2 ERROR...  too many timesteps'
          info = 11
       endif
c
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
       return 
c
       end
