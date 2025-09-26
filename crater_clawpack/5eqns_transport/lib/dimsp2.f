c
c ---------------------------------------------------------------
c
      subroutine dimsp2(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &           q,qold,aux,auxold,xlower,ylower,dx,dy,told,dt,
     &           method,mthlim,mthbc,
     &           cfl,cflv,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &           aux1,aux2,aux3,work,mwork,rpn2,rpt2)
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension cflv(4)
      dimension qadd(1-mbc:maxm+mbc,meqn)
      dimension fadd(1-mbc:maxm+mbc,meqn)
      dimension gadd(1-mbc:maxm+mbc,meqn,2)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension aux1(1-mbc:maxm+mbc,maux)
      dimension aux2(1-mbc:maxm+mbc,maux)
      dimension aux3(1-mbc:maxm+mbc,maux)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension method(10),mthlim(mwaves),mthbc(4)
      dimension work(mwork)
c
c     # Take one time step, updating q, using dimensional
c     # splitting. Two choices are available:
c     #
c     # method(3) = -1   gives Godunov splitting:
c     #    time step dt in x-direction
c     #    time step dt in y-direction
c
c     # method(3) = -2   gives Strang splitting
c     #    time step dt/2 in x-direction
c     #    time step dt   in y-direction
c     #    time step dt/2 in x-direction
c
c     # Godunov splitting is recommended over Strang splitting normally
c     # since it typically works as well, is faster, and boundary
c     # conditions are handled properly.
c
c     # If method(3) = -1, take a full time step in x.
c     # If method(3) = -2, take a half time step in x.
c
      dt2 = dt/2.d0
c
      if (method(3) .eq. -2)then
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &         q,qold,aux,auxold,xlower,ylower,dx,dy,told,dt2,
     &         method,mthlim,mthbc,cflx,
     &         qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &         aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      else
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &         q,qold,aux,auxold,xlower,ylower,dx,dy,told,dt,
     &         method,mthlim,mthbc,cflx,
     &         qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &         aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      endif
c
      if (cflx .gt. cflv(1)) then
c        # Abort if the Courant number was too large in x-sweep
         cfl = cflx
         return
      endif
c
c     # take full step in y-direction
c
      call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux,told,dt,mthbc)
c
      call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &     q,q,aux,auxold,xlower,ylower,dx,dy,told,dt,
     &     method,mthlim,mthbc,cfly,
     &     qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &     aux1,aux2,aux3,work,mwork,rpn2,rpt2,2)
c
      cfl = dmax1(cflx,cfly)
c
c     # Finally, take a half time step in the x-direction
c     # if Strang splitting is used.  NOTE: boundary conditions may
c     # not be set properly for this sweep.
c
      if (method(3) .eq. -2)then
          if (cfly .gt. cflv(1)) then
c             # Abort if the Courant number was too large in y-sweep
              cfl = cfly
              return
           endif
c
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &         q,q,aux,auxold,xlower,ylower,dx,dy,told,dt2,
     &         method,mthlim,mthbc,cflx,
     &         qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &         aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
c
          cfl = dmax1(cfl,cflx)
      endif
      return
c
      end
