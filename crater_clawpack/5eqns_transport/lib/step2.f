c
c --------------------------------------------------------------
c
      subroutine step2(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &           q,qold,aux,auxold,xlower,ylower,dx,dy,told,dt,
     &           method,mthlim,mthbc,
     &           cfl,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &           aux1,aux2,aux3,work,mwork,rpn2,rpt2)
      implicit double precision (a-h,o-z)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      external rpn2,rpt2
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension qold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension auxold(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension q1d(1-mbc:maxm+mbc,meqn)
      dimension qadd(1-mbc:maxm+mbc,meqn)
      dimension fadd(1-mbc:maxm+mbc,meqn)
      dimension gadd(1-mbc:maxm+mbc,meqn,2)
      dimension aux1(1-mbc:maxm+mbc,maux)
      dimension aux2(1-mbc:maxm+mbc,maux)
      dimension aux3(1-mbc:maxm+mbc,maux)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension method(10),mthlim(mwaves),mthbc(4)
      dimension work(mwork)
c
c     # Take one time step, updating q.
c     # On entry, q and qold should be identical and give the
c     #    initial data for this step
c     # On exit, q returns values at the end of the time step.
c     #    qold is unchanged.
c
c     # qadd is used to return increments to q from flux2
c     # fadd and gadd are used to return flux increments from flux2.
c     # See the flux2 documentation for more information.
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave  = 1
      i0s     = i0wave+(maxm+2*mbc)*meqn*mwaves
      i0amdq  = i0s+(maxm+2*mbc)*mwaves
      i0apdq  = i0amdq+(maxm+2*mbc)*meqn
      i0cqxx  = i0apdq+(maxm+2*mbc)*meqn
      i0bmadq = i0cqxx+(maxm+2*mbc)*meqn
      i0bpadq = i0bmadq+(maxm+2*mbc)*meqn
      i0sigma = i0bpadq+(maxm+2*mbc)*meqn
      iused   = i0sigma+(maxm+2*mbc)*meqn*mwaves-1
c
      if (iused .gt. mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop 
      endif
c
      mcapa = method(6)
c
      cfl  = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      if (mcapa.eq.0) then
c         # no capa array
          do i=1-mbc,maxm+mbc
             dtdx1d(i) = dtdx
             dtdy1d(i) = dtdy
             enddo
      endif
c
c     # perform x-sweeps
      do j=0,my+1
c        # copy data along a slice into 1d arrays
         do m=1,meqn
            do i=1-mbc,mx+mbc
               q1d(i,m) = qold(i,j,m)
               enddo
            enddo    
c
         if (mcapa .gt. 0)  then
             do i=1-mbc,mx+mbc
                dtdx1d(i) = dtdx/aux(i,j,mcapa)
                enddo    
         endif
c
         if (maux .gt. 0)  then
             do ma=1,maux
                do i=1-mbc,mx+mbc
                   aux1(i,ma) = aux(i,j-1,ma)
                   aux2(i,ma) = aux(i,j  ,ma)
                   aux3(i,ma) = aux(i,j+1,ma)
                   enddo
                enddo 
         endif
c
c        # store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,mwaves,maux,mbc,mx,
     &        q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &        qadd,fadd,gadd,cfl1d,
     &        work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &        work(i0sigma),work(i0cqxx),work(i0bmadq),work(i0bpadq),
     &        rpn2,rpt2)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update q by flux differencing.
c        # (rather than maintaining arrays f and g for the total fluxes,
c        # the modifications are used immediately to update q
c        # in order to save storage.)
c
         if (mcapa.eq.0) then
c            # no capa array
c            # standard flux differencing:
             do m=1,meqn
                do i=1,mx
                   q(i,j,m) = q(i,j,m)+qadd(i,m)-
     &                        dtdx*(fadd(i+1,m)-fadd(i,m))-
     &                        dtdy*(gadd(i,m,2)-gadd(i,m,1))
                   q(i,j-1,m) = q(i,j-1,m)-dtdy*gadd(i,m,1)
                   q(i,j+1,m) = q(i,j+1,m)+dtdy*gadd(i,m,2)
                   enddo
                enddo
         else
c            # with capa array  
             do m=1,meqn
                do i=1,mx
                   q(i,j,m) = q(i,j,m)+qadd(i,m)-
     &                        (dtdx*(fadd(i+1,m)-fadd(i,m))+
     &                         dtdy*(gadd(i,m,2)-gadd(i,m,1)))/
     &                         aux(i,j,mcapa)
                   q(i,j-1,m) = q(i,j-1,m)-dtdy*gadd(i,m,1)/
     &                          aux(i,j-1,mcapa)
                   q(i,j+1,m) = q(i,j+1,m)+dtdy*gadd(i,m,2)/
     &                             aux(i,j+1,mcapa)
                   enddo
                enddo
         endif
         enddo
c
c     # perform y sweeps
      do i=0,mx+1
c        # copy data along a slice into 1d arrays
         do m=1,meqn
            do j=1-mbc,my+mbc
               q1d(j,m) = qold(i,j,m)
               enddo
            enddo
c
         if (mcapa.gt.0)  then
             do j=1-mbc,my+mbc
                dtdy1d(j) = dtdy/aux(i,j,mcapa)
                enddo
         endif
c
         if (maux .gt. 0)  then
             do ma=1,maux
               do j=1-mbc,my+mbc
                  aux1(j,ma) = aux(i-1,j,ma)
                  aux2(j,ma) = aux(i,  j,ma)
                  aux3(j,ma) = aux(i+1,j,ma)
                  enddo
               enddo
         endif
c
c        # store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,mwaves,maux,mbc,my,
     &        q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,
     &        qadd,fadd,gadd,cfl1d,
     &        work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &        work(i0sigma),work(i0cqxx),work(i0bmadq),work(i0bpadq),
     &        rpn2,rpt2)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update q by flux differencing.
c        # Note that the roles of fadd and gadd are reversed for
c        # the y-sweeps -- fadd is the modification to g-fluxes and
c        # gadd is the modification to f-fluxes to the left and right.
c
         if (mcapa.eq.0) then
c            # no capa array
c            # standard flux differencing
             do m=1,meqn
                do j=1,my
                   q(i,j,m) = q(i,j,m)+(qadd(j,m)-
     &                        dtdy*(fadd(j+1,m)-fadd(j,m))-
     &                        dtdx*(gadd(j,m,2)-gadd(j,m,1)))
                   q(i-1,j,m) = q(i-1,j,m)-dtdx*gadd(j,m,1)
                   q(i+1,j,m) = q(i+1,j,m)+dtdx*gadd(j,m,2)
                   enddo
                enddo
         else
c            # with capa array  
             do m=1,meqn
                do j=1,my
                   q(i,j,m) = q(i,j,m)+qadd(j,m)-
     &                        (dtdy*(fadd(j+1,m)-fadd(j,m))+
     &                        dtdx*(gadd(j,m,2)-gadd(j,m,1)))/
     &                        aux(i,j,mcapa)
                   q(i-1,j,m) = q(i-1,j,m)-dtdx*gadd(j,m,1)/
     &                          aux(i-1,j,mcapa)
                   q(i+1,j,m) = q(i+1,j,m)+dtdx*gadd(j,m,2)/
     &                             aux(i+1,j,mcapa)
                   enddo
                enddo
         endif
         enddo
      return
c
      end
