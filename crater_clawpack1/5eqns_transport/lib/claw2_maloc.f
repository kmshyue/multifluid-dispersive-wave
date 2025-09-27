c
c ---------------------------------------------------------------------
c
      subroutine claw_maloc(mx,my,maxmx,maxmy,meqn,mwaves,mbc,
     &           maux,mwork,method)
      implicit double precision (a-h,o-z)
      dimension method(10)
c
c     # check that enough storage has been allocated:
c
      maxm = max0(maxmx,maxmy)
c
c     # in claw2
      mwork1 = 0
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+2*(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)
      mwork1 = mwork1+(maxm+2*mbc)
      mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
      mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*maux
      mwork1 = mwork1+(maxm+2*mbc)*maux
      mwork1 = mwork1+(maxm+2*mbc)*maux
      mwork1 = mwork1+(maxm+2*mbc)*maux  
c
c     # dq
      mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
c
c     # in step2
      mwork1 = mwork1+(maxm+2*mbc)*meqn*mwaves
      mwork1 = mwork1+(maxm+2*mbc)*mwaves
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn
      mwork1 = mwork1+(maxm+2*mbc)*meqn*mwaves
c
      if (method(8) .eq. 1 .or.
     &    method(8) .eq. 12) then
c         # additional storage for moving mesh in b4step2
          mwork1= mwork1+
     &            (maxmx+2*mbc+1)*(maxmy+2*mbc+1)*2
          mwork1= mwork1+
     &            (maxmx+2*mbc+1)*(maxmy+2*mbc+1)*2
          mwork1= mwork1+
     &            (maxmx+2*mbc+1)*(maxmy+2*mbc+1)*2
          mwork1= mwork1+(maxm+2*mbc)*meqn
          mwork1= mwork1+(maxm+2*mbc)*meqn
          mwork1= mwork1+(maxm+2*mbc)*meqn
          mwork1= mwork1+(maxm+2*mbc)*meqn
      endif
c
      if (method(8) .eq. 2 .or.
     &    method(8) .eq. 12) then
c         # additional storage for interface sharpening in a4step2
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
      endif
c
      if (method(8) .eq. 3) then
c         # additional storage for anti-diffusion 
c         # interface sharpening in a4step2
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)*2
          mwork1 = mwork1+
     &             (maxmx+2*mbc)*(maxmy+2*mbc)*2
          mwork1 = mwork1+(maxm+2*mbc)
      endif
c
      if (method(9) .eq. 1) then
c         # in an2step (diffusive step for drift flux)
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
      elseif (method(9) .eq. 4) then
c         # memory in a4step1 (for Helmholtz solver step)
c
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)
      endif
c
      if (method(10) .eq. 1) then
c         # in an2step (CSF step for surface tension)
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*2
      elseif (method(10) .eq. 2) then
c         # in b4step2 (VSF step for surface tension)
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
      elseif (method(10) .eq. 3) then
c         # in a4step2 (split VSF step for surface tension)
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)*meqn
      elseif (method(10) .eq. 4) then
c         # in a4step2 (split VSF step for surface tension)
c         # anti-diffusion on volume fraction included
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
      elseif (method(10) .eq. 5) then
c         # in a4step2 (split VSF step for surface tension) 
c
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)
          mwork1 = mwork1+(maxm+2*mbc)*meqn
c
c         # in an2step (split viscous step for diffusion)
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
          mwork1 = mwork1+(maxm+2*mbc)*meqn
c
c         # in step2_vis_rkc
          mwork1 = mwork1+(maxmx+2*mbc)*(maxmy+2*mbc)*meqn
      endif
c
      if (mx .gt. maxmx .or. my .gt. maxmy .or.
     &    mwork .lt. mwork1) then
c
          write(6,*) 'mwork=',mwork,', mwork1=',mwork1
c
c         # insufficient storage
          maxmx1 = max0(mx,maxmx)
          maxmy1 = max0(my,maxmy)
          maxm1  = max0(maxmx1,maxmy1)
c
c         # in claw2
          mwork1 = 0
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+2*(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)
          mwork1 = mwork1+(maxm1+2*mbc)
          mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
          mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*maux
          mwork1 = mwork1+(maxm1+2*mbc)*maux
          mwork1 = mwork1+(maxm1+2*mbc)*maux
          mwork1 = mwork1+(maxm1+2*mbc)*maux
c
c         # dq
          mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
c
c         # in step2
          mwork1 = mwork1+(maxm1+2*mbc)*meqn*mwaves
          mwork1 = mwork1+(maxm1+2*mbc)*mwaves
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn
          mwork1 = mwork1+(maxm1+2*mbc)*meqn*mwaves
c
          if (method(8) .eq. 1 .or.
     &        method(8) .eq. 12) then
c             # additional storage for moving mesh in b4step2
              mwork1= mwork1+
     &                (maxmx1+2*mbc+1)*(maxmy1+2*mbc+1)*2
              mwork1= mwork1+
     &                (maxmx1+2*mbc+1)*(maxmy1+2*mbc+1)*2
              mwork1= mwork1+
     &                (maxmx1+2*mbc+1)*(maxmy1+2*mbc+1)*2
              mwork1= mwork1+(maxm1+2*mbc)*meqn
              mwork1= mwork1+(maxm1+2*mbc)*meqn
              mwork1= mwork1+(maxm1+2*mbc)*meqn
              mwork1= mwork1+(maxm1+2*mbc)*meqn
          endif
c
          if (method(8) .eq. 2 .or.
     &        method(8) .eq. 12) then
c             # additional storage for interface sharpening in a4step2
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+
     &                 (maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
          endif
c
          if (method(8) .eq. 3) then
c             # additional storage for anti-diffusion 
c             # interface sharpening in a4step2
              mwork1 = mwork1+
     &                 (maxmx+2*mbc)*(maxmy+2*mbc)*2
              mwork1 = mwork1+
     &                 (maxmx+2*mbc)*(maxmy+2*mbc)*2
              mwork1 = mwork1+(maxm+2*mbc)
          endif
c
          if (method(9) .eq. 1) then
c             # in an2step (diffusive step for drift flux)
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
          elseif (method(9) .eq. 4) then
c             # memory in a4step2 (for Helmholtz solver step)
c
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)
          endif
c
          if (method(10) .eq. 1) then
c             # in an2step (CSF step for surface tension)
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*2
          elseif (method(10) .eq. 2) then
c             # in b4step2 (VSF step for surface tension)
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
          elseif (method(10) .eq. 3) then
c             # in a4step2 (split VSF step for surface tension)
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
          elseif (method(10) .eq. 4) then
c             # in a4step2 (split VSF step for surface tension)
c             # anti-diffusion on volume fraction included
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
          elseif (method(10) .eq. 5) then
c             # in a4step2 (split VSF step for surface tension)
c
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
c
c             # in an2step (split viscous step for diffusion)
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
              mwork1 = mwork1+(maxm1+2*mbc)*meqn
c
c             # in step2_vis_rkc
              mwork1 = mwork1+(maxmx1+2*mbc)*(maxmy1+2*mbc)*meqn
          endif
c
          write(6,*) ' '
          write(6,*) '*** ERROR *** Insufficient storage allocated'
          write(6,*) 'Recompile after increasing values in driver.f:'
          write(6,611) maxmx1
          write(6,612) maxmy1
          write(6,613) mwork1
c
          stop
      endif
      return
c
 611  format(/,'parameter (maxmx = ',i5,')')
 612  format('parameter (maxmy = ',i5,')')
 613  format('parameter (mwork = ',i9,')',/)
      end
