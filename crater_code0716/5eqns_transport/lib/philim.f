c
c --------------------------------------------------------------
c
      double precision function philim(a,b,meth)
      implicit double precision(a-h,o-z)
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
c     # NOTE: This routine is obsolete.  Instead of using limiter.f,
c     # which calls philim.f for every wave, it is more efficient to 
c     # use inlinelimiter.f, which eliminates all these function calls
c     # to philim.  If you wish to change the limiter function and are
c     # using inlinelimiter.f, the formulas must be changed in that routine.
c
      r = b/a
c
      if (meth .eq. 0) then
c         # no (LW) limiter
          philim = 1.d0
      elseif (meth .eq. 1) then
c         # minmod
          philim = dmax1(0.d0,dmin1(1.d0,r))
      elseif (meth .eq. 2) then
c         # superbee
          philim = dmax1(0.d0,dmin1(1.d0,2.d0*r), 
     &                        dmin1(2.d0,r))
      elseif (meth .eq. 3) then
c         # van Leer
          philim = (r+dabs(r))/(1.d0+dabs(r))
      elseif (meth .eq. 4) then
c         # monotinized centered 
          c = (1.d0+r)/2.d0
          philim = dmax1(0.d0,dmin1(c,2.d0,2.d0*r))
      elseif (meth .eq. 5) then
c         # Beam-Warming
          philim = r
      elseif (meth .eq. 6) then
          write(6,*) 'error in philim: meth=',meth
          stop
      endif
      return
c
      end
