c
c -------------------------------------------------------------
c
      subroutine limiter(maxm,meqn,mwaves,mbc,mx,
     &           sigma,wave,s,mthlim)
      implicit double precision(a-h,o-z)
      dimension mthlim(mwaves)
      dimension sigma(1-mbc:maxm+mbc,meqn,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension work1(20),work2(20)
c
c     # Apply a limiter to the waves
c     # component by component
c
      do mw=1,mwaves
         do i=0,mx+1
            do m=1,meqn
               sigma(i,m,mw) = 0.d0
               enddo
c
            if (s(i,mw) .gt. 0.d0) then
                do m=1,meqn
                   work2(m) = wave(i-1,m,mw)         
                   enddo
            else
                do m=1,meqn
                   work2(m) = wave(i+1,m,mw)
                   enddo
            endif
c
            do m=1,meqn
               work1(m) = wave(i,m,mw)
               enddo
c
c           # limit over each component of waves separately
            do m=1,meqn
               if (dabs(work1(m)) .gt. 1.d-20) then
                   sigma(i,m,mw) = 
     &              philim(work1(m),work2(m),mthlim(mw))*
     &                        work1(m)
               endif
               enddo
            enddo   
         enddo   
      return
c   
      end
