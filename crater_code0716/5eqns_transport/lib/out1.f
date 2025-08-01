c
c --------------------------------------------------------------
c
      subroutine out1(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &           dx,dy,q,t,iframe,aux,maux)
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension qloc(20),qp(20)
      dimension qave(20)
      character*10 fname
c
c     # a one-dimensional array of data
c
      fname = 'fort.cxxxx'
      nstp = iframe
      do ipos=10,7,-1
         idigit = mod(nstp,10)
         fname(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp/10
         enddo    
c
      open(unit=50,file=fname,status='unknown',form='formatted')
c
c     # compute depth average 
      do i=1,mx
         uave = 0.0d0
         have = 0.0d0
c
         do j=1,my
c           # radial velocity
            uave = uave+(1.0d0-q(i,j,6))*dy*
     &             q(i,j,3)/(q(i,j,1)+q(i,j,2))
c
c           # water height
            have = have+(1.0d0-q(i,j,6))*dy
            enddo
c
         if (dabs(uave) .lt. 1.d-40) uave = 0.0d0
c
         xc = xlower+dble(i-1)*dx+0.5d0*dx
         write(50,1005) xc,uave/have,have        
         enddo           
c
      close(unit=50)
      return
c
 1005 format(3e18.10)
      end      

