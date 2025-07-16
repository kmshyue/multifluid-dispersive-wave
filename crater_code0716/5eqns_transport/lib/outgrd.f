c
c ------------------------------------------------------------
c
      subroutine outgrd(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &           dx,dy,iframe,aux,maux)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      character*10 fname
c
c     # output physical grid system
c
      if (iframe .gt. 0) return
c
c     # cell center
      fname = 'fort.gxxxx'
      nstp = iframe
      do ipos=10,7,-1
         idigit = mod(nstp,10)
         fname(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp/10
         enddo
c
      open(unit=50,file=fname,status='unknown',form='formatted')
c
      write(50,*) mx,my
c
      do j=1,my
         yc = ylower+dble(j-1)*dy+0.5d0*dy
         do i=1,mx
            xc = xlower+dble(i-1)*dx+0.5d0*dx
c
            write(50,1001) xc,yc
            enddo
         enddo
c
      close(unit=50)
      return
c
 1001 format(2f26.12)
      end
