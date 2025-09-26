c
c -----------------------------------------------------------
c
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,
     &           dx,dy,maux,aux)
      implicit double precision (a-h,o-z)
      double precision kappa_max,kappa_min,kappa_avg
      double precision normals(2,2)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension quad(0:1,0:1,2),slength(2)
c
c     # set underlying uniform grid (mapped or not)
c
      mcapa = 7

      if (maux .lt. mcapa) then
         write(6,*) 'setaux2 : maux must be >= 7 when using '
         stop
      endif
c
      do j=1-mbc,my+mbc
         ye = ylower+dble(j-1)*dy
         do i=1-mbc,mx+mbc
            xe = xlower+dble(i-1)*dx
c
            do icell=0,1
               xcorner = xe+dble(icell)*dx
               do jcell= 0,1
                  ycorner = ye+dble(jcell)*dy
                  call mapc2p(xcorner,ycorner,xp,yp)
                  quad(icell,jcell,1) = xp
                  quad(icell,jcell,2) = yp
               enddo
            enddo

c           # call same routine as above, this time only using
c           # normals information
            call compute_info2(quad,normals,slength,area,
     &           xc,yc)

c          # Store information in aux array.  Of course, this can be
c          # stored in any arrangement desired.
c          # Normal at x face - first row in normals matrix
           aux(i,j,1) = normals(1,1)
           aux(i,j,2) = normals(1,2)
           aux(i,j,3) = slength(1)/dy

c          # Normal at y face - second row in normals matrix.
           aux(i,j,4) = normals(2,1)
           aux(i,j,5) = normals(2,2)
           aux(i,j,6) = slength(2)/dx

c          # Capacity
           if (area .le. 0.d0) then
              write(6,*) 'setquadinfo : capacity <= 0'
              write(6,*) 'area = ', area
              write(6,*) 'i,j = ', i,j
              write(6,'(A,A,A,A)') 'This usually means that there ',
     &              'is a problem with your mapping function. ',
     &              'Check to be sure your mapping is valid in ',
     &              'ghost cell regions. '
           endif
c
           aux(i,j,mcapa) = dabs(area)/(dx*dy)
c
c           aux(i,j,8) = xc/area
c           aux(i,j,9) = yc/area
c           aux(i,j,10) = quad(0,0,1)     
c           aux(i,j,11) = quad(0,0,2)    
        enddo
      enddo
c
      kappa_max = 0
      kappa_min = 100
      kappa_avg = 0;
      do j=1,my
         do i=1,mx
            kappa_max = max(aux(i,j,mcapa),kappa_max)
            kappa_min = min(aux(i,j,mcapa),kappa_min)
            kappa_avg = kappa_avg+aux(i,j,mcapa)
            enddo
         enddo
c
      write(6,'(A,F16.8)') 'Max kappa : ', kappa_max
      write(6,'(A,F16.8)') 'Min kappa : ', kappa_min
      write(6,'(A,F16.8)') 'Avg kappa : ', kappa_avg/(mx*my)
      write(6,'(A,F16.8)') 'Ratio     : ', kappa_max/kappa_min
      write(6,*) ' '
      return
c
      end
c
c -------------------------------------------------------------------
c
      subroutine compute_info2(quad,normals,slength,area,xc,yc)
      implicit double precision (a-h,o-z)
      double precision normals(2,2)
      dimension quad(0:1,0:1,2),slength(2)
      dimension v(2),w(2)
      dimension xpcorn(5),ypcorn(5)
      dimension subxc(3),subyc(3),subar(3)
c
c     # This is where the area, length and normals are computed
c
c     # compute normals to left edge
c     # \xi constant
c
c     # x_\eta
      v(1) = quad(0,1,1)-quad(0,0,1)
c
c     # y_\eta
      v(2) = quad(0,1,2)-quad(0,0,2)
c
      vnorm = dsqrt(v(1)*v(1)+v(2)*v(2))
      slength(1) = vnorm
      if (vnorm .eq. 0.d0) then
          normals(1,1) = 1.d0
          normals(1,2) = 0.d0
      else
c         # normal should point in positive x dir.
c
          normals(1,1) = v(2)/vnorm  
          normals(1,2) = -v(1)/vnorm
      endif
c
c     # compute normals to bottom edge
c     # \eta constant
c
c     # x_\xi
      w(1) = quad(1,0,1)-quad(0,0,1)
c
c     # y_\xi
      w(2) = quad(1,0,2)-quad(0,0,2)
c
      wnorm = dsqrt(w(1)*w(1)+w(2)*w(2))
      slength(2) = wnorm
      if (wnorm .eq. 0.d0) then
          normals(2,1) = 1.d0
          normals(2,2) = 0.d0
      else
c         # normal should point in positive y dir.
c
          normals(2,1) = -w(2)/wnorm 
          normals(2,2) = w(1)/wnorm
      endif
c
      xpcorn(1) = quad(0,0,1)
      xpcorn(2) = quad(0,1,1)
      xpcorn(3) = quad(1,1,1)
      xpcorn(4) = quad(1,0,1)
      xpcorn(5) = quad(0,0,1)
c
      ypcorn(1) = quad(0,0,2)
      ypcorn(2) = quad(0,1,2)
      ypcorn(3) = quad(1,1,2)
      ypcorn(4) = quad(1,0,2)
      ypcorn(5) = quad(0,0,2)
c
      xpcorn(5) = xpcorn(1)
      ypcorn(5) = ypcorn(1)
      do ip=3,5
         subxc(ip-2) = (xpcorn(ip-1)+xpcorn(ip)+xpcorn(1))/3.d0
         subyc(ip-2) = (ypcorn(ip-1)+ypcorn(ip)+ypcorn(1))/3.d0
         subar(ip-2) = 0.5d0*((ypcorn(1)+ypcorn(ip-1))*
     &                        (xpcorn(ip-1)-xpcorn(1))+
     &                        (ypcorn(ip-1)+ypcorn(ip))*
     &                        (xpcorn(ip)-xpcorn(ip-1))+
     &                        (ypcorn(ip)+ypcorn(1))*
     &                        (xpcorn(1)-xpcorn(ip)))
         enddo
c
      area = 0.d0
      xc   = 0.d0
      yc   = 0.d0
      do ip=1,3
         area = area+subar(ip)
         xc   = xc+subar(ip)*subxc(ip)
         yc   = yc+subar(ip)*subyc(ip)
         enddo
      return
c
      end
c
c -------------------------------------------------------------
c
      subroutine get_aux_locations_n(ixy,mcapa,locrot,locarea)
      implicit double precision (a-h,o-z)
c
c     # Get aux locations for normal solve
      mcapa   = 7
      locrot  = (ixy-1)*3+1
      locarea = locrot+2
      return
c
      end
c
c -----------------------------------------------------------------------
c
      subroutine get_aux_locations_t(ixy,mcapa,locrot,locarea)
      implicit double precision (a-h,o-z)

c     # Get aux locations for transverse solves.

c     # Get transverse direction
      icoor = 3-ixy
c
      call get_aux_locations_n(icoor,mcapa,locrot,locarea)
      return
c
      end
c
c --------------------------------------------------------------
c
      subroutine compute_tangent(rot)
      implicit double precision (a-h,o-z)
      dimension rot(4)
c
      rot(3) = -rot(2)
      rot(4) = rot(1)
      return
c
      end
c
c ---------------------------------------------------------------------
c
      subroutine rotate2(rot,velcomps)
      implicit double precision (a-h,o-z)
      dimension rot(4),velcomps(2)
c
      v1 = velcomps(1)
      v2 = velcomps(2)
c
      velcomps(1) = rot(1)*v1+rot(2)*v2
      velcomps(2) = rot(3)*v1+rot(4)*v2
      return
c
      end
c
c ----------------------------------------------------------
c
      subroutine rotate2_tr(rot,velcomps)
      implicit double precision (a-h,o-z)
      dimension velcomps(2),rot(4)
c
      v1 = velcomps(1)
      v2 = velcomps(2)
c
      velcomps(1) = rot(1)*v1+rot(3)*v2
      velcomps(2) = rot(2)*v1+rot(4)*v2
      return
c
      end
